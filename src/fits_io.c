#include <fitsio.h>
#include "dadafits_internal.h"

fitsfile *output[NSYNS_MAX];

float fits_offset[NCHANNELS * NPOLS];
float fits_scale[NCHANNELS * NPOLS];
float fits_weights[NCHANNELS];
float fits_freqs[NCHANNELS];

// Fit column ID's, looked up from the template.
// If not present in the template, it set it to '-1' and it is not written
int col_data = 17;
int col_freqs = 13;
int col_offset = 15;
int col_scale = 16;
int col_weights = 14;
int col_offs_sub = 2;

/**
 * pretty print the fits error to the log, and close down cleanly
 * @param {int} status The status code returned by the (failed) fits call
 */
void fits_error_and_exit(int status) {
  if (runlog) {
    fits_report_error(runlog, status);
  }
  fits_report_error(stdout, status);

  close_fits();
  exit(EXIT_FAILURE);
}

/**
 * Close all opened fits files
 */
void close_fits() {
  int tab, status;

  LOG("Closing files\n");
  for (tab=0; tab<NTABS_MAX; tab++) {
    if (output[tab]) {
      // ignore errors on closing files; cfitsio 3.37 reports junk error codes
      fits_close_file(output[tab], &status);

      // FUTURE VERSION:
      // if (fits_close_file (output[tab], &status)) {
      //   if (runlog) fits_report_error(runlog, status);
      //   fits_report_error(stdout, status);
      // }
    }
  }
}

/**
 * Retrun the column id for the parameter with the given name,
 * -1 when not found
 */
int dadafits_find_column(char *name, fitsfile *file) {
  int column = -1;
  int status = 0;

  fits_get_colnum(file, 0, name, &column, &status);
  if (status != 0) {
    if (status == 219) {
      // not found
      column = -1;
    } else {
      // different error, abort
      fits_error_and_exit(status);
    }
  }
  LOG("%15s to column: %i\n", name, column);
  return column;
}

/**
 * Write a row of data to a FITS BINTABLE.SUBINT
 *
 * Optionally uses the global arrays: 'weights', 'scale', 'offset', and 'packed'
 *
 * @param {const int} tab                Tied array beam index used to select output file
 * @param {const int} channels           The number of channels to use
 * @param {const int} pols               The number of polarizations to use
 * @param {const int} rowid              Row number in the SUBINT table, corresponds to ringbuffer page number + 1
 * @param {const int} rowlength          Size of a data row
 * @param {const unsigned char *} data   Row to write
 */
void write_fits(const int tab, const int channels, const int pols, const long rowid, const int rowlength, unsigned char *data) {
  int status;
  fitsfile *fptr = output[tab];

  // From the cfitsio documentation:
  // Note that it is *not* necessary to insert rows in a table before writing data to those rows (indeed, it
  // would be inefficient to do so). Instead, one may simply write data to any row of the table, whether
  // that row of data already exists or not.
  //
  // status = 0;
  // if (fits_insert_rows(fptr, rowid, 1, &status)) {
  //   fits_error_and_exit(status);
  // }

  double offs_sub = (double) rowid * 1.024; // OFFS_SUB is in seconds since start of run, but may not be zero

  if (col_offs_sub >= 0) {
    status = 0;
    if (fits_write_col(fptr, TDOUBLE, col_offs_sub, rowid, 1, 1, &offs_sub, &status)) {
      fits_error_and_exit(status);
    }
  }

  if (col_freqs >= 0) {
    status = 0;
    if (fits_write_col(fptr, TFLOAT, col_freqs, rowid, 1, channels, fits_freqs, &status)) {
      fits_error_and_exit(status);
    }
  }

  if (col_weights >= 0) {
    status = 0;
    if (fits_write_col(fptr, TFLOAT, col_weights, rowid, 1, channels, fits_weights, &status)) {
      fits_error_and_exit(status);
    }
  }

  if (col_offset >= 0) {
    status = 0;
    if (fits_write_col(fptr, TFLOAT, col_offset, rowid, 1, channels * pols, fits_offset, &status)) {
      fits_error_and_exit(status);
    }
  }

  if (col_scale >= 0) {
    status = 0;
    if (fits_write_col(fptr, TFLOAT, col_scale, rowid, 1, channels * pols, fits_scale, &status)) {
      fits_error_and_exit(status);
    }
  }

  status = 0;
  if (fits_write_col(fptr, TBYTE,  col_data, rowid, 1, rowlength, data, &status)) {
    fits_error_and_exit(status);
  }
}

/**
 * Initialize the CFITSIO library
 * @param {char *} template_dir     Directory containing FITS templates
 * @param {char *} template_file    FITS template to use for creating initial file.
 * @param {char *} output_directory Directoy where output FITS files can be written.
 * @param {int} ntabs               Number of beams
 * @param {int} mode                0: one file per tab, 1: one file per selected synthesized beam
 * @param {float } min_frequency    Center of lowest frequency band of observation
 * @param {float } bandwidth        Bandwith per channel, after optional downsampling
 */
void dadafits_fits_init (const char *template_dir, const char *template_file, char *output_directory,
    int ntabs, int mode, float min_frequency, float bandwidth) {
  int status;
  float version;
  fits_get_version(&version);
  LOG("Using FITS library version %f\n", version);

  // set filename prefix according to mode:
  char *prefix;
  char *tab_prefix = "tab";
  char *synthesized_beam_prefix = "syn";
  if (mode == 0) {
    // one file per tab, name after tab
    prefix = tab_prefix;
  } else {
    // one file per synthesized beam, name after synthesize beam
    prefix = synthesized_beam_prefix;
  }

  int t;
  for (t=0; t<NSYNS_MAX; t++) {
    char fname[256];
    fitsfile *fptr;

    if (mode == 0 && t >= ntabs) {
      // when one file per tab, stop after ntab files
      break;
    } else if (mode == 1 && ! synthesized_beam_selected[t]) {
      // when one file per synthesized beam, skip when not selected
      continue;
    }

    if (output_directory) {
      snprintf(fname, 256, "%s/%s%c.fits(%s/%s)", output_directory, prefix, 'A'+t, template_dir, template_file);
    } else {
      snprintf(fname, 256, "%s%c.fits(%s/%s)", prefix, 'A'+t, template_dir, template_file);
    }
    LOG("Writing %s %02i to file %s\n", prefix, t, fname);

    status = 0; if (fits_create_file(&fptr, fname, &status)) fits_error_and_exit(status);
    status = 0; if (fits_movabs_hdu(fptr, 1, NULL, &status)) fits_error_and_exit(status);
    status = 0; if (fits_write_date(fptr, &status))          fits_error_and_exit(status);
    status = 0; if (fits_write_chksum(fptr, &status))        fits_error_and_exit(status);

    status = 0; if (fits_movabs_hdu(fptr, 2, NULL, &status)) fits_error_and_exit(status);

    output[t] = fptr;
  }

  // Set scaling, weights, and offsets to neutral values
  int i;
  for (i=0; i<NCHANNELS * NPOLS; i++) {
    fits_offset[i] = 0.0;
    fits_scale[i] = 1.0;
  }

  for (i=0; i<NCHANNELS; i++) {
    fits_weights[i] = 1.0;
    fits_freqs[i] = min_frequency + i * bandwidth;
  }
}
