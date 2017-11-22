/**
 * program: dadafits
 *          Written for the AA-Alert project, ASTRON
 *
 * Purpose: connect to a ring buffer and create FITS output per TAB
 *          Depending on science case and mode, reduce time and frequency resolution to 1 bit
 *          Fits files are created using templates
 *
 * Science case 3, mode 0:
 *          template: sc34_1bit_I_reduced.txt
 *
 *          A ringbuffer page is interpreted as an array of Stokes I:
 *          [NTABS, NCHANNELS, padded_size] = [12, 1536, > 12500]
 *
 *          The code reduces (by summation) from 12500 to 500 timesteps
 *          and from 1536 to 384 channels.
 *          Time dimension padding is required by other programes (GPU pipeline)
 *          that connects to the same ringbuffer.
 *
 * Science case 3, mode 1:
 *          template: sc3_8bit_IQUV_full.txt
 *
 *          A ringbuffer page is interpreted as an interleaved array of Stokes IQUV:
 *          [tab][channel_offset][sequence_number][packet]
 *
 *          tab             := ranges from 0 to NTABS
 *          channel_offset  := ranges from 0 to NCHANNELS/4 (384)
 *          sequence_number := ranges from 0 to 25
 *
 *          with a packet: [t0 .. t499][c0 .. c3][IQUV] total of 500*4*4=8000 bytes
 *          t = tn + sequence_number * 25
 *          c = cn + channel_offset * 4
 *
 * Science case 3, mode 2:
 *          template: sc34_1bit_I_reduced.txt
 *
 *          A ringbuffer page is interpreted as an array of Stokes I:
 *          [NTABS, NCHANNELS, padded_size] = [1, 1536, > 12500]
 *
 *          The code reduces (by summation) from 12500 to 500 timesteps
 *          and from 1536 to 384 channels.
 *          Time dimension padding is required by other programes (GPU pipeline)
 *          that connects to the same ringbuffer.
 *
 * Science case 4, mode 0:
 *          template: sc34_1bit_I_reduced.txt
 *
 *          A ringbuffer page is interpreted as an array of Stokes I:
 *          [NTABS, NCHANNELS, padded_size] = [12, 1536, > 25000]
 *
 *          The code reduces (by summation) from 25000 to 500 timesteps
 *          and from 1536 to 384 channels.
 *          Time dimension padding is required by other programes (GPU pipeline)
 *          that connects to the same ringbuffer.
 *
 * Science case 4, mode 1:
 *          template: sc4_8bit_IQUV_full.txt
 *
 *          A ringbuffer page is interpreted as an interleaved array of Stokes IQUV:
 *          [tab][channel_offset][sequence_number][packet]
 *
 *          tab             := ranges from 0 to NTABS
 *          channel_offset  := ranges from 0 to NCHANNELS/4 (384)
 *          sequence_number := ranges from 0 to 50
 *
 *          with a packet: [t0 .. t499][c0 .. c3][IQUV] total of 500*4*4=8000 bytes
 *          t = tn + sequence_number * 50
 *          c = cn + channel_offset * 4
 *
 * Science case 4, mode 0:
 *          template: sc34_1bit_I_reduced.txt
 *
 *          A ringbuffer page is interpreted as an array of Stokes I:
 *          [NTABS, NCHANNELS, padded_size] = [1, 1536, > 25000]
 *
 *          The code reduces (by summation) from 25000 to 500 timesteps
 *          and from 1536 to 384 channels.
 *          Time dimension padding is required by other programes (GPU pipeline)
 *          that connects to the same ringbuffer.
 *
 * Author: Jisk Attema, Netherlands eScience Center
 * Licencse: Apache v2.0
 */

#include <string.h>
#include <stdio.h>
#include <fitsio.h>
#include <getopt.h>
#include <fitsio.h>
#include <math.h>
#include <signal.h>
#include <fenv.h>
#include <errno.h>

#include "dada_hdu.h"
#include "ascii_header.h"

FILE *runlog = NULL;
#define LOG(...) {fprintf(stdout, __VA_ARGS__); fprintf(runlog, __VA_ARGS__); fflush(stdout);}

#define NTABS_MAX 12
#define NCHANNELS 1536
#define NPOLS 4

// 80 microsecond -> 2 milisecond
#define SC3_NTIMES 12500
#define SC3_DOWNSAMPLE_TIME 25

// 40 microsecond -> 2 milisecond
#define SC4_NTIMES 25000
#define SC4_DOWNSAMPLE_TIME 50

// the same for both science case 3 and 4
#define NTIMES_LOW 500
#define NCHANNELS_LOW (NCHANNELS / 4)

// Variables read from ring buffer header
float min_frequency = 1492;
float channel_bandwidth = 0.1953125;

// Variables used for FITS output
fitsfile *output[NTABS_MAX];

// Fit column ID's, looked up from the template.
// If not present in the template, it is set to '-1' and not written
int col_data = -1;
int col_freqs = -1;
int col_offset = -1;
int col_scale = -1;
int col_weights = -1;
int col_offs_sub = -1;

unsigned int downsampled[NCHANNELS_LOW * NTIMES_LOW];
unsigned char packed[NCHANNELS_LOW * NTIMES_LOW / 8];
unsigned char transposed[NTABS_MAX * NCHANNELS * NPOLS * SC4_NTIMES]; // big enough for largest case
float offset[NCHANNELS];
float scale[NCHANNELS];
float weights[NCHANNELS];
float freqs[NCHANNELS];

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
 * Downsample timeseries by summation over time and frequency
 *
 * Not optimized, but such large sums are memory bound, not CPU bound
 * On modern laptop (dell xps13) this program stays below a few percent CPU usage.
 *
 * Data is copied from the buffer to the global 'downsampled' array
 *
 * Also, we are using 8bit integers, which are much faster than SSE/AVX operations on floats.
 * See for some interesting reading https://github.com/jodavies/dot-product
 *
 * @param {int} tab                                     TAB index
 * @param {uchar[NTABS, NCHANNELS, padded_size]} buffer Ring buffer page to downsample
 * @param {int} padded_size                             Size of fastest dimension, as timeseries are padded for optimal memory layout on GPU
 * @param {int} science_case                            Science case; either 3 (12500 samples) or 4 (25000 samples) per batch of 1.024 seconds
 */
void downsample_sc3(const int tab, const unsigned char *buffer, const int padded_size) {
  unsigned int *temp1 = downsampled;
  int dc; // downsampled channel
  int dt; // downsampled time
  int t; // full time

  for (dc=0; dc < NCHANNELS_LOW; dc++) {
    // pointer to next sample in the four channels
    unsigned const char *s0 = &buffer[tab * NCHANNELS * padded_size + ((dc << 2) + 0) * padded_size];
    unsigned const char *s1 = &buffer[tab * NCHANNELS * padded_size + ((dc << 2) + 1) * padded_size];
    unsigned const char *s2 = &buffer[tab * NCHANNELS * padded_size + ((dc << 2) + 2) * padded_size];
    unsigned const char *s3 = &buffer[tab * NCHANNELS * padded_size + ((dc << 2) + 3) * padded_size];

    for (dt=0; dt < NTIMES_LOW; dt++) {
      // partial sums (per channel)
      unsigned int ps0 = 0;
      unsigned int ps1 = 0;
      unsigned int ps2 = 0;
      unsigned int ps3 = 0;

      for (t=0; t < SC3_DOWNSAMPLE_TIME; t++) {
        ps0 += *s0++;
        ps1 += *s1++;
        ps2 += *s2++;
        ps3 += *s3++;
      }
      *temp1++ = ps0 + ps1 + ps2 + ps3;
    }
  }
}

void downsample_sc4(const int tab, const unsigned char *buffer, const int padded_size) {
  unsigned int *temp1 = downsampled;
  int dc; // downsampled channel
  int dt; // downsampled time
  int t; // full time

  for (dc=0; dc < NCHANNELS_LOW; dc++) {
    // pointer to next sample in the four channels
    unsigned const char *s0 = &buffer[tab * NCHANNELS * padded_size + ((dc << 2) + 0) * padded_size];
    unsigned const char *s1 = &buffer[tab * NCHANNELS * padded_size + ((dc << 2) + 1) * padded_size];
    unsigned const char *s2 = &buffer[tab * NCHANNELS * padded_size + ((dc << 2) + 2) * padded_size];
    unsigned const char *s3 = &buffer[tab * NCHANNELS * padded_size + ((dc << 2) + 3) * padded_size];

    for (dt=0; dt < NTIMES_LOW; dt++) {
      // partial sums (per channel)
      unsigned int ps0 = 0;
      unsigned int ps1 = 0;
      unsigned int ps2 = 0;
      unsigned int ps3 = 0;

      for (t=0; t < SC4_DOWNSAMPLE_TIME; t++) {
        ps0 += *s0++;
        ps1 += *s1++;
        ps2 += *s2++;
        ps3 += *s3++;
      }
      *temp1++ = ps0 + ps1 + ps2 + ps3;
    }
  }
}

/**
 * Pack series of 8-bit StokesI to 1-bit
 *   NBIN*NCHAN*NPOL*NSBLK => 1 x 384 x 1 x 500 bits equals or 24000 bytes
 */
void pack_sc34() {
  unsigned int *temp1;
  unsigned char *temp2;

  // DEBUG NaNs:
  errno = 0;
  feclearexcept(FE_ALL_EXCEPT);

  int dc;
  for (dc = 0; dc < NCHANNELS_LOW; dc++) {

    // First pass: calculate average(=offset) and stdev(=scale)
    temp1 = &downsampled[dc * NTIMES_LOW];

    // Sum:
    // total sum is over 25000 samples times 4 frequencies is 100,000 samples (50,000 for sc3)
    // maxium value is 100000 * 255 = 5,500,000
    // an unsigned int should hold 4,294,967,295, so no overflow here
    //
    // Sos:
    // maxium value of a downsampled sample is 255 * 50 * 4 = 51,000
    // of these, we are summing 500 squares, so maximum value is
    //             51000 * 51000 * 500 = 1,300,500,000,000
    // unsigned long long holds 18,446,744,073,709,551,615, so ok.
    //
    unsigned int sum = 0;
    unsigned long long sos = 0;

    int dt;
    for (dt=0; dt < NTIMES_LOW; dt++) {
      sos += (*temp1) * (*temp1);
      sum += *temp1++;
    }
    offset[dc] = sum / (1.0 * NTIMES_LOW);
    scale[dc]  = sqrtf(sos / (1.0 * NTIMES_LOW) - offset[dc] * offset[dc]);

    // Set cutoff to 1 stdev above average
    unsigned int cutoff = offset[dc] + scale[dc];

    // Second pass: convert to 1 bit
    temp1 = &downsampled[dc * NTIMES_LOW];
    for (dt=0; dt < NTIMES_LOW; dt++) {
      *temp1++ = *temp1 > cutoff ? 1 : 0;
    }
  }

  // Third pass: pack bits in bytes
  temp1 = downsampled;
  temp2 = packed;
  int count;
  for (count=0; count < NCHANNELS_LOW * NTIMES_LOW / 8; count++) {
    *temp2  = *temp1++ ? 1 << 7 : 0;
    *temp2 += *temp1++ ? 1 << 6 : 0;
    *temp2 += *temp1++ ? 1 << 5 : 0;
    *temp2 += *temp1++ ? 1 << 4 : 0;
    *temp2 += *temp1++ ? 1 << 3 : 0;
    *temp2 += *temp1++ ? 1 << 2 : 0;
    *temp2 += *temp1++ ? 1 << 1 : 0;
    *temp2 += *temp1++ ? 1      : 0;
    temp2++;
  }

  // DEBUG NaNs:
  int except = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  if (errno || except) {
    LOG("Error in packing data: errn=%i, fetestexcept=", errno);
    switch (except) {
      case FE_INVALID: printf("FE_INVALID\n"); break;
      case FE_DIVBYZERO: printf("FE_DIVBYZERO\n"); break;
      case FE_OVERFLOW: printf("FE_OVERFLOW\n"); break;
      case FE_UNDERFLOW: printf("FE_UNDERFLOW\n"); break;
      default: printf("-\n");
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
 * Initialize the CFITSIO library
 * @param {char *} template_file    FITS template to use for creating initial file.
 * @param {char *} output_directory Directoy where output FITS files can be written.
 * @param {int} ntabs               Number of beams
 * @param {float } bandwidth        Bandwith per channel, after optional downsampling
 */
void dadafits_fits_init (char *template_file, char *output_directory, int ntabs, float bandwidth) {
  int status;
  float version;
  fits_get_version(&version);
  LOG("Using FITS library version %f\n", version);

  int t;
  for (t=0; t<ntabs; t++) {
    char fname[256];
    fitsfile *fptr;

    if (output_directory) {
      snprintf(fname, 256, "%s/beam%c.fits(%s)", output_directory, 'A'+t, template_file);
    } else {
      snprintf(fname, 256, "beam%c.fits(%s)", 'A'+t, template_file);
    }
    LOG("Writing tab %02i to file %s\n", t, fname);

    status = 0; if (fits_create_file(&fptr, fname, &status)) fits_error_and_exit(status);
    status = 0; if (fits_movabs_hdu(fptr, 1, NULL, &status)) fits_error_and_exit(status);
    status = 0; if (fits_write_date(fptr, &status))          fits_error_and_exit(status);
    status = 0; if (fits_write_chksum(fptr, &status))        fits_error_and_exit(status);

    status = 0; if (fits_movabs_hdu(fptr, 2, NULL, &status)) fits_error_and_exit(status);

    output[t] = fptr;
  }

  // Set scaling, weights, and offsets to neutral values
  int i;
  for (i=0; i<NCHANNELS; i++) {
    offset[i] = 0.0;
    scale[i] = 1.0;
    weights[i] = 1.0;
    freqs[i] = min_frequency + i * bandwidth;
  }

  // Get required column ID's
  col_freqs = dadafits_find_column("DAT_FREQ", output[0]);
  col_data = dadafits_find_column("DATA", output[0]);
  col_offset = dadafits_find_column("DAT_OFFS", output[0]);
  col_scale = dadafits_find_column("DAT_SCL", output[0]);
  col_weights = dadafits_find_column("DAT_WTS", output[0]);
  col_offs_sub = dadafits_find_column("OFFS_SUB", output[0]);
}

/**
 * Write a row of 8-bits Stokes IQUV data a FITS BINTABLE.SUBINT
 *
 * Uses the array 'transposed'
 *
 * @param {int} tab                Tight array beam index used to select output file
 * @param {long} rowid             Row number in the SUBINT table, corresponds to ringbuffer page number + 1
 * @param {int} rowlength          Size of a data row
 */
void write_fits_IQUV(int tab, long rowid, int rowlength) {
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
    if (fits_write_col(fptr, TFLOAT, col_freqs, rowid, 1, NCHANNELS, freqs, &status)) {
      fits_error_and_exit(status);
    }
  }

  status = 0;
  if (fits_write_col(fptr, TBYTE, col_data, rowid, 1, rowlength, transposed, &status)) {
    fits_error_and_exit(status);
  }
}

/**
 * Write a row of packed, 1 bit Stokes I to a FITS BINTABLE.SUBINT
 *
 * Uses the global arrays: 'scale', 'offset', and 'packed'
 *
 * @param {int} tab                Tight array beam index used to select output file
 * @param {int} rowid              Row number in the SUBINT table, corresponds to ringbuffer page number + 1
 */
void write_fits_packedI(int tab, long rowid) {
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
    if (fits_write_col(fptr, TFLOAT, col_freqs, rowid, 1, NCHANNELS_LOW, freqs, &status)) {
      fits_error_and_exit(status);
    }
  }

  if (col_weights >= 0) {
    status = 0;
    if (fits_write_col(fptr, TFLOAT, col_weights, rowid, 1, NCHANNELS_LOW, weights, &status)) {
      fits_error_and_exit(status);
    }
  }

  if (col_offset >= 0) {
    status = 0;
    if (fits_write_col(fptr, TFLOAT, col_offset, rowid, 1, NCHANNELS_LOW, offset, &status)) {
      fits_error_and_exit(status);
    }
  }

  if (col_scale >= 0) {
    status = 0;
    if (fits_write_col(fptr, TFLOAT, col_scale, rowid, 1, NCHANNELS_LOW, scale, &status)) {
      fits_error_and_exit(status);
    }
  }

  status = 0;
  if (fits_write_col(fptr, TBYTE,  col_data, rowid, 1, NCHANNELS_LOW * NTIMES_LOW / 8, packed, &status)) {
    fits_error_and_exit(status);
  }
}

/**
 * Open a connection to the ringbuffer
 *
 * @param {char *} key String containing the shared memory key as hexadecimal number
 * @returns {hdu *} A connected HDU
 */
dada_hdu_t *init_ringbuffer(char *key) {
  uint64_t nbufs;

  multilog_t* multilog = NULL; // TODO: See if this is used in anyway by dada

  // create hdu
  dada_hdu_t *hdu = dada_hdu_create (multilog);

  // init key
  key_t shmkey;
  sscanf(key, "%x", &shmkey);
  dada_hdu_set_key(hdu, shmkey);
  LOG("dadafits SHMKEY: %s\n", key);

  // connect
  if (dada_hdu_connect (hdu) < 0) {
    LOG("ERROR in dada_hdu_connect\n");
    exit(EXIT_FAILURE);
  }

  // Make data buffers readable
  if (dada_hdu_lock_read(hdu) < 0) {
    LOG("ERROR in dada_hdu_open_view\n");
    exit(EXIT_FAILURE);
  }

  // get write address
  char *header;
  uint64_t bufsz;
  LOG("dadafits reading header");
  header = ipcbuf_get_next_read (hdu->header_block, &bufsz);
  if (! header || ! bufsz) {
    LOG("ERROR. Get next header block error\n");
    exit(EXIT_FAILURE);
  }

  // parse header
  unsigned int uintValue;
  float floatValue[2];
  ascii_header_get(header, "MIN_FREQUENCY", "%f", &min_frequency);
  ascii_header_get(header, "CHANNEL_BANDWIDTH", "%f", &channel_bandwidth);

  // tell the ringbuffer the header has been read
  if (ipcbuf_mark_cleared(hdu->header_block) < 0) {
    LOG("ERROR. Cannot mark the header as cleared\n");
    exit(EXIT_FAILURE);
  }

  LOG("psrdada HEADER:\n%s\n", header);

  return hdu;
}

/**
 * Print commandline options
 */
void printOptions() {
  printf("usage: dadafits -k <hexadecimal key> -l <logfile> -c <science_case> -m <science_mode> -b <padded_size> -t <template> -d <output_directory>\n");
  printf("e.g. dadafits -k dada -l log.txt -c 3 -m 0 -b 25088 -t /full/path/template.txt -d /output/directory\n");
  return;
}

/**
 * Parse commandline
 */
void parseOptions(int argc, char *argv[], char **key, int *padded_size, char **logfile, int *science_case, int *science_mode, char **template_file, char **output_directory) {
  int c;

  int setk=0, setb=0, setl=0, setc=0, setm=0, sett=0, setd=0;
  while((c=getopt(argc,argv,"k:b:l:c:m:t:d:"))!=-1) {
    switch(c) {
      // -t <template_file>
      case('t'):
        *template_file = strdup(optarg);
        sett=1;
        break;

      // OPTIONAL: -d <output_directory>
      // DEFAULT: CWD
      case('d'):
        *output_directory = strdup(optarg);
        setd=1;
        break;

      // -k <hexadecimal_key>
      case('k'):
        *key = strdup(optarg);
        setk=1;
        break;

      // -b padded_size (bytes)
      case('b'):
        *padded_size = atoi(optarg);
        setb=1;
        break;

      // -l log file
      case('l'):
        *logfile = strdup(optarg);
        setl=1;
        break;

      // -c science_case
      case('c'):
        *science_case = atoi(optarg);
        setc = 1;
        break;

      // -m science_mode
      case('m'):
        *science_mode = atoi(optarg);
        setm = 1;
        break;

      default:
        printOptions();
        exit(0);
    }
  }

  // Required arguments
  if (!setk || !setl || !setb || !setc || !setm || !sett) {
    printOptions();
    exit(EXIT_FAILURE);
  }
}
/**
 * Deinterleave (transpose) an IQUV ring buffer page to the ordering needed for FITS files
 * Note that this is probably a slow function, and is not meant to be run real-time
 * Suggested use is:
 *   1. realtime: ringbuffer -> [trigger] -> dada_dbdisk
 *   2. offline: dada_dbdisk -> ringbuffer -> dadafits
 *
 *  @param {const unsigned char *} page    Ringbuffer page with interleaved data
 *  @param {int} ntabs                     Number of tabs
 *  @param {int} sequence_length           Number of packets per
 */
void deinterleave (const unsigned char *page, int ntabs, int sequence_length) {
  // ring buffer page contains matrix:
  //   [tab][channel_offset][sequence_number][8000]
  //
  // tab             : ranges from 0 to (1 or 12) mode TAB / IAB
  // channel_offset  : ranges from 0 to NCHANNELS/4 (=1536/4=384)
  // sequence_number : ranges from 0 to 25 or 50 (sc3 or sc4)
  //
  // the 8000 bytes are the original packets, containing stokes IQUV:
  //   [t0 .. t499][c0 .. c3][the 4 components IQUV]
  //
  //   t0, .., t499 = sequence_number * 500 + tx
  //   c0, c1, c2, c3 = curr_channel + 0, 1, 2, 3
  //
  // Transposed buffer will contain:
  // (NBIN,NCHAN,NPOL,NSBLK) = (1,1536,4,12500) or (1,1536,4,25000); sc3 or sc4

  // Tranpose by linearly processing original packets from the page
  const unsigned char *packet = page;

  // and find the matching address in the transposed buffer
  int tab = 0;
  for (tab = 0; tab < ntabs; tab++) {
    int channel_offset = 0;
    for (channel_offset = 0; channel_offset < NCHANNELS; channel_offset+=4) {
      int sequence_number = 0;
      for (sequence_number = 0; sequence_number < sequence_length; sequence_number++) {
        // process packet
        int tn,cn,pn;
        for (tn = 0; tn < 500; tn++) { // 500 samples per packet
          for (cn = 0; cn < 4; cn++) { // 4 channels per packet
            for (pn = 0; pn < NPOLS; pn++) {
              transposed[
                ((tab * NCHANNELS +
                  cn + channel_offset) * NPOLS +
                 pn) * sequence_length * 500 +
                  tn + sequence_number * 500
              ] = *packet++;
            }
          }
        }
      }
    }
  }
}

int main (int argc, char *argv[]) {
  char *key;
  int padded_size;
  char *logfile;
  char *template_file;
  char *output_directory = NULL; // defaults to CWD
  int science_case;
  int science_mode;
  int ntabs;
  int nchannels;
  int sequence_length;

  // parse commandline
  parseOptions(argc, argv, &key, &padded_size, &logfile, &science_case, &science_mode, &template_file, &output_directory);

  // set up logging
  if (logfile) {
    runlog = fopen(logfile, "w");
    if (! runlog) {
      LOG("ERROR opening logfile: %s\n", logfile);
      exit(EXIT_FAILURE);
    }
    LOG("Logging to logfile: %s\n", logfile);
    free (logfile);
  }

  LOG("dadafits version: " VERSION "\n");

  switch (science_mode) {
    case 0: // IQUV + TAB to deinterleave
      nchannels = NCHANNELS_LOW;
      ntabs = 12;
      break;
    case 1: // I + TAB to be compressed and downsampled
      nchannels = NCHANNELS;
      ntabs = 12;
      break;
    case 2: // I + IAB to be compressed and downsampled
      nchannels = NCHANNELS_LOW;
      ntabs = 1;
      break;
    case 3: // IQUV + IAB to deinterleave
      nchannels = NCHANNELS;
      ntabs = 1;
      break;
    default:
      fprintf(stderr, "Illegal science mode %i\n", science_mode);
      exit(EXIT_FAILURE);
  }
  switch (science_case) {
    case 3:
      sequence_length = 25;
      if (padded_size < 12500)
        fprintf(stderr, "Error: padded_size too small, should be at least 12500 for science case 3\n");
      break;
    case 4:
      sequence_length = 50;
      if (padded_size < 25000)
        fprintf(stderr, "Error: padded_size too small, should be at least 25000 for science case 4\n");
      break;
    default:
      fprintf(stderr, "Illegal science case %i\n", science_case);
  }

  int quit = 0;
  long page_count = 0;

  // Trap Ctr-C to properly close fits files on exit
  signal(SIGTERM, fits_error_and_exit);

#ifdef DRY_RUN
  // do 10 iterations with random data, ignore ringbuffer
  int mysize = NTABS_MAX * NCHANNELS * padded_size * 4; // IQUV
  char *page = malloc(mysize);

  dadafits_fits_init(template_file, output_directory, ntabs,
      channel_bandwidth * NCHANNELS / nchannels);

  int g_seed = 1234;
  inline unsigned char fastrand() {
    g_seed = (214013*g_seed+2531011);
    // return (g_seed>>16)&0x7FFF;
    return (g_seed>>16)&0xFF;
  }

  while(page_count < 10) {
    int kkk;
    for(kkk=0; kkk<mysize; kkk++) {
      page[kkk] = fastrand();
    }
#else
  // normal operation with ringbuffer
  char *page = NULL;

  // must init ringbuffer before fits, as this reads parameters
  // like channel_bandwidth from ring buffer header
  dada_hdu_t *ringbuffer = init_ringbuffer(key);
  ipcbuf_t *data_block = (ipcbuf_t *) ringbuffer->data_block;
  ipcio_t *ipc = ringbuffer->data_block;
  uint64_t bufsz = ipc->curbufsz;

  dadafits_fits_init(template_file, output_directory, ntabs,
      channel_bandwidth * NCHANNELS / nchannels);

  while(!quit && !ipcbuf_eod(data_block)) {
    page = ipcbuf_get_next_read(data_block, &bufsz);
#endif

    int tab; // tight array beam

    if (! page) {
      quit = 1;
    } else {
      switch (science_mode) {
        // stokesI data to compress, downsample, and write
        case 0:
        case 2:
          for (tab = 0; tab < ntabs; tab++) {
            // move data from the page to the downsampled array
            if (science_case == 3) {
              downsample_sc3(tab, page, padded_size);
            } else if (science_case == 4) {
              downsample_sc4(tab, page, padded_size);
            } else {
              exit(EXIT_FAILURE);
            }
            // pack data from the downsampled array to the packed array,
            // and set scale and offset arrays with used values
            pack_sc34();

            // write data from the packed, scale, and offset arrays
            write_fits_packedI(tab, page_count + 1); // pagecount starts at 0, but FITS rowid at 1
          }
        break;

        // stokesIQUV data to write
        case 1:
        case 3:
          // transpose data from page to transposed buffer
          deinterleave(page, ntabs, sequence_length);

          for (tab = 0; tab < ntabs; tab++) {
            // write data from transposed buffer
            write_fits_IQUV(tab,
                page_count + 1, // pagecount starts at 0, but FITS rowid at 1
                sequence_length * 500 * 4 * NCHANNELS);
          }
          break;

        default:
          // should not happen
          fprintf(stderr, "Illegal science mode %i\n", science_mode);
          quit = 1;
          break;
      }
#ifndef DRY_RUN
      ipcbuf_mark_cleared((ipcbuf_t *) ipc);
#endif
      page_count++;
    }
  }

#ifndef DRY_RUN
  if (ipcbuf_eod(data_block)) {
    LOG("End of data received\n");
  }

  dada_hdu_unlock_read(ringbuffer);
  dada_hdu_disconnect(ringbuffer);
#endif

  LOG("Read %li pages\n", page_count);

  close_fits();
}
