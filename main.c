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
 * Science case 4, mode 2:
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

#include "dada_hdu.h"
#include "ascii_header.h"

FILE *runlog = NULL;
#define LOG(...) {fprintf(stdout, __VA_ARGS__); fprintf(runlog, __VA_ARGS__); fflush(stdout);}

#define NTABS_MAX 12
#define NCHANNELS 1536
#define NCHANNELS_LOW (NCHANNELS / 4)

// 80 microsecond -> 2 milisecond
#define SC3_NTIMES 12500
#define SC3_DOWNSAMPLE_TIME 25

// 40 microsecond -> 2 milisecond
#define SC4_NTIMES 25000
#define SC4_DOWNSAMPLE_TIME 50

// the same for both science case 3 and 4
#define NTIMES_LOW 500

fitsfile *output[NTABS_MAX];
unsigned int downsampled[NCHANNELS_LOW * NTIMES_LOW];
unsigned char packed[NCHANNELS_LOW * NTIMES_LOW / 8];
float offset[NCHANNELS_LOW];
float scale[NCHANNELS_LOW];

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
void fits_error_and_exit(int status, int line) {
  if (runlog) {
    fprintf(runlog, "Fits error at line %i: ", line);
    fits_report_error(runlog, status);
  }
  fprintf(stdout, "Fits error at line %i: ", line);
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

  int dc;
  for (dc = 0; dc < NCHANNELS_LOW; dc++) {

    // First pass: calculate average(=offset) and stdev(=scale)
    temp1 = &downsampled[dc * NTIMES_LOW];
    unsigned int sum = 0;
    unsigned int sos = 0;

    int dt;
    // for (dt=dc*NTIMES_LOW; dt < (dc+1)*NTIMES_LOW; dt++) {
    for (dt=0; dt < NTIMES_LOW; dt++) {
      sos += (*temp1) * (*temp1);
      sum += *temp1++;
    }
    offset[dc] = sum / (1.0 * NTIMES_LOW);
    scale[dc] = sqrt((sos / (1.0 * NTIMES_LOW)) - offset[dc] * offset[dc]);

    // Set cutoff to 1 stdev above average
    int cutoff = offset[dc] + scale[dc];

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
}

/**
 * Initialize the CFITSIO library
 * @param {char *} template_file    FITS template to use for creating initial file.
 * @param {char *} output_directory Directoy where output FITS files can be written.
 * @param {int} ntabs               Number of beams
 */
void dadafits_fits_init (char *template_file, char *output_directory, int ntabs) {
  float version;
  fits_get_version(&version);
  LOG("Using FITS library version %f\n", version);

  int t;
  for (t=0; t<ntabs; t++) {
    int status;
    char fname[256];
    fitsfile *fptr;

    if (output_directory) {
      snprintf(fname, 256, "%s/beam%c.fits(%s)", output_directory, 'A'+t, template_file);
    } else {
      snprintf(fname, 256, "beam%c.fits(%s)", 'A'+t, template_file);
    }
    LOG("Writing tab %02i to file %s\n", t, fname);

    status = 0; if (fits_create_file(&fptr, fname, &status)) fits_error_and_exit(status, __LINE__);
    status = 0; if (fits_movabs_hdu(fptr, 1, NULL, &status)) fits_error_and_exit(status, __LINE__);
    status = 0; if (fits_write_date(fptr, &status))          fits_error_and_exit(status, __LINE__);
    status = 0; if (fits_write_chksum(fptr, &status))        fits_error_and_exit(status, __LINE__);

    status = 0; if (fits_movabs_hdu(fptr, 2, NULL, &status)) fits_error_and_exit(status, __LINE__);

    output[t] = fptr;
  }
}

/**
 * Write a row of packed, 1 bit Stokes I to a FITS BINTABLE.SUBINT
 *
 * Uses the global arrays: 'scale', 'offset', and 'packed'
 *
 * @param {int} tab                Tight array beam index used to select output file
 * @param {int} rowid              Row number in the SUBINT table, corresponds to ringbuffer page number
 */
void write_fits_packed(int tab, int rowid) {
  int status;
  fitsfile *fptr = output[tab];

  // fits_write_col parameters:
  // type
  // column
  // firstrow
  // firstelem
  // nelement
  // *array
  // *status

  status = 0;
  if (fits_insert_rows(fptr, rowid, 1, &status)) {
    fits_error_and_exit(status, __LINE__);
  }

  double offs_sub = rowid + 1;
  status = 0;
  if (fits_write_col(fptr, TDOUBLE, 2, rowid + 1, 1, 1, &offs_sub, &status)) {
    fits_error_and_exit(status, __LINE__);
  }

  status = 0;
  if (fits_write_col(fptr, TFLOAT, 15, rowid + 1, 1, NCHANNELS_LOW, &offset, &status)) {
    fits_error_and_exit(status, __LINE__);
  }

  status = 0;
  if (fits_write_col(fptr, TFLOAT, 16, rowid + 1, 1, NCHANNELS_LOW, &scale, &status)) {
    fits_error_and_exit(status, __LINE__);
  }

  status = 0;
  if (fits_write_col(fptr, TBYTE,  17, rowid + 1, 1, NCHANNELS_LOW * NTIMES_LOW / 8, packed, &status)) {
    fits_error_and_exit(status, __LINE__);
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
  ascii_header_get(header, "SAMPLES_PER_BATCH", "%d", &uintValue);
  ascii_header_get(header, "CHANNELS", "%d", &uintValue);
  ascii_header_get(header, "MIN_FREQUENCY", "%f", &floatValue[0]);
  ascii_header_get(header, "CHANNEL_BANDWIDTH", "%f", &floatValue[1]);

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

int main (int argc, char *argv[]) {
  char *key;
  int padded_size;
  char *logfile;
  char *template_file;
  char *output_directory = NULL; // defaults to CWD
  int science_case;
  int science_mode;
  int ntabs;

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

  int nchannels;
  int row_size = 0;

  if (science_mode == 0) { // I + TAB to be compressed and downsampled
    row_size = NCHANNELS_LOW * NTIMES_LOW / 8;
    nchannels = NCHANNELS_LOW;
    ntabs = 12;
  } else if (science_mode == 1) { // IQUV + TAB
    exit(EXIT_FAILURE);
  } else if (science_mode == 2) { // I + IAB
    row_size = NCHANNELS_LOW * NTIMES_LOW / 8;
    nchannels = NCHANNELS_LOW;
    ntabs = 1;
  } else if (science_mode == 3) {
    exit(EXIT_FAILURE);
  } else {
    fprintf(stderr, "Illegal science mode %i\n", science_mode);
    exit(EXIT_FAILURE);
  }

  unsigned char *data = malloc(row_size);

  dadafits_fits_init(template_file, output_directory, ntabs);

  // dada_hdu_t *ringbuffer = init_ringbuffer(key);
  // ipcbuf_t *data_block = (ipcbuf_t *) ringbuffer->data_block;
  // ipcio_t *ipc = ringbuffer->data_block;

  int quit = 0;
  // uint64_t bufsz = ipc->curbufsz;

  int page_count = 0;
  // char *page = NULL;
  
    int mysize = NTABS_MAX * NCHANNELS * padded_size;
    char *page = malloc(mysize);
    int kkk;
    for(kkk=0; kkk<mysize; kkk++) {
      page[kkk] = round(10.0 * rand()/RAND_MAX );
    }
  

  // while(!quit && !ipcbuf_eod(data_block)) {
  while(!quit) {
    int tab; // tight array beam

    // Trap Ctr-C and kill commands to properly close fits files on exit
    signal(SIGTERM, fits_error_and_exit);

    // page = ipcbuf_get_next_read(data_block, &bufsz);
    if (page_count == 100) {
      page = NULL;
    }

    if (! page) {
      quit = 1;
    } else {
      switch (science_mode) {
        case 0: // stokesI data to compress and downsample
        case 2: // stokesI data to compress and downsample
          for (tab = 0; tab < ntabs; tab++) {
            if (science_case == 3) {
              downsample_sc3(tab, page, padded_size); // moves data from the page to the downsampled array
            } else if (science_case == 4) {
              downsample_sc4(tab, page, padded_size); // moves data from the page to the downsampled array
            } else {
              exit(EXIT_FAILURE);
            }
            pack_sc34(); // moves data from the downsampled array to the packed array, and sets scale and offset array
            write_fits_packed(tab, page_count); // writes data from the packed, scale, and offset array
          }
        break;

        case 1: // stokes IQUV to dump
          // TODO
          exit(EXIT_FAILURE);
          break;

        default:
          // should not happen
          fprintf(stderr, "Illegal science mode %i\n", science_mode);
          quit = 1;
          break;
      }
      // ipcbuf_mark_cleared((ipcbuf_t *) ipc);
      page_count++;
    }
  }

  // if (ipcbuf_eod(data_block)) {
  //   LOG("End of data received\n");
  // }

  close_fits();

  // dada_hdu_unlock_read(ringbuffer);
  // dada_hdu_disconnect(ringbuffer);
  LOG("Read %i pages\n", page_count);
}
