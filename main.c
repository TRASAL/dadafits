#include <string.h>
#include <stdio.h>
#include <fitsio.h>
#include <getopt.h>

#include "dada_hdu.h"
#include "ascii_header.h"

// #define LOG(...) {fprintf(stdout, __VA_ARGS__); fflush(stdout);}
#define LOG(...) {fprintf(stdout, __VA_ARGS__); fprintf(runlog, __VA_ARGS__); fflush(stdout);}

#define NTABS 12
#define NCHANNELS 1536
#define NTIMES 25000
#define DOWNSAMPLE_TIME 25

FILE *runlog = NULL;

/**
 * Downsample timeseries by summation
 *
 * @param {int} index TAB index
 * @param {ushort[NTABS,NCHANNELS, padded_size]} buffer Ring buffer page to downsample
 * @param {int} padded_size Size of fastest dimension, as timeseries are padded for optimal memory layout on GPU
 */
void downsample(int index, unsigned char *buffer, int padded_size) {
  int tab[(NCHANNELS/4) * (25000/DOWNSAMPLE_TIME)];

  int dc, dt;
  for (dc=0; dc < NCHANNELS / 4; dc++) {
    int ps0, ps1, ps2, ps3; // partial sums (per channel)
    unsigned char *s0; // channel 0
    unsigned char *s1; // channel 1
    unsigned char *s2; // channel 2
    unsigned char *s3; // channel 3

    s0 = &buffer[index * NCHANNELS * padded_size + ((dc << 2) + 0) * padded_size];
    s1 = &buffer[index * NCHANNELS * padded_size + ((dc << 2) + 1) * padded_size];
    s2 = &buffer[index * NCHANNELS * padded_size + ((dc << 2) + 2) * padded_size];
    s3 = &buffer[index * NCHANNELS * padded_size + ((dc << 2) + 3) * padded_size];

    for (dt=0; dt < NTIMES / DOWNSAMPLE_TIME; dt++) {
      ps0 = 0;
      ps1 = 0;
      ps2 = 0;
      ps3 = 0;

      int t;
      for (t=0; t < DOWNSAMPLE_TIME; t++) {
        ps0 += *s0++;
        ps1 += *s1++;
        ps2 += *s2++;
        ps3 += *s3++;
      }
      tab[dc * (NTIMES  / DOWNSAMPLE_TIME) + dt] = ps0 + ps1 + ps2 + ps3;
    }
  }
}

/**
 * Initialize the CFITSIO library
 */
void dadafits_fits_init () {
  float version;

  fits_get_version(&version);
  LOG("Using FITS library version %f\n", version);
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

  LOG("psrdada HEADER: %s\n", header);

  return hdu;
}

/**
 * Print commandline options
 */
void printOptions() {
  printf("usage: dadafits -k <hexadecimal key> -l <logfile> -b <padded_size>\n");
  printf("e.g. dadafits -k dada -l log.txt -b 25088\n");
  return;
}

/**
 * Parse commandline
 */
void parseOptions(int argc, char *argv[], char **key, int *padded_size, char **logfile) {
  int c;

  int setk=0, setb=0, setl=0;
  while((c=getopt(argc,argv,"k:b:l:"))!=-1) {
    switch(c) {
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

      default:
        printOptions();
        exit(0);
    }
  }

  // All arguments are required
  if (!setk || !setl || !setb) {
    printOptions();
    exit(EXIT_FAILURE);
  }
}


int main (int argc, char *argv[]) {
  char *key;
  int padded_size;
  char *logfile;

  // parse commandline
  parseOptions(argc, argv, &key, &padded_size, &logfile);

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

  dadafits_fits_init ();

  dada_hdu_t *ringbuffer = init_ringbuffer(key);
  ipcbuf_t *data_block = (ipcbuf_t *) ringbuffer->data_block;
  ipcio_t *ipc = ringbuffer->data_block;

  int quit = 0;
  uint64_t bufsz = ipc->curbufsz;

  int page_count = 0;
  char *page = NULL;
  while(!quit && !ipcbuf_eod(data_block)) {
    int tab; // tight array beam

    page = ipcbuf_get_next_read(data_block, &bufsz);
    if (! page) {
      quit = 1;
    } else {
      for (tab = 0; tab < NTABS; tab++) {
        downsample(tab, page, padded_size);
      }
      ipcbuf_mark_cleared((ipcbuf_t *) ipc);
      page_count++;
    }
  }

  if (ipcbuf_eod(data_block)) {
    LOG("End of data received\n");
  }

  dada_hdu_unlock_read(ringbuffer);
  dada_hdu_disconnect(ringbuffer);
  LOG("Read %i pages\n", page_count);
}
