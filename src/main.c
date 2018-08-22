/**
 * program: dadafits
 *          Written for the AA-Alert project, ASTRON
 *
 * Purpose: connect to a ring buffer and create FITS output per TAB
 *          Depending on science case and mode, reduce time and frequency resolution to 1 bit
 *          Fits files are created using templates
 *
 * Science case 3, mode 0
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
 *          template: sc3_IQUV.txt
 *
 *          A ringbuffer page is interpreted as an interleaved array of Stokes IQUV:
 *          [tab][channel_offset][sequence_number][packet]
 *
 *          tab             := ranges from 0 to NTABS (12)
 *          channel_offset  := ranges from 0 to NCHANNELS/4 (384)
 *          sequence_number := ranges from 0 to 25
 *
 *          with a packet: [t0 .. t499][c0 .. c3][IQUV] total of 500*4*4=8000 bytes
 *          t = tn + sequence_number * 500
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
 * Science case 3, mode 3:
 *          template: sc3_IQUV.txt
 *
 *          A ringbuffer page is interpreted as an interleaved array of Stokes IQUV:
 *          [tab][channel_offset][sequence_number][packet]
 *
 *          tab             := ranges from 0 to NTABS (1)
 *          channel_offset  := ranges from 0 to NCHANNELS/4 (384)
 *          sequence_number := ranges from 0 to 25
 *
 *          with a packet: [t0 .. t499][c0 .. c3][IQUV] total of 500*4*4=8000 bytes
 *          t = tn + sequence_number * 500
 *          c = cn + channel_offset * 4
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
 *          template: sc4_IQUV.txt
 *
 *          A ringbuffer page is interpreted as an interleaved array of Stokes IQUV:
 *          [tab][channel_offset][sequence_number][packet]
 *
 *          tab             := ranges from 0 to NTABS (12)
 *          channel_offset  := ranges from 0 to NCHANNELS/4 (384)
 *          sequence_number := ranges from 0 to 50
 *
 *          with a packet: [t0 .. t499][c0 .. c3][IQUV] total of 500*4*4=8000 bytes
 *          t = tn + sequence_number * 500
 *          c = cn + channel_offset * 4
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
 * Science case 4, mode 3:
 *          template: sc4_IQUV.txt
 *
 *          A ringbuffer page is interpreted as an interleaved array of Stokes IQUV:
 *          [tab][channel_offset][sequence_number][packet]
 *
 *          tab             := ranges from 0 to NTABS (1)
 *          channel_offset  := ranges from 0 to NCHANNELS/4 (384)
 *          sequence_number := ranges from 0 to 50
 *
 *          with a packet: [t0 .. t499][c0 .. c3][IQUV] total of 500*4*4=8000 bytes
 *          t = tn + sequence_number * 500
 *          c = cn + channel_offset * 4
 *
 * Author: Jisk Attema, Netherlands eScience Center
 * Licencse: Apache v2.0
 */

#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <signal.h>

#include "dada_hdu.h"
#include "ascii_header.h"
#include "config.h"

#include "dadafits_internal.h"
int padded_size;
int science_case;
int science_mode;

FILE *runlog = NULL;

const char *science_modes[] = {"I+TAB", "IQUV+TAB", "I+IAB", "IQUV+IAB"};
const char *template_case3mode13 = "sc3_IQUV.txt";
const char *template_case34mode02 = "sc34_1bit_I_reduced.txt";
const char *template_case4mode13 = "sc4_IQUV.txt";

// Variables read from ring buffer header
float min_frequency;
float bandwidth = 300;
char ra_hms[256];
char dec_hms[256];
float scanlen;
float center_frequency;
char parset[24567];
char source_name[256];
char utc_start[256];
double mjd_start;
double lst_start;
float az_start;
float za_start;

// Variables set from commandline
int make_synthesized_beams = 0;

// Memory buffers
unsigned int downsampled[NCHANNELS_LOW * NTIMES_LOW];
unsigned char packed[NCHANNELS_LOW * NTIMES_LOW / 8];
unsigned char *transposed = NULL; // Stokes IQUV buffer of approx 2 GB, allocated only when necessary
unsigned char *synthesized = NULL; // Stokes IQUV for a single synthesized beam

/**
 * Open a connection to the ringbuffer
 *
 * @param {char *} key String containing the shared memory key as hexadecimal number
 * @returns {hdu *} A connected HDU
 */
dada_hdu_t *init_ringbuffer(char *key) {
  uint64_t nbufs;
  int header_incomplete = 0;

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
  if (ascii_header_get(header, "MIN_FREQUENCY", "%f", &min_frequency) == -1) {
    LOG("ERROR. MIN_FREQUENCY not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "BW", "%f", &bandwidth) == -1) {
    LOG("ERROR. BW not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "PADDED_SIZE", "%i", &padded_size) == -1) {
    LOG("ERROR. PADDED_SIZE not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "SCIENCE_CASE", "%i", &science_case) == -1) {
    LOG("ERROR. SCIENCE_CASE not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "SCIENCE_MODE", "%i", &science_mode) == -1) {
    LOG("ERROR. SCIENCE_MODE not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "RA_HMS", "%s", ra_hms) == -1) {
    LOG("ERROR. RA not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "DEC_HMS", "%s", dec_hms) == -1) {
    LOG("ERROR. DEC not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "SCANLEN", "%f", &scanlen) == -1) {
    LOG("ERROR. SCANLEN not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "FREQ", "%f", &center_frequency) == -1) {
    LOG("ERROR. FREQ not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "SOURCE", "%s", source_name) == -1) {
    LOG("ERROR. SOURCE not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "UTC_START", "%s", utc_start) == -1) {
    LOG("ERROR. UTC_START not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "MJD_START", "%lf", &mjd_start) == -1) {
    LOG("ERROR. MJD_START not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "LST_START", "%lf", &lst_start) == -1) {
    LOG("ERROR. LST_START not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "AZ_START", "%f", &az_start) == -1) {
    LOG("ERROR. AZ_START not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "ZA_START", "%f", &za_start) == -1) {
    LOG("ERROR. ZA_START not set in dada buffer\n");
    header_incomplete = 1;
  }
  if (ascii_header_get(header, "PARSET", "%s", parset) == -1) {
    LOG("ERROR. PARSET not set in dada buffer\n");
    header_incomplete = 1;
  }

  // tell the ringbuffer the header has been read
  if (ipcbuf_mark_cleared(hdu->header_block) < 0) {
    LOG("ERROR. Cannot mark the header as cleared\n");
    exit(EXIT_FAILURE);
  }

  LOG("psrdada HEADER:\n%s\n", header);

  if (header_incomplete) {
    exit(EXIT_FAILURE);
  }

  return hdu;
}

/**
 * Print commandline options
 */
void printOptions() {
  printf("usage: dadafits -k <hexadecimal key> -l <logfile> -t <template> -d <output_directory> -S <synthesized beam table> -s <synthesize these beams>\n");
  printf("e.g. dadafits -k dada -l log.txt -c 3 -m 0 -b 25088 -t /full/path/template.txt -S table.txt -s 0,1,4-8 -d /output/directory\n");
  return;
}

/**
 * Parse commandline
 */
void parseOptions(int argc, char *argv[], char **key, char **logfile, char **template_dir, char **table_name, char **sb_selection, char **output_directory) {
  int c;

  int setk=0, setl=0;
  while((c=getopt(argc,argv,"k:l:t:d:s:S:"))!=-1) {
    switch(c) {
      // OPTIONAL: -d <output_directory>
      // DEFAULT: CWD
      case('d'):
        *output_directory = strdup(optarg);
        break;

      // OPTIONAL: -t <template_dir>
      // DEFAULT: CWD/templates
      case('t'):
        *template_dir = strdup(optarg);
        break;

      // -k <hexadecimal_key>
      case('k'):
        *key = strdup(optarg);
        setk=1;
        break;

      // -l log file
      case('l'):
        *logfile = strdup(optarg);
        setl=1;
        break;

      // OPTIONAL: -S synthesized beam table
      case('S'):
        *table_name = strdup(optarg);
        break;

      // OPTIONAL: -s synthesized beam selection
      case('s'):
        *sb_selection = strdup(optarg);
        break;

      default:
        printOptions();
        fprintf(stderr, "Unknown option: %c\n",  c);
        exit(EXIT_FAILURE);
    }
  }

  // Required arguments
  if (!setk || !setl) {
    printOptions();
    exit(EXIT_FAILURE);
  }
}

int main (int argc, char *argv[]) {
  char *key;
  char *logfile;
  const char *template_file = NULL;
  char *template_dir = "templates";
  char *table_name = NULL; // optional argument
  char *sb_selection = NULL; // optional argument, defaults to all beams
  char *output_directory = NULL; // defaults to CWD
  int ntabs;
  int nchannels; // for FITS outputfile (so after optional compression)
  int ntimes; // for FITS outputfile (so after optional compression)
  int npols; // for FITS outputfile (so after optional compression)
  int sequence_length;

  // parse commandline
  parseOptions(argc, argv, &key, &logfile, &template_dir, &table_name, &sb_selection, &output_directory);

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

  // must init ringbuffer before fits, as this reads parameters
  // like bandwidth from ring buffer header
  dada_hdu_t *ringbuffer = init_ringbuffer(key);
  ipcbuf_t *data_block = (ipcbuf_t *) ringbuffer->data_block;
  ipcio_t *ipc = ringbuffer->data_block;
  uint64_t bufsz = ipc->curbufsz;

  LOG("dadafits version: " VERSION "\n");

  if (table_name) {
    LOG("Writing synthesized beams\n");
    make_synthesized_beams = 1;
    read_synthesized_beam_table(table_name);
    parse_synthesized_beam_selection(sb_selection);
  } else {
    LOG("Writing TABs (not synthesized beams)\n");
    make_synthesized_beams = 0;
  }

  switch (science_mode) {
    case 0: // I + TAB to be compressed and downsampled
      ntimes = NTIMES_LOW;
      nchannels = NCHANNELS_LOW;
      ntabs = 12;
      npols = 1;

      // adjust min_frequency for downsampling:
      // before |  x  |     |    |    |
      // after  |  x        X         | small 'x' should be large 'X' : add 1.5 of the original channels
      min_frequency = min_frequency + (1.5 * bandwidth / ((float) NCHANNELS));

      if (make_synthesized_beams) {
        LOG("Cannot write synthesized beams for compressed I+TAB\n");
        exit(EXIT_FAILURE);
      }
      template_file = template_case34mode02;
      break;
    case 1: // IQUV + TAB to deinterleave
      ntimes = science_case == 3 ? 12500 : 25000;
      nchannels = NCHANNELS;
      ntabs = 12;
      npols = 4;
      break;
    case 2: // I + IAB to be compressed and downsampled
      ntimes = NTIMES_LOW;
      nchannels = NCHANNELS_LOW;
      ntabs = 1;
      npols = 1;

      // adjust min_frequency for downsampling:
      // before |  x  |     |    |    |
      // after  |  x        X         | small 'x' should be large 'X' : add 1.5 of the original channels
      min_frequency = min_frequency + (1.5 * bandwidth / ((float) NCHANNELS));

      if (make_synthesized_beams) {
        LOG("Cannot write synthesized beams for compressed I+IAB\n");
        exit(EXIT_FAILURE);
      }
      template_file = template_case34mode02;
      break;
    case 3: // IQUV + IAB to deinterleave
      ntimes = science_case == 3 ? 12500 : 25000;
      nchannels = NCHANNELS;
      ntabs = 1;
      npols = 4;
      break;
    default:
      LOG("Illegal science mode %i\n", science_mode);
      exit(EXIT_FAILURE);
  }

  switch (science_case) {
    case 3:
      sequence_length = 25;
      if (padded_size < 12500) {
        LOG("Error: padded_size too small, should be at least 12500 for science case 3\n");
        exit(EXIT_FAILURE);
      }
      if (! template_file) {
        template_file = template_case3mode13;
      }
      break;
    case 4:
      sequence_length = 50;
      if (padded_size < 25000) {
        LOG("Error: padded_size too small, should be at least 25000 for science case 4\n");
        exit(EXIT_FAILURE);
      }
      if (! template_file) {
        template_file = template_case4mode13;
      }
      break;
    default:
      LOG("Illegal science case %i\n", science_case);
      exit(EXIT_FAILURE);
  }

  LOG("Science mode: %i [ %s ]\n", science_mode, science_modes[science_mode]);
  LOG("Science case: %i\n", science_case);
  LOG("Template: %s\n", template_file);

  LOG("Output to FITS tabs: %i, channels: %i, polarizations: %i, samples: %i\n", ntabs, nchannels, npols, ntimes);
  dadafits_fits_init(template_dir, template_file, output_directory,
      ntabs, make_synthesized_beams, scanlen, center_frequency, bandwidth, min_frequency, bandwidth / nchannels,
      ra_hms, dec_hms, source_name, utc_start, mjd_start, lst_start, parset);

  if (science_mode == 1 || science_mode == 3) {
    LOG("Allocating Stokes IQUV transpose buffer (%i,%i,%i,%i)\n", ntabs, NCHANNELS, NPOLS, ntimes);
    transposed = malloc(ntabs * NCHANNELS * NPOLS * ntimes * sizeof(char));
    if (transposed == NULL) {
      LOG("Could not allocate stokes IQUV transpose matrix\n");
      exit(EXIT_FAILURE);
    }
  }

  if (make_synthesized_beams) {
    LOG("Allocating Stokes IQUV synthesized beam buffer (1,%i,%i,%i)\n", NCHANNELS, NPOLS, ntimes);
    synthesized = malloc(1 * NCHANNELS * NPOLS * ntimes * sizeof(char));
    if (synthesized == NULL) {
      LOG("Could not allocate stokes IQUV synthesized beam buffer\n");
      exit(EXIT_FAILURE);
    }
  }

  int quit = 0;
  long page_count = 0;
  char *page = NULL;

  // Trap Ctr-C to properly close fits files on exit
  signal(SIGTERM, fits_error_and_exit);

  while(!quit && !ipcbuf_eod(data_block)) {
    page = ipcbuf_get_next_read(data_block, &bufsz);

    int tab; // Tied array beam
    int sb;  // synthesized beam

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
              downsample_sc3(&page[tab * NCHANNELS * padded_size], padded_size, downsampled);
            } else if (science_case == 4) {
              downsample_sc4(&page[tab * NCHANNELS * padded_size], padded_size, downsampled);
            } else {
              exit(EXIT_FAILURE);
            }
            // pack data from the downsampled array to the packed array,
            // and set scale and offset arrays with used values
            pack_sc34(downsampled, packed);

            // NOTE: Use hardcoded values instead of the variables ntimes, nchannels, npols
            // because at this point in the program they can only have these values, and this could
            // possibly allow some more optimizations
            write_fits(
              tab,
              NCHANNELS_LOW,
              1, // only Stokes I
              page_count + 1, // page_count starts at 0, but FITS rowid at 1
              NCHANNELS_LOW * NTIMES_LOW / 8,
              packed, // write data from the packed array to file, also uses scale, weights, and offset arrays
              az_start, za_start
            );
          }
          break;

        // stokesIQUV data to (optionally synthesize) and write
        case 1:
        case 3:
          // transpose data from page to transposed buffer
          deinterleave(page, ntabs, sequence_length, transposed);

          if (make_synthesized_beams) {
            // synthesize beams
            //
            // Input: transposed buffer   [TABS, CHANNELS, POLS, TIMES]
            // Output: synthesized buffer [CHANNELS, POLS, TIMES]
            
            for (sb = 0; sb < synthesized_beam_count; sb++) {
              if (synthesized_beam_selected[sb]) {
                int band;

                // a subband contains 1536/32=48 frequencies from a TAB
                for (band = 0; band < NSUBBANDS; band++) {
                  // find the TAB for this subband, and check validity
                  tab = synthesized_beam_table[sb][band]; // TODO: subband from table -> tab?
                  if (tab == SUBBAND_UNSET || tab > ntabs) {
                    LOG("Error: illegal subband index %i in synthesized beam %i\n", tab, sb);
                    exit(EXIT_FAILURE);
                  }

                  memcpy(
                    &synthesized[band * FREQS_PER_SUBBAND * NPOLS * ntimes],
                    &transposed[(tab * NCHANNELS + band * FREQS_PER_SUBBAND) * NPOLS * ntimes],
                    FREQS_PER_SUBBAND * NPOLS * ntimes
                  );
                }

                // write data from synthesized buffer
                write_fits(
                  sb,
                  NCHANNELS,
                  NPOLS, // full Stokes IQUV
                  page_count + 1, // page_count starts at 0, but FITS rowid at 1
                  NCHANNELS * NPOLS * ntimes,
                  synthesized,
                  az_start, za_start
                );
              }
            }
          } else {
            LOG("TABs %i\n", ntabs);
            // do not synthesize, but use TABs
            for (tab = 0; tab < ntabs; tab++) {
              // write data from transposed buffer, also uses scale, weights, and offset arrays (but set to neutral values)
              write_fits(
                tab,
                NCHANNELS,
                NPOLS, // full Stokes IQUV
                page_count + 1, // page_count starts at 0, but FITS rowid at 1
                NCHANNELS * NPOLS * ntimes,
                &transposed[tab * NCHANNELS * NPOLS * ntimes],
                az_start, za_start
              );
            }
          }
          break;

        default:
          // should not happen
          LOG("Illegal science mode %i\n", science_mode);
          quit = 1;
          break;
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

  LOG("Read %li pages\n", page_count);

  close_fits();
}
