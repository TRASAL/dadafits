/**
 * Functions to:
 *  parse a Synthesized Beam table
 *  parse a selection string
 */
#include <stdlib.h>
#include <string.h>

#include "dadafits_internal.h"

int synthesized_beam_table[NSYNS_MAX][NSUBBANDS];
int synthesized_beam_selected[NSYNS_MAX];
int synthesized_beam_count; // number of SBs in the table

/**
 * Read the synthesized beam table
 */
int read_synthesized_beam_table(char *fname) {
  int subband_index, syn_index;

  LOG("Reading synthesized beam table '%s'\n", fname);

  // Permitted whitespace
  char delim[4];
  delim[0] = ' '; // space
  delim[1] = '\t'; // tab
  delim[2] = '\n'; // newline If not included, empty lines in the txt file result in  a single '\n' token
  delim[3] = '\0';

  FILE *table = fopen(fname, "r");
  if (! table) {
    LOG("Error: Cannot read file: '%s'\n", fname);
    return 1;
  }

  syn_index = 0;

#define LINELENGTH 512
  char line[LINELENGTH];
  while (fgets(line,LINELENGTH,table)) {
    // remove comments
    char *hash = index(line, '#');
    if (hash) {
      *hash = '\0';
    }

    subband_index = 0;

    char *saveptr;
    char *key = strtok_r(line, delim, &saveptr);
    while (key) {
      if (subband_index == NSUBBANDS) {
        LOG("Error: Too many subbands (more than %i), increase NSUBBANDS\n", NSUBBANDS);
        exit(EXIT_FAILURE);
      }

      synthesized_beam_table[syn_index][subband_index] = atoi(key);

      // Convert indexing from [-NTABS_MAX/2, NTABS_MAX/2  - 1] := [-6,5] to [0,NTABS_MAX] := [0,11]
      synthesized_beam_table[syn_index][subband_index] += NTABS_MAX/2;
      if (synthesized_beam_table[syn_index][subband_index] < 0 || synthesized_beam_table[syn_index][subband_index] >= NTABS_MAX) {
        LOG("Error: illegal TAB entry '%s' at %i for synthesized beam %i\n", key, subband_index, syn_index);
        exit(EXIT_FAILURE);
      }


      key = strtok_r(NULL, delim, &saveptr); // next token
      subband_index++;
    }

    if (subband_index != 0) { // did we read anything from this line?
      if (subband_index != NSUBBANDS) {
        LOG("Error: wrong number of subbands (%i) for beam %i\n", subband_index, syn_index);
        exit(EXIT_FAILURE);
      }

      // go to the next row
      syn_index++;
      if (syn_index == NSYNS_MAX) {
        LOG("Too many synthesized beams (more than %i), increase NSYNS_MAX\n", NSYNS_MAX);
        exit(EXIT_FAILURE);
      }
    }
  }

  synthesized_beam_count = syn_index;
  LOG("Read %i synthesized beams\n", synthesized_beam_count);

  // clear the remaining rows of the table
  while (syn_index < NSYNS_MAX) {
    subband_index = 0;
    while (subband_index < NSUBBANDS) {
      synthesized_beam_table[syn_index][subband_index++] = SUBBAND_UNSET;
    }
    syn_index++;
  }
}

/**
 * Parse commandline for synthesized beams
 * 
 * Read a table via the -S option
 * Setup up a selection of synthesized beams to process using -s
 *     
 *  Selections are a comma separated list of
 *   a) single beam numbers
 *   b) beam ranges as start/end where both limits are inclusive
 *                     
 *  NOTE:
 *   * Repeating a beam has no effect
 *   * Beams start counting at -
 *                                      
 *  So something like this: -s 0,3/4,6/7,3
 *  would result in beams 0, 3, 4, 6, and 7 to be processed
 *                                                   
 * Program terminates on unparseable/illegal beam selections
 */
void parse_synthesized_beam_selection (char *selection) {
  int sb;

  for (sb=0; sb<synthesized_beam_count; sb++) {
    if (! selection) {
      // without selection, process all
      synthesized_beam_selected[sb] = 1;
    } else {
      // by default, do not process
      synthesized_beam_selected[sb] = 0;
    }
  }

  if (! selection) {
    LOG("Processing all synthesized beams\n");
    return;
  }

  LOG("Sythesized beam list: ");

  // split string on commas
  char delim[2];
  delim[0]=',';
  delim[1]='\0';
  char *saveptr;

  char *key = strtok_r(selection, delim, &saveptr);
  while (key) {
    char *key_start = key;
    char *key_end = index(key, '/');
  
    int s = atoi(key_start);
    int e = 0;
    if (key_end) {
        e = atoi(&key_end[1]); // skip leading dash, worst case: key='3/' the string key_end is now empty ('\0')
    }

    if (key_start && key_end) {
      if (s>=synthesized_beam_count || e>=synthesized_beam_count || e < s || s<0 || e<0) {
        LOG("Error: Invalid range: '%s'\n", key);
        exit(EXIT_FAILURE);
      }
      for (sb = s; sb <= e; sb++) {
        synthesized_beam_selected[sb] = 1;
        LOG(" %i", sb);
      }
    } else {
      if (s<0 || s >=synthesized_beam_count) {
        LOG("Error: Invalid beam: '%s'\n", key);
        exit(EXIT_FAILURE);
      }
      synthesized_beam_selected[s] = 1;
      LOG(" %i", s);
    }
    key = strtok_r(NULL, delim, &saveptr); // next token
  }
  LOG("\n");
}
