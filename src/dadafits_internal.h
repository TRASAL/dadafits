#ifndef __HAVE_DADAFITS_INTERNAL_H__
#define __HAVE_DADAFITS_INTERNAL_H__

#include <stdio.h>

extern FILE *runlog;
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

// The synthesized beams table
#define NSYNS_MAX 256
#define NSUBBANDS 32
#define SUBBAND_UNSET 9999
#define FREQS_PER_SUBBAND 48

// Global parameter definintions
extern int science_case;
extern int science_mode;
extern int padded_size;

extern float fits_offset[NCHANNELS * NPOLS];
extern float fits_scale[NCHANNELS * NPOLS];
extern float fits_weights[NCHANNELS];
extern float fits_freqs[NCHANNELS];

extern int synthesized_beam_table[NSYNS_MAX][NSUBBANDS];
extern int synthesized_beam_selected[NSYNS_MAX];
extern int synthesized_beam_count; // number of SBs in the table

// Function definitions

// from downsample.c
extern void downsample_sc3(const unsigned char *buffer, const int padded_size, unsigned int downsampled[NCHANNELS_LOW * NTIMES_LOW]);
extern void downsample_sc4(const unsigned char *buffer, const int padded_size, unsigned int downsampled[NCHANNELS_LOW * NTIMES_LOW]);

// from sb_util.c
extern int read_synthesized_beam_table(char *fname);
extern void parse_synthesized_beam_selection (char *selection);

// from fits_io.c
extern void dadafits_fits_init (const char *template_dir, const char *template_file, const char *output_directory,
    const int ntabs, const int mode, const float min_frequency, const float channelwidth, char *ra_hms, char *dec_hms,
    char *source_name, const char *utc_start, const double mjd_start, double lst_start);
extern void write_fits(const int tab, const int channels, const int pols, const long rowid, const int rowlength, unsigned char *data, const float telaz, const float telza);
extern void close_fits();
extern void fits_error_and_exit(int status); // needed for trapping C-c

// from manipulate.c
extern void deinterleave (const unsigned char *page, const int ntabs, const int sequence_length, unsigned char *transposed);
extern void pack_sc34(unsigned int downsampled[NCHANNELS_LOW * NTIMES_LOW], unsigned char packed[NCHANNELS_LOW * NTIMES_LOW/8]);

#endif
