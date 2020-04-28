#include "dadafits_internal.h"

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
 * @param {uchar[NCHANNELS, padded_size]} buffer        Buffer page to downsample
 * @param {int} padded_size                             Size of fastest dimension, as timeseries are padded for optimal memory layout on GPU
 * @param {uint[NCHANNELS_LOW, NTIMES_LOW]} downsampled Output array holding downsampled data
 */
void downsample_sc3(const unsigned char *buffer, const int padded_size, unsigned int downsampled[NCHANNELS_LOW * NTIMES_LOW]) {
  unsigned int *temp1 = downsampled;
  int dc; // downsampled channel
  int dt; // downsampled time
  int t; // full time

  for (dc=0; dc < NCHANNELS_LOW; dc++) {
    // pointer to next sample in the two channels
    unsigned const char *s0 = &buffer[((dc << 1) + 0) * padded_size];
    unsigned const char *s1 = &buffer[((dc << 1) + 1) * padded_size];

    for (dt=0; dt < NTIMES_LOW; dt++) {
      // partial sums (per channel)
      unsigned int ps0 = 0;
      unsigned int ps1 = 0;

      for (t=0; t < SC3_DOWNSAMPLE_TIME; t++) {
        ps0 += *s0++;
        ps1 += *s1++;
      }
      *temp1++ = ps0 + ps1;
    }
  }
}

void downsample_sc4(const unsigned char *buffer, const int padded_size, unsigned int downsampled[NCHANNELS_LOW * NTIMES_LOW]) {
  unsigned int *temp1 = downsampled;
  int dc; // downsampled channel
  int dt; // downsampled time
  int t; // full time

  for (dc=0; dc < NCHANNELS_LOW; dc++) {
    // pointer to next sample in the two channels
    unsigned const char *s0 = &buffer[((dc << 1) + 0) * padded_size];
    unsigned const char *s1 = &buffer[((dc << 1) + 1) * padded_size];

    for (dt=0; dt < NTIMES_LOW; dt++) {
      // partial sums (per channel)
      unsigned int ps0 = 0;
      unsigned int ps1 = 0;

      for (t=0; t < SC4_DOWNSAMPLE_TIME; t++) {
        ps0 += *s0++;
        ps1 += *s1++;
      }
      *temp1++ = ps0 + ps1;
    }
  }
}
