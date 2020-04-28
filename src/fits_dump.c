/**
 * program: fits_dump
 *          Written for the AA-Alert project, ASTRON
 *
 * Purpose: help debugging / reading 1-bit compressed FITS files
 *
 * Author: Jisk Attema, Netherlands eScience Center
 * Licencse: Apache v2.0
 */

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include "dadafits_internal.h"

void print_table(fitsfile *fptr) {
  int status, tstatus;

  long nrows;
  fits_get_num_rows(fptr, &nrows, &status);

  int ncols, c;
  fits_get_num_cols(fptr, &ncols, &status);

  printf("Table contains %li rows and %i columns\n", nrows, ncols);

  char colname[256];
  int typecode;
  long repeat, width;

  fits_get_colname(fptr, CASEINSEN, "*", colname, &c, &status);
  while(status != COL_NOT_FOUND) {
    fits_get_coltype(fptr, c, &typecode, &repeat, &width, &tstatus);
    switch (typecode) {
      case TSTRING:     printf("%03i %10s TSTRING %li %li\n", c, colname, repeat, width); break;
      case TSHORT:      printf("%03i %10s TSHORT %li %li\n", c, colname, repeat, width); break;
      case TLONG:       printf("%03i %10s TLONG %li %li\n", c, colname, repeat, width); break;
      case TFLOAT:      printf("%03i %10s TFLOAT %li %li\n", c, colname, repeat, width); break;
      case TDOUBLE:     printf("%03i %10s TDOUBLE %li %li\n", c, colname, repeat, width); break;
      case TLOGICAL:    printf("%03i %10s TLOGICAL %li %li\n", c, colname, repeat, width); break;
      case TBIT:        printf("%03i %10s TBIT %li %li\n", c, colname, repeat, width); break;
      case TBYTE:       printf("%03i %10s TBYTE %li %li\n", c, colname, repeat, width); break;
      case TCOMPLEX:    printf("%03i %10s TCOMPLEX %li %li\n", c, colname, repeat, width); break;
      case TDBLCOMPLEX: printf("%03i %10s TDBLCOMPLEX %li %li\n", c, colname, repeat, width); break;
      // case TINT32BIT:   printf("%03i TINT32BIT %li %li\n", c, colname, repeat, width); break;
      default:          printf("%03i %10s Unknown %li %li\n", c, colname, repeat, width); break;
    }
    fits_get_colname(fptr, CASEINSEN, "*", colname, &c, &status);
  }

}

void fitsinfo(fitsfile *fptr) {
  int status = 0;
  fits_movabs_hdu(fptr, 2, NULL, &status); if (status) fits_report_error(stdout, status);

  int cf, cw, cs, co, cos, cd;
  status = 0; fits_get_colnum(fptr, 0, "OFFS_SUB", &cos, &status); if (status) fits_report_error(stdout, status);
  status = 0; fits_get_colnum(fptr, 0, "DAT_FREQ", &cf, &status); if (status) fits_report_error(stdout, status);
  status = 0; fits_get_colnum(fptr, 0, "DAT_WTS", &cw, &status); if (status) fits_report_error(stdout, status);
  status = 0; fits_get_colnum(fptr, 0, "DAT_OFFS", &co, &status); if (status) fits_report_error(stdout, status);
  status = 0; fits_get_colnum(fptr, 0, "DAT_SCL", &cs, &status); if (status) fits_report_error(stdout, status);
  status = 0; fits_get_colnum(fptr, 0, "DATA", &cd, &status); if (status) fits_report_error(stdout, status);

  double offs_sub;
  float freqs[NCHANNELS_LOW], wts[NCHANNELS_LOW], scale[NCHANNELS_LOW], offs[NCHANNELS_LOW];
  unsigned char data[NCHANNELS_LOW * NTIMES_LOW / 8]; // array of 1 bit datapoints of [channel=NCHANNELS_LOW, time=NTIMES_LOW]

  long nrows;
  fits_get_num_rows(fptr, &nrows, &status);

  for (long row = 1; row <= nrows; row++) {
    int channel, time, channeltime;

    status = 0; fits_read_col(fptr, TDOUBLE, cos, row, 1L, 1L, 0, &offs_sub, NULL, &status);
    status = 0; fits_read_col(fptr, TFLOAT, cf, row, 1L, NCHANNELS_LOW, 0, freqs, NULL, &status);
    status = 0; fits_read_col(fptr, TFLOAT, cw, row, 1L, NCHANNELS_LOW, 0, wts,   NULL, &status);
    status = 0; fits_read_col(fptr, TFLOAT, cs, row, 1L, NCHANNELS_LOW, 0, scale, NULL, &status);
    status = 0; fits_read_col(fptr, TFLOAT, co, row, 1L, NCHANNELS_LOW, 0, offs,  NULL, &status);
    status = 0; fits_read_col(fptr, TBYTE, cd, 1L, 1L, NCHANNELS_LOW * NTIMES_LOW / 8, 0, data, NULL, &status);
    if (status) fits_report_error(stdout, status);

    printf("\n\n# Frequencies, weights, scales and offsets for subint: %li, subint offset %e\n", row, offs_sub);

    for (channel=0; channel<NCHANNELS_LOW; channel++){
      printf("%i %lf %lf %lf %lf\n", channel, freqs[channel], wts[channel], scale[channel], offs[channel]);
    }

    printf("\n\n# Data, 1 bit, shown per byte [channel, time]:\n");

    for (channeltime=0; channeltime<NCHANNELS_LOW*NTIMES_LOW/8; channeltime++) {
      channel = (channeltime * 8) / NTIMES_LOW;
      time = (channeltime * 8) % NTIMES_LOW;
      if (channeltime % 20 == 0) {
        printf( "\n[channel= % 6i, time=  % 6i] ", channel, time);
      }
      printf( " % 4i", data[channeltime]);
    }
    printf("\n");
  }
}

void print_hdus(fitsfile *fptr) {
  int hdunum;
  int status;

  status = 0;
  fits_get_num_hdus(fptr, &hdunum, &status);
  fits_report_error(stdout, status);

  printf("File contains %i HDUs\n", hdunum);

  int h;
  for (h=1; h <= hdunum; h++) {
    int type;
    fits_movabs_hdu(fptr, h, &type, &status);
    switch(type) {
      case IMAGE_HDU:
        printf("%i is an IMAGE_HDU\n", h);
        break;
      case ASCII_TBL:
        printf("%i is an ASCII_TBL\n", h);
        print_table(fptr);
        break;
      case BINARY_TBL:
        printf("%i is an BINARY_TBL\n", h);
        print_table(fptr);
        break;
      default:
        printf("%i is an unknown type\n", h);
        break;
    }
  }
}

int main(int argc, char *argv[]) {
  int status;
  fitsfile *fptr;

  printf("Opening: '%s'\n", argv[1]);
  status = 0;
  fits_open_file(&fptr, argv[1], READONLY, &status);

  print_hdus(fptr);
  fitsinfo(fptr);

  fits_close_file (fptr, &status);
}
