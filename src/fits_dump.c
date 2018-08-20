#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0')
// printf( BYTE_TO_BINARY_PATTERN "\n", BYTE_TO_BINARY(data[i]));

#include <string.h>
#include <stdio.h>
#include "fitsio.h"

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
  float freqs[384], wts[384], scale[384], offs[384];
  unsigned char data[384 * 500 / 8]; // array of 1 bit datapoints of [channel=384, time=500]

  long nrows;
  fits_get_num_rows(fptr, &nrows, &status);

  for (long row = 1; row <= nrows; row++) {
    int channel, time, channeltime;

    status = 0; fits_read_col(fptr, TDOUBLE, cos, row, 1L, 1L, 0, &offs_sub, NULL, &status);
    status = 0; fits_read_col(fptr, TFLOAT, cf, row, 1L, 384L, 0, freqs, NULL, &status);
    status = 0; fits_read_col(fptr, TFLOAT, cw, row, 1L, 384L, 0, wts,   NULL, &status);
    status = 0; fits_read_col(fptr, TFLOAT, cs, row, 1L, 384L, 0, scale, NULL, &status);
    status = 0; fits_read_col(fptr, TFLOAT, co, row, 1L, 384L, 0, offs,  NULL, &status);
    status = 0; fits_read_col(fptr, TBYTE, cd, 1L, 1L, 384 * 500 / 8, 0, data, NULL, &status);
    if (status) fits_report_error(stdout, status);

    printf("\n\n# Frequencies, weights, scales and offsets for subint: %li, subint offset %e\n", row, offs_sub);

    for (channel=0; channel<384; channel++){
      printf("%i %lf %lf %lf %lf\n", channel, freqs[channel], wts[channel], scale[channel], offs[channel]);
    }

    printf("\n\n# Data, 1 bit, shown per byte [channel, time]:\n");

    for (channeltime=0; channeltime<384*500/8; channeltime++) {
      channel = (channeltime * 8) / 500;
      time = (channeltime * 8) % 500;
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

  fitsinfo(fptr);

  fits_close_file (fptr, &status);
}

int write_example(int argc, char *argv[]) {
  int status;
  fitsfile *fptr;

  fits_create_file(&fptr, "example.fits.gz(template.txt)", &status);
 
  fits_movabs_hdu(fptr, 1, NULL, &status);
  fits_write_date(fptr, &status);
  fits_write_chksum(fptr, &status);
 
  fits_movabs_hdu(fptr, 2, NULL, &status);
  fits_insert_rows(fptr, 0, 1, &status);
 
  char bits[48000];
  fits_write_col(fptr,
      TBYTE,        // type
      17,           // colnum
      1,            // firstrow
      1,            // firstelem
      48000,        // nelement
      &bits,        // *array
      &status       // *status
  );
 
  fits_write_chksum(fptr, &status);
  fits_close_file (fptr, &status);
}

/*
  fits_write_col(fptr, TDOUBLE, 2, rowid, 1, 1, offs_sub, &status);
  fits_write_col(fptr, TDOUBLE, 3, rowid, 1, 1, lst_sub,  &status);
  fits_write_col(fptr, TDOUBLE, 4, rowid, 1, 1, ra_sub,   &status);
  fits_write_col(fptr, TDOUBLE, 5, rowid, 1, 1, dec_sub,  &status);
  fits_write_col(fptr, TDOUBLE, 6, rowid, 1, 1, glon_sub, &status);
  fits_write_col(fptr, TDOUBLE, 7, rowid, 1, 1, glat_sub, &status);

  fits_write_col(fptr, TFLOAT,  8, rowid, 1, 1, fd_ang,  &status);
  fits_write_col(fptr, TFLOAT,  9, rowid, 1, 1, pos_ang, &status);
  fits_write_col(fptr, TFLOAT, 10, rowid, 1, 1, par_ang, &status);
  fits_write_col(fptr, TFLOAT, 11, rowid, 1, 1, tel_az,  &status);
  fits_write_col(fptr, TFLOAT, 12, rowid, 1, 1, tel_zen, &status);

  fits_write_col(fptr, TFLOAT, 13, rowid, 1, 384, dat_freq, &status);
  fits_write_col(fptr, TFLOAT, 14, rowid, 1, 384, dat_wts,  &status);
  fits_write_col(fptr, TFLOAT, 15, rowid, 1, 384, dat_offs, &status);
  fits_write_col(fptr, TFLOAT, 16, rowid, 1, 384, dat_scl,  &status);

  fits_write_col(fptr, TBYTE, 17, rowid, 1, 48000, data, &status);
*/

