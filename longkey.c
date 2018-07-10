#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#define MAX_LONGSTRING 20000

int main(int argc, char *argv[]) {
  int status;
  fitsfile *fptr;
  char longstring[MAX_LONGSTRING];
  int c;

  // create a very long string
  for (c=0; c<MAX_LONGSTRING-2; c++) {
    longstring[c] = (c % ('Z'-'A')) + 'A';
  }
  longstring[MAX_LONGSTRING-1] = '\0';

  fits_create_file(&fptr, "example.fits(template.txt)", &status);
 
  fits_movabs_hdu(fptr, 1, NULL, &status);
  fits_write_date(fptr, &status);
  fits_write_key_longwarn (fptr, &status);
  fits_write_key_longstr (fptr, "parset", longstring, NULL, &status);
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

