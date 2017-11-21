SOURCE_ROOT ?= $(HOME)
BIN_DIR := $(SOURCE_ROOT)/install/bin

VERSION := $(shell git rev-parse HEAD )

# http://psrdada.sourceforge.net/
PSRDADA  ?= $(SOURCE_ROOT)/src/psrdada

FITS_CFLAGS ?= `pkg-config --cflags cfitsio`
FITS_LIBS ?= `pkg-config --libs cfitsio`

#OPTIMIZATION ?= -Ofast -march=native
OPTIMIZATION := -DDRY_RUN

INCLUDES := -I"$(PSRDADA)/src/"
DADA_DEPS := $(PSRDADA)/src/dada_hdu.o $(PSRDADA)/src/ipcbuf.o $(PSRDADA)/src/ipcio.o $(PSRDADA)/src/ipcutil.o $(PSRDADA)/src/ascii_header.o $(PSRDADA)/src/multilog.o $(PSRDADA)/src/tmutil.o $(PSRDADA)/src/fileread.o $(PSRDADA)/src/filesize.o

dadafits: main.c
	gcc $(FITS_CFLAGS) -o dadafits main.c $(DADA_DEPS) $(FITS_LIBS) -lm -I"$(PSRDADA)/src" $(OPTIMIZATION) -DVERSION='"$(VERSION)"'

fits_example: fits_example.c
	gcc $(FITS_CFLAGS) $(OPTIMIZATION) -o fits_example fits_example.c $(FITS_LIBS)
