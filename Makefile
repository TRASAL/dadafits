SOURCE_ROOT ?= $(HOME)
BIN_DIR := $(SOURCE_ROOT)/install/bin

VERSION := $(shell git rev-parse HEAD )

# http://psrdada.sourceforge.net/
PSRDADA  ?= $(SOURCE_ROOT)/src/psrdada

FITS_CFLAGS ?= `pkg-config --cflags cfitsio`
FITS_LIBS ?= `pkg-config --libs cfitsio`

#OPTIMIZATION ?= -Ofast -march=native
OPTIMIZATION := -O0 -march=native -g -ggdb

INCLUDES := -I"$(PSRDADA)/src/"
DADA_DEPS := $(PSRDADA)/src/dada_hdu.o $(PSRDADA)/src/ipcbuf.o $(PSRDADA)/src/ipcio.o $(PSRDADA)/src/ipcutil.o $(PSRDADA)/src/ascii_header.o $(PSRDADA)/src/multilog.o $(PSRDADA)/src/tmutil.o $(PSRDADA)/src/fileread.o $(PSRDADA)/src/filesize.o

SOURCES := src/main.c src/downsample.c src/sb_util.c src/fits_io.c src/manipulate.c
HEADERS := src/dadafits_internal.h

dadafits: $(SOURCES) $(HEADERS)
	gcc $(FITS_CFLAGS) -o dadafits $(SOURCES) $(DADA_DEPS) $(FITS_LIBS) -lm -I"$(PSRDADA)/src" $(OPTIMIZATION) -DVERSION='"$(VERSION)"'

fits_dump: fits_dump.c
	gcc $(FITS_CFLAGS) $(OPTIMIZATION) -o fits_dump fits_dump.c $(FITS_LIBS)

longkey: longkey.c
	gcc $(FITS_CFLAGS) $(OPTIMIZATION) -o longkey longkey.c $(FITS_LIBS)
