SOURCE_ROOT ?= $(HOME)
BIN_DIR := $(SOURCE_ROOT)/install/bin

VERSION := $(shell git rev-parse HEAD )

# http://psrdada.sourceforge.net/
PSRDADA  := $(SOURCE_ROOT)/src/psrdada

OPTIMIZATION := -O1
# OPTIMIZATION := -Ofast -march=native 

INCLUDES := -I"$(PSRDADA)/src/"
DADA_DEPS := $(PSRDADA)/src/dada_hdu.o $(PSRDADA)/src/ipcbuf.o $(PSRDADA)/src/ipcio.o $(PSRDADA)/src/ipcutil.o $(PSRDADA)/src/ascii_header.o $(PSRDADA)/src/multilog.o $(PSRDADA)/src/tmutil.o $(PSRDADA)/src/fileread.o $(PSRDADA)/src/filesize.o

dadafits: main.c
	gcc -o dadafits main.c `pkg-config --cflags --libs cfitsio` $(DADA_DEPS) -I"$(PSRDADA)/src" $(OPTIMIZATION) -DVERSION='"$(VERSION)"' 
fits_example:
	gcc -o fits_example fits_example.c `pkg-config --cflags --libs cfitsio`
