# dadafilterbank

Connect to a [PSRdada](http://psrdada.sourceforge.net/) ringbuffer, optionally downsample and compress, and write out the data in [FITS](https://fits.gsfc.nasa.gov/fits_home.html) format.

This program is part of the data handling pipeline for the AA-ALERT project.
See [dadatrigger](https://github.com/AA-ALERT/dadatrigger) for an introduction and dataflow schema.

# Usage

```bash
 $ dadafits -k <hexadecimal key> -l <logfile> -c <science_case> -m <science_mode> -b <padded_size> -t <template> -d <output_directory> -S <synthesized beam table> -s <synthesize these beams>
```

Command line arguments:
 * *-k* Set the (hexadecimal) key to connect to the ringbuffer.
 * *-l* Absolute path to a logfile (to be overwritten)
 * *-n* Prefix for the fitlerbank output files
 * *-c* move to ringbuffer header
 * *-m* move to ringbuffer header
 * *-b* move to ringbuffer header
 * *-t* determined by case/mode?
 * *-d*
 * *-S* 
 * *-s*

# Modes of operation

## As part of the real time pipeline 

These modes are for archiving data; the program can be run as part of the realtime pipeline.

### Science modes

The program implements different modes:
- mode 0: Stokes I + IAB (coherent beams, so only one tied array beam)
- mode 2: Stokes I + TAB (12 tied array beams)

In these modes data is also:
* integrated over time to reduce sample rate to 500 samples per 1.024 seconds
* summed over frequencies to reduce total number of frequencies to 384

## Science cases

The data input rate is set per science case.
Supported cases:
- case 3: 12500 samples per second
- case 4: 25000 samples per second

## As part of an event-based postprocessing step

These modes are for analysing event data, and are not optimized for real time use.

### Science modes

The program implements different modes:
- mode 1: Stokes IQUV + IAB
- mode 3: Stokes IQUV + TAB

## Science cases

The data input rate is set per science case.
Supported cases:
- case 3: 12500 samples per second
- case 4: 25000 samples per second

# The ringbuffer

## Header block

Metadata is read from the PSRdada header block.
Note that some of the metadata available in the header block is ignored, due to code constraints and optimizations.
For values that should be present see the table below.

|header key| description | notes | units |
|----------|-------------|-------|-------|
| MIN\_FREQUENCY | lowest frequency                           |                              | |
| BW             | Bandwidth of a frequency channel           |                              | |
| PADDED\_SIZE   | Length of the fastest dimension of the data array |                       | |
| SCIENCE\_CASE  | Mode of operation of ARTS, determines data rate   |  TODO                 | |
| SCIENCE\_MODE  | Mode of operation of ARTS, determines data layout |     TODO              | |

## Data block

A ringbuffer page is interpreted as an array of Stokes I: [NTABS, NCHANNELS, padded\_size]
Array padding along the fastest dimension is implemented to facilitate memory copies.

# FITS output files

# Building

To connect to the PSRDada ring buffer, we depend on some object files that can be obtained when compiling PSRDada.
The location of these files is assumed to be in the **PSRDADA** directory.
Alternatively, set **SOURCE\_ROOT** such that the files are in **SOURCE\_ROOT/src/psrdada**.

Building is then done using the Makefile:
```bash
  make
```

# Contributers

Jisk Attema, Netherlands eScience Center

# NOTES

