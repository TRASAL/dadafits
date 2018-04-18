# dadafits

Connect to a [PSRdada](http://psrdada.sourceforge.net/) ringbuffer, optionally downsample and compress, and write out the data in [FITS](https://fits.gsfc.nasa.gov/fits_home.html) format.

This program is part of the data handling pipeline for the AA-ALERT project.
See [dadatrigger](https://github.com/AA-ALERT/dadatrigger) for an introduction and dataflow schema.

# Usage

```bash
 $ dadafits -k <hexadecimal key> -l <logfile> -t <template_directory> -d <output_directory> -S <synthesized beam table> -s <synthesize these beams>
```

Command line arguments:
 * *-k* Set the (hexadecimal) key to connect to the ringbuffer.
 * *-l* Absolute path to a logfile (to be overwritten)
 * *-n* Prefix for the fitlerbank output files
 * *-t* Template directory (defaults to the directory **templates** in the current working directory)
 * *-d* Output directory
 * *-S* Synthesized beam table
 * *-s* Selection of synthesized beams

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

### Science cases

The data input rate is set per science case.
Supported cases:
- case 3: 12500 samples per second
- case 4: 25000 samples per second

# The ringbuffer

## Header block

Metadata is read from the PSRdada header block.
Note that some of the metadata available in the header block is ignored, due to code constraints and optimizations.
For values that should be present see the table below.

|header key      | type    | units | description | notes |
|----------------|---------|-------|-------------|-------|
| MIN\_FREQUENCY | double  | Mhz                 | Center of lowest frequency band of observation    |  |
| BW             | double  | Mhz                 | Total bandwidth of observation                    |  |
| PADDED\_SIZE   | int     | bytes               | Length of the fastest dimension of the data array |  |
| SCIENCE\_CASE  | int     | 1                   | Mode of operation of ARTS, determines data rate   | Must be 3 or 4 |
| SCIENCE\_MODE  | int     | 1                   | Mode of operation of ARTS, determines data layout | Either 1,2,3, or 4 |
| RA\_HMS        | string  | HH:MM:SS.ssss       | Right ascension                                   | maps to RA |
| DEC\_HMS       | string  |+HH:MM:SS.ssss       | Declination                                       | maps to DEC |
| SOURCE         | string  | text                | Source name                                       | maps to SRC\_NAME |
| UTC\_START     | char    | YYYY-MM-DDTHH:MM:SS | Human readable timestamp of the start of the observation. | The program will silently modify the separators to conform to FITS standard. However, whitespace characters as in '2018-04-18 14:40:10' will not work | maps to DATE-OBS |
| MJD\_START     | double  | days since epoch    | Modified Julian Date                              | maps to STT\_IMJD and STT\_SMJD |
| LST\_START     | double  | degrees             | Local siderial time                               | maps to STT\_LST |
| AZ\_START      | float   | degrees             | Azimuth angle of telescope                        | set per row in the SUBINT binary table as TEL\_AZ, assumed constant over the run |
| ZA\_START      | float   | degrees             | Zenith angle of telescope                         | set per row in the SUBINT binary table as TEL\_ZEN, assumed constant over the run |

## Data block

For modes 0 and 2 (ie Stokes I data), a ringbuffer page is interpreted as an array of Stokes I: [NTABS, NCHANNELS, padded\_size]
Array padding along the fastest dimension is implemented to facilitate memory copies.

For modes 1 and 3 (ie Stokes IQUV), a ringbuffer page is an interleaved array: [tab, channel\_offset, sequence\_number, packet]
Where:
- tab ranges from 0 to 0 or 11 (modes IAB or TAB)
- channel\_offset ranges from 0 upto 383 (NCHANNELS/4 - 1)
- sequence\_number ranges from 0 upto 24
- packet is a direct copy of a UDP datapacket coming from the network, making up 8000 bytes

The packet itself is an array: [time, channel, polarization]
where:
- time runs from 0 to 499, to get actual time, *sequence\_number * 500* should be added
- channel runs from 0 to 4, to get actual channel, *channel\_offset * 4* should be added
- polarization stands for the 4 Stokes components, IQUV.

# FITS output files

Output files are created in the directory specified on the commandline.
A template is used for the FITS file and is selected based on science case and mode.
Templates are searched for in the **template** directory in the current working directory; or its location can be specified as a command line argument.

Data is stored one beam per file.

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
Leon Oostrum, ASTRON / UvA

# NOTES

