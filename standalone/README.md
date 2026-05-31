# csdr-eti Standalone Tool

This directory contains a standalone command-line tool version of the csdr-eti DAB ETI decoder.

## Overview

The standalone tool reads raw IQ data (float 32bit) from a file or stdin and outputs ETI (Encoded Transport Interface) data to a file or stdout. Unlike the csdr module, this tool can be used independently without the full csdr pipeline.

## Building

### Prerequisites

- csdr >= 0.18
- FFTW3 library
- C++17 compatible compiler
- CMake >= 3.0

### Compilation

```bash
cd standalone
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
sudo make install
```

### Usage

Note: the common input format *u8* is not accepted, so we need to convert it with sox to *float32*

```bash
rtl_sdr -f 215072000 -s 2048000 -w 1700000 - | sox -r 2048000 -t u8 /dev/stdin -t f32 /dev/stdout | csdr-eti-standalone | dablin_gtk
```
(for 10D with 1,7MHz bandwidth and output to `dablin_gtk`)
