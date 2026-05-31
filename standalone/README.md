# csdr-eti Standalone Tool

This directory contains a standalone command-line tool version of the csdr-eti DAB ETI decoder.

## Overview

The standalone tool reads raw IQ data from a file or stdin and outputs ETI (Encoded Transport Interface) data to a file or stdout. Unlike the csdr module, this tool can be used independently without the full csdr pipeline.

## Building

### Prerequisites

- csdr >= 0.18
- FFTW3 library
- C++17 compatible compiler
- CMake >= 3.0

### Compilation

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
sudo make install
