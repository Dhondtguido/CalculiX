# CalculiX Binary

This directory contains the compiled CalculiX binary used for benchmark testing.

- **Version**: 2.21
- **Size**: 12MB
- **Compiled**: Ubuntu 24.04, gfortran
- **Dependencies**: SPOOLES, ARPACK, OpenBLAS

## Usage
```bash
./ccx_2.21 input_file
```

## Rebuilding
To rebuild from source:
```bash
cd /work/src
make
```
