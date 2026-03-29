# Mesosim

Mesosim is a Monte Carlo simulation tool for mesoscale modeling of atomic systems governed by bond-breaking physics, inspired by LAMMPS and SPPARKS. It supports flexible input scripts, various crystal structures, and multiple simulation flavors (KMC, MC). Mesosim is designed for research and educational use in computational materials science.

## Features
- Flexible input script format (similar to LAMMPS/SPPARKS)
- Supports FCC, BCC, and SC crystal structures
- Sheet, cluster, and file-based geometry initialization
- N-ary atomic systems
- System evolution via diffusion and evaporation/dissolution governed by a broken-bond model
- M-nearest neighbor shell energy specification
- CSV and XYZ data recording formats
- Multiple simulation flavors: Kinetic Monte Carlo (KMC), Monte Carlo (MC)

## Input File Format
Input scripts are read line-by-line. Repeated commands are overwritten by the most recent instance. Simulation starts after the full input file is read. For a full list of commands and details, see `docs/Commands.md`.

Example input:
```
systemsize 20 20 20
struct FCC
geometry cluster 8
atomtype A B
composition 0.5 0.5
dissolution true false
nnlevels 1
nne 0.1 0.2 0.3 0.1 0.2 0.3
run time 1000
flavor KMC
```

### Key Commands
- `systemsize NX NY NZ` — Set system size
- `struct FCC|BCC|SC` — Crystal structure
- `geometry (sheet N|cluster R|file path)` — Geometry initialization
- `atomtype A B C ...` — Atom types
- `composition xA xB xC ...` — Atomic fractions
- `dissolution true|false ...` — Atom dissolution flags
- `nnlevels N` — Number of neighbor shells
- `nne eAA eAB ...` — Neighbor energies
- `run time|iteration value` — Simulation end condition
- `flavor KMC|MC` — Simulation algorithm
- `potential initial [ramp max]` — Electric potential sweep
- `seed default|random|int` — Random seed
- `output [type] filename ...` — Output file name

## Building Mesosim

Mesosim uses a flexible Makefile supporting multiple build types and traditional Make variables to control compilation.

### Build Types
- **release** (default): Optimized build for production
- **debug**: Debug build with symbols and no optimization
- **test**: Build with unit tests
- **dbtest**: Debug build with unit tests
- **docs**: Build documentation tools

### Building
To build a specific type, use:
```
make [release|debug|test|dbtest|docs]
```
For example, to build the release version (default):
```
make release
```
Or simply:
```
make
```

### Overriding Compilation Variables
You can override compile variables by passing them on the command line:
```
make CC=clang CFLAGS="-O2 -Wall" SHELL=/bin/bash
```
This will use `clang` as the compiler, set `CFLAGS` to `-O2 -Wall`, and use `/bin/bash` as the shell.

### Clean
To remove all build outputs:
```
make clean
```

## Running Simulations
After building, run Mesosim with your input file:
```
./build/release/mesosim input/your_input.in
```
(Output location and executable name may vary by build type and platform.)

## Directory Structure
- src/ — Source code
- include/ — Header files
- build/ — Build outputs (debug, release, test, etc.)
- docs/ — Documentation
- test/ — Unit and integration test files

## Documentation
- docs/Commands.md — Input file format and command details

## License
[MIT](https://choosealicense.com/licenses/mit/)
