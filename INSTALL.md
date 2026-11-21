# Build dependencies

## Requirements

| Dependency | Purpose |
|---|---|
| C compiler (`cc`, `clang`, `gcc`) | builds the C11 sources |
| [SHTns](https://nschaeff.bitbucket.io/shtns/) (with FFTW3) | spherical-harmonic transforms |
| NetCDF-C | reads surface boundary data, writes monthly output |
| LAPACK | matrix inversion (the Accelerate framework on macOS) |
| Open MPI (`mpicc`, `mpirun`) | MPI builds and runs |
| `pkg-config` | locates libraries and link flags |

## Installing on macOS

PlaSiC's horizontal scalar and vector spherical-harmonic transforms are performed by
[SHTns](https://nschaeff.bitbucket.io/shtns/), which internally uses FFTW3. This
project has been tested with SHTns 3.7.5.

Below is the complete procedure for installing SHTns into a persistent user
directory on macOS. The installation prefix is `$HOME/.local/shtns-3.7.5`.

### 1. Install the base dependencies

Install [Homebrew](https://brew.sh/) first, then FFTW3, NetCDF-C and `pkg-config`:

```bash
brew install fftw netcdf pkg-config
```

For MPI builds and runs, also install Open MPI:

```bash
brew install open-mpi
```

### 2. Set the version, prefix and work directory

```bash
SHTNS_VERSION="3.7.5"
SHTNS_PREFIX="$HOME/.local/shtns-${SHTNS_VERSION}"
SHTNS_WORK="$HOME/Downloads/shtns-${SHTNS_VERSION}-build"

mkdir -p "$SHTNS_WORK"
cd "$SHTNS_WORK"
```

`SHTNS_PREFIX` is the final installation directory; `SHTNS_WORK` only holds the
downloaded source and intermediate build files.

### 3. Download SHTns 3.7.5

```bash
curl -fL \
  "https://files.pythonhosted.org/packages/35/60/241291c0a737d269efd0ae02a488172721f1ad4b226729613cf8d12d3a4f/shtns-3.7.5.tar.gz" \
  -o "shtns-3.7.5.tar.gz"
```

### 4. Unpack and configure

```bash
tar -xzf "shtns-3.7.5.tar.gz"
cd "shtns-3.7.5"

FFTW_PREFIX="$(brew --prefix fftw)"

CPPFLAGS="-I${FFTW_PREFIX}/include" \
LDFLAGS="-L${FFTW_PREFIX}/lib" \
./configure \
  --prefix="$SHTNS_PREFIX" \
  --disable-openmp
```

In PlaSiC's teaching-oriented MPI path, every process performs its own complete
SHTns transform. SHTns' internal OpenMP is therefore disabled here, so that MPI
processes and OpenMP threads do not compete for the same CPUs.

### 5. Compile SHTns

```bash
make -j"$(sysctl -n hw.logicalcpu)"
```

A successful build produces the static library:

```bash
ls -lh libshtns.a
```

### 6. Install into the persistent prefix

In a non-GPU build, `make install` from the SHTns 3.7.5 source package may still try
to copy CUDA headers that are not present. PlaSiC only needs `shtns.h` and
`libshtns.a`, so install those two files directly:

```bash
mkdir -p "$SHTNS_PREFIX/include"
mkdir -p "$SHTNS_PREFIX/lib"

install -m 0644 shtns.h \
  "$SHTNS_PREFIX/include/shtns.h"

install -m 0644 libshtns.a \
  "$SHTNS_PREFIX/lib/libshtns.a"
```

Check the result:

```bash
ls -lh \
  "$SHTNS_PREFIX/include/shtns.h" \
  "$SHTNS_PREFIX/lib/libshtns.a"
```

### 7. Let PlaSiC find SHTns permanently

Add to `~/.zshrc`:

```bash
export SHTNS_PREFIX="$HOME/.local/shtns-3.7.5"
```

Then reload the shell configuration:

```bash
source ~/.zshrc
```

Verify the variable:

```bash
echo "$SHTNS_PREFIX"
```

which should print something like:

```text
/Users/your-name/.local/shtns-3.7.5
```

PlaSiC's Makefile reads this variable, so later builds do not need it passed on the
command line. If you prefer not to edit `~/.zshrc`, specify it for each build
instead:

```bash
make -C src SHTNS_PREFIX="$HOME/.local/shtns-3.7.5"
```

### 8. Build and test PlaSiC

From the project root, run a strict-warning build first:

```bash
make -C src strict
```

Build the serial model:

```bash
make -C src MPI=0
```

Build the two-process MPI model:

```bash
make -C src MPI=1 NPRO=2
```

The Makefile locates FFTW3 with `pkg-config` when possible; `SHTNS_CFLAGS`/`SHTNS_LIBS`
and `FFTW_CFLAGS`/`FFTW_LIBS` can also be overridden directly. Because `libshtns.a`
is linked statically, running PlaSiC does not require setting `DYLD_LIBRARY_PATH`.
