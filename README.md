
Flick is a radiative-transfer modelling framework for generating the
optical properties of the Earth’s atmosphere and ocean. It includes
Monte Carlo methods for radiance simulations and provides input and
output functionality for the AccuRT model. Flick uses the AccuRT
kernel to solve the radiative-transfer equation. Because AccuRT is
commercial software, it is not included in this release; these
capabilities require a separate AccuRT licence and installation.


#### Clone source files from GitHub with

 - *git clone https://github.com/uiboptics/flick.git*

#### For macOS and Linux, install and run with 

 - *cd flick* - to enter the flick directory

 - *make* - to compile, build, and test

 - *make install* - to make the already-built `flick` command and
   Python module available from any directory. For a current-user
   installation, use `make install BINDIR="$HOME/.local/bin"`. For a
   system-wide installation, use `sudo make install
   BINDIR=/usr/local/bin PYTHON_INSTALL_FLAGS=`. On Linux
   distributions that use an externally managed Python (PEP 668), the
   current-user installation automatically passes pip's required
   `--break-system-packages` opt-in together with `--user`; it does
   not write to the system Python directories.

 - *make uninstall* - to remove the installed launcher and Python
   module (use the same `BINDIR` passed to `make install`, if any)

 - *flick help* - for further documentation and usage examples.

You will need a c++ compiler (e.g. clang++ or g++) already installed.
On most Linux variants "sudo apt install build-essential" should
provide the needed compiler and tools.

#### Dependencies
 
The linear algebra library
[*Eigen*](https://gitlab.com/libeigen/eigen.git) (MPL 2.0) is
automatically cloned into the flick/external directory during
installation.

#### Name and license

Flick is named after the flickering of a flame and the unit
occasionally used for spectral radiance. It is distributed under the
MIT license, https://docs.openmc.org/en/stable/license.html.
