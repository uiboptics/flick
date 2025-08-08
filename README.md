
Flick is a radiative transfer modeling framework that currently supports generation of optical properties for the Earth’s atmosphere and ocean. It also includes Monte Carlo implementations for radiance simulations. In addition, Flick provides input–output functionalities for users of AccuRT.      

#### Clone from GitHub with

 - *git clone https://github.com/uiboptics/flick.git*

#### For macOS and Linux, install and run with 

 - *cd flick* - to enter the flick directory

 - *make* - to compile and build

 - *flick help* - for further documentation and usage examples.

You will need a c++ compiler (e.g. clang++ or g++) already installed.

#### Dependencies
 
The linear algebra library [*Eigen*](https://gitlab.com/libeigen/eigen.git) (MPL 2.0) is automatically cloned into the flick/external directory during installation.

#### Name and license

Flick is named after the flickering of a flame and the unit occasionally used for spectral radiance. It is distributed under the MIT license, https://docs.openmc.org/en/stable/license.html.
