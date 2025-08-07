
Flick is a general radiative transfer modeling framework  

## 📦 Dependencies

This project depends on the following lools and libraries:

### Required:
- **C++20-compatible compiler** – e.g., g++ or clang++
- [Eigen](https://gitlab.com/libeigen/eigen.git) (MPL 2.0) – Included as a Git submodule


### Optional:
- **Python + matplotlib** – For optional plotting scripts

---

## 🔧 Building the Project

### Step 1: Clone the repository with submodules

```bash
 git clone --recurse-submodules https://github.com/uiboptics/flick.git  

### Step 2: Enter directory and compile
```bash
 cd flick
 make






Clone from GitHub with:  
 git clone https://github.com/uiboptics/flick.git  

For macOS and Linux:

 enter the flick directory  

 run 'make' to compile, build, and test  

 run 'flick help' for further documentation and usage examples  

You will need to have a c++ compiler (e.g. clang++ or g++) already installed.  


## Name

Flick is named after the flickering of a flame and the unit
occasionally used for spectral radiance, and it is distributed under
the MIT license.  
https://docs.openmc.org/en/stable/license.html


## Third-Party Libraries

This project includes the [Eigen](http://eigen.tuxfamily.org) C++
template library for linear algebra.

Eigen is included in the `external/eigen/` directory, and is licensed
under the **Mozilla Public License v2.0 (MPL-2.0)**.

See the full license in `external/eigen/COPYING.MPL2`.

Your use of this project must comply with both the **MIT License**
(for this code) and the **MPL-2.0** (for Eigen).
