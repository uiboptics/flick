"""
Compute the spectral albedo including the configuration parameters set in
make_config.py.

See also README.txt

Unless otherwise specified, SI (MKS) base units are used, with angles given in
degrees.
"""

import numpy as np
from make_config import c, sys, os, flick

""" See flick_tmp/config for parameter documentation """

c.set('aerosol_od', 0.1)
c.set('cloud_liquid', 1e-4)
c.set("snow_ice", 1)
c.set("snow_radius", 1e-3)
c.set("snow_impurity_names", "EIK1")
c.set("snow_impurity_scaling_factors", 1)
c.set("source_zenith_angle", 60)
c.set("wavelengths",np.linspace(300e-9, 900e-9, 8))

albedo = flick.run("accurt flick_tmp/config")
albedo[:,0] *= 1e9 # To nm
flick.save_two_columns(albedo,"albedo.txt")

