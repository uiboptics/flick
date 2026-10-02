"""Compute the spectral albedo including the configurations set in
make_config.py.

See also README.txt

Unless otherwise specified, SI (MKS) units are used, with angles given in
degrees.

"""

import os
import numpy as np
from make_config import f

""" See flick_tmp/config for parameter documentation """
f.set('aerosol_od', 0.1)
f.set('cloud_liquid', 1e-4)
f.set("snow_ice", 1)
f.set("snow_radius", 1e-3)
f.set("snow_impurity_names", "EIK1")
f.set("snow_impurity_scaling_factors", 1)

n_wl = 7 # Increase to improve spectral resolution
wavelengths = np.linspace(300e-9, 900e-9, n_wl) 
solar_zenith_angle = 60  
albedo = f.spectrum(wavelengths, solar_zenith_angle)

albedo[:,0] *= 1e9 # [nm]
if not os.path.exists('output'):
    os.makedirs('output')
np.savetxt('output/computed_albedo.txt', albedo, fmt=['%6.2f ','%8.3e'])
