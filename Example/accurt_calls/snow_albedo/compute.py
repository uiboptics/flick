"""
Compute the spectral albedo including the configuration parameters set in
make_config.py.

See also README.txt

Unless otherwise specified, SI (MKS) base units are used, with angles given in
degrees.
"""

import numpy as np
import make_config
from make_config import c, sys, os, flick

""" See flick_tmp/config for parameter documentation """

c.set('aerosol_od', 0.1)
c.set('cloud_liquid', 1e-4)
c.set("snow_ice", 1)
c.set("snow_radius", 1e-3)
c.set("snow_impurity_names", "EIK1")
c.set("snow_impurity_scaling_factors", 1)
c.set("source_zenith_angle", 85)
c.set("wavelengths",np.linspace(300e-9, 900e-9, 2))
""" simulte horizontal sensor face """
"""
albedo = flick.run("accurt flick_tmp/config")
albedo[:,0] *= 1e9 # To nm
flick.save_two_columns(albedo,"albedo.txt")
"""

simulate_tilted_sensor_face = True

if simulate_tilted_sensor_face:
    c.set_streams(16)
    c.set("detector_radiance_distribution_override",[20,20])
    c.set('reference_detector_height', 120e3)

    c.set('detector_orientation','up')
    c.set("detector_conical_radiation_override",[0,0,90])
    F_d = flick.run("accurt flick_tmp/config")
"""
    c.set('detector_orientation','down')
    c.set("detector_conical_radiation_override",[-3,0,90])
    F_u = flick.run("accurt flick_tmp/config")

    albedo = F_u/F_d
    albedo[:,0] *= 1e9 # To nm
    flick.save_two_columns(albedo,"albedo_tilted.txt")
"""
