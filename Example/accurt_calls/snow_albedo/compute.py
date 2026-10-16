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

c.set('aerosol_od', 0.05)
c.set('cloud_liquid', 0)
c.set("snow_ice", 1)
c.set("snow_radius", 1e-3)
c.set("snow_impurity_names", "EIK1")
c.set("snow_impurity_scaling_factors", 0.1)
c.set("source_zenith_angle", 80)
c.set("wavelengths",np.linspace(300e-9, 900e-9, 2))
c.set_streams(16)

simulate_horizontal = True
simulate_tilted = False

def run():
    return flick.run("accurt flick_tmp/config")

if simulate_horizontal:
    albedo = run()
    albedo[:, 0] *= 1e9  # m to nm
    flick.save_two_columns(albedo, "albedo_horizontal.txt")

if simulate_tilted:
    c.set("detector_radiance_distribution_override", [20, 15])
    c.set("reference_detector_height", 120e3)

    c.set("detector_orientation", "up")
    c.set("detector_conical_radiation_override", [3, 180, 90])
    downward = run()

    c.set("detector_orientation", "down")
    c.set("detector_conical_radiation_override", [177, 0, 90])
    upward = run()

    albedo = np.column_stack((
        downward[:, 0] * 1e9,  # Wavelength [nm]
        upward[:, 1] / downward[:,1]
    ))
    flick.save_two_columns(albedo, "albedo_tilted.txt")
