"""
Compute conical irradiance for a detector with varying acceptance
half-angles and polar orientation angles relative to the Sun.

See also README.txt

Unless otherwise specified, SI (MKS) base units are used, with angles given in
degrees.

This script works only with link to the external AccuRT kernel.

"""

import numpy as np
import make_config
from make_config import c, sys, os, flick

""" See flick_tmp/config for parameter documentation """

c.set('aerosol_od', 0.08)
c.set('cloud_liquid', 0)
c.set('detector_height', 0.1)
c.set('cdom_440', 0.1)
c.set('chl_concentration', 1e-6)
c.set('nap_concentration', 1e-3)
c.set("gases",["o3", "o2", "h2o"])
c.set("pure_water_volume_fraction",0)

c.set("detector_radiance_distribution_override",[30,30])
c.set("use_AccuRT_kernel_if_exists","true")

""" Sun at one astronomical unit distance and solar zenith angle 60
degrees """
time_pos = [2000, 4, 4, 9, 32, 20, 60.42, 5.27]
c.set('toa_solar_multiplication',time_pos+[10e-9])

wl0 = 440e-9
c.set("wavelengths",[wl0]) #[wl0-3e-9, wl0, wl0+3e-9])

c.set_streams(16)

compute_all_angles = True

if compute_all_angles:
    orientation_angles = [0]#[0, 90, 180]
    acceptance_half_angles = np.linspace(1, 90, 2)
else:
    orientation_angles = [90]
    acceptance_half_angles = [10, 90] 
    
for oa in orientation_angles:
    F = []
    for aha in acceptance_half_angles:
        c.set("detector_conical_radiation_override", [oa, 90, aha])
        irradiance = flick.run("accurt flick_tmp/config")[1]*1e-9 # [W/m2/nm]  
        F.append([aha, irradiance])

    flick.save_two_columns(F, f"{oa}.txt")
    
