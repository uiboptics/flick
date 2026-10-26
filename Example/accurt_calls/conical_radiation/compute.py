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
c.set('detector_height', -5)
c.set('cdom_440', 0.1)
c.set('chl_concentration', 1e-6)
c.set('nap_concentration', 1e-3)
c.set("gases",["o3", "o2", "h2o"])
c.set("use_AccuRT_kernel_if_exists","true")

""" Sun at one astronomical unit distance and solar zenith angle 60
degrees """
time_pos = [2000, 4, 4, 9, 32, 20, 60.42, 5.27]
c.set('toa_solar_multiplication',time_pos+[10e-9])

wl0 = 440e-9
c.set("wavelengths",[wl0])

c.set_streams(24)

compute_all_angles = False

orientation_angles = [0, 90, 180]
if compute_all_angles:
    acceptance_half_angles = np.linspace(1, 90, 30)
else:
    acceptance_half_angles = [10, 90] 

c.remove_parallel_directories()
conical = np.empty((len(acceptance_half_angles), len(orientation_angles)))
for i, acceptance_half_angle in enumerate(acceptance_half_angles):
    for j, orientation_angle in enumerate(orientation_angles):
        c.copy_config_to_parallel(j)
        c.set("detector_conical_radiation_override", [orientation_angle,
                                                      90,
                                                      acceptance_half_angle])
    parallel = flick.run("accurt --parallel flick_tmp")
    if len(parallel) != len(orientation_angles):
        raise RuntimeError(
            f"Expected {len(orientation_angles)} result matrices, "
            f"but AccuRT returned {len(parallel)}"
        )
    for j, result in enumerate(parallel):
        conical[i, j] = result[0, 1] * 1e-9  # [W/m2/nm]

for j, orientation_angle in enumerate(orientation_angles):
    output = np.column_stack((acceptance_half_angles, conical[:, j]))
    flick.save_two_columns(output, f"{orientation_angle:03d}.txt")
