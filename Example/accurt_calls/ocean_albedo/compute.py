"""
Compute the spectrally integrated surface albedo of the ocean for
varying volume fractions of nano bubbles.

See also README.txt

Unless otherwise specified, SI (MKS) base units are used, with angles given in
degrees.

This script works only with link to the external AccuRT kernel.

"""

import numpy as np
import time
import make_config
from make_config import c, sys, os, flick

""" See flick_tmp/config for parameter documentation """

c.set('aerosol_od', 0.08)
c.set('cloud_liquid', 0)
c.set('detector_height', 1)
c.set('cdom_440', 0.1)
c.set('chl_concentration', 1e-6)
c.set('nap_concentration', 1e-3)
c.set("gases",["o3", "o2", "h2o"])
c.set("use_AccuRT_kernel_if_exists","true")
c.set("ti_names","small_marine_bubbles")
c.set("gas_spectral_region","solar")

""" Sun at one astronomical unit distance and solar zenith angle 60
degrees """
time_pos = [2000, 4, 4, 9, 32, 20, 60.42, 5.27]
c.set('toa_solar_multiplication',time_pos+[10e-9])

wavelengths = np.geomspace(300e-9, 2500e-9, 30)
n_cores = max(1, (os.cpu_count() or 1) - 1)
print(f'Using {n_cores} CPU cores')
bubble_volume_fractions = np.geomspace(1e-8,1e-1,n_cores)
c.set("wavelengths", wavelengths)

c.set_streams(8)

def integrated_irradiance(direction):
    c.remove_parallel_directories()
    c.set("detector_orientation", direction)

    for i, bvf in enumerate(bubble_volume_fractions):
        c.copy_config_to_parallel(i)
        c.set("ti_volume_fractions", bvf)
    spectra = flick.run("accurt --parallel flick_tmp")
    integrated = np.empty(len(bubble_volume_fractions))
    for i, s in enumerate(spectra):
        integrated[i] = np.trapezoid(s[:, 1], s[:, 0])
    return integrated

calculation_start = time.perf_counter()
albedo = integrated_irradiance("down") / integrated_irradiance("up")
calculation_duration = time.perf_counter() - calculation_start
print(f"Calculation of {n_cores} volume fractions completed in {calculation_duration:.2f} seconds")

output = np.column_stack((bubble_volume_fractions, albedo))
flick.save_two_columns(output, "albedo_versus_volume_fraction.txt")

