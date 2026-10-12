"""
Generate a Flick reference spectrum of downward plane irradiance
just above the ocean surface, in W/m²/nm.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

full_spectrum = False
save_data = False
overwrite_current_standard = False

# Sun at one astronomical unit distance and solar zenith angle 60 degrees
time_pos = [2000, 4, 4, 9, 32, 20, 60.42, 5.27, 10e-9]

# Number of wavelengths used to compute transmittance of atmosphere
if full_spectrum:
    n_wls = 230
else:
    n_wls = 15

# Run Flick with AccuRT
c = flick.accurt_config()
c.set_streams(8)
c.set('ozone', 0.003)
c.set('pressure',101325)
c.set('temperature',288.15)
c.set('aerosol_od', 0.08)
c.set('cloud_liquid', 0)
c.set('detector_type', 'plane_irradiance')
c.set('detector_height', 0.1)
c.set('cdom_440', 0.1)
c.set('chl_concentration', 1e-6)
c.set('nap_concentration', 1e-3)

wl = np.linspace(290e-9, 1040e-9, n_wls)
c.set("wavelengths",wl)
c.set('toa_solar_multiplication',time_pos)
E = flick.run("accurt flick_tmp/config")

if save_data:
    file_name = 'surface_reference.txt' 
    flick.save_two_columns(radiation,file_name)
    
    # Insert header
    with open(file_name, 'r') as file:
        original_content = file.read()
    with open('header.txt', 'r') as file:
        header = file.read()
    with open(file_name, 'w') as file:
        file.write(header + '\n' + original_content)

if overwrite_current_standard:
    shutil.copy(file_name,'../../../radiator/.')

# Show spectrum
Er = flick.run('radiator surface_reference')
plt.plot(Er[:,0],Er[:,1],label='Current reference')
plt.plot(E[:,0],E[:,1],label='Suggested reference')
plt.plot(wl,np.interp(wl,E[:,0],E[:,1]),'.',label='calculated wavelengths')
plt.xlabel('Wavelength [m]')
plt.ylabel(r'Irradiance [W$\,$m$^{-2}\,$m$^{-1}$]')
plt.grid()
plt.legend()

if __name__ == "__main__":
    plt.show()
