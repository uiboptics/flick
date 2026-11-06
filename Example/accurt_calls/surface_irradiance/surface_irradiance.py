"""
Calculate and plot the downward UV irradiance spectrum at the
surface for a specified spectral bandwidth. The simulation does not
include an ocean.

See also `README.md`.

"""
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from make_config import c
import os
import flick

time_pos = [2022, 5, 20, 11, 24, 0, 59.88, 5.7]
c.set('toa_solar_multiplication',time_pos + [3e-9])
c.set('wavelengths',np.linspace(280e-9,400e-9,10))
c.set("aerosol_od", 0)
c.set("cloud_liquid", 0)
c.set("ozone", 0.003)
c.set("water_vapor", 10)
c.set("detector_height",1)
c.set("pure_water_volume_fraction",0)
c.set_streams(8)

E = flick.run("accurt flick_tmp/config")

fig, ax = plt.subplots()
fig.set_size_inches(5,4)
ax.plot(E[:,0]*1e9,E[:,1]*1e-9)
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('Surface irradiance [W m$^{-2}$ nm$^{-1}$]')
ax.grid()

if __name__ == "__main__":
    plt.show()
