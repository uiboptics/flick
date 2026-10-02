"""
Compute scalar irradiance, plane irradiance, and nadir radiance as
measured by three different Ramses sensors, all oriented vertically
upwards above snow-covered sea ice.

See also README.txt

Unless otherwise specified, SI (MKS) units are used, with angles given in
degrees.
"""

import os
import numpy as np
from dataclasses import dataclass
from make_config import f
import flick

def set_parameters(detector_type):
    """ See flick_tmp/config for parameter documentation """

    f.set('aerosol_od', 0.1)
    f.set('cloud_liquid', 1e-4)

    f.set("snow_ice", 1)
    f.set("snow_radius", 1e-3)
    f.set("snow_impurity_names", "EIK1")
    f.set("snow_impurity_scaling_factors", 1)

    f.set('ice_depths',1)
    f.set('ice_bubble_fraction',[0.005])
    f.set('ice_brine_fraction',[0.02])

    f.set('cdom_440',0.1)        
    f.set('chl_concentration',0.1e-6)        
    f.set('nap_concentration',0.1e-3)  
   
    f.set('detector_type', detector_type)
    if detector_type == 'radiance':
        f.set_n_angles(round(16**1.6))
    else:
        f.set_n_angles(round(8**1.6))

@dataclass
class spacetime_point:
    time_point_utc: str
    latitude: float
    longitude: float

def change_units(detector_type, value):
    value[:,0] *= 1e9 # wavelength [nm]
    value[:,1] *= 1e-9 # radiation [per nm]
    return value
    
def save(radiation, name):
    r = change_units(d,radiation)
    if not os.path.exists('output'):
        os.makedirs('output')
    np.savetxt(f'output/{name}.txt', r, fmt=['%6.2e ','%8.3e'])

  
if __name__ == "__main__":
    n_wl = 7 # Increase to improve spectral accuracy
    wavelengths = np.linspace(320e-9, 940e-9, n_wl)
    wl_band_width = 10e-9 # Ramses band width

    # More spacetimes may be added here. Default is North Pole only. 
    spacetimes = [
        spacetime_point('2026 08 29 08 10 0',90,0),
        #spacetime_point('2026 08 30 08 23 0',88+51.6/60,48+4.4/60)
    ]
    detector_types = ['plane_irradiance','scalar_irradiance','radiance']
    for s in spacetimes:
        for d in detector_types:
            name = f'{s.time_point_utc[0:10]} {d}'
            print(f'running for {name}')
            set_parameters(d)
            radiation = f.spectrum(wavelengths, wl_band_width, s.time_point_utc,
                                   s.latitude, s.longitude)
            save(radiation,name)
