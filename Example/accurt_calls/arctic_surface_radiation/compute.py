"""
Compute scalar irradiance, plane irradiance, and nadir radiance as
measured by three different Ramses sensors, all oriented vertically
upwards above snow-covered sea ice.

See also README.txt

Unless otherwise specified, SI (MKS) base units are used, with angles given in
degrees.
"""

import numpy as np
from dataclasses import dataclass
from make_config import c, sys, os, flick

""" See flick_tmp/config for parameter documentation """

c.set('aerosol_od', 0.1)
c.set('cloud_liquid', 1e-4)

c.set("snow_ice", 1)
c.set("snow_radius", 1e-3)
c.set("snow_impurity_names", "EIK1")
c.set("snow_impurity_scaling_factors", 0.1)

c.set('ice_depths',1)
c.set('ice_bubble_fraction',[0.005])
c.set('ice_brine_fraction',[0.02])

c.set('cdom_440',0.1)        
c.set('chl_concentration',0.1e-6)        
c.set('nap_concentration',0.1e-3)  

def set_detector(detector_type):
    c.set('detector_type', detector_type)
    if detector_type == 'radiance':
        c.set_streams(16)
    else:
        c.set_streams(8)

@dataclass
class spacetime_point:
    time_point_utc: str
    latitude: float
    longitude: float

    #def outstring(self):
    #    return f'{self.time_point_utc} {self.latitude} {self.longitude}'

    #def __str__(self):
    #    return self.outstring()

def change_units(detector_type, value):
    value[:,0] *= 1e9 # wavelength [nm]
    value[:,1] *= 1e-9 # radiation [per nm]
    return value
    
def save(radiation, name):
    r = change_units(d,radiation)
    flick.save_two_columns(radiation,f"{name}.txt")
  
if __name__ == "__main__":
    n_wl = 7 # Increase to improve spectral accuracy
    wavelengths = np.linspace(320e-9, 940e-9, n_wl)
    wl_band_width = 10e-9 # Ramses band width

    # More spacetimes may be added here. Default is North Pole only. 
    spacetimes = [
        [2026, 8, 29, 8, 10, 0, 90, 0]
        #, [2026, 8, 30, 8, 23, 0, 88+51.6/60, 48+4.4/60]
                  ]
    detector_types = ['plane_irradiance','scalar_irradiance','radiance']
    for s in spacetimes:
        for d in detector_types:
            name = f"{" ".join(map(str, s))} {d}" 
            print(f'running for: {name}')
            set_detector(d)
            c.set('toa_solar_multiplication',s + [wl_band_width])
            radiation = flick.run("accurt flick_tmp/config")
            save(radiation,name)
