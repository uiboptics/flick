"""
Compute the surface scalar irradiance above snow-covered sea ice,
the downwelling plane irradiance, the nadir radiance, and the radiance
distribution beneath the snow-covered sea ice.

See the flick_tmp/config file that will be
generated after the first run for documentation on all variables that
may be set with the 'set' function used in this script. SI-units (mks)
and degrees are used unless otherwise specified.

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

path = os.environ['FLICK_PATH']+"/Example/accurt_calls/sea_ice"
os.chdir(path)

if not os.path.exists('output'):
    os.makedirs('output')
    
class sea_ice:
    detector_height = -5 # [m]
    wl_low = 400e-9 # [m]
    wl_high = 600e-9 # [m]
    n_wl =  10
    wl_band_width = 10e-9 # [m]
    time_point_utc = '2025 08 01 12 00' # yyyy mm dd hh mm
    latitude = 85 # [degree]
    longitude = 5 # [degree]
    ocean_depth_grid = [0, 0.01, 0.01001, 100]
    
    def __init__(self, flick_radiation_object):
        self.f = flick_radiation_object
        self.f.set("snow_ice", 0.01)
        self.f.set("snow_radius", 1e-3)        
        self.f.set('ice_depths',2)
        self.f.set('ice_bubble_fraction',[0.01, 0.01])
        self.f.set('ice_brine_fraction',[0.01, 0.01])
        self.set_derived_parameters()

    def set_derived_parameters(self):
        self.f.set('bottom_depth',self.ocean_depth_grid[-1])
        self.f.set('concentration_relative_depths',self.absolute_to_relative(self.ocean_depth_grid))
        self.f.set('concentration_scaling_factors',np.ones(len(self.ocean_depth_grid)))
 
    def set(self,parameter_name,value):
        self.f.set(parameter_name,value)

    def to_W_per_m2_nm(self, radiation_values):
        return self.f.to_W_per_m2_nm(radiation_values)

    def to_mW_per_m2_nm_sr(self, radiation_values):
        return self.f.to_mW_per_m2_nm_sr(radiation_values)
    
    def absolute_to_relative(self, depths):
        b = self.f.get('bottom_depth')
        return np.array(depths)/b
    
    def wavelength(self):
        return flick.atmosphere_wavelengths('config_for_wavelengths',
                                            self.wl_low, self.wl_high, self.n_wl)
        
    def radiation(self):
        self.set_derived_parameters()
        return self.f.spectrum(self.wavelength(), self.wl_band_width,
                               self.time_point_utc+' 0', self.latitude, self.longitude)

def plane_downward_irradiance(height):
    si = sea_ice(flick.ocean_downward_plane_irradiance())
    si.set('detector_height', height)
    Ed = si.to_W_per_m2_nm(si.radiation())
    return Ed[:,0], Ed[:,1]

def surface_downward_scalar_irradiance():
    si = sea_ice(flick.ocean_downward_plane_irradiance())
    si.set('detector_type','scalar_irradiance')
    si.set('detector_height', 1)
    Eds = si.to_W_per_m2_nm(si.radiation())
    return Eds[:,0], Eds[:,1]

def nadir_radiance(height):
    si = sea_ice(flick.ocean_nadir_radiance())
    si.set('detector_height', height)
    Lu = si.to_mW_per_m2_nm_sr(si.radiation())
    return Lu[:,0], Lu[:,1]

def save(file_name, x, y):
    f = open(file_name,'w')
    if not f:
        raise IOError(file_name)
    for i in range(len(x)):
        f.write(f"{x[i]:#.{5}g}\t{y[i]:#.{4}g}\n")
    f.close()

def plot_irradiance():
    depth = 9
    x1, y1 = plane_downward_irradiance(height=-depth)
    x2, y2 = surface_downward_scalar_irradiance()
    plt.plot(x1,y1,label=f'plane at {depth:#.{3}g} m depth')
    plt.plot(x2,y2,label='scalar at surface')
    plt.legend()
    plt.grid()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel(r'Downward irradiance [W m$^{-2}$ nm$^{-1}$]')
    save('output/plane_irradiance.txt',x1,y1)
    save('output/surface_scalar_irradiance.txt',x2,y2)

def plot_radiance():
    depth = 9.37
    x, y = nadir_radiance(height=-9.37)
    plt.plot(x,y,label=f'at {depth:#.{3}g} m depth')
    plt.legend()
    plt.grid()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel(r'Nadir radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
    save('output/nadir_radiance.txt',x,y)

def plot_all():
    plt.figure(figsize=(9, 4))
    plt.subplot(1,2,1)
    plot_irradiance()
    plt.subplot(1,2,2)
    plot_radiance()
    plt.tight_layout()
    
if __name__ == "__main__":
    plot_all()
    plt.show()
    
