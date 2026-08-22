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
import os
import sys

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), ".matplotlib_cache"))

import matplotlib.pyplot as plt
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

path = os.environ['FLICK_PATH']+"/Example/accurt_calls/sea_ice"
os.chdir(path)

if not os.path.exists('output'):
    os.makedirs('output')
    
class sea_ice:
    detector_height = 1 # [m]
    wl_low = 310e-9 # [m]
    wl_high = 1200e-9 # [m]
    n_wl =  80
    wl_band_width = 10e-9 # [m]
    time_point_utc = '2026 08 30 12 00' # yyyy mm dd hh mm
    latitude = 90 # [degree]
    longitude = 90 # [degree]
    ocean_depth_grid = [0, 0.5, 0.5001, 500] # [m]
    def __init__(self, flick_radiation_object, run_info):
        self.f = flick_radiation_object
        self.f.set('cloud_liquid', 1e-4)
        if run_info == 'snow_on_ice':
            self.f.set("snow_ice", 0.01)
            self.f.set("snow_radius", 1e-3)        
            self.f.set('ice_depths',2)
            self.f.set('ice_bubble_fraction',[0.005, 0.005])
            self.f.set('ice_brine_fraction',[0.02, 0.02])
        self.f.set('gases',['o3','o2','h2o'])
        self.f.set('gas_spectral_region','solar')
        self.f.set('cdom_440',0.1)        
        self.f.set('chl_concentration',0.1e-6)        
        self.f.set('nap_concentration',0.1e-3)        
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
        return flick.atmosphere_wavelengths('./flick_tmp/config',                                            self.wl_low, self.wl_high, self.n_wl)
        
    def radiation(self):
        self.set_derived_parameters()
        return self.f.spectrum(self.wavelength(), self.wl_band_width,
                               self.time_point_utc+' 0', self.latitude, self.longitude)

def downward_plane_irradiance(height, run_info):
    si = sea_ice(flick.ocean_downward_plane_irradiance(), run_info)
    si.set('detector_type','plane_irradiance')
    si.set('detector_orientation','up')
    si.set('detector_height', height)
    Ed = si.to_W_per_m2_nm(si.radiation())
    return Ed[:,0], Ed[:,1]

def downward_scalar_irradiance(height, run_info):
    si = sea_ice(flick.ocean_downward_plane_irradiance(), run_info)
    si.set('detector_type','scalar_irradiance')
    si.set('detector_orientation','up')
    si.set('detector_height', height)
    Eds = si.to_W_per_m2_nm(si.radiation())
    return Eds[:,0], Eds[:,1]

def nadir_radiance(height, run_info):
    si = sea_ice(flick.ocean_nadir_radiance(), run_info)
    si.set('detector_height', height)
    Lu = si.to_mW_per_m2_nm_sr(si.radiation())
    return Lu[:,0], Lu[:,1]

def zenith_radiance(height, run_info):
    si = sea_ice(flick.ocean_nadir_radiance(), run_info)
    si.set('detector_height', height)
    si.set('detector_orientation','up')
    Ld = si.to_mW_per_m2_nm_sr(si.radiation())
    return Ld[:,0], Ld[:,1]

def save(file_name, x, y):
    f = open(file_name,'w')
    if not f:
        raise IOError(file_name)
    for i in range(len(x)):
        f.write(f"{x[i]:#.{5}g}\t{y[i]:#.{4}g}\n")
    f.close()

def plot_and_save_all(run_info, height):
    x1, y1 = downward_plane_irradiance(height, run_info)
    x2, y2 = downward_scalar_irradiance(height, run_info)
    x3, y3 = zenith_radiance(height, run_info)
    y3 *= 1e-3 # From mW to W
    run_info_file = run_info.replace(' ', '_')
    title = (f"{run_info}; Position: {sea_ice.latitude:g} deg N, {sea_ice.longitude:g} deg E; "
             f"Time: {sea_ice.time_point_utc} UTC")
    output_file_base = (f"output/{run_info_file}_lat_{sea_ice.latitude:g}_lon_{sea_ice.longitude:g}_"
                        f"{sea_ice.time_point_utc.replace(' ', '_')}")
    plot_file_name = f"{output_file_base}.png"
    
    plt.figure(figsize=(9, 4))
    plt.suptitle(title, y=0.94)
    plt.subplot(1,2,1)
    plt.plot(x1,y1*2,label=r'plane $\times$ 2')
    plt.plot(x2,y2,label='scalar')
    plt.ylabel(r'Irradiance [W m$^{-2}$ nm$^{-1}$]')
    plt.xlabel('Wavelength [nm]')
    plt.legend()
    plt.grid()
    plt.subplot(1,2,2)
    plt.plot(x3,y3*2*np.pi,label=r'zenith radiance $\times 2\pi$')
    plt.legend()
    plt.grid()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel(r'Zenith radiance [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.savefig(plot_file_name, dpi=300)
    
    save(f"{output_file_base}_plane_irradiance.txt",x1,y1)
    save(f"{output_file_base}_scalar_irradiance.txt",x2,y2)
    save(f"{output_file_base}_radiance.txt",x3,y3)

    
if __name__ == "__main__":
    height = 1
    for run_info in ['no_ice', 'snow_on_ice']:
        plot_and_save_all(run_info, height)
    plt.show()
    
