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
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
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
    #time_point_utc = '2026 08 30 12 00' # yyyy mm dd hh mm
    time_point_utc = '2026 08 29 08 00' # yyyy mm dd hh mm
    #latitude = 90 # [degree]
    #longitude = 90 # [degree]
    #latitude = 86 # [degree]
    #longitude = 104 # [degree]
    latitude = 90 # [degree]
    longitude = 0 # [degree]
    
    ocean_depth_grid = [0, 0.5, 0.5001, 500] # [m]

    @staticmethod
    def parse_run_info(run_info):
        run_conditions = [condition.strip() for condition in run_info.split(',')]
        if len(run_conditions) != 2:
            raise ValueError('run_info must be "<clear/cloudy>, <open water/snow and ice>"')
        return run_conditions

    def __init__(self, flick_radiation_object, run_info):
        self.f = flick_radiation_object
        cloud_condition, ice_condition = self.parse_run_info(run_info)

        if cloud_condition == 'clear':
            print('clear...')
            self.f.set('cloud_liquid', 0)
        elif cloud_condition == 'cloudy':
            print('cloudy...')
            self.f.set('cloud_liquid', 1e-4) 
        else:
            raise ValueError('wrong clear/cloudy condition') 
        if ice_condition == 'open water':
            print('open water...')
            self.f.set("snow_ice", 0)
            self.f.set('ice_depths',0)
        elif ice_condition == 'snow and ice':
            print('snow and ice...')
            self.f.set("snow_ice", 0.01)
            self.f.set("snow_radius", 1e-3)        
            self.f.set('ice_depths',2)
            self.f.set('ice_bubble_fraction',[0.005, 0.005])
            self.f.set('ice_brine_fraction',[0.02, 0.02])
        else:
            raise ValueError('wrong snow, ice, open water condition') 
        self.f.set('gases',['o3','o2','h2o'])
        self.f.set('gas_spectral_region','solar')
        self.f.set('cdom_440',0.1)        
        self.f.set('chl_concentration',0.1e-6)        
        self.f.set('nap_concentration',0.1e-3)        
        self.set_derived_parameters()

    def normalized_radiance_distribution(self, wl):
        r = self.f.values(wl,source_zenith_angle=78)
        r_max = np.nanmax(r)
        if r_max > 0:
            r = r/r_max
        r = np.nan_to_num(r)
        polar_axis = np.argmin(r.shape)
        if r.shape[polar_axis] > 2:
            if polar_axis == 0:
                r[0, :] = np.mean(r[1, :])
                r[-1, :] = np.mean(r[-2, :])
            else:
                r[:, 0] = np.mean(r[:, 1])
                r[:, -1] = np.mean(r[:, -2])
        return r

    def plot_radiance_distribution(self, wl, ax=None):
        r = self.normalized_radiance_distribution(wl)
        x,y,z = self.f.radiance_surface(r)
        color_values = r
        if color_values.shape != x.shape and color_values.T.shape == x.shape:
            color_values = color_values.T
        color_values = 0.2 + 0.8*color_values
        colors = cm.inferno(color_values)
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 6), subplot_kw=dict(projection='3d'))
            fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=colors,
                        linewidth=0, antialiased=True, shade=False)
        axis_radius = 1.08*np.nanmax(np.abs([x, y, z]))
        ax.set_xlim(-axis_radius, axis_radius)
        ax.set_ylim(-axis_radius, axis_radius)
        ax.set_zlim(-axis_radius, axis_radius)
        ax.set_box_aspect((1, 1, 1), zoom=1.25)
        ax.view_init(elev=25, azim=-45)
        ax.set_xlabel('x', color='black', labelpad=-2)
        ax.set_ylabel('y', color='black', labelpad=-2)
        ax.set_zlabel('z', color='black', labelpad=-2)
        ax.set_title((f"Radiance distribution at {wl*1e9:g} nm\n"
                      r'[W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'), color='black', pad=2)
        ax.tick_params(colors='black', labelsize=7, pad=-2)
        ax.xaxis.label.set_color('black')
        ax.yaxis.label.set_color('black')
        ax.zaxis.label.set_color('black')
        ax.xaxis.pane.set_facecolor((1, 1, 1, 0))
        ax.yaxis.pane.set_facecolor((1, 1, 1, 0))
        ax.zaxis.pane.set_facecolor((1, 1, 1, 0))
        ax.xaxis.pane.set_edgecolor((0, 0, 0, 0.25))
        ax.yaxis.pane.set_edgecolor((0, 0, 0, 0.25))
        ax.zaxis.pane.set_edgecolor((0, 0, 0, 0.25))
        ax.xaxis._axinfo['grid']['color'] = (0, 0, 0, 0.12)
        ax.yaxis._axinfo['grid']['color'] = (0, 0, 0, 0.12)
        ax.zaxis._axinfo['grid']['color'] = (0, 0, 0, 0.12)

    def plot_radiance_profiles(self, wl, ax=None):
        r = self.normalized_radiance_distribution(wl)
        polar_angles = np.degrees(self.f.polar_angles())
        azimuth_angles = np.degrees(self.f.azimuth_angles())

        if r.shape == (len(azimuth_angles), len(polar_angles)):
            r = r.T

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 4))
            fig.patch.set_facecolor('white')

        norm = Normalize(vmin=np.min(azimuth_angles), vmax=np.max(azimuth_angles))
        cmap = cm.viridis
        for azimuth_angle, profile in zip(azimuth_angles, r.T):
            ax.plot(polar_angles, profile, color=cmap(norm(azimuth_angle)),
                    linewidth=0.65, alpha=0.55)

        ax.set_xlabel('Polar angle [degree]')
        ax.set_ylabel('r')
        ax.set_title(f"Normalized radiance radius at {wl*1e9:g} nm")
        ax.set_xlim(polar_angles[0], polar_angles[-1])
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.35)
        colorbar = ax.figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                                      pad=0.02, fraction=0.05)
        colorbar.set_label('Azimuth angle [degree]')
        
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


def plot_radiance_distribution(height, run_info):
    n_polar = 50
    n_azimuth = 150
    wl = 500e-9
    si = sea_ice(flick.radiance_distribution(n_polar,n_azimuth), run_info)
    si.set('detector_height', height)
    si.plot_radiance_distribution(wl)
    
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
    run_info_file = run_info.replace(',', '')
    title = (f"{run_info}; Position: {sea_ice.latitude:g} deg N, {sea_ice.longitude:g} deg E; "
             f"Time: {sea_ice.time_point_utc} UTC")
    output_file_base = (f"output/{run_info_file}_lat_{sea_ice.latitude:g}_lon_{sea_ice.longitude:g}_"
                        f"{sea_ice.time_point_utc.replace(' ', '_')}")
    plot_file_name = f"{output_file_base}.png"
    
    fig = plt.figure(figsize=(10.5, 8.2))
    fig.suptitle(title, y=0.97)
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0], projection='3d')
    ax4 = fig.add_subplot(gs[1, 1])

    ax1.plot(x1,y1*2,label=r'plane $\times$ 2')
    ax1.plot(x2,y2,label='scalar')
    ax1.set_ylabel(r'Irradiance [W m$^{-2}$ nm$^{-1}$]')
    ax1.set_xlabel('Wavelength [nm]')
    ax1.legend()
    ax1.grid()

    ax2.plot(x3,y3*2*np.pi,label=r'zenith radiance $\times 2\pi$')
    ax2.legend()
    ax2.grid()
    ax2.set_xlabel('Wavelength [nm]')
    ax2.set_ylabel(r'Zenith radiance [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')

    n_polar = 50
    n_azimuth = 150
    wl = 500e-9
    si = sea_ice(flick.radiance_distribution(n_polar,n_azimuth), run_info)
    si.set('detector_height', height)
    si.plot_radiance_distribution(wl, ax3)
    si.plot_radiance_profiles(wl, ax4)

    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.08, top=0.91, wspace=0.28, hspace=0.34)
    plt.savefig(plot_file_name, dpi=300, bbox_inches='tight', pad_inches=0.02)
    
    save(f"{output_file_base}_plane_irradiance.txt",x1,y1)
    save(f"{output_file_base}_scalar_irradiance.txt",x2,y2)
    save(f"{output_file_base}_radiance.txt",x3,y3)

    
if __name__ == "__main__":
    height = 1
    ri = ['cloudy, snow and ice', 'clear, snow and ice', 'cloudy, open water', 'clear, open water']
    #ri = ['clear, snow and ice']
    #ri = ['clear, open water']
    for run_info in ri:
        print('plot_and_save...')
        plot_and_save_all(run_info, height)
    plt.show()
    
