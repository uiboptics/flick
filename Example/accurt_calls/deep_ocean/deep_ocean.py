"""
See the flick_tmp/config file, which will be
generated after the first run, for documentation on all variables that
may be set with the set function used in this logo script. SI-units and
degrees are used unless otherwise specified.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
from matplotlib import cm
from matplotlib.colors import LightSource
import flick

def xyz_surface(radiance, n_polar, n_azimuth):
    r = radiance
    theta = np.linspace(0, np.pi, n_polar)
    phi = np.linspace(-np.pi, np.pi, n_azimuth)
    rt = r.transpose() 
    t, p = np.meshgrid(theta, phi)
    x = rt * np.sin(t) * np.cos(p)
    y = rt * np.sin(t) * np.sin(p)
    z = rt * np.cos(t)
    return x,y,z

class deep_ocean:
    def __init__(self, depth, sza):
        self.depth = depth
        self.F_0 = 1
        self.config = flick.basic_radiation()
        self.config._generate_config('toa_reflectance', './flick_tmp')
        self.config.set('detector_height', -depth)
        self.config.set_n_angles(100)
        self.config.set("n_heights", 8)
        self.config.set("pressure",0.001)
        self.config.set("source_zenith_angle", sza)
        self.config.set("aerosol_od", 0)
        self.config.set("cloud_liquid", 0)
        self.config.set("nap_concentration", 0)
        self.config.set("chl_concentration", 0)
        self.config.set("cdom_440", 0.004)
        self.config.set("bottom_depth", depth)
        self.config.set('detector_wavelengths', [430e-9])

    def run(self):
        return flick.run('accurt ./flick_tmp/config')
        
    def normalized_irradiance_down(self):
        self.config.set('detector_type','plane_irradiance')
        self.config.set('detector_orientation','up')
        return self.run()[1]/self.F_0
    
    def normalized_irradiance_up(self):
        self.config.set('detector_type','plane_irradiance')
        self.config.set('detector_orientation','down')
        return self.run()[1]/self.F_0 

    def normalized_radiance(self, n_polar, n_azimuth):
        self.config.set("detector_radiance_distribution_override",
                        [n_polar, n_azimuth])
        return self.run()/self.F_0

    def plot_normalized_radiance(self, n_polar, n_azimuth):
        r = self.normalized_radiance(n_polar, n_azimuth)
        x,y,z = xyz_surface(r, n_polar, n_azimuth)
        ls = LightSource(azdeg=0, altdeg=55)
        rgb = ls.shade(np.transpose(r), cmap=cm.hot, vert_exag=0, blend_mode='soft')
        fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
        ax.set_aspect('equal')
        surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
                       linewidth=1, antialiased=True, shade=False,alpha=1)
        ax.set_title(r'Normalized radiance [sr$^{-1}$]'+
                     f' at {self.depth} m depth')
        
if __name__ == "__main__":
    do = deep_ocean(2000,sza=0)
    print(f'Ocean transmittance: {do.normalized_irradiance_down():#.4g}')
    do.plot_normalized_radiance(150,50)
    plt.show()
