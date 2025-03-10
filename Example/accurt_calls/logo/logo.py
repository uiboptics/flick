"""
Make the Flick logo. See the flick_tmp/config file, which will be
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

recalculate = True
black_background = True
n_polar = 50
n_azimuth = 150
wl = 300e-9
f = flick.radiance_distribution(n_polar,n_azimuth)

f.set_n_angles(100)
f.set("n_heights", 8)
f.set("detector_height", -1)
f.set("source_zenith_angle", 0)
f.set("aerosol_od", 0)
f.set("cloud_liquid", 0)
f.set("nap_concentration", 0)
f.set("chl_concentration", 0)
f.set("cdom_440", 0)

if (recalculate):
    r = f.values(wl,source_zenith_angle=0)
    np.save("flick_tmp/logo",r)
else:
    r = np.load("flick_tmp/logo.npy")

r = r/r.max()    
x,y,z = f.radiance_surface(r)

if black_background:
    plt.style.use('dark_background')
ls = LightSource(azdeg=0, altdeg=55)
rgb = ls.shade(np.transpose(r), cmap=cm.hot, vert_exag=0, blend_mode='soft')
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
                       linewidth=0, antialiased=False, shade=False,alpha=1)

font = {'fontname':'Arial'}
ax.text(0,0.77,-0.35,"flick",fontsize=180,fontweight='heavy',**font)

ax.plot(0,1.93,0.67,'o',color='black',markerfacecolor=[0.81,0,0],
        markeredgewidth=0,
        markersize=35,zorder=10)
ax.view_init(elev=0, azim=0, roll=0)

ax.set_aspect('equal')
fig.set_size_inches(9,3.5,forward=True)
plt.subplots_adjust(left=-0.3,
                    bottom=-0.3,
                    right=0.9,
                    top=1.28)
plt.axis('off')

if __name__ == "__main__":
    plt.show()
    if black_background:
        fig.savefig("logo_black.png",dpi=500)
        fig.savefig("logo_black.jpg",dpi=500)
    else:
        fig.savefig("logo_white.png",dpi=500)
        fig.savefig("logo_white.jpg",dpi=500)
