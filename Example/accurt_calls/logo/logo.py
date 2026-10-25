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
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Ellipse
from matplotlib.textpath import TextPath
import flick

from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

black_background = True
n_polar = 90
n_azimuth = 150
wl = 300e-9

def make_logo(recalculate=True):
    if recalculate:
        c = flick.accurt_config()
        c.set_streams(32)
        c.set("n_heights", 8)
        c.set("detector_height", -10)
        c.set("source_zenith_angle", 60)
        c.set("aerosol_od", 0.1)
        c.set("cloud_liquid", 0)
        c.set("nap_concentration", 0)
        c.set("chl_concentration", 0)
        c.set("cdom_440", 0)
        c.set("wavelengths", wl)
        c.set("detector_radiance_distribution_override", [n_polar, n_azimuth])
        r = flick.run(f"accurt {c.filename}")
        np.save("flick_tmp/logo", r)
    else:
        r = np.load("flick_tmp/logo.npy")

    theta = np.linspace(0, np.pi, n_polar)
    phi = np.linspace(-np.pi, np.pi, n_azimuth)
    r = (r / r.max()).transpose()
    t, p = np.meshgrid(theta+np.pi, phi)
    x = r * np.sin(t) * np.cos(p)
    y = r * np.sin(t) * np.sin(p)
    z = r * np.cos(t)

    if black_background:
        plt.style.use('dark_background')
    ls = LightSource(azdeg=0, altdeg=95)
    rgb = ls.shade(r, cmap=cm.hot, vert_exag=0, blend_mode='soft')
    fig = plt.figure(figsize=(9, 3.5))
    # Keep the radiance figure in its own area so it cannot overlap the name.
    ax = fig.add_axes([0.02, 0.05, 0.34, 0.9], projection='3d')
    ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
                    linewidth=0, antialiased=False, shade=False, alpha=1)

    font_size = 150
    font = FontProperties(family='Arial', weight='heavy')
    text_x, text_y = 0.3, 0.29
    fig.text(text_x, text_y, "flick", fontsize=font_size,
             fontproperties=font, color='white' if black_background else 'black')

    # Locate the centre of the i from the actual font metrics.  Drawing this in
    # figure coordinates keeps the red circle over the dot at every DPI.
    prefix_width = TextPath((0, 0), "Fl", size=font_size,
                            prop=font).get_extents().xmax
    i_width = TextPath((0, 0), "i", size=font_size,
                       prop=font).get_extents().xmax
    dot_x = text_x + (prefix_width + i_width / 2) / (72 * fig.get_figwidth())-0.04
    dot_y = text_y + 0.735 * font_size / (72 * fig.get_figheight())-0.045
    dot_diameter = 30 / 72
    fig.add_artist(Ellipse((dot_x, dot_y),
                           width=dot_diameter / fig.get_figwidth(),
                           height=dot_diameter / fig.get_figheight(),
                           transform=fig.transFigure, facecolor=[0.81, 0, 0],
                           edgecolor='none', zorder=10))
    ax.view_init(elev=0, azim=-150, roll=0)
    ax.set_aspect('equal')
    ax.set_axis_off()
    return fig

if __name__ == "__main__":
    fig = make_logo()
    if black_background:
        fig.savefig("logo_black.png",dpi=500)
        fig.savefig("logo_black.jpg",dpi=500)
    else:
        fig.savefig("logo_white.png",dpi=500)
        fig.savefig("logo_white.jpg",dpi=500)
    plt.show()
