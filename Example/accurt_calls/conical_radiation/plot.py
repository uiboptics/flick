"""
This script only works if `compute.py` has been run first.

See also README.txt
"""
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, ScalarFormatter
import numpy as np
from pathlib import Path

fig, ax = plt.subplots()

for output_file in sorted(Path("output").glob("*.txt")):
    acceptance_halfangle, conical_irradiance = np.loadtxt(output_file, unpack=True)
    sold_angle = 2*np.pi*(1-np.cos(acceptance_halfangle/180*np.pi))
    ax.plot(
        acceptance_halfangle,
        conical_irradiance/sold_angle*1e3,
        label=output_file.stem,
    )

ax.grid()
ax.set_yscale("log")
ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.set_xlim(-1, 91)
ax.set_title('Solid-angle-averaged conical irradiance')
ax.set_ylabel(r'Irradiance per solid angle [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
ax.set_xlabel(r'Acceptance half-angle [deg]')
ax.legend(title='orientation [deg]')
fig.savefig('output/conical_irradiance.pdf')

if __name__ == "__main__":
    plt.show()


