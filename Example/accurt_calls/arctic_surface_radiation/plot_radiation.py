"""
This script only works if `compute_radiation.py` has been run first.

See also README.txt
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

sys.path.append(os.environ['FLICK_PATH'] + "/python_script")
import flick

os.chdir(Path(__file__).resolve().parent)

fig, ax = plt.subplots()

directory = Path("output")

for file in sorted(directory.glob("*.txt")):
    if file.is_file():
        x, y = np.loadtxt(file, unpack=True)
        ax.plot(x, y, label=file.stem)

ax.legend()
ax.grid()
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'Radiation [W m$^{-2}$ nm$^{-1}$] or [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')

fig.savefig('output/computed_radiation.pdf')

if __name__ == "__main__":
    plt.show()
    

    

