"""
This script only works if `compute_albedo.py` has been run first.

See also README.txt
"""
import numpy as np
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

name = 'computed_albedo'
data = flick.table(f'output/{name}.txt')
wls = data[:,0]
albedo = data[:,1]
fig, ax = plt.subplots()
ax.plot(wls,albedo)
ax.grid()
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'Snow albedo')
ax.set_ylim([0.6,1.02])
fig.savefig('output/{name}.pdf')

if __name__ == "__main__":
    plt.show()


