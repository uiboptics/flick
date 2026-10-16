"""
This script only works if `compute.py` has been run first.

See also README.txt
"""
import matplotlib.pyplot as plt
import numpy as np
from make_config import sys, os, flick

fig, ax = plt.subplots()

tilted_path = 'output/albedo_horizontal.txt'
if os.path.exists(tilted_path):
    a_tilted = np.loadtxt(tilted_path, ndmin=2)
    ax.plot(a_tilted[:,0], a_tilted[:,1],'-o', label='Horizontal')
    ax.legend()

tilted_path = 'output/albedo_tilted.txt'
if os.path.exists(tilted_path):
    a_tilted = np.loadtxt(tilted_path, ndmin=2)
    ax.plot(a_tilted[:,0], a_tilted[:,1],'-o', label='Tilted')
    ax.legend()

ax.grid()
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'Snow albedo')
ax.set_ylim([0.5,1.1])
fig.savefig('output/albedo.pdf')

if __name__ == "__main__":
    plt.show()


