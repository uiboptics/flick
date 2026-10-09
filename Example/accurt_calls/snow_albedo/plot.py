"""
This script only works if `compute.py` has been run first.

See also README.txt
"""
import matplotlib.pyplot as plt
from make_config import sys, os, flick

name = 'albedo'
a = flick.table(f'output/{name}.txt')
fig, ax = plt.subplots()
ax.plot(a[:,0],a[:,1])
ax.grid()
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'Snow albedo')
ax.set_ylim([0.6,1.02])
fig.savefig('output/{name}.pdf')

if __name__ == "__main__":
    plt.show()


