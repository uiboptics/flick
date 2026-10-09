"""
This script only works if `compute.py` has been run first.

See also README.txt
"""
import matplotlib.pyplot as plt
from make_config import sys, os, flick

name = 'albedo'
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


