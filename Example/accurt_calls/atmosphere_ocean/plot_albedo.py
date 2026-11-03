"""
Plot ocean albedo data.
Run compute_albedo.py first.
"""
import numpy as np
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import os
import flick
flick.setup_script(__file__)


data = flick.table('output/computed_albedo.txt')
wls = data[:,0]
albedo = data[:,1]
fig, ax = plt.subplots()
ax.plot(wls,albedo)
ax.grid()
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'Ocean albedo')
fig.savefig('output/computed_albedo.pdf')

if __name__ == "__main__":
    plt.show()


