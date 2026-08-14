""" Run to make albedo data for an ocean with bubbles, then
plot with plot_albedo.py

"""
import numpy as np
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

path = os.environ['FLICK_PATH']+"/Example/accurt_calls/atmosphere_ocean"
os.chdir(path)

f = flick.relative_radiation()

f.set('detector_type','plane_irradiance')
f.set('detector_orientation','down')
f.set('detector_height', 1)
f.set('reference_detector_height', 1)

f.set('aerosol_od', 0.1)
f.set('cloud_liquid', 0)
f.set('temperature', 273+15)
f.set('salinity', 33)
f.set('water_vapor',20)
f.set('mp_names', 'input/ECOSENS_HF22_D1')
f.set('mp_concentrations', 0.01)
f.set('mcdom_names', 'input/ECOSENS_HF22_D1')
f.set('mcdom_scaling_factors', 1)
f.set('detector_height', 1)
f.set('bottom_depth', 100)

f.set('bubble_volume_fraction',0.0001)
f.set('bubble_calculator','parameterized_mie')
#f.set('bubble_calculator','full_mie')
f.set('bubble_radius',1e-6)
f.set('bubble_sigma',0.0)

f.set('concentration_relative_depths', [0,0.01,0.01001,1])
f.set('concentration_scaling_factors', [1,1,1,1])
f.set('concentration_exception_names', 'bubbles')
f.set('concentration_exception_scaling_factors', [1,1,0,0])
       
wl = np.linspace(310e-9, 700e-9, 3)
albedo = f.spectrum(wl, source_zenith_angle = 20)
albedo[:,0] *= 1e9

if not os.path.exists('output'):
    os.makedirs('output')
np.savetxt('output/computed_albedo.txt', albedo, fmt=['%6.2f ','%8.3e'])


