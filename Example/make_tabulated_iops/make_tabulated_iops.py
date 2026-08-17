"""
Example of how to make water bubbles input data to the tabulated_iops
material.  Customize by modifying parameters at the bottom of this
script.

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.environ['FLICK_PATH']+"/python_script")
import flick

def refractive_index_host(wl):
    T = 293
    S = 33
    sub_command = f"iop refractive_index {wl} {wl} 1"
    real = flick.run(f"{sub_command} pure_water {T} {S}")[1]
    sub_command = f"iop absorption_length {wl} {wl} 1"
    abs_length = flick.run(f"{sub_command} pure_water {T} {S}")
    abs_coef = 1/abs_length[1]
    imag = abs_coef * wl / (4 * np.pi);
    return f'{real}+{imag}i'

def refractive_index_sphere():
    return '1'

def mie_sub_command(wl):
    h = refractive_index_host(wl)
    s = refractive_index_sphere()
    return f'mie {wl} {h} {s} {r_median} {sigma} {percent_accuracy}'

def scattering_elements(row,col,wl):
    c = f'{mie_sub_command(wl)} scattering_matrix_element {row} {col} {n_angles}'
    return flick.run(c)

def normalized_scattering_elements(row,col,wl):
    return scattering_elements(row,col,wl)/scattering_cross_section(wl)

def absorption_coefficient_per_vf(wl):
    return flick.run(f'{mie_sub_command(wl)} absorption_coefficient_per_vf')

def absorption_cross_section(wl):
    return flick.run(f'{mie_sub_command(wl)} absorption_cross_section')

def scattering_coefficient_per_vf(wl):
    return flick.run(f'{mie_sub_command(wl)} scattering_coefficient_per_vf')

def scattering_cross_section(wl):
    return flick.run(f'{mie_sub_command(wl)} scattering_cross_section')

def stream_normalized_scattering_elements(wls,row,col):
    s = """/* 
Normalized scattering matrix element for the material described
in ab.txt. Element 1,1 is the phase function with 4 pi integral equal
to one. This file is for the element given by the file name,
e.g. 's_11' indicates element 1,1.

First data row: wavelengths [nm]

Second data row: scattering angles [degrees]

Then: normalized scattering values [1/sr] for each wavelength (rows)
and each angle (columns).  
*/

"""
    wls_nm = wls*1e9
    s += " ".join(f"{x:#.4g}" for x in wls_nm)+'\n\n'
    for i in range(len(wls)):
        se = normalized_scattering_elements(row,col,wls[i])
        if i == 0:
            angles_deg = se[:,0]/np.pi*180;
            s += " ".join(f"{x:#.4g}" for x in angles_deg)+'\n\n'            
        s += " ".join(f"{x:#.2e}" for x in se[:,1])+'\n'    
    return s

def write_all_normalized_scattering_elements(wls, include_elements):
    for el in include_elements:
        row = el[0]
        col = el[1]
        with open(f's_{row+1}{col+1}.txt', 'w') as f1:
            print(stream_normalized_scattering_elements(wls,row,col),file=f1)

def write_include_elements(include_elements):
    header = """/*         
Only the following elements will be loaded into the normalized
scattering matrix. Other elements will be zero.

First column: scattering matrix row number, counting from 1
Second column: scattering matrix column number, counting form 1
*/
"""
    with open(f'include_elements.txt', 'w') as f:
        print(header,file=f)
        for el in include_elements:
            row = el[0]
            col = el[1]
            print(f'{row+1} {col+1}',file=f)

def write_ab(wls):
    header = """/*
Spherical bubbles in sea water with S = 33 PSU and T = 293 K.

Frist column: wavelengths [nm]
Second column: absorption coefficient per volume fraction [1/m]
Third column: scattering coefficient per volume fraction [1/m]
*/
"""
    with open(f'ab.txt', 'w') as f:
        print(header,file=f)
        for wl in wls:
            a = absorption_coefficient_per_vf(wl)
            b = scattering_coefficient_per_vf(wl)
            wl_nm = wl*1e9
            print(f'{wl_nm:#.4g}\t{a:#.2e}\t{b:#.2e}',file=f)

if __name__ == "__main__":
    n_angles = 5
    include_elements = [[0,0], [1,1]]
    r_median = 1e-7
    sigma = 0.1
    percent_accuracy = 1
    wls = np.linspace(300e-9, 700e-9, 4)

    write_ab(wls)
    write_all_normalized_scattering_elements(wls, include_elements)
    write_include_elements(include_elements)

    
