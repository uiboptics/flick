"""
Compute top-of-atmosphere (TOA) radiances at 21 spectral bands as measured
by the OLCI radiometer onboard the Sentinel-3 satellite.

See also README.txt

Unless otherwise specified, SI base units are used, with angles given in
degrees.

"""

"""This script is meant for HTEK301. Teaming up with your favorite AI
is encorraged to improve it as much as possible within the short
timefame of the optics module.
"""

import numpy as np
from types import SimpleNamespace
import flick

flick.setup_script(__file__)

from make_config import c

"""Example Sentinel-3 OLCI data from Hardangerfjord (Hamre et al., 2026)."""


"""UTC date and time followed by geographical position: 
   [year, month, day, hour, minute, second, latitude, longitude].
"""
time_pos = [2022, 5, 20, 10, 51, 0, 59.87, 5.68]


"""OLCI observation angles [degrees]:
   [observation zenith angle (OZA), observation azimuth angle (OAA)].
"""
obs_ang = [14.0, -68.0]


"""Number of wavelengths sampled within each OLCI band.
   May be increased until the calculated band-averaged radiances have
   converged sufficiently for the intended application. Note, however,
   that increasing this value can substantially increase the
   computational cost.
"""
n_wls = 1


def sentinel_wavelength():
    """Nominal central wavelengths of the 21 OLCI bands [nm]."""
    return [400.0, 412.5, 442.5, 490.0, 510.0, 560.0, 620.0,
            665.0, 673.8, 681.3, 708.8, 753.8, 761.3, 764.4,
            767.5, 778.8, 865.0, 885.0, 900.0, 940.0, 1020.0]


def sentinel_radiance():
    """Example OLCI TOA radiance spectrum for Hardangerfjord.
       Units are mW m-2 nm-1 sr-1. Data from Hamre et al. (2026).
    """
    return [75.76, 81.17, 82.73, 81.42, 79.10, 71.51, 31.02,
            21.75, 20.23, 19.31, 13.15, 6.914, 2.669, 3.919,
            5.910, 5.929, 3.430, 3.030, 2.267, 1.115, 1.985]


def sentinel_deviation():
    """Geometric standard deviation within a 3-by-3 pixel area.
       The returned values are dimensionless multiplicative factors.
    """
    return [1.002, 1.002, 1.002, 1.003, 1.003, 1.003, 1.003, 1.004,
            1.004, 1.005, 1.004, 1.008, 1.017, 1.008, 1.011, 1.005,
            1.009, 1.008, 1.013, 1.028, 1.040]


def set_fixed_parameter():
    """See flick_tmp/config for complete parameter documentation"""
    c.set('cloud_liquid', 0)
    c.set('water_vapor',19)
    c.set('pressure',101e3)
    c.set('temperature',293)
    c.set('ozone',0.003)
    c.set('toa_solar_multiplication',time_pos+[0])
    c.set('detector_sentinel3_override',[n_wls,obs_ang[0],obs_ang[1]])
    c.set_streams(12)

class parameter:
    """Use SI units, for example kg m-3 for concentrations."""

    name = ["aerosol_od","cdom_440","chl_concentration","nap_concentration"]
    lower_limit = np.array([0.02, 0.02, 0.1e-6, 0.2e-3])
    upper_limit = np.array([1.0, 1.0, 10e-6, 22e-3])
    value = (upper_limit + lower_limit) / 2

    rng = np.random.default_rng()

    @classmethod
    def set_random(cls):
        cls.set_value(cls.rng.uniform(cls.lower_limit, cls.upper_limit))
        return cls

    @classmethod
    def clip(cls, n, value):
        return np.clip(value,cls.lower_limit[n],cls.upper_limit[n])

    @classmethod
    def set_value(cls, new_value):
        new_value = np.asarray(new_value, dtype=float)
        cls.value = np.clip(new_value,cls.lower_limit,cls.upper_limit)
        return cls
            
    
def modeled_radiance(parameter):
    """Calculate the modeled Sentinel-3 OLCI radiance spectrum.

    Possible project extensions:

    Can the cost-function minimum be improved by converting some of the
    fixed parameters into fitted parameters? For example water
    vapor and ozone? Increasing the number of fitted parameters may,
    however, make the optimization slower and could make the solution
    less well constrained.

    This function is computationally expensive because it runs the
    full radiative-transfer Flick-AccuRT model. Training a
    neural-network surrogate on a large collection of spectra
    calculated using randomly sampled parameters could be investigated
    as a possible acceleration method. However, developing and
    validating such a surrogate may be far beyond the scope of the
    HTEK301 project, even with AI support :)
    """
    for i in range(len(parameter.name)):
        c.set(parameter.name[i], parameter.value[i])
    return flick.run("accurt flick_tmp/config")


class cost_function:
    parameter = parameter
    wavelengths = sentinel_wavelength()
    measured = np.asarray(sentinel_radiance())
    modeled = None
    value = np.inf

    @classmethod
    def simple_least_square(cls, parameter):
        """Calculate the mean squared difference between the spectra.

        Possible project extensions:

        This is a simple least-squares cost function. Could the fitting be
        improved by using relative residuals, a logarithmic cost function,
        a spectral-angle metric, or correlated measurement uncertainties?
        
        Another possibility is a wavelength-integrated squared residual.
        Would that provide a meaningful way to account for the varying
        spectral separation between the OLCI bands?
        """
        cls.parameter = parameter
        cls.modeled = np.asarray(modeled_radiance(cls.parameter))[:,1]
        print(cls.modeled)
        cls.value = np.mean((cls.modeled - cls.measured) ** 2)
        return cls

    @classmethod
    def copy(cls):
        return SimpleNamespace(
            parameter=SimpleNamespace(
                name=cls.parameter.name.copy(),
                value=cls.parameter.value.copy(),
            ),
            modeled=cls.modeled.copy(),
            value=float(cls.value),
        )

    
class optimizer:
    """Optimize the model parameters by minimizing the cost function."""

    best = None

    @classmethod
    def store_if_best(cls, candidate):
        if cls.best is None or candidate.value < cls.best.value:
            cls.best = candidate.copy()

    @classmethod
    def lucky_guess(cls, n_iterations=5):
        """Perform a simple random-search inversion.

        This is probably the simplest and slowest optimization algorithm.
        It repeatedly guesses parameter values and checks whether the
        resulting cost is smaller than the best cost found so far.

        Possible project extensions:

        Could the optimization be performed more efficiently using a more
        sophisticated method? Possibilities include steepest descent
        with finite-difference derivatives, Powell's method, L-BFGS-B,
        or another suitable optimization method.
        """
        for _ in range(n_iterations):
            candidate = parameter.set_random()
            cls.store_if_best(cost_function.simple_least_square(candidate))
        return cls.best
    

if __name__ == "__main__":
    set_fixed_parameter()
    best = optimizer.lucky_guess(5)
    
#output = np.column_stack((bubble_volume_fractions, albedo))
#flick.save_two_columns(output, "albedo_versus_volume_fraction.txt")

