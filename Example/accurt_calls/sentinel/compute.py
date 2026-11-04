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
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from types import SimpleNamespace
from time import perf_counter
import textwrap
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
    upper_limit = np.array([1.0, 1.0, 15e-6, 25e-3])
    value = (upper_limit + lower_limit) / 2
    machup = np.array([0.14, 0.14, 6.8e-6, 21e-3-6.8e-6/0.01])
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
        cls.value = np.mean((cls.modeled - cls.measured) ** 2)
        return cls

    @classmethod
    def my_super_version(cls, parameter):
        """ better method could be implement here"""
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

    
class optimization_plot:
    """Display the live spectrum and summary for an optimization run."""

    figure = None
    spectrum_axis = None
    summary_axis = None

    @staticmethod
    def wrap_table_text(value, width):
        """Wrap table labels at spaces and underscores."""
        wrappable = str(value).replace("_", "_ ")
        lines = textwrap.wrap(
            wrappable,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
        return "\n".join(lines).replace("_ ", "_")

    @classmethod
    def format_spectrum_axis(cls):
        """Show every OLCI band centre and both major and minor grids."""
        wavelengths = cost_function.wavelengths
        labels = [f"{wavelength:g}" for wavelength in wavelengths]

        cls.spectrum_axis.set_xticks(wavelengths, labels=labels)
        cls.spectrum_axis.tick_params(axis="x", labelrotation=90, labelsize=8)
        cls.spectrum_axis.xaxis.set_minor_locator(
            LogLocator(base=10, subs=np.arange(2, 10) * 0.1)
        )
        cls.spectrum_axis.yaxis.set_minor_locator(
            LogLocator(base=10, subs=np.arange(2, 10) * 0.1)
        )
        cls.spectrum_axis.yaxis.set_major_formatter(
            FuncFormatter(
                lambda value, _: np.format_float_positional(value, trim="-")
            )
        )
        cls.spectrum_axis.xaxis.set_minor_formatter(NullFormatter())
        cls.spectrum_axis.yaxis.set_minor_formatter(NullFormatter())
        cls.spectrum_axis.margins(x=0.02)
        cls.spectrum_axis.grid(True, which="major", alpha=0.35)
        cls.spectrum_axis.grid(
            True, which="minor", linestyle=":", linewidth=0.6, alpha=0.3
        )

    @classmethod
    def initialize(cls):
        """Create the live spectrum and optimization-summary display."""
        if cls.figure is not None and plt.fignum_exists(cls.figure.number):
            plt.close(cls.figure)
        plt.ion()
        cls.figure, (cls.spectrum_axis, cls.summary_axis) = plt.subplots(
            1,
            2,
            figsize=(12, 5),
            gridspec_kw={"width_ratios": (3, 2)},
        )
        cls.figure.canvas.manager.set_window_title("Sentinel spectrum optimization")

        cls.spectrum_axis.plot(
            cost_function.wavelengths,
            cost_function.measured,
            "o-",
            label="Measured",
        )
        cls.spectrum_axis.set(
            xlabel="Wavelength [nm]",
            ylabel=r"Spectral radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]",
            title="Measured and best spectra",
            xscale="log",
            yscale="log",
        )
        cls.format_spectrum_axis()
        cls.spectrum_axis.legend()

        cls.summary_axis.set_title("Current best parameters")
        cls.summary_axis.axis("off")
        cls.summary_axis.text(
            0.5,
            0.5,
            "Waiting for the first evaluation...",
            ha="center",
            va="center",
            transform=cls.summary_axis.transAxes,
        )
        cls.figure.tight_layout()
        cls.figure.canvas.draw()
        cls.figure.canvas.flush_events()

    @classmethod
    def update(cls, optimization):
        """Refresh the figure after one cost-function evaluation."""
        if cls.figure is None or not plt.fignum_exists(cls.figure.number):
            cls.initialize()

        cls.spectrum_axis.clear()
        for number, previous_best in enumerate(optimization.best_history[:-1]):
            cls.spectrum_axis.plot(
                cost_function.wavelengths,
                previous_best.modeled,
                "o-",
                color="0.85",
                alpha=0.8,
                label="Previous best" if number == 0 else None,
            )
        cls.spectrum_axis.plot(
            cost_function.wavelengths,
            cost_function.measured,
            "o-",
            label="Measured",
        )
        cls.spectrum_axis.plot(
            cost_function.wavelengths,
            optimization.best.modeled,
            "o-",
            label="Current best",
        )
        cls.spectrum_axis.set(
            xlabel="Wavelength [nm]",
            ylabel=r"Spectral radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]",
            title="Measured and best spectra",
            xscale="log",
            yscale="log",
        )
        cls.format_spectrum_axis()
        cls.spectrum_axis.legend()

        cls.summary_axis.clear()
        cls.summary_axis.set_title("Current best parameters")
        cls.summary_axis.axis("off")
        rows = [
            [
                cls.wrap_table_text(name, 12),
                f"{lower:#.1e}",
                f"{value:#.1e}",
                f"{upper:#.1e}",
                f"{value / machup:#.2g}" if np.isfinite(machup) and machup != 0 else "N/A",
            ]
            for name, lower, value, upper, machup in zip(
                parameter.name,
                parameter.lower_limit,
                optimization.best.parameter.value,
                parameter.upper_limit,
                parameter.machup,
            )
        ]
        table = cls.summary_axis.table(
            cellText=rows,
            colLabels=[
                cls.wrap_table_text(label, 10)
                for label in (
                    "Parameter",
                    "Lower limit",
                    "Best value",
                    "Upper limit",
                    "Best / matchup",
                )
            ],
            cellLoc="center",
            colLoc="center",
            loc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(0.9, 1.4)
        for row in range(len(rows) + 1):
            line_count = max(
                table[row, column].get_text().get_text().count("\n") + 1
                for column in range(5)
            )
            for column in range(5):
                cell = table[row, column]
                cell.set_height(cell.get_height() * line_count)
        cls.summary_axis.text(
            0.5,
            0.12,
            (
                f"Best cost-function value: {optimization.best.value:#.1e}\n"
                f"Iterations completed: {len(optimization.cost_history)}\n"
                f"Total time used: {perf_counter() - optimization.start_time:#.1e} s"
            ),
            ha="center",
            va="center",
            fontsize=11,
            transform=cls.summary_axis.transAxes,
        )

        cls.figure.tight_layout()
        cls.figure.canvas.draw_idle()
        cls.figure.canvas.flush_events()
        plt.pause(0.001)


class optimizer:
    """Optimize the model parameters by minimizing the cost function."""

    best = None
    best_history = []
    cost_history = []
    start_time = None
    plot = optimization_plot

    @classmethod
    def store_if_best(cls, candidate):
        if cls.best is None or candidate.value < cls.best.value:
            cls.best = candidate.copy()
            cls.best_history.append(cls.best)
    
    @classmethod
    def initialize(cls):
        cls.best = None
        cls.best_history = []
        cls.cost_history = []
        cls.start_time = perf_counter()
        cls.plot.initialize()
        
    @classmethod
    def update(cls, evaluated):
        cls.cost_history.append(evaluated.value)
        cls.store_if_best(evaluated)
        cls.plot.update(cls)
        
    @classmethod
    def lucky_guess(cls, n_iterations=10):
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
        cls.initialize()
        for _ in range(n_iterations):
            candidate = parameter.set_random()
            evaluated = cost_function.simple_least_square(candidate)
            cls.update(evaluated)
        return cls.best

    @classmethod
    def my_version_not_based_on_luck(cls):
        cls.initialize()
        tolerance = 1; # find suitable value from cost function
        happy = False
        while not happy:
            """
            better optimization ideas could be implemented here. AI could help.
            candidate = find_parameters_in_a_clever_way()
            """
            evaluated = cost_function.my_super_version(candidate)
            cls.update(evaluated)
            if cls.best.value < tolerance:
                happy = True
        return cls.best
    

if __name__ == "__main__":
    set_fixed_parameter()
    best = optimizer.lucky_guess(8)

    p = best.parameter.value
    output = np.column_stack((np.arange(1, len(p) + 1), p))
    flick.save_two_columns(output, "best_parameter_values.txt")
    optimizer.plot.figure.savefig("output/optimization_result.png", dpi=300)

    plt.ioff()
    plt.show()
