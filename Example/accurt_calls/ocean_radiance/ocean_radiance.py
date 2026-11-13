"""
Plot the two hemispheres of an ocean radiance field as polar
contours.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
import numpy as np
import flick
flick.setup_script(__file__)

n_polar = 90
n_azimuth = 150
center_wavelength = 515e-9
band_width = 10e-9

def to_3d_array(radiance_2d):
    # Reshape stacked wavelength blocks to wavelength, polar, azimuth.
    radiance = np.asarray(radiance_2d).squeeze()
    wavelength_count = radiance.shape[0] // n_polar
    return radiance.reshape(wavelength_count, n_polar, n_azimuth)

def calculate_radiance():
    config = flick.accurt_config()
    config.set_streams(32)
    config.set("n_heights", 10)  # Atmosphere layers
    config.set("detector_height", -5)

    config.set("aerosol_od", 0.06)  # Typical Arctic value
    config.set("cloud_liquid", 0)  # Clear blue sky

    # Turbid Beaufort Sea material stored in marine-particle IOP tables.
    config.set("mp_names", "PB1819_K1c_s")
    config.set("mp_concentrations", 2e-3)  # From metadata
    config.set("cdom_440", 0.1)

    # Time and position of IOP measurements in the Beaufort Sea.
    time_pos = [2019, 8, 11, 1, 12, 0, 70.60, -149.04]
    config.set("toa_solar_multiplication", time_pos + [band_width])
    
    wavelength_low = center_wavelength - band_width / 2
    wavelength_high = center_wavelength + band_width / 2
    config.set("wavelengths", np.linspace(wavelength_low, wavelength_high, 3))
    config.set("detector_radiance_distribution_override",[n_polar, n_azimuth])

    solar_zenith_angle = flick.solar_zenith_angle(time_pos)
    print(f"Solar zenith angle: {solar_zenith_angle:.1f} deg")

    radiance = flick.run(f"accurt {config.filename}")
    return to_3d_array(radiance)
    
    
class plot:
    def format_radiance(self, value, position=None):
        del position
        return f"{value:.2g}"

    def close_azimuth(self, radiance):
        """Repeat the first azimuth at 360 degrees to remove the polar seam."""
        return np.concatenate((radiance, radiance[:1, :]), axis=0)

    def make_plot(self):
        """Return a two-panel polar contour plot of the angular light field."""
        raw = calculate_radiance()
        middle_wavelength = raw.shape[0] // 2
        radiance = raw[middle_wavelength].T * 1e-6  # W/m to mW/nm

        polar = np.linspace(0.0, np.pi, n_polar)
        azimuth = np.linspace(-np.pi, np.pi, n_azimuth)
        azimuth_closed = np.append(azimuth, azimuth[0] + 2 * np.pi)
        radiance_closed = self.close_azimuth(radiance)

        upward = polar <= np.pi / 2
        downward = polar >= np.pi / 2
        panels = (
            (
                "Downwelling radiance",
                np.pi - polar[downward][::-1],
                radiance_closed[:, downward][:, ::-1],
            ),
            (
                "Upwelling radiance",
                polar[upward],
                radiance_closed[:, upward],
            ),
        )

        right_positive = panels[1][2][panels[1][2] > 0]
        if right_positive.size == 0:
            raise ValueError("Upwelling radiance contains no positive values")
        shared_boundary = right_positive.max()

        left_positive = panels[0][2][panels[0][2] > 0]
        if left_positive.size == 0:
            raise ValueError("Downwelling radiance contains no positive values")
        if shared_boundary >= left_positive.max():
            raise ValueError(
                "The upwelling maximum must be below the downwelling maximum"
            )

        color_limits = (
            (shared_boundary, left_positive.max()),
            (right_positive.min(), shared_boundary),
        )
        color_maps = ("hot", "hot")

        figure, axes = plt.subplots(
            1, 2, figsize=(10.5, 4.8), subplot_kw={"projection": "polar"}
        )
        figure.subplots_adjust(
            left=0.06, right=0.94, bottom=0.12, top=0.88, wspace=0.42
        )

        for axis, (label, radius, values), (value_min, value_max), color_map in zip(
            axes, panels, color_limits, color_maps
        ):
            theta_grid, radius_grid = np.meshgrid(
                azimuth_closed, radius, indexing="ij"
            )
            positive_values = values[values > 0]
            if positive_values.size == 0:
                raise ValueError(f"{label} contains no positive radiance values")

            levels = np.geomspace(value_min, value_max, 13)
            plotted_values = np.ma.clip(
                np.ma.masked_less_equal(values, 0), value_min, value_max
            )
            contour = axis.contourf(
                theta_grid,
                radius_grid,
                plotted_values,
                levels=levels,
                norm=LogNorm(vmin=value_min, vmax=value_max),
                cmap=color_map,
            )
            axis.set_theta_zero_location("E")
            axis.set_theta_direction(1)
            axis.set_ylim(0, np.pi / 2)
            axis.set_thetagrids(np.arange(0, 360, 45))
            axis.set_rgrids(
                np.deg2rad([30, 60, 90]),
                labels=[r"30$^\circ$", r"60$^\circ$", r"90$^\circ$"],
            )
            axis.set_rlabel_position(25)
            axis.grid(color="0.65", linewidth=0.8)
            axis.set_title(label, pad=20)

            colorbar = figure.colorbar(contour, ax=axis, pad=0.10, shrink=0.88)
            colorbar.formatter = FuncFormatter(self.format_radiance)
            colorbar.update_ticks()
            colorbar.set_label(
                r"Radiance [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]"
            )
        return figure


if __name__ == "__main__":
    plot = plot()
    figure = plot.make_plot()
    figure.savefig("ocean_radiance.png", dpi=500, bbox_inches="tight")
    plt.show()
