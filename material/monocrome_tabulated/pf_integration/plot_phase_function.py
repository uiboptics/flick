#!/usr/bin/env python3
"""Plot a tabulated phase function in mu and in a stretched coordinate u."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_phase_function(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load angle (degrees) and phase function from the repository's text format.

    The angle occupies the first eight characters.  Reading it as fixed-width is
    important because, for large phase values, there is no separating space.
    """
    angles: list[float] = []
    phase: list[float] = []

    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            try:
                fields = line.split()
                if len(fields) == 2:
                    angle_text, phase_text = fields
                else:
                    # Legacy tables use an eight-character angle field and
                    # sometimes omit whitespace before a large phase value.
                    angle_text, phase_text = line[:8], line[8:]
                angles.append(float(angle_text))
                phase.append(float(phase_text))
            except ValueError as error:
                raise ValueError(
                    f"Could not parse {path} line {line_number}: {line.rstrip()!r}"
                ) from error

    if len(angles) < 2:
        raise ValueError(f"{path} must contain at least two tabulated points")

    return np.asarray(angles), np.asarray(phase)


def phase_in_half_angle_coordinate(
    theta_deg: np.ndarray, phase: np.ndarray, q: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return the original samples sorted in u, where u=sin(theta/2)^(2/q)."""
    if q <= 0:
        raise ValueError("q must be positive")
    if theta_deg.shape != phase.shape:
        raise ValueError("theta_deg and phase must have the same shape")

    mu = np.cos(np.deg2rad(theta_deg))
    # Clipping prevents tiny round-off excursions outside [0, 1].
    x = np.clip((1.0 - mu) / 2.0, 0.0, 1.0)
    u = x ** (1.0 / q)
    order = np.argsort(u)
    return u[order], phase[order]


def half_angle_quadrature(
    theta_deg: np.ndarray,
    phase: np.ndarray,
    q: float = 4.0,
    points: int = 200,
    method: str = "analytic",
) -> tuple[float, float, float]:
    """Compute normalization N, forward-peak-suppressed D, and g=1-D/N.

    The phase function is interpolated linearly in the transformed coordinate
    ``u=((1-mu)/2)**(1/q)``.  ``method="analytic"`` integrates that
    piecewise-linear interpolant exactly.  The alternatives are ``"gauss"``
    for Gauss--Legendre quadrature and ``"uniform"`` for a uniform-u
    trapezoidal rule.  If the table does not reach u=0 or u=1, the nearest
    tabulated value is used over the uncovered interval.
    """
    if points < 2:
        raise ValueError("points must be at least 2")

    u_data, phase_u = phase_in_half_angle_coordinate(theta_deg, phase, q)

    if method == "analytic":
        # Add the endpoints with constant extrapolation, then integrate
        # p(u)=slope*u+intercept exactly on every tabulated interval.
        if u_data[0] > 0.0:
            u_data = np.concatenate(([0.0], u_data))
            phase_u = np.concatenate(([phase_u[0]], phase_u))
        if u_data[-1] < 1.0:
            u_data = np.concatenate((u_data, [1.0]))
            phase_u = np.concatenate((phase_u, [phase_u[-1]]))

        def integrate_power_times_linear(power: float) -> float:
            lower = u_data[:-1]
            upper = u_data[1:]
            slope = np.diff(phase_u) / np.diff(u_data)
            intercept = phase_u[:-1] - slope * lower
            slope_part = slope * (upper ** (power + 2.0) - lower ** (power + 2.0))
            constant_part = intercept * (
                upper ** (power + 1.0) - lower ** (power + 1.0)
            )
            return float(
                np.sum(slope_part / (power + 2.0) + constant_part / (power + 1.0))
            )

        normalization_integral = integrate_power_times_linear(q - 1.0)
        deficit_integral = integrate_power_times_linear(2.0 * q - 1.0)
    elif method == "gauss":
        nodes, weights = np.polynomial.legendre.leggauss(points)
        u_eval = 0.5 * (nodes + 1.0)
        weights = 0.5 * weights  # Map weights from [-1, 1] to [0, 1].
    elif method == "uniform":
        u_eval = np.linspace(0.0, 1.0, points)
        weights = None
    else:
        raise ValueError("method must be 'analytic', 'gauss', or 'uniform'")

    if method != "analytic":
        p_eval = np.interp(
            u_eval, u_data, phase_u, left=phase_u[0], right=phase_u[-1]
        )
        normalization_integrand = u_eval ** (q - 1.0) * p_eval
        deficit_integrand = u_eval ** (2.0 * q - 1.0) * p_eval

        if weights is None:
            normalization_integral = np.trapezoid(normalization_integrand, u_eval)
            deficit_integral = np.trapezoid(deficit_integrand, u_eval)
        else:
            normalization_integral = np.dot(weights, normalization_integrand)
            deficit_integral = np.dot(weights, deficit_integrand)

    normalization = 4.0 * np.pi * q * normalization_integral
    deficit = 8.0 * np.pi * q * deficit_integral
    if not np.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("phase function has a non-positive or invalid normalization")

    g = 1.0 - deficit / normalization
    return float(normalization), float(deficit), float(g)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot a phase function in mu=cos(theta) and transformed u."
    )
    parser.add_argument(
        "file",
        nargs="?",
        type=Path,
        default=Path("../phase_functions/vsfPart1_bp0.1.txt"),
        help="phase-function table (default: ../phase_functions/vsfPart1_bp0.1.txt)",
    )
    parser.add_argument(
        "--q",
        type=float,
        choices=(2.0, 3,4,5.0),
        default=4.0,
        help="half-angle stretching exponent (default: 2)",
    )
    parser.add_argument(
        "--points",
        type=int,
        default=300,
        help="number of points in each interpolated curve (default: 2000)",
    )
    parser.add_argument(
        "--linear-y", action="store_true", help="use a linear instead of logarithmic y-axis"
    )
    parser.add_argument(
        "--quadrature",
        choices=("analytic", "gauss", "uniform"),
        default="analytic",
        help="integration rule in u (default: exact analytic integration of the interpolant)",
    )
    parser.add_argument("--save", type=Path, help="save the figure instead of opening it")
    args = parser.parse_args()

    if args.q <= 0:
        parser.error("--q must be positive")
    if args.points < 2:
        parser.error("--points must be at least 2")

    theta_deg, phase = load_phase_function(args.file)
    mu = np.cos(np.deg2rad(theta_deg))

    normalization, deficit, g = half_angle_quadrature(
        theta_deg, phase, q=args.q, points=args.points, method=args.quadrature
    )
    print(f"quadrature = {args.quadrature}, q = {args.q:g}, points = {args.points}")
    print(f"N = {normalization:.12g}")
    print(f"D = 1-g numerator = {deficit:.12g}")
    print(f"g = 1-D/N = {g:.12g}")

    # np.interp requires increasing coordinates.  Angle increases while mu
    # decreases, so sort rather than assuming a particular input ordering.
    mu_order = np.argsort(mu)
    mu_data = mu[mu_order]
    phase_mu = phase[mu_order]
    mu_grid = np.linspace(mu_data[0], mu_data[-1], args.points)
    ordinary_interp = np.interp(mu_grid, mu_data, phase_mu)

    # Transform the original tabulated coordinates before interpolating.
    # u increases with angle (and decreases with mu), hence another explicit sort.
    u_data, phase_u = phase_in_half_angle_coordinate(theta_deg, phase, args.q)
    u_grid = np.linspace(0.0, 1.0, args.points)
    transformed_interp = np.interp(
        u_grid, u_data, phase_u, left=phase_u[0], right=phase_u[-1]
    )
    normalization_integrand = (
        4.0 * np.pi * args.q * u_grid ** (args.q - 1.0) * transformed_interp
    )
    deficit_integrand = (
        8.0 * np.pi * args.q * u_grid ** (2.0 * args.q - 1.0) * transformed_interp
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)

    axes[0, 0].plot(mu_grid, ordinary_interp, label="linear interpolation in $\\mu$")
    axes[0, 0].scatter(mu, phase, s=15, color="black", zorder=3, label="tabulated points")
    axes[0, 0].set_xlabel(r"$\mu=\cos(\theta)$")
    axes[0, 0].set_ylabel("Phase function")
    axes[0, 0].set_title("Original coordinate")

    axes[0, 1].plot(
        u_grid, transformed_interp, color="tab:orange", label="linear interpolation in $u$"
    )
    axes[0, 1].scatter(
        u_data, phase_u, s=15, color="black", zorder=3, label="transformed points"
    )
    axes[0, 1].set_xlabel(r"$u=((1-\mu)/2)^{1/q}$")
    axes[0, 1].set_ylabel("Phase function")
    axes[0, 1].set_title(fr"Transformed coordinate ($q={args.q:g}$)")

    axes[1, 0].plot(u_grid, normalization_integrand, color="tab:green")
    axes[1, 0].fill_between(u_grid, normalization_integrand, color="tab:green", alpha=0.15)
    axes[1, 0].set_xlabel("$u$")
    axes[1, 0].set_ylabel(r"$4\pi q\,u^{q-1}p(1-2u^q)$")
    axes[1, 0].set_title(fr"Normalization integrand ($N={normalization:.6g}$)")

    axes[1, 1].plot(u_grid, deficit_integrand, color="tab:red")
    axes[1, 1].fill_between(u_grid, deficit_integrand, color="tab:red", alpha=0.15)
    axes[1, 1].set_xlabel("$u$")
    axes[1, 1].set_ylabel(r"$8\pi q\,u^{2q-1}p(1-2u^q)$")
    axes[1, 1].set_title(fr"Asymmetry-deficit integrand ($D={deficit:.6g}$, $g={g:.6g}$)")

    for axis in axes.flat:
        axis.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend()
    axes[0, 1].legend()
    if not args.linear_y:
        axes[0, 0].set_yscale("log")
        axes[0, 1].set_yscale("log")

    fig.suptitle(args.file.name)
    if args.save:
        fig.savefig(args.save, dpi=180)
        print(f"Saved {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
