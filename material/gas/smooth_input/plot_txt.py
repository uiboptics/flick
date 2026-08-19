#!/usr/bin/env python3
from __future__ import annotations

"""Plot a two-column wavelength/cross-section text file."""
"""python3 plot_txt.py --nm --logy -o no2_uv_vis_T0_P1.txt"""

import argparse
import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib_cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(SCRIPT_DIR / ".cache"))


DEFAULT_FILE = "no2_uv_vis_T0_P1.txt"


def read_two_column_file(path: Path) -> tuple[list[float], list[float]]:
    wavelengths: list[float] = []
    cross_sections: list[float] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(("#", "//", "/*", "*", "*/")):
                continue

            parts = stripped.split()
            if len(parts) < 2:
                continue

            try:
                wavelength = float(parts[0])
                cross_section = float(parts[1])
            except ValueError:
                continue

            wavelengths.append(wavelength)
            cross_sections.append(cross_section)

    if not wavelengths:
        raise ValueError(f"No numeric two-column data found in {path}")

    return wavelengths, cross_sections


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot wavelength and cross-section columns from a text file."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=DEFAULT_FILE,
        help=f"text file to plot (default: {DEFAULT_FILE})",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="save the figure to this path instead of opening an interactive window",
    )
    parser.add_argument(
        "--nm",
        action="store_true",
        help="convert wavelengths from meters to nanometers on the x-axis",
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help="use a logarithmic scale for the y-axis",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.output:
        import matplotlib

        matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    input_path = Path(args.input)
    wavelengths, cross_sections = read_two_column_file(input_path)

    if args.nm:
        wavelengths = [value * 1e9 for value in wavelengths]
        x_label = "Wavelength [nm]"
    else:
        x_label = "Wavelength [m]"

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(wavelengths, cross_sections, linewidth=1.2)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Cross section")
    ax.set_title(input_path.name)
    ax.grid(True, alpha=0.3)

    if args.logy:
        ax.set_yscale("log")

    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=200)
    else:
        plt.show()


if __name__ == "__main__":
    main()
