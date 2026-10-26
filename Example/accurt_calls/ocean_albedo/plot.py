"""Plot ocean albedo versus bubble volume fraction.

Run ``compute.py`` first to generate the input data.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "output" / "albedo_versus_volume_fraction.txt"
OUTPUT_FILE = BASE_DIR / "output" / "albedo_versus_volume_fraction.pdf"


def main() -> None:
    data = np.loadtxt(INPUT_FILE, ndmin=2)
    volume_fraction = data[:, 0]
    albedo = data[:, 1]

    fig, ax = plt.subplots()
    ax.plot(volume_fraction, albedo, "o-")
    ax.set_xlabel("Nano-bubble volume fraction")
    ax.set_ylabel("Albedo")
    ax.set_title("Ocean broadband plane irradiance albedo")
    ax.set_xscale("log")
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks(np.arange(0.0, 1.01, 0.1))
    ax.yaxis.set_major_formatter("{x:.1f}")
    ax.minorticks_on()
    ax.grid(True, which="both")
    fig.tight_layout()
    fig.savefig(OUTPUT_FILE)
    plt.show()


if __name__ == "__main__":
    main()
