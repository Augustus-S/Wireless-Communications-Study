"""Zadoff-Chu sequence properties and cyclic-shift detection example.

Run from the repository root:

    python3 zadoff-chu-sequence/examples/zc_sequence_example.py

The script writes plots to:

    zadoff-chu-sequence/examples/output/zc_sequence_properties.png
"""

from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache"))

import matplotlib.pyplot as plt
import numpy as np


def zadoff_chu_sequence(length: int, root: int) -> np.ndarray:
    """Generate a Zadoff-Chu sequence.

    For odd lengths, this uses x[n] = exp(-j*pi*u*n*(n+1)/N).
    For even lengths, this uses x[n] = exp(-j*pi*u*n*n/N).
    The root and length must be coprime for the ideal CAZAC property.
    """

    if length <= 0:
        raise ValueError("length must be positive")
    if math.gcd(root, length) != 1:
        raise ValueError("root and length must be coprime")

    n = np.arange(length)
    if length % 2 == 0:
        phase = -np.pi * root * n * n / length
    else:
        phase = -np.pi * root * n * (n + 1) / length
    return np.exp(1j * phase)


def cyclic_correlation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Compute cyclic correlation between equal-length complex sequences."""

    if a.shape != b.shape:
        raise ValueError("sequences must have the same length")

    length = a.size
    result = np.empty(length, dtype=np.complex128)
    for shift in range(length):
        result[shift] = np.vdot(np.roll(b, shift), a)
    return result


def main() -> None:
    length = 63
    root = 25
    other_root = 29
    true_shift = 17
    rng = np.random.default_rng(7)

    zc = zadoff_chu_sequence(length, root)
    zc_other = zadoff_chu_sequence(length, other_root)

    magnitude = np.abs(zc)
    auto_corr = cyclic_correlation(zc, zc) / length
    cross_corr = cyclic_correlation(zc, zc_other) / length

    spectrum = np.fft.fftshift(np.fft.fft(zc)) / np.sqrt(length)
    spectrum_magnitude = np.abs(spectrum)

    received = np.roll(zc, true_shift)
    noise_power = 0.15
    noise = noise_power * (rng.standard_normal(length) + 1j * rng.standard_normal(length))
    received_noisy = received + noise
    detection_corr = cyclic_correlation(received_noisy, zc) / length
    estimated_shift = int(np.argmax(np.abs(detection_corr)))

    nonzero_shifts = np.arange(1, length)
    max_magnitude_error = float(np.max(np.abs(magnitude - 1.0)))
    max_auto_sidelobe = float(np.max(np.abs(auto_corr[nonzero_shifts])))
    max_cross_corr = float(np.max(np.abs(cross_corr)))
    detection_peak = float(np.max(np.abs(detection_corr)))

    print("Zadoff-Chu sequence example")
    print(f"length N              = {length}")
    print(f"root u                = {root}")
    print(f"other root            = {other_root}")
    print(f"max | |x[n]| - 1 |    = {max_magnitude_error:.3e}")
    print(f"max autocorr sidelobe = {max_auto_sidelobe:.3e}")
    print(f"max cross-correlation = {max_cross_corr:.3e}")
    print(f"true cyclic shift     = {true_shift}")
    print(f"estimated shift       = {estimated_shift}")
    print(f"detection peak        = {detection_peak:.3f}")

    output_dir = Path(__file__).resolve().parent / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "zc_sequence_properties.png"

    shifts = np.arange(length)
    freq_bins = np.arange(-length // 2, length // 2 + 1)
    if freq_bins.size > length:
        freq_bins = freq_bins[:length]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    axes[0, 0].plot(np.real(zc), np.imag(zc), "o", markersize=4)
    unit_circle = np.exp(1j * np.linspace(0.0, 2.0 * np.pi, 400))
    axes[0, 0].plot(np.real(unit_circle), np.imag(unit_circle), "--", linewidth=1)
    axes[0, 0].set_title("ZC Sequence Samples on Unit Circle")
    axes[0, 0].set_xlabel("Real")
    axes[0, 0].set_ylabel("Imag")
    axes[0, 0].axis("equal")
    axes[0, 0].grid(True)

    axes[0, 1].stem(shifts, np.abs(auto_corr), basefmt=" ")
    axes[0, 1].set_title("Normalized Cyclic Autocorrelation")
    axes[0, 1].set_xlabel("Cyclic shift")
    axes[0, 1].set_ylabel("Magnitude")
    axes[0, 1].grid(True)

    axes[1, 0].plot(freq_bins, spectrum_magnitude, linewidth=1.5)
    axes[1, 0].set_title("DFT Magnitude")
    axes[1, 0].set_xlabel("Shifted DFT bin")
    axes[1, 0].set_ylabel("Magnitude")
    axes[1, 0].grid(True)

    axes[1, 1].stem(shifts, np.abs(detection_corr), basefmt=" ")
    axes[1, 1].axvline(true_shift, color="tab:red", linestyle="--", label="true shift")
    axes[1, 1].axvline(estimated_shift, color="tab:green", linestyle=":", label="estimated")
    axes[1, 1].set_title("Cyclic Shift Detection with Noise")
    axes[1, 1].set_xlabel("Candidate shift")
    axes[1, 1].set_ylabel("Correlation magnitude")
    axes[1, 1].legend()
    axes[1, 1].grid(True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    print(f"saved plot            = {output_path}")


if __name__ == "__main__":
    main()
