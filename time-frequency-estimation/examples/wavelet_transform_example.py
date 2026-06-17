"""Continuous wavelet transform example using a Morlet wavelet."""

from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", str(Path("outputs/.matplotlib")))
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


def make_multiscale_signal(fs: float, duration: float) -> tuple[np.ndarray, np.ndarray]:
    """Create a signal with a low-frequency component and a short high-frequency burst."""
    t = np.arange(0.0, duration, 1.0 / fs)
    low = np.sin(2.0 * np.pi * 20.0 * t)
    burst_window = np.exp(-0.5 * ((t - 1.1) / 0.08) ** 2)
    burst = 1.4 * burst_window * np.sin(2.0 * np.pi * 140.0 * t)
    pulse = signal.gausspulse(t - 1.55, fc=90.0, bw=0.35)
    return t, low + burst + pulse


def morlet_cwt(
    x: np.ndarray,
    fs: float,
    frequencies: np.ndarray,
    omega0: float = 6.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute a Morlet CWT without requiring PyWavelets."""
    scales = omega0 * fs / (2.0 * np.pi * frequencies)
    coefficients = np.empty((len(frequencies), len(x)), dtype=np.complex128)

    for row, scale in enumerate(scales):
        half_width = int(np.ceil(5.0 * scale))
        n = np.arange(-half_width, half_width + 1)
        wavelet = np.pi ** (-0.25) * np.exp(1j * omega0 * n / scale)
        wavelet *= np.exp(-0.5 * (n / scale) ** 2) / np.sqrt(scale)
        coefficients[row] = signal.fftconvolve(x, np.conj(wavelet[::-1]), mode="same")

    return coefficients, scales


def main() -> None:
    fs = 1_000.0
    duration = 2.0
    t, x = make_multiscale_signal(fs, duration)
    frequencies = np.linspace(10.0, 220.0, 160)
    coefficients, scales = morlet_cwt(x, fs, frequencies)
    power = np.abs(coefficients) ** 2

    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), constrained_layout=True)
    axes[0].plot(t, x, linewidth=0.9)
    axes[0].set_title("Multiscale signal")
    axes[0].set_xlabel("Time [s]")
    axes[0].set_ylabel("Amplitude")

    mesh = axes[1].pcolormesh(t, frequencies, power, shading="auto")
    axes[1].set_title("Morlet CWT scalogram")
    axes[1].set_xlabel("Time [s]")
    axes[1].set_ylabel("Pseudo-frequency [Hz]")
    fig.colorbar(mesh, ax=axes[1], label="Power")

    output_path = output_dir / "wavelet_transform_example.png"
    fig.savefig(output_path, dpi=160)
    print(f"Saved {output_path}")
    print(f"Scale range: {scales.min():.2f} to {scales.max():.2f} samples")


if __name__ == "__main__":
    main()
