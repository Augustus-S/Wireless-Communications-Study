"""STFT time-frequency estimation example."""

from pathlib import Path
import os

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT_DIR / ".matplotlib"))
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


def make_chirp(fs: float, duration: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create a linear chirp with known instantaneous frequency."""
    t = np.arange(0.0, duration, 1.0 / fs)
    f0 = 40.0
    f1 = 220.0
    k = (f1 - f0) / duration
    phase = 2.0 * np.pi * (f0 * t + 0.5 * k * t**2)
    x = np.cos(phase)
    true_if = f0 + k * t
    return t, x, true_if


def compute_stft(
    x: np.ndarray,
    fs: float,
    window_length: int = 256,
    overlap: int = 192,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute STFT and return frequency bins, time bins, and power."""
    frequencies, times, spectrum = signal.stft(
        x,
        fs=fs,
        window="hann",
        nperseg=window_length,
        noverlap=overlap,
        nfft=512,
        boundary=None,
    )
    power = np.abs(spectrum) ** 2
    return frequencies, times, power


def estimate_ridge(frequencies: np.ndarray, power: np.ndarray) -> np.ndarray:
    """Estimate a simple frequency ridge by selecting max power per time bin."""
    ridge_indices = np.argmax(power, axis=0)
    return frequencies[ridge_indices]


def main() -> None:
    fs = 1_000.0
    duration = 2.0
    t, x, true_if = make_chirp(fs, duration)
    frequencies, times, power = compute_stft(x, fs)
    ridge = estimate_ridge(frequencies, power)

    OUTPUT_DIR.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), constrained_layout=True)
    axes[0].plot(t, x, linewidth=0.9)
    axes[0].set_title("Linear chirp signal")
    axes[0].set_xlabel("Time [s]")
    axes[0].set_ylabel("Amplitude")

    mesh = axes[1].pcolormesh(times, frequencies, 10 * np.log10(power + 1e-12), shading="auto")
    axes[1].plot(times, ridge, "w--", linewidth=1.4, label="STFT ridge")
    axes[1].plot(t, true_if, "r", linewidth=1.0, alpha=0.8, label="True IF")
    axes[1].set_ylim(0, 260)
    axes[1].set_title("STFT spectrogram")
    axes[1].set_xlabel("Time [s]")
    axes[1].set_ylabel("Frequency [Hz]")
    axes[1].legend(loc="upper left")
    fig.colorbar(mesh, ax=axes[1], label="Power [dB]")

    output_path = OUTPUT_DIR / "stft_example.png"
    fig.savefig(output_path, dpi=160)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
