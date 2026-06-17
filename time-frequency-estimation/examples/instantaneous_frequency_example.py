"""Instantaneous frequency estimation from analytic phase."""

from pathlib import Path
import os

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT_DIR / ".matplotlib"))
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


def make_noisy_chirp(
    fs: float,
    duration: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create a noisy linear chirp and its true instantaneous frequency."""
    t = np.arange(0.0, duration, 1.0 / fs)
    f0 = 35.0
    f1 = 180.0
    k = (f1 - f0) / duration
    phase = 2.0 * np.pi * (f0 * t + 0.5 * k * t**2)
    clean = np.cos(phase)
    noise = 0.18 * rng.standard_normal(len(t))
    return t, clean + noise, f0 + k * t


def estimate_instantaneous_frequency(
    x: np.ndarray,
    fs: float,
    smoothing_length: int = 31,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate instantaneous frequency from the derivative of analytic phase."""
    z = signal.hilbert(x)
    phase = np.unwrap(np.angle(z))
    instantaneous_frequency = np.gradient(phase) * fs / (2.0 * np.pi)

    if smoothing_length % 2 == 0:
        smoothing_length += 1
    smoothed = signal.medfilt(instantaneous_frequency, kernel_size=smoothing_length)
    return phase, instantaneous_frequency, smoothed


def main() -> None:
    fs = 1_000.0
    duration = 2.0
    rng = np.random.default_rng(7)
    t, x, true_if = make_noisy_chirp(fs, duration, rng)
    phase, raw_if, smoothed_if = estimate_instantaneous_frequency(x, fs)

    OUTPUT_DIR.mkdir(exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(9, 7), constrained_layout=True)
    axes[0].plot(t, x, linewidth=0.8)
    axes[0].set_title("Noisy chirp")
    axes[0].set_xlabel("Time [s]")
    axes[0].set_ylabel("Amplitude")

    axes[1].plot(t, phase, linewidth=0.9)
    axes[1].set_title("Unwrapped phase")
    axes[1].set_xlabel("Time [s]")
    axes[1].set_ylabel("Phase [rad]")

    axes[2].plot(t, raw_if, color="0.75", linewidth=0.8, label="Raw estimate")
    axes[2].plot(t, smoothed_if, "b", linewidth=1.2, label="Median filtered")
    axes[2].plot(t, true_if, "r--", linewidth=1.2, label="True IF")
    axes[2].set_ylim(0, 230)
    axes[2].set_title("Instantaneous frequency")
    axes[2].set_xlabel("Time [s]")
    axes[2].set_ylabel("Frequency [Hz]")
    axes[2].legend(loc="upper left")

    output_path = OUTPUT_DIR / "instantaneous_frequency_example.png"
    fig.savefig(output_path, dpi=160)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
