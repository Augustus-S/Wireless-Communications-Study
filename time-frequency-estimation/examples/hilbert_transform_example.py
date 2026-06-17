"""Hilbert transform and analytic signal example."""

from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", str(Path("outputs/.matplotlib")))
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


def make_am_signal(fs: float, duration: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create an amplitude-modulated narrowband signal."""
    t = np.arange(0.0, duration, 1.0 / fs)
    envelope = 1.0 + 0.45 * np.sin(2.0 * np.pi * 3.0 * t)
    carrier = np.cos(2.0 * np.pi * 80.0 * t)
    return t, envelope * carrier, envelope


def analytic_signal(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return analytic signal, amplitude envelope, and unwrapped phase."""
    z = signal.hilbert(x)
    envelope = np.abs(z)
    phase = np.unwrap(np.angle(z))
    return z, envelope, phase


def main() -> None:
    fs = 1_000.0
    duration = 1.0
    t, x, true_envelope = make_am_signal(fs, duration)
    z, envelope, phase = analytic_signal(x)

    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(9, 7), constrained_layout=True)
    axes[0].plot(t, x, linewidth=0.8, label="Real signal")
    axes[0].plot(t, envelope, "r", linewidth=1.2, label="Estimated envelope")
    axes[0].plot(t, -envelope, "r", linewidth=1.2)
    axes[0].plot(t, true_envelope, "k--", linewidth=1.0, label="True envelope")
    axes[0].set_title("Hilbert envelope")
    axes[0].set_xlabel("Time [s]")
    axes[0].set_ylabel("Amplitude")
    axes[0].legend(loc="upper right")

    axes[1].plot(t, phase, linewidth=1.0)
    axes[1].set_title("Unwrapped analytic phase")
    axes[1].set_xlabel("Time [s]")
    axes[1].set_ylabel("Phase [rad]")

    axes[2].plot(np.real(z), np.imag(z), linewidth=0.8)
    axes[2].set_title("Analytic signal in complex plane")
    axes[2].set_xlabel("Real part")
    axes[2].set_ylabel("Hilbert imaginary part")
    axes[2].axis("equal")

    output_path = output_dir / "hilbert_transform_example.png"
    fig.savefig(output_path, dpi=160)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
