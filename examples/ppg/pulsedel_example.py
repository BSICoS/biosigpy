"""Delineate the apex, basal, mid-amplitude and onset points of PPG pulses."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt, lfilter

from biosigpy.ppg import pulsedel, pulsedet
from biosigpy.tools import lpdfilter


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"


def main():
    # Load PPG signal from fixtures
    ppg_data = np.genfromtxt(
        EXAMPLES_ROOT / "fixtures/ppg/ppg_signals.csv", delimiter=",", names=True,
    )
    fs = 1000

    # Use only the first 2 minutes of the signal for demonstration
    signal = ppg_data["sig"][:2 * 60 * fs]
    t = ppg_data["t"][:2 * 60 * fs]

    # Apply high-pass filter to remove baseline drift
    b_high, a_high = butter(4, 0.5 / (fs / 2), btype="high")
    padlen = 3 * (max(len(b_high), len(a_high)) - 1)
    signal = filtfilt(b_high, a_high, signal, padlen=padlen)

    # Apply LPD (Low-Pass Differentiator) filter
    fp_lpd = 7.8        # Pass-band frequency (Hz)
    fc_lpd = 8          # Cut-off frequency (Hz)
    order_lpd = 100     # Filter order (samples)

    # Generate LPD filter and apply it
    b, delay = lpdfilter(fs, fc_lpd, pass_frequency=fp_lpd, order=order_lpd)
    delay = int(delay)
    dppg = lfilter(b, [1], signal)
    dppg = np.r_[dppg[delay:], np.zeros(delay)]

    # Reuse the aligned derivative for pulse detection and onset delineation.
    t_d, _ = pulsedet(dppg, fs)

    # window_a searches forward from t_d for the apex (default 0.25 s).
    # window_b searches backward from t_d for the basal point (default 0.15 s).
    # points="ABMO" (default) returns apex, basal, mid-amplitude and onset times.
    # Only onset uses dppg. If requested without it, t_o contains one NaN per pulse.
    t_a, t_b, t_m, t_o = pulsedel(
        signal, fs, t_d, dppg=dppg, window_a=0.25, window_b=0.15, points="ABMO"
    )

    # To request a subset, list the letters in the desired output order:
    # t_a, t_b = pulsedel(signal, fs, t_d, points="AB")
    # t_a = pulsedel(signal, fs, t_d, points="A")  # One letter returns one array.

    # Show the same fiducials on both signals.
    figure, axes = plt.subplots(2, 1, sharex=True)
    points = [
        (t_d, "o", (.47, .67, .19), r"$t_D$"),
        (t_a, "v", (.00, .45, .74), r"$t_A$"),
        (t_b, "^", (.00, .45, .74), r"$t_B$"),
        (t_m, "s", (.85, .33, .10), r"$t_M$"),
        (t_o, "d", (.49, .18, .56), r"$t_O$"),
    ]
    for axis, samples, label in zip(axes, (signal, dppg), ("Original PPG", "LPD-Filtered PPG")):
        axis.plot(t, samples, "k", linewidth=1, label=label)
        for times, marker, color, name in points:
            times = times[np.isfinite(times)]
            indices = np.floor(times * fs + .5).astype(int)
            axis.plot(times, samples[indices], marker, color=color, label=name)
        axis.set_xlabel("Time (s)")
        axis.set_ylabel("Amplitude")
        axis.grid(True)
        axis.legend()
    axes[0].set_title("PPG Delineation")
    axes[1].axvline(0, color="k", linestyle=":")
    figure.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
