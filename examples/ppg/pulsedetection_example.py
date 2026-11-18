"""Detect PPG pulses and demonstrate optional recovery after a gain change.

Run from the repository root with:
    python examples/ppg/pulsedetection_example.py
"""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt, lfilter

from biosigpy.ppg import pulsedetection
from biosigpy.tools import lpdfilter


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"
if str(EXAMPLES_ROOT) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_ROOT))

from _support.figures import add_save_figure_argument, show_or_save_figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_save_figure_argument(parser)
    args = parser.parse_args()

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
    signal_filtered = lfilter(b, [1], signal)
    signal_filtered = np.r_[signal_filtered[delay:], np.zeros(delay)]

    # Run pulse detection on filtered signal
    t_d, threshold = pulsedetection(signal_filtered, fs)
    detection_samples = np.floor(t_d * fs + .5).astype(int)

    # Plot results
    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 7))
    marker_color = (.47, .67, .19)
    axes[0].plot(t, signal, "k", linewidth=1, label="Original PPG")
    axes[0].plot(t_d, signal[detection_samples], "o", color=marker_color, label=r"$t_D$")
    axes[0].set_title("PPG Detection")
    axes[1].plot(t, signal_filtered, "k", linewidth=1, label="LPD-Filtered PPG")
    axes[1].plot(t, threshold, linewidth=1, label="Adaptive Threshold")
    axes[1].plot(t_d, signal_filtered[detection_samples], "o",
                 color=marker_color, label=r"$t_D$")
    for axis in axes:
        axis.set_xlabel("Time (s)")
        axis.set_ylabel("Amplitude")
        axis.grid(True)
        axis.legend()
    figure.tight_layout()

    # Optional feature: reinitialize after a device gain change
    # reset_after is disabled by default (None): the detector keeps its adaptive
    # state even during long periods without detections. Set a finite duration in
    # seconds to allow reinitialization after that time without a detected pulse.
    #
    # One possible application is a device that changes its acquisition gain.
    # If the gain decreases substantially, the previous threshold floor can remain
    # above the new pulses and prevent the detector from recovering on its own.
    # Here we simulate an interruption followed by a gain reduction to demonstrate
    # how enabling reset_after can help. This is an illustrative setting, not a
    # universal reset duration for all devices or recordings.
    #
    # The pulse morphology comes from real PPG: concatenate three fixture copies
    # to obtain 15 minutes. The interruption and gain change are introduced below;
    # they are not events observed in the original recording. Filter first so the
    # inserted NaNs do not propagate through the preprocessing filters.
    recovery_fs = fs
    recovery_signal = np.tile(ppg_data["sig"], 3)
    recovery_signal = filtfilt(b_high, a_high, recovery_signal, padlen=padlen)
    recovery_signal = lfilter(b, [1], recovery_signal)
    recovery_signal = np.r_[recovery_signal[delay:], np.zeros(delay)]
    recovery_time = np.arange(recovery_signal.size) / recovery_fs
    reset_after = 5

    # Control: run the same recorded signal with the optional feature disabled
    # and enabled, before introducing the simulated interruption and gain change.
    t_d_original_without_reset, _ = pulsedetection(recovery_signal, recovery_fs)
    t_d_original_with_reset, _ = pulsedetection(
        recovery_signal, recovery_fs, reset_after=reset_after
    )

    # Simulate 10 seconds of unavailable data, followed by a device gain of 0.05
    # times its previous value. Scaling the filtered waveform preserves the real
    # pulse shapes while isolating the amplitude change from filter transients.
    gap_start = ppg_data.size / recovery_fs
    gap_end = gap_start + 10
    amplitude_factor = .05
    recovery_signal[(recovery_time >= gap_start) & (recovery_time < gap_end)] = np.nan
    recovery_signal[recovery_time >= gap_end] *= amplitude_factor

    # Compare the identical modified input in both conditions. Without reset,
    # the old threshold may remain too high. With reset_after = 5 s, the missing
    # span has already exceeded the timeout when data returns, so the detector
    # reinitializes using the available signal at its new amplitude.
    t_d_without_reset, threshold_without_reset = pulsedetection(recovery_signal, recovery_fs)
    t_d_with_reset, threshold_with_reset = pulsedetection(
        recovery_signal, recovery_fs, reset_after=reset_after
    )

    # Left panels show the interruption and gain change; right panels zoom into
    # recovery. Each column uses matching axes for a direct comparison.
    recovery_figure, recovery_axes = plt.subplots(
        2, 2, sharex="col", sharey="col", figsize=(14, 7)
    )
    conditions = [
        (t_d_without_reset, threshold_without_reset, "Reset disabled"),
        (t_d_with_reset, threshold_with_reset, f"reset_after = {reset_after:g} s"),
    ]
    for row, (detection_times, recovery_threshold, condition_title) in enumerate(conditions):
        for view in range(2):
            limits = (gap_start - 8, gap_end + 15) if view == 0 else (gap_end, gap_end + 12)
            visible = (recovery_time >= limits[0]) & (recovery_time <= limits[1])
            visible_times = detection_times[
                (detection_times >= limits[0]) & (detection_times <= limits[1])
            ]
            axis = recovery_axes[row, view]
            axis.plot(recovery_time[visible], recovery_signal[visible], "k",
                      label="Recorded pulse morphology")
            axis.plot(recovery_time[visible], recovery_threshold[visible], "r",
                      label="Adaptive threshold")
            indices = np.floor(visible_times * recovery_fs + .5).astype(int)
            axis.plot(visible_times, recovery_signal[indices], "bo", label=r"$t_D$")
            axis.set_xlim(limits)
            if view == 0:
                axis.set_title(condition_title + " - interruption and amplitude change")
                axis.axvline(gap_start, color="gray", linestyle="--")
                axis.axvline(gap_end, color="gray", linestyle="--")
                axis.text(gap_start, .95, "Simulated gap", transform=axis.get_xaxis_transform(),
                          va="top", rotation=90)
                axis.text(gap_end, .95, "Amplitude x 0.05", transform=axis.get_xaxis_transform(),
                          va="top", rotation=90)
            else:
                axis.set_title(condition_title + " - recovery detail")
            axis.set_xlabel("Time (s)")
            axis.set_ylabel("Amplitude (a.u.)")
            axis.grid(True)
    recovery_axes[0, 0].legend(loc="best")
    recovery_figure.tight_layout()

    reset_path = None
    if args.save_figure is not None:
        reset_path = args.save_figure.with_name(
            args.save_figure.stem + "_reset" + args.save_figure.suffix
        )
    show_or_save_figure(figure, args.save_figure)
    show_or_save_figure(recovery_figure, reset_path)


if __name__ == "__main__":
    main()
