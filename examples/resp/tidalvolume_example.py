"""Estimate respiratory amplitude and illustrate gap-limited interpolation on a real recording."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, detrend, filtfilt

from biosigpy.resp import tidalvolume


def main():
    fixture = Path(__file__).resolve().parents[1] / "fixtures/ecg/medicom_mtd_ecg_respiration.csv"
    data = np.genfromtxt(fixture, delimiter=",", names=True)
    fs = 256
    t = np.arange(data.size) / fs
    b, a = butter(4, 1, fs=fs)
    resp = filtfilt(b, a, data["respiration"], padlen=3 * (len(b) - 1))

    # Estimate the peak-to-following-valley amplitude of each complete breath.
    tdvol, upper, lower = tidalvolume(resp, fs, mindist=1)
    fig, (signal_ax, volume_ax) = plt.subplots(2, 1, sharex=True, layout="constrained")
    signal_ax.plot(t, detrend(resp), "k", linewidth=1, label="Respiration")
    signal_ax.plot(t, upper, "r", label="Upper envelope")
    signal_ax.plot(t, lower, "b", label="Lower envelope")
    signal_ax.legend()
    signal_ax.set_title("Respiration signal with upper and lower envelopes")
    volume_ax.plot(t, tdvol, "k", linewidth=1)
    volume_ax.set(xlabel="Time (s)", title="Tidal volume")

    # Remove samples after filtering. A missing interval interrupts detection.
    # Each curve applies MaxGap to its own peaks, valleys or cycle midpoints.
    max_gap = 15
    gaps = [(26, 28), (50, 70)]
    missing_resp = resp.copy()
    for start, end in gaps:
        missing_resp[(t >= start) & (t <= end)] = np.nan
    limited, upper, lower = tidalvolume(missing_resp, fs, mindist=1, max_gap=max_gap)

    fig, (signal_ax, volume_ax) = plt.subplots(
        2, 1, sharex=True, figsize=(10, 5.1), layout="constrained"
    )
    # Plot the observed samples on the same detrended scale as the envelopes.
    valid = np.isfinite(missing_resp)
    trend = np.polyfit(t[valid], missing_resp[valid], 1)
    signal_ax.plot(t, missing_resp - np.polyval(trend, t), color="#4d5459", linewidth=1,
                   label="Respiration")
    signal_ax.plot(t, upper, "r", linewidth=1.3, label="Upper envelope")
    signal_ax.plot(t, lower, "b", linewidth=1.3, label="Lower envelope")
    signal_ax.legend(loc="upper left", ncols=3)
    signal_ax.set(ylabel="Respiration (a.u.)", title=f"Respiration and envelopes with MaxGap = {max_gap} s")
    volume_ax.plot(t, limited, color="#146d88", linewidth=1.8)
    volume_ax.set(xlabel="Time (s)", ylabel="Peak-to-valley amplitude (a.u.)", xlim=(10, 88))
    labels = ["Short gap (interpolated)", "Long gap (no interpolation)"]
    for (start, end), label in zip(gaps, labels):
        for ax in (signal_ax, volume_ax):
            ax.axvspan(start, end, color="#b3a68f", alpha=.21, linewidth=0, zorder=0)
        volume_ax.text((start + end) / 2, 1.025, label, transform=volume_ax.get_xaxis_transform(),
                       ha="center", va="bottom", fontsize=9)
    plt.show()


if __name__ == "__main__":
    main()
