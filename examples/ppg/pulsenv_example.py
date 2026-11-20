"""Estimate pulse envelopes and compare short and long gaps on a real PPG."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import lfilter

from biosigpy.ppg import pulsedel, pulsedet, pulsenv
from biosigpy.tools import lpdfilter


def main():
    # Load 20 seconds of PPG and detect pulses on the aligned derivative.
    fixture = Path(__file__).resolve().parents[1] / "fixtures/ppg/ppg_signals.csv"
    data = np.genfromtxt(fixture, delimiter=",", names=True, max_rows=20000)
    t, ppg = data["t"], data["sig"]
    fs = 1000
    b, delay = lpdfilter(fs, 8, pass_frequency=7.8, order=100)
    delay = int(delay)
    dppg = lfilter(b, [1], ppg)
    dppg = np.r_[dppg[delay:], np.zeros(delay)]
    t_d, _ = pulsedet(dppg, fs)

    lower, upper = pulsenv(ppg, fs, t_d)
    fig, ax = plt.subplots()
    ax.plot(t, ppg, "k", linewidth=1, label="PPG")
    ax.plot(t, lower, "b", label="Lower envelope")
    ax.plot(t, upper, "r", label="Upper envelope")
    samples = np.floor(t_d * fs + .5).astype(int)
    ax.plot(t_d, ppg[samples], "o", color=(.47, .67, .19), label=r"$t_D$")
    ax.set(xlabel="Time (s)", ylabel="Amplitude", title="PPG envelopes")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()

    # Remove observed samples after detection to demonstrate optional MaxGap.
    # The shaded intervals mark missing PPG. MaxGap is measured separately
    # between consecutive basal points and between consecutive apex points.
    max_gap = 2.0
    gaps = [(t_d[5] - .22, t_d[5] + .24), (t_d[11] - .22, t_d[14] + .24)]
    missing_ppg = ppg.copy()
    for start, end in gaps:
        missing_ppg[(t >= start) & (t <= end)] = np.nan
    t_a, t_b = pulsedel(missing_ppg, fs, t_d, points="AB")
    lower, upper = pulsenv(missing_ppg, fs, t_a=t_a, t_b=t_b, max_gap=max_gap)

    fig, ax = plt.subplots(figsize=(10, 3.7), layout="constrained")
    ax.plot(t, missing_ppg, color="#62696d", linewidth=1, label="PPG")
    for times, envelope, color, label in [
        (t_b, lower, "#146d88", "Lower envelope"),
        (t_a, upper, "#ad441f", "Upper envelope"),
    ]:
        ax.plot(t, envelope, color=color, linewidth=1.7, label=label)
        anchors = np.floor(times[np.isfinite(times)] * fs + .5).astype(int)
        ax.scatter(t[anchors], missing_ppg[anchors], s=14, color=color, zorder=3)
    for (start, end), label in zip(gaps, ["Short gap (interpolated)", "Long gap (no interpolation)"]):
        ax.axvspan(start, end, color="#b3a68f", alpha=.21, linewidth=0)
        ax.text((start + end) / 2, 1.025, label, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9)
    ax.set(xlabel="Time (s)", ylabel="PPG (a.u.)", xlim=(t_d[2] - .16, t_d[18] + .4))
    ax.spines[["right", "top"]].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.22), ncols=3, frameon=False)
    fig.suptitle("PPG envelopes with MaxGap = 2 s")
    plt.show()


if __name__ == "__main__":
    main()
