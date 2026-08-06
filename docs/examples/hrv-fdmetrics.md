# HRV Frequency-domain Metrics

Script: `examples/hrv/fdmetrics_example.py`

This example reproduces the Biosigmat workflow with the local Medicom MTD
fixtures. It reconstructs the TVIPFM modulating signal, aligns detrended
respiration to the same 4 Hz grid, and estimates their spectra. It reports
conventional metrics with limited and unlimited HF bands, then applies OSP
and reports the respiration-related frequency-domain metrics.

Run it from the repository root:

```bash
python examples/hrv/fdmetrics_example.py
```
