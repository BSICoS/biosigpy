# HRV Event Preprocessing

Script: `examples/hrv/preprocessing_example.py`

This synthetic example inserts a false-positive event and applies `removefp`.
It plots the event timestamps and interval series before and after removing the
false-positive detection. The numerical function has no GUI side effects;
plotting belongs only to this example.

Run it from the repository root:

```bash
python examples/hrv/preprocessing_example.py
```

To save the figure instead of opening an interactive window:

```bash
python examples/hrv/preprocessing_example.py --save-figure hrv-preprocessing.png
```
