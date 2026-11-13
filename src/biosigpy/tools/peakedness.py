"""Measure concentration and relative peak height in spectral bands."""

import numpy as np
from numpy.typing import ArrayLike

from biosigpy.tools._validation import as_positive_real_scalar, as_real_vector


def peakedness(pxx: ArrayLike, f: ArrayLike, reference_freq: float | None = None,
               window: float = 0.125) -> tuple[float | np.ndarray, float | np.ndarray]:
    """Return power concentration (PKL) and relative maximum (AKL) in percent.

    Parameters
    ----------
    pxx : array_like
        Nonempty, real, finite, nonnegative spectrum or frequency-by-spectrum
        matrix. Row and column vectors are treated as a single spectrum.
    f : array_like
        Real, finite, strictly increasing frequency vector in Hz, matching pxx.
    reference_freq : float or None, optional
        Center frequency in Hz. None selects each spectrum's global maximum;
        ties select the first frequency. No fixed frequency range is applied.
    window : float, optional
        Finite positive width of the wide band in Hz, default 0.125. The narrow
        band has 40 percent of this width. Band endpoints are included.

    Returns
    -------
    pkl : float or ndarray
        100 times the sum in the narrow band divided by the sum in the wide
        band.
    akl : float or ndarray
        100 times the maximum in the wide band divided by the global maximum.
        In automatic mode this is 100 for every nonzero spectrum. Zero spectra
        and bands without samples return zero for the corresponding measure.
        Vector inputs return scalars; matrices return one value per column.
    """
    # Parse and validate inputs
    values = np.asarray(pxx)
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("pxx must be real numeric data")
    if values.ndim > 2 or values.size == 0:
        raise ValueError("pxx must be a nonempty spectrum or matrix")
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("pxx must be finite and nonnegative")
    values = np.asarray(values, dtype=float)
    frequencies = as_real_vector(f, name="f")
    if frequencies.size == 0 or not np.all(np.isfinite(frequencies)) or np.any(np.diff(frequencies) <= 0):
        raise ValueError("f must be finite, nonempty and strictly increasing")
    width = as_positive_real_scalar(window, name="window")
    use_adaptive_method = reference_freq is None
    if not use_adaptive_method:
        reference = np.asarray(reference_freq)
        if not np.issubdtype(reference.dtype, np.number) or np.iscomplexobj(reference):
            raise TypeError("reference_freq must be real numeric data")
        if reference.ndim != 0 or not np.isfinite(reference):
            raise ValueError("reference_freq must be a finite scalar or None")
        reference_freq = float(reference)

    # Ensure pxx is a matrix with spectra as columns
    vector_input = values.ndim < 2 or 1 in values.shape
    if vector_input:
        values = values.reshape(-1, 1)
    if values.shape[0] != frequencies.size:
        raise ValueError("first dimension of pxx must match length of f")

    # Pre-allocate output arrays
    num_spectra = values.shape[1]
    pkl = np.zeros(num_spectra)
    akl = np.zeros(num_spectra)

    # Calculate peakedness for each spectrum
    for i in range(num_spectra):
        current_pxx = values[:, i]
        if np.all(current_pxx == 0):
            continue

        # Determine center frequency for window definition
        if use_adaptive_method:
            center_freq = frequencies[np.argmax(current_pxx)]
        else:
            center_freq = reference_freq
        window_narrow = 0.4 * width

        # Define search windows around center frequency
        search_window_wide = ((frequencies >= center_freq - width / 2)
                              & (frequencies <= center_freq + width / 2))
        search_window_narrow = ((frequencies >= center_freq - window_narrow / 2)
                                & (frequencies <= center_freq + window_narrow / 2))
        if not np.any(search_window_wide):
            continue

        # Calculate power concentration peakedness
        power_in_narrow_window = np.sum(current_pxx[search_window_narrow])
        power_in_wide_window = np.sum(current_pxx[search_window_wide])
        if power_in_wide_window > 0:
            pkl[i] = 100 * power_in_narrow_window / power_in_wide_window

        # Calculate absolute maximum peakedness
        max_in_window = np.max(current_pxx[search_window_wide])
        global_max = np.max(current_pxx)
        if global_max > 0:
            akl[i] = 100 * max_in_window / global_max

        # Ensure peakedness values are within 0-100% range
        pkl[i] = np.clip(pkl[i], 0, 100)
        akl[i] = np.clip(akl[i], 0, 100)

    if vector_input:
        return float(pkl[0]), float(akl[0])
    return pkl, akl
