"""Low-pass differentiator FIR filter design."""

import numpy as np
from scipy import special

from biosigpy.tools._validation import as_integer_scalar, as_positive_real_scalar


def lpd_filter(
    sampling_frequency: float,
    stop_frequency: float,
    pass_frequency: float | None = None,
    order: int | None = None,
) -> tuple[np.ndarray, float]:
    """Design a Biosiglib-compatible low-pass differentiator FIR filter.

    Parameters
    ----------
    sampling_frequency : float
        Positive sampling frequency in Hz.
    stop_frequency : float
        Positive stop frequency in Hz. The value must be less than the
        Nyquist frequency.
    pass_frequency : float or None, optional
        Positive pass frequency in Hz. When omitted, ``stop_frequency - 0.2``
        is used. The pass frequency must be less than ``stop_frequency``.
    order : int or None, optional
        Explicit positive filter order. Automatic order selection is currently
        unsupported.

    Returns
    -------
    tuple[numpy.ndarray, float]
        ``(b, delay)`` where ``b`` is a one-dimensional array of FIR
        coefficients and ``delay`` is the filter delay in samples.

    Raises
    ------
    TypeError
        If scalar inputs are not real numeric values.
    ValueError
        If frequencies or order are outside their accepted ranges.
    NotImplementedError
        If ``order`` is ``None``.

    Notes
    -----
    This function implements the Biosiglib ``tools.lpd_filter`` specification
    and follows the MATLAB/Biosigmat-compatible low-pass differentiator design.

    Explicit ``order`` design is supported. Automatic order selection is not
    currently implemented. The effective order is rounded up to an even value.
    Coefficients are scaled by ``sampling_frequency / (2*pi)``, and the
    returned delay is ``filter_order / 2`` for the effective filter order.

    Examples
    --------
    >>> from biosigpy.tools import lpd_filter
    >>> b, delay = lpd_filter(256.0, 12.0, order=4)
    >>> b.shape
    (5,)
    >>> delay
    2.0
    """

    sampling_frequency = as_positive_real_scalar(
        sampling_frequency, name="sampling_frequency"
    )
    stop_frequency = as_positive_real_scalar(stop_frequency, name="stop_frequency")

    if pass_frequency is None:
        pass_frequency = stop_frequency - 0.2
    pass_frequency = as_positive_real_scalar(pass_frequency, name="pass_frequency")

    if stop_frequency >= sampling_frequency / 2.0:
        raise ValueError("stop_frequency must be less than Nyquist frequency")
    if pass_frequency >= stop_frequency:
        raise ValueError("pass_frequency must be less than stop_frequency")

    if order is None:
        raise NotImplementedError("automatic lpd_filter order is not supported")
    filter_order = as_integer_scalar(order, name="order")
    if filter_order <= 0:
        raise ValueError("order must be positive")
    filter_order = filter_order + (filter_order % 2)

    w_pass = pass_frequency / (sampling_frequency / 2.0)
    w_stop = stop_frequency / (sampling_frequency / 2.0)
    b = _firls_differentiator(filter_order, w_pass, w_stop)
    b = b * sampling_frequency / (2.0 * np.pi)
    delay = filter_order / 2.0

    return b, delay


def _firls_differentiator(
    filter_order: int, w_pass: float, w_stop: float
) -> np.ndarray:
    # Equivalent to MATLAB:
    # firls(filter_order, [0 w_pass w_stop 1], [0 pi*w_pass 0 0],
    #       [1 1], 'differentiator')
    freq = np.array([0.0, w_pass, w_stop, 1.0], dtype=np.float64)
    amp = np.array([0.0, np.pi * w_pass, 0.0, 0.0], dtype=np.float64)
    weight = np.ones(2, dtype=np.float64)

    frequency = freq / 2.0
    wt = np.sqrt(weight)
    filter_length = filter_order + 1
    half_order = (filter_length - 1) // 2
    m = np.arange(1, half_order + 1, dtype=np.float64)
    k = m.reshape(-1, 1)
    i1, i2 = _init_matrices(m)

    rhs = np.zeros((m.size, 1), dtype=np.float64)
    gram = np.zeros((m.size, m.size), dtype=np.float64)
    do_weight = (
        np.abs(amp[0::2]) + np.abs(amp[1::2])
    ) > 0.0

    for band_index, band_start in enumerate(range(0, frequency.size, 2)):
        f_start = float(frequency[band_start])
        f_end = float(frequency[band_start + 1])
        amp_start = float(amp[band_start])
        amp_end = float(amp[band_start + 1])
        weight_squared = wt[band_index] ** 2

        if do_weight[band_index]:
            if f_start == 0.0:
                f_start = 1e-5
            slope = (amp_end - amp_start) / (f_end - f_start)
            intercept = amp_start - slope * f_start
            rhs += _weighted_rhs(
                k, f_start, f_end, slope, intercept
            ) * weight_squared
            gram += _weighted_gram(i1, i2, f_start, f_end) * weight_squared
        else:
            slope = (amp_end - amp_start) / (f_end - f_start)
            intercept = amp_start - slope * f_start
            rhs += _unweighted_rhs(
                k, f_start, f_end, slope, intercept
            ) * weight_squared
            gram += _unweighted_gram(i1, i2, f_start, f_end) * weight_squared

    a = np.linalg.solve(gram, rhs).ravel()
    return -0.5 * np.concatenate((a[::-1], np.array([0.0]), -a))


def _weighted_rhs(
    k: np.ndarray, f_start: float, f_end: float, slope: float, intercept: float
) -> np.ndarray:
    snint = _sineint(2.0 * np.pi * k * f_end) - _sineint(
        2.0 * np.pi * k * f_start
    )
    csint = _cosint(2.0 * np.pi * k * f_end) - _cosint(
        2.0 * np.pi * k * f_start
    )
    return slope * snint + intercept * 2.0 * np.pi * k * (
        -_sinc(2.0 * k * f_end) + _sinc(2.0 * k * f_start) + csint
    )


def _weighted_gram(
    i1: np.ndarray, i2: np.ndarray, f_start: float, f_end: float
) -> np.ndarray:
    term_end = -0.5 * (
        np.cos(2.0 * np.pi * f_end * (-i2)) / f_end
        - 2.0 * _sineint(2.0 * np.pi * f_end * (-i2)) * np.pi * i2
        - np.cos(2.0 * np.pi * f_end * i1) / f_end
        - 2.0 * _sineint(2.0 * np.pi * f_end * i1) * np.pi * i1
    )
    term_start = -0.5 * (
        np.cos(2.0 * np.pi * f_start * (-i2)) / f_start
        - 2.0 * _sineint(2.0 * np.pi * f_start * (-i2)) * np.pi * i2
        - np.cos(2.0 * np.pi * f_start * i1) / f_start
        - 2.0 * _sineint(2.0 * np.pi * f_start * i1) * np.pi * i1
    )
    return -(term_end - term_start)


def _unweighted_rhs(
    k: np.ndarray, f_start: float, f_end: float, slope: float, intercept: float
) -> np.ndarray:
    rhs = slope / (4.0 * np.pi**2.0) * (
        np.sin(2.0 * np.pi * k * f_end)
        - np.sin(2.0 * np.pi * k * f_start)
    ) / (k * k)
    rhs += (
        (slope * f_start + intercept) * np.cos(2.0 * np.pi * k * f_start)
        - (slope * f_end + intercept) * np.cos(2.0 * np.pi * k * f_end)
    ) / (2.0 * np.pi * k)
    return rhs


def _unweighted_gram(
    i1: np.ndarray, i2: np.ndarray, f_start: float, f_end: float
) -> np.ndarray:
    return 0.5 * f_end * (
        _sinc(2.0 * i1 * f_end) - _sinc(2.0 * i2 * f_end)
    ) - 0.5 * f_start * (
        _sinc(2.0 * i1 * f_start) - _sinc(2.0 * i2 * f_start)
    )


def _init_matrices(m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    k = m.reshape(-1, 1)
    x = np.tile(k, (1, m.size))
    y = np.tile(m, (m.size, 1))
    return x + y, x - y


def _sineint(x: np.ndarray) -> np.ndarray:
    return special.sici(x)[0]


def _cosint(x: np.ndarray) -> np.ndarray:
    return special.sici(x)[1]


def _sinc(x: np.ndarray) -> np.ndarray:
    return np.sinc(x)
