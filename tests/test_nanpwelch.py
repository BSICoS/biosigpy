"""Shared PSD normalization, gap handling and channel alignment."""

import warnings
import numpy as np
import pytest
from biosigpy.tools import nanpwelch, NanPwelchWarning
from conformance import assert_expected_outputs, assert_expected_warnings, case_id, cases_for_specification, load_input


@pytest.mark.parametrize("case_definition",cases_for_specification("tools.nanpwelch"),ids=case_id)
def test_conformance(case_definition):
    args=[load_input(case_definition,key) for key in
          ("signal","window","overlap","nfft","sampling_frequency")]
    # JSON NaN markers are strings; numerical fixtures already contain NaNs.
    args[0] = np.asarray(args[0], dtype=float)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        psd,f,windows=nanpwelch(*args,**case_definition["parameters"])
    actual_warnings=[item.message for item in captured if isinstance(item.message,NanPwelchWarning)]
    assert len(actual_warnings)==len(captured)
    assert_expected_warnings(actual_warnings,case_definition)
    if isinstance(windows,list):
        flat=np.concatenate([item.reshape(-1) for item in windows])
    else:
        flat=windows.reshape(-1)
    assert_expected_outputs(dict(psd=psd.reshape(-1,order="F"),frequencies=f,window_psd=flat),case_definition)


@pytest.mark.parametrize("overrides",[{"window":0},{"window":[0,0]}, {"overlap":4},
    {"nfft":0}, {"nfft":2.5}, {"sampling_frequency":np.inf}, {"max_gap":-1}])
def test_invalid_parameters(overrides):
    args=dict(signal=np.ones(8),window=4,overlap=0,nfft=4,sampling_frequency=4)
    args.update(overrides)
    with pytest.raises((ValueError,TypeError)):
        nanpwelch(**args)



def test_matrix_channel_collections_and_optional_gap():
    signal=np.column_stack((np.ones(8),np.full(8,np.nan),np.full(8,2)))
    for max_gap in (None, [], 0):
        psd,f,windows=nanpwelch(signal,4,0,4,4,max_gap)
        assert psd.shape==(3,3) and len(windows)==3
        assert windows[0].shape==(2,3) and windows[1].size==0
        assert windows[2].shape==(2,3)
        np.testing.assert_array_equal(psd[:,1],np.full(3,np.nan))
