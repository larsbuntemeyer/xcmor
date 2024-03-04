import pytest
import xarray as xr

from ..mapping import dtype_map
from ..rules import rules


def test_rules():
    da = xr.DataArray(
        data=[1, 4, 2, 9],
        dims=("x"),
        coords={"x": [0, 1, 2, 3]},
    )
    da.attrs["type"] = "float64"
    da.name = "variable"
    da.attrs = {
        "coordinates": "longitude",
        "type": "real",
        "out_name": "tas",
        "valid_min": 0,
        "valid_max": 10,
        "standard_name": "air_temperature",
    }
    for attr in da.attrs:
        if hasattr(rules, attr):
            da = getattr(rules, attr)(da)

    assert da.name == "tas"
    assert da.dtype == dtype_map["real"]

    da.attrs = {
        "standard_name": "air_temp",
    }
    with pytest.raises(Exception) as e_info:
        rules.standard_name(da)
        assert e_info == f"{da.standard_name} is not a valid standard name"

    da.attrs = {
        "valid_min": 0,
    }
    da[:] = [2.0, 4.0, 2.0, 9.0]
    with pytest.raises(Exception) as e_info:
        rules.valid_min(da)
        assert e_info == f"{da.name} is violating valid_min: {da.valid_min}"
