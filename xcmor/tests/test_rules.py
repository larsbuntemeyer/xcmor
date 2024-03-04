import pytest
import xarray as xr

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

    da.attrs = {
        "standard_name": "air_temp",
    }
    with pytest.raises(Exception) as e_info:
        rules.standard_name(da)
        assert e_info == f"{da.standard_name} is not a valid standard name"
