import pytest
import xarray as xr

from ..mapping import dtype_map
from ..rules import rules


def test_rules():
    varname = "foo"
    da = xr.DataArray(
        data=[1, 4, 2, 9],
        dims=("x"),
        coords={"x": [0, 1, 2, 3]},
    )
    da.attrs["type"] = "float64"
    da.name = varname
    da.attrs = {
        "coordinates": "longitude",
        "type": "real",
        "out_name": "tas",
        "valid_min": 0,
        "valid_max": 10,
        "standard_name": "air_temperature",
    }
    ds = da.to_dataset()
    for attr in ds[varname].attrs:
        if hasattr(rules, attr):
            ds = getattr(rules, attr)(ds, varname)

    assert ds[varname].name == "foo"
    assert ds[varname].dtype == dtype_map["real"]

    ds[varname].attrs = {
        "standard_name": "air_temp",
    }
    with pytest.raises(Exception) as e_info:
        rules.standard_name(ds, varname)
        assert e_info == f"{ds[varname].standard_name} is not a valid standard name"

    ds[varname].attrs = {
        "valid_min": 0,
    }
    ds[varname][:] = [2.0, 4.0, 2.0, 9.0]
    with pytest.raises(Exception) as e_info:
        rules.valid_min(ds, varname)
        assert (
            e_info
            == f"{ds[varname].name} is violating valid_min: {ds[varname].valid_min}"
        )
