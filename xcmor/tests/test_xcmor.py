import numpy as np
import xarray as xr
from cf_xarray.datasets import rotds

from ..datasets import plev_ds, reg_ds
from ..xcmor import (
    Cmorizer,
    _add_var_attrs,
    _get_lon_lat,
    _get_x_y_coords,
    _is_curvilinear,
    cmorize,
)
from .tables import coords, dataset, mip_amon

expected_var_attrs = [
    "standard_name",
    "units",
    "cell_methods",
    "cell_measures",
    "long_name",
]


def test_get_x_y_coords():
    result = _get_x_y_coords(rotds)
    xr.testing.assert_equal(result[0], rotds.rlon)
    xr.testing.assert_equal(result[1], rotds.rlat)

    result = _get_x_y_coords(reg_ds)
    expect_lon = reg_ds.cf["lon"]
    expect_lat = reg_ds.cf["lat"]
    expect_lon.attrs = {
        "units": "degrees_east",
        "standard_name": "longitude",
        "axis": "X",
    }
    expect_lat.attrs = {
        "units": "degrees_north",
        "standard_name": "latitude",
        "axis": "Y",
    }
    xr.testing.assert_equal(result[0], expect_lon)
    xr.testing.assert_equal(result[1], expect_lat)


def test_get_lon_lat():
    result = _get_lon_lat(rotds)
    xr.testing.assert_equal(result[0], rotds.cf["longitude"])
    xr.testing.assert_equal(result[1], rotds.cf["latitude"])

    result = _get_lon_lat(reg_ds)
    xr.testing.assert_equal(result[0], reg_ds.cf["lon"])
    xr.testing.assert_equal(result[1], reg_ds.cf["lat"])


def test_curvilinear():
    assert _is_curvilinear(rotds) is True
    assert _is_curvilinear(reg_ds) is False


def test_add_variable_attrs():
    mip_table = mip_amon
    ds = reg_ds[["temperature"]].rename(temperature="ta")
    result = _add_var_attrs(ds, mip_table)
    assert result.ta.units == mip_table["ta"]["units"]
    assert result.ta.standard_name == mip_table["ta"]["standard_name"]
    assert result.frequency == mip_table["ta"]["frequency"]


def test_cmorize_minimal():
    ds = reg_ds.copy()
    mip_table = mip_amon
    ds_out = cmorize(
        ds.rename({"temperature": "ta", "precipitation": "pr"}),
        mip_table=mip_table,
    )

    return ds_out

    expected_global_attrs = ["frequency"]

    for var in ds_out.data_vars:
        da = ds_out[var]
        # ensure cmorizer does not change values
        np.testing.assert_allclose(da, ds[var])
        # test for expected attributes
        for k in expected_var_attrs:
            assert da.attrs[k] == mip_table[var][k]
        for k in expected_global_attrs:
            assert ds_out.attrs[k] == mip_table[var][k]


def test_cmorize():
    ds = plev_ds
    mip_table = mip_amon
    ds_out = cmorize(
        ds.rename({"temperature": "ta"}),
        mip_table=mip_table,
        coords_table=coords,
        dataset_table=dataset,
    )

    expected_global_attrs = ["frequency"]

    for var in ds_out.data_vars:
        da = ds_out[var]
        # ensure cmorizer does not change values
        np.testing.assert_allclose(da, ds["temperature"])
        # test for expected attributes
        for k in expected_var_attrs:
            assert da.attrs[k] == mip_table[var][k]
        for k in expected_global_attrs:
            assert ds_out.attrs[k] == mip_table[var][k]


def test_cmorizer():
    cmorizer = Cmorizer()
    assert cmorizer.project == "CMIP6"
    assert (
        cmorizer.tables.get_url("Amon")
        == "https://raw.githubusercontent.com/PCMDI/cmip6-cmor-tables/master/Tables/CMIP6_Amon.json"
    )

    cmorizer = Cmorizer(project="CORDEX")
    assert cmorizer.project == "CORDEX"
    assert (
        cmorizer.tables.get_url("mon")
        == "https://raw.githubusercontent.com/WCRP-CORDEX/cordex-cmip6-cmor-tables/main/Tables/CORDEX_mon.json"
    )


def test_cmorizer_cmorize():
    ds = reg_ds.rename(temperature="tas")
    cmor = Cmorizer(project="CMIP6")
    mip_table = "Amon"
    ds_out = cmor.cmorize(ds.tas, mip_table, cmor.tables["input_example"])

    expected_global_attrs = ["frequency"]

    for var in ds_out.data_vars:
        da = ds_out[var]
        # ensure cmorizer does not change values
        np.testing.assert_allclose(da, ds[var])
        # test for expected attributes
        for k in expected_var_attrs:
            assert da.attrs[k] == cmor.tables[mip_table]["variable_entry"][var][k]
        for k in expected_global_attrs:
            assert ds_out.attrs[k] == cmor.tables[mip_table]["variable_entry"][var][k]

    assert ds_out.to_netcdf("test.nc") is None
