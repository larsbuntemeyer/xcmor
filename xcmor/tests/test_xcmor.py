import numpy as np
import xarray as xr
from cf_xarray.datasets import rotds

from ..datasets import plev_ds, reg_ds
from ..mapping import dtype_map
from ..xcmor import (
    Cmorizer,
    _add_var_attrs,
    _get_lon_lat,
    _get_x_y_coords,
    _interpret_var_dims,
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


def test_interpret_var_dims():
    ds = xr.Dataset(
        data_vars=dict(
            temp=(["x"], [1, 4, 2, 9], {"dimensions": "longitude"}),
        ),
        coords=dict(
            x=(["x"], [1, 2, 3, 4]),
        ),
        attrs=dict(description="Weather related data."),
    )
    ds_out = _interpret_var_dims(
        ds.cf.guess_coord_axis(verbose=True), coords["axis_entry"]
    )
    assert coords["axis_entry"]["longitude"]["out_name"] in ds_out.coords

    ds.temp.attrs = {"dimensions": "longitude height2m"}
    ds_out = _interpret_var_dims(
        ds.cf.guess_coord_axis(verbose=True), coords["axis_entry"]
    )
    assert coords["axis_entry"]["longitude"]["out_name"] in ds_out.coords
    assert coords["axis_entry"]["height2m"]["out_name"] in ds_out.coords


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
        # ensure cmorizer does not change values but datatype
        np.testing.assert_allclose(da, ds["temperature"])
        # test for expected attributes
        for k in expected_var_attrs:
            assert da.attrs[k] == mip_table[var][k]
        for k in expected_global_attrs:
            assert ds_out.attrs[k] == mip_table[var][k]

    # test cmorization of dataset wiht two variables and a mapping table
    ds_out = cmorize(
        reg_ds,
        mip_table=mip_amon,
        coords_table=coords,
        dataset_table=dataset,
        mapping_table={"temperature": "tas", "precipitation": "pr"},
    )

    assert "tas" in ds_out
    assert "pr" in ds_out

    for var in ds_out:
        da = ds_out[var]
        assert da.dtype == dtype_map["real"]
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

    cmorizer = Cmorizer(project="CORDEX-CMIP6")
    assert cmorizer.project == "CORDEX-CMIP6"
    assert (
        cmorizer.tables.get_url("mon")
        == "https://raw.githubusercontent.com/WCRP-CORDEX/cordex-cmip6-cmor-tables/main/Tables/CORDEX-CMIP6_mon.json"
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
