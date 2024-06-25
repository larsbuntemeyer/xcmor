import numpy as np
import xarray as xr
from cf_xarray.datasets import rotds

from ..datasets import plev_ds, reg_ds, temp_ds
from ..mapping import dtype_map
from ..xcmor import (
    Cmorizer,
    _add_var_attrs,
    _get_lon_lat_coords,
    _get_x_y_coords,
    _guess_dims_attr,
    #    _interpret_var_dims,
    _is_curvilinear,
    _transpose,
    _units_convert,
    cmorize,
)
from . import requires_pint_xarray
from .tables import coords, dataset, mip_amon

expected_var_attrs = [
    "standard_name",
    "units",
    "cell_methods",
    "cell_measures",
    "long_name",
]

# def test_encode_time():

#     expected = "days since 2014-09-06T00:00:00"
#     time_out = _encode_time(reg_ds, cf_units="days since ?")
#     assert time_out.encoding['units'] == expected


@requires_pint_xarray
def test_units_convert():
    ds = reg_ds.copy()
    da = ds.temperature
    da.attrs["original_units"] = "degC"
    da.attrs["units"] = "K"
    da_conv = _units_convert(da)
    expected_data = np.array(
        [
            [
                [302.26241877, 291.35125767, 295.97990387],
                [306.07714559, 303.09046392, 280.33177696],
            ],
            [
                [295.75070734, 286.93914233, 287.32424919],
                [291.43478802, 289.30234857, 299.78418806],
            ],
        ]
    )
    np.testing.assert_allclose(da_conv, expected_data)
    assert "original data with units degC converted to K" in da_conv.history


def test_transpose():
    ds_out = _transpose(reg_ds)
    assert list(ds_out.temperature.dims) == ["time", "y", "x"]


def test_guess_dims_attr():
    assert _guess_dims_attr(plev_ds) == ["longitude", "latitude", "lev", "time"]
    assert _guess_dims_attr(reg_ds) == ["longitude", "latitude", "time"]
    assert _guess_dims_attr(temp_ds) == ["longitude", "latitude", "time"]


# def test_interpret_var_dims():
#    ds = xr.Dataset(
#        data_vars=dict(
#            temp=(["x"], [1, 4, 2, 9], {"dimensions": "longitude"}),
#        ),
#        coords=dict(
#            x=(["x"], [1, 2, 3, 4]),
#        ),
#        attrs=dict(description="Weather related data."),
#    )
#    ds_out = _interpret_var_dims(
#        ds.cf.guess_coord_axis(verbose=True), coords["axis_entry"]
#    )
#    assert coords["axis_entry"]["longitude"]["out_name"] in ds_out.coords
#
#    ds.temp.attrs = {"dimensions": "longitude height2m"}
#    ds_out = _interpret_var_dims(
#        ds.cf.guess_coord_axis(verbose=True), coords["axis_entry"]
#    )
#    assert coords["axis_entry"]["longitude"]["out_name"] in ds_out.coords
#    assert coords["axis_entry"]["height2m"]["out_name"] in ds_out.coords


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


def test_get_lon_lat_coords():
    result = _get_lon_lat_coords(rotds)
    xr.testing.assert_equal(result[0], rotds.cf["longitude"])
    xr.testing.assert_equal(result[1], rotds.cf["latitude"])

    result = _get_lon_lat_coords(reg_ds)
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


def test_cmorize_auto():
    ds = reg_ds.copy()
    ds_out = cmorize(
        ds.rename({"temperature": "tas", "precipitation": "pr"}),
    )

    return ds_out


def test_cmorize_minimal():
    ds = reg_ds.copy()
    mip_table = mip_amon
    mapping_table = {"temperature": "tas", "precipitation": "pr"}
    ds_out = cmorize(
        ds.rename(mapping_table),
        mip_table=mip_table,
    )

    expected_global_attrs = ["frequency"]

    for var in ds.data_vars:
        cf_name = mapping_table.get(var)
        assert cf_name in ds_out
        cf_var = ds_out[cf_name]
        # ensure cmorizer does not change values
        # np.testing.assert_allclose(cf_var, ds[var])
        # test for expected attributes
        for k in expected_var_attrs:
            assert cf_var.attrs[k] == mip_table[cf_name][k]
        for k in expected_global_attrs:
            assert ds_out.attrs[k] == mip_table[cf_name][k]


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
        # np.testing.assert_allclose(da, ds["temperature"])
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


def test_cmorizer_cmorize_cmip6():
    ds = reg_ds.rename(temperature="tas")
    cmor = Cmorizer(project="CMIP6")
    mip_table = "Amon"
    ds_out = cmor.cmorize(ds.tas, mip_table, cmor.tables["input_example"])

    expected_global_attrs = ["frequency"]

    for var in ds_out.data_vars:
        da = ds_out[var]
        # ensure cmorizer does not change values
        # np.testing.assert_allclose(da, ds[var])
        # test for expected attributes
        for k in expected_var_attrs:
            assert da.attrs[k] == cmor.tables[mip_table]["variable_entry"][var][k]
        for k in expected_global_attrs:
            assert ds_out.attrs[k] == cmor.tables[mip_table]["variable_entry"][var][k]

    assert ds_out.to_netcdf("test.nc") is None
