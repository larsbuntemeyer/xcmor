from ..datasets import reg_ds
from ..xcmor import add_variable_attrs, cmorize

mip_amon = {
    "ta": {
        "frequency": "mon",
        "modeling_realm": "atmos",
        "standard_name": "air_temperature",
        "units": "K",
        "cell_methods": "time: mean",
        "cell_measures": "area: areacella",
        "long_name": "Air Temperature",
        "comment": "Air Temperature",
        "dimensions": "longitude latitude plev19 time",
        "out_name": "ta",
        "type": "real",
        "positive": "",
        "valid_min": "",
        "valid_max": "",
        "ok_min_mean_abs": "",
        "ok_max_mean_abs": "",
    },
    "tas": {
        "frequency": "mon",
        "modeling_realm": "atmos",
        "standard_name": "air_temperature",
        "units": "K",
        "cell_methods": "area: time: mean",
        "cell_measures": "area: areacella",
        "long_name": "Near-Surface Air Temperature",
        "comment": "near-surface (usually, 2 meter) air temperature",
        "dimensions": "longitude latitude time height2m",
        "out_name": "tas",
        "type": "real",
        "positive": "",
        "valid_min": "",
        "valid_max": "",
        "ok_min_mean_abs": "",
        "ok_max_mean_abs": "",
    },
    "pr": {
        "frequency": "mon",
        "modeling_realm": "atmos",
        "standard_name": "precipitation_flux",
        "units": "kg m-2 s-1",
        "cell_methods": "area: time: mean",
        "cell_measures": "area: areacella",
        "long_name": "Precipitation",
        "comment": "includes both liquid and solid phases",
        "dimensions": "longitude latitude time",
        "out_name": "pr",
        "type": "real",
        "positive": "",
        "valid_min": "",
        "valid_max": "",
        "ok_min_mean_abs": "",
        "ok_max_mean_abs": "",
    },
}


def test_add_variable_attrs():
    mip_table = mip_amon
    ds = reg_ds[["temperature"]].rename(temperature="ta")
    result = add_variable_attrs(ds, mip_table)
    assert result.ta.units == mip_table["ta"]["units"]
    assert result.ta.standard_name == mip_table["ta"]["standard_name"]
    assert result.frequency == mip_table["ta"]["frequency"]


def test_cmorize():
    ds = reg_ds.copy()
    mip_table = mip_amon
    ds_out = cmorize(
        ds.rename({"temperature": "ta", "precipitation": "pr"}), mip_table=mip_table
    )

    expected_var_attrs = [
        "standard_name",
        "units",
        "cell_methods",
        "cell_measures",
        "long_name",
    ]
    expected_global_attrs = ["frequency"]

    for var in ds_out.data_vars:
        da = ds_out[var]
        for k in expected_var_attrs:
            assert da.attrs[k] == mip_table[var][k]
        for k in expected_global_attrs:
            assert ds_out.attrs[k] == mip_table[var][k]

    return ds_out
