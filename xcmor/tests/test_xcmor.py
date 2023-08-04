from ..datasets import reg_ds
from ..xcmor import add_variable_attrs, cmorize
from .tables import mip_amon


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
