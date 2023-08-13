import numpy as np

from ..datasets import plev_ds, reg_ds
from ..xcmor import Cmorizer, _add_var_attrs, cmorize
from .tables import coords, dataset, mip_amon


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
