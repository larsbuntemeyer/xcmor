from datetime import date

from xarray import DataArray

# from .utils import filter_table_by_value


def cmorize(ds, mip_table=None, coords_table=None, dataset_table=None, mapping=None):
    ds = ds.copy()

    if mapping is None:
        mapping = {}

    if isinstance(ds, DataArray):
        ds = ds.to_dataset()

    ds = ds.rename({v: mapping.get(v) for v in ds})

    ds = add_variable_attrs(ds, mip_table["variable_entry"])

    for var in ds.data_vars:
        da = apply_dimensions(ds[var], coords_table)
        ds[var] = da
        ds = ds.assign_coords(da.coords)

    ds = add_global_attributes(ds, dataset_table)

    ds = add_version_attribute(ds)

    return ds


def add_variable_attrs(ds, mip_table):
    """add variable attributes"""

    for v in ds.data_vars:
        ds[v].attrs = mip_table[v]

    return ds


def apply_dimensions(da, coords_table):
    """apply dimensions from coordinates table"""
    da = da.copy()

    coords_table = coords_table.get("axis_entry") or coords_table

    dims = da.attrs.get("dimensions")
    print(dims)

    if dims:
        dims = {d: coords_table[d] for d in dims.split()}

    for d, v in dims.items():
        keys = ["out_name", "standard_name"]
        out_name = v["out_name"]
        for k in keys:
            if v[k] in da.coords:
                da.coords[v[k]].attrs = v
                da = da.rename({v[k]: out_name})
                break

    return da


def add_version_attribute(ds):
    """add version attribute"""
    now = date.today().strftime("%Y%m%d")
    ds.attrs["version"] = now

    return ds


def add_global_attributes(ds, dataset_table):
    ds.attrs = {k: v for k, v in dataset_table.items() if not k.startswith("#")}

    return ds
