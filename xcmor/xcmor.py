import collections
from datetime import date
from warnings import warn

from xarray import DataArray

# from .utils import filter_table_by_value


def cmorize(
    ds,
    mip_table=None,
    coords_table=None,
    dataset_table=None,
    cv_table=None,
    mapping=None,
):
    """Lazy cmorization"""

    ds = ds.copy()

    coords_table = coords_table.get("axis_entry") or coords_table

    if mapping is None:
        mapping = {}

    if isinstance(ds, DataArray):
        ds = ds.to_dataset()

    ds = ds.rename({v: mapping.get(v) for v in ds})

    ds = add_variable_attrs(ds, mip_table["variable_entry"] or mip_table)

    for var in ds.data_vars:
        dims = ds[var].attrs.get("dimensions")
        dims = {d: coords_table[d] for d in dims.split()}
        ds = apply_dimensions(ds, dims, coords_table)

    if dataset_table:
        ds = update_global_attributes(ds, dataset_table)

    ds = add_version_attribute(ds)

    if mip_table.get("Header"):
        ds = add_header_attributes(ds, mip_table.get("Header"), cv_table)

    if cv_table:
        ds = add_derived_attributes(ds, cv_table)
        check_cv(ds, cv_table)

    # sort attributes
    ds.attrs = collections.OrderedDict(sorted(ds.attrs.items()))

    return ds


def add_variable_attrs(ds, mip_table):
    """add variable attributes"""

    for v in ds.data_vars:
        ds[v].attrs = mip_table[v]
        ds.attrs["variable_id"] = v

    return ds


def apply_dimensions(da, dims, coords_table):
    """apply dimensions from coordinates table"""

    for d, v in dims.items():
        if d in da.coords:
            da = add_coordinate(da, d, v)
            continue
        keys = ["out_name", "standard_name"]
        for k in keys:
            if v[k] in da.coords:
                da = add_coordinate(da, v[k], v)
                break

    return da


def add_coordinate(da, d, axis_entry):
    out_name = axis_entry["out_name"]
    da = da.rename({d: out_name})
    da.coords[out_name].attrs = axis_entry
    dims = da.coords[out_name].dims
    if len(dims) == 1:
        da = da.swap_dims({dims[0]: out_name})

    return da


def add_version_attribute(ds):
    """add version attribute"""
    now = date.today().strftime("%Y%m%d")
    ds.attrs["version"] = now

    return ds


def update_global_attributes(ds, dataset_table):
    ds.attrs.update({k: v for k, v in dataset_table.items() if not k.startswith("#")})

    return ds


def check_cv(ds, cv_table):
    cv = cv_table.get("CV") or cv_table

    req_attrs = cv["required_global_attributes"]

    for attr in req_attrs:
        cv_values = cv.get(attr)
        v = ds.attrs.get(attr)
        if not v:
            warn(f"{attr} not found")
        elif cv_values and v not in list(cv_values):
            warn(f"value '{v}' for '{attr}' not in {list(cv_values)}")


def add_derived_attributes(ds, cv_table):
    cv = cv_table.get("CV") or cv_table

    req_attrs = cv["required_global_attributes"]

    for attr in req_attrs:
        actual_value = ds.attrs.get(attr)
        cv_values = cv.get(attr)
        if isinstance(cv_values, dict) and actual_value in cv_values.keys():
            ds = add_derived_attribute(ds, attr, cv_values.get(actual_value))

    return ds


def add_derived_attribute(ds, attr, cv_values):
    if isinstance(cv_values, str) and attr.endswith("_id"):
        v = cv_values
        k = attr.replace("_id", "")
        warn(f"for attribute '{k}' --> add value '{v}'")
        ds.attrs[k] = v
        return ds

    if isinstance(cv_values, str):
        v = cv_values
        k = attr + "_description"
        warn(f"for attribute '{k}' --> add value '{v}'")
        ds.attrs[k] = v
        return ds

    for k, v in cv_values.items():
        actual_value = ds.attrs.get(k)
        if isinstance(v, list) and actual_value:
            if actual_value not in v:
                warn(f"actual_value '{actual_value}' for {k} not in {v}")
        elif isinstance(v, str) and actual_value is None:
            warn(f"for attribute '{k}' --> add value '{v}'")
            ds.attrs[k] = v
        elif isinstance(v, str) and actual_value:
            if actual_value != v:
                warn(
                    f"attribute '{k}' is set to '{actual_value}' but CV has '{v}' derived from '{attr}'!"
                )
            ds.attrs[k] = v

    return ds


def add_header_attributes(ds, header, cv_table=None):
    default_header_attrs = ["table_id", "realm", "product", "mip_era", "Conventions"]

    if cv_table:
        cv = cv_table.get("CV") or cv_table
        header_attrs = [
            a for a in header.keys() if a in cv["required_global_attributes"]
        ]
    else:
        header_attrs = default_header_attrs

    if "table_id" in header_attrs:
        header["table_id"] = header["table_id"].split()[-1]

    ds.attrs.update({k: header[k] for k in header_attrs})

    return ds
