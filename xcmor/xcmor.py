import collections
from datetime import date

# from warnings import warn
import cf_xarray as cfxr  # noqa
import numpy as np
from xarray import DataArray

from .log import get_logger
from .mapping import dtype_map
from .resources import get_project_tables

logger = get_logger(__name__)


def _guess_X_Y_coords(ds):
    """Guess linear X and Y coordinates"""
    ds = ds.cf.guess_coord_axis()
    X = None
    Y = None
    if all([c in ds.cf.coords for c in ["X", "Y"]]):
        X = ds.cf["X"].name
        Y = ds.cf["Y"].name
    elif all([c in ds.cf.coords for c in ["longitude", "latitude"]]):
        lon = ds.cf["longitude"]
        lat = ds.cf["latitude"]
        if lon.ndim == 1 and lat.ndim == 1:
            X = lon.name
            Y = lat.name
    return X, Y


def _add_var_attrs(ds, mip_table):
    """add variable attributes"""

    for v in ds.data_vars:
        ds[v].attrs = mip_table[v]
        ds.attrs["variable_id"] = v
        if mip_table[v].get("frequency"):
            ds.attrs["frequency"] = mip_table[v].get("frequency")
            del ds[v].attrs["frequency"]
        if mip_table[v].get("modeling_realm"):
            ds.attrs["realm"] = mip_table[v].get("modeling_realm")
            del ds[v].attrs["modeling_realm"]

    return ds


def _interpret_var_attrs(ds, mip_table):
    """Apply variable attributes found in the mip table.

    This will interpret attributes found in the mip table, e.g.,
    valid_min, valid_max, convert dtypes, etc...
    Once attributes were interpreted they are removed from the
    variables attributes dictionary.

    """

    for v in ds.data_vars:
        attrs = ds[v].attrs
        if ds[v].dtype != dtype_map[attrs["type"]]:
            logger.warning(
                f"converting {v} from {ds[v].dtype} to {dtype_map[attrs['type']]}"
            )
            ds[v] = ds[v].astype(dtype_map[attrs["type"]])
        del ds[v].attrs["type"]

        ds.rename({v: attrs.get("out_name") or v})
        del ds[v].attrs["out_name"]

        valid_min, valid_max = attrs.get("valid_min"), attrs.get("valid_max")
        del ds[v].attrs["valid_min"]
        del ds[v].attrs["valid_max"]
        if valid_min:
            assert ds[v].min() >= valid_min

        if valid_max:
            assert ds[v].max() <= valid_max
            del ds[v].attrs["valid_max"]

    return ds


def _interpret_var_dims(ds, coords_table):
    """Interpret variable dimensions attribute.

    This will look up the dimensions defined for variables
    in the mip table and update coordinates acoording to
    meta data in the coordinates table.

    """
    for var in ds.data_vars:
        dims = ds[var].attrs.get("dimensions")
        dims = {d: coords_table[d] for d in dims.split()}
        ds = _apply_dims(ds, dims, coords_table)
        # add coordinates attribute, e.g., for 0D coordinate variables
        # e.g., height2m, etc...
        coordinates = " ".join(
            [d["out_name"] for d in dims.values() if d["out_name"] not in ds.indexes]
        )

        if coordinates:
            print(f"coordinates: {coordinates}")
            ds[var].attrs["coordinates"] = coordinates

        # del ds[var].attrs["dimensions"]
    return ds


def _apply_dims(da, dims, coords_table):
    """Apply dimensions from coordinates table"""

    for d, v in dims.items():
        # we find the coordinate already by its correct cf name
        if d in da.coords:
            da = _add_coord(da, d, v)
            continue

        # search for a coordinate by attributes (using cf_xarray)
        keys = ["out_name", "standard_name", "axis"]
        for k in keys:
            if v[k] in da.cf.coords or v[k] in da.coords:
                # print(f"found {v[k]} by {k}")
                da = _add_coord(da, v[k], v)
                break

        # seems to be a scalar coordinate that we need to create
        if v["out_name"] not in da.coords:
            logger.info(f"adding coordinate: {d}")
            value = float(v["value"])
            dtype = v["type"]
            coord = DataArray(value).astype(dtype)
            da = da.assign_coords({v["out_name"]: coord})
            da.coords[v["out_name"]].attrs = v

    return da


def _add_coord(da, d, axis_entry):
    """Add coordinate attributes from coordinates table"""
    out_name = axis_entry["out_name"]
    da = da.cf.rename({d: out_name})
    da.coords[out_name].attrs = axis_entry
    dims = da.coords[out_name].dims
    if len(dims) == 1:
        da = da.swap_dims({dims[0]: out_name})

    dtype = axis_entry.get("type")
    if dtype:
        da[out_name] = da[out_name].astype(dtype)

    requested = axis_entry.get("requested")
    if requested:
        requested = list(map(float, requested))
        # print(f"requested: {requested}")
        # print(f"values: {da[out_name].values}")
        assert np.allclose(da[out_name].values, requested)

    return da


def _add_version_attr(ds):
    """add version attribute"""
    now = date.today().strftime("%Y%m%d")
    ds.attrs["version"] = now

    return ds


def _update_global_attrs(ds, dataset_table):
    ds.attrs.update({k: v for k, v in dataset_table.items() if not k.startswith("#")})

    return ds


def _check_cv(ds, cv_table):
    cv = cv_table.get("CV") or cv_table

    req_attrs = cv["required_global_attributes"]

    for attr in req_attrs:
        cv_values = cv.get(attr)
        v = ds.attrs.get(attr)
        if not v:
            logger.warn(f"{attr} not found")
        elif cv_values and v not in list(cv_values):
            logger.warn(f"value '{v[0:50]}...' for '{attr}' not in {list(cv_values)}")


def _add_derived_attrs(ds, cv_table):
    """Add derived global attributes from CV

    Attributes in the CV table might contain derived
    attributes that we add automatically.


    """
    cv = cv_table.get("CV") or cv_table

    req_attrs = cv["required_global_attributes"]

    for attr in req_attrs:
        actual_value = ds.attrs.get(attr)
        cv_values = cv.get(attr)
        if isinstance(cv_values, dict) and actual_value in cv_values.keys():
            ds = _add_derived_attr(ds, attr, cv_values.get(actual_value))

    return ds


def _add_derived_attr(ds, attr, cv_values):
    if isinstance(cv_values, str) and attr.endswith("_id"):
        # for all attributes that end with "*_id", and that
        # have a description in the CV, we add another attribute
        # containing that description, e.g. institution_id
        v = cv_values
        k = attr.replace("_id", "")
        logger.info(f"for attribute '{k}' --> add value '{v}'")
        ds.attrs[k] = v
        return ds

    if isinstance(cv_values, str):
        # for all attributes that have a description in the CV
        # we add another attribute ending on "*_info" that
        # adds the description, .e.g, frequency
        v = cv_values
        k = attr + "_info"
        logger.info(f"for attribute '{k}' --> add value '{v}'")
        ds.attrs[k] = v
        return ds

    for k, v in cv_values.items():
        # if cv_values is a dict,
        actual_value = ds.attrs.get(k)
        if isinstance(v, list) and actual_value:
            if actual_value not in v:
                logger.warn(
                    f"actual_value '{actual_value}' for {k} not in list of expected values: {v}"
                )
        elif isinstance(v, str) and actual_value is None:
            message = (
                f"for attribute '{k}' --> add value '{v}' because attribute '{attr}' in CV "
                f"asks for attribute '{k}' to be set to '{v}'"
            )
            logger.info(message)
            ds.attrs[k] = v
        elif isinstance(v, str) and actual_value:
            if actual_value != v:
                logger.warn(
                    f"attribute '{k}' is set to '{actual_value}' but CV has '{v}' derived from '{attr}'!"
                )
            ds.attrs[k] = v

    return ds


def _add_header_attrs(ds, header, cv_table=None):
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


def cmorize(
    ds,
    mip_table=None,
    coords_table=None,
    dataset_table=None,
    cv_table=None,
    mapping=None,
):
    """Lazy cmorization.

    Cmorizes an xarray Dataset or DataArray object. The cmorizations tries
    to follow the approach of the original `cmor <https://github.com/PCMDI/cmor>`_
    library in adding, manipulating and interpreting dataseta attributes and
    cmor table vocabulary.

    Parameters
    ----------
    ds : DataArray, Dataset
        Dataset that should be cmorized.
    mip_table : dict, str
        The MIP table, can either be a dictionary or a path to a cmor table
        in json format.
    coords_table : dict, str
        The cmor coordinates table, can either be a dictionary or a path to a cmor table
        in json format.
    dataset_table : dict, str
        The input dataset cmor table, can either be a dictionary or a path to a cmor table
        in json format.
    cv_table: dict, str
        The controlled vocabulary table, can either be a dictionary or a path to a cmor table
        in json format.
    mapping: dict
        The mapping table mapping input variable names to cmor table axis entry keys.

    Returns
    -------
    Cmorized Dataset.

    """

    ds = ds.copy()

    ds = ds.cf.guess_coord_axis(verbose=False)

    if coords_table is None:
        coords_table = {}
    else:
        coords_table = coords_table.get("axis_entry") or coords_table

    if mapping is None:
        mapping = {}

    if isinstance(ds, DataArray):
        ds = ds.to_dataset()

    ds = ds.rename({v: (mapping.get(v) or v) for v in ds})

    ds = _add_var_attrs(ds, mip_table.get("variable_entry") or mip_table)

    ds = _interpret_var_attrs(ds, mip_table.get("variable_entry") or mip_table)

    if coords_table:
        ds = _interpret_var_dims(ds, coords_table)

    if dataset_table:
        ds = _update_global_attrs(ds, dataset_table)

    ds = _add_version_attr(ds)

    if mip_table.get("Header"):
        ds = _add_header_attrs(ds, mip_table.get("Header"), cv_table)

    if cv_table:
        ds = _add_derived_attrs(ds, cv_table)
        _check_cv(ds, cv_table)

    # sort attributes
    ds.attrs = collections.OrderedDict(sorted(ds.attrs.items()))

    return ds


class Cmorizer:
    def __init__(self, project=None, url=None, template=None):
        """
        Cmorizer class supporting preconfigured MIPs.

        Parameters
        ----------
        project : str, optional
            Pre-configures MIP, e.g.,

            - CMIP6
            - CORDEX

        url : str, optional
            Base URL or directory of cmor tables.

        template : str, optional
            CMOR talbe naming template, e.g.::

                CMIP6_{table_id}.json

            e.g. CMIP6_Amon.json

        Returns
        -------
        cmorizer : Cmorizer object.

        """
        self._init_tables(project, url, template)

    def _init_tables(self, project, url, template):
        if project is None and url is None:
            self.project = "CMIP6"
        else:
            self.project = project
        self.tables = get_project_tables(url, self.project, template)

    @property
    def required(self):
        """List required global attributes."""
        return self.tables.cv["CV"].get("required_global_attributes")

    def cmorize(self, ds, mip_table, dataset_table, mapping=None):
        """Lazy cmorization.

        Cmorizes an xarray Dataset or DataArray object. The cmorizations tries
        to follow the approach of the original `cmor <https://github.com/PCMDI/cmor>`_
        library in adding, manipulating and interpreting dataseta attributes and
        cmor table vocabulary.

        Parameters
        ----------
        ds : DataArray, Dataset
            Dataset that should be cmorized.
        mip_table : dict, str
            The MIP table, can either be a dictionary or a path to a cmor table
            in json format or a table_id from the MIP.
        dataset_table : dict, str
            The input dataset cmor table, can either be a dictionary or a path to a cmor table
            in json format.
        mapping : dict
            The mapping table mapping input variable names to cmor table axis entry keys.

        Returns
        -------
        Cmorized Dataset.

        Examples
        --------

        >>> from xcmor.datasets import reg_ds
        >>> from xcmor import Cmorizer
        >>>
        >>> cmor = Cmorizer()
        >>> ds_out = cmor.cmorize(
        ...     reg_ds.rename(temperature="tas").tas,
        ...     "Amon",
        ...     cmor.tables["input_example"],
        ... )

        """

        if not isinstance(mip_table, dict):
            mip_table = self.tables[mip_table]

        return cmorize(
            ds,
            mip_table=mip_table,
            dataset_table=dataset_table,
            coords_table=self.tables.coords,
            cv_table=self.tables.cv,
            mapping=mapping,
        )
