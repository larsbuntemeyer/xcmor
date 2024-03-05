import collections
from datetime import date

# from warnings import warn
import cf_xarray as cfxr  # noqa
import xarray as xr
from xarray import DataArray

from .log import get_logger
from .resources import get_project_tables
from .rules import rules
from .tests.tables import coords as coords_default
from .utils import cf_table

logger = get_logger(__name__)


def _encode_time(ds, cf_units=None):
    """Encode time units and calendar"""
    time = ds.cf["time"]
    cf_units = cf_units or time.attrs.get("units") or time.encoding.get("units")

    if cf_units is None:
        cf_units = "days since ?"
    else:
        del ds.time.attrs["units"]

    start_format = "%Y-%m-%dT%H:%M:%S"

    # check if time is datetime-like, maybe there is a better way?
    # decode times if not datetime-like
    try:
        time.dt
    except Exception:
        cf_units = cf_units.replace("?", "1950")
        logger.warning(
            f"time axis does not seem to be datetime-like, encoding with units '{cf_units}'"
        )

        ds.time.attrs["units"] = cf_units  # .replace("?", "1950")
        ds = xr.decode_cf(ds, decode_times=True, decode_coords=False)
        return ds.time

    start_str = f"{time[0].dt.strftime(start_format).item()}"
    units = cf_units.replace("?", start_str)

    logger.debug(f"setting time units: {units}")
    time.encoding["units"] = units
    return ds.cf["time"]


def _units_convert(da, cf_units, format=None):
    """Use pint_xarray to convert units"""
    import pint_xarray  # noqa
    from cf_xarray.units import units  # noqa

    if format is None:
        format = "cf"
    if units.Unit(da.units) != units.Unit(cf_units):
        logger.warn(
            f"converting units {da.units} from input data to CF units {cf_units}"
        )
        da_quant = da.pint.quantify()
        da = da_quant.pint.to(cf_units).pint.dequantify(format=format)
    return da


def _get_x_y_coords(obj):
    """Guess linear X and Y coordinates"""
    obj = obj.cf.guess_coord_axis()
    X = None
    Y = None
    # cfxr finds the X and Y coordinates right away
    if "X" in obj.cf.coords and "Y" in obj.cf.coords:
        X = obj.cf["X"]
        Y = obj.cf["Y"]
    # cfxr finds longitude and latitude, let's check if they are 1D
    elif "longitude" in obj.cf.coords and "latitude" in obj.cf.coords:
        lon = obj.cf["longitude"]
        lat = obj.cf["latitude"]
        if lon.ndim == 1 and lat.ndim == 1:
            X = lon
            Y = lat
    # ensure the attributes to make CF conform
    if X is not None and Y is not None:
        X.attrs["axis"] = "X"
        Y.attrs["axis"] = "Y"
    return X, Y


def _get_lon_lat(obj):
    """Return lon and lat extracted from ds

    Use cf_xarray to identify longitude and latitude coordinates.
    Might be 1D or 2D coordinates.

    """
    obj = obj.copy().cf.guess_coord_axis()
    try:
        lon = obj.cf["longitude"]
        lat = obj.cf["latitude"]
    except KeyError:
        raise KeyError("could not identify longitude/latitude")

    return lon, lat


def _is_curvilinear(obj):
    """Check for curvilinear

    Pretty naive definition here, curvilinear for us here simply
    means if longitude and latitude are not 1D coordinates.

    """
    lon, lat = _get_lon_lat(obj)
    return lon.ndim > 1 and lat.ndim > 1


def _guess_dims_attr(obj):
    """Try to guess dimensions attribute"""
    obj = obj.copy().cf.guess_coord_axis()
    dimensions = []
    try:
        lon, lat = _get_lon_lat(obj)
        dimensions.extend(["longitude", "latitude"])
    except KeyError:
        logger.warning(
            f"Could not guess longitude and latitude coordinates from {list(obj.coords)}"
        )
    if "Z" in obj.cf.coords:
        dimensions.append(obj.cf.coords["Z"].name)
    if "time" in obj.cf.coords:
        dimensions.append("time")
    return dimensions


def _add_var_attrs(ds, mip_table):
    """add variable attributes"""

    for var in ds.data_vars:
        ds[var].attrs = {k: v for k, v in mip_table[var].items() if v}
        ds.attrs["variable_id"] = var
        if mip_table[var].get("frequency"):
            ds.attrs["frequency"] = mip_table[var].get("frequency")
            del ds[var].attrs["frequency"]
        if mip_table[var].get("modeling_realm"):
            ds.attrs["realm"] = mip_table[var].get("modeling_realm")
            del ds[var].attrs["modeling_realm"]

    return ds


def _interpret_var_attrs(ds, mip_table):
    """Apply variable attributes found in the mip table.

    This will interpret attributes found in the mip table, e.g.,
    valid_min, valid_max, convert dtypes, etc...
    Once attributes were interpreted they are removed from the
    variables attributes dictionary.

    """

    for v in ds.data_vars:
        da = ds[v]
        for attr in da.attrs.copy():
            if hasattr(rules, attr):
                da = getattr(rules, attr)(da)
        ds = ds.assign({da.name: da})  # .drop_vars(v)

    return ds


def _interpret_coord_attrs(ds, time_units=None):
    """Apply coordinates attributes.

    This will interpret attributes found in the mip table, e.g.,
    valid_min, valid_max, convert dtypes, etc...
    Once attributes were interpreted they are removed from the
    variables attributes dictionary.

    """

    for v in ds.coords:
        da = ds.coords[v]
        # print(v, da)
        # logger.info(v)
        for attr in da.attrs.copy():
            # logger.info(attr)
            if hasattr(rules, attr):
                da = getattr(rules, attr)(da)
        ds = ds.assign_coords({da.name: da})  # .drop_vars(v)

    if "time" in ds.cf.coords:
        ds = ds.cf.assign_coords(time=_encode_time(ds, time_units))

    return ds


def _find_coord_key(da, axis_entry):
    """find datarray coordinate by cf attributes from coordinates table"""

    keys = ["out_name", "standard_name", "axis"]
    for k in keys:
        if axis_entry[k] in da.cf.coords or axis_entry[k] in da.coords:
            # print(f"found {v[k]} by {k}")
            return axis_entry[k]
    return None


def _add_coord_attrs(da, axis_entry):
    """Add coordinate attributes from coordinates table"""

    out_name = axis_entry["out_name"]
    coord_key = out_name

    if coord_key not in da.coords:
        coord_key = _find_coord_key(da, axis_entry)

    if coord_key is None:
        # we could not find the coordinate in the dataset
        logger.info(f"adding coordinate: {out_name}")
        value = float(axis_entry["value"])
        da = da.assign_coords({out_name: DataArray(value)})
    else:
        # rename coord key to actual coordinate out_name
        da = da.cf.rename({coord_key: out_name})

    # add required attributes
    da.coords[out_name].attrs = {k: v for k, v in axis_entry.items() if v}

    dims = da.coords[out_name].dims

    # this is a coordinate variable (not auxilliary), swap dims
    if len(dims) == 1:
        da = da.swap_dims({dims[0]: out_name})

    return da


def _apply_dims(da, dims):
    """Apply dimensions from coordinates table

    Parameters
    ----------
    da : DataArray, Dataset
        DataArray of which coordinates should be cmorized.
    dims : dict
        Dictionary with dimension names a keys and cmor coordinate
        table entries as values.

    """

    # d is a cmor coordinate table key, v is the coordinates table entry
    for d, v in dims.items():
        if v:
            da = _add_coord_attrs(da, v)
        else:
            logger.warning(f"found no coordinate attributes for coordinate '{d}'")
        # we find the coordinate already by its correct cf out_name
        # if v["out_name"] in da.coords:
        #     da = _add_coord_attrs(da, d, v)
        #     continue

        # # search for a coordinate by attributes (using cf_xarray)
        # keys = ["out_name", "standard_name", "axis"]
        # for k in keys:
        #     if v[k] in da.cf.coords or v[k] in da.coords:
        #         # print(f"found {v[k]} by {k}")
        #         da = _add_coord_attrs(da, v[k], v)
        #         break

        # # seems to be a scalar coordinate that we need to create
        # if v["out_name"] not in da.coords:
        #     logger.info(f"adding coordinate: {d}")
        #     value = float(v["value"])
        #     coord = DataArray(value)
        #     dtype = v["type"]
        #     coord = DataArray(value).astype(dtype)
        #     da = da.assign_coords({v["out_name"]: coord})
        #     da.coords[v["out_name"]].attrs = v

    return da


def _interpret_var_dims(ds, coords_table, drop=True):
    """Interpret variable dimensions attribute.

    This will look up the dimensions defined for variables
    in the mip table and update coordinates acoording to
    meta data in the coordinates table.

    """
    all_dims = []

    for var in ds.data_vars:
        dims = ds[var].attrs.get("dimensions")
        if not dims:
            dims = _guess_dims_attr(ds[var])
            logger.debug(f"dims of {var}: {dims}")
        else:
            del ds[var].attrs["dimensions"]
            dims = dims.split()

        dims = {d: coords_table.get(d) or {"out_name": d} for d in dims}

        ds = _apply_dims(ds, dims)
        # add coordinates attribute, e.g., for 0D coordinate variables
        # e.g., height2m, if not a dataarray index.
        coordinates = " ".join(
            [
                d["out_name"]
                for d in dims.values()
                if d["out_name"] not in ds[var].indexes and d["out_name"] in ds.coords
            ]
        )

        if coordinates:
            ds[var].attrs["coordinates"] = coordinates

        all_dims.extend([v["out_name"] for v in dims.values()])

    logger.debug(f"added coordinates: {list(all_dims)}")
    # drop unneccessary coordinates
    if drop is True:
        drops = [c for c in ds.coords if c not in list(all_dims)]
        logger.debug(f"dropping coordinates: {drops}")
        ds = ds.drop(drops)

    return ds


def _add_version_attr(ds):
    """add version attribute"""
    now = date.today().strftime("%Y%m%d")
    ds.attrs["version"] = now

    return ds


def _update_global_attrs(ds, dataset_table):
    ds.attrs.update(
        {
            k: v
            for k, v in dataset_table.items()
            if (not k.startswith("#") and not k.startswith("_"))
        }
    )

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
        logger.debug(f"for attribute '{k}' --> add value '{v}'")
        ds.attrs[k] = v
        return ds

    if isinstance(cv_values, str):
        # for all attributes that have a description in the CV
        # we add another attribute ending on "*_info" that
        # adds the description, .e.g, frequency
        v = cv_values
        k = attr + "_info"
        logger.debug(f"for attribute '{k}' --> add value '{v}'")
        ds.attrs[k] = v
        return ds

    for k, v in cv_values.items():
        # if cv_values is a dict,
        actual_value = ds.attrs.get(k)
        if isinstance(v, list) and actual_value:
            if actual_value not in v:
                logger.warn(
                    f"attribute '{attr}' has value '{ds.attrs.get(attr)}' but attribute '{k}' has value '{actual_value}' which is not in the list of expected values: {v}"
                )
        elif isinstance(v, str) and actual_value is None:
            message = f"attribute '{attr}' has value '{ds.attrs.get(attr)}' and requires attribute '{k}' to be set to '{v}'"
            logger.info(message)
            ds.attrs[k] = v
        elif isinstance(v, str) and actual_value:
            if actual_value != v:
                logger.warn(
                    f"attribute '{attr}' has value '{ds.attrs.get(attr)}' but attribute '{k}' is set to '{actual_value}' but CV requires '{v}'!"
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
    grids_table=None,
    mapping_table=None,
    guess=True,
    time_units=None,
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
    grids_table: dict, str
        The grids table, can either be a dictionary or a path to a cmor table
        in json format.
    mapping_table: dict
        The mapping table maps input variable names to mip table variable keys.
    time_units: str
        Time units for NetCDF encoding. Default is ``days since`` the beginning of the
        time interval.

    Returns
    -------
    Cmorized Dataset.

    """
    guess = True

    ds = ds.copy()

    # ensure dataset
    if isinstance(ds, DataArray):
        ds = ds.to_dataset()

    if mip_table is None:
        logger.debug("using default cf variable table")
        mip_table = cf_table().to_dict(orient="index")

    if coords_table is None:
        logger.debug("using default coords table")
        coords_table = coords_default

    if guess is True:
        ds = ds.cf.guess_coord_axis(verbose=False)

    # ensure grid mappings and bounds in coords, not in data_vars
    ds = xr.decode_cf(ds, decode_coords="all")

    if mapping_table is not None:
        ds = ds.rename_vars({v: (mapping_table.get(v) or v) for v in ds})

    # add variable attributes from mip table entries
    ds = _add_var_attrs(ds, mip_table.get("variable_entry") or mip_table)

    # interprets variable attributes
    ds = _interpret_var_attrs(ds, mip_table.get("variable_entry") or mip_table)

    if coords_table:
        ds = _interpret_var_dims(ds, coords_table.get("axis_entry") or coords_table)
        ds = _interpret_coord_attrs(ds, time_units)

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

    def cmorize(
        self, ds, mip_table, dataset_table, mapping_table=None, time_units=None
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
            in json format or a table_id from the MIP.
        dataset_table : dict, str
            The input dataset cmor table, can either be a dictionary or a path to a cmor table
            in json format.
        mapping_table: dict
            The mapping table maps input variable names to mip table variable keys.
        time_units: str
            Time units for NetCDF encoding. Default is ``days since`` the beginning of the
            time interval.

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
            mapping_table=mapping_table,
            time_units=time_units,
        )
