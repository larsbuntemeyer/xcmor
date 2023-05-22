import json
import re
import tempfile
from warnings import warn

import xarray as xr
import yaml


def filter_table_by_value(table, key, value):
    return {k: v for k, v in table.items() if v.get(key) == value}


def parse_cell_methods(cm_string):
    # https://stackoverflow.com/questions/52340963/how-to-insert-a-newline-character-before-a-words-that-contains-a-colon
    ys = re.sub(r"(\w+):", r"\n\1:", cm_string).strip()
    d = yaml.safe_load(ys)

    if "area" in d and d.get("area") is None:
        d["area"] = d["time"]

    return d


def _encode_time(time):
    """encode xarray time axis into cf values

    see https://github.com/pydata/xarray/issues/4412

    """
    return xr.conventions.encode_cf_variable(time.variable)


def _read_table(table):
    return read_json(table)


def read_json(filename):
    with open(filename) as f:
        data = json.load(f)
    return data


def write_json(filename, data):
    with open(filename, "w") as fp:
        json.dump(data, fp, indent=4)
    return filename


def _get_cfvarinfo(out_name, table):
    """Returns variable entry from cmor table"""
    if isinstance(table, str):
        table = _read_table(table)
    info = table["variable_entry"].get(out_name, None)
    if info is None:
        raise Exception(
            "{} not found in table {}".format(out_name, get_table_id(table))
        )
    return info


def get_table_id(table):
    """parse the table_id from a cmor table header"""
    separator = " "
    table_id = table["Header"].get("table_id", None)
    if table_id is None:
        raise Exception("no table_id in Header")
    if separator in table_id:
        return table_id.split(separator)[1]
    return table_id


def _tmp_table(table, format="json"):
    """creates a temporay table json file"""
    _, filename = tempfile.mkstemp()
    warn(f"writing temporary table to {filename}")
    if format == "json":
        return write_json(filename, table)


def _get_time_cell_method(cf_varname, table):
    return _strip_time_cell_method(_get_cfvarinfo(cf_varname, table))


def _strip_time_cell_method(cfvarinfo):
    try:
        return cfvarinfo["cell_methods"].split("time:")[1].strip()
    except Exception:
        return None
