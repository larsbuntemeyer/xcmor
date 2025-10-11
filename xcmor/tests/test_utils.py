import pytest

from ..resources import ProjectTables, cmip6, cordex, retrieve_cmor_table
from ..utils import (
    _tmp_table,
    filter_table_by_value,
    parse_cell_methods,
    posix_to_python_regex,
    read_tables,
    table_to_dataframe,
)
from . import requires_pooch
from .tables import mip_amon


@pytest.mark.parametrize(
    ["cell_methods", "expected"],
    (
        ("time: mean", {"time": "mean"}),
        ("area: mean time: mean", {"area": "mean", "time": "mean"}),
        ("area: time: mean", {"area": "mean", "time": "mean"}),
    ),
)
def test_parse_cell_methods(cell_methods, expected):
    parsed = parse_cell_methods(cell_methods)
    assert parsed == expected


def test_filter_by_value():
    filtered = filter_table_by_value(mip_amon, "standard_name", "air_temperature")
    expected = ["ta", "tas"]
    assert list(filtered.keys()) == expected
    for k, v in filtered.items():
        assert v == mip_amon[k]

    filtered = filter_table_by_value(mip_amon, "standard_name", "precipitation_flux")
    expected = ["pr"]
    assert list(filtered.keys()) == expected
    for k, v in filtered.items():
        assert v == mip_amon[k]


def test_table_to_dataframe():
    df = table_to_dataframe(mip_amon)
    assert df.index.to_list() == list(mip_amon.keys())
    df = table_to_dataframe(mip_amon, index_name="variable_key")
    assert df.variable_key.to_list() == list(mip_amon.keys())
    assert df.standard_name.to_list() == [v["standard_name"] for v in mip_amon.values()]
    assert df.units.to_list() == [v["units"] for v in mip_amon.values()]


def test_read_tables_decorator():
    table = {"x": 1, "y": 2}
    table_fname = _tmp_table(table)

    @read_tables(tables=["mip_table", "coords_table"])
    def cmorize(mip_table=None, coords_table=None):
        """Cmorizer function"""
        return mip_table, coords_table

    assert cmorize(mip_table=table, coords_table=table_fname) == (
        {"x": 1, "y": 2},
        {"x": 1, "y": 2},
    )


@requires_pooch
@pytest.mark.parametrize(
    ["table_id", "project"],
    (
        ("Amon", "CMIP6"),
        ("day", "CMIP6"),
    ),
)
def test_table_retrieve(table_id, project):
    table = retrieve_cmor_table(table_id, project)
    assert table["Header"]["table_id"] == f"Table {table_id}"


def test_url_tables():
    tables = ProjectTables("table/dir", template="MIP_{table_id}.json")
    assert tables.get_url("coordinate") == "table/dir/MIP_coordinate.json"
    tables = ProjectTables("table/dir", template="MIP_{table_id}.yaml")
    assert tables.get_url("coordinate") == "table/dir/MIP_coordinate.yaml"


@requires_pooch
def test_project_tables():
    assert isinstance(cmip6["coordinate"], dict)
    assert cmip6.grids["Header"]["table_id"] == "Table grids"
    assert isinstance(cmip6.coords, dict)
    assert isinstance(cmip6.terms, dict)
    assert isinstance(cmip6["Amon"], dict)
    assert isinstance(cordex["mon"], dict)


@pytest.mark.parametrize(
    "posix_pattern,expected_py,valid,invalid",
    [
        (
            r"r[[:digit:]]\{1,\}i[[:digit:]]\{1,\}p[[:digit:]]\{1,\}f[[:digit:]]\{1,\}$",
            r"r[0-9]{1,}i[0-9]{1,}p[0-9]{1,}f[0-9]{1,}$",
            ["r1i1p1f1", "r12i3p45f6"],
            ["xr1i1p1f1", "r1i1p1f", "r1i1p1f1x"],
        ),
        (
            r"[[:alpha:]]\{3\}[[:digit:]]\{2,4\}",
            r"[A-Za-z]{3}[0-9]{2,4}",
            ["Abc12", "abc1234"],
            ["ab12", "AB12", "abcd123"],
        ),
        (
            r"a[[:digit:]]\{1,\}b",
            r"a[0-9]{1,}b",
            ["a1b", "a123b"],
            ["ab", "a_b"],
        ),
        (
            r"foo\([[:digit:]]\{2\}\)",
            r"foo([0-9]{2})",
            ["foo12", "foo99"],
            ["foo1", "foo123", "foo2"],
        ),
        (
            # pattern using an unsupported collating symbol should pass through unchanged
            r"[.ch.]",
            r"[.ch.]",
            [".", "c", "h"],
            ["ch"],
        ),
    ],
)
def test_posix_to_python_regex(posix_pattern, expected_py, valid, invalid):
    """Test conversion of POSIX BRE subset to Python regex and matching semantics."""
    converted = posix_to_python_regex(posix_pattern)
    assert converted == expected_py
    # compile; if pattern unchanged (e.g., unsupported features) it may not be valid Python
    # In that case we skip match assertions.
    import re as _re

    try:
        regex_obj = _re.compile(converted)
    except Exception:
        # If compilation fails we only assert conversion outcome
        return

    for s in valid:
        assert regex_obj.fullmatch(s), f"Should match: {s} with {converted}"
    for s in invalid:
        assert not regex_obj.fullmatch(s), f"Should NOT match: {s} with {converted}"


def test_posix_to_python_regex_idempotent_simple():
    # Patterns already Python style should stay the same
    pat = r"abc[0-9]{2,}"  # no POSIX classes
    assert posix_to_python_regex(pat) == pat
