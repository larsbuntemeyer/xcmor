import pytest

from ..resources import ProjectTables, cmip6, cordex, retrieve_cmor_table
from ..utils import filter_table_by_value, parse_cell_methods, table_to_dataframe
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
    assert isinstance(cmip6["Amon"], dict)
    assert isinstance(cordex["mon"], dict)
