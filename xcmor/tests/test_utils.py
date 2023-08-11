import pytest

from ..utils import filter_table_by_value, parse_cell_methods, table_to_dataframe
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
