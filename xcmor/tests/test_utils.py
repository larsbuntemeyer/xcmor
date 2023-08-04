import pytest

from ..utils import filter_table_by_value, parse_cell_methods
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
