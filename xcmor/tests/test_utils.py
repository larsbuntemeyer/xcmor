import pytest

from ..utils import parse_cell_methods


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
