import importlib

import pytest


def _importorskip(modname):
    try:
        importlib.import_module(modname)
        has = True
    except ImportError:
        has = False
    func = pytest.mark.skipif(not has, reason=f"requires {modname}")
    return has, func


has_pooch, requires_pooch = _importorskip("pooch")
has_pint_xarray, requires_pint_xarray = _importorskip("pint_xarray")
