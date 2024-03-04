import numpy as np
import pandas as pd
import xarray as xr


def create_regular_ds():
    np.random.seed(0)
    temperature = 15 + 8 * np.random.randn(2, 2, 3)
    precipitation = 10 * np.random.rand(2, 2, 3)
    lon = [-99.83 + 180.0, -99.32 + 180.0]
    lat = [42.25, 42.21]
    time = pd.date_range("2014-09-06", periods=3)
    reference_time = pd.Timestamp("2014-09-05")

    ds = xr.Dataset(
        data_vars=dict(
            temperature=(["x", "y", "time"], temperature),
            precipitation=(["x", "y", "time"], precipitation),
        ),
        coords=dict(
            lon=(["x"], lon),
            lat=(["y"], lat),
            time=time,
            reference_time=reference_time,
        ),
        attrs=dict(description="Weather related data."),
    )

    return ds


def create_curvilinear_ds():
    np.random.seed(0)
    temperature = 15 + 8 * np.random.randn(2, 2, 3)
    precipitation = 10 * np.random.rand(2, 2, 3)
    lon = [[-99.83, -99.32], [-99.79, -99.23]]
    lat = [[42.25, 42.21], [42.63, 42.59]]
    time = pd.date_range("2014-09-06", periods=3)
    reference_time = pd.Timestamp("2014-09-05")

    ds = xr.Dataset(
        data_vars=dict(
            temperature=(["x", "y", "time"], temperature),
            precipitation=(["x", "y", "time"], precipitation),
        ),
        coords=dict(
            lon=(["x", "y"], lon),
            lat=(["x", "y"], lat),
            time=time,
            reference_time=reference_time,
        ),
        attrs=dict(description="Weather related data."),
    )

    return ds


def create_plev_ds():
    ta = 10.0 * np.random.random_sample((2, 19, 3, 4)) + 250.0
    lat = np.array([10, 20, 30])

    lon = np.array([0, 90, 180, 270])

    time = np.array([15.5, 45])

    lev = np.array(
        [
            100000,
            92500,
            85000,
            70000,
            60000,
            50000,
            40000,
            30000,
            25000,
            20000,
            15000,
            10000,
            7000,
            5000,
            3000,
            2000,
            1000,
            500,
            100,
        ]
    )

    ds = xr.Dataset(
        data_vars=dict(
            temperature=(["time", "z", "y", "x"], ta),
        ),
        coords=dict(
            lon=(["x"], lon),
            lat=(["y"], lat),
            lev=(["z"], lev),
            time=time,
        ),
        attrs=dict(description="Weather related data."),
    )

    return ds


temp_ds = create_curvilinear_ds()
reg_ds = create_regular_ds()
plev_ds = create_plev_ds()
