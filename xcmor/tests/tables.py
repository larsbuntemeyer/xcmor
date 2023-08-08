mip_amon = {
    "ta": {
        "frequency": "mon",
        "modeling_realm": "atmos",
        "standard_name": "air_temperature",
        "units": "K",
        "cell_methods": "time: mean",
        "cell_measures": "area: areacella",
        "long_name": "Air Temperature",
        "comment": "Air Temperature",
        "dimensions": "longitude latitude plev19 time",
        "out_name": "ta",
        "type": "real",
        "positive": "",
        "valid_min": "",
        "valid_max": "",
        "ok_min_mean_abs": "",
        "ok_max_mean_abs": "",
    },
    "tas": {
        "frequency": "mon",
        "modeling_realm": "atmos",
        "standard_name": "air_temperature",
        "units": "K",
        "cell_methods": "area: time: mean",
        "cell_measures": "area: areacella",
        "long_name": "Near-Surface Air Temperature",
        "comment": "near-surface (usually, 2 meter) air temperature",
        "dimensions": "longitude latitude time height2m",
        "out_name": "tas",
        "type": "real",
        "positive": "",
        "valid_min": "",
        "valid_max": "",
        "ok_min_mean_abs": "",
        "ok_max_mean_abs": "",
    },
    "pr": {
        "frequency": "mon",
        "modeling_realm": "atmos",
        "standard_name": "precipitation_flux",
        "units": "kg m-2 s-1",
        "cell_methods": "area: time: mean",
        "cell_measures": "area: areacella",
        "long_name": "Precipitation",
        "comment": "includes both liquid and solid phases",
        "dimensions": "longitude latitude time",
        "out_name": "pr",
        "type": "real",
        "positive": "",
        "valid_min": "",
        "valid_max": "",
        "ok_min_mean_abs": "",
        "ok_max_mean_abs": "",
    },
}


coords = {
    "axis_entry": {
        "gridlatitude": {
            "standard_name": "grid_latitude",
            "units": "degrees",
            "axis": "Y",
            "long_name": "Grid Latitude",
            "climatology": "",
            "formula": "",
            "must_have_bounds": "yes",
            "out_name": "rlat",
            "positive": "",
            "requested": "",
            "requested_bounds": "",
            "stored_direction": "increasing",
            "tolerance": "",
            "type": "double",
            "valid_max": "90.0",
            "valid_min": "-90.0",
            "value": "",
            "z_bounds_factors": "",
            "z_factors": "",
            "bounds_values": "",
            "generic_level_name": "",
        },
        "height100m": {
            "standard_name": "height",
            "units": "m",
            "axis": "Z",
            "long_name": "height",
            "climatology": "",
            "formula": "",
            "must_have_bounds": "no",
            "out_name": "height",
            "positive": "up",
            "requested": "",
            "requested_bounds": "",
            "stored_direction": "increasing",
            "tolerance": "",
            "type": "double",
            "valid_max": "120.0",
            "valid_min": "80.0",
            "value": "100.",
            "z_bounds_factors": "",
            "z_factors": "",
            "bounds_values": "",
            "generic_level_name": "",
        },
        "height10m": {
            "standard_name": "height",
            "units": "m",
            "axis": "Z",
            "long_name": "height",
            "climatology": "",
            "formula": "",
            "must_have_bounds": "no",
            "out_name": "height",
            "positive": "up",
            "requested": "",
            "requested_bounds": "",
            "stored_direction": "increasing",
            "tolerance": "",
            "type": "double",
            "valid_max": "30.0",
            "valid_min": "1.0",
            "value": "10.",
            "z_bounds_factors": "",
            "z_factors": "",
            "bounds_values": "",
            "generic_level_name": "",
        },
        "height2m": {
            "standard_name": "height",
            "units": "m",
            "axis": "Z",
            "long_name": "height",
            "climatology": "",
            "formula": "",
            "must_have_bounds": "no",
            "out_name": "height",
            "positive": "up",
            "requested": "",
            "requested_bounds": "",
            "stored_direction": "increasing",
            "tolerance": "",
            "type": "double",
            "valid_max": "10.0",
            "valid_min": "1.0",
            "value": "2.",
            "z_bounds_factors": "",
            "z_factors": "",
            "bounds_values": "",
            "generic_level_name": "",
        },
        "latitude": {
            "standard_name": "latitude",
            "units": "degrees_north",
            "axis": "Y",
            "long_name": "Latitude",
            "climatology": "",
            "formula": "",
            "must_have_bounds": "yes",
            "out_name": "lat",
            "positive": "",
            "requested": "",
            "requested_bounds": "",
            "stored_direction": "increasing",
            "tolerance": "",
            "type": "double",
            "valid_max": "90.0",
            "valid_min": "-90.0",
            "value": "",
            "z_bounds_factors": "",
            "z_factors": "",
            "bounds_values": "",
            "generic_level_name": "",
        },
        "longitude": {
            "standard_name": "longitude",
            "units": "degrees_east",
            "axis": "X",
            "long_name": "Longitude",
            "climatology": "",
            "formula": "",
            "must_have_bounds": "yes",
            "out_name": "lon",
            "positive": "",
            "requested": "",
            "requested_bounds": "",
            "stored_direction": "increasing",
            "tolerance": "",
            "type": "double",
            "valid_max": "360.0",
            "valid_min": "0.0",
            "value": "",
            "z_bounds_factors": "",
            "z_factors": "",
            "bounds_values": "",
            "generic_level_name": "",
        },
    },
}
