import numpy as np

dtype_map = {
    "real": np.dtype("float32"),
    "double": np.dtype("float64"),
}


freq_map = {"1hr": "1H", "3hr": "3H", "6hr": "6H", "day": "D", "mon": "M", "year": "Y"}

da_attrs = [
    "standard_name",
    "long_name",
    "comment",
    "units",
    "cell_methods",
    "cell_measures",
    "dimensions",
]
# da_keywords = ["out_name", "type", "dimensions"]
