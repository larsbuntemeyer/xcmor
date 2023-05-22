def freq_map(cf_freq):
    """map cf frequency to pandas"""
    if cf_freq[-2:] == "hr":
        return cf_freq[0] + "H"
    else:
        return cf_freq[0].upper()


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
