import json
from os import path as op

urls = {
    "CMIP6": "https://raw.githubusercontent.com/PCMDI/cmip6-cmor-tables/master/Tables"
}


def retrieve_cmor_table(table_id, project="CMIP6"):
    import pooch

    suffix = ".json"

    url = op.join(urls[project], f"{project}_{table_id}{suffix}")

    filepath = pooch.retrieve(url, known_hash=None)

    with open(filepath) as f:
        return json.load(f)
