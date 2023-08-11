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


class ProjectTables:
    def __init__(self, table_url, project=None, template=None, suffix=".json"):
        self.project = project
        self.table_url = table_url
        self.suffix = suffix
        if template is None:
            self.template = f"{project}_" + "{table_id}"
        else:
            self.template = template

    def _get_url(self, table_id):
        return op.join(
            self.table_url, self.template.format(table_id=table_id) + self.suffix
        )

    def _retrieve(self, url):
        import pooch

        return pooch.retrieve(url, known_hash=None)

    @property
    def coords(self):
        return self._retrieve(self._get_url("coordinate"))

    @property
    def grids(self):
        return self._retrieve(self._get_url("grids"))

    @property
    def terms(self):
        return self._retrieve(self._get_url("formula_terms"))

    def __getattr__(self, table):
        return self._retrieve(self._get_url(table))

    def __getitem__(self, key):
        return self.__getattr__(key)


cmip6 = ProjectTables(urls["CMIP6"], "CMIP6")
