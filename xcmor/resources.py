import json
from os import path as op
from urllib.parse import urlparse

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
    def __init__(self, url, project=None, template=None, suffix=None):
        self.project = project
        self.url = url
        if suffix is None:
            self.suffix = ".json"
        else:
            self.suffix = suffix
        if template is None:
            self.template = f"{project}_" + "{table_id}" + ("" or self.suffix)
            print(self.template)
        else:
            self.template = template
        self.url_parsed = urlparse(url)

    def get_url(self, table_id):
        return op.join(self.url, self.template.format(table_id=table_id))

    def _retrieve(self, url):
        if self.url_parsed.scheme:
            import pooch

            return pooch.retrieve(url, known_hash=None)
        else:
            return url

    @property
    def coords(self):
        # return self._retrieve(self.get_url("coordinate"))
        return self.__getitem__("coordinate")

    @property
    def grids(self):
        return self.__getitem__("grids")

    @property
    def terms(self):
        return self.__getitem__("formula_terms")

    def __getitem__(self, key):
        with open(self._retrieve(self.get_url(key))) as f:
            return json.load(f)


cmip6 = ProjectTables(urls["CMIP6"], "CMIP6")


def get_project_tables(url=None, project=None, template=None, suffix=None):
    if project is None:
        project = "CMIP6"
    if url is None:
        url = urls[project]
    return ProjectTables(url, project, template, suffix)
