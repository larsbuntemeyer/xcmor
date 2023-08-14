from pprint import pprint

from xcmor import Cmorizer
from xcmor.datasets import reg_ds

cmor = Cmorizer()
ds_out = cmor.cmorize(
    reg_ds.rename(temperature="tas").tas, "Amon", cmor.tables["input_example"]
)

print(ds_out)
pprint(ds_out.tas.attrs)
pprint(ds_out.attrs)
ds_out.to_netcdf("test.nc")
