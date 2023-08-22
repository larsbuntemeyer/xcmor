---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  name: python3
---

```{eval-rst}
.. currentmodule:: xarray
```

# What does *cmorization* mean?

The process of rewriting climate model output (cmor) to a common intercomparable data model is usually referred to as *cmorization*.
The goal is to produce [CF-compliant](https://cfconventions.org/) inter-comparable climate model output for use in
model intercomparison projects (MIPs) like, e.g., [CMIP6](https://wcrp-cmip.org/cmip-phase-6-cmip6).
Usually, the [cmor](https://cmor.llnl.gov/) library is used which provides APIs for C, Fortran and python and
is also linked to, e.g., the [cmor operator](https://code.mpimet.mpg.de/projects/cdo/wiki/CDO_CMOR_Operator)
of the climate data operator tool (CDO).

Data and metadata conventions for a MIP are usually stored in the form of json tables that contain entries for each
variable and its metadata. Metadata is attached in the form of NetCDF variable or global attributes. Additionally,
data might have to be manipulated (e.g., unit conversion or conversion of precision)

## Tables

```{code-cell}
from xcmor.resources import get_project_tables

tables = get_project_tables(project="CMIP6")
tables.coords['axis_entry'].keys()
```
