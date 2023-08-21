# What does *cmorization* mean?

The process of rewriting climate model output (cmor) to a common intercomparable data model is usually referred to as *cmorization*.
The goal is to produce [CF-compliant](https://cfconventions.org/) inter-comparable climate model output for use in
model intercomparison projects like, e.g., [CMIP6](https://wcrp-cmip.org/cmip-phase-6-cmip6).
Usually, the [cmor](https://cmor.llnl.gov/) library is used which provides APIs for C, Fortran and python and
is also linked to, e.g., the [cmor operator](https://code.mpimet.mpg.de/projects/cdo/wiki/CDO_CMOR_Operator)
of the climate data operator tool (CDO).
However, the cmorization process usually requires an actual rewrite of the data
