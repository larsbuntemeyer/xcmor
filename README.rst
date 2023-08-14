xcmor
=====

!:construction: UNDER CONSTRUCTION!

.. image:: https://github.com/larsbuntemeyer/xcmor/actions/workflows/ci.yaml/badge.svg
   :target: https://github.com/larsbuntemeyer/xcmor/actions/workflows/ci.yaml
   :alt: CI

.. image:: https://codecov.io/gh/larsbuntemeyer/xcmor/branch/main/graph/badge.svg?token=YW7PBUTMZ6
   :target: https://codecov.io/gh/larsbuntemeyer/xcmor

.. image:: https://results.pre-commit.ci/badge/github/larsbuntemeyer/xcmor/main.svg
   :target: https://results.pre-commit.ci/latest/github/larsbuntemeyer/xcmor/main
   :alt: pre-commit.ci status

.. image:: https://readthedocs.org/projects/xcmor/badge/?version=latest
    :target: https://xcmor.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Lazy in-memory cmorization with `xarray <https://docs.xarray.dev>`_.

Motivation
----------
This package aims at implementing a climate model output rewriter (cmor) tool based on xarray.
Actually, the process of *cmorization* mostly deals with making climate model output data
compliant with the Climate and Weather Forecast meta data conventions
(`CF-conventions <https://cfconventions.org/>`_) and, in most cases, does not really touch the data
itself. That makes is ideal for handling this process with xarray data structures since
they allow easy manipulation of meta data using python dictionaries.

While the `original cmor library <https://github.com/PCMDI/cmor>`_ offers a wide variety
of APIs for different programming languages (including python) to rewrite climate model output
to NetCDF files, xcmor focuses more on lazy cmorization in memory without neccessarily
having to actually rewrite the dataset to a filesystem. However, this is of course possible
using xarrays `versatile IO features <https://docs.xarray.dev/en/stable/user-guide/io.html>`_.
