Welcome to ``xcmor``
========================

``xcmor`` allows for lazy in-memory cmorization with `xarray <https://docs.xarray.dev/>`_.

Motivation
----------
This package aims at implementing a climate model output rewriter (cmor) tool based on xarray.
Actually, the process of *cmorization* mostly deals with making climate model output data
compliant with the Climate and Weather Forecast meta data conventions
(`CF-conventions <https://cfconventions.org/>`_) and, in most cases, does not really touch the data
itself. This makes it ideal for handling cmorization with xarray data structures since
they allow easy manipulation of meta data using python dictionaries.

While the `original cmor library <https://github.com/PCMDI/cmor>`_ offers a wide variety
of APIs for different programming languages (including python) to rewrite climate model output
to NetCDF files, xcmor focuses more on lazy cmorization in memory without neccessarily
having to actually rewrite the dataset to a filesystem. However, this is, of course, also possible
using xarrays `versatile IO features <https://docs.xarray.dev/en/stable/user-guide/io.html>`_.

Installing
----------

``xcmor`` can be installed using ``pip``

    >>> pip install xcmor


or using ``conda``

    >>> conda install -c conda-forge xcmor


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: In-depth Examples

   examples/introduction

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User Guide

   quickstart
   API Reference <api>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For contributors

   Whats New <whats-new>
   GitHub repository <https://github.com/xarray-contrib/cf-xarray>
