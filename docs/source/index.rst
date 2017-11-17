.. sen1mosaic documentation master file, created by
   sphinx-quickstart on Mon Nov  6 12:18:40 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to sen1mosaic's documentation!
======================================


This is a set of tools to aid in the production of large-scale seasonal mosaic products from Sentinel-1 data.

Building large-scale mosaics of Sentinel-1 data for land cover mapping is difficult, with existing tools still under-development and hard to use. The goal of these tools to streamline this processing chain with a set of straightforward command line tools.

This repository contains three command-line based scripts to perform the following tasks:

* Downloading Sentinel-1 data from the `Copernicus Open Access Hub <https://scihub.copernicus.eu/>`_ for a particular latitude/longitude window, specifying date ranges and ascending/descending orbits. This is based on the `Sentinelsat <https://github.com/sentinelsat/sentinelsat/>`_ utility.
* Executing `SNAP <http://step.esa.int/main/toolboxes/snap/>`_ graph processing tool to calibrate, filter and perform terrain-correction on Sentinel-1 GRD images.
* Mosaicking pre-processed Sentinel-1 files into larger GeoTIFF files that are suitable for image classification.

How do I get set up?
--------------------

These tools are written in Python for use in Linux. You will need to have first successfully installed the following:

* `Sentinelhub <https://github.com/sinergise/sentinelhub>`_: A library for searching and downloading Sentinel products.
* `SNAP <http://step.esa.int/main/toolboxes/snap/>`_: Pre-processing tools for Sentinel-1 data.

The modules used in these scripts are all available in Anaconda Python.

Who do I talk to?
-----------------

Written and maintained by Samuel Bowers (`sam.bowers@ed.ac.uk <mailto:sam.bowers@ed.ac.uk>`_).

Contents:

.. toctree::
   :numbered:
   :maxdepth: 2

   setup.rst
   command_line.rst
   sen1mosaic.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
