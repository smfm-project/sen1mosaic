# README #

### What is this repository for? ###

Building large-scale mosaics of Sentinel-1 data for land cover mapping is difficult, with existing tools still under-development and frequently confusing.

This is a set of tools to aid in the production of large-scale seasonal mosaic products from Sentinel-1 data. The goal is to streamline this processing chain with a set of straightforward command line tools.

This repository contains three command-line based scripts to perform the following tasks:


* Downloading Sentinel-1 data from the [Copernicus Open Access Hub](https://scihub.copernicus.eu/) for a particular latitude/longitude window, specifying date ranges and ascending/descending orbits. This is based on the [Sentinelsat](https://github.com/sentinelsat/sentinelsat) utility.
* Executing [SNAP](http://step.esa.int/main/toolboxes/snap/) graph processing tool to calibrate, filter and perform terrain-correction on Sentinel-1 GRD images.
* Mosaicking pre-processed Sentinel-1 files into larger GeoTIFF files that are suitable for image classification.

### How do I get set up? ###

These tools are written in Python for use in Linux. You will need to have first successfully installed:

* [sentinelhub](https://github.com/sinergise/sentinelhub): A library for searching and downloading Sentinel-2 products.
* [SNAP](http://step.esa.int/main/toolboxes/snap/): Tools for pre-processing data from Sentinel-1.

The modules used in these scripts are all available in [Anaconda](https://www.anaconda.com/download/) Python.

### How does it work? ###

Full documentation is hosted at: [http://sen1mosaic.readthedocs.io/](http://sen1mosaic.readthedocs.io/).

### Who do I talk to? ###

Written and maintained by Samuel Bowers ([sam.bowers@ed.ac.uk](mailto:sam.bowers@ed.ac.uk)).
