
Command line tools
==================

The most straightforward way of using sen1mosaic it to call its various stages from the Linux command line. Here the functionality of each of the three commands is explained. In the next section we show how it can be used by example.

.. note:: Remember, each command line script has a help flag, which can be examined when in doubt.


Downloading Sentinel-1 data
---------------------------

Data from Sentinel-1 are available from the `Copernicus Open Access Data Hub <https://scihub.copernicus.eu/>`_, which has a graphical interface to download scenes from selected areas. Whilst useful for smaller areas, generating mosaics at national scales requires a volume of data which makes this extremely labour intensive.

The alternative is to download data using the `API Hub <https://scihub.copernicus.eu/twiki/do/view/SciHubWebPortal/APIHubDescription>`_. This system allows users to search for files using conditions on the command line, and automatically download files. To interface with the API hub, we use an open source utility called `Sentinelsat <https://sentinelsat.readthedocs.io/en/v0.12/>`_. This operates both as a command line tool, and as a Python API, which we use here. You will need to sign up for an account at `Scihub <https://scihub.copernicus.eu/>`_.

``download.py`` is a program to interface with Sentinelsat to download Sentinel-1 files, specifying a particular latitude/longitude ranges, dates and orbital directions.

Help for ``0_download.py`` can be viewed by typing ``s1m download --help``:

.. code-block:: console
    
    HELP TEXT

For example, to download all data for the August-September 2017 for the longitudes of 34 to 35 degrees and latitudes of -19 to -18 degrees (around Gorongosa National Park in Mozambique), specifying an output location, use the following command:

.. code-block:: console
    
    s1m download -u user.name -p supersecret -a 34 -19 35 -18 -s 20170801 -e 20170930 -o /path/to/S1_data/

.. note:: If you already have access to Sentinel-1 GRD IW data, you can skip straight to the next section. This may be the case if you're using a cloud platform where Sentinel-1 data archives are stored at the same location as servers.


Pre-processing Sentinel-1 data
------------------------------

Once you have Sentinel-1 (GRD IW) data, the next step is to calibrate, filter, and perform geometric correction on the data.

``preprocess.py`` takes a list of Sentinel-1 .zip files as input, and inititates a series of SNAP processing chains.

Help for ``preprocess.py`` can be viewed by typing ``s1m preprocess --help``:

.. code-block:: console
    
    HELP TEXT

For example, to run ``preprocess.py`` on a set of Sentinel-1 GRD IW .zip files in a directory (specifying an output and a temporary files directory), use the following command:

.. code-block:: console
    
    s1m preprocess -o /path/to/S1_data/ -t /scratch/ /path/to/S1_data/S1*_IW_GRDH_.zip

    
Processing to GeoTiff tiles
---------------------------

The final step is to process Sentinel-1 data into a user-specified tiling grid. This script takes Sentinel-1 .dim files as input, selects the tiles that fall within the specified spatial extent, and mosaics available data into single-band GeoTiff files for easy use in classification systems.

``2_tile.py`` takes a directory containing level 3A .SAFE files, an output image extent (xmin, ymin, xmax, ymax) and projection EPSG code as input.

Help for ``2_tile.py`` can be viewed by typing ``python /path/to/sen1mosaic/2_tile.py --help``:

.. code-block:: console

    HELP TEXT

For example, to run ``2_tile.py`` in the directory ``/path/to/S1_data/`` which contains pre-processed Sentinel-1 files to create a 200 x 200 km output tile in the UTM36S projection at 20 m resolution, input:

.. code-block:: console
    
    python /path/to/sen1mosaic/2_tile.py -te 600000 7900000 800000 8100000 -e 32736 -r 20 /path/to/S1_data/*.dim

To do the same operation, but specifying an output directory and a name to prepend to outputs from this tile, input:

.. code-block:: console
    
    python /path/to/sen1mosaic/2_tile.py -te 600000 7900000 800000 8100000 -e 32736 -r 20 -o /path/to/output/ -n tile_1 /path/to/S1_data/*.dim





