
Command line tools
==================

The most straightforward way of using sen1mosaic it to call its various stages from the Linux command line. Here the functionality of each of the three commands is explained. In the next section we show how it can be used by example.

.. note:: Remember, each command line script has a help flag, which can point you in the right direction when in doubt.


Downloading Sentinel-1 data
---------------------------

Data from Sentinel-1 are available from the `Copernicus Open Access Data Hub <https://scihub.copernicus.eu/>`_, which has a graphical interface to download scenes from selected areas. Whilst useful for smaller areas, generating mosaics at national scales requires a volume of data which makes this extremely labour intensive.

The alternative is to download data using the `API Hub <https://scihub.copernicus.eu/twiki/do/view/SciHubWebPortal/APIHubDescription>`_. This system allows users to search for files using conditions on the command line, and automatically download files. To interface with the API hub, we use an open source utility called `Sentinelsat <https://sentinelsat.readthedocs.io/en/v0.12/>`_. This operates both as a command line tool, and as a Python API, which we use here. You will need to sign up for an account at `Scihub <https://scihub.copernicus.eu/>`_.

``download.py`` is a program to interface with Sentinelsat to download Sentinel-1 files, specifying a particular latitude/longitude ranges, dates and orbital directions.

Help for ``download.py`` can be viewed by typing ``s1m download --help``:

.. code-block:: console
    
    usage: download.py [-h] -u USER -p PASS -a LONMIN LATMIN LONMAX LATMAX -s
                    YYYYMMDD -e YYYYMMDD [-o PATH] [-d DIR]

    Download Sentinel-1 IW GRD data from the Copernicus Open Access Hub,
    specifying a latitude/longitude range, date ranges and ascending/descending
    orbits. Files that are already present in the destination directory won't be
    re-downloaded.

    Required arguments:
    -u USER, --user USER  Scihub username
    -p PASS, --password PASS
                            Scihub password
    -a LONMIN LATMIN LONMAX LATMAX, --search_area LONMIN LATMIN LONMAX LATMAX
                            Extent of search area, in format <lonmin latmin lonmax
                            latmax>.
    -s YYYYMMDD, --start YYYYMMDD
                            Start date for search in format YYYYMMDD.
    -e YYYYMMDD, --end YYYYMMDD
                            End date for search in format YYYYMMDD.

    Optional arguments:
    -o PATH, --output_dir PATH
                            Specify an output directory. Defaults to the current
                            working directory.
    -d DIR, --direction DIR
                            Specify <ASCENDING> or <DESCENDING> if only a single
                            orbital direction should be downloaded. Defaults to
                            downloading both.

For example, to download all data for the August-September 2017 for the longitudes of 34 to 35 degrees and latitudes of -19 to -18 degrees (around Gorongosa National Park in Mozambique), specifying an output location, use the following command:

.. code-block:: console
    
    s1m download -u user.name -p supersecret -a 34 -19 35 -18 -s 20170801 -e 20170930 -o /path/to/S1_data/

.. note:: sen1mosaic is only compatible with Sentinel-1 data in Ground Range Detected (GRD) Interferometic Wide swath (IW) mode. If you already have access to Sentinel-1 GRD IW data, you can skip straight to the next section. This may be the case if you're using a cloud platform where Sentinel-1 data archives are stored at the same location as servers.


Pre-processing Sentinel-1 data
------------------------------

Once you have Sentinel-1 (GRD IW) data, the next step is to calibrate, filter, and perform geometric correction on the data.

``preprocess.py`` takes a list of Sentinel-1 .zip files as input, and inititates a series of SNAP processing chains.

Help for ``preprocess.py`` can be viewed by typing ``s1m preprocess --help``:

.. code-block:: console
    
    usage: preprocess.py [-h] [-o PATH] [-t PATH] [-m N] [-l N] [-f] [-s] [-r]
                        [-v] [-p N]
                        [S1_FILES [S1_FILES ...]]

    Pre-process Sentinel-1 IW GRD data from the Copernicus Open Access Hub to
    radiometric/terrain corrected images.

    Optional arguments:
    S1_FILES              Input files. Specify a valid S1 input file (.zip),
                            multiple files through wildcards, or a directory.
                            Defaults to processing all S1 files in current working
                            directory.
    -o PATH, --output_dir PATH
                            Output directory for processed files. Defaults to
                            current working directory.
    -t PATH, --temp_dir PATH
                            Output directory for intermediate files. Defaults to
                            current working directory.
    -m N, --max_scenes N  Maximum number of scenes from an overpass to
                            reconstitute and process together. Higher values
                            result in fewer output files with fewer artefacts at
                            scene boundaries, but require more RAM. Defaults to 3
                            scenes.
    -l N, --multilook N   Multilooking reduces image noise by degrading output
                            resolution from ~10 x 10 m by a factor. Defaults to 2
                            (~20 x 20 m output).
    -f, --speckle_filter  Apply a speckle filter (Refined Lee) to output images.
    -s, --short           Perform a more rapid processing chain, ommitting some
                            nonessential preprocessing steps.
    -r, --remove          Delete input files after processing is complete.
    -v, --verbose         Print script progress.
    -p N, --processes N   Specify a maximum number of tiles to process in
                            paralell. Note: more processes will require more
                            resources. Defaults to 1.

For example, to run ``preprocess.py`` on a set of Sentinel-1 GRD IW .zip files in a directory (specifying an output and a temporary files directory), use the following command:

.. code-block:: console
    
    s1m preprocess -o /path/to/S1_data/ -t /scratch/ /path/to/S1_data/


Processing to GeoTiff tiles
---------------------------

The final step is to process Sentinel-1 data into a user-specified tiling grid. This script takes Sentinel-1 .dim files or a directory containing .dim files as input, selects the tiles that fall within the specified spatial extent, and mosaics available data into single-band GeoTiff files for easy use in classification systems.

``mosaic.py`` takes input .dim files and generates an output image with a specifed extent (xmin, ymin, xmax, ymax) and projection EPSG code as input. The script outputs a mean, minimum, maximum, and standard deviation of Sentinel-1 backscatter for each available polarisation.

Help for ``mosaic.py`` can be viewed by typing ``s1m mosaic --help``:

.. code-block:: console
    
    usage: mosaic.py [-h] [-te XMIN YMIN XMAX YMAX] [-e EPSG] [-r RES] [-o PATH]
                    [-n NAME] [-p POL] [-v]
                    [S1_FILES [S1_FILES ...]]

    Collate preprocessed Sentinel-1 data into mosaicked tiles. This script mosaics
    Sentinel-1 data into a customisable grid square, based on specified UTM
    coordinate bounds. Files are output as GeoTiffs of mean, min, max, and
    standard deviation of each available backscatter.

    required arguments:
    -te XMIN YMIN XMAX YMAX, --target_extent XMIN YMIN XMAX YMAX
                            Extent of output image tile, in format <xmin, ymin,
                            xmax, ymax>.
    -e EPSG, --epsg EPSG  EPSG code for output image tile CRS. This must be UTM.
                            Find the EPSG code of your output CRS as https://www
                            .epsg-registry.org/.

    optional arguments:
    S1_FILES              Input files from preprocess.py. Specify a valid S1
                            input file (.dim), multiple files through wildcards,
                            or a directory. Defaults to processing all S1 files in
                            current working directory.
    -r RES, --resolution RES
                            Output resolution in metres. Defaults to 20 m.
    -o PATH, --output_dir PATH
                            Output directory. If nothing specified, downloads will
                            output to the present working directory, given a
                            standard filename.
    -n NAME, --output_name NAME
                            Optionally specify a string to precede output
                            filename.
    -p POL, --pol POL     Specify a single polarisation ('VV' or 'VH') or
                            'both'. Defaults to processing both.
    -v, --verbose         Print script progress.


For example, to run ``mosaic.py`` in the directory ``/path/to/S1_data/`` which contains pre-processed Sentinel-1 files to create a 200 x 200 km output tile in the UTM36S projection at 20 m resolution, input:

.. code-block:: console
    
    s1m mosaic -te 600000 7900000 800000 8100000 -e 32736 -r 20 /path/to/S1_data

To do the same operation, but specifying an output directory and a name to prepend to outputs from this tile, input:

.. code-block:: console
    
    s1m mosaic -te 600000 7900000 800000 8100000 -e 32736 -r 20 -o /path/to/output/ -n tile_1 /path/to/S1_data/





