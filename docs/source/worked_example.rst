.. _worked_example_commandline:

Worked example on the command line
==================================

Here we'll show you by example how the sen1mosaic processing chain works in practice. We will focus on an example from Zambezia Province of Mozambique, with the aim of creating a composite GeoTiff mosaic product for the province. The extent of Zambezia Province extends roughly from **35 - 40** degrees longitude and **-19 - -15** degrees latitude. We will generate a mosaic for the late dry season of 2016 (**October and Novemeber**) of **2016**, in anticipation of multiple seasonally-specific mosaics improving classification accuracy.

Preparation
-----------

First ensure that you've followed :ref:`setup` successfully.

Open a terminal, and use ``cd`` to navigate to the location you'd like to store data.

.. code-block:: console
    
    cd /home/user/DATA
    mkdir worked_example
    cd worked_example

Use mkdir to make a directory for the region to be downloaded, and navigate to that directory.

.. code-block:: console
    
    mkdir zambezia_data
    cd zambezia_data

Downloading data
----------------

The first step is to download Sentinel-1 GRD IW data from the `Copernicus Open Access Data Hub <https://scihub.copernicus.eu/>`_.

For this we use the ``download.py`` tool. We will need to specify a Scihub username and password (sign up for an account at `Scihub <https://scihub.copernicus.eu/>`_), the extent to download (in degrees lat/lon), and a start and end date in the format YYYYMMDD.

These options can be encoded as follows:

.. code-block:: console
    
    s1m download -u user.name -p supersecret -a 35 -19 40 -15 -s 20170501 -e 20170630

As we didn't specify the option ``-o`` (``--output``), data will output to the current working directory. The data format will be ``.zip``, which can be read directly by sen1mosiac without extraction.

Wait for all files to finish downloading before proceeding to the next step. By the time the processing is complete, your ``zambezia_data/`` directory should contain the following files (show files in the currenty working directory with the command ``ls``).

.. code-block:: console
    
    S1B_IW_GRDH_1SDV_20161119T155925_20161119T155950_003029_005263_2736.zip
    S1B_IW_GRDH_1SSV_20161025T030720_20161025T030745_002657_0047DF_EDDE.zip
    S1B_IW_GRDH_1SSV_20161123T031538_20161123T031603_003080_0053D4_DCC0.zip
    ...


Preprocessing Sentinel-1 GRD data
---------------------------------

The next step is to perform the preprocessing steps that convert raw Sentinel-1 data to usable terrain corrected images. We do this with the graph processing tool (gpt) bundled with SNAP.

To perform atmospheric correction and cloud masking we call the tool ``preprocess.py``. We need to specify Sentinel-1 input files, a directory containing Sentinel-1 .zip files, or a single Sentinel-1 .zip file.

To process all Sentinel-1 input files, we can submit the following line:

.. code-block:: console

    s1m preprocess ...

This command will loop through each Sentinel-2 level 1C file and process them one at a time. You might alternatively want to specify a single level 1C .SAFE file, and run several commands similtaneously. Bear in mind that this will require access to a large quantity of memory.

Here we didn't specify the options ``-o`` (``--output_dir``) and ``--g`` (``--gipp``), which can be used to output data to a location other than the directory containing input files, or the ``-r`` (``--remove``) option, which would delete Sentinel-2 level 1C data once data is finished processing.

Wait for all files to be processed to level 2A before proceeding. If you run ``ls`` again, your ``36KWA/`` directory should now contain a new set of files:

.. code-block:: console
    
    ...
    S2A_MSIL2A_20170506T074241_N0205_R049_T36KWA_20170506T075325.SAFE
    S2A_MSIL2A_20170516T072621_N0205_R049_T36KWA_20170516T075513.SAFE
    S2A_MSIL2A_20170519T075221_N0205_R092_T36KWA_20170519T080547.SAFE
    S2A_MSIL2A_20170526T074241_N0205_R049_T36KWA_20170526T074901.SAFE
    S2A_MSIL2A_20170529T073611_N0205_R092_T36KWA_20170529T075550.SAFE
    S2A_MSIL2A_20170605T072621_N0205_R049_T36KWA_20170605T075534.SAFE
    S2A_MSIL2A_20170608T075211_N0205_R092_T36KWA_20170608T080546.SAFE
    S2A_MSIL2A_20170628T075211_N0205_R092_T36KWA_20170628T080542.SAFE


Generating a cloud-free composite image
---------------------------------------

Each of these Sentinel-2 level 2A images is now atmospherically corrected, but each still contains areas of cloud. The goal of this step is to combine the cloud-free pixels of each image to generate a single cloud-free composite image. We do this with the ESA program ``sen2three``.

To perform this step we call the tool ``L3A.py``. We need to specify the directory that contains Sentinel-2 level 2A input files. Note: the code will not run if the directory contains level 2A files from multiple tiles.

To run the process, we need to submit the following line:

.. code-block:: console

    python /path/to/sen2mosaic/L3A.py /home/user/DATA/worked_example/36KWA/

Here we didn't specify the ``-r`` (``--remove``) option, which would delete Sentinel-2 level 2A data once data is finished processing.

.. warning: sen2three requires access to a lot of memory. If this is an issue, consider inputting a smaller number of level 2A fies.

Wait for sen2three to finish processing (which may take several hours). If you run ``ls`` again, your ``36KWA/`` directory should now contain a new level-3 file:

.. code-block:: console
    
    S2A_MSIL03_20170506T074241_N0205_R049_T36KWA_20160101T000000.SAFE
    
Repeat for other tiles
----------------------

The download, atmospheric correction and compositing stages need to be repeated for each tile of interest.

Now it's your turn! ``cd`` to the 36KWB folder, and generate a Sentinel-2 level-3 image using the methods we've just employed for tile 36KWA.

Generating a mosaic for classification
--------------------------------------

Once you have multiple level 3A files, the final step is to mosaic these into a larger tiling system in preparation for image classification. Whilst it is possible to classify the level 3A tiles directly, the .SAFE file format is difficult to work with, and tiles might not be the size you might prefer to work with. We recommend a grid of tiles that's approximately equal to the area of four Sentinel-2 tiles (~200,000 x 200,000 m). We call this the (unofficial) level 3B product, which is output in the easy to work with GeoTiff format.

Here we only have two tiles (36KWA and 36KWB), so we'll just perform a small-scale demonstration, generating an output with the limits **500,000 - 600,000** m Eastings and **7,550,000 - 7,650,000** m Northings (**UTM 36S**).

To perform this step we call the tool ``L3B.py``. We need to specify the location of all input files (with wildcards), the exent of the output image and the EPSG code describing the output coordinate reference system. We'll also give output data a name to identfy this tile.

First cd to the directory containing all Sentinel-2 level 3 data.

.. code-block:: console
    
    cd /home/user/DATA/worked_example/

To run ``L3B.py``, 
    
.. code-block:: console
    
    python /path/to/sen2mosaic/L3B.py -te 500000 7550000 600000 7650000 -e 32736 -n worked_example 36KW*/*_MSIL03_*.SAFE

Here we didn't specify the ``-o`` (``--output_dir``) option, meaning that results will be output to the current working directory. Once processing is complte, you can use ``ls`` to view the newly created output files:

.. code-block:: console
    
    ...

    
Viewing data
------------

In addition to a GeoTiff file for each Sentinel-2 band, ``L3B.py`` outputs two 3-band GDAL virtual dataset files (``.vrt``). These are labelled ``_RGB.vrt`` and ``_NIR.vrt``, and can be opened in QGIS to show a true colour and false colour composite (NIR, Red, Green) image.

[INSERT IMAGE]

See also
--------

This example required a lot of manual typing. We can achieve further automation through Python. To see an example of how to run achieve the same results in Python, see the scripts in the sectiomn :ref:`worked_example_python`.

