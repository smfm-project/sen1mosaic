Setup instructions
==================

Requirements
------------

This toolset is written for use in Linux.

You will need access to a PC or server with at least:

* 8 GB of RAM to run SNAP.
* 8+ GB of RAM to combine images into a mosaic tile (depending on resolution/extent).

Installing Anaconda Python
--------------------------

These tools are written in Python. We recommend the Anaconda distribution of Python, which contains all the modules necessary to run these scripts.

To install Anaconda Python, open a terminal window, change directory to the location you'd like to install Anaconda Python, and run the following commands:

.. code-block:: console
    
    wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh
    bash Anaconda2-4.2.0-Linux-x86_64.sh


Once complete, you'll need to add this version of Python to your .bashrc file as follows:

.. code-block:: console
    
    # Substitute root for the path to your system's installation and .bashrc file.
    echo 'export PATH="~/anaconda2/bin:$PATH"' >> ~/.bashrc
    exec -l $SHELL


If this has functioned, on executing ``python`` in a terminal window, you should ssee the following:

.. code-block:: console

    Python 2.7.12 |Anaconda custom (64-bit)| (default, Jul  2 2016, 17:42:40) 
    [GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    Anaconda is brought to you by Continuum Analytics.
    Please check out: http://continuum.io/thanks and https://anaconda.org
    >>> 


Installing SNAP
---------------

SNAP is an ESA toolset used for (amongst other things) pre-processing data from Sentinel-1. SNAP for Linux can be downloaded from http://step.esa.int/main/download/.

To install SNAP, open a terminal window, change directory to the location you'd like to download SNAP, and run the following commands:

.. code-block:: console

    wget http://step.esa.int/downloads/5.0/installers/esa-snap_sentinel_unix_5_0.sh
    bash esa-snap_sentinel_unix_5_0.sh
    
...and follow the instructions. The default installation instructions should work fine with sen1mosaic.

For further details and up-to-date installation instructions, see the `SNAP website <http://step.esa.int/main/toolboxes/snap/>`_.


Installing sentinelsat
----------------------

Sentinelsat is the toolset used to access data from the Sentinel-1 archive at the `Copernicus Open Access Data Hub <https://scihub.copernicus.eu/>`_.

Up-to-date installation instructions can be found `here <https://pypi.python.org/pypi/sentinelsat>`_.

At the time of writing, the installation process is as follows:

.. code-block:: console

    pip install sentinelsat


Installing sen1mosaic
---------------------

sen1mosaic can be downloaded to a machine from its `repository<https://bitbucket.org/sambowers/sen1mosaic>`_ . To do this, open a terminal window and input:

.. code-block:: console

    git clone https://sambowers@bitbucket.org/sambowers/sen1mosaic.git


Where do I get help?
--------------------

For help installing SNAP, it's best to refer to the `ESA STEP forum <http://forum.step.esa.int/>`_. For assistance in setting up and using sen1mosaic, email `sam.bowers@ed.ac.uk <mailto:sam.bowers@ed.ac.uk>`_.

