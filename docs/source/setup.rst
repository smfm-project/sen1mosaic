Setup instructions
==================

Requirements
------------

This toolset is written for use in Linux.

You will need access to a PC or server with at least:

* Python 3
* 8 GB of RAM to run SNAP.
* 8+ GB of RAM to combine images into a mosaic tile (depending on resolution/extent).

Installing Anaconda Python
--------------------------

These tools are written in Python. We recommend the Anaconda distribution of Python, which contains all the modules necessary to run these scripts.

To install Anaconda Python, open a terminal window, change directory to the location you'd like to install Anaconda Python, and run the following commands:

.. code-block:: console
    
    wget https://repo.anaconda.com/archive/Anaconda2-5.1.0-Linux-x86_64.sh
    chmod +x Anaconda2-5.1.0-Linux-x86_64.sh 
    ./Anaconda2-5.1.0-Linux-x86_64.sh 

If this has functioned, on executing ``python`` in a terminal window, you should ssee the following:

.. code-block:: console
    
    Python 2.7.14 |Anaconda, Inc.| (default, Dec  7 2017, 17:05:42) 
    [GCC 7.2.0] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> 

Setting up your Anaconda environment
------------------------------------

.. note:: The Anaconda environment required for sen1mosaic and sen2mosaic is identical. If you already have a sen2mosaic environment set up, it can be used in place of a new environment.

To ensure you are working with the appropriate version of Python as well as the correct modules, we recommend that you create an Anaconda virtual environment set up for running ``sen1mosaic``. This is done by running the following commands in your terminal or the Anaconda prompt (recommended procedure):

.. code-block:: console
    
    conda create -n sen1mosaic -c conda-forge python=3.7 scipy pandas psutil scikit-image gdal opencv

Activate the ``sen1mosaic`` environment whenever opening a new terminal window by running this command:

.. code-block:: console
    
    conda activate sen1mosaic

Installing SNAP
---------------

SNAP is an ESA toolset used for (amongst other things) pre-processing data from Sentinel-1. SNAP for Linux can be downloaded from http://step.esa.int/main/download/.

To install SNAP, open a terminal window, change directory to the location you'd like to download SNAP, and run the following commands:

.. code-block:: console

    wget http://step.esa.int/downloads/6.0/installers/esa-snap_sentinel_unix_6_0.sh
    bash esa-snap_sentinel_unix_6_0.sh
    
...and follow the instructions. The default installation instructions should work fine with sen1mosaic.

Sen1mosaic uses the 'graph processing tool' to operate SNAP commands from the command line. To give access system-wide to the graph processing tool, you'll need to add an alias to your ``.bashrc`` file as follows:

.. code-block:: console
    
    echo 'alias gpt=~/snap/bin/gpt' >> ~/.bashrc

It's a good idea to increase the memory allocation to SNAP. This is controlled by the text file ``~/snap/bin/gpt.vmoptions``. This can be done with following line:

.. code-block:: console
    
    echo '-Xmx8G' >> ~/snap/bin/gpt.vmoptions

Some SNAP operations are currently having trouble with the latest Sentinel-1 data (after March 2018). This can be fixed by installing updated through the SNAP GUI (``Help >> Check for Updates``), or with the following line in the terminal:

.. code-block:: console
    
    snap --nosplash --nogui --modules --update-all

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

To install sen1mosaic, navigate to the sen1mosaic directory and run the following within your sen1mosaic environment

.. code-block:: console
    
    python setup.py install

To avoid having to reference the full path of the Python scripts in sen1mosaic, it's a good idea add the following line to your ``.bashrc`` file:

.. code-block:: console

    echo "alias s1m='_s1m() { python ~/sen1mosaic/sen1mosaic/\"\$1\".py \$(shift; echo \"\$@\") ;}; _s1m'" >> ~/.bashrc

Installing sen2mosaic
---------------------

sen1mosaic makes use of some of the functons of sen2mosaic. To install sen2mosaic:

.. code-block:: console

    git clone https://sambowers@bitbucket.org/sambowers/sen2mosaic.git

To install sen2mosaic, navigate to the sen2mosaic directory and run the following within your sen2mosaic environment

.. code-block:: console
    
    python setup.py install

Is there a Dockerfile?
----------------------

Coming soon!

Where do I get help?
--------------------

For help installing SNAP, it's best to refer to the `ESA STEP forum <http://forum.step.esa.int/>`_. For assistance in setting up and using sen1mosaic, email `sam.bowers@ed.ac.uk <mailto:sam.bowers@ed.ac.uk>`_.

