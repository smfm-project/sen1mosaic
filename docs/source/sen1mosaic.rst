Using sen1mosaic in Python
==========================

This is harder than the command line, but you may be interested in importing the sen1mosaic functions into Python in order to customise the processing chain.

To make sen1mosaic accesible in Python, edit your ``.bashrc`` file (located at ``~/.bashrc``) to contain the line:

.. code-block:: console
    
    export PYTHONPATH=$PYTHONPATH:/path/to/sen1mosaic/

You should now be able to import each of the four modules in Python as follows:

.. code-block:: python
    
    import sen1mosaic.download
    import sen1mosaic.preprocess
    import sen1mosaic.mosaic

Help for each function can be accessed interactively from Python. For example:

.. code-block:: python
    
    >>> help(sen1mosaic.download.connectToAPI)
            Help on function connectToAPI in module sen1mosaic.download:

            connectToAPI(username, password)
            Connect to the SciHub API with sentinelsat.
            
            Args:
                username: Scihub username. Sign up at https://scihub.copernicus.eu/.
                password: Scihub password.

On this page each of the functions from the download, preprocess, and mosaic modules are documented. Note that the ``main()`` function in each is what is driven by the command line tools, so in addition to it's component parts you can call the entire processing chain from Python.

Download module
---------------

.. automodule:: sen1mosaic.download
    :members:
    :undoc-members:
    :show-inheritance:

Preprocessing module
--------------------

.. automodule:: sen1mosaic.preprocess
    :members:
    :undoc-members:
    :show-inheritance:

Mosaicking module
-----------------

.. automodule:: sen1mosaic.mosaic
    :members:
    :undoc-members:
    :show-inheritance:

