Installation
============

The pyCrust package requires `pyshtools` (>=4.7.1). This package can be installed using either `conda` or `pip`::

    conda install -c conda-forge pyshtools  # Linux and macOS only
    pip install pyshtools

After `pyshtools` is installed, you can install the pyCrust module (excluding the example scripts) using the command ::

    pip install pycrust

**Working with the example scripts**

To access the example scripts, you must download the entire pyCrust
repository from GitHub. The easiest way to do this is by cloning the repo::

    git clone https://github.com/MarkWieczorek/pycrust.git
    cd pycrust

If pyCrust is not already installed from `conda` or `pip`, it can be installed from the downloaded source using one of the two commands::

    pip install .

will install pyCrust in the active Python environment lib folder, whereas ::

    pip install -e .

will install the files in the current working directory and link them to the system Python directory. The second method is preferred if you plan on modifying the pyCrust source code.

To execute a script, it is only necessary to enter the `examples` directory and to run the file using the python command ::

    cd examples
    python pyCrust_Moon.py

.. note::
    Depending on how your system is set up, it might be necessary to use
    explicitly ``python3`` and ``pip3`` instead of ``python`` and ``pip`` in
    the above commands.
