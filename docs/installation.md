# Installation

## Dependencies

The `ctplanet` package requires `pyshtools` (>=4.13.1), which can be installed using either `conda`

```
conda install -c conda-forge pyshtools
```

or `pip`

```
pip install pyshtools
```

## Install using pip

Install the `ctplanet` module using pip

```
pip install ctplanet
```

## Install from source

First, clone the repo on your computer and cd to the new directory

```
git clone https://github.com/MarkWieczorek/ctplanet.git
cd ctplanet
```

To install `ctplanet` in the active Python environment lib folder, use

```
pip install .
```

To instead install the files in the current working directory and link them to the system Python directory, use

```
pip install -e .
```

The second method is preferred if you plan on modifying the `ctplanet` source code.

## Working with the example scripts

To access the example scripts, you must download the `ctplanet`
repository from GitHub. The easiest way to do this is by cloning the repo

```
git clone https://github.com/MarkWieczorek/ctplanet.git
```

To execute a script, it is only necessary to enter the `examples` directory and run the file using the python command

```
cd ctplanet/examples
python Moon-Crust.py
```

:::{note}
Depending on how your system is set up, it might be necessary to use
explicitly ``python3`` and ``pip3`` instead of ``python`` and ``pip`` in
the above commands.
:::
