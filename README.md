The Georgia Tech Tokamak Transport Code (GT3). Also known as GTEDGE2.

**Installation**

**Using pip**

`$ pip install GT3`

**Installing from GitHub**

- **Master branch**

Clone GT3 master branch from github

`$ git clone https://github.com/gt-frc/gt3.git`

- **Other branches**

You can clone another branch from github as follows:

`$ git clone -b <branch> https://github.com/gt-frc/gt3.git`

**Installing Branches from GitHub via pip**

To use pip to install a development version of GT3, simply use

`$ pip install git+https://github.com/gt-frc/gt3@development`

**Neutrals Calculations**

GT3 utilizes NeutPy (https://github.com/gt-frc/neutpy) for the calculation of 
neutral particles recycling from the wall. Neutpy is an optional dependency. NeutPy
also requires the Triangle 2D meshing package (see the NeutPy github for details).

To install NeutPy:

`$ pip install neutpy`

**Usage**

GT3 can be run in 2 different ways: via a config file and profiles or via a test data class.

- **From Files**

The primary method of running GT3 is via the gt3 class

```python
from GT3 import gt3
```

If `neutpy` is installed, you'll see a message indicating that it is being imported. Otherwise, a warning 
will be given. If you will be running `neutpy`, the main configuration file needs to be placed in the
CWD. A sample configuration file can be found on the neutpy github as `neutpy.conf`

GT3 is instantiated in this manner by providing an input file.

```python
myPlasma = gt3(inputFile="myInputFile")
```

See the `/inputs` directory for example input files (generally called `togt3_d3d_shotid_timeid`).
The `inputFile` argument takes a file relative to the current working directory (CWD). An input file will
include some plasma parameters and meshing information. In addition, the locations of 1D and 2D
profile data are entered into this file. The locations also must be relative to the CWD. See the
`inputs` directory for an example with DIII-D shot 144977.3000.

GT3 includes various modules that provide information and calculations about the plasma. Running `gt3()`
with the `mode` argument will run different sets of modules (see `gt3.py` for a list of modes).
Modules can also be run using the various arguments (`run_IOL()`, `run_NBI()`, `run_radial_transport()`, etc.).

To run the full radial transport code, use the `mode='radialtrans'` argument in `gt3()` or run 
`myPlasma.run_radial_transport()`.

- **From A Class**

*Documentation Coming Soon*