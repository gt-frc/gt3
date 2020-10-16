The Georgia Tech Tokamak Transport Code (GT3)

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

