The Georgia Tech Tokamak Transport Code (GT3)

**Installation**


- **Master branch**

Clone GT3 master branch from github

`$ git clone https://github.com/gt-frc/gt3.git`

- **Other branches**

You can clone another branch from github as follows:

`$ git clone -b <branch> https://github.com/gt-frc/gt3.git`

Enter GT3

`$ cd gt3`

Setup your virtual environment

`$ virtualenv venv`

Activate it

`$ source venv/bin/activate`

Install packages

`$ pip install -r requirements.txt`

Install NeutPy

`$ cd Neutrals`

`$ git clone https://github.com/gt-frc/neutpy.git`