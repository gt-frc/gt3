#!/usr/bin/python

try:
    import neutpy
    print("importing neutpy")
except:
    print("""
    #############################################
    NeutPy not found. You will be unable to run neutrals calculations without NeutPy.
    
    To install, use 
    
    $ pip install neutpy
    #############################################""")

from .gt3 import *
from ._version import __version__