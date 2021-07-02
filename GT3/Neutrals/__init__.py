#!/usr/bin/python
import pkg_resources
from packaging import version

try:
    from neutpy import neutrals
except ModuleNotFoundError:
    raise ModuleNotFoundError("Neutpy is not installed or could not be loaded. Neutrals data will be unavailable.")
try:
    if version.parse(pkg_resources.get_distribution('numpy').version) < version.parse('0.0.5'):
        raise ImportError("NeutPy v0.0.5 or newer required")
except pkg_resources.DistributionNotFound as e:
    raise

from .neutrals import *