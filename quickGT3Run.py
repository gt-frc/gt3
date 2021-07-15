#!/usr/bin/python

import GT3

shot="togt3_d3d_163477_1800"
plasma = GT3.gt3(inputFile="inputs/"+shot)
plasma.run_radial_transport()