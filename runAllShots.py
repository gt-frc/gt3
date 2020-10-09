#!/usr/bin/python

from GT3 import gt3

shots = ["inputs/togt3_d3d_118888_1525",
         "inputs/togt3_d3d_118888_1570",
         "inputs/togt3_d3d_123301_2800",
         "inputs/togt3_d3d_123302_2810",
         "inputs/togt3_d3d_144977_3000",
         "inputs/togt3_d3d_164436_3720",
         "inputs/togt3_d3d_164436_3740",
         "inputs/togt3_d3d_170672_1900"]

runs = []

for shot in shots:
    plasma = gt3(inputFile=shot)
    plasma.run_radial_transport()
    runs.append(plasma)