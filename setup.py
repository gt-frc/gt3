#!/usr/bin/python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
deps = ['matplotlib',
        'multiprocess',
        'scipy',
        'numpy',
        'pandas',
        'pathos',
        'Shapely',
        'contours',
        'enum34',
        'PyYAML',
        'texttable',
        'deprecation']

setuptools.setup(
    name="GT3",
    version=open("GT3/_version.py").readlines()[-1].split()[-1].strip("\"'"),
    author="Maxwell D. Hill, Jonathan J. Roveto, Nicholas Piper",
    install_requires=deps,
    include_package_data=True,
    author_email="max.hill@pm.me, veto1024@gmail.com, doom@gatech.edu",
    description="GT3 - The Georgia Tech Tokamak Transport codebase",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gt-frc/gt3/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    python_requires='>=3.8',
)

if __name__ == '__main__':
    pass
