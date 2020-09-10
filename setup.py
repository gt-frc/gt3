#!/usr/bin/python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
deps = ['matplotlib', 'multiprocess', 'scipy', 'numpy', 'pandas', 'pathos', 'Shapely', 'contours', 'enum34', 'PyYAML']

setuptools.setup(
    name="GT3",
    version="0.0.1",
    author="Maxwell D. Hill, Jonathan J. Roveto, Nicholas Piper",
    install_requires=deps,
    author_email="max.hill@pm.me, veto1024@gmail.com, doom@gatech.edu",
    description="GT3 - The Georgia Tech Tokamak Transport codebase",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gt-frc/gt3/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    python_requires='>=2.7',
)

if __name__ == '__main__':
    pass
