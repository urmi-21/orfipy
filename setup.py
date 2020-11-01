#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:46:36 2020

@author: usingh
"""
import setuptools
#import os
import sys
#from distutils.core import setup
from distutils.extension import Extension

#try cython
try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
except:
    print("You don't have Cython installed. Please run 'pip install cython' to install cython first.")
    sys.exit(1)

#exit if python 2
if sys.version_info.major != 3:
    raise EnvironmentError("""orfipy requires python 3.5 or higher. Please upgrade your python.""")
    
#read description
with open("README.md", "r") as fh:
    long_description = fh.read()
    
    
extensions=[
    Extension("orfipy_core",
              ["orfipy/orfipy_core.pyx"]) 
]

setuptools.setup(
    name="orfipy",
    version=open("orfipy/version.py").readlines()[-1].split('=')[-1].strip().replace('"',''),
    author="Urminder Singh",
    author_email="usingh@iastate.edu",
    description="orfipy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/urmi-21/orfipy",
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
            "": []
            },
    scripts=[],
    entry_points={
            'console_scripts': [
                    'orfipy = orfipy.__main__:main'                    
                    ]
            },
    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
    tests_require=["pytest"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.5',
    cmdclass = {"build_ext": build_ext},
    ext_modules = cythonize(extensions)
)
