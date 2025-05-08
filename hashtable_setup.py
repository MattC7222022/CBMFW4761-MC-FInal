#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 21:23:47 2025

@author: mattc
"""

from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("MC_hash_table.pyx", language_level=3)
)

"""
run this to compule the cython
python hashtable_setup.py build_ext --inplace
"""