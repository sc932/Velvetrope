#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


example_module = Extension('_clubgen_c',
                           sources=['clubgen_c_wrap.c', 'clubgen_c.c'],
                           )

setup (name = 'clubgen_c',
       version = '0.16',
       author      = "Scott Clark",
       description = """Velvetrope C module""",
       ext_modules = [example_module],
       py_modules = ["clubgen_c"],
       )