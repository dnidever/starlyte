import os
import sys

from setuptools import setup, find_packages, Extension
from setuptools.command.build_py import build_py as _build_py

setup(name='popstarlight',
      version='1.0.0',
      author = 'David Nidever',
      author_email='dnidever@montana.edu',
      url = 'https://github.com/dnidever/popstarlight',
      package_dir = {'': 'python'},
      packages = ["popstarlight"]
)
