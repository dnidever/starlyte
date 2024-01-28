import os
import sys

from setuptools import setup, find_packages, Extension
from setuptools.command.build_py import build_py as _build_py

setup(name='starlyte',
      version='1.0.3',
      author = 'David Nidever',
      author_email='dnidever@montana.edu',
      url = 'https://github.com/dnidever/starlyte',
      package_dir = {'': 'python'},
      packages = ["starlyte"],
      install_requires=['numpy','astropy(>=4.0)','scipy','dlnpyutils(>=1.0.3)','ferre','gdown','chronos'],
      include_package_data=True
)
