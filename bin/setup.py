from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("phyloHmm_cython.pyx")
)