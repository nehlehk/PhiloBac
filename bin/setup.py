from setuptools import setup,Extension
from Cython.Build import cythonize
import numpy


# setup(
#     ext_modules=cythonize("PhiloBacter.pyx")
#    ,include_path='/home/nehleh/anaconda3/lib/python3.8/site-packages/numpy/core/include/numpy'
#    ,include_dirs=[numpy.get_include()]
#    ,compiler_directives={'language_level' : "3"}
# )

cython_directives = {'language_level': "3"}

setup(
   ext_modules=cythonize(
                    Extension(
                             name='PhiloBacter'
                            ,sources=["PhiloBacter.pyx"]
                            ,include_dirs=[numpy.get_include()]
                            # ,language='c'
                               )
                            # ,include_path=["/home/nehleh/anaconda3/lib/python3.8/site-packages/numpy/core/include/numpy"]
                            ,include_path=["/home/nehleh/anaconda3/envs/PhiloBacteria/lib/python3.8/site-packages/numpy/core/include/numpy"]
                            # ,compiler_directives={'language_level' : "3"}
                            ,compiler_directives= cython_directives
                            )
)

