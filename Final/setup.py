from setuptools import Extension, setup
import numpy

module = Extension(
    "symnmfmodule",
    sources=['symnmfmodule.c', 'symnmf.c', 'vector.c', 'matrix.c'],
    include_dirs=[numpy.get_include()]
)

setup(
    name='symnmfmodule',
     version='1.0',
     description='Python wrapper for SymNMF algorithm implementation',
     ext_modules=[module]
)
