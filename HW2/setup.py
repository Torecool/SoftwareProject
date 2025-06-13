from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['kmeansmodule.c', 'kmeans.c', 'vector.c'])
setup(
    name='mykmeanssp',
     version='1.0',
     description='Python wrapper for Kmeans algorithm implementation',
     ext_modules=[module]
)
