from setuptools import Extension, setup

module = Extension("mykmeanspp", sources=['kmeansmodule.c', 'kmeans.c', 'vector.c'])
setup(
    name='mykmeanspp',
     version='1.0',
     description='Python wrapper for K-Means algorithm implementation',
     ext_modules=[module]
)
