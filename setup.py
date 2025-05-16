from setuptools import setup, Extension
import pybind11

module = Extension(
    's_curve',
    sources=['src/s_curve_py.cpp', 'src/s_curve.cpp'],
    include_dirs=[pybind11.get_include(), 'include'],
    language='c++',
    extra_compile_args=['-std=c++11'],
)

setup(
    name='s_curve',
    version='0.1',
    ext_modules=[module],
)