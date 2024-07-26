from setuptools import setup, Extension
import pybind11
import sysconfig
import subprocess

# Get the extension suffix (e.g., .so or .pyd)
extension_suffix = sysconfig.get_config_var('EXT_SUFFIX')

# Get the Pybind11 includes using `python3 -m pybind11 --includes`
pybind11_includes = subprocess.check_output(
    ['python3', '-m', 'pybind11', '--includes']
).decode('utf-8').strip().split()

# Define the extension module
extension_mod = Extension(
    'cooperative_tasep_lib',  # Name of the generated shared library
    sources=['cooperative_tasep.cpp', 'binding_cooperative_tasep.cpp'],
    include_dirs=[pybind11.get_include(), '/path/to/pybind11/include'],
    extra_compile_args=['-O3', '-Wall', '-std=c++17', '-fPIC'],
    language='c++',
)

# Setup script
setup(
    name='cooperative_tasep_lib',
    version='0.1',
    ext_modules=[extension_mod],
)
