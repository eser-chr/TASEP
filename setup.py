# from setuptools import setup, Extension
# import pybind11
# import sysconfig
# import subprocess

# # Get the extension suffix (e.g., .so or .pyd)
# extension_suffix = sysconfig.get_config_var('EXT_SUFFIX')

# # Get the Pybind11 includes using `python3 -m pybind11 --includes`
# pybind11_includes = subprocess.check_output(
#     ['python3', '-m', 'pybind11', '--includes']
# ).decode('utf-8').strip().split()

# # Define the extension module
# extension_mod = Extension(
#     'cooperative_tasep_lib',  # Name of the generated shared library
#     sources=['cooperative_tasep.cpp', 'binding_cooperative_tasep.cpp'],
#     include_dirs=[pybind11.get_include(), '/path/to/pybind11/include'],
#     extra_compile_args=['-O3', '-Wall', '-std=c++17', '-fPIC'],
#     language='c++',
# )

# # Setup script
# setup(
#     name='cooperative_tasep_lib',
#     version='0.1',
#     ext_modules=[extension_mod],
# )


import os
import re
import sys
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import platform
from distutils.version import LooseVersion

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''), self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

setup(
    name='tasep',
    version='0.1',
    author='Christoforos Eseroglou',
    author_email='chriseseroglou@gmail.com',
    description='TASEP simulation module',
    long_description='',
    ext_modules=[CMakeExtension('tasep', sourcedir='src')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
