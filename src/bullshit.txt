cmake_minimum_required(VERSION 3.10)
project(TASEP)

set(CMAKE_CXX_STANDARD 11)

# Locate pybind11
find_package(pybind11 REQUIRED)

# Create the shared library
pybind11_add_module(tasep MODULE
    binding_cooperative_tasep.cpp
    cooperative_tasep.cpp
)
