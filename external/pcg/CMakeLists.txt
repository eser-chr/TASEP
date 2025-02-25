cmake_minimum_required(VERSION 3.16)
project(pcg LANGUAGES C CXX)

add_library(pcg INTERFACE)

target_include_directories(
    pcg INTERFACE
    include/
)

target_compile_features(pcg INTERFACE cxx_std_20)

#For other projects
export(TARGETS pcg FILE pcgConfig.cmake)