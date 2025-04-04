cmake_minimum_required(VERSION 3.11)
project(tasep LANGUAGES C CXX)

# -------------------------------
# 1) If you’re using scikit-build, you typically do NOT need to manually set
#    CMAKE_PREFIX_PATH to your conda environment. scikit-build sets up a
#    correct environment. If you do want to force it for some reason, do so
#    conditionally. But usually it’s best to skip this entirely.

# if(DEFINED ENV{CONDA_PREFIX})
#     set(CMAKE_PREFIX_PATH
#         "$ENV{CONDA_PREFIX}"
#         ${CMAKE_PREFIX_PATH}
#     )
# endif()
# -------------------------------

# 2) Find pybind11
find_package(pybind11 REQUIRED)

# 3) Add your subdirectories
add_subdirectory(external/pcg)
add_subdirectory(external/bucket)

# 4) Add the pybind11 module
add_library(
    tasep MODULE
    src/Bindings.cpp
    src/abstract.cpp
    src/basic_iteration.cpp
    src/kins_only.cpp
    src/neighbors.cpp
    src/nneighbors.cpp
    src/profile.cpp
)

target_include_directories(tasep PRIVATE include)
target_link_libraries(tasep PRIVATE pybind11::module bucket pcg)
target_compile_definitions(tasep PRIVATE ENABLE_CHECKS TIME_ME)
target_compile_options(tasep PRIVATE -O3 -march=native)

# Adjust suffix for Windows vs. Unix-likes
if(WIN32)
    set_target_properties(tasep PROPERTIES SUFFIX ".pyd")
else()
    set_target_properties(tasep PROPERTIES PREFIX "" SUFFIX ".so")
endif()

# 5) Use a relative install DESTINATION so scikit-build can package it
install(TARGETS tasep
    LIBRARY DESTINATION .
    RUNTIME DESTINATION .
)
