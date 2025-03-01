cmake_minimum_required(VERSION 3.11)
project(tasep LANGUAGES C CXX)

# Find Python interpreter and extract version information
find_package(PythonInterp REQUIRED)
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
    OUTPUT_VARIABLE PYTHON_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Check if CONDA_PREFIX is set (meaning we're in a Conda environment)
if(DEFINED ENV{CONDA_PREFIX})
    # Use the Conda environment's path to locate pybind11
    set(CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}/lib/python${PYTHON_VERSION}/site-packages/pybind11/share/cmake/pybind11" ${CMAKE_PREFIX_PATH})
    message(${CMAKE_PREFIX_PATH})
endif()

# Find pybind11
find_package(pybind11 REQUIRED)

add_subdirectory(external/pcg)
add_subdirectory(external/bucket)

# Specify the source files and create the library
add_library(
    tasep MODULE
    src/Bindings.cpp
    src/abstract.cpp
    src/basic_iteration.cpp
    src/phase.cpp
    src/kins_only.cpp
    src/neighbors.cpp
    src/nneighbors.cpp
)

# add_library(
#     neighbors MODULE
#     src/abstract.cpp
#     src/neighbors.cpp
#     python/neighbors.cpp
# )
# add_library(
#     tasep_profile MODULE
#     src/abstract.cpp
#     src/phase.cpp
#     python/profile.cpp
# )

# Include directories
target_include_directories(tasep PRIVATE include)
# target_include_directories(neighbors PRIVATE include)
# target_include_directories(tasep_profile PRIVATE include)


# Link pybind11
target_link_libraries(tasep PRIVATE pybind11::module)
target_link_libraries(tasep PRIVATE bucket)
target_link_libraries(tasep PRIVATE pcg)

# target_link_libraries(neighbors PRIVATE pybind11::module)
# target_link_libraries(neighbors PRIVATE bucket)

# target_link_libraries(tasep_profile PRIVATE pybind11::module)
# target_link_libraries(tasep_profile PRIVATE bucket)

target_compile_definitions(tasep PRIVATE ENABLE_CHECKS TIME_ME)
target_compile_options(tasep PRIVATE -O3 -march=native)
# target_compile_options(neighbors PRIVATE -O3)
# target_compile_options(tasep_profile PRIVATE -O3)



# Set module properties dynamically based on the platform
if(WIN32)
    set_target_properties(tasep PROPERTIES SUFFIX ".pyd")
else()
    set_target_properties(tasep PROPERTIES PREFIX "" SUFFIX ".so")

endif()

# if(WIN32)
#     set_target_properties(tasep neighbors PROPERTIES SUFFIX ".pyd")
# else()
#     set_target_properties(tasep neighbors PROPERTIES PREFIX "" SUFFIX ".so")
# endif()
# if(WIN32)
#     set_target_properties(tasep neighbors tasep_profile PROPERTIES SUFFIX ".pyd")
# else()
#     set_target_properties(tasep neighbors tasep_profile PROPERTIES PREFIX "" SUFFIX ".so")
# endif()


# Install the module to the Python site-packages directory dynamically
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)



message(STATUS "Python site-packages: ${PYTHON_SITE_PACKAGES}")


install(TARGETS tasep
    LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}
)
# install(TARGETS neighbors
#     LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}
# )
# install(TARGETS tasep_profile
#     LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}
# )

