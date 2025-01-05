# cmake_minimum_required(VERSION 3.1)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# # Project name
# project(tasep)

# # Find Python interpreter
# find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
# include_directories(${Python3_INCLUDE_DIRS})


# # Dynamically locate pybind11
# execute_process(
#     COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --cmakedir
#     OUTPUT_VARIABLE pybind11_DIR
#     OUTPUT_STRIP_TRAILING_WHITESPACE
# )
# list(APPEND CMAKE_PREFIX_PATH ${pybind11_DIR})
# find_package(pybind11 REQUIRED)

# # Add library
# add_library(
#     tasep MODULE
#     src/cooperative_tasep.cpp
#     src/binding_cooperative_tasep.cpp
#     src/basic.cpp
#     src/basiciterations.cpp
# )

# target_include_directories(tasep PRIVATE include)
# target_link_libraries(tasep PRIVATE ${Python3_LIBRARIES})
# target_link_libraries(tasep PRIVATE pybind11::module)

# # Set module properties for shared libraries
# if (WIN32)
#     set_target_properties(tasep PROPERTIES SUFFIX ".pyd")
# else()
#     set_target_properties(tasep PROPERTIES PREFIX "" SUFFIX ".so")
# endif()

# # Install to Python site-packages
# execute_process(
#     COMMAND ${PYTHON_EXECUTABLE} -c "import site; print(site.getsitepackages()[0])"
#     OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
#     OUTPUT_STRIP_TRAILING_WHITESPACE
# )
# message(STATUS "Python site-packages detected: ${PYTHON_SITE_PACKAGES}")

# # Determine the Python site-packages path
# execute_process(
#     COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('purelib'))"
#     OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
#     OUTPUT_STRIP_TRAILING_WHITESPACE
# )

# message(STATUS "Python site-packages detected: ${PYTHON_SITE_PACKAGES}")

# set(PYTHON_SITE_PACKAGES "/home/chris/Desktop/TASEP_C/venv-tasep/lib/python3.12/site-packages")


# # Install the module to the Python site-packages directory
# install(TARGETS tasep
#     LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}
# )



# FOR WORKSTATION

cmake_minimum_required(VERSION 3.1)
project(tasep LANGUAGES C CXX)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

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

add_subdirectory(external/bucket)

# Specify the source files and create the library
add_library(
    tasep MODULE
    src/Bindings.cpp
    src/abstract.cpp
    src/full.cpp
    src/kins_only.cpp
    src/neighbors.cpp
)

# Include directories
target_include_directories(tasep PRIVATE include)
# target_include_directories(tasep PRIVATE include)

# Link pybind11
target_link_libraries(tasep PRIVATE pybind11::module)
target_link_libraries(tasep PRIVATE bucket)
# target_compile_definitions(tasep PRIVATE DEBUGB ENABLE_CHECKS TIME_ME)
target_compile_options(tasep PRIVATE -O3)



# Set module properties dynamically based on the platform
if(WIN32)
    set_target_properties(tasep PROPERTIES SUFFIX ".pyd")
else()
    set_target_properties(tasep PROPERTIES PREFIX "" SUFFIX ".so")

endif()

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


