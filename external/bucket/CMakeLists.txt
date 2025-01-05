cmake_minimum_required(VERSION 3.16)
project(bucket LANGUAGES C CXX)

add_library(bucket INTERFACE)
# target_include_directories(
#     bucket INTERFACE
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#     $<INSTALL_INTERFACE:include>
# )

target_include_directories(
    bucket INTERFACE
    include/
)

target_compile_features(bucket INTERFACE cxx_std_20)
# install(TARGETS bucket EXPORT bucket_Targets)
# install(DIRECTORY include/ DESTINATION include)


#For other projects
export(TARGETS bucket FILE bucketConfig.cmake)