cmake_minimum_required(VERSION 3.10)
project(Com_toolbox)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-w)

include_directories(${CMAKE_SOURCE_DIR}/dependencies/external/matplotlib-cpp)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
include_directories(${NUMPY_INCLUDE_DIR})

file(GLOB INCLUDE_DIR ${PROJECT_SOURCE_DIR}/dependencies/include)
file(GLOB LIB_DIR ${PROJECT_SOURCE_DIR}/dependencies/lib)
file(GLOB TEST_DIR ${PROJECT_SOURCE_DIR}/dependencies/Test)

include_directories(${INCLUDE_DIR} ${TEST_DIR} ${LIB_DIR})

file(GLOB SRC_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB LIB_SOURCES "${PROJECT_SOURCE_DIR}/dependencies/lib/*.cpp")
file(GLOB TEST_SOURCES "${PROJECT_SOURCE_DIR}/dependencies/Test/*.cpp")

# Create an executable target
add_executable(Com_toolbox ${SRC_SOURCES} ${TEST_SOURCES} ${LIB_SOURCES})


target_link_libraries(${PROJECT_NAME} Python3::Python Python3::Module)

