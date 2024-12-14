cmake_minimum_required(VERSION 3.5.0)
project(lab_9 VERSION 0.1.0 LANGUAGES C CXX)


add_executable(lab_9 main.cpp)
find_package(ZLIB REQUIRED) 
target_link_libraries(lab_9 ${ZLIB_LIBRARIES})

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
