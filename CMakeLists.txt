cmake_minimum_required(VERSION 2.8)

project(TB_GRA Fortran)

set(CMAKE_Fortran_COMPILER      "mpif90")
#set(CMAKE_Fortran_COMPILER      "mpiifort") 

set(CMAKE_Fortran_FLAGS   "-O0")
#set(CMAKE_Fortran_FLAGS   "-O0 -fcheck=bounds -llapack")
#set(CMAKE_Fortran_FLAGS   "-O3 -ipo -xHOST")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

set(TARGET_NAME   "tb_gra")

add_subdirectory(src)

