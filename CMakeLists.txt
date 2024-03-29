#MFCFD CMAKE SCRIPT

cmake_minimum_required(VERSION 3.10)
project(mfcfd)
ENABLE_LANGUAGE(Fortran)
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if( NOT "${MFCFD}" MATCHES "CUDA" )
#------------------------------------------------------------------------------#
# PETSc support
#------------------------------------------------------------------------------#
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
set(LIBS ${LIBS} ${PETSC_LIBRARIES})

else()
#------------------------------------------------------------------------------#
# CUDA support
#------------------------------------------------------------------------------#
ENABLE_LANGUAGE(CUDA)
include_directories(${CUDA_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CUDA_LIBRARIES})

endif()

#------------------------------------------------------------------------------#
# HDF5 support
#------------------------------------------------------------------------------#

if( NOT "${MFCFD}" MATCHES "CUJOINT" )
find_package(HDF5 COMPONENTS C CXX Fortran REQUIRED)
link_directories( ${HDF5_LIBRARY_DIRS})
include_directories(${HDF5_Fortran_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_Fortran_LIBRARIES})
endif()

if( "${MFCFD}" MATCHES "PRIMAL" )
message("Building Primal MeshFree Solver")
add_compile_options(-march=native -O3 -ffree-line-length-512)
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_primal/*.F90")
add_executable(execname ${MFCFD_SOURCES})
target_link_libraries(execname ${LIBS})

elseif( "${MFCFD}" MATCHES "TANGENT" )
message("Building Tangent Meshfree Solver")
add_compile_options(-march=native -O3 -ffree-line-length-512)
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_forward/*.F90")
FILE(GLOB MFCFD_SENSE_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_forward/sensitivity_files/*.F90")
add_executable(execname ${MFCFD_SOURCES} ${MFCFD_SENSE_SOURCES})
target_link_libraries(execname ${LIBS})

elseif( "${MFCFD}" MATCHES "ADJOINT" )
message("Building Adjoint Meshfree Solver")

add_compile_options(-march=native -O3)
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_adjoint/*.F90")
FILE(GLOB MFCFD_SENSE_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_adjoint/sensitivity_files/*.F90")
set(MFCFD_AD_SOURCES ${CMAKE_SOURCE_DIR}/src_mpi_adjoint/ADFirstAidKit/adBuffer.f ${CMAKE_SOURCE_DIR}/src_mpi_adjoint/ADFirstAidKit/adStack.c ${CMAKE_SOURCE_DIR}/src_mpi_adjoint/revolve.f)
#message(${MFCFD_AD_SOURCES})
set_source_files_properties(${MFCFD_SOURCES} PROPERTIES COMPILE_FLAGS -ffree-line-length-512)
set_source_files_properties(${MFCFD_SENSE_SOURCES} PROPERTIES COMPILE_FLAGS -ffree-line-length-512)
set_source_files_properties(${CMAKE_SOURCE_DIR}/src_mpi_adjoint/revolve.f PROPERTIES COMPILE_FLAGS -fno-range-check)
set_source_files_properties(${CMAKE_SOURCE_DIR}/src_mpi_adjoint/ADFirstAidKit/adBuffer.f PROPERTIES COMPILE_FLAGS -fno-range-check)
add_executable(execname ${MFCFD_SOURCES} ${MFCFD_SENSE_SOURCES} ${MFCFD_AD_SOURCES})
target_link_libraries(execname ${LIBS})

elseif( "${MFCFD}" MATCHES "CUDA" )
message("Building Cuda Meshfree Solver")
#add_compile_options(-Mcuda=nordc,lineinfo -Mcuda=keepptx,ptxinfo)
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_cuda/*.F90")
SET(MFCFD_CUDA_SOURCES ${CMAKE_SOURCE_DIR}/src_cuda/device_data_structure.cuf ${CMAKE_SOURCE_DIR}/src_cuda/q_variables.cuf ${CMAKE_SOURCE_DIR}/src_cuda/flux_residual.cuf ${CMAKE_SOURCE_DIR}/src_cuda/state_update.cuf ${CMAKE_SOURCE_DIR}/src_cuda/q_lskum.cuf)
file(WRITE cuda.F90.in "")
foreach(MFCFD_CUDA_FILE ${MFCFD_CUDA_SOURCES})
  file(READ ${MFCFD_CUDA_FILE} CONTENTS)
  file(APPEND cuda.F90.in "${CONTENTS}")
endforeach()
configure_file(cuda.F90.in cuda.F90 COPYONLY)
set(CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcuda=nordc,lineinfo  -Minfo")
add_executable(execname ${MFCFD_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/cuda.F90)
target_link_libraries(execname ${CUDA_cusparse_LIBRARY} ${LIBS})


elseif( "${MFCFD}" MATCHES "CUJOINT" )
ENABLE_LANGUAGE(CUDA)
include_directories(${CUDA_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CUDA_LIBRARIES})
message("Building Cuda Adjoint Meshfree Solver")
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/*.F90")
FILE(GLOB MFCFD_SENSE_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/sensitivity_files/*.F90")
set(MFCFD_AD_SOURCES ${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/ADFirstAidKit/adBuffer.f ${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/ADFirstAidKit/adStack.c ${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/revolve.f)
SET(MFCFD_CUDA_SOURCES ${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/sensitivity_files/device_data_structure.cuf ${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/sensitivity_files/q_lskum_chkpts_cuda.cuf)
file(WRITE cuda.F90.in "")
foreach(MFCFD_CUDA_FILE ${MFCFD_CUDA_SOURCES})
  file(READ ${MFCFD_CUDA_FILE} CONTENTS)
  file(APPEND cuda.F90.in "${CONTENTS}")
endforeach()
configure_file(cuda.F90.in cuda.F90 COPYONLY)
set(CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcuda=nordc,lineinfo  -Minfo")
#message(${MFCFD_AD_SOURCES})
# set_source_files_properties(${MFCFD_SOURCES} PROPERTIES COMPILE_FLAGS -ffree-line-length-512)
# set_source_files_properties(${MFCFD_SENSE_SOURCES} PROPERTIES COMPILE_FLAGS -ffree-line-length-512)
# set_source_files_properties(${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/revolve.f PROPERTIES COMPILE_FLAGS -fno-range-check)
# set_source_files_properties(${CMAKE_SOURCE_DIR}/src_mpi_adjoint_cuda/ADFirstAidKit/adBuffer.f PROPERTIES COMPILE_FLAGS -fno-range-check)
add_executable(execname ${MFCFD_SOURCES} ${MFCFD_SENSE_SOURCES} ${MFCFD_AD_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/cuda.F90)
target_link_libraries(execname ${LIBS})

else()
message("Building Default Primal MeshFree Solver")
FILE(GLOB MFCFD_SOURCES "${CMAKE_SOURCE_DIR}/src_mpi_primal/*.F90")
add_executable(execname ${MFCFD_SOURCES})
target_link_libraries(execname ${LIBS})

endif()
