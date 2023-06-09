cmake_minimum_required(VERSION 3.16.3)
project(cfem)

set(CMAKE_CXX_STANDARD 17)

###############################################
### Set your PETSc/MPI path here or bashrc  ###
### The only things to modify is the        ###
### following two lines(PETSC/MPI_DIR)      ###
###############################################
SET(WINDOWS "False" CACHE STRING "WINDOWS")

if(WINDOWS)
    message("We are running on windows system with cygwin/mingw ...")
    message(WARNING "Please use Cygwin or Mingw for the compiling !")
    if(EXISTS $ENV{PETSC_DIR})
        set(PETSC_DIR $ENV{PETSC_DIR})
        message("PETSc dir is: ${PETSC_DIR}")
        set(PETSC_INC "${PETSC_DIR}/include")
        set(PETSC_LIB "${PETSC_DIR}/lib/libpetsc.dll")
        include_directories("${PETSC_DIR}/include")
        link_libraries("${PETSC_LIB}")
        message("    PETSc include dir is: ${PETSC_INC}")
        message("    PETSc library dir is: ${PETSC_LIB}")
    else()
        message(FATAL_ERROR "PETSC_DIR is not defined in your cygwin/mingw terminal, please set the path correctly !")
    endif()
    if(EXISTS $ENV{MPI_DIR})
        message(FATAL_ERROR "MPI_DIR is not supported in cygwin/mingw terminal, please remove it !")
    endif()
elseif(UNIX)
    message ("We are running on linux system ...")
    ########################################
    ### for PETSc
    ########################################
    if(EXISTS $ENV{MPI_DIR})
        set(MPI_DIR $ENV{MPI_DIR})
        message("MPI dir is: ${MPI_DIR}")
        include_directories("${MPI_DIR}/include")
        find_library(MPI_LIB NAMES mpi mpicxx HINTS ${MPI_DIR}/lib ${MPI_DIR}/lib64)
        if(NOT MPI_LIB)
            message(FATAL_ERROR " MPI lib dir can\'t be found!")
        endif()
        message("    MPI include dir is: ${MPI_DIR}/include")
    else()
        find_path(MPI_INC NAMES mpi.h  
                  PATHS /usr/lib/x86_64-linux-gnu/openmpi/include 
                  /lib64/mpi/gcc/openmpi4/include 
                  /usr/lib64/petsc/*/linux-gnu-c-opt/include/petsc/mpiuni)
        find_library(MPI_LIB NAMES mpi 
                     PATHS /usr/lib 
                     /usr/lib64 
                     /usr/lib64/mpi/gcc/openmpi4/lib64 
                     /opt/openmpi/lib)
        message (WARNING "MPI location (MPI_DIR) is not defined in your PATH, cfem will try to use the default one installed in your system. If your mpi path is not correct, please modifiy MPI_INC and MPI_LIB in your CMakeList.txt")
        if(NOT MPI_INC)
            message(FATAL_ERROR " MPI include dir can\'t be found!")
        endif()
        if(NOT MPI_LIB)
            message(FATAL_ERROR " MPI lib dir can\'t be found!")
        endif()
        include_directories("${MPI_INC}")
        message("    MPI include dir is: ${MPI_INC}")
    endif()
    link_directories("${MPI_LIB}")
    message("    MPI library dir is: ${MPI_LIB}")
    ########################################
    ### for PETSc
    ########################################
    if(EXISTS $ENV{PETSC_DIR})
        set(PETSC_DIR $ENV{PETSC_DIR})
        message("PETSC dir is: ${PETSC_DIR}")
        include_directories("${PETSC_DIR}/include")
        find_library(PETSC_LIB NAMES petsc HINTS ${PETSC_DIR}/lib ${PETSC_DIR}/lib64)
        if(NOT PETSC_LIB)
            message(FATAL_ERROR " PETSc lib dir can\'t be found!")
        endif()
        message("    PETSc include dir is: ${PETSC_DIR}/include")
    else()
        message (WARNING "PETSc location (PETSC_DIR) is not defined in your PATH, cfem will try to use the default one installed in your system. If your PETSC_DIR path is not correct, please modifiy PETSC_INC and PETSC_LIB in your CMakeList.txt")
        find_path(PETSC_INC NAMES petsc.h 
                  PATHS /usr/lib/petsc 
                  /usr/lib/petscdir/* 
                  /opt/petsc/linux-c-opt 
                  /opt/local/lib/petsc 
                  /usr/lib64/petsc/*/linux-gnu-c-opt 
                  /usr/lib64/mpi/gcc/openmpi4/lib64/ 
                  /usr/lib64/mpi/gcc/openmpi4/lib64/petsc/*/linux-gnu-c-opt/lib)
        find_library(PETSC_LIB NAMES petsc PATHS /usr)
        if(NOT PETSC_INC)
            message(FATAL_ERROR " PETSc include dir can\'t be found!")
        endif()
        if(NOT PETSC_LIB)
            message(FATAL_ERROR " PETSc lib dir can\'t be found!")
        endif()
        message(" PETSc include dir is: ${PETSC_INC}")
        include_directories("${PETSC_INC}")
    endif()
    message("    PETSc library dir is: ${PETSC_LIB}")
    link_directories("${PETSC_LIB}")
endif()


###############################################
# For external packages                     ###
###############################################
### for eigen
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/eigen")
### for json
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external")

###############################################
### set debug or release mode               ###
###############################################
if (CMAKE_BUILD_TYPE STREQUAL "")
    # user should use -DCMAKE_BUILD_TYPE=Release[Debug] option
    set (CMAKE_BUILD_TYPE "Depdebug")
    # set (CMAKE_BUILD_TYPE "Release")
endif ()

###############################################
### set LTO for cfem                       ###
###############################################

# LTO is Link-Time Optimization.
# This block check current compiler supports whether the IPO/LTO technology.
# added by wwj(bbsy789@126.com) on 2022.12.10.

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT currentCompilerLTOsupported OUTPUT error)
endif ()

###############################################
### For linux platform                      ###
###############################################
if(UNIX)
    if (CMAKE_BUILD_TYPE STREQUAL "Depdebug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O0 -fsanitize=address -g -fno-omit-frame-pointer -std=c++14")
    elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O1 -g -std=c++14")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -g -march=native -mtune=native -DNDEBUG -std=c++14")
    else()
        message (FATAL_ERROR "Unknown cmake build type (CMAKE_BUILD_TYPE)")
    endif()
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W1 /arch:AVX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /openmp")
endif()

message("cfem will be compiled in ${CMAKE_BUILD_TYPE} mode !")
message("CXX compiler options are: ${CMAKE_CXX_FLAGS}")
message("The executable file 'cfem' will be generated in ${CMAKE_CURRENT_SOURCE_DIR}/bin folder")
###############################################
### Do not edit the following two lines !!! ###
###############################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set(src src/main.cpp)
# #############################################################
# ### For input system                                      ###
# #############################################################
set(inc ${inc} include/InputSystem/InputSystem.h)
set(src ${src} src/InputSystem/InputSystem.cpp)
#############################################################
### For message printer utils                             ###
#############################################################
set(inc ${inc} include/Utils/MessagePrinter.h include/Utils/MessageColor.h)
set(src ${src} src/Utils/MessagePrinter.cpp)
#############################################################
### For string utils                                      ###
#############################################################
set(inc ${inc} include/Utils/StringUtils.h)
set(src ${src} src/Utils/StringUtils.cpp)

#############################################################
### For json utils                                        ###
#############################################################
set(inc ${inc} include/Utils/JsonUtils.h)
set(src ${src} src/Utils/JsonUtils.cpp)

#############################################################
### For timer                                             ###
#############################################################
set(inc ${inc} include/Utils/Timer.h)
set(src ${src} src/Utils/Timer.cpp)

#############################################################
### For mathematic utils                                  ###
#############################################################
set(inc ${inc} include/MathUtils/Vector2d.h)
set(src ${src} src/MathUtils/Vector2d.cpp)
set(inc ${inc} include/MathUtils/Vector3d.h)
set(src ${src} src/MathUtils/Vector3d.cpp)
set(inc ${inc} include/MathUtils/VectorXd.h)
set(src ${src} src/MathUtils/VectorXd.cpp)
set(inc ${inc} include/MathUtils/MatrixXd.h)
set(src ${src} src/MathUtils/MatrixXd.cpp)
### for 2D Viogt tensor
set(inc ${inc} include/MathUtils/ViogtRank2Tensor2D.h)
set(src ${src} src/MathUtils/ViogtRank2Tensor2D.cpp)
set(inc ${inc} include/MathUtils/ViogtRank4Tensor2D.h)
set(src ${src} src/MathUtils/ViogtRank4Tensor2D.cpp)
set(inc ${inc} include/MathUtils/TensorConst2D.h)
set(src ${src} src/MathUtils/TensorConst2D.cpp)
### for rank-2 tensor
set(inc ${inc} include/MathUtils/Rank2Tensor3d.h)
set(src ${src} src/MathUtils/Rank2Tensor3d.cpp)
set(inc ${inc} include/MathUtils/Rank2Tensor2d.h)
set(src ${src} src/MathUtils/Rank2Tensor2d.cpp)
### for rank-4 tensor
set(inc ${inc} include/MathUtils/Rank4Tensor3d.h)
set(src ${src} src/MathUtils/Rank4Tensor3d.cpp)
# ### for general math funs
set(inc ${inc} include/MathUtils/MathFuns.h)
set(src ${src} src/MathUtils/MathFuns.cpp)
set(inc ${inc} include/MathUtils/PetsExtension.h)
set(src ${src} src/MathUtils/PetsExtension.cpp)
#############################################################
### For MeshSystem                                        ###
#############################################################
set(inc ${inc} include/MeshSystem/NodeVarInfo.h)
set(src ${src} src/MeshSystem/NodeVarInfo.cpp)
set(inc ${inc} include/MeshSystem/SetManager.h)
set(src ${src} src/MeshSystem/SetManager.cpp)
set(inc ${inc} include/MeshSystem/MeshSystem.h)
set(src ${src} src/MeshSystem/MeshSystem.cpp)
set(inc ${inc} include/MeshSystem/StructuredMesh2D.h)
set(src ${src} src/MeshSystem/StructuredMesh2D.cpp)
#############################################################
### For MaterialSystem                                    ###
#############################################################
set(inc ${inc} include/MaterialSystem/ElasticConst.h)
set(inc ${inc} include/MaterialSystem/Material.h)
set(inc ${inc} include/MaterialSystem/Material2D.h)
set(inc ${inc} include/MaterialSystem/LinearElasticMat2D.h)
set(src ${src} src/MaterialSystem/ElasticConst.cpp)
set(src ${src} src/MaterialSystem/LinearElasticMat2D.cpp)
set(inc ${inc} include/MaterialSystem/ElmtVarInfo.h)
set(src ${src} src/MaterialSystem/ElmtVarInfo.cpp)
set(inc ${inc} include/MaterialSystem/NeoHookeanAbq2d.h)
set(src ${src} src/MaterialSystem/NeoHookeanAbq2d.cpp)
#############################################################
### For ElementSystem                                     ###
#############################################################
set(inc ${inc} include/ElementSystem/ElementSystem.h)
set(src ${src} src/ElementSystem/ElementSystem.cpp)
set(inc ${inc} include/ElementSystem/Element/element.h)
set(src ${src} src/ElementSystem/Element/element.cpp)
set(inc ${inc} include/ElementSystem/Element/Element2D/CPE4R.h)
set(src ${src} src/ElementSystem/Element/Element2D/CPE4R.cpp)
set(inc ${inc} include/ElementSystem/Shpfun/Shpfun2D.h)
set(inc ${inc} include/ElementSystem/Shpfun/ShpfunQuad4.h)
set(src ${src} src/ElementSystem/Shpfun/ShpfunQuad4.cpp)
#############################################################
### For BCSSystem                                         ###
#############################################################
set(inc ${inc} include/BCsSystem/BCsSystem.h)
set(inc ${inc} include/BCsSystem/BCsSysStructured2d.h)
set(src ${src} src/BCsSystem/BCsSysStructured2d.cpp)
#############################################################
### For LoadController                                    ###
#############################################################
set(inc ${inc} include/LoadController/LoadController.h)
set(src ${src} src/LoadController/LoadController.cpp)
#############################################################
### For SolutionSystem                                    ###
#############################################################
set(inc ${inc} include/SolutionSystem/SolutionCtx.h)
set(inc ${inc} include/SolutionSystem/SolutionSystem.h)
set(src ${src} src/SolutionSystem/SolutionSystem.cpp)
set(inc ${inc} include/SolutionSystem/ArcLengthSolver.h)
set(src ${src} src/SolutionSystem/ArcLengthSolver.cpp)
#############################################################
### For system item inition                               ###
#############################################################
set(inc ${inc} include/Init/SystemInit.h)
set(src ${src} src/Init/SystemInit.cpp)
#############################################################
### For Post process system                               ###
#############################################################
set(inc ${inc} include/PostProcessSystem/OutputVarInfo.h)
set(src ${src} src/PostProcessSystem/OutputVarInfo.cpp)
set(inc ${inc} include/PostProcessSystem/PostProcessSystem.h)
set(src ${src} src/PostProcessSystem/PostProcessSystem.cpp)
set(inc ${inc} include/PostProcessSystem/PostStructured2d.h)
set(src ${src} src/PostProcessSystem/PostStructured2d.cpp)
set(src ${src} src/PostProcessSystem/OutputStructured2d.cpp)

add_executable(cfem ${inc} ${src})
target_link_libraries(cfem PUBLIC ${MPI_LIB})
target_link_libraries(cfem PUBLIC ${PETSC_LIB})

###############################################
### set LTO for cfem                        ###
###############################################

# LTO is Link-Time Optimization.
# This block check current compiler supports whether the IPO/LTO technology.
# added by wwj(bbsy789@126.com) on 2022.12.10

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(currentCompilerLTOsupported)
    message(STATUS "IPO/LTO enabled")
    set_property(TARGET cfem PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
    message(STATUS "IPO/LTO not supported: <${error}>")
    endif()
endif()
