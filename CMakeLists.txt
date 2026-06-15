cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(synq VERSION 0.0.1)


include(FetchContent)
# --------------------
# google test
# --------------------
find_package(GTest CONFIG REQUIRED)


# --------------------
# Eigen
# --------------------
find_package(Eigen3 REQUIRED NO_MODULE)

# --------------------
# MKL(WIN32) or Lapack(UNIX)
# --------------------
if (WIN32)
    add_definitions(-DUSE_MKL)
    find_package(MKL CONFIG REQUIRED)
else()
    find_package(LAPACK REQUIRED)
    find_library(LAPACKE_LIBRARY lapacke REQUIRED)
    find_path(LAPACKE_INCLUDE_DIR lapacke.h REQUIRED)
endif()

enable_testing()

add_library(synq
        src/qspUcrNode.cpp
        include/qspUcrNode.h
        include/csd.h
        src/csd.cpp
        src/one_qubit_gate.cpp
        include/one_qubit_gate.h
        src/unitary.cpp
        include/unitary.h
        src/ctrl_qubit_gate.cpp
        include/ctrl_qubit_gate.h
        src/ucrzNode.cpp
        include/ucrzNode.h
        src/ucryNode.cpp
        include/ucryNode.h
        src/iast.cpp
        include/iast.h
        src/nodeVisitor.cpp
        include/nodeVisitor.h
        src/randomUnitary.cpp
        include/randomUnitary.h
        src/qsd.cpp
        include/qsd.h
        src/ryNode.cpp
        include/ryNode.h
        include/rzNode.h
        src/rzNode.cpp
        include/ucRotationNode.h
        src/ucRotationNode.cpp
        src/csdNode.cpp
        include/csdNode.h
        src/qsdNode.cpp
        include/qsdNode.h
        src/unitaryGateNode.cpp
        include/unitaryGateNode.h
        src/twoQubitsUnitary.cpp
        include/twoQubitsUnitary.h
        src/diagonalGateNode.cpp
        include/diagonalGateNode.h
)

add_executable(Test
        test/test_one_qubit.cpp
        test/test_mcgates.cpp
        test/test_qsp.cpp
        test/test_csd.cpp
        test/test_rnd_unitary.cpp
        test/test_qsd.cpp
        test/test_unitary.cpp
        test/test_diagonal.cpp)

if(WIN32)
    target_link_libraries(synq 
        Eigen3::Eigen
        MKL::MKL)
else()
    target_link_libraries(synq 
        Eigen3::Eigen
        ${BLAS_LIBRARIES}
        ${LAPACKE_LIBRARY})
endif()

target_link_libraries(Test
        GTest::gtest_main
        synq)

include(GoogleTest)
gtest_discover_tests(Test)


