cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(synq VERSION 0.0.1)


include(FetchContent)
# --------------------
# google test
# --------------------
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


# --------------------
# Eigen
# --------------------
# Eigen (header-only)
FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
)

FetchContent_MakeAvailable(eigen)


# --------------------
# Lapack
# --------------------
find_package(LAPACK REQUIRED)
find_library(LAPACKE_LIBRARY lapacke REQUIRED)
find_path(LAPACKE_INCLUDE_DIR lapacke.h REQUIRED)

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

target_link_libraries(synq Eigen3::Eigen)

target_link_libraries(
        Test
        GTest::gtest_main
        synq
)

target_link_libraries(Test
        ${BLAS_LIBRARIES}
        ${LAPACKE_LIBRARY})

include(GoogleTest)



