cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED On)

IF (WIN32)
include_directories(c:/include)
ENDIF()

project(SynthQ VERSION 0.0.1)
project(Test)
# GoogleTest requires at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

if (${APPLE})
        set(LAPACKE_DIR "/opt/homebrew/Cellar/lapack/3.12.1/")
        include_directories(${LAPACKE_DIR}/include
                /opt/homebrew/include/eigen3)
        link_directories(${LAPACKE_DIR}/lib)
endif()

enable_testing()

add_library(SynthQ
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
        )

add_executable(Test
        test/test_one_qubit.cpp
        test/test_mcgates.cpp
        test/test_qsp.cpp
        test/test_csd.cpp)

target_link_libraries(
        Test
        GTest::gtest_main
        SynthQ
        lapacke
)

include(GoogleTest)



