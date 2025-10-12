cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED On)

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


enable_testing()
add_library(SynthQ
        src/one_qubit_gate.cpp
        include/one_qubit_gate.h
        src/unitary.cpp
        include/unitary.h
        src/ctrl_qubit_gate.cpp
        include/ctrl_qubit_gate.h
        src/ucrzNode.cpp
        include/ucrzNode.h
        src/nodeVisitor.cpp
        include/nodeVisitor.h
        src/qsd.cpp
        include/qsd.h
        csd/csd.cpp
        csd/csd.h
        src/csdNode.cpp
        include/csdNode.h)

target_link_libraries(SynthQ -llapacke)

add_executable(Test
        test/test_one_qubit.cpp
        test/test_mcgates.cpp
        test/test_qsd.cpp)

target_link_libraries(
        Test
        GTest::gtest_main
        SynthQ
)

include(GoogleTest)





