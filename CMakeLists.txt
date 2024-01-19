cmake_minimum_required(VERSION 3.27)
project(cpp_verkle)

set(CMAKE_CXX_STANDARD 23)

add_library(cpp_verkle SHARED library.cpp
        FrE/FrE.cpp
        FrE/FrE.h
        FpE/FpE.cpp
        FpE/FpE.h
        AffinePoint/AffinePoint.cpp
        AffinePoint/AffinePoint.h
)
