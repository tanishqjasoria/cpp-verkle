cmake_minimum_required(VERSION 3.27)
project(cpp_verkle)

set(CMAKE_CXX_STANDARD 23)

add_library(cpp_verkle SHARED library.cpp
        Fields/FrE.cpp
        Fields/FrE.h
        Fields/FpE.cpp
        Fields/FpE.h
        Curve/AffinePoint.cpp
        Curve/AffinePoint.h
        Curve/ExtendedPoint.cpp
        Curve/ExtendedPoint.h
        Fields/FrE.cpp
        Fields/FrE.h
        Fields/FpE.cpp
        Fields/FpE.h
        Fields/utils.cpp
        Fields/utils.cpp
        Curve/curve.cpp
        Curve/Banderwagon.cpp
        Curve/Banderwagon.h
)
