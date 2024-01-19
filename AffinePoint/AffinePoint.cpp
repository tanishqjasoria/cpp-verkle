//
// Created by eurus on 1/19/24.
//

#include "AffinePoint.h"
#include "../FpE/FpE.cpp"

struct AffinePoint {
    union {
        struct{
            FpE X;
            FpE Y;
        };
        char bytes[64]; // Size of the union
    };
};