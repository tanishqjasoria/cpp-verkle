//
// Created by eurus on 1/19/24.
//

#include "ExtendedPoint.h"
#include "../Fields/FpE.cpp"
#include "../Fields/FrE.cpp"
#include "AffinePoint.cpp"

struct ExtendedPoint {
    union {
        struct{
            FpE X;
            FpE Y;
            FpE Z;
        };
        char bytes[96]; // Size of the union
    };
};


static ExtendedPoint FromAffine(AffinePoint point)
{
    return {point.X, point.Y, One};
}


