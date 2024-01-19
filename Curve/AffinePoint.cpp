//
// Created by eurus on 1/19/24.
//

#include "AffinePoint.h"
#include "curve.cpp"
#include "../Fields/FpE.cpp"
#include "../Fields/FrE.cpp"


struct AffinePoint {
    union {
        struct{
            FpE X;
            FpE Y;
        };
        char bytes[64]; // Size of the union
    };
};

const AffinePoint Identity = {Zero, One};

static AffinePoint Neg(const AffinePoint& p)
{
    return {Negative(p.X), p.Y};
}

static AffinePoint Add(const AffinePoint& p, const AffinePoint& q)
{
    FpE x1Y2 = p.X * q.Y;
    FpE y1X2 = p.Y * q.X;
    FpE x1X2A = p.X * q.X * A;
    FpE y1Y2 = p.Y * q.Y;

    FpE x1X2Y1Y2D = x1Y2 * y1X2 * D;

    FpE xNum = x1Y2 + y1X2;

    FpE xDen = One + x1X2Y1Y2D;

    FpE yNum = y1Y2 - x1X2A;

    FpE yDen = One - x1X2Y1Y2D;

    FpE x = xNum / xDen;

    FpE y = yNum / yDen;

    return {x,y};
}

static AffinePoint Sub(const AffinePoint& p, const AffinePoint& q)
{
    AffinePoint neg = Neg(q);
    return Add(p, neg);
}

static AffinePoint Double(const AffinePoint& p)
{
    FpE xSq = p.X * p.X;
    FpE ySq = p.Y * p.Y;

    FpE xY = p.X * p.Y;
    FpE xSqA = xSq * A;

    FpE xSqYSqD = xSq * ySq * D;

    FpE xNum = xY + xY;

    FpE xDen = One + xSqYSqD;

    FpE yNum = ySq - xSqA;

    FpE yDen = One - xSqYSqD;

    FpE x = xNum / xDen;
    FpE y = yNum / yDen;

    return {x, y};
}

static bool IsOnCurve(const AffinePoint& p)
{
    FpE xSq = p.X * p.X;
    FpE ySq = p.Y * p.Y;

    FpE dxySq = xSq * ySq * D;
    FpE aXSq = A * xSq;

    FpE one = One;

    FpE rhs = one + dxySq;
    FpE lhs = aXSq + ySq;

    return Equals(lhs, rhs);
}

static AffinePoint ScalarMultiplication(const AffinePoint& point, const FrE& scalarMont)
{
    AffinePoint result = Identity;
    FrE scalar = FromMontgomery(scalarMont);
    AffinePoint pointToUse = point;

    for (int i = 0; i < scalar.BitLen(); ++i) {
        if(scalar.Bit(i)) result = Add(result, pointToUse);
        pointToUse = Double(pointToUse);
    }
    return result;
}

static FpE GetYCoordinate(const FpE& x, bool returnPositiveY)
{
    FpE one = One;
    FpE num = x * x;
    FpE den = (num * D) - one;
    num = (num * A) - one;

    FpE y = num / den;

    FpE z = {};
    if (!Sqrt(y, &z)) return {};

    bool isLargest = LexicographicallyLargest(z);

    return isLargest == returnPositiveY ? z : Negative(z);
}