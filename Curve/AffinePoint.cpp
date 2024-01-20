//
// Created by eurus on 1/19/24.
//

#include "curve.cpp"

struct AffinePoint {
    union {
        struct{
            ElementFp::FpE X;
            ElementFp::FpE Y;
        };
        char bytes[64]; // Size of the union
    };
};

const AffinePoint AffinePointIdentity = {ElementFp::Zero, ElementFp::One};

static AffinePoint Neg(const AffinePoint& p)
{
    return {Negative(p.X), p.Y};
}

static AffinePoint Add(const AffinePoint& p, const AffinePoint& q)
{
    ElementFp::FpE x1Y2 = p.X * q.Y;
    ElementFp::FpE y1X2 = p.Y * q.X;
    ElementFp::FpE x1X2A = p.X * q.X * A;
    ElementFp::FpE y1Y2 = p.Y * q.Y;

    ElementFp::FpE x1X2Y1Y2D = x1Y2 * y1X2 * D;

    ElementFp::FpE xNum = x1Y2 + y1X2;

    ElementFp::FpE xDen = ElementFp::One + x1X2Y1Y2D;

    ElementFp::FpE yNum = y1Y2 - x1X2A;

    ElementFp::FpE yDen = ElementFp::One - x1X2Y1Y2D;

    ElementFp::FpE x = xNum / xDen;

    ElementFp::FpE y = yNum / yDen;

    return {x,y};
}

static AffinePoint Sub(const AffinePoint& p, const AffinePoint& q)
{
    AffinePoint neg = Neg(q);
    return Add(p, neg);
}

static AffinePoint Double(const AffinePoint& p)
{
    ElementFp::FpE xSq = p.X * p.X;
    ElementFp::FpE ySq = p.Y * p.Y;

    ElementFp::FpE xY = p.X * p.Y;
    ElementFp::FpE xSqA = xSq * A;

    ElementFp::FpE xSqYSqD = xSq * ySq * D;

    ElementFp::FpE xNum = xY + xY;

    ElementFp::FpE xDen = ElementFp::One + xSqYSqD;

    ElementFp::FpE yNum = ySq - xSqA;

    ElementFp::FpE yDen = ElementFp::One - xSqYSqD;

    ElementFp::FpE x = xNum / xDen;
    ElementFp::FpE y = yNum / yDen;

    return {x, y};
}

static bool IsOnCurve(const AffinePoint& p)
{
    ElementFp::FpE xSq = p.X * p.X;
    ElementFp::FpE ySq = p.Y * p.Y;

    ElementFp::FpE dxySq = xSq * ySq * D;
    ElementFp::FpE aXSq = A * xSq;

    ElementFp::FpE one = ElementFp::One;

    ElementFp::FpE rhs = one + dxySq;
    ElementFp::FpE lhs = aXSq + ySq;

    return Equals(lhs, rhs);
}

static AffinePoint ScalarMultiplication(const AffinePoint& point, const ElementFr::FrE& scalarMont)
{
    AffinePoint result = AffinePointIdentity;
    ElementFr::FrE scalar = ElementFr::FromMontgomery(scalarMont);
    AffinePoint pointToUse = point;

    for (int i = 0; i < scalar.BitLen(); ++i) {
        if(scalar.Bit(i)) result = Add(result, pointToUse);
        pointToUse = Double(pointToUse);
    }
    return result;
}

static ElementFp::FpE GetYCoordinate(const ElementFp::FpE& x, bool returnPositiveY)
{
    ElementFp::FpE one = ElementFp::One;
    ElementFp::FpE num = x * x;
    ElementFp::FpE den = (num * D) - one;
    num = (num * A) - one;

    ElementFp::FpE y = num / den;

    ElementFp::FpE z = {};
    if (!Sqrt(y, &z)) throw;

    bool isLargest = ElementFp::LexicographicallyLargest(z);

    return isLargest == returnPositiveY ? z : Negative(z);
}