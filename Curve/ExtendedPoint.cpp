//
// Created by eurus on 1/19/24.
//

#include "ExtendedPoint.h"
#include "AffinePoint.cpp"

struct ExtendedPoint {
    union {
        struct{
            ElementFp::FpE X;
            ElementFp::FpE Y;
            ElementFp::FpE Z;
        };
        char bytes[96]; // Size of the union
    };

    bool IsZero() const
    {
        return X.IsZero() && Equals(Y, Z) && !Y.IsZero();
    }

    AffinePoint ToAffine()
    {
        if (IsZero()) return AffinePointIdentity;
        if (Z.IsZero()) throw;
        if (Equals(Z, ElementFp::One)) return {X, Y};

        ElementFp::FpE zInv = ElementFp::FpE::Inverse(Z);
        ElementFp::FpE xAff = X * zInv;
        ElementFp::FpE yAff = Y * zInv;

        return {xAff, yAff};
    }
};


static ExtendedPoint FromAffine(const AffinePoint& point)
{
    return {point.X, point.Y, ElementFp::One};
}

const ExtendedPoint ExtendedPointIdentity = {AffinePointIdentity.X, AffinePointIdentity.Y, ElementFp::One};

static bool Equals(const ExtendedPoint& p, const ExtendedPoint& q)
{
    if (p.IsZero()) return q.IsZero();
    if (q.IsZero()) return false;

    return Equals(p.X * q.Z, p.Z * q.X) && Equals(p.Y * q.Z, q.Y * p.Z);
}

static ExtendedPoint Neg(const ExtendedPoint& p)
{
    return {Negative(p.X), p.Y, p.Z};
}

// https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html
static ExtendedPoint Add(const ExtendedPoint& p, const ExtendedPoint& q)
{
    ElementFp::FpE x1 = p.X;
    ElementFp::FpE y1 = p.Y;
    ElementFp::FpE z1 = p.Z;

    ElementFp::FpE x2 = q.X;
    ElementFp::FpE y2 = q.Y;
    ElementFp::FpE z2 = q.Z;

    ElementFp::FpE a = z1 * z2;
    ElementFp::FpE b = a * a;

    ElementFp::FpE c = x1 * x2;

    ElementFp::FpE d = y1 * y2;

    ElementFp::FpE e = D * c * d;

    ElementFp::FpE f = b - e;
    ElementFp::FpE g = b + e;

    ElementFp::FpE x3 = a * f * (((x1 + y1) * (x2 + y2)) - c - d);
    ElementFp::FpE y3 = a * g * (d - (A * c));
    ElementFp::FpE z3 = f * g;

    return {x3, y3, z3};
}

static ExtendedPoint Sub(ExtendedPoint p, ExtendedPoint q)
{
    return Add(p, Neg(q));
}

static ExtendedPoint Double(ExtendedPoint p)
{
    ElementFp::FpE x1 = p.X;
    ElementFp::FpE y1 = p.Y;
    ElementFp::FpE z1 = p.Z;

    ElementFp::FpE b = (x1 + y1) * (x1 + y1);
    ElementFp::FpE c = x1 * x1;
    ElementFp::FpE d = y1 * y1;

    ElementFp::FpE e = A * c;
    ElementFp::FpE f = e + d;
    ElementFp::FpE h = z1 * z1;
    ElementFp::FpE j = f - (h + h);

    ElementFp::FpE x3 = (b - c - d) * j;
    ElementFp::FpE y3 = f * (e - d);
    ElementFp::FpE z3 = f * j;
    return {x3, y3, z3};
}

static ExtendedPoint ScalarMultiplication(ExtendedPoint point, ElementFr::FrE scalarMont) {
    ExtendedPoint result = ExtendedPointIdentity;

    ElementFr::FrE scalar = ElementFr::FromMontgomery(scalarMont);

    int len = scalar.BitLen();
    for (int i = len; i >= 0; i--) {
        result = Double(result);
        if (scalar.Bit(i)) result = Add(result, point);
    }
    return result;
}