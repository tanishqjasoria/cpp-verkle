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

    bool IsZero() const
    {
        return X.IsZero() && Equals(Y, Z) && !Y.IsZero();
    }

    AffinePoint ToAffine()
    {
        if (IsZero()) return AffinePointIdentity;
        if (Z.IsZero()) throw;
        if (Equals(Z, One)) return {X, Y};

        FpE zInv = Inverse(Z);
        FpE xAff = X * zInv;
        FpE yAff = Y * zInv;

        return {xAff, yAff};
    }
};


static ExtendedPoint FromAffine(const AffinePoint& point)
{
    return {point.X, point.Y, One};
}

const ExtendedPoint ExtendedPointIdentity = {AffinePointIdentity.X, AffinePointIdentity.Y, One};

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
    FpE x1 = p.X;
    FpE y1 = p.Y;
    FpE z1 = p.Z;

    FpE x2 = q.X;
    FpE y2 = q.Y;
    FpE z2 = q.Z;

    FpE a = z1 * z2;
    FpE b = a * a;

    FpE c = x1 * x2;

    FpE d = y1 * y2;

    FpE e = D * c * d;

    FpE f = b - e;
    FpE g = b + e;

    FpE x3 = a * f * (((x1 + y1) * (x2 + y2)) - c - d);
    FpE y3 = a * g * (d - (A * c));
    FpE z3 = f * g;

    return {x3, y3, z3};
}

static ExtendedPoint Sub(ExtendedPoint p, ExtendedPoint q)
{
    return Add(p, Neg(q));
}

static ExtendedPoint Double(ExtendedPoint p)
{
    FpE x1 = p.X;
    FpE y1 = p.Y;
    FpE z1 = p.Z;

    FpE b = (x1 + y1) * (x1 + y1);
    FpE c = x1 * x1;
    FpE d = y1 * y1;

    FpE e = A * c;
    FpE f = e + d;
    FpE h = z1 * z1;
    FpE j = f - (h + h);

    FpE x3 = (b - c - d) * j;
    FpE y3 = f * (e - d);
    FpE z3 = f * j;
    return {x3, y3, z3};
}

static ExtendedPoint ScalarMultiplication(ExtendedPoint point, FrE scalarMont) {
    ExtendedPoint result = ExtendedPointIdentity;

    FrE scalar = FromMontgomery(scalarMont);

    int len = scalar.BitLen();
    for (int i = len; i >= 0; i--) {
        result = Double(result);
        if (scalar.Bit(i)) result = Add(result, point);
    }
    return result;
}