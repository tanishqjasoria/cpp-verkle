//
// Created by eurus on 1/19/24.
//

#include "Banderwagon.h"
#include "AffinePoint.cpp"

struct Banderwagon {
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

static int SubgroupCheck(const FpE& x)
{
    FpE res = x * x;
    res = res * A;
    res = Negative(res);
    res = res + One;

    return Legendre(res);
}

static Banderwagon FromBytes(uint8_t* bytes, bool isBigEndian, bool subgroupCheck)
{
    FpE x = New(bytes, isBigEndian);

    FpE y = GetYCoordinate(x, true);

    if (!subgroupCheck) return {x, y, One};
    if(SubgroupCheck(x) != 1) throw;
    return {x, y, One};
}


static Banderwagon FromAffine(const AffinePoint& point)
{
    return {point.X, point.Y, One};
}

const Banderwagon BanderwagonIdentity = {AffinePointIdentity.X, AffinePointIdentity.Y, One};

static bool Equals(const Banderwagon& p, const Banderwagon& q)
{
    if (p.IsZero()) return q.IsZero();
    if (q.IsZero()) return false;

    return Equals(p.X * q.Z, p.Z * q.X) && Equals(p.Y * q.Z, q.Y * p.Z);
}

static Banderwagon Neg(const Banderwagon& p)
{
    return {Negative(p.X), p.Y, p.Z};
}

// https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html
static Banderwagon Add(const Banderwagon& p, const Banderwagon& q)
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

static Banderwagon Sub(Banderwagon p, Banderwagon q)
{
    return Add(p, Neg(q));
}

static Banderwagon Double(Banderwagon p)
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

static Banderwagon ScalarMultiplication(Banderwagon point, FrE scalarMont) {
    Banderwagon result = BanderwagonIdentity;

    FrE scalar = FromMontgomery(scalarMont);

    int len = scalar.BitLen();
    for (int i = len; i >= 0; i--) {
        result = Double(result);
        if (scalar.Bit(i)) result = Add(result, point);
    }
    return result;
}