//
// Created by eurus on 1/19/24.
//

#include "FpE.h"
#include <stdio.h>
#include <iostream>
#include <x86intrin.h>
#include "utils.h"


namespace ElementFp
{

    struct FpE {
        union {
            struct {
                unsigned long u0;
                unsigned long u1;
                unsigned long u2;
                unsigned long u3;
            };
            char bytes[32]; // Size of the union
        };

        static FpE AddMod(const FpE& a, const FpE& b);
        static FpE SubMod(const FpE& a, const FpE& b);
        static FpE MultiplyMod(const FpE& x, const FpE& y);
        static FpE Inverse(const FpE& x);

        FpE operator*(const FpE& obj) const
        {
            return MultiplyMod(*this, obj);
        }

        FpE operator+(const FpE& obj) const
        {
            return AddMod(*this, obj);
        }

        FpE operator-(const FpE& obj) const
        {
            return SubMod(*this, obj);
        }

        FpE operator/(const FpE& obj) const
        {
            return MultiplyMod(*this, Inverse(obj));
        }

        int BitLen() const
        {
            return u3 != 0
                   ? 192 + Utils::Len64(u3)
                   : u2 != 0
                     ? 128 + Utils::Len64(u2)
                     : u1 != 0
                       ? 64 + Utils::Len64(u1)
                       : Utils::Len64(u0);
        }

        bool Bit(int n) const
        {
            int bucket = n / 64 % 4;
            int position = n % 64;
            ulong bucketValue;
            switch (bucket) {
                case 0:
                    bucketValue = this->u0;
                    break;
                case 1:
                    bucketValue = this->u1;
                    break;
                case 2:
                    bucketValue = this->u2;
                    break;
                case 3:
                    bucketValue = this->u3;
                    break;
            }
            return (bucketValue & ((ulong)1 << position)) != 0;
        }

        bool IsZero() const
        {
            return (u0 | u1 | u2 | u3) == 0;
        }

    };


//struct FpE {
//    unsigned long u0;
//    unsigned long u1;
//    unsigned long u2;
//    unsigned long u3;
//};

    const ulong QInvNeg = 18446744069414584319UL;
    const ulong Q0 = 18446744069414584321UL;
    const ulong Q1 = 6034159408538082302UL;
    const ulong Q2 = 3691218898639771653UL;
    const ulong Q3 = 8353516859464449352UL;
    const FpE qElement = FpE{Q0, Q1, Q2, Q3};

    const ulong R0 = 14526898881837571181UL;
    const ulong R1 = 3129137299524312099UL;
    const ulong R2 = 419701826671360399UL;
    const ulong R3 = 524908885293268753UL;
    const FpE rSquare = FpE{R0, R1, R2, R3};

    const FpE Zero = {0};

    const ulong One0 = 8589934590UL;
    const ulong One1 = 6378425256633387010UL;
    const ulong One2 = 11064306276430008309UL;
    const ulong One3 = 1739710354780652911UL;
    const FpE One = {One0, One1, One2, One3};

    const ulong QM0 = 9223372034707292161UL;
    const ulong QM1 = 12240451741123816959UL;
    const ulong QM2 = 1845609449319885826UL;
    const ulong QM3 = 4176758429732224676UL;
    const FpE qMinOne = {QM0, QM1, QM2, QM3};

    const ulong G0 = 11289237133041595516UL;
    const ulong G1 = 2081200955273736677UL;
    const ulong G2 = 967625415375836421UL;
    const ulong G3 = 4543825880697944938UL;
    const FpE gResidue = {G0, G1, G2, G3};

    const ulong BS0 = 9223141137265459199UL;
    const ulong BS1 = 347036667491570177UL;
    const ulong BS2 = 10722717374829358084UL;
    const ulong BS3 = 972477353UL;
    const FpE _bSqrtExponentElement = {BS0, BS1, BS2, BS3};

    const ulong BL0 = 9223372034707292160UL;
    const ulong BL1 = 12240451741123816959UL;
    const ulong BL2 = 1845609449319885826UL;
    const ulong BL3 = 4176758429732224676UL;
    const FpE _bLegendreExponentElement = {BL0, BL1, BL2, BL3};


    const ulong SqrtR = 32;

    static bool LessThan(const FpE& a, const FpE& b)
    {
        if (a.u3 != b.u3)
            return a.u3 < b.u3;
        if (a.u2 != b.u2)
            return a.u2 < b.u2;
        if (a.u1 != b.u1)
            return a.u1 < b.u1;
        return a.u0 < b.u0;
    }

    static FpE RightShiftByOne(const FpE& x) {
        return FpE{
                (x.u0 >> 1) | (x.u1 << 63),
                (x.u1 >> 1) | (x.u2 << 63),
                (x.u2 >> 1) | (x.u3 << 63),
                x.u3 >> 1
        };
    }

    static bool Equals(const FpE& lhs, const FpE& rhs)
    {
        return lhs.u0 == rhs.u0 && lhs.u1 == rhs.u1 && lhs.u2 == rhs.u2 && lhs.u3 == rhs.u3;
    }

    static unsigned long AddWithCarry(unsigned long a, unsigned long b, unsigned long *carry) {
        unsigned long res = a + b + *carry;
        *carry = ((a & b) | ((a | b) & ~res)) >> 63;
        return res;
    }

    static unsigned long SubtractWithBorrow(unsigned long a, unsigned long b, unsigned long *borrow) {
        unsigned long res = a - b - *borrow;
        *borrow = ((~a & b) | (~(a ^ b) & res)) >> 63;
        return res;
    }

    static ulong BigMul(ulong a, ulong b, ulong* low)
    {
        uint num1 = (uint) a;
        uint num2 = (uint) (a >> 32);
        uint num3 = (uint) b;
        uint num4 = (uint) (b >> 32);
        ulong num5 = (ulong) num1 * (ulong) num3;
        ulong num6 = (ulong) num2 * (ulong) num3 + (num5 >> 32);
        ulong num7 = (ulong) num1 * (ulong) num4 + (ulong) (uint) num6;
        *low = num7 << 32 | (ulong) (uint) num5;
        return (ulong) num2 * (ulong) num4 + (num6 >> 32) + (num7 >> 32);
    }


    static unsigned long MAdd0(unsigned long a, unsigned long b, unsigned long c)
    {
        ulong carry = 0;
        ulong lo;
        ulong hi = BigMul(a,b,&lo);
        lo = AddWithCarry(lo, c, &carry);
        hi = hi + carry;
        return hi;
    }

    static unsigned long MAdd1(unsigned long a, unsigned long b, unsigned long c, unsigned long *lo)
    {
        ulong carry = 0;
        ulong hi = BigMul(a,b,lo);
        *lo = AddWithCarry(*lo, c, &carry);
        hi = hi + carry;
        return hi;
    }

    static unsigned long MAdd2(unsigned long a, unsigned long b, unsigned long c, unsigned long d, unsigned long *lo)
    {
        ulong carry = 0;
        ulong hi = BigMul(a,b,lo);
        c = AddWithCarry(c, d, &carry);
        hi = hi + carry;
        carry = 0;
        *lo = AddWithCarry(*lo, c, &carry);
        hi = hi + carry;
        return hi;
    }

    static unsigned long MAdd3(unsigned long a, unsigned long b, unsigned long c, unsigned long d, unsigned long e, unsigned long *lo)
    {
        ulong carry = 0;
        ulong hi = BigMul(a,b,lo);
        c = AddWithCarry(c, d, &carry);
        hi = hi + carry;
        carry = 0;
        *lo = AddWithCarry(*lo, c, &carry);
        hi = hi + e + carry;
        return hi;
    }
    static bool AddOverflow(const FpE& a, const FpE& b, FpE* c)
    {
        FpE p{};
        unsigned long carry = 0;
        p.u0 = AddWithCarry(a.u0, b.u0, &carry);
        p.u1 = AddWithCarry(a.u1, b.u1, &carry);
        p.u2 = AddWithCarry(a.u2, b.u2, &carry);
        p.u3 = AddWithCarry(a.u3, b.u3, &carry);
        *c = p;
        return carry != 0;
    }

    static bool SubtractUnderflow(const FpE& a, const FpE& b, FpE* c)
    {
        FpE p{};
        unsigned long borrow = 0;
        p.u0 = SubtractWithBorrow(a.u0, b.u0, &borrow);
        p.u1 = SubtractWithBorrow(a.u1, b.u1, &borrow);
        p.u2 = SubtractWithBorrow(a.u2, b.u2, &borrow);
        p.u3 = SubtractWithBorrow(a.u3, b.u3, &borrow);
        *c = p;
        return borrow != 0;
    }

    FpE FpE::AddMod(const FpE& a, const FpE& b)
    {
        FpE c{};
        AddOverflow(a,b,&c);

        if(!LessThan(c, qElement))
        {
            SubtractUnderflow(c,qElement,&c);
        }
        return c;
    }

    FpE FpE::SubMod(const FpE& a, const FpE& b)
    {
        FpE c{};
        if (SubtractUnderflow(a, b, &c))
        {
            AddOverflow(qElement, c, &c);
        }
        return c;
    }

    FpE FpE::MultiplyMod(const FpE& x, const FpE& y)
    {
        ulong t[4] = {0,0,0,0};
        ulong c[3] = {0,0,0};
        ulong z[4] = {0,0,0,0};

        // round 0
        c[1] = BigMul(x.u0, y.u0, & c[0]);
        ulong m = c[0] * QInvNeg;
        c[2] = MAdd0(m, Q0, c[0]);
        c[1] = MAdd1(x.u0, y.u1, c[1], & c[0]);
        c[2] = MAdd2(m, Q1, c[2], c[0], & t[0]);
        c[1] = MAdd1(x.u0, y.u2, c[1], & c[0]);
        c[2] = MAdd2(m, Q2, c[2], c[0], & t[1]);
        c[1] = MAdd1(x.u0, y.u3, c[1], & c[0]);
        t[3] = MAdd3(m, Q3, c[0], c[2], c[1], & t[2]);

        // round 1
        c[1] = MAdd1(x.u1, y.u0, t[0], & c[0]);
        m = c[0] * QInvNeg;
        c[2] = MAdd0(m, Q0, c[0]);
        c[1] = MAdd2(x.u1, y.u1, c[1], t[1], & c[0]);
        c[2] = MAdd2(m, Q1, c[2], c[0], & t[0]);
        c[1] = MAdd2(x.u1, y.u2, c[1], t[2], & c[0]);
        c[2] = MAdd2(m, Q2, c[2], c[0], & t[1]);
        c[1] = MAdd2(x.u1, y.u3, c[1], t[3], & c[0]);
        t[3] = MAdd3(m, Q3, c[0], c[2], c[1], & t[2]);

        // round 2
        c[1] = MAdd1(x.u2, y.u0, t[0], & c[0]);
        m = c[0] * QInvNeg;
        c[2] = MAdd0(m, Q0, c[0]);
        c[1] = MAdd2(x.u2, y.u1, c[1], t[1], & c[0]);
        c[2] = MAdd2(m, Q1, c[2], c[0], & t[0]);
        c[1] = MAdd2(x.u2, y.u2, c[1], t[2], & c[0]);
        c[2] = MAdd2(m, Q2, c[2], c[0], & t[1]);
        c[1] = MAdd2(x.u2, y.u3, c[1], t[3], & c[0]);
        t[3] = MAdd3(m, Q3, c[0], c[2], c[1], & t[2]);

        // round 3
        c[1] = MAdd1(x.u3, y.u0, t[0], & c[0]);
        m = c[0] * QInvNeg;
        c[2] = MAdd0(m, Q0, c[0]);
        c[1] = MAdd2(x.u3, y.u1, c[1], t[1], & c[0]);
        c[2] = MAdd2(m, Q1, c[2], c[0], & z[0]);
        c[1] = MAdd2(x.u3, y.u2, c[1], t[2], & c[0]);
        c[2] = MAdd2(m, Q2, c[2], c[0], & z[1]);
        c[1] = MAdd2(x.u3, y.u3, c[1], t[3], & c[0]);
        z[3] = MAdd3(m, Q3, c[0], c[2], c[1], & z[2]);


        FpE res = {z[0],z[1],z[2],z[3]};
        if (LessThan(qElement, res)) SubtractUnderflow(res, qElement, & res);
        return res;
    }

    static FpE SetElement(ulong u0 = 0, ulong u1 = 0, ulong u2 = 0, ulong u3 = 0)
    {
        FpE x = {u0,u1,u2,u3};
        return FpE::MultiplyMod(x, rSquare);
    }


    static FpE Mod(const FpE& x) {
        FpE modElem = x;
        if (LessThan(qElement, x))
        {
            if (SubtractUnderflow(x, qElement, & modElem)) throw;
        }
        return modElem;
    }

    static FpE Negative(const FpE& a)
    {
        return FpE::SubMod(FpE{0,0,0,0}, a);
    }

    static FpE New(const uint8_t* bytes, bool isBigEndian = false)
    {
        FpE x = FpE{0,0,0,0};
        if (isBigEndian) {
            x.u3 = (static_cast<uint64_t>(bytes[0]) << 56) | (static_cast<uint64_t>(bytes[1]) << 48) |
                   (static_cast<uint64_t>(bytes[2]) << 40) | (static_cast<uint64_t>(bytes[3]) << 32) |
                   (static_cast<uint64_t>(bytes[4]) << 24) | (static_cast<uint64_t>(bytes[5]) << 16) |
                   (static_cast<uint64_t>(bytes[6]) << 8) | static_cast<uint64_t>(bytes[7]);

            x.u2 = (static_cast<uint64_t>(bytes[8]) << 56) | (static_cast<uint64_t>(bytes[9]) << 48) |
                   (static_cast<uint64_t>(bytes[10]) << 40) | (static_cast<uint64_t>(bytes[11]) << 32) |
                   (static_cast<uint64_t>(bytes[12]) << 24) | (static_cast<uint64_t>(bytes[13]) << 16) |
                   (static_cast<uint64_t>(bytes[14]) << 8) | static_cast<uint64_t>(bytes[15]);

            x.u1 = (static_cast<uint64_t>(bytes[16]) << 56) | (static_cast<uint64_t>(bytes[17]) << 48) |
                   (static_cast<uint64_t>(bytes[18]) << 40) | (static_cast<uint64_t>(bytes[19]) << 32) |
                   (static_cast<uint64_t>(bytes[20]) << 24) | (static_cast<uint64_t>(bytes[21]) << 16) |
                   (static_cast<uint64_t>(bytes[22]) << 8) | static_cast<uint64_t>(bytes[23]);

            x.u0 = (static_cast<uint64_t>(bytes[24]) << 56) | (static_cast<uint64_t>(bytes[25]) << 48) |
                   (static_cast<uint64_t>(bytes[26]) << 40) | (static_cast<uint64_t>(bytes[27]) << 32) |
                   (static_cast<uint64_t>(bytes[28]) << 24) | (static_cast<uint64_t>(bytes[29]) << 16) |
                   (static_cast<uint64_t>(bytes[30]) << 8) | static_cast<uint64_t>(bytes[31]);
        } else {
            x.u0 = (static_cast<uint64_t>(bytes[0]) << 56) | (static_cast<uint64_t>(bytes[1]) << 48) |
                   (static_cast<uint64_t>(bytes[2]) << 40) | (static_cast<uint64_t>(bytes[3]) << 32) |
                   (static_cast<uint64_t>(bytes[4]) << 24) | (static_cast<uint64_t>(bytes[5]) << 16) |
                   (static_cast<uint64_t>(bytes[6]) << 8) | static_cast<uint64_t>(bytes[7]);

            x.u1 = (static_cast<uint64_t>(bytes[8]) << 56) | (static_cast<uint64_t>(bytes[9]) << 48) |
                   (static_cast<uint64_t>(bytes[10]) << 40) | (static_cast<uint64_t>(bytes[11]) << 32) |
                   (static_cast<uint64_t>(bytes[12]) << 24) | (static_cast<uint64_t>(bytes[13]) << 16) |
                   (static_cast<uint64_t>(bytes[14]) << 8) | static_cast<uint64_t>(bytes[15]);

            x.u2 = (static_cast<uint64_t>(bytes[16]) << 56) | (static_cast<uint64_t>(bytes[17]) << 48) |
                   (static_cast<uint64_t>(bytes[18]) << 40) | (static_cast<uint64_t>(bytes[19]) << 32) |
                   (static_cast<uint64_t>(bytes[20]) << 24) | (static_cast<uint64_t>(bytes[21]) << 16) |
                   (static_cast<uint64_t>(bytes[22]) << 8) | static_cast<uint64_t>(bytes[23]);

            x.u3 = (static_cast<uint64_t>(bytes[24]) << 56) | (static_cast<uint64_t>(bytes[25]) << 48) |
                   (static_cast<uint64_t>(bytes[26]) << 40) | (static_cast<uint64_t>(bytes[27]) << 32) |
                   (static_cast<uint64_t>(bytes[28]) << 24) | (static_cast<uint64_t>(bytes[29]) << 16) |
                   (static_cast<uint64_t>(bytes[30]) << 8) | static_cast<uint64_t>(bytes[31]);
        }
        return x;
    }

    static FpE FromBytesReduced(const uint8_t* bytes, bool isBigEndian = false)
    {
        FpE res = New(bytes, isBigEndian);
        res = Mod(res);
        return FpE::MultiplyMod(res, rSquare);
    }

    static FpE FromMontgomery(const FpE& x)
    {
        ulong z[4] = {x.u0,x.u1,x.u2,x.u3};

        ulong m = z[0] * QInvNeg;
        ulong c = MAdd0(m, Q0, z[0]);
        c = MAdd2(m, Q1, z[1], c, & z[0]);
        c = MAdd2(m, Q2, z[2], c, & z[1]);
        c = MAdd2(m, Q3, z[3], c, & z[2]);
        z[3] = c;

        m = z[0] * QInvNeg;
        c = MAdd0(m, Q0, z[0]);
        c = MAdd2(m, Q1, z[1], c, & z[0]);
        c = MAdd2(m, Q2, z[2], c, & z[1]);
        c = MAdd2(m, Q3, z[3], c, & z[2]);
        z[3] = c;

        m = z[0] * QInvNeg;
        c = MAdd0(m, Q0, z[0]);
        c = MAdd2(m, Q1, z[1], c, & z[0]);
        c = MAdd2(m, Q2, z[2], c, & z[1]);
        c = MAdd2(m, Q3, z[3], c, & z[2]);
        z[3] = c;

        m = z[0] * QInvNeg;
        c = MAdd0(m, Q0, z[0]);
        c = MAdd2(m, Q1, z[1], c, & z[0]);
        c = MAdd2(m, Q2, z[2], c, & z[1]);
        c = MAdd2(m, Q3, z[3], c, & z[2]);
        z[3] = c;

        FpE res = {z[0],z[1],z[2],z[3]};
        if (LessThan(qElement, res)) SubtractUnderflow(res, qElement, & res);
        return res;
    }








    FpE FpE::Inverse(const FpE& x)
    {
        // initialize u = q
        FpE u = qElement;
        // initialize s = r^2
        FpE s = rSquare;
        FpE r = {0,0,0,0};
        FpE v = x;


        while (true)
        {
            while ((v.u0 & 1) == 0)
            {
                v = RightShiftByOne(v);
                if ((s.u0 & 1) == 1) AddOverflow(s, qElement, & s);

                s = RightShiftByOne(s);
            }

            while ((u.u0 & 1) == 0)
            {
                u = RightShiftByOne(u);
                if ((r.u0 & 1) == 1) AddOverflow(r, qElement, & r);

                r = RightShiftByOne(r);
            }

            if (!LessThan(v, u))
            {
                SubtractUnderflow(v, u, & v);
                s = SubMod(s, r);
            }
            else
            {
                SubtractUnderflow(u, v, & u);
                r = SubMod(r, s);
            }


            if (u.u0 == 1 && (u.u3 | u.u2 | u.u1) == 0)
            {
                return r;
            }

            if (v.u0 == 1 && (v.u3 | v.u2 | v.u1) == 0)
            {
                return s;
            }
        }
    }

    static FpE Exp(const FpE& b, const FpE& e)
    {
        FpE result = One;
        FpE bs = b;
        int len = e.BitLen();
        for (int i = 0; i < len; i++)
        {
            if (e.Bit(i)) result = FpE::MultiplyMod(result, bs);

            bs = FpE::MultiplyMod(bs, bs);
        }
        return result;
    }

    static int Legendre(const FpE& z)
    {
        FpE res = Exp(z, _bLegendreExponentElement);

        if (res.IsZero()) return 0;
        if (Equals(res, One)) return 1;

        return -1;
    }

    static bool Sqrt(const FpE& x, FpE* z)
    {
        FpE w = Exp(x, _bSqrtExponentElement);

        FpE y = x * w;
        FpE b = w * y;

        ulong r = SqrtR;

        FpE t = b;

        for (ulong i = 0; i < r-1; ++i) {
            t = t * t;
        }

        if(t.IsZero())
        {
            *z = Zero;
            return true;
        }

        if(!Equals(t, One))
        {
            *z = Zero;
            return false;
        }

        FpE g = gResidue;

        while (true)
        {
            ulong m = 0;
            t = b;

            while (!Equals(t, One))
            {
                t = t * t;
                m++;
            }

            if (m == 0)
            {
                *z = y;
                return true;
            }

            // t = g^(2^(r-m-1)) (mod q)
            int ge = (int)(r - m - 1);
            t = g;
            while (ge > 0)
            {
                t = t * t;
                ge--;
            }

            g = t * t;
            y = t * y;
            b = b * g;

            r = m;
        }

    }

    static bool LexicographicallyLargest(FpE x)
    {
        FpE mont = FromMontgomery(x);
        FpE y = {};
        return !SubtractUnderflow(mont, qMinOne, &y);
    }

};

