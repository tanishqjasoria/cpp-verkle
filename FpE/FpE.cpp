//
// Created by eurus on 1/19/24.
//

#include "FpE.h"
#include <stdio.h>
#include <iostream>
#include <x86intrin.h>

typedef unsigned long ulong;

const ulong QInvNeg = 18446744069414584319UL;
const ulong Q0 = 18446744069414584321UL;
const ulong Q1 = 6034159408538082302UL;
const ulong Q2 = 3691218898639771653UL;
const ulong Q3 = 8353516859464449352UL;

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
};


//struct FpE {
//    unsigned long u0;
//    unsigned long u1;
//    unsigned long u2;
//    unsigned long u3;
//};

const FpE qElement = FpE{Q0, Q1, Q2, Q3};

const ulong R0 = 14526898881837571181UL;
const ulong R1 = 3129137299524312099UL;
const ulong R2 = 419701826671360399UL;
const ulong R3 = 524908885293268753UL;
const FpE rSquare = FpE{R0, R1, R2, R3};

static FpE SetElement(ulong u0 = 0, ulong u1 = 0, ulong u2 = 0, ulong u3 = 0)
{
    FpE x = {u0,u1,u2,u3};
    return MultiplyMod(x, rSquare);
}

ulong BigMul(ulong a, ulong b, ulong* low)
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

unsigned long AddWithCarry(unsigned long a, unsigned long b, unsigned long *carry) {
    unsigned long res = a + b + *carry;
    *carry = ((a & b) | ((a | b) & ~res)) >> 63;
    return res;
}

unsigned long MAdd0(unsigned long a, unsigned long b, unsigned long c)
{
    ulong carry = 0;
    ulong lo;
    ulong hi = BigMul(a,b,&lo);
    lo = AddWithCarry(lo, c, &carry);
    hi = hi + carry;
    return hi;
}

unsigned long MAdd1(unsigned long a, unsigned long b, unsigned long c, unsigned long *lo)
{
    ulong carry = 0;
    ulong hi = BigMul(a,b,lo);
    *lo = AddWithCarry(*lo, c, &carry);
    hi = hi + carry;
    return hi;
}

unsigned long MAdd2(unsigned long a, unsigned long b, unsigned long c, unsigned long d, unsigned long *lo)
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

unsigned long MAdd3(unsigned long a, unsigned long b, unsigned long c, unsigned long d, unsigned long e, unsigned long *lo)
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

unsigned long SubtractWithBorrow(unsigned long a, unsigned long b, unsigned long *borrow) {
    unsigned long res = a - b - *borrow;
    *borrow = ((~a & b) | (~(a ^ b) & res)) >> 63;
    return res;
}

bool LessThan(const FpE& a, const FpE& b)
{
    if (a.u3 != b.u3)
        return a.u3 < b.u3;
    if (a.u2 != b.u2)
        return a.u2 < b.u2;
    if (a.u1 != b.u1)
        return a.u1 < b.u1;
    return a.u0 < b.u0;
}

FpE RightShiftByOne(const FpE& x) {
    return FpE{
            (x.u0 >> 1) | (x.u1 << 63),
            (x.u1 >> 1) | (x.u2 << 63),
            (x.u2 >> 1) | (x.u3 << 63),
            x.u3 >> 1
    };
}


FpE Negative(const FpE& a)
{
    return SubMod(FpE{0,0,0,0}, a);
}

bool AddOverflow(const FpE& a, const FpE& b, FpE* c)
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

bool SubtractUnderflow(const FpE& a, const FpE& b, FpE* c)
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

FpE AddMod(const FpE& a, const FpE& b)
{
    FpE c{};
    AddOverflow(a,b,&c);

    if(!LessThan(c, qElement))
    {
        SubtractUnderflow(c,qElement,&c);
    }
    return c;
}

FpE SubMod(const FpE& a, const FpE& b)
{
    FpE c{};
    if (SubtractUnderflow(a, b, &c))
    {
        AddOverflow(qElement, c, &c);
    }
    return c;
}

FpE MultiplyMod(const FpE& x, const FpE& y)
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

FpE Inverse(const FpE& x)
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
