//
// Created by eurus on 1/19/24.
//

#ifndef CPP_VERKLE_FRE_H
#define CPP_VERKLE_FRE_H

struct FrE;

extern "C" {
    static bool AddOverflow(const FrE& a, const FrE& b, FrE* c);
static bool SubtractUnderflow(const FrE& a, const FrE& b, FrE* c);
static  FrE AddMod(const FrE& a, const FrE& b);
static  FrE SubMod(const FrE& a, const FrE& b);
static FrE MultiplyMod(const FrE& x, const FrE& y) ;
static  FrE Inverse(FrE x);
}


#endif //CPP_VERKLE_FRE_H
