//
// Created by eurus on 1/19/24.
//

#ifndef CPP_VERKLE_FPE_H
#define CPP_VERKLE_FPE_H


struct FpE;

extern "C" {
    static bool AddOverflow(const FpE& a, const FpE& b, FpE* c);
    static bool SubtractUnderflow(const FpE& a, const FpE& b, FpE* c);
    static FpE AddMod(const FpE& a, const FpE& b);
    static FpE SubMod(const FpE& a, const FpE& b);
    static FpE MultiplyMod(const FpE& x, const FpE& y) ;
    static FpE Inverse(const FpE& x);
}


#endif //CPP_VERKLE_FPE_H
