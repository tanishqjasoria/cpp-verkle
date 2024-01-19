//
// Created by eurus on 1/19/24.
//

#ifndef CPP_VERKLE_FPE_H
#define CPP_VERKLE_FPE_H


struct FpE;

extern "C" {
    bool AddOverflow(const FpE& a, const FpE& b, FpE* c);
    bool SubtractUnderflow(const FpE& a, const FpE& b, FpE* c);
    FpE AddMod(const FpE& a, const FpE& b);
    FpE SubMod(const FpE& a, const FpE& b);
    FpE MultiplyMod(const FpE& x, const FpE& y) ;
    FpE Inverse(FpE x);
}


#endif //CPP_VERKLE_FPE_H
