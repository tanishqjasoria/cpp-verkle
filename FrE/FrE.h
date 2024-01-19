//
// Created by eurus on 1/19/24.
//

#ifndef CPP_VERKLE_FRE_H
#define CPP_VERKLE_FRE_H

struct FrE;

extern "C" {
    bool AddOverflow(const FrE& a, const FrE& b, FrE* c);
    bool SubtractUnderflow(const FrE& a, const FrE& b, FrE* c);
    FrE AddMod(const FrE& a, const FrE& b, const FrE& qElement);
    FrE SubMod(const FrE& a, const FrE& b, const FrE& qElement);
    FrE MultiplyMod(const FrE& x, const FrE& y) ;
    FrE Inverse(FrE x);
}


#endif //CPP_VERKLE_FRE_H
