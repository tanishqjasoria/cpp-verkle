//
// Created by eurus on 1/19/24.
//

#ifndef CPP_VERKLE_BANDERWAGON_H
#define CPP_VERKLE_BANDERWAGON_H

#include <cstdint>
#include "../Fields/FrE.h"
struct Banderwagon;

extern "C"
{
    static Banderwagon Double(Banderwagon p);
    static Banderwagon FromBytes(uint8_t* bytes, bool isBigEndian = true, bool subgroupCheck = true);
    static Banderwagon Add(const Banderwagon& p, const Banderwagon& q);
    static Banderwagon Sub(Banderwagon p, Banderwagon q);
    static Banderwagon ScalarMultiplication(Banderwagon point, FrE scalarMont);
};


#endif //CPP_VERKLE_BANDERWAGON_H
