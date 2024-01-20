//
// Created by eurus on 1/20/24.
//

#ifndef CPP_VERKLE_UTILS_H
#define CPP_VERKLE_UTILS_H
#ifdef _MSC_VER // MSVC
#include <intrin.h>
#pragma intrinsic(_BitScanReverse)
#else // GCC or Clang
#include <climits>
#endif

typedef unsigned long ulong;

namespace Utils
{
    int LeadingZeroCount(ulong x);

    int Len64(ulong x);
}
#endif //CPP_VERKLE_UTILS_H
