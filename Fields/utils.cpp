//
// Created by eurus on 1/19/24.
//

#ifdef _MSC_VER // MSVC
#include <intrin.h>
#pragma intrinsic(_BitScanReverse)
#else // GCC or Clang
#include <climits>
#endif

typedef unsigned long ulong;

static int LeadingZeroCount(ulong x)
{
#ifdef _MSC_VER // MSVC
    unsigned long index;
    if (_BitScanReverse(&index, u3) != 0)
        return static_cast<int>(CHAR_BIT * sizeof(unsigned long) - 1 - index);
#else // GCC or Clang
    if (x != 0) return __builtin_clzl(x);
#endif
    return 64;
}

static int Len64(ulong x)
{
    return 64 - LeadingZeroCount(x);
}