//
// Created by eurus on 1/19/24.
//
#pragma once
#include "../Fields/FpE.cpp"
#include "../Fields/FrE.cpp"

const ElementFp::FpE A = Negative(ElementFp::SetElement(5));
const ElementFp::FpE D = {12167860994669987632UL, 4043113551995129031UL, 6052647550941614584UL, 3904213385886034240UL};

const uint8_t NumX[32] = {
        24, 174, 82, 162, 102, 24, 231, 225,
        101, 132, 153, 173, 34, 192, 121, 43,
        243, 66, 190, 123, 119, 17, 55,116,
        197, 52, 11, 44, 204, 50, 193, 41
};

const uint8_t NumY[32] = {
        102, 65, 151, 204, 182, 103, 49, 94,
        96, 100, 228, 238, 129, 173, 140,53,
        134, 213, 220, 186, 80, 139, 125,21,
        15,62, 18, 218, 158, 102, 108, 42
};

const ElementFp::FpE YTe = ElementFp::FromBytesReduced(NumX);
const ElementFp::FpE XTe = ElementFp::FromBytesReduced(NumY);