#include "hex.hpp"
#include <iostream>

int main() {
    hex::HexCoord<double> origin = hex::HexCoord<double>(1, -1, 0);

    std::cout << "Location:  " << origin.toString_Length(2, true, true) << "\n";
    std::cout << "Taxicab:   " << origin.taxicabDistance() << "\nEuclidean: " << origin.euclideanDistance() << "\n\n";

    for (unsigned char i = 0; i < 6; i++) {
        std::cout << nstrings::toString_Length(i * 60, 3, true) << ": " << origin.getDiagonal(i).toString_Length(2, true, true) << "\n";
        std::cout << nstrings::toString_Length(i * 60 + 30, 3, true) << ": " << origin.getAdjacent(i).toString_Length(2, true, true) << "\n";
    }
    origin.move(HEX_DIR_090);
    for (unsigned char i = 0; i < 12; i++) {
        std::cout << nstrings::toString_Length(i * 30, 3, true) << ": " << origin.getNeighbor(i).toString_Length(2, true, true) << "\n";
    }

    double num = 123.456;
    int lzeros = 7, tzeros = 6, length = 10;
    std::cout << "\nNum:            " << num << "\nLeading Zeros:  " << lzeros << "\nTrailing Zeros: " << tzeros << "\nLength:         " << length << "\n\n";
    std::cout << "base:       " << std::to_string(num) << "\n";
    std::cout << "regular:    " << nstrings::toString(num, true) << "\n";
    std::cout << "lead/trail: " << nstrings::toString_Places(num, lzeros, tzeros, false) << "\n";
    std::cout << "lead/trail: " << nstrings::toString_Places(num, lzeros, tzeros, true) << "\n";
    std::cout << "length:     " << nstrings::toString_Length(num, length, true) << "\n";
    std::cout << "length:     " << nstrings::toString_Length(num, length, false) << "\n\n";
    std::wcout << "base:       " << std::to_wstring(num) << "\n";
    std::wcout << "regular:    " << nstrings::toWideString(num) << "\n";
    std::wcout << "lead/trail: " << nstrings::toWideString_Places(num, lzeros, tzeros, false) << "\n";
    std::wcout << "lead/trail: " << nstrings::toWideString_Places(num, lzeros, tzeros, true) << "\n";
    std::wcout << "length:     " << nstrings::toWideString_Length(num, length, true) << "\n";
    std::wcout << "length:     " << nstrings::toWideString_Length(num, length, false) << "\n\n";
}
