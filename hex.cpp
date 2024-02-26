#include "hex.hpp"
#include <iostream>

int main() {
    hex::HexCoord<double> origin = hex::HexCoord<double>(1, -1, 0);

    std::cout << "Location:  " << origin.toString_Length(2, true, true) << "\n";
    std::cout << "Taxicab:   " << origin.taxicabDistance() << "\nEuclidean: " << origin.euclideanDistance() << "\n\n";

    for (unsigned char i = 0; i < 6; i++) {
        std::cout << str::toString_Length(i * 60, 3, true) << ": " << origin.getIntermediate(i).toString_Length(2, true, true) << "\n";
        std::cout << str::toString_Length(i * 60 + 30, 3, true) << ": " << origin.getCardinal(i).toString_Length(2, true, true) << "\n";
    }

    double num = 123.456;
    int lzeros = 7, tzeros = 6, length = 10;
    std::cout << "\nNum:            " << num << "\nLeading Zeros:  " << lzeros << "\nTrailing Zeros: " << tzeros << "\nLength:         " << length << "\n\n";
    std::cout << "base:       " << std::to_string(num) << "\n";
    std::cout << "regular:    " << str::toString(num, true) << "\n";
    std::cout << "lead/trail: " << str::toString_Places(num, lzeros, tzeros, false) << "\n";
    std::cout << "lead/trail: " << str::toString_Places(num, lzeros, tzeros, true) << "\n";
    std::cout << "length:     " << str::toString_Length(num, length, true) << "\n";
    std::cout << "length:     " << str::toString_Length(num, length, false) << "\n\n";
    std::wcout << "base:       " << std::to_wstring(num) << "\n";
    std::wcout << "regular:    " << str::toWideString(num) << "\n";
    std::wcout << "lead/trail: " << str::toWideString_Places(num, lzeros, tzeros, false) << "\n";
    std::wcout << "lead/trail: " << str::toWideString_Places(num, lzeros, tzeros, true) << "\n";
    std::wcout << "length:     " << str::toWideString_Length(num, length, true) << "\n";
    std::wcout << "length:     " << str::toWideString_Length(num, length, false) << "\n\n";
}
