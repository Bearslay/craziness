#include "hex.hpp"
#include "ncursespp.hpp"

std::wstring hexagon[6] = {
    L"  ╶──────╴  ",
    L" ╱  +000  ╲ ",
    L"╱   +000   ╲",
    L"╲   +000   ╱",
    L" ╲        ╱ ",
    L"  ╶──────╴  "
};

void printHex(hex::HexCoord<int> &coord) {
    for (unsigned char i = 0; i < 6; i++) {
        for (unsigned char j = 0; j < 12; j++) {
            if (hexagon[i][j] != L' ') {
                npp::mwin.wchar(npp::mwin.gdimy() / 2 - 3 + i + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 6 + coord.getX_npp(12, 2) + j, hexagon[i][j]);
            }
        }
        npp::mwin.wstr(npp::mwin.gdimy() / 2 - 2 + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 2 + coord.getX_npp(12, 2), str::toWideString_Length(coord.getI(), 4, true, true));
        npp::mwin.wstr(npp::mwin.gdimy() / 2 - 1 + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 2 + coord.getX_npp(12, 2), str::toWideString_Length(coord.getJ(), 4, true, true));
        npp::mwin.wstr(npp::mwin.gdimy() / 2 + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 2 + coord.getX_npp(12, 2), str::toWideString_Length(coord.getK(), 4, true, true));
    }
}

int main() {
    npp::init();

    hex::HexCoord<int> origin = hex::HexCoord<int>(0, 0, 0);
    printHex(origin);

    hex::HexCoord<int> a = hex::HexCoord<int>(3, -2, -1);
    printHex(a);
    a.rotate();
    printHex(a);
    a.rotate(-2);
    printHex(a);
    a.rotate(3);
    printHex(a);
    a.rotate();
    printHex(a);
    a.rotate();
    printHex(a);

    npp::mwin.gchar();

    hex::HexCoord<int> b = hex::HexCoord<int>(2, -2, 0);
    printHex(b);
    npp::mwin.gchar();
    b.reflectI();
    printHex(b);
    npp::mwin.gchar();
    b.reflectI();
    b.reflectJ();
    printHex(b);
    npp::mwin.gchar();
    b.reflectJ();
    b.reflectK();
    printHex(b);

    npp::mwin.gchar();
    npp::mwin.reset();

    hex::HexCoord<int> coord = hex::HexCoord<int>(0, 0, 0);

    // Algorithm for finding all hexagons within a radius around a central hexagon
    int radius = 3;
    for (int i = -radius; i <= radius; i++) {
        for (int j = std::max(-radius, -i - radius); j <= std::min(radius, -i + radius); j++) {
            hex::HexCoord<int> tempHex = coord + hex::HexCoord<int>(i, j, -i - j);
            printHex(tempHex);
        }
    }
    npp::mwin.gchar();
    npp::mwin.reset();



    npp::mwin.gchar();
    npp::end();
}
