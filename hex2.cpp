#include "hex.hpp"
#include "ncursespp.hpp"

std::wstring hexagon[6] = {
    L"  ╶──────╴  ",
    L" ╱  +000  ╲ ",
    L"╱   +000   ╲",
    L"╲   +000   ╱",
    L" ╲  ____  ╱ ",
    L"  ╶──────╴  "
};

void printHex(hex::HexCoord<int> coord, bool adjacent) {
    for (unsigned char i = 0; i < 6; i++) {
        for (unsigned char j = 0; j < 12; j++) {
            if (hexagon[i][j] != L' ') {
                npp::mwin.wchar(npp::mwin.gdimy() / 2 - 3 + i + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 6 + coord.getX_npp(12, 2) + j, hexagon[i][j]);
            }
        }
        npp::mwin.wstr(npp::mwin.gdimy() / 2 - 2 + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 2 + coord.getX_npp(12, 2), str::toWideString_Length(coord.getI(), 4, true, true));
        npp::mwin.wstr(npp::mwin.gdimy() / 2 - 1 + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 2 + coord.getX_npp(12, 2), str::toWideString_Length(coord.getJ(), 4, true, true));
        npp::mwin.wstr(npp::mwin.gdimy() / 2 + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 2 + coord.getX_npp(12, 2), str::toWideString_Length(coord.getK(), 4, true, true));
        npp::mwin.wstr(npp::mwin.gdimy() / 2 + 1 + coord.getY_npp(6), npp::mwin.gdimx() / 2 - 2 + coord.getX_npp(12, 2), adjacent ? L"ADJT" : L"DIAG");
    }
}

int main() {
    npp::init();

    hex::HexCoord<int> coord = hex::HexCoord<int>(1, -2, 1);

    printHex(coord, true);
    npp::mwin.gchar();

    char radius = 3;
    for (char i = -radius; i <= radius; i++) {
        for (char j = -radius - i + (i > 0 ? i : 0); j <= radius - (i > 0 ? i : 0); j++) {
            printHex(coord + hex::HexCoord<int>(i, j, -i - j), true);
        }
    }
    npp::mwin.gchar();
    npp::mwin.reset();

    printHex(coord, true);
    for (unsigned char i = 0; i < 6; i++) {
        printHex(coord.getAdjacent(i), true);
        printHex(coord.getDiagonal(i), false);
        printHex(coord.getAdjacent(i, 2), true);
        printHex(coord.getDiagonal(i, 2), false);
        printHex(coord.getAdjacent(i, 3), true);
        printHex(coord.getAdjacent(i, 4), true);
    }
    npp::mwin.gchar();

    npp::end();
}
