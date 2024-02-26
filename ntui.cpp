#include "ntui.hpp"
#include <iostream>

int main() {
    ntui::MathVector<int> vector = ntui::MathVector<int>(1, 2, 3);

    std::cout << vector.geti() << vector.getj() << vector.getk() << "\n";
}
