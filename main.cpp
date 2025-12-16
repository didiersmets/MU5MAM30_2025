#include <iostream>
#include "vec3.h"
#include "stiffness.h"

int main() {
    Vec3d AB = {1,0,0};
    Vec3d AC = {0,1,0};

    double S[6];
    stiffness(AB, AC, S);

    std::cout << "Matriz de rigidez (6 entradas):\n";
    std::cout << "S00 = " << S[0] << "\n";
    std::cout << "S11 = " << S[1] << "\n";
    std::cout << "S22 = " << S[2] << "\n"; 
    std::cout << "S01 = " << S[3] << "\n";
    std::cout << "S12 = " << S[4] << "\n";
    std::cout << "S20 = " << S[5] << "\n";

    return 0;
}
