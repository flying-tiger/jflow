#include "yaml-cpp/yaml.h"
#include <iostream>

int main(void) {
    auto lineup = YAML::Load("{1B: Prince Fielder, 2B: Rickie Weeks, LF: Ryan Braun}");
    std::cout << lineup;
    std::cin.ignore();
    return 0;
}
