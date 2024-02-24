// #include "ntui.hpp"

#include <unordered_map>
#include <vector>
#include <utility>
#include <iostream>

template <typename TypeKey, typename TypeElement> std::vector<TypeKey> getKeys(std::unordered_map<TypeKey, TypeElement> input) {
    std::vector<TypeKey> output;
    for (std::pair<TypeKey, TypeElement> i : input) {
        output.emplace_back(i.first);
    }
    return output;
}
template <typename TypeKey, typename TypeElement> std::vector<TypeElement> getElements(std::unordered_map<TypeKey, TypeElement> input) {
    std::vector<TypeKey> output;
    for (std::pair<TypeKey, TypeElement> i : input) {
        output.emplace_back(i.second);
    }
    return output;
}
template <typename TypeKey, typename TypeElement> std::vector<std::pair<TypeKey, TypeElement>> mapToVector(std::unordered_map<TypeKey, TypeElement> input) {
    std::vector<std::pair<TypeKey, TypeElement>> output;
    for (std::pair<TypeKey, TypeElement> i : input) {
        output.emplace_back(i);
    }
    return output;
}

int main() {
    std::unordered_map<int, int> testing = {{1012, 25}, {-1234, 361}, {124545, 235}};
    std::vector<int> keys = getKeys(testing);
    std::vector<int> elements = getElements(testing);
    std::vector<std::pair<int, int>> full = mapToVector(testing);

    for (int i = 0; i < testing.size(); i++) {
        std::cout << keys[i] << ": " << elements[i] << " | " << full[i].first << ": " << full[i].second << "\n";
    }
}
