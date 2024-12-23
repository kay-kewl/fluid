#include "utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

static inline std::string trim_leading_whitespace(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r");
    return (start == std::string::npos) ? "" : str.substr(start);
}

namespace utils {
    std::vector<std::string> readFieldFromFile(const char* filename) {
        std::vector<std::string> lines;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening file: " + std::string(filename));
        }
        std::string line;
        while (std::getline(file, line)) {
            lines.push_back(line);
        }
        file.close();
        return lines;
    }
}
