#include <iostream>
#include <string>
#include <stdexcept>
#include <chrono>
#include <sstream>
#include <fstream>
#include <vector>

#include "simulator.h"
#include "config.h"
#include "utils.h"
#include "macros.h"

int main(int argc, char* argv[]) {
    const char* filename = "../data/default.txt";
    const char* p_type_str = "FIXED(32,16)";
    const char* v_type_str = "FIXED(32,16)";
    const char* vf_type_str = "FIXED(32,16)";
    size_t steps = 10000;
    size_t checkpoint_interval = 1;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::cout << argv[i] << std::endl;
        if (arg == "--p-type" && i + 1 < argc) {
            p_type_str = argv[++i];
        } else if (arg == "--v-type" && i + 1 < argc) {
            v_type_str = argv[++i];
        } else if (arg == "--v-flow-type" && i + 1 < argc) {
            vf_type_str = argv[++i];
        } else if (arg == "--file" && i + 1 < argc) {
            filename = argv[++i];
        } else if (arg == "--steps" && i + 1 < argc) {
            steps = std::stoi(argv[++i]);
        } else if (arg == "--checkpoint" && i + 1 < argc) {
            checkpoint_interval = std::stoi(argv[++i]);
        }
    }

    try {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<std::string> field_data_input = utils::readFieldFromFile(filename);

        std::unique_ptr<FluidSimulatorBase> simulator = createSimulatorInstance(
            field_data_input, p_type_str, v_type_str, vf_type_str
        );
        simulator->run(steps, checkpoint_interval);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Simulation took " << duration.count() << " ms\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
