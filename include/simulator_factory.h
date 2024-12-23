#pragma once
#include <memory>
#include <vector>
#include <string>
#include "simulator.h"

struct SimulatorFactory {
    template <typename PType, typename VType, typename VFType, size_t N = 0, size_t K = 0>
    static std::unique_ptr<FluidSimulatorBase> create(const std::vector<std::string>& field_data_input) {
        return std::make_unique<FluidSimulator<PType, VType, VFType, N, K>>(field_data_input);
    }
};
