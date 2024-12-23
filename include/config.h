#pragma once
#include <memory>
#include <string>
#include <stdexcept>
#include <tuple>
#include <vector>
#include "fixed.h"

using SupportedTypes = std::tuple<
    float,
    double,
    Fixed<32,16>,
    Fixed<64,32>,
    FastFixed<16,8>,
    FastFixed<32,16>
>;

template<typename T>
struct is_valid_simulator_type : std::false_type {};

template<>
struct is_valid_simulator_type<float> : std::true_type {};

template<>
struct is_valid_simulator_type<double> : std::true_type {};

template<size_t N, size_t K>
struct is_valid_simulator_type<Fixed<N,K>> : std::true_type {};

template<size_t N, size_t K>
struct is_valid_simulator_type<FastFixed<N,K>> : std::true_type {};

struct TypeInfo {
    std::string base_type;
    size_t N{0};
    size_t K{0};
};

TypeInfo parse_type_info(const std::string& type_str) {
    if (type_str == "FLOAT") return {"float", 0, 0};
    if (type_str == "DOUBLE") return {"double", 0, 0};
    
    size_t pos = type_str.find('(');
    if (pos != std::string::npos) {
        std::string base = type_str.substr(0, pos);
        size_t comma_pos = type_str.find(',', pos);
        size_t close_pos = type_str.find(')', comma_pos);
        
        if (comma_pos != std::string::npos && close_pos != std::string::npos) {
            size_t N = std::stoul(type_str.substr(pos + 1, comma_pos - pos - 1));
            size_t K = std::stoul(type_str.substr(comma_pos + 1, close_pos - comma_pos - 1));
            return {base, N, K};
        }
    }
    throw std::runtime_error("Invalid type format: " + type_str);
}

template<typename T> 
bool matches_type_info(const TypeInfo& info) {
    if constexpr (std::is_same_v<T, float>) {
        return info.base_type == "FLOAT" || info.base_type == "float";
    }
    else if constexpr (std::is_same_v<T, double>) {
        return info.base_type == "DOUBLE" || info.base_type == "double"; 
    }
    else if constexpr (is_valid_simulator_type<T>::value) {
        if constexpr (std::is_base_of_v<Fixed<T::Bits, T::Fraction>, T>) {
            return (info.base_type == "FIXED" || info.base_type == "Fixed") && 
                   info.N == T::Bits && 
                   info.K == T::Fraction;
        }
        else if constexpr (std::is_base_of_v<FastFixed<T::Bits, T::Fraction>, T>) {
            return (info.base_type == "FAST_FIXED" || info.base_type == "FastFixed") &&
                   info.N == T::Bits && 
                   info.K == T::Fraction;
        }
    }
    return false;
}

template<typename Tuple, size_t I = 0>
typename std::tuple_element_t<0, Tuple> find_matching_type(const TypeInfo& info) {
    if constexpr (I >= std::tuple_size_v<Tuple>) {
        throw std::runtime_error("No matching type found for: " + info.base_type);
    }
    else {
        using T = std::tuple_element_t<I, Tuple>;
        if (matches_type_info<T>(info)) {
            return T{};
        }
        return find_matching_type<Tuple, I + 1>(info);
    }
}

inline std::unique_ptr<FluidSimulatorBase> createSimulatorInstance(
    const std::vector<std::string>& field_data_input,
    const char* p_type_str,
    const char* v_type_str,
    const char* vf_type_str
) {
    try {
        auto p_info = parse_type_info(p_type_str);
        auto v_info = parse_type_info(v_type_str);
        auto vf_info = parse_type_info(vf_type_str);

        using PType = std::remove_cvref_t<decltype(find_matching_type<SupportedTypes>(p_info))>;
        using VType = std::remove_cvref_t<decltype(find_matching_type<SupportedTypes>(v_info))>;
        using VFType = std::remove_cvref_t<decltype(find_matching_type<SupportedTypes>(vf_info))>;

        static_assert(is_valid_simulator_type<PType>::value, "Invalid pressure type");
        static_assert(is_valid_simulator_type<VType>::value, "Invalid velocity type");
        static_assert(is_valid_simulator_type<VFType>::value, "Invalid velocity field type");

        return std::make_unique<FluidSimulator<PType, VType, VFType>>(field_data_input);
    }
    catch (const std::exception& e) {
        throw std::runtime_error(std::string("Failed to create simulator: ") + e.what());
    }
}