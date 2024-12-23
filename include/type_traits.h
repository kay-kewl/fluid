#pragma once
#include <type_traits>

template<typename T>
struct size_of {
    static constexpr size_t value = sizeof(T);
};