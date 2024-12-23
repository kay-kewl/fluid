#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

template<typename T>
class VectorField {
public:
    using value_type = T;
    using size_type = std::size_t;
    using delta_array = std::array<std::pair<int, int>, 4>;

    VectorField() = default;
    VectorField(const VectorField&) = default;
    VectorField(VectorField&&) noexcept = default;
    VectorField& operator=(const VectorField&) = default;
    VectorField& operator=(VectorField&&) noexcept = default;
    ~VectorField() = default;

    void init(size_type rows, size_type cols) {
        v.resize(rows);
        for(auto& row : v) {
            row.resize(cols);
            for(auto& cell : row) {
                cell.fill(T(0));
            }
        }
    }

    T& add(size_type x, size_type y, int dx, int dy, T dv,
           const delta_array& deltas) {
        assert(is_valid_position(x, y));
        
        auto it = std::find(deltas.begin(), deltas.end(), std::make_pair(dx, dy));
        if (it == deltas.end()) {
            std::cout << "Invalid delta: (" << dx << "," << dy << ")\n";
            std::cout << "Available deltas:\n";
            for (const auto& [dx_, dy_] : deltas) {
                std::cout << "(" << dx_ << "," << dy_ << ") ";
            }
            std::cout << "\n";
            throw std::runtime_error("Invalid delta values");
        }
        
        size_t i = std::distance(deltas.begin(), it);
        assert(i < deltas.size());
        return v[x][y][i] += dv;
    }

    T& get(size_type x, size_type y, int dx, int dy, const delta_array& deltas) {
        size_t i = std::distance(deltas.begin(), std::find(deltas.begin(), deltas.end(), std::make_pair(dx, dy)));
        return v[x][y][i];
    }

    void reset() {
        for(auto& row : v) {
            for(auto& cell : row) {
                cell.fill(T(0));
            }
        }
    }

    const auto& operator[](size_type i) const { 
        assert(i < rows());
        return v[i]; 
    }
    
    auto& operator[](size_type i) { 
        assert(i < rows());
        return v[i]; 
    }

    static bool is_valid_delta(int dx, int dy, const delta_array& deltas) {
        return std::find(deltas.begin(), deltas.end(), 
                        std::make_pair(dx, dy)) != deltas.end();
    }

    const T& at(size_type x, size_type y, size_type i) const {
        assert(is_valid_position(x, y) && i < 4);
        return v[x][y][i];
    }

    T& at(size_type x, size_type y, size_type i) {
        assert(is_valid_position(x, y) && i < 4);
        return v[x][y][i];
    }

    void swap(VectorField& other) noexcept {
        v.swap(other.v);
    }

    size_type rows() const { return v.size(); }
    size_type cols() const { return v.empty() ? 0 : v[0].size(); }
    bool empty() const { return v.empty(); }
    bool is_valid_position(size_type x, size_type y) const {
        return x < rows() && y < cols();
    }

    std::array<T, 4> get_array(size_t x, size_t y) const {
        return v[x][y];
    }
    
    void set_array(size_t x, size_t y, const std::array<T, 4>& arr) {
        v[x][y] = arr;
    }

private:
    std::vector<std::vector<std::array<T, 4>>> v;
};

template<typename T, size_t N, size_t K>
class StaticVectorField {
public:
    using value_type = T;
    using size_type = std::size_t;
    using delta_array = std::array<std::pair<int, int>, 4>;
    
    StaticVectorField() { reset(); }
    StaticVectorField(const StaticVectorField&) = default;
    StaticVectorField(StaticVectorField&&) noexcept = default;
    StaticVectorField& operator=(const StaticVectorField&) = default;
    StaticVectorField& operator=(StaticVectorField&&) noexcept = default;
    ~StaticVectorField() = default;

    void init(size_t rows, size_t cols) {
        v = {};
    }

    T& add(size_type x, size_type y, int dx, int dy, T dv,
           const delta_array& deltas) {
        assert(is_valid_position(x, y));
        size_t i = std::distance(deltas.begin(),
            std::find(deltas.begin(), deltas.end(), std::make_pair(dx, dy)));
        assert(i < deltas.size());
        return v[x][y][i] += dv;
    }

    T& get(size_type x, size_type y, int dx, int dy,
           const delta_array& deltas) {
        assert(is_valid_position(x, y));
        size_t i = std::distance(deltas.begin(),
            std::find(deltas.begin(), deltas.end(), std::make_pair(dx, dy)));
        assert(i < deltas.size());
        return v[x][y][i];
    }

    void reset() {
        for(auto& plane : v) {
            for(auto& row : plane) {
                row.fill(T(0));
            }
        }
    }

    static bool is_valid_delta(int dx, int dy, const delta_array& deltas) {
        return std::find(deltas.begin(), deltas.end(), 
                        std::make_pair(dx, dy)) != deltas.end();
    }

    const T& at(size_type x, size_type y, size_type i) const {
        assert(is_valid_position(x, y) && i < 4);
        return v[x][y][i];
    }

    T& at(size_type x, size_type y, size_type i) {
        assert(is_valid_position(x, y) && i < 4);
        return v[x][y][i];
    }

    void swap(StaticVectorField& other) noexcept {
        std::swap(v, other.v);
    }

    const auto& operator[](size_type i) const { 
        assert(i < N);
        return v[i]; 
    }
    
    auto& operator[](size_type i) { 
        assert(i < N);
        return v[i]; 
    }

    std::array<T, 4> get_array(size_t x, size_t y) const {
        return v[x][y];
    }
    
    void set_array(size_t x, size_t y, const std::array<T, 4>& arr) {
        v[x][y] = arr;
    }

    static constexpr size_type rows() { return N; }
    static constexpr size_type cols() { return K; }
    static constexpr bool empty() { return N == 0 || K == 0; }
    static constexpr bool is_valid_position(size_type x, size_type y) {
        return x < N && y < K;
    }

private:
    std::array<std::array<std::array<T, 4>, K>, N> v;
};