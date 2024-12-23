#pragma once

#include <vector>
#include <array>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <tuple>
#include "fixed.h"
#include "vector_field.h"

class FluidSimulatorBase {
public:
    virtual ~FluidSimulatorBase() = default;
    virtual void run(size_t steps, size_t checkpoint_interval) = 0;
    virtual void load_state(const char* filename) = 0;
    virtual void save_state(const char* filename) = 0;
};

template<typename PType, typename VType, typename VFType, size_t N = 0, size_t K = 0>
class FluidSimulator : public FluidSimulatorBase {
public:
    explicit FluidSimulator(const std::vector<std::string>& field_data_input);
    void run(size_t steps, size_t checkpoint_interval) override;
    void load_state(const char* filename) override;
    void save_state(const char* filename) override;

private:
    struct ParticleParams {
        char type;
        PType cur_p;
        std::array<VFType, 4> v;

        void swap_with(FluidSimulator* sim, size_t x, size_t y) {
            if (sim->use_static) {
                std::swap(sim->field[x][y], type);
                std::swap(sim->static_p[x][y], cur_p);
                v = sim->static_velocity.get_array(x, y);
                sim->static_velocity.set_array(x, y, v);
            } else {
                std::swap(sim->field_data[x][y], type);
                std::swap(sim->p[x][y], cur_p);
                v = sim->velocity.get_array(x, y);
                sim->velocity.set_array(x, y, v);
            }
        }
    };

    std::vector<std::string> field_data;
    std::vector<std::vector<PType>> p;
    std::vector<std::vector<PType>> old_p;

    VectorField<VFType> velocity;
    VectorField<VFType> velocity_flow;
    std::vector<std::vector<int>> last_use;

    std::mt19937 rnd;

    size_t rows{0}, cols{0};
    size_t UT{0};
    std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    std::vector<PType> rho;
    PType g{0};

    std::array<std::array<char, K + 1>, N> field{};
    std::array<std::array<PType, K>, N> static_p{};
    std::array<std::array<PType, K>, N> static_old_p{};
    StaticVectorField<VFType, N, K> static_velocity;
    StaticVectorField<VFType, N, K> static_velocity_flow;
    std::array<std::array<int, K>, N> static_last_use{};
    bool use_static{N > 0 && K > 0};

    void initialize_field();
    PType random01() noexcept;

    std::tuple<PType, bool, std::pair<int, int>>
    propagate_flow(int x, int y, PType lim);
    void propagate_stop(int x, int y, bool force = false);
    PType move_prob(int x, int y);
    bool propagate_move(int x, int y, bool is_first, int depth = 0) {
        const int MAX_DEPTH = 1000;
        last_use[x][y] = UT - is_first;
        if (depth > MAX_DEPTH) {
            std::cerr << "Max recursion depth reached at (" << x << ", " << y << ")\n";
            return false;
        }

        if (use_static) {
            return static_propagate_move(x, y, is_first);
        }

        bool ret = false;
        int target_x = -1, target_y = -1;
        do {
            std::array<VFType, 4> thresholds{};
            std::array<VFType, 4> velocities{};
            VFType sum = VFType(0);

            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (nx >= static_cast<int>(rows) || ny >= static_cast<int>(cols) || nx < 0 || ny < 0) {
                    continue;
                }
                if ((field_data[nx][ny] == '#' || last_use[nx][ny] == UT)) {
                    velocities[i] = VFType(0);
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy, deltas);
                if (v < 0) {
                    thresholds[i] = sum;
                    continue;
                }
                sum += v;
                thresholds[i] = sum;
            }

            if (sum == VFType(0)) {
                break;
            }

            VFType r = random01() * sum;
            size_t dir = 0;
            for (size_t i = 0; i < thresholds.size(); i++) {
                if (thresholds[i] > r) {
                    dir = i;
                    break;
                }
            }
            
            auto [dx, dy] = deltas[dir];

            target_x = x + dx;
            target_y = y + dy;
            
            if (target_x < 0 || target_y < 0 || target_x >= static_cast<int>(rows) || target_y >= static_cast<int>(cols)) {
                continue;
            }

            ret = (last_use[target_x][target_y] == UT - 1 || propagate_move(target_x, target_y, false, depth + 1));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (nx > -1 && ny > -1 && nx < rows && ny < cols && (field_data[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy, deltas) < 0)) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(this, x, y);
                pp.swap_with(this, target_x, target_y);
                pp.swap_with(this, x, y);
            }
        }
        return ret;
    }

    std::tuple<PType, bool, std::pair<size_t, size_t>>
    static_propagate_flow(size_t x, size_t y, PType lim);
    void static_propagate_stop(size_t x, size_t y, bool force = false);
    PType static_move_prob(size_t x, size_t y);
    bool static_propagate_move(size_t x, size_t y, bool is_first);
};

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
FluidSimulator<PType, VType, VFType, N, K>::FluidSimulator(
    const std::vector<std::string>& field_data_input)
    : field_data(field_data_input), rnd(1337) {
    initialize_field();
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
void FluidSimulator<PType, VType, VFType, N, K>::initialize_field() {
    std::cout << "Field data contains " << field_data.size() << " lines:\n";
    for (size_t i = 0; i < field_data.size(); i++) {
        std::cout << "Line " << i << ": " << field_data[i] << "\n";
    }
    std::cout << "---End of field data---\n";
    std::stringstream ss(field_data[0]);
    ss >> rows >> cols;
    
    g = std::stod(field_data[1]);

    if (rows > N && N != 0) {
        throw std::runtime_error("Invalid rows");
    }
    if (cols > K && K != 0) {
        throw std::runtime_error("Invalid cols");
    }
    if (rows == 0 || cols == 0) {
        throw std::runtime_error("Invalid rows or cols");
    }

    rho.resize(256, PType(0.01));

    std::vector<std::string> field_lines(field_data.begin() + 2, field_data.begin() + 2 + rows);

    for (size_t i = 2 + rows; i < field_data.size(); i++) {
        std::string line = field_data[i];
        if (line.empty()) continue;
        
        std::stringstream rho_ss(line);
        char symbol;
        std::string equals;
        double value;
        
        if (rho_ss >> symbol >> equals >> value) {
            rho[static_cast<size_t>(symbol)] = PType(value);
        }
    }

    if (use_static) {
        for (size_t x = 0; x < rows; ++x) {
            for (size_t y = 0; y < cols; ++y) {
                field[x][y] = field_lines[x][y];
            }
        }
    } else {
        field_data = field_lines; 
        p.resize(rows, std::vector<PType>(cols, PType(0)));
        old_p.resize(rows, std::vector<PType>(cols, PType(0)));
        velocity.init(rows, cols);
        velocity_flow.init(rows, cols);
        last_use.resize(rows, std::vector<int>(cols, 0));
    }

    std::cout << "\n=== Current Simulator State ===\n";
    std::cout << "Dimensions: " << rows << "x" << cols << "\n";
    std::cout << "Gravity: " << g << "\n";
    
    std::cout << "\nField Layout:\n";
    if (use_static) {
        for (size_t x = 0; x < rows; ++x) {
            for (size_t y = 0; y < cols; ++y) {
                std::cout << field[x][y];
            }
            std::cout << "\n";
        }
    } else {
        for (const auto& row : field_data) {
            std::cout << row << "\n";
        }
    }

    std::cout << "\nDensity Values:\n";
    for (size_t i = 0; i < rho.size(); ++i) {
        if (rho[i] != PType(0.01)) { 
            std::cout << "'" << (char)i << "': " << rho[i] << "\n";
        }
    }

    std::cout << "\nCurrent Pressures:\n";
    for (size_t x = 0; x < rows; ++x) {
        for (size_t y = 0; y < cols; ++y) {
            if (use_static) {
                std::cout << static_p[x][y] << " ";
            } else {
                std::cout << p[x][y] << " ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "===========================\n\n";
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
PType FluidSimulator<PType, VType, VFType, N, K>::random01() noexcept {
    static std::uniform_real_distribution<VFType> dist(0.0, 1.0);
    return dist(rnd);
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
std::tuple<PType, bool, std::pair<int, int>>
FluidSimulator<PType, VType, VFType, N, K>::propagate_flow(int x, int y, PType lim) {
    last_use[x][y] = UT - 1;
    
    if (use_static) {
        return static_propagate_flow(x, y, lim);
    }

    if (x < 0 || y < 0 || x >= rows || y >= cols || field_data[x][y] == '#') {
        return {PType(0), false, {0, 0}};
    }

    PType ret = PType(0);

    for (const auto& [dx, dy] : deltas) {
        int nx = x + dx;
        int ny = y + dy;
        
        if (nx < 0 || ny < 0 || nx >= rows || ny >= cols || field_data[nx][ny] == '#') {
            continue;
        }

        if (field_data[nx][ny] != '#' && last_use[nx][ny] < UT) {
            auto cap = velocity.get(x, y, dx, dy, deltas);
            auto flow = velocity_flow.get(x, y, dx, dy, deltas);
            if (flow == cap) {
                continue;
            }
            // assert(v >= velocity_flow.get(x, y, dx, dy));
            auto vp = std::min(lim, cap - flow);
            if (last_use[nx][ny] == UT - 1) {
                velocity_flow.add(x, y, dx, dy, vp, deltas);
                last_use[x][y] = UT;
                // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                return {vp, 1, {nx, ny}};
            }
            auto [t, prop, end] = propagate_flow(nx, ny, vp);
            ret += t;
            if (prop) {
                velocity_flow.add(x, y, dx, dy, t, deltas);
                last_use[x][y] = UT;
                // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                return {t, prop && end != std::pair(x, y), end};
            }
        }
    }
    last_use[x][y] = UT;

    return {ret, 0, {0, 0}};
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
void FluidSimulator<PType, VType, VFType, N, K>::propagate_stop(int x, int y, bool force) {
    if (use_static) {
        static_propagate_stop(x, y, force);
        return;
    }

    if (!force) {
        bool stop = true;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (nx >= 0 && ny >= 0 && nx < static_cast<int>(rows) && ny < static_cast<int>(cols) &&
                field_data[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && 
                velocity.get(x, y, dx, dy, deltas) > VFType(0)) {
                stop = false;
                break;
            }
        }
        if (!stop) {
            return;
        }
    }

    last_use[x][y] = UT;
    for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (nx < 0 || ny < 0 || nx >= rows || ny >= cols) {
            continue;
        }
        if (field_data[nx][ny] == '#' || last_use[nx][ny] == UT || 
            velocity.get(x, y, dx, dy, deltas) > VFType(0)) {
            continue; 
        }
        propagate_stop(nx, ny);
    }
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
PType FluidSimulator<PType, VType, VFType, N, K>::move_prob(int x, int y) {
    if (use_static) {
        return static_move_prob(x, y);
    }

    PType sum = PType(0);
    for (const auto& [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field_data[nx][ny] == '#' || last_use[nx][ny] == UT) {
            continue;
        }
        VFType v = velocity.get(x, y, dx, dy, deltas);
        if (v >= VFType(0)) {
            sum += v;
        }
    }
    return sum;
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
void FluidSimulator<PType, VType, VFType, N, K>::run(size_t steps, size_t checkpoint_interval) {
    std::vector<std::vector<int>> dirs(rows, std::vector<int>(cols, 0));

    if (!use_static) {
        for (size_t x = 0; x < rows; ++x) {
            for (size_t y = 0; y < cols; ++y) {
                if (field_data[x][y] == '#') continue;
                for (auto [dx, dy] : deltas) {
                    if (x + dx >= 0 && x + dx < rows && y + dy >= 0 && y + dy < cols) {
                        dirs[x][y] += (field_data[x + dx][y + dy] != '#');
                    }
                }
            }
        }
    } else {
        for (size_t x = 0; x < rows; ++x) {
            for (size_t y = 0; y < cols; ++y) {
                if (field[x][y] == '#') continue;
                for (auto [dx, dy] : deltas) {
                    if (x + dx >= 0 && x + dx < rows && y + dy >= 0 && y + dy < cols) {
                        dirs[x][y] += (field[x + dx][y + dy] != '#');
                    }
                }
            }
        }
    }

    for (size_t step = 0; step < steps; ++step) {
        std::cout << "Starting step " << step + 1 << "\n";

        PType total_delta_p = PType(0);

        if (!use_static) {
            std::cout << "Applying gravity...\n";

            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < cols; ++y) {
                    if (field_data[x][y] == '#') continue;
                    if (field_data[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, g, deltas);
                }
            }

            old_p = p;

            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < cols; ++y) {
                    if (field_data[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (nx < 0 || ny < 0 || nx >= rows || ny >= cols) {continue;}
                        if (field_data[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            auto delta_p = old_p[x][y] - old_p[nx][ny];
                            auto force = delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy, deltas);
                            if (contr * rho[(int) field_data[nx][ny]] >= force) {
                                contr -= force / rho[(int) field_data[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field_data[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, force / rho[(int) field_data[x][y]], deltas);
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
                        }
                    }
                }
            }

            velocity_flow = {};
            velocity_flow.init(rows, cols);

            bool prop = false;
            do {
                UT += 2;
                prop = false;
                for (size_t x = 0; x < rows; ++x) {
                    for (size_t y = 0; y < cols; ++y) {
                        if (field_data[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1);
                            if (t > 0) {
                                prop = true;
                            }
                        }
                    }
                }
            } while (prop);

            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < cols; ++y) {
                    if (field_data[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (nx < 0 || ny < 0 || nx >= rows || ny >= cols) {continue;}
                        auto old_v = velocity.get(x, y, dx, dy, deltas);
                        auto new_v = velocity_flow.get(x, y, dx, dy, deltas);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy, deltas) = new_v;
                            auto force = (old_v - new_v) * rho[(int) field_data[x][y]];
                            if (field_data[x][y] == '.')
                                force *= 0.8;
                            if (field_data[nx][ny] == '#') {
                                p[x][y] += force / dirs[x][y];
                                total_delta_p += force / dirs[x][y];
                            } else {
                                p[nx][ny] += force / dirs[nx][ny];
                                total_delta_p += force / dirs[nx][ny];
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;

            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < cols; ++y) {
                    if (field_data[x][y] != '#' && last_use[x][y] != UT) {
                        auto pr = random01();
                        auto pr1 = move_prob(x, y);
                        if (pr < pr1) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop) {
                std::cout << "Tick " << step++ << ":\n";
                for (size_t x = 0; x < rows; ++x) {
                    std::cout << field_data[x] << "\n";
                }
            }

        } else {
            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < cols; ++y) {
                    if (field[x][y] == '#') continue;
                    if (x + 1 < rows && field[x + 1][y] != '#') {
                        static_velocity.add(x, y, 1, 0, g, deltas);
                    }
                }
            }

            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < cols; ++y) {
                    static_old_p[x][y] = static_p[x][y];
                }
            }

            static_velocity_flow = {};

            bool prop;
            do {
                UT += 2;
                prop = false;
                for (size_t x = 0; x < rows; ++x) {
                    for (size_t y = 0; y < cols; ++y) {
                        if (field[x][y] != '#' && static_last_use[x][y] != UT) {
                            if (random01() < static_move_prob(x, y)) {
                                prop = true;
                                static_propagate_move(x, y, true);
                            } else {
                                static_propagate_stop(x, y, true);
                            }
                        }
                    }
                }
            } while (prop);

            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < cols; ++y) {
                    if (field[x][y] == '#') continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = static_velocity.get(x, y, dx, dy, deltas);
                        auto new_v = static_velocity_flow.get(x, y, dx, dy, deltas);
                        if (old_v > 0) {
                            static_velocity.add(x, y, dx, dy, new_v - old_v, deltas);
                            auto force = (old_v - new_v) * rho[(int)field[x][y]];
                            if (field[x][y] == '.') force *= PType(0.8);

                            size_t nx = x + dx, ny = y + dy;
                            if (nx < rows && ny < cols && field[nx][ny] == '#') {
                                static_p[x][y] += force / PType(dirs[x][y]);
                                total_delta_p += force / PType(dirs[x][y]);
                            } else if (nx < rows && cols < cols) {
                                static_p[nx][ny] += force / PType(dirs[nx][ny]);
                                total_delta_p += force / PType(dirs[nx][ny]);
                            }
                        }
                    }
                }
            }
        }

    }
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
std::tuple<PType, bool, std::pair<size_t, size_t>>
FluidSimulator<PType, VType, VFType, N, K>::static_propagate_flow(size_t x, size_t y, PType lim) {
    if (x >= rows || y >= cols || field[x][y] == '#') {
        return {PType(0), false, {0, 0}};
    }

    PType sum_flow = PType(0);
    bool has_flow = false;
    std::pair<size_t, size_t> max_flow_pos{0, 0};
    PType max_flow = PType(0);

    for (const auto& [dx, dy] : deltas) {
        size_t nx = x + dx;
        size_t ny = y + dy;
        
        if (nx >= rows || ny >= cols || field[nx][ny] == '#') {
            continue;
        }

        PType dp = static_p[x][y] - static_p[nx][ny];
        if (dp <= PType(0)) {
            continue;
        }

        PType flow = std::min(dp, lim);
        sum_flow += flow;
        has_flow = true;

        if (flow > max_flow) {
            max_flow = flow;
            max_flow_pos = {nx, ny};
        }

        static_velocity_flow.add(x, y, dx, dy, flow, deltas);
    }

    return {sum_flow, has_flow, max_flow_pos};
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
void FluidSimulator<PType, VType, VFType, N, K>::static_propagate_stop(size_t x, size_t y, bool force) {
    if (x >= rows || y >= cols || field[x][y] == '#') {
        return;
    }

    bool should_stop = force;
    if (!force) {
        PType v_sum = PType(0);
        for (const auto& [dx, dy] : deltas) {
            v_sum += abs(static_velocity.get(x, y, dx, dy, deltas));
        }
        should_stop = v_sum < PType(0.1);
    }

    if (should_stop) {
        for (const auto& [dx, dy] : deltas) {
            static_velocity.add(x, y, dx, dy, PType(0), deltas);
        }
    }
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
PType FluidSimulator<PType, VType, VFType, N, K>::static_move_prob(size_t x, size_t y) {
    if (x >= rows || y >= cols || field[x][y] == '#') {
        return PType(0);
    }

    PType sum = PType(0);
    for (const auto& [dx, dy] : deltas) {
        VFType v = static_velocity.get(x, y, dx, dy, deltas);
        if (v > VFType(0)) {
            sum += v;
        }
    }
    return sum; 
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
bool FluidSimulator<PType, VType, VFType, N, K>::static_propagate_move(size_t x, size_t y, bool is_first) {
    if (x >= rows || y >= cols || field[x][y] == '#') {
        return false;
    }

    if (!is_first && static_last_use[x][y] == UT) {
        return false;
    }
    static_last_use[x][y] = UT - is_first;

    std::array<VFType, 4> thresholds;
    VFType sum = VFType(0);
    
    for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        size_t nx = x + dx, ny = y + dy;
        
        if (nx >= rows || ny >= cols || field[nx][ny] == '#' || static_last_use[nx][ny] == UT) {
            thresholds[i] = sum;
            continue;
        }
        
        VFType v = static_velocity.get(x, y, dx, dy, deltas);
        if (v <= VFType(0)) {
            thresholds[i] = sum;
            continue;
        }
        
        sum += v;
        thresholds[i] = sum;
    }

    if (sum == VFType(0)) {
        return false;
    }

    VFType r = random01() * sum;
    size_t dir = std::upper_bound(thresholds.begin(), thresholds.end(), r) - thresholds.begin();
    
    auto [dx, dy] = deltas[dir];
    size_t nx = x + dx, ny = y + dy;

    bool moved = (static_last_use[nx][ny] == UT - 1) || static_propagate_move(nx, ny, false);
    
    if (moved && !is_first) {
        std::swap(field[x][y], field[nx][ny]);
        std::swap(static_p[x][y], static_p[nx][ny]);
        
        auto tmp_vel = static_velocity.get_array(x, y);
        static_velocity.set_array(x, y, static_velocity.get_array(nx, ny));
        static_velocity.set_array(nx, ny, tmp_vel);
    }

    return moved;
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
void FluidSimulator<PType, VType, VFType, N, K>::save_state(const char* filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening file for saving state");
    }
    file << rows << " " << cols << std::endl;
    file << g << std::endl;
    
    if (!use_static) {
        for (size_t x = 0; x < rows; x++) {
            file << field_data[x] << std::endl;
        }
    } else {
        for (size_t x = 0; x < rows; x++) {
            for (size_t y = 0; y < cols; y++) {
                file << field[x][y];
            }
            file << std::endl;
        }
    }

    double default_rho = 0.01;
    for (size_t i = 0; i < rho.size(); ++i) {
        if (rho[i] != static_cast<PType>(default_rho)) {
            file << static_cast<char>(i) << " = " << rho[i] << std::endl;
        }
    }

    file.close();
}

template<typename PType, typename VType, typename VFType, size_t N, size_t K>
void FluidSimulator<PType, VType, VFType, N, K>::load_state(const char* filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Failed to open file for reading");
    }

    size_t new_rows, new_cols;
    file >> new_rows >> new_cols;
    file >> g;

    if (use_static) {
        if (new_rows > rows || new_cols > cols) {
            throw std::runtime_error("Saved dimensions exceed static allocation");
        }
        for (size_t i = 0; i < new_rows; i++) {
            std::string line;
            file >> line;
            for (size_t j = 0; j < new_cols; j++) {
                field[i][j] = line[j];
            }
        }

        for (size_t i = 0; i < new_rows; i++) {
            for (size_t j = 0; j < new_cols; j++) {
                file >> static_p[i][j];
                file >> static_old_p[i][j];
            }
        }

        static_velocity.init(new_rows, new_cols);
        static_velocity_flow.init(new_rows, new_cols);

        for (size_t i = 0; i < new_rows; i++) {
            for (size_t j = 0; j < new_cols; j++) {
                for (size_t k = 0; k < 4; k++) {
                    VFType val;
                    file >> val;
                    static_velocity.add(i, j, deltas[k].first, deltas[k].second, val, deltas);
                }
            }
        }

        file >> UT;
        rows = new_rows;
        cols = new_cols;
    } else {
        field_data.resize(new_rows);
        for (size_t i = 0; i < new_rows; i++) {
            file >> field_data[i];
        }

        p.resize(new_rows, std::vector<PType>(new_cols));
        old_p.resize(new_rows, std::vector<PType>(new_cols));
        
        for (size_t i = 0; i < new_rows; i++) {
            for (size_t j = 0; j < new_cols; j++) {
                file >> p[i][j];
                file >> old_p[i][j];
            }
        }

        velocity.init(new_rows, new_cols);
        velocity_flow.init(new_rows, new_cols);
        last_use.resize(new_rows, std::vector<int>(new_cols));
        
        for (size_t i = 0; i < new_rows; i++) {
            for (size_t j = 0; j < new_cols; j++) {
                for (size_t k = 0; k < 4; k++) {
                    VFType val;
                    file >> val;
                    velocity.add(i, j, deltas[k].first, deltas[k].second, val, deltas);
                }
            }
        }

        file >> UT;
        rows = new_rows;
        cols = new_cols;
    }

    double default_rho = 0.01;
    rho.resize(256, PType(default_rho));
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        char ch;
        std::string eq;
        double value;
        if (iss >> ch >> eq >> value) {
            rho[static_cast<size_t>(ch)] = PType(value);
        }
    }
    file.close();
}

