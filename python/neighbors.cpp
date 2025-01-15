#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <functional>
#include <future>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "derived.hpp"

namespace py = pybind11;
double calc_Tequil(double kon, double koff) { return 1.0 / (kon + koff); }

double calc_Tsim(double kon, double koff) { return 10.0 / (kon + koff); }
template <typename T>
py::array_t<T> vector_to_numpy(std::vector<T> &&vec) {
    py::array_t<T> result({vec.size()}, {sizeof(T)});
    auto result_buffer = result.request();
    T *result_ptr = static_cast<T *>(result_buffer.ptr);
    std::move(vec.begin(), vec.end(), result_ptr);
    return result;
}

template <typename T>
py::array_t<T> vector_to_numpy(std::vector<T> &&vec, size_t rows, size_t cols) {
    py::array_t<T> result({rows, cols});
    auto result_buffer = result.request();
    T *result_ptr = static_cast<T *>(result_buffer.ptr);
    std::move(vec.begin(), vec.end(), result_ptr);

    return result;
}

// auto D_neighbors(int L, int ITERS, double kon, double koff, double kstep, double q, double kq) {
//     fastTasep::Neighbors<double> sim(L, ITERS, kon, koff, kstep, q, kq);
//     sim.simulation();
//     return std::make_tuple(std::move(sim.NEIGHBORS), std::move(sim.TIMES),
//     std::move(sim.ACTION));
// }

std::vector<uint16_t> single_sim(int L, int ITERS, double kon, double koff, double kstep, double q,
                                 double kq) {
    double T_equil = calc_Tequil(kon, koff);
    double T_sim = calc_Tsim(kon, koff);
    fastTasep::Neighbors<double> sim(L, ITERS, kon, koff, kstep, q, kq);
    sim.simulation();

    // auto [nn, times, action] = D_neighbors(L, ITERS, kon, koff, kstep, q, kq);

    if (sim.TIMES.back() < T_equil) {
        std::cerr << sim.TIMES.back() << " < " << T_equil << std::endl;
        throw std::runtime_error("Not enough ITERS to capture the equilibrium");
    }
    if (sim.TIMES.back() < T_sim) {
        std::cerr << sim.TIMES.back() << " < " << T_sim << std::endl;
        std::cerr << "Warning: Not enough simulation time. Increase ITERS" << std::endl;
    }

    auto idx = std::lower_bound(sim.TIMES.begin(), sim.TIMES.end(), T_equil) - sim.TIMES.begin();
    int EBidx = std::accumulate(sim.ACTION.begin(), sim.ACTION.begin() + idx, 0,
                                [](int acc, int val) { return acc + (val == 0 ? 1 : 0); });

    return std::vector<uint16_t>(std::make_move_iterator(sim.NEIGHBORS.begin() + 2 * EBidx),
                                 std::make_move_iterator(sim.NEIGHBORS.end()));
}

// void many_sims_parallel(int L, int ITERS, double kon, double koff, double kstep, double q,
//                         double kq, int SAMPLES = 10) {
//     std::vector<std::future<std::vector<uint16_t>>> futures;
//     for (int i = 0; i < SAMPLES; ++i) {
//         futures.push_back(
//             std::async(std::launch::async, single_sim, L, ITERS, kon, koff, kstep, q, kq));
//     }

//     std::vector<uint16_t> results;
//     for (auto &fut : futures) {
//         const std::vector<uint16_t> &&partial_result = fut.get();
//         results.insert(results.end(), std::make_move_iterator(partial_result.begin()),
//                        std::make_move_iterator(partial_result.end()));
//     }
//     // return vector_to_numpy(std::move(results));
// }

class Many_Neighbor_Sims {
   public:
    Many_Neighbor_Sims(int L, int ITERS, double kon, double koff, double kstep, double q, double kq,
                       int SAMPLES = 10)
        : L(L), ITERS(ITERS), kon(kon), koff(koff), kstep(kstep), q(q), kq(kq), SAMPLES(SAMPLES) {
        results.reserve(ITERS * 2 * SAMPLES);
    }
    ~Many_Neighbor_Sims() { std::cout << "Destructor called" << std::endl; }
    void run() {
        for (int i = 0; i < SAMPLES; ++i) {
            futures.push_back(
                std::async(std::launch::async, single_sim, L, ITERS, kon, koff, kstep, q, kq));
        }

        // std::vector<uint16_t> results;
        for (auto &fut : futures) {
            const std::vector<uint16_t> &&partial_result = fut.get();
            results.insert(results.end(), std::make_move_iterator(partial_result.begin()),
                           std::make_move_iterator(partial_result.end()));
        }
    }

   public:
    std::vector<uint16_t> results;

   private:
    std::vector<std::future<std::vector<uint16_t>>> futures;
    int L, ITERS, SAMPLES;
    double kon, koff, kstep, q, kq;
};

py::array_t<uint16_t> py_many_sims_parallel(int L, int ITERS, double kon, double koff, double kstep,
                                            double q, double kq, int SAMPLES = 10) {
    std::vector<uint16_t> final(2 * ITERS * SAMPLES);
    {
        Many_Neighbor_Sims sim(L, ITERS, kon, koff, kstep, q, kq, SAMPLES);
        sim.run();
        for (size_t i = 0; i < sim.results.size(); i++) {
            final[i] = sim.results[i];
        }
    }
    return vector_to_numpy(std::move(final));
}

PYBIND11_MODULE(neighbors, m) {
    m.def("many_sims_parallel", &py_many_sims_parallel, "Run many simulations in parallel");
}
