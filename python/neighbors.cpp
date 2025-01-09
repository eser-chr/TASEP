#include <algorithm>
#include <functional>
#include <future>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
// #include "utils.hpp"
#include "derived.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
double calc_Tequil(double kon, double koff) {
    return 1.0 / (kon + koff);
}

double calc_Tsim(double kon, double koff) {
    return 10.0 / (kon + koff);
}
template <typename T>
py::array_t<T> vector_to_numpy(const std::vector<T> &&vec) {
  py::array_t<T> result({vec.size()}, {sizeof(T)});
  auto result_buffer = result.request();
  T *result_ptr = static_cast<T *>(result_buffer.ptr);
  std::move(vec.begin(), vec.end(), result_ptr);
  return result;
}

template <typename T>
py::array_t<T> vector_to_numpy(const std::vector<T> &&vec, size_t rows,
                               size_t cols) {
  py::array_t<T> result({rows, cols});
  auto result_buffer = result.request();
  T *result_ptr = static_cast<T *>(result_buffer.ptr);
  std::move(vec.begin(), vec.end(), result_ptr);

  return result;
}

auto D_neighbors(int L, int ITERS, double kon, double koff, double kstep, double q, double kq) {
    fastTasep::Neighbors<double> sim(L, ITERS, kon, koff, kstep, q, kq);
    sim.simulation();
    return std::make_tuple(
        std::move(sim.NEIGHBORS),
        std::move(sim.TIMES),
        std::move(sim.ACTION),
        std::move(sim.SIDE)
    );
}

std::vector<uint16_t> single_sim(int L, int ITERS, double kon, double koff, double kstep,
                                 double q, double kq) {
    double T_equil = calc_Tequil(kon, koff);
    double T_sim = calc_Tsim(kon, koff);

    auto [nn, times, action, sides] = D_neighbors(L, ITERS, kon, koff, kstep, q, kq);

    if (times.back() < T_equil) {
        std::cerr << times.back() << " < " << T_equil << std::endl;
        throw std::runtime_error("Not enough ITERS to capture the equilibrium");
    }
    if (times.back() < T_sim) {
        std::cerr << times.back() << " < " << T_sim << std::endl;
        std::cerr << "Warning: Not enough simulation time. Increase ITERS" << std::endl;
    }

    auto idx = std::lower_bound(times.begin(), times.end(), T_equil) - times.begin();
    int EBidx = std::accumulate(
        action.begin(), action.begin() + idx, 0,
        [](int acc, int val) { return acc + (val == 0 ? 1 : 0); }
    );

    return std::vector<uint16_t>(nn.begin() + 2 * EBidx, nn.end());
}

auto many_sims_parallel(int L, int ITERS, double kon, double koff, double kstep,
                             double q, double kq, int SAMPLES = 10) {
    std::vector<std::future<std::vector<uint16_t>>> futures;
    for (int i = 0; i < SAMPLES; ++i) {
        futures.push_back(std::async(std::launch::async, single_sim,
                                     L, ITERS, kon, koff, kstep, q, kq));
    }

    std::vector<uint16_t> results;
    for (auto& fut : futures) {
        auto partial_result = fut.get();
        results.insert(results.end(), partial_result.begin(), partial_result.end());
    }

    return vector_to_numpy(std::move(results));
}

PYBIND11_MODULE(neighbors, m) {
    m.def("many_sims_parallel", &many_sims_parallel, "Run many simulations in parallel");
}
