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

std::vector<double> D_profile(int L, int ITERS, double kon, double koff, double kstep, double q,
                              double kq) {
    fastTasep::Profile<double> sim(L, ITERS, kon, koff, kstep, q, kq);
    sim.simulation();
    std::vector<double> data(L);
    double temp = static_cast<double>(L * (sim.t_final - sim.t_equil));
    for (size_t i = 0; i < L; i++) data[i] = static_cast<double>(sim.DATA[i]) / temp;
    return data;
}

auto many_sims_parallel(int L, int ITERS, double kon, double koff, double kstep, double q,
                        double kq, int SAMPLES = 10) {
    std::cout << "Running " << SAMPLES << " simulations in parallel" << std::endl;
    std::vector<std::future<std::vector<double>>> futures;
    for (int i = 0; i < SAMPLES; ++i) {
        futures.push_back(
            std::async(std::launch::async, D_profile, L, ITERS, kon, koff, kstep, q, kq));
    }

    std::vector<double> results(L, 0.0);
    for (auto &fut : futures) {
        auto partial_result = fut.get();
        for (size_t i = 0; i < L; i++) results[i] += partial_result[i] / SAMPLES;
    }

    return vector_to_numpy(std::move(results));
}

PYBIND11_MODULE(tasep_profile, m) {
    m.def("many_sims_parallel", &many_sims_parallel, "Run many simulations in parallel");
}
