#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <tuple>
#include <vector>
#include "cooperative_tasep.h"

using namespace std;
namespace py = pybind11;

py::array_t<int> matrix_to_numpy(const std::vector<std::vector<int>>& vec) {
    // Get the dimensions
    size_t rows = vec.size();
    size_t cols = vec.empty() ? 0 : vec[0].size();

    // Allocate a buffer to hold the data
    py::array_t<int> result({rows, cols});
    auto result_buffer = result.request();
    int* result_ptr = static_cast<int*>(result_buffer.ptr);

    // Copy the data from the std::vector to the buffer
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result_ptr[i * cols + j] = vec[i][j];
        }
    }

    return result;
}

py::array_t<int> vector_to_numpy(const std::vector<int>& vec) {
    // Allocate a buffer to hold the data
    py::array_t<int> result(vec.size());
    auto result_buffer = result.request();
    int* result_ptr = static_cast<int*>(result_buffer.ptr);

    // Copy the data from the std::vector to the buffer
    for (size_t i = 0; i < vec.size(); ++i) {
        result_ptr[i] = vec[i];
    }

    return result;
}


py::array_t<double> vector_to_numpy(const std::vector<double>& vec) {
    // Allocate a buffer to hold the data
    py::array_t<double> result(vec.size());
    auto result_buffer = result.request();
    double* result_ptr = static_cast<double*>(result_buffer.ptr);

    // Copy the data from the std::vector to the buffer
    for (size_t i = 0; i < vec.size(); ++i) {
        result_ptr[i] = vec[i];
    }

    return result;
}

// std::tuple< std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>, std::vector<double> > 
// line_simulation(int L, double T, double kon, double koff, double kstep, double kq, double q);


py::tuple cooperative_tasep_sim(int L, double T, double kon, double koff, double kstep, double kq, double q){
    CooperativeTasep sim(L, T, kon, koff, kstep, q, kq);
    sim.simulation();
    auto result = sim.get_results();
    return py::make_tuple(
        matrix_to_numpy(get<0>(result)),
        vector_to_numpy(get<1>(result))
    );
}
py::tuple specific_cooperative_tasep_sim(int L, double T, double kon, double koff, double kstep, double kq, double q){
    SpecificCooperativeTasep sim(L, T, kon, koff, kstep, q, kq);
    sim.simulation();
    auto result = sim.get_results();
    return py::make_tuple(
        matrix_to_numpy(get<0>(result)),
        matrix_to_numpy(get<1>(result)),
        matrix_to_numpy(get<2>(result)),
        vector_to_numpy(get<3>(result)),
        vector_to_numpy(get<4>(result)),
        vector_to_numpy(get<5>(result))
    );
}


// //CYCLE--------------------
// py::tuple py_cyclic_sim(int L, double T, double kon, double koff, double kstep, double kq, double q, int initial_density) {
//     auto result = cyclic_simulation(L, T, kon, koff, kstep, kq, q, initial_density);
    
//     return py::make_tuple(
//         matrix_to_numpy(get<0>(result)),
//         matrix_to_numpy(get<1>(result)),
//         vector_to_numpy(get<2>(result)), 
//         vector_to_numpy(get<3>(result)),
//         vector_to_numpy(get<4>(result))
//         );
// }








PYBIND11_MODULE(cooperative_tasep_lib, m) {
    m.def("sim", &cooperative_tasep_sim, "A function to run the simulation");
    m.def("ssim", &specific_cooperative_tasep_sim, "A function to run the simulation");
    // m.def("_cyclic_sim", &py_cyclic_sim, "A function to run the simulation");
}
