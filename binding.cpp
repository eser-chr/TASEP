#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <tuple>
#include <vector>

using namespace std;
namespace py = pybind11;

tuple<vector<vector<int>>, vector<double>, vector<int>, vector<double>> 
cyclic_simulation(int L, double T, double kon, double koff, double kstep, double kq, double q, int initial_density);


std::tuple< std::vector<std::vector<int>>, std::vector<double>, std::vector<int>, std::vector<double> > 
line_simulation(int L, double T, double kon, double koff, double kstep, double kq, double q);

template <typename T>
py::array_t<T> convert_to_array(std::vector<T> & vec){
    auto res = py::array_t<T>(vec.size());
    auto buffer = res.request();
    auto ptr = (T*) buffer.ptr;

    for(int i =0; i<vec.size(); i++){
        ptr[i] = vec[i];
    }
    return res;
}





py::tuple py_cyclic_sim(int L, double T, double kon, double koff, double kstep, double kq, double q, int initial_density) {
    auto result = cyclic_simulation(L, T, kon, koff, kstep, kq, q, initial_density);
    
    return py::make_tuple(
        get<0>(result),
        get<1>(result), 
        get<2>(result), 
        get<3>(result)
        );
}

py::tuple py_line_sim(int L, double T, double kon, double koff, double kstep, double kq, double q) {
    auto result = line_simulation(L, T, kon, koff, kstep, kq, q);
    
    return py::make_tuple(
        get<0>(result),
        get<1>(result), 
        get<2>(result), 
        get<3>(result)
        );
}






PYBIND11_MODULE(binding, m) {
    m.def("_line_sim", &py_line_sim, "A function to run the simulation");
    m.def("_cyclic_sim", &py_cyclic_sim, "A function to run the simulation");
}
