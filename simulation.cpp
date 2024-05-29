#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tuple>

#include "main.cpp"


// PYBIND11_MODULE(myfe, m)
// {
//   cout << "Loading myfe library" << endl;

//   ngcomp::ExportFESpace<ngcomp::MyFESpace>(m, "MyFESpace", true)
//     ;
// }    


// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>

// #include "simulation.cpp" // Include the C++ source file containing the function definition

namespace py = pybind11;

PYBIND11_MODULE(simulation, m) {
    m.doc() = "Simulation module";
    m.def("new_simulation", &simulation, "A function which performs simulation");
}
