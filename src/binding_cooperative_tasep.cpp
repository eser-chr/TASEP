#include "cooperative_tasep.h"
#include "timer.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace std;
namespace py = pybind11;

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

// py::array_t<double> vector_to_numpy(const std::vector<double>& vec) {
//     // Allocate a buffer to hold the data
//     py::array_t<double> result(vec.size());
//     auto result_buffer = result.request();
//     double* result_ptr = static_cast<double*>(result_buffer.ptr);

//     // Copy the data from the std::vector to the buffer
//     for (size_t i = 0; i < vec.size(); ++i) {
//         result_ptr[i] = vec[i];
//     }

//     return result;
// }

// std::tuple< std::vector<std::vector<int>>, std::vector<std::vector<int>>,
// std::vector<double>, std::vector<int>, std::vector<double> >
// line_simulation(int L, double T, double kon, double koff, double kstep,
// double kq, double q);

// py::tuple cooperative_tasep_sim(int L, double T, double kon, double koff,
// double kstep, double kq, double q){
//     CooperativeTasep sim(L, T, kon, koff, kstep, q, kq);
//     sim.simulation();
//     auto result = sim.get_results();
//     return py::make_tuple(
//         matrix_to_numpy(get<0>(result)),
//         vector_to_numpy(get<1>(result))
//     );
// }
// py::tuple specific_cooperative_tasep_sim(int L, double T, double kon, double
// koff, double kstep, double q, double kq){
//     SpecificCooperativeTasep sim(L, T, kon, koff, kstep, q, kq);
//     sim.simulation();
//     auto result = sim.get_results();
//     return py::make_tuple(
//         matrix_to_numpy(get<0>(result)),
//         matrix_to_numpy(get<1>(result)),
//         matrix_to_numpy(get<2>(result)),
//         vector_to_numpy(get<3>(result)),
//         vector_to_numpy(get<4>(result)),
//         vector_to_numpy(get<5>(result)),
//         vector_to_numpy(get<6>(result))
//     );
// }

// py::tuple basic_tasep_sim(int L, double T, double kon, double koff,
//                           double kstep, double q, double kq, bool trajectory,
//                           bool details, double period) {
//   tasep::Basic sim(L, T, kon, koff, kstep, q, kq, trajectory, details,
//   period); sim.simulation(); if (trajectory && details) {
//     auto traj = sim.get_trajectory();
//     auto details = sim.get_details();
//     auto &data = get<0>(traj);

//     // for(int i =0; i<100; i++){std::cout<<dts[i]<<" ";}
//     return py::make_tuple(
//         vector_to_numpy(data, data.size() / L, L),
//         vector_to_numpy(get<1>(traj)), vector_to_numpy(get<0>(details)),
//         vector_to_numpy(get<1>(details)), vector_to_numpy(get<2>(details)));
//   } else if (trajectory && !details) {
//     auto traj = sim.get_trajectory();
//     auto data = get<0>(traj);
//     return py::make_tuple(vector_to_numpy(data, data.size() / L, L),
//                           vector_to_numpy(get<1>(traj)));
//   } else if (!trajectory && details) {
//     auto details = sim.get_details();
//     return py::make_tuple(vector_to_numpy(get<0>(details)),
//                           vector_to_numpy(get<1>(details)),
//                           vector_to_numpy(get<2>(details)));

//   } else {
//     throw std::runtime_error("You need some analytics");
//   }
// }
#ifdef TIME_ME
py::tuple iter_sim(int L, int ITERS, double kon, double koff, double kstep,
                   double q, double kq, bool verbose = false) {
  double t1, t2, t3;
  MyTimer timer1{};
  tasep::BasicIteration sim(L, ITERS, kon, koff, kstep, q, kq);
  t1 = timer1.get();

  if (verbose) {
    sim.printme();
  }
  MyTimer timer2{};
  sim.simulation();
  t2 = timer2.get();

  MyTimer timer3{};
  const auto &to_rtn =
      py::make_tuple(vector_to_numpy(std::move(sim.DATA), ITERS, L),
                     vector_to_numpy(std::move(sim.TIMES)),
                     vector_to_numpy(std::move(sim.ACTION)),
                     vector_to_numpy(std::move(sim.SIDE)));
  t3 = timer3.get();

  std::cout << "times: " << t1 << "," << t2 << "," << t3 << std::endl;
  return to_rtn;
}
#else
py::tuple iter_sim(int L, int ITERS, double kon, double koff, double kstep,
                   double q, double kq, bool verbose = false) {
  tasep::BasicIteration sim(L, ITERS, kon, koff, kstep, q, kq);
  if (verbose) sim.printme();
  sim.simulation();
  return py::make_tuple(
    vector_to_numpy(std::move(sim.DATA), ITERS, L),
                        vector_to_numpy(std::move(sim.TIMES)),
                        vector_to_numpy(std::move(sim.ACTION)),
                        vector_to_numpy(std::move(sim.SIDE))
                        );
}
#endif

PYBIND11_MODULE(tasep, m) {
  // m.def("sim", &cooperative_tasep_sim, "A function to run the simulation");
  // m.def("ssim", &specific_cooperative_tasep_sim, "A function to run the
  // simulation");
  // m.def("bsim", &basic_tasep_sim, "A function to run the simulation");
  m.def("itersim", &iter_sim, "A function to run the simulation");
}
