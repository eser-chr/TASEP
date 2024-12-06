#include "cooperative_tasep.h"
#include "count_kins.h"
#include "neighbors.h"
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

// -----------------------------------------
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
  return py::make_tuple(vector_to_numpy(std::move(sim.DATA), ITERS, L),
                        vector_to_numpy(std::move(sim.TIMES)),
                        vector_to_numpy(std::move(sim.ACTION)),
                        vector_to_numpy(std::move(sim.SIDE)));
}
#endif

// -----------------------------------------------------
// -----------------------------------------------------

#ifdef TIME_ME
py::tuple Dkins_time_sim(int L, int ITERS, double kon, double koff,
                         double kstep, double q, double kq,
                         bool verbose = false) {
  double t1, t2, t3;
  MyTimer timer1{};
  tasep::CountKins<double> sim(L, ITERS, kon, koff, kstep, q, kq);
  t1 = timer1.get();

  // if (verbose) {
  //   sim.printme();
  // }
  MyTimer timer2{};
  sim.simulation();
  t2 = timer2.get();

  MyTimer timer3{};
  const auto &to_rtn = py::make_tuple(vector_to_numpy(std::move(sim.KINS)),
                                      vector_to_numpy(std::move(sim.TIMES)));
  t3 = timer3.get();

  std::cout << "times: " << t1 << "," << t2 << "," << t3 << std::endl;
  return to_rtn;
}

py::tuple Fkins_time_sim(int L, int ITERS, float kon, float koff, float kstep,
                         float q, float kq, bool verbose = false) {
  double t1, t2, t3;
  MyTimer timer1{};
  tasep::CountKins<float> sim(L, ITERS, kon, koff, kstep, q, kq);
  t1 = timer1.get();

  // if (verbose) {
  //   sim.printme();
  // }
  MyTimer timer2{};
  sim.simulation();
  t2 = timer2.get();

  MyTimer timer3{};
  const auto &to_rtn = py::make_tuple(vector_to_numpy(std::move(sim.KINS)),
                                      vector_to_numpy(std::move(sim.TIMES)));
  t3 = timer3.get();

  std::cout << "times: " << t1 << "," << t2 << "," << t3 << std::endl;
  return to_rtn;
}

#else
py::tuple Dkins_time_sim(int L, int ITERS, double kon, double koff,
                         double kstep, double q, double kq,
                         bool verbose = false) {
  tasep::CountKins<double> sim(L, ITERS, kon, koff, kstep, q, kq);
  if (verbose) sim.printme();
  sim.simulation();
  return py::make_tuple(vector_to_numpy(std::move(sim.KINS)),
                        vector_to_numpy(std::move(sim.TIMES)));
}

py::tuple Fkins_time_sim(int L, int ITERS, float kon, float koff, float kstep,
                         double q, float kq) {
  tasep::CountKins<float> sim(L, ITERS, kon, koff, kstep, q, kq);
  sim.simulation();
  return py::make_tuple(vector_to_numpy(std::move(sim.KINS)),
                        vector_to_numpy(std::move(sim.TIMES)));
}
#endif

// -----------------------------------------------------
// NEIGHBORS

#ifdef TIME_ME
py::tuple Dneighbors_sim(int L, int ITERS, double kon, double koff,
                         double kstep, double q, double kq,
                         bool verbose = false) {
  double t1, t2, t3;
  MyTimer timer1{};
  tasep::Neighbors<double> sim(L, ITERS, kon, koff, kstep, q, kq);
  t1 = timer1.get();

  // if (verbose) {
  //   sim.printme();
  // }
  MyTimer timer2{};
  sim.simulation();
  t2 = timer2.get();

  MyTimer timer3{};
  const auto &to_rtn = py::make_tuple(vector_to_numpy(std::move(sim.NEIGHBORS)),
                                      vector_to_numpy(std::move(sim.TIMES)),
                                      vector_to_numpy(std::move(sim.ACTION)),
                                      vector_to_numpy(std::move(sim.SIDE)));
  t3 = timer3.get();

  std::cout << "times: " << t1 << "," << t2 << "," << t3 << std::endl;
  return to_rtn;
}

py::tuple Fneighbors_sim(int L, int ITERS, float kon, float koff, float kstep,
                         float q, float kq, bool verbose = false) {
  double t1, t2, t3;
  MyTimer timer1{};
  tasep::Neighbors<float> sim(L, ITERS, kon, koff, kstep, q, kq);
  t1 = timer1.get();

  // if (verbose) {
  //   sim.printme();
  // }
  MyTimer timer2{};
  sim.simulation();
  t2 = timer2.get();

  MyTimer timer3{};
  const auto &to_rtn = py::make_tuple(vector_to_numpy(std::move(sim.NEIGHBORS)),
                                      vector_to_numpy(std::move(sim.TIMES)),
                                      vector_to_numpy(std::move(sim.ACTION)),
                                      vector_to_numpy(std::move(sim.SIDE)));
  t3 = timer3.get();

  std::cout << "times: " << t1 << "," << t2 << "," << t3 << std::endl;
  return to_rtn;
}

#else
py::tuple Dneighbors_sim(int L, int ITERS, double kon, double koff,
                         double kstep, double q, double kq,
                         bool verbose = false) {
  tasep::Neighbors<double> sim(L, ITERS, kon, koff, kstep, q, kq);
  if (verbose) sim.printme();
  sim.simulation();
  return py::make_tuple(vector_to_numpy(std::move(sim.NEIGHBORS)),
                        vector_to_numpy(std::move(sim.TIMES)),
                        vector_to_numpy(std::move(sim.ACTION)),
                        vector_to_numpy(std::move(sim.SIDE)));
}

py::tuple Fneighbors_sim(int L, int ITERS, float kon, float koff, float kstep,
                         double q, float kq) {
  tasep::Neighbors<float> sim(L, ITERS, kon, koff, kstep, q, kq);
  sim.simulation();
  return py::make_tuple(vector_to_numpy(std::move(sim.NEIGHBORS)),
                        vector_to_numpy(std::move(sim.TIMES)),
                        vector_to_numpy(std::move(sim.ACTION)),
                        vector_to_numpy(std::move(sim.SIDE)));
}
#endif

//===========================================================================
//===========================================================================
//===========================================================================
#ifdef TIME_ME
py::tuple Dnneighbors_sim(int L, int ITERS, double kon, double koff,
                          double kstep, double q, double kq,
                          bool verbose = false) {
  double t1, t2, t3;
  MyTimer timer1{};
  tasep::NNeighbors<double> sim(L, ITERS, kon, koff, kstep, q, kq);
  t1 = timer1.get();

  // if (verbose) {
  //   sim.printme();
  // }
  MyTimer timer2{};
  sim.simulation();
  t2 = timer2.get();

  MyTimer timer3{};
  const auto &to_rtn = py::make_tuple(vector_to_numpy(std::move(sim.NEIGHBORS)),
                                      vector_to_numpy(std::move(sim.TIMES)),
                                      vector_to_numpy(std::move(sim.ACTION)),
                                      vector_to_numpy(std::move(sim.SIDE)));
  t3 = timer3.get();

  std::cout << "times: " << t1 << "," << t2 << "," << t3 << std::endl;
  return to_rtn;
}

py::tuple Fnneighbors_sim(int L, int ITERS, float kon, float koff, float kstep,
                          float q, float kq, bool verbose = false) {
  double t1, t2, t3;
  MyTimer timer1{};
  tasep::NNeighbors<float> sim(L, ITERS, kon, koff, kstep, q, kq);
  t1 = timer1.get();

  // if (verbose) {
  //   sim.printme();
  // }
  MyTimer timer2{};
  sim.simulation();
  t2 = timer2.get();

  MyTimer timer3{};
  const auto &to_rtn = py::make_tuple(vector_to_numpy(std::move(sim.NEIGHBORS)),
                                      vector_to_numpy(std::move(sim.TIMES)),
                                      vector_to_numpy(std::move(sim.ACTION)),
                                      vector_to_numpy(std::move(sim.SIDE)));
  t3 = timer3.get();

  std::cout << "times: " << t1 << "," << t2 << "," << t3 << std::endl;
  return to_rtn;
}

#else
py::tuple Dnneighbors_sim(int L, int ITERS, double kon, double koff,
                          double kstep, double q, double kq,
                          bool verbose = false) {
  tasep::NNeighbors<double> sim(L, ITERS, kon, koff, kstep, q, kq);
  if (verbose) sim.printme();
  sim.simulation();
  return py::make_tuple(vector_to_numpy(std::move(sim.NEIGHBORS)),
                        vector_to_numpy(std::move(sim.TIMES)),
                        vector_to_numpy(std::move(sim.ACTION)),
                        vector_to_numpy(std::move(sim.SIDE)));
}

py::tuple Fnneighbors_sim(int L, int ITERS, float kon, float koff, float kstep,
                          double q, float kq) {
  tasep::NNeighbors<float> sim(L, ITERS, kon, koff, kstep, q, kq);
  sim.simulation();
  return py::make_tuple(vector_to_numpy(std::move(sim.NEIGHBORS)),
                        vector_to_numpy(std::move(sim.TIMES)),
                        vector_to_numpy(std::move(sim.ACTION)),
                        vector_to_numpy(std::move(sim.SIDE)));
}
#endif

//--------------=================================
//--------------=================================
//--------------=================================

PYBIND11_MODULE(tasep, m) {
  m.def("itersim", &iter_sim, "A function to run the simulation");
  m.def("kins_time", &Dkins_time_sim, "Returns only the total number of kins");
  m.def("Fkins_time", &Fkins_time_sim, "Returns only the total number of kins");
  m.def("neighbors", &Dneighbors_sim, "Returns only the total number of kins");
  m.def("Fneighbors", &Fneighbors_sim, "Returns only the total number of kins");
  m.def("nneighbors", &Dnneighbors_sim,
        "Returns only the total number of kins");
  m.def("Fnneighbors", &Fnneighbors_sim,
        "Returns only the total number of kins");
}
