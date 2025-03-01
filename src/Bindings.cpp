#include "derived.hpp"
#include "timer.hpp"
#include "utils.hpp"

#ifdef TIME_ME

template <template <class> class SimType, typename T>
py::tuple runsim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    MyTimer timer;

    SimType<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    timer.add_lap();
    sim.simulation();
    timer.add_lap();

    py::tuple to_rtn = sim.export_python();
    timer.add_lap();
    timer.print_times();

    return to_rtn;
}

#else //-------------------------------------------------

template <template <class> class SimType, typename T>
py::tuple runsim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    SimType<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    sim.simulation();
    return sim.export_python();
}

#endif

//=================================

PYBIND11_MODULE(tasep, m) {
    m.def("D_kins_time", &runsim<fastTasep::CountKins, double>,
          "Returns only the total number of kins");
    m.def("F_kins_time", &runsim<fastTasep::CountKins, float>,
          "Returns only the total number of kins");
    m.def("D_neighbors", &runsim<fastTasep::Neighbors, double>,
          "Returns only the total number of kins");
    m.def("F_neighbors", &runsim<fastTasep::Neighbors, float>,
          "Returns only the total number of kins");
    m.def("D_nneighbors", &runsim<fastTasep::NearestNeighbor, double>,
          "Returns only the total number of kins");
    m.def("F_nneighbors", &runsim<fastTasep::NearestNeighbor, float>,
          "Returns only the total number of kins");
    m.def("D_new_itersim", &runsim<fastTasep::BasicIteration, double>,
          "Returns only the total number of kins");
    m.def("F_new_itersim", &runsim<fastTasep::BasicIteration, float>,
          "Returns only the total number of kins");
    m.def("D_profile", &runsim<fastTasep::Profile, double>,
          "Returns only the total number of kins");
    m.def("F_profile", &runsim<fastTasep::Profile, float>,
          "Returns only the total number of kins");
}
