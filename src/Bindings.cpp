#include "derived.hpp"
#include "timer.hpp"
#include "utils.hpp"

// If TIME_ME is defined, include timing functionality in the simulation.
#ifdef TIME_ME

/**
 * @brief Runs the simulation with timing enabled.
 *
 * This function initializes the simulation, runs it, and measures the time taken
 * for different parts of the simulation. It returns the results as a Python tuple.
 *
 * @tparam SimType Template parameter for the simulation type.
 * @tparam T Data type for the simulation parameters.
 * @param L Length of the simulation.
 * @param ITERS Number of iterations for the simulation.
 * @param kon Rate constant for binding.
 * @param koff Rate constant for unbinding.
 * @param kstep Step size for the rate constant.
 * @param q Parameter q for the simulation.
 * @param kq Parameter kq for the simulation.
 * @return py::tuple Results of the simulation.
 */
template <template <class> class SimType, typename T>
py::tuple runsim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    // Implementation with timing
}

#else  //-------------------------------------------------

/**
 * @brief Runs the simulation without timing.
 *
 * This function initializes the simulation, runs it, and returns the results
 * as a Python tuple.
 *
 * @tparam SimType Template parameter for the simulation type.
 * @tparam T Data type for the simulation parameters.
 * @param L Length of the simulation.
 * @param ITERS Number of iterations for the simulation.
 * @param kon Rate constant for binding.
 * @param koff Rate constant for unbinding.
 * @param kstep Step size for the rate constant.
 * @param q Parameter q for the simulation.
 * @param kq Parameter kq for the simulation.
 * @return py::tuple Results of the simulation.
 */
template <template <class> class SimType, typename T>
py::tuple runsim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    // Implementation without timing
}

#endif

//=================================

/**
 * @brief Binds the simulation functions to the Python module.
 *
 * This function defines the Python module 'tasep' and binds various simulation
 * functions to it. Each function runs a specific type of simulation and returns
 * the results as a Python tuple.
 *
 * @param m The Python module to which the functions are bound.
 */
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
    m.def("D_new_itersim", &runsim<fastTasep::BasicIteration, double>,
          "Returns only the total number of kins");
    m.def("F_new_itersim", &runsim<fastTasep::BasicIteration, float>,
          "Returns only the total number of kins");
    m.def("D_profile", &runsim<fastTasep::Profile, double>,
          "Returns only the total number of kins");
    m.def("F_profile", &runsim<fastTasep::Profile, float>, "Returns only the total number of kins");
}
