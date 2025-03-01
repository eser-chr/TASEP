
/**
 * @file profile.cpp
 * @brief Implementation of the Profile class for the fastTasep namespace.
 * 
 * This file contains the implementation of the Profile class template, which
 * inherits from the AbstractIteration class. The Profile class is used to
 * manage and process the trajectory data of a TASEP (Totally Asymmetric Simple
 * Exclusion Process) simulation.
 */

#include "derived.hpp"

/**
 * @brief Constructor for the Profile class.
 * 
 * @tparam T The data type for the profile (e.g., float, double).
 * @param L The length of the profile.
 * @param ITERS The number of iterations.
 * @param kon The rate constant for binding.
 * @param koff The rate constant for unbinding.
 * @param kstep The step rate constant.
 * @param q The hopping rate.
 * @param kq The rate constant for hopping.
 * 
 * Initializes the Profile object with the given parameters and allocates
 * memory for the DATA vector.
 */
template <typename T>
fastTasep::Profile<T>::Profile(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq), DATA(L) {}

/**
 * @brief Appends the current trajectory to the profile data.
 * 
 * This function updates the DATA vector with the current state of the grid
 * if the current time (_time) is greater than or equal to the equilibration
 * time (Tequil). The grid values are accumulated over time.
 */
template <typename T>
void fastTasep::Profile<T>::append_trajectory() {
    if (this->_time >= this->Tequil) {
        for (size_t i = this->l_ghost; i < this->L + this->r_ghost; ++i) {
            DATA[i] += static_cast<double>(this->grid[i]) * this->dt;
        }
        ttotal += this->dt;
    }
}

/**
 * @brief Exports the profile data to a Python tuple.
 * 
 * @tparam T The data type for the profile (e.g., float, double).
 * @return A Python tuple containing the profile data as a NumPy array.
 * 
 * This function normalizes the DATA vector by dividing each element by the
 * total time (ttotal) and then converts the vector to a NumPy array, which
 * is returned as a Python tuple.
 */
template <typename T>
py::tuple fastTasep::Profile<T>::export_python() {
    for (auto &d : DATA) d /= ttotal;
    return py::make_tuple(vector_to_numpy(std::move(DATA)));
}

// Explicit instantiation
template class fastTasep::Profile<double>;
template class fastTasep::Profile<float>;