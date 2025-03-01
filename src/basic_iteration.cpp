
/**
 * @file basic_iteration.cpp
 * @brief Implementation of the BasicIteration class for the TASEP simulation.
 */

#include <cstring>
#include "derived.hpp"

/**
 * @brief Constructor for the BasicIteration class.
 * 
 * @tparam T The data type for the iteration (e.g., float, double).
 * @param L The length of the grid.
 * @param ITERS The number of iterations.
 * @param kon The rate of particle attachment.
 * @param koff The rate of particle detachment.
 * @param kstep The step rate.
 * @param q The hopping rate.
 * @param kq The rate of hopping with a particle.
 */
template <typename T>
fastTasep::BasicIteration<T>::BasicIteration(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq) {
    DATA.resize(L * ITERS);
    TIMES.resize(ITERS);
}

/**
 * @brief Appends the current trajectory to the DATA vector.
 * 
 * This function copies the current state of the grid (excluding ghost cells)
 * into the DATA vector and records the current time in the TIMES vector.
 */
template <typename T>
void fastTasep::BasicIteration<T>::append_trajectory() {
    std::copy(this->grid.begin() + this->l_ghost, this->grid.end() - this->r_ghost,
              DATA.begin() + this->_iter * this->L);

    TIMES[this->_iter] = this->_time;
}

/**
 * @brief Exports the trajectory data and times to Python.
 * 
 * @return A Python tuple containing the trajectory data and times as NumPy arrays.
 */
template <typename T>
py::tuple fastTasep::BasicIteration<T>::export_python() {
    return py::make_tuple(vector_to_numpy(std::move(DATA), this->ITERS, this->L),
                          vector_to_numpy(std::move(TIMES)));
}

// Explicit template instantiation for float and double types.
template class fastTasep::BasicIteration<float>;
template class fastTasep::BasicIteration<double>;