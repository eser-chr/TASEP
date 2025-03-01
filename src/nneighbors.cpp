
/**
 * @file nneighbors.cpp
 * @brief Implementation of the NearestNeighbor class for the fastTasep namespace.
 */

#include "derived.hpp"

/**
 * @class NearestNeighbor
 * @brief A class that extends AbstractIteration to implement nearest neighbor interactions.
 * 
 * @tparam T The type of the rate constants (e.g., float, double).
 */

/**
 * @brief Constructor for the NearestNeighbor class.
 * 
 * @param L Length of the grid.
 * @param ITERS Number of iterations.
 * @param kon Binding rate constant.
 * @param koff Unbinding rate constant.
 * @param kstep Step rate constant.
 * @param q Probability factor.
 * @param kq Rate constant for q.
 */
template <typename T>
fastTasep::NearestNeighbor<T>::NearestNeighbor(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq), NEIGHBORS(L) {}

/**
 * @brief Binds a particle to the grid and updates the nearest neighbor counts.
 * 
 * @param side The position on the grid where the particle is to be bound.
 */
template <typename T>
void fastTasep::NearestNeighbor<T>::bind(int side) {
    // Implementation details...
}

/**
 * @brief Appends the current state to the trajectory.
 * 
 * This function is currently a placeholder and does not perform any operations.
 */
template <typename T>
void fastTasep::NearestNeighbor<T>::append_trajectory() {}

/**
 * @brief Exports the nearest neighbor data to a Python tuple.
 * 
 * @return A Python tuple containing the nearest neighbor data as a NumPy array.
 */
template <typename T>
py::tuple fastTasep::NearestNeighbor<T>::export_python() {
    // Implementation details...
}

// Explicit template instantiation for double and float types.
template class fastTasep::NearestNeighbor<double>;
template class fastTasep::NearestNeighbor<float>;