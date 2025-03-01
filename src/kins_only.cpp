/**
 * @class CountKins
 * @brief A class that extends AbstractIteration to count kinetic events in a TASEP simulation.
 *
 * @tparam T The data type for the kinetic rates and other parameters (e.g., float, double).
 */

#include "derived.hpp"

/**
 * @brief Constructor for CountKins.
 *
 * @param L The length of the lattice.
 * @param ITERS The number of iterations.
 * @param kon The rate of binding.
 * @param koff The rate of unbinding.
 * @param kstep The rate of stepping.
 * @param q The hopping rate.
 * @param kq The rate of kinetic events.
 */
template <typename T>
fastTasep::CountKins<T>::CountKins(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);

/**
 * @brief Appends the current trajectory data to the KINS and TIMES vectors.
 */
template <typename T>
void fastTasep::CountKins<T>::append_trajectory();

/**
 * @brief Handles the binding event and increments the total kinetic events counter.
 *
 * @param side The side on which the binding occurs.
 */
template <typename T>
void fastTasep::CountKins<T>::bind(int side);

/**
 * @brief Handles the unbinding event and decrements the total kinetic events counter.
 *
 * @param side The side on which the unbinding occurs.
 */
template <typename T>
void fastTasep::CountKins<T>::unbind(int side);

/**
 * @brief Handles the stepping event and decrements the total kinetic events counter if stepping off
 * the lattice.
 *
 * @param side The side on which the stepping occurs.
 */
template <typename T>
void fastTasep::CountKins<T>::step(int side);

/**
 * @brief Exports the kinetic events and times to Python as numpy arrays.
 *
 * @return A Python tuple containing the KINS and TIMES vectors as numpy arrays.
 */
template <typename T>
py::tuple fastTasep::CountKins<T>::export_python();

/**
 * @brief Explicit template instantiation for CountKins with double type.
 */
template class fastTasep::CountKins<double>;

/**
 * @brief Explicit template instantiation for CountKins with float type.
 */
template class fastTasep::CountKins<float>;
