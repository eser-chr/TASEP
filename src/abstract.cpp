#include "abstract.hpp"

template <typename T>
/**
 * @brief Constructs an AbstractIteration object with the given parameters.
 * 
 * @tparam T The type of the rate constants and propensities.
 * @param L The length of the grid.
 * @param ITERS The number of iterations.
 * @param kon The binding rate constant.
 * @param koff The unbinding rate constant.
 * @param kstep The stepping rate constant.
 * @param q The hopping rate constant.
 * @param kq The queue rate constant.
 * 
 * This constructor initializes the grid and propensities, and sets up the 
 * random number generator and the bucket manager. The grid is resized to 
 * accommodate ghost cells, and the propensities are initialized with the 
 * binding rate constant for the appropriate positions.
 */



/**
 * @brief Runs the entire simulation for the specified number of iterations.
 * 
 * This function repeatedly performs iterations until the specified number of
 * iterations (ITERS) is reached.
 */
fastTasep::AbstractIteration<T>::AbstractIteration(int L, int ITERS, T kon, T koff, T kstep, T q,
                                                   T kq)
    : L(L),
      ITERS(ITERS),
      kon(kon),
      koff(koff),
      kstep(kstep),
      q(q),
      kq(kq),
      _rng(std::make_unique<PCGRNG<T>>()) {
    grid.resize(L + ghost, 0);

    int total_size = 4 * (L + ghost);
    COLS = std::ceil(std::sqrt(total_size) / std::sqrt(2));
    ROWS = total_size / COLS + 2;
    Nactions_per_col = static_cast<double>(N_actions) / COLS;

    propensities.resize(ROWS * COLS, 0.0);

    for (int i = l_ghost * N_actions + BIND; i < (l_ghost + L) * N_actions; i += N_actions) {
        propensities[i] = kon;
    }

    _manager = std::make_unique<bucket<T>>(ROWS, COLS, propensities);
}



/**
 * @brief Prints the rate constants and hopping rate constant.
 */
template <typename T>
void fastTasep::AbstractIteration<T>::printme() {
    std::cout << "kon: " << kon << "\n";
    std::cout << "koff: " << koff << "\n";
    std::cout << "kstep: " << kstep << "\n";
    std::cout << "kq: " << kq << "\n";
    std::cout << "q: " << q << std::endl;
}

/**
 * @brief Binds a particle to the grid at the specified side.
 * 
 * @param side The position on the grid where the particle will bind.
 * 
 * This function updates the propensities and grid state to reflect the binding
 * of a particle at the specified side.
 */
template <typename T>
void fastTasep::AbstractIteration<T>::bind(int side) {
    ASSERT(grid[side] == 0, "wrong binding");
    int temp = side * N_actions;
    propensities[temp + STEP - N_actions] = 0.0;
    propensities[temp + BIND] = 0.0;
    propensities[temp + UNBIND] = koff;
    propensities[temp + STEP] = (T)kstep * (1 - grid[side + 1]);
    propensities[temp + DEACTIVATE] = 0.0;
    grid[side] = static_cast<uint8_t>(1);
}

/**
 * @brief Unbinds a particle from the grid at the specified side.
 * 
 * @param side The position on the grid where the particle will unbind.
 * 
 * This function updates the propensities and grid state to reflect the unbinding
 * of a particle at the specified side.
 */
template <typename T>
void fastTasep::AbstractIteration<T>::unbind(int side) {
    ASSERT(grid[side] == 1, "Wrong unbinding");
    int temp = side * N_actions;
    propensities[temp + STEP - N_actions] = (T)kstep * grid[side - 1];
    propensities[temp + BIND] = kon * q;
    propensities[temp + UNBIND] = 0.0;
    propensities[temp + STEP] = 0.0;
    propensities[temp + DEACTIVATE] = kq;
    grid[side] = static_cast<uint8_t>(0);
}

/**
 * @brief Moves a particle one step forward on the grid.
 * 
 * @param side The position on the grid where the particle will step.
 * 
 * This function updates the propensities and grid state to reflect the stepping
 * of a particle from the specified side to the next position.
 */
template <typename T>
void fastTasep::AbstractIteration<T>::step(int side) {
    int temp = side * N_actions;

    propensities[temp + STEP - N_actions] = (T)kstep * grid[side - 1];

    propensities[temp + BIND] = kon * q;
    propensities[temp + UNBIND] = 0;
    propensities[temp + STEP] = 0.0;
    propensities[temp + DEACTIVATE] = kq;

    propensities[temp + BIND + N_actions] = 0.0;
    propensities[temp + UNBIND + N_actions] = koff;
    propensities[temp + STEP + N_actions] = (T)kstep * (1 - grid[side + 2]);
    propensities[temp + DEACTIVATE + N_actions] = 0.0;
    grid[side] = static_cast<uint8_t>(0);
    grid[side + 1] = static_cast<uint8_t>(1);
}

/**
 * @brief Deactivates a particle at the specified side.
 * 
 * @param side The position on the grid where the particle will be deactivated.
 * 
 * This function updates the propensities and grid state to reflect the deactivation
 * of a particle at the specified side.
 */
template <typename T>
void fastTasep::AbstractIteration<T>::deactivate(int side) {
    ASSERT(grid[side] == 0, "Wrong step");
    int temp = side * N_actions;
    propensities[temp + BIND] = kon;
    propensities[temp + UNBIND] = 0;
    propensities[temp + STEP] = 0.0;
    propensities[temp + DEACTIVATE] = 0;
}

/**
 * @brief Fixes the boundary conditions of the grid.
 * 
 * This function ensures that the propensities and grid state at the boundaries
 * are correctly set to zero.
 */
template <typename T>
void fastTasep::AbstractIteration<T>::fixBoundaries() {
    std::fill(propensities.begin(), propensities.begin() + l_ghost * N_actions, 0.0);
    std::fill(propensities.begin() + (L + l_ghost) * N_actions,
              propensities.begin() + (L + ghost) * N_actions, 0.0);
    grid[0] = static_cast<uint8_t>(0);
    grid[L + 1] = static_cast<uint8_t>(0);
    ASSERT(grid[L + 2] == 0, "Boundaries were affected");
}

/**
 * @brief Fixes the cumulative sum of propensities after an update.
 * 
 * @param side The position on the grid where the update occurred.
 * 
 * This function updates the cumulative sum of propensities to reflect changes
 * in the grid state.
 */
template <typename T>
void fastTasep::AbstractIteration<T>::fixCumsum(int side) {
    int I = static_cast<int>((side - 1) * Nactions_per_col);
    int J = static_cast<int>((side + 2) * Nactions_per_col);
    _manager->update_sum_at_row(I);
    if (I != J) _manager->update_sum_at_row(J);
    _manager->refresh_cumsum();
}


/**
 * @brief Performs a single iteration of the simulation.
 * 
 * This function executes one iteration of the simulation, updating the grid state,
 * propensities, and cumulative sum based on the selected action.
 */
template <typename T>
void fastTasep::AbstractIteration<T>::iteration() {
    T back = _manager->_p_cum_sums.back();
    T r = _rng->random();
    dt = -(1.0 / back) * log(r);
    _time += dt;

    r = _rng->random() * back;
    _index = _manager->find_upper_bound(r);

    _action = _index & 3;
    _side = _index >> 2;

    ASSERT(_side != 0, "_side is 0");
    ASSERT(_side < L + 1, "side is bigger than the microtubule");

    executeAction(_action, _side);
    fixBoundaries();
    fixCumsum(_side);
}

template <typename T>
void fastTasep::AbstractIteration<T>::simulation() {
    while (_iter < ITERS) {
        append_trajectory();
        iteration();
        _iter++;
    }
}

template class fastTasep::AbstractIteration<float>;
template class fastTasep::AbstractIteration<double>;