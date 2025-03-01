#include "abstract.hpp"

template <typename T>
fastTasep::AbstractIteration<T>::AbstractIteration(int L, int ITERS, T kon, T koff, T kstep, T q,
                                                   T kq)
    : L(L),
      ITERS(ITERS),
      kon(kon),
      koff(koff),
      kstep(kstep),
      q(q),
      kq(kq),
      _rng(std::make_unique<PCGRNG<T>>())
{
    int total_size = 4 * (L + ghost);
    COLS = std::ceil(std::sqrt(total_size) / std::sqrt(2));
    ROWS = total_size / COLS + 2;
    Nactions_per_col = static_cast<double>(N_actions) / COLS;

    grid.resize(L + ghost, 0);
    propensities.resize(ROWS * COLS, 0.0);

    for (int i = l_ghost * N_actions + BIND; i < (l_ghost + L) * N_actions; i += N_actions) {
        propensities[i] = kon;
    }

    _manager = std::make_unique<bucket<T>>(ROWS, COLS, propensities);
}

// --------------------------------------------------------------------------

template <typename T>
void fastTasep::AbstractIteration<T>::printme() {
    std::cout << "kon: " << kon << "\n";
    std::cout << "koff: " << koff << "\n";
    std::cout << "kstep: " << kstep << "\n";
    std::cout << "kq: " << kq << "\n";
    std::cout << "q: " << q << std::endl;
}
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
template <typename T>
void fastTasep::AbstractIteration<T>::deactivate(int side) {
    ASSERT(grid[side] == 0, "Wrong step");
    int temp = side * N_actions;
    propensities[temp + BIND] = kon;
    propensities[temp + UNBIND] = 0;
    propensities[temp + STEP] = 0.0;
    propensities[temp + DEACTIVATE] = 0;
}
template <typename T>
void fastTasep::AbstractIteration<T>::fixBoundaries() {
    std::fill(propensities.begin(), propensities.begin() + l_ghost * N_actions, 0.0);
    std::fill(propensities.begin() + (L + l_ghost) * N_actions,
              propensities.begin() + (L + ghost) * N_actions, 0.0);
    grid[0] = static_cast<uint8_t>(0);
    grid[L + 1] = static_cast<uint8_t>(0);
    ASSERT(grid[L + 2] == 0, "Boundaries were affected");
};
template <typename T>
void fastTasep::AbstractIteration<T>::fixCumsum(int side) {
    int I = static_cast<int>((side - 1) * Nactions_per_col);
    int J = static_cast<int>((side + 2) * Nactions_per_col);
    _manager->update_sum_at_row(I);
    if (I != J) _manager->update_sum_at_row(J);
    _manager->refresh_cumsum();
};

template <typename T>
void fastTasep::AbstractIteration<T>::iteration() {
    T back = _manager->_p_cum_sums.back();
    T r = _rng->random();
    _time -= (1.0 / back) * log(r);

    r = _rng->random() * back;
    _index = _manager->find_upper_bound(r);

    _action = _index & 3;
    _side = _index >> 2;

    ASSERT(_side != 0, "_side is 0");
    ASSERT(_side < L + 1, "side is bigger than the microtubule");

    executeAction(_action, _side);
    fixBoundaries();
    fixCumsum(_side);
};

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