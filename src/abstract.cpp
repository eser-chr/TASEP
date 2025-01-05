#include "abstract.hpp"

template <typename T>
fastTasep::AbstractIteration<T>::AbstractIteration(int L, int ITERS, T kon, T koff,
                                                              T kstep, T q, T kq)
    : L(L),
      ITERS(ITERS),
      kon(kon),
      koff(koff),
      kstep(kstep),
      q(q),
      kq(kq),
      gen(std::random_device{}()),
      dis(0.0, 1.0) {

    int total_size = 4 * (L + ghost);
    COLS = std::ceil(std::sqrt(total_size) / std::sqrt(2));
    ROWS = total_size / COLS + 2;

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
    temp = side * N_actions;
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
    temp = side * N_actions;
    propensities[temp + STEP - N_actions] = (T)kstep * grid[side - 1];
    propensities[temp + BIND] = kon * q;
    propensities[temp + UNBIND] = 0.0;
    propensities[temp + STEP] = 0.0;
    propensities[temp + DEACTIVATE] = kq;
    grid[side] = 0;
}
template <typename T>
void fastTasep::AbstractIteration<T>::step(int side) {
    temp = side * N_actions;

    propensities[temp + STEP - N_actions] = (T)kstep * grid[side - 1];

    propensities[temp + BIND] = kon * q;
    propensities[temp + UNBIND] = 0;
    propensities[temp + STEP] = 0.0;
    propensities[temp + DEACTIVATE] = kq;

    propensities[temp + BIND + N_actions] = 0.0;
    propensities[temp + UNBIND + N_actions] = koff;
    propensities[temp + STEP + N_actions] = (T)kstep * (1 - grid[side + 2]);
    propensities[temp + DEACTIVATE + N_actions] = 0.0;
    grid[side] = 0;
    grid[side + 1] = 1;
}
template <typename T>
void fastTasep::AbstractIteration<T>::deactivate(int side) {
    ASSERT(grid[side] == 0, "Wrong step");
    temp = side * N_actions;
    propensities[temp + BIND] = kon;
    propensities[temp + UNBIND] = 0;
    propensities[temp + STEP] = 0.0;
    propensities[temp + DEACTIVATE] = 0;
}
template <typename T>
void fastTasep::AbstractIteration<T>::fix_boundaries() {
    std::fill(propensities.begin(), propensities.begin() + l_ghost * N_actions, 0.0);
    std::fill(propensities.begin() + (L + l_ghost) * N_actions,
              propensities.begin() + (L + ghost) * N_actions, 0.0);
    grid[0] = 0;
    grid[L + 1] = 0;
    ASSERT(grid[L + 2] == 0, "Boundaries were affected");
};
template <typename T>
void fastTasep::AbstractIteration<T>::fix_cumsum(int side) {
    int I = (side - 1) * N_actions / COLS;
    int J = (side + 2) * N_actions / COLS;
    _manager->update_sum_at_row(I);
    if (I != J) _manager->update_sum_at_row(J);
    _manager->refresh_cumsum();
    // _manager->update_cumsum();
};

template <typename T>
void fastTasep::AbstractIteration<T>::iteration() {
    r1 = dis(gen);
    r2 = dis(gen) * _manager->_p_cum_sums.back();
    dt = (1.0 / _manager->_p_cum_sums.back()) * log(1 / r1);
    _index = _manager->find_upper_bound(r2);


    _action = _index & 3;
    _side = _index >> 2;

    ASSERT(_side != 0, "_side is 0");
    ASSERT(_side < L + 1, "side is bigger than the microtubule");

    switch (_action) {
        case 0:
            bind(_side);
            break;
        case 1:
            unbind(_side);
            break;
        case 2:
            step(_side);
            break;
        case 3:
            deactivate(_side);
            break;
        default:
            throw std::runtime_error("An error occurred in switch statement");
    }

    fix_boundaries();
    fix_cumsum(_side);
};

template <typename T>
void fastTasep::AbstractIteration<T>::simulation() {
    while (_iter < ITERS) {
        append_trajectory();
        iteration();
        time += dt;
        _iter++;
    }
}


template class fastTasep::AbstractIteration<float>;
template class fastTasep::AbstractIteration<double>;