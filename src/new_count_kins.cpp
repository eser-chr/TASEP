#include "count_kins.h"
#include "timer.hpp"

#include <cstdlib> // for abort()
#include <iostream>

#ifdef DEBUGB
#define ASSERT(condition, message)                                             \
  do {                                                                         \
    if (!(condition)) {                                                        \
      std::cout << "Assertion failed: " << (message) << std::endl;             \
    }                                                                          \
  } while (false)

// } //     if (!(condition)) {                                          \
    //         std::cerr << "Assertion failed: (" #condition ") in "    \
    //                   << __FILE__ << ", line " << __LINE__ << ": "   \
    //                   << message << std::endl;                       \
    //         std::abort();                                            \
    //     }                                                            \


#else
#define ASSERT(condition, message)
#endif

template <typename T>
fastTasep::CountKins<T>::CountKins(int L, int ITERS, T kon, T koff, T kstep,
                                   T q, T kq)
    : L(L), ITERS(ITERS), kon(kon), koff(koff), kstep(kstep), q(q), kq(kq),
      gen(std::random_device{}()), dis(0.0, 1.0) {

  // std::cout << " I am in the constructor" << std::endl;

  COLS = 40;
  ROWS = 2 + ((L + ghost) * N_actions / COLS);

  // assert((L + l_ghost) % 10 == 0);

  grid.resize(L + ghost, 0);
  propensities.resize(ROWS * COLS, 0.0);
  // sum_propensities.resize(N_actions * (L + ghost), 0.0);

  sum_of_rows.resize(ROWS, 0.0);
  cumsum_of_rows.resize(ROWS, 0.0);
  cumsum_of_rows[0] = 0.0;

  KINS.resize(ITERS);
  TIMES.resize(ITERS);

  for (int i = l_ghost * N_actions + BIND; i < (l_ghost + L) * N_actions;
       i += N_actions) {
    propensities[i] = kon;
  }

#ifdef DEBUGB
  for (int kl = 0; kl < 20; kl++) {
    std::cout << propensities[kl] << " ";
  }
  std::cout << std::endl;
  // std::partial_sum(propensities.begin(), propensities.end(),
  //                  sum_propensities.begin());
  std::cout << "BEfore cumm" << std::endl;
#endif

  for (int row = 0; row < ROWS; row++) {
    sum_of_rows[row] =
        std::accumulate(propensities.begin() + row * COLS,
                        propensities.begin() + (row + 1) * COLS, 0.0);
  }

  for (int row = 1; row < ROWS; row++) {
    cumsum_of_rows[row] = cumsum_of_rows[row - 1] + sum_of_rows[row - 1];
  }
#ifdef DEBUGB

  std::cout << _iter << " " << ITERS << std::endl;
  std::cout << "Done from constructor" << std::endl;
  for (int kl = 0; kl < ROWS; kl++) {
    std::cout << sum_of_rows[kl] << " ";
  }
  std::cout << std::endl;
  for (int kl = 0; kl < ROWS; kl++) {
    std::cout << cumsum_of_rows[kl] << " ";
  }
  std::cout << std::endl;
#endif
}

// --------------------------------------------------------------------------

template <typename T> void fastTasep::CountKins<T>::printme() {
  std::cout << "kon: " << kon << "\n";
  std::cout << "koff: " << koff << "\n";
  std::cout << "kstep: " << kstep << "\n";
  std::cout << "kq: " << kq << "\n";
  std::cout << "q: " << q << std::endl;
}
template <typename T> void fastTasep::CountKins<T>::bind(int side) {

  ASSERT(grid[side] == 0, "wrong binding");
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = 0.0;
  propensities[temp + BIND] = 0;
  propensities[temp + UNBIND] = koff;
  propensities[temp + STEP] = (double)kstep * (1 - grid[side + 1]);
  propensities[temp + DEACTIVATE] = 0.0;
  grid[side] = 1;
  TOTAL_KINS++;
}
template <typename T> void fastTasep::CountKins<T>::unbind(int side) {
  ASSERT(grid[side] == 1, "Wrong unbinding");
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = (double)kstep * grid[side - 1];
  propensities[temp + BIND] = kon * q;
  propensities[temp + UNBIND] = 0.0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = kq;
  grid[side] = 0;
  TOTAL_KINS--;
}
template <typename T> void fastTasep::CountKins<T>::step(int side) {
  temp = side * N_actions;

  propensities[temp + STEP - N_actions] = (double)kstep * grid[side - 1];

  propensities[temp + BIND] = kon * q;
  propensities[temp + UNBIND] = 0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = kq;

  propensities[temp + BIND + N_actions] = 0.0;
  propensities[temp + UNBIND + N_actions] = koff;
  propensities[temp + STEP + N_actions] = (double)kstep * (1 - grid[side + 2]);
  propensities[temp + DEACTIVATE + N_actions] = 0.0;
  grid[side] = 0;
  grid[side + 1] = 1;
  if (side == L) TOTAL_KINS--;
}
template <typename T> void fastTasep::CountKins<T>::deactivate(int side) {
  ASSERT(grid[side] == 0, "Wrong step");
  temp = side * N_actions;
  propensities[temp + BIND] = kon;
  propensities[temp + UNBIND] = 0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = 0;
}
template <typename T> void fastTasep::CountKins<T>::fix_boundaries() {
  std::fill(propensities.begin(), propensities.begin() + l_ghost * N_actions,
            0.0);
  std::fill(propensities.begin() + (L + l_ghost) * N_actions,
            propensities.begin() + (L + ghost) * N_actions, 0.0);
  grid[0] = 0;
  grid[L + 1] = 0;
  ASSERT(grid[L + 2] == 0, "Boundaries were affected");
};
template <typename T> void fastTasep::CountKins<T>::fix_cumsum(int side) {
  int first_el_index =
      (side - 1) *
      N_actions; // find the ROW of the first element that might have changed
  int I = first_el_index / COLS;

  sum_of_rows[I] = std::accumulate(propensities.begin() + I * COLS,
                                   propensities.begin() + (I + 1) * COLS, 0.0);
  sum_of_rows[I + 1] =
      std::accumulate(propensities.begin() + (I + 1) * COLS,
                      propensities.begin() + (I + 2) * COLS, 0.0);

  for (int temp = I; temp < ROWS; temp++) { // remove temp
    cumsum_of_rows[temp + 1] = cumsum_of_rows[temp] + sum_of_rows[temp];
  }
  // std::cout << std::endl;
};

template <typename T> void fastTasep::CountKins<T>::set_index() {
  _INDEX = std::distance(
      cumsum_of_rows.begin(),
      std::upper_bound(cumsum_of_rows.begin(), cumsum_of_rows.end(), r2));

  _index = (_INDEX - 1) * COLS;
  // std::cout << _index << "    <---INDEX BEFORE" << std::endl;
  T temp_sum = cumsum_of_rows[_INDEX - 1];
  // std::cout << "temp_sum, r2  " << temp_sum << "," << r2 << std::endl;
  while (temp_sum < r2) {
    // std::cout << "_index, temp_sum, r2 " << _index << ", " << temp_sum << ",
    // "
    //           << r2 << std::endl;
    temp_sum += propensities[_index];
    _index++;
  }
  _index--;
};

template <typename T> void fastTasep::CountKins<T>::iteration() {
  r1 = dis(gen);
  r2 = dis(gen) * cumsum_of_rows.back();

  dt = (1.0 / cumsum_of_rows.back()) * log(1 / r1);
  // _index = 0;

  set_index();

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
template <typename T> void fastTasep::CountKins<T>::append_trajectory() {

  // uint16_t sum = std::accumulate(grid.begin(), grid.end(), 0);
  KINS[_iter] = TOTAL_KINS;
  TIMES[_iter] = time;
};
template <typename T> void fastTasep::CountKins<T>::simulation() {
  while (_iter < ITERS) {
    append_trajectory();
    iteration();
    time += dt;
    _iter++;
  }
}

template class fastTasep::CountKins<float>;
template class fastTasep::CountKins<double>;