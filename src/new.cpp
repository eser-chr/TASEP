#include "new.h"
#include "timer.hpp"

#ifdef DEBUG
#define ASSERT(condition) assert(condition)
#else
#define ASSERT(condition)
#endif

template <typename T>
fastTasep::BasicIteration<T>::BasicIteration(int L, int ITERS, T kon, T koff,
                                             T kstep, T q, T kq)
    : L(L), ITERS(ITERS), kon(kon), koff(koff), kstep(kstep), q(q), kq(kq),
      gen(std::random_device{}()), dis(0.0, 1.0) {

L =    1006;
K = 64;     

  grid.resize(L + ghost, 0);
  propensities.resize(N_actions * (L + ghost), 0.0);
  sum_propensities.resize(N_actions * (L + ghost), 0.0);
  DATA.resize(L * ITERS);
  TIMES.resize(ITERS);
  ACTION.resize(ITERS);
  SIDE.resize(ITERS);

  for (int i = l_ghost * N_actions + BIND;
       i < propensities.size() - r_ghost * N_actions; i += N_actions) {
    propensities[i] = kon;
  }

#ifdef DEBUG
  for (int kl = 0; kl < 20; kl++) {
    std::cout << propensities[kl] << " ";
  }
  std::cout << std::endl;
#endif

  std::partial_sum(propensities.begin(), propensities.end(),
                   sum_propensities.begin());

#ifdef DEBUG
  for (int kl = 0; kl < 20; kl++) {
    std::cout << sum_propensities[kl] << " ";
  }
  std::cout << std::endl;
#endif
}

// --------------------------------------------------------------------------

template <typename T> void fastTasep::BasicIteration<T>::printme() {
  std::cout << "kon: " << kon << "\n";
  std::cout << "koff: " << koff << "\n";
  std::cout << "kstep: " << kstep << "\n";
  std::cout << "kq: " << kq << "\n";
  std::cout << "q: " << q << std::endl;
}
template <typename T> void fastTasep::BasicIteration<T>::bind(int side) {

  ASSERT(grid[side] == 0);
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = 0.0;
  propensities[temp + BIND] = 0;
  propensities[temp + UNBIND] = koff;
  propensities[temp + STEP] = (double)kstep * (1 - grid[side + 1]);
  propensities[temp + DEACTIVATE] = 0.0;
  grid[side] = 1;
}
template <typename T> void fastTasep::BasicIteration<T>::unbind(int side) {
  ASSERT(grid[side] == 1);
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = (double)kstep * grid[side - 1];
  propensities[temp + BIND] = kon * q;
  propensities[temp + UNBIND] = 0.0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = kq;
  grid[side] = 0;
}
template <typename T> void fastTasep::BasicIteration<T>::step(int side) {
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
}
template <typename T> void fastTasep::BasicIteration<T>::deactivate(int side) {
  ASSERT(grid[side] == 0);
  temp = side * N_actions;
  propensities[temp + BIND] = kon;
  propensities[temp + UNBIND] = 0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = 0;
}
template <typename T> void fastTasep::BasicIteration<T>::fix_boundaries() {
  std::fill(propensities.begin(), propensities.begin() + l_ghost * N_actions,
            0.0);
  std::fill(propensities.end() - r_ghost * N_actions, propensities.end(), 0.0);
  grid[0] = 0;
  grid[L + 1] = 0;
  ASSERT(grid[L + 2] == 0);
};
template <typename T> void fastTasep::BasicIteration<T>::fix_cumsum(int side) {
  int I = (side - 1) * N_actions;
  int Row = I / K;
  T R2 = sum_propensities[Row + 2];

  sum_propensities[Row + 1] =
      std::accumulate(propensities.size() + Row * K,
                      propensities.size() + (Row + 1) * K, sum_propensities[Row]);
  
  sum_propensities[Row + 2] =
      std::accumulate(propensities.size() + Row * K,
                      propensities.size() + (Row + 1) * K, sum_propensities[Row+1]);

Row= Row+2;
T diff = sum_propensities[Row]-R2;
  for(;Row<LEN_SUM_PROPENSITIES; Row++){
    sum_propensities[Row]+=diff;
  }
};

template <typename T> void fastTasep::BasicIteration<T>::set_index() {
  _INDEX = std::distance(
      sum_propensities.begin(),
      std::lower_bound(sum_propensities.begin(), sum_propensities.end(), r2));

  _index = _INDEX * K;
  temp_sum = sum_propensities[_INDEX];
  while (temp_sum < r2) {
    temp_sum += propensities[_index];
    _index++;
  }
};

template <typename T> void fastTasep::BasicIteration<T>::iteration() {
  r1 = dis(gen);
  r2 = dis(gen) * sum_propensities.back();
  dt = (1.0 / sum_propensities.back()) * log(1 / r1);

  set_index();

  _action = _index & 3;
  _side = _index >> 2;

  ASSERT(side != 0);
  ASSERT(side != L);
  ASSERT(side != L + 1);

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
template <typename T> void fastTasep::BasicIteration<T>::append_trajectory() {
  std::copy(grid.begin() + l_ghost, grid.end() - r_ghost,
            DATA.begin() + _iter * L);
  TIMES[_iter] = time;
  ACTION[_iter] = static_cast<u_int8_t>(_action);
  SIDE[_iter] = static_cast<u_int16_t>(_side - l_ghost);
};
template <typename T> void fastTasep::BasicIteration<T>::simulation() {
  while (_iter < ITERS) {
    append_trajectory();
    iteration();
    time += dt;
    _iter++;
  }
}

template class fastTasep::BasicIteration<float>;
template class fastTasep::BasicIteration<double>;