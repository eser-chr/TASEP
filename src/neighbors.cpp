#include "neighbors.h"
#include "timer.hpp"

#ifdef DEBUG
#define ASSERT(condition) assert(condition)
#else
#define ASSERT(condition)
#endif
template <typename T>
tasep::NeighborsBase<T>::NeighborsBase(int L, int ITERS, T kon, T koff, T kstep,
                                       T q, T kq)
    : L(L), ITERS(ITERS), kon(kon), koff(koff), kstep(kstep), q(q), kq(kq),
      gen(std::random_device{}()), dis(0.0, 1.0) {
  grid.resize(L + ghost, 0);
  propensities.resize(N_actions * (L + ghost), 0.0);
  sum_propensities.resize(N_actions * (L + ghost), 0.0);

  TIMES.resize(ITERS);
  ACTION.resize(ITERS);
  SIDE.resize(ITERS);
  for (int i = l_ghost * N_actions + BIND;
       i < propensities.size() - r_ghost * N_actions; i += N_actions) {
    propensities[i] = kon;
  }
  std::partial_sum(propensities.begin(), propensities.end(),
                   sum_propensities.begin());
}

template <typename T>
tasep::Neighbors<T>::Neighbors(int L, int ITERS, T kon, T koff, T kstep, T q,
                               T kq)
    : NeighborsBase<T>(L, ITERS, kon, koff, kstep, q, kq) {
  this->NEIGHBORS.reserve(2 * ITERS);
}

template <typename T>
tasep::NNeighbors<T>::NNeighbors(int L, int ITERS, T kon, T koff, T kstep, T q,
                                 T kq)
    : NeighborsBase<T>(L, ITERS, kon, koff, kstep, q, kq) {
  this->NEIGHBORS.reserve(ITERS);
}

// --------------------------------------------------------------------------

template <typename T> void tasep::NeighborsBase<T>::printme() {
  std::cout << "kon: " << kon << "\n";
  std::cout << "koff: " << koff << "\n";
  std::cout << "kstep: " << kstep << "\n";
  std::cout << "kq: " << kq << "\n";
  std::cout << "q: " << q << std::endl;
}

template <typename T> void tasep::NeighborsBase<T>::bind(int side) {
  ASSERT(grid[side] == 0);
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = 0.0;
  propensities[temp + BIND] = 0;
  propensities[temp + UNBIND] = koff;
  propensities[temp + STEP] = (double)kstep * (1 - grid[side + 1]);
  propensities[temp + DEACTIVATE] = 0.0;
  grid[side] = 1;
  TOTAL_KINS++;
}

template <typename T> void tasep::Neighbors<T>::bind(int side) {
NeighborsBase<T>::bind(side);
  RNN = 0;
  LNN = 0;
  lnn = static_cast<int16_t>(side - 1);
  rnn = static_cast<int16_t>(side + 1);

  while (lnn > this->l_ghost) {
    if (this->grid[lnn]) {
      LNN = lnn - side;
      break;
    }
    lnn--;
  }

  while (rnn < this->L + this->r_ghost) {
    if (this->grid[rnn]) {
      RNN = rnn - side;
      break;
    }
    rnn++;
  }
  
  this->NEIGHBORS.push_back(LNN);
  this->NEIGHBORS.push_back(RNN);
}

template <typename T>
void tasep::NNeighbors<T>::bind(int side) { 
NeighborsBase<T>::bind(side);
  NN = 0;
  lnn = side - 1;
  rnn = side + 1;

  while (lnn > this->l_ghost && rnn < this->L+ this->r_ghost) {
    if (this->grid[lnn]) {
      NN = lnn - side;
      break;
    };
    if (this->grid[rnn]) {
      NN = rnn - side;
      break;
    };

    lnn--;
    rnn++;
  }
  this->NEIGHBORS.push_back(NN);
}

template <typename T> void tasep::NeighborsBase<T>::unbind(int side) {
  ASSERT(grid[side] == 1);
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = (double)kstep * grid[side - 1];
  propensities[temp + BIND] = kon * q;
  propensities[temp + UNBIND] = 0.0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = kq;
  grid[side] = 0;
  TOTAL_KINS--;
}

template <typename T> void tasep::NeighborsBase<T>::step(int side) {
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

template <typename T> void tasep::NeighborsBase<T>::deactivate(int side) {
  ASSERT(grid[side] == 0);
  temp = side * N_actions;
  propensities[temp + BIND] = kon;
  propensities[temp + UNBIND] = 0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = 0;
}

template <typename T> void tasep::NeighborsBase<T>::fix_boundaries() {
  std::fill(propensities.begin(), propensities.begin() + l_ghost * N_actions,
            0.0);
  std::fill(propensities.end() - r_ghost * N_actions, propensities.end(), 0.0);
  grid[0] = 0;
  grid[L + 1] = 0;
  ASSERT(grid[L + 2] == 0);
};

template <typename T> void tasep::NeighborsBase<T>::fix_cumsum(int side) {
  for (int i = (side - 1) * N_actions; i < propensities.size(); i++)
    sum_propensities[i] = sum_propensities[i - 1] + propensities[i];
};

template <typename T> void tasep::NeighborsBase<T>::iteration() {
  r1 = dis(gen);
  r2 = dis(gen) * sum_propensities.back();
  dt = (1.0 / sum_propensities.back()) * log(1 / r1);

  _index = std::distance(
      sum_propensities.begin(),
      std::upper_bound(sum_propensities.begin(), sum_propensities.end(), r2));

  _action = _index & 3;
  _side = _index >> 2;

  ASSERT(_side != 0);
  ASSERT(_side != L);
  ASSERT(_side != L + 1);

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

template <typename T> void tasep::NeighborsBase<T>::append_trajectory() {

  TIMES[_iter] = time;
  ACTION[_iter] = _action;
  SIDE[_iter] = _side;
};

template <typename T> void tasep::NeighborsBase<T>::simulation() {
  while (_iter < ITERS) {
    append_trajectory();
    iteration();
    time += dt;
    _iter++;
  }
}

template class tasep::NeighborsBase<double>;
template class tasep::NeighborsBase<float>;
template class tasep::Neighbors<double>;
template class tasep::Neighbors<float>;
template class tasep::NNeighbors<double>;
template class tasep::NNeighbors<float>;
