#include "cooperative_tasep.h"
#include "timer.hpp"

#ifdef DEBUG
#define ASSERT(condition) assert(condition)
#else
#define ASSERT(condition)
#endif
using DATATYPE = double;
tasep::BasicIteration::BasicIteration(int L, int ITERS, DATATYPE kon,
                                      DATATYPE koff, DATATYPE kstep, DATATYPE q,
                                      DATATYPE kq)
    : L(L), ITERS(ITERS), kon(kon), koff(koff), kstep(kstep), q(q), kq(kq),
      gen(std::random_device{}()), dis(0.0, 1.0) {


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

void tasep::BasicIteration::printme() {
  std::cout << "kon: " << kon << "\n";
  std::cout << "koff: " << koff << "\n";
  std::cout << "kstep: " << kstep << "\n";
  std::cout << "kq: " << kq << "\n";
  std::cout << "q: " << q << std::endl;
}

void tasep::BasicIteration::bind(int side) {

  ASSERT(grid[side] == 0);
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = 0.0;
  propensities[temp + BIND] = 0;
  propensities[temp + UNBIND] = koff;
  propensities[temp + STEP] = (double)kstep * (1 - grid[side + 1]);
  propensities[temp + DEACTIVATE] = 0.0;
  grid[side] = 1;
}

void tasep::BasicIteration::unbind(int side) {
  ASSERT(grid[side] == 1);
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = (double)kstep * grid[side - 1];
  propensities[temp + BIND] = kon * q;
  propensities[temp + UNBIND] = 0.0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = kq;
  grid[side] = 0;
}

void tasep::BasicIteration::step(int side) {
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

void tasep::BasicIteration::deactivate(int side) {
  ASSERT(grid[side] == 0);
  temp = side * N_actions;
  propensities[temp + BIND] = kon;
  propensities[temp + UNBIND] = 0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = 0;
}

void tasep::BasicIteration::fix_boundaries() {
  std::fill(propensities.begin(), propensities.begin() + l_ghost * N_actions,
            0.0);
  std::fill(propensities.end() - r_ghost * N_actions, propensities.end(), 0.0);
  grid[0] = 0;
  grid[L + 1] = 0;
  ASSERT(grid[L + 2] == 0);
};

void tasep::BasicIteration::fix_cumsum(int side) {
  for (int i = (side - 1) * N_actions; i < propensities.size(); i++)
    sum_propensities[i] = sum_propensities[i - 1] + propensities[i];
};

void tasep::BasicIteration::iteration() {
  r1 = dis(gen);
  r2 = dis(gen) * sum_propensities.back();
  dt = (1.0 / sum_propensities.back()) * log(1 / r1);

  _index = std::distance(
      sum_propensities.begin(),
      std::upper_bound(sum_propensities.begin(), sum_propensities.end(), r2));

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

void tasep::BasicIteration::append_trajectory() {
  std::copy(grid.begin() + l_ghost, grid.end() - r_ghost,
            DATA.begin() + _iter * L);
  TIMES[_iter] = time;
  ACTION[_iter] = static_cast<u_int8_t>(_action);
  SIDE[_iter] = static_cast<u_int16_t>(_side - l_ghost);
};

void tasep::BasicIteration::simulation() {
  while (_iter < ITERS) {
    append_trajectory();
    iteration();
    time += dt;
    _iter++;
  }
}

// std::tuple<std::vector<uint8_t>, std::vector<DATATYPE>,
// std::vector<u_int8_t>,
//            std::vector<uint16_t>>
// tasep::BasicIteration::get_trajectory() {
//   std::cout << std::flush;
//   std::vector<u_int8_t> action_subvector(
//       std::make_move_iterator(ACTION.begin() + 1),
//       std::make_move_iterator(ACTION.end()));
//   std::vector<u_int16_t> side_subvector(
//       std::make_move_iterator(SIDE.begin() + 1),
//       std::make_move_iterator(SIDE.end()));

//   return std::make_tuple(DATA, TIMES, action_subvector, side_subvector);
// };
