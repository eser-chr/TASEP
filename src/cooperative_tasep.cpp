#include "cooperative_tasep.h"

CooperativeTasep::CooperativeTasep(int L, int T, double kon, double koff,
                                   double kstep, double q, double kq)
    : L(L), T(T), kon(kon), koff(koff), kstep(kstep), q(q), kq(kq),
      gen(std::random_device{}()), dis(0.0, 1.0) {


  propensities.resize(N_actions * (L + ghost), 0);
  sum_propensities.resize(N_actions * (L + ghost));
  std::fill(propensities.begin() + l_ghost, propensities.begin() + L - r_ghost,
            kon);
  std::partial_sum(propensities.begin(), propensities.end(),
                   sum_propensities.begin());

  grid.resize(L + ghost, 0);


}

SpecificCooperativeTasep::SpecificCooperativeTasep(int L, int T, double kon,
                                                   double koff, double kstep,
                                                   double q, double kq)
    : CooperativeTasep(L, T, kon, koff, kstep, q, kq) {}

void CooperativeTasep::bind(int side) {
  propensities[__BIND + side] = 0;

  propensities[__UNBIND + side] = koff;

  propensities[__STEP + side] = (double)kstep * (!grid[l_ghost + side + 1]);
  propensities[__STEP + side - 1] = 0.0;

  propensities[__DEACTIVATE + side] = 0.0;

  grid[l_ghost + side] = 1;
}

void CooperativeTasep::unbind(int side) {
  propensities[__BIND + side] = kon * q;

  propensities[__UNBIND + side] = 0.0;

  propensities[__STEP + side] = 0.0;
  propensities[__STEP + side - 1] = (double)kstep * grid[l_ghost + side - 1];

  propensities[__DEACTIVATE + side] = kq;

  grid[l_ghost + side] = 0;
}

void CooperativeTasep::step(int side) {
  propensities[__BIND + side] = kon * q;
  propensities[__BIND + side + 1] = 0.0;

  propensities[__UNBIND + side] = 0;
  propensities[__UNBIND + side + 1] = koff;

  propensities[__STEP + side - 1] = (double)kstep * grid[l_ghost + side - 1];
  propensities[__STEP + side] = 0.0;
  propensities[__STEP + side + 1] = (double)kstep * (!grid[l_ghost + side + 2]);

  propensities[__DEACTIVATE + side] = kq;
  propensities[__DEACTIVATE + side + 1] = 0.0;

  grid[l_ghost + side] = 0;
  grid[l_ghost + side + 1] = 1;
}
void CooperativeTasep::deactivate(int side) {
  propensities[__BIND + side] = kon;

  propensities[__UNBIND + side] = 0;

  propensities[__STEP + side - 1] = (double)kstep * grid[l_ghost + side - 1];
  propensities[__STEP + side] = 0.0;

  propensities[__DEACTIVATE + side] = 0;
}

void CooperativeTasep::fix_boundaries() {
  for (int left_ghost = 0; left_ghost < l_ghost; left_ghost++) {
    propensities[__BIND - l_ghost + left_ghost] = 0;
    propensities[__UNBIND - l_ghost + left_ghost] = 0;
    propensities[__STEP - l_ghost + left_ghost] = 0;
    propensities[__DEACTIVATE - l_ghost + left_ghost] = 0;
  }

  for (int right_ghost = 0; right_ghost < r_ghost; right_ghost++) {
    propensities[__BIND + L + right_ghost] = 0;
    propensities[__UNBIND + L + right_ghost] = 0;
    propensities[__STEP + L + right_ghost] = 0;
    propensities[__DEACTIVATE + L + right_ghost] = 0;
  }
}

void CooperativeTasep::fix_cumsum() {
  std::partial_sum(propensities.begin(), propensities.end(),
                   sum_propensities.begin());
}

void CooperativeTasep::iteration() { // Possible std::async
  r1 = dis(gen);
  dt = (1.0 / sum_propensities.back()) * log(1 / r1);
  time += dt;

  r2 = dis(gen) * sum_propensities.back();

  _index = std::distance(
      sum_propensities.begin(),
      std::upper_bound(sum_propensities.begin(), sum_propensities.end(), r2));


  // result = div(_index, (L + ghost));
  // _action = result.quot;
  // _side = result.rem;
  _action = _index / (L + ghost);
  _side = _index % (L + ghost);

  // _action = _index & 3;
  // _side = _index >>2 ;



  // Thread 1 -----

  switch (_action) {
  case 0:
    bind(_side - l_ghost);
    break;
  case 1:
    unbind(_side - l_ghost);
    break;
  case 2:
    step(_side - l_ghost);
    break;
  case 3:
    deactivate(_side - l_ghost);
    break;
  default:
    std::cout << _action << " " << _side << "\n";
    for (int i = 10; i < 16; i++)
      std::cout << grid[i] << " ";
    std::cout << std::endl;

    for (int i = 10; i < 16; i++)
      std::cout << propensities[i] << " ";
    std::cout << std::endl;

    for (int i = 10; i < 16; i++)
      std::cout << propensities[L + ghost + i] << " ";
    std::cout << std::endl;
    for (int i = 10; i < 16; i++)
      std::cout << propensities[2 * (L + ghost) + i] << " ";
    std::cout << std::endl;
    for (int i = 10; i < 16; i++)
      std::cout << propensities[3 * (L + ghost) + i] << " ";
    std::cout << std::endl;

    throw std::runtime_error("An error occurred in switch statement");
  }

  fix_boundaries();
  // Thread 2
  fix_cumsum();
}

void CooperativeTasep::append_to_data() { // Need to resize at the constructor
                                          // to avoid rellocation
  DATA.emplace_back(L);
  std::transform(grid.begin() + __BIND, grid.begin() + __BIND + L,
                 DATA.back().begin(),
                 [](bool grid_val) { return (int)grid_val; });   // Consider initialize grid as int.

  TIMES.push_back(time);
}

void CooperativeTasep::simulation() {
  while (time < T) {
    iteration();
    if (next_write_time < time) {
      append_to_data();
      next_write_time += writing_period;
    }
  }
}

void SpecificCooperativeTasep::append_to_data() {
  CooperativeTasep::append_to_data();
  ACTIVATION.emplace_back(L);
  std::transform(propensities.begin() + __DEACTIVATE,
                 propensities.begin() + __DEACTIVATE + L,
                 ACTIVATION.back().begin(),
                 [](int propensity) { return propensity == 0 ? 0 : 1; });
}

void SpecificCooperativeTasep::bind(int side) {
  CooperativeTasep::bind(side);

  left_ptr = -1;
  right_ptr = 1;

  while (side + left_ptr > l_ghost) {
    if (grid[l_ghost + side + left_ptr]) {
      break;
    }
    left_ptr--;
  }
  // std::cout<<side<<" ";
  while (side + right_ptr < L + l_ghost) {
    if (grid[l_ghost + side + right_ptr]) {
      break;
    }
    right_ptr++;
  }

  if (side + left_ptr != l_ghost &&
      side + right_ptr !=
          L + l_ghost) { // Count only real neighbors and not bdrs
    std::vector<int> temp = {left_ptr, right_ptr};
    nearest_neighbor.push_back(temp);
  }
}

void SpecificCooperativeTasep::simulation() {
  while (time < T) {
    iteration();
    res.push_back(_action);
    sides.push_back(_side);
    dts.push_back(dt);
    if (next_write_time < time) {
      append_to_data();
      next_write_time += writing_period;
    }
  }
}

std::tuple<std::vector<std::vector<int>>, std::vector<double>>
CooperativeTasep::get_results() {
  std::cout << std::flush;
  return std::make_tuple(DATA, TIMES);
  // return std::make_tuple(DATA, ACTIVATION, TIMES, res, dts);
}

std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>,
           std::vector<std::vector<int>>, std::vector<double>, std::vector<int>,
           std::vector<int>, std::vector<double>>
SpecificCooperativeTasep::get_results() {
  std::cout << std::flush;
  return std::make_tuple(DATA, ACTIVATION, nearest_neighbor, TIMES, res, sides,
                         dts);
}
