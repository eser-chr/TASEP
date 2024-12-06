#include "cooperative_tasep.h"

template <typename T>
tasep::Basic<T>::Basic(int L, T TIME, T kon, T koff, T kstep, T q, T kq,
                       bool trajectory, bool details, T period = 0.1)
    : L(L), TIME(TIME), kon(kon), koff(koff), kstep(kstep), q(q), kq(kq),
      trajectory(trajectory), details(details), period(period), time(0.0),
      gen(std::random_device{}()), dis(0.0, 1.0) {

  grid.resize(L + ghost, 0);
  propensities.resize(N_actions * (L + ghost), 0);
  sum_propensities.resize(N_actions * (L + ghost));

  for (int i = l_ghost * N_actions + BIND;
       i < propensities.size() - r_ghost * N_actions; i += N_actions) {
    propensities[i] = kon;
  }

  std::partial_sum(propensities.begin(), propensities.end(),
                   sum_propensities.begin());

  if (!trajectory && !details) {
    throw std::runtime_error("You have to specify analytics!");
  }

  if (trajectory) {
    DATA.resize(L * ceil(TIME / period));
    TIMES.reserve(ceil(TIME / period));
    // analysis_function = &Basic<T><T>::append_trajectory;
  }
  if (details) {
    res_actions.reserve(10000);
    res_sides.reserve(10000);
    res_dts.reserve(10000); // Heuristic choice.
    // analysis_function = &Basic<T>::append_details;
  }
  if (trajectory && !details) {
    analysis_function = &Basic<T>::append_trajectory;
  } else if (!trajectory && details) {
    analysis_function = &Basic<T>::append_details;
  } else if (trajectory && details) {
    analysis_function = &Basic<T>::append_all;
  }
}

template <typename T> void tasep::Basic<T>::bind(int side) {
  // if (grid[side]) { throw std::runtime_error("BIND"); }
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = 0.0;
  propensities[temp + BIND] = 0;
  propensities[temp + UNBIND] = koff;
  propensities[temp + STEP] = (double)kstep * (!grid[side + 1]);
  propensities[temp + DEACTIVATE] = 0.0;
  grid[side] = 1;
}

template <typename T> void tasep::Basic<T>::unbind(int side) {
  // if (!grid[side]) { throw std::runtime_error("UNBIND"); }
  temp = side * N_actions;
  propensities[temp + STEP - N_actions] = (double)kstep * grid[side - 1];
  propensities[temp + BIND] = kon * q;
  propensities[temp + UNBIND] = 0.0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = kq;
  grid[side] = 0;
}

template <typename T> void tasep::Basic<T>::step(int side) {
  temp = side * N_actions;
  // if (!grid[side]) { throw std::runtime_error("STEP2"); }
  // if (grid[side+1]) { throw std::runtime_error("STEP1"); }
  //  DYNAMIC_ASSERT(!grid[side + 1], "step");
  //  DYNAMIC_ASSERT(grid[side], "step2");
  propensities[temp + STEP - N_actions] = (double)kstep * grid[side - 1];

  propensities[temp + BIND] = kon * q;
  propensities[temp + UNBIND] = 0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = kq;

  propensities[temp + BIND + N_actions] = 0.0;
  propensities[temp + UNBIND + N_actions] = koff;
  propensities[temp + STEP + N_actions] = (double)kstep * (!grid[side + 2]);
  propensities[temp + DEACTIVATE + N_actions] = 0.0;
  grid[side] = 0;
  grid[side + 1] = 1;
}

template <typename T> void tasep::Basic<T>::deactivate(int side) {
  // DYNAMIC_ASSERT(!grid[side], "deactivate");
  // if (grid[side]) { throw std::runtime_error("DEACTIVATE"); }
  temp = side * N_actions;
  propensities[temp + BIND] = kon;
  propensities[temp + UNBIND] = 0;
  propensities[temp + STEP] = 0.0;
  propensities[temp + DEACTIVATE] = 0;
  propensities[temp + STEP - N_actions] = (double)kstep * grid[side - 1];
}

template <typename T> void tasep::Basic<T>::fix_boundaries() {
  std::fill(propensities.begin(), propensities.begin() + l_ghost * N_actions,
            0.0);
  std::fill(propensities.end() - r_ghost * N_actions, propensities.end(), 0.0);
  grid[0] = 0;
  grid[L + 1] = 0;
  grid[L + 2] = 0;
};

template <typename T> void tasep::Basic<T>::fix_cumsum(int side) {
  for (int i = (side - 1) * N_actions; i < propensities.size(); i++)
    sum_propensities[i] = sum_propensities[i - 1] + propensities[i];

  // DYNAMIC_ASSERT(sum_propensities.back()>0, "cumsum");
};

template <typename T> void tasep::Basic<T>::iteration() {
  r1 = dis(gen);
  r2 = dis(gen) * sum_propensities.back();
  dt = (1.0 / sum_propensities.back()) * log(1 / r1);

  _index = std::distance(
      sum_propensities.begin(),
      std::upper_bound(sum_propensities.begin(), sum_propensities.end(), r2));

  // DYNAMIC_ASSERT(_index > 4, "index");
  // DYNAMIC_ASSERT(_index< 4*(L+3)-8, "index");  /// To be removed!!

  _action = _index & 3;
  _side = _index >> 2; // Thats because N_actions is power of 2, e.g 4 = 2^2

  // if(_side<1 || _side > L+1){throw std::runtime_error("Side is off");}

  // Thread 1 -----

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
  // Thread 2
  fix_cumsum(_side);
};

template <typename T> void tasep::Basic<T>::append_details() {
  res_actions.push_back(_action);
  res_dts.push_back(dt);
  res_sides.push_back(_side);
};
template <typename T> void tasep::Basic<T>::append_trajectory() {
  if (next_write_time < time) {
    if (counter * L < DATA.size()) {

      std::transform(grid.begin() + l_ghost, grid.begin() + l_ghost + L,
                     DATA.begin() + counter * L, [](bool grid_val) {
                       return (uint8_t)grid_val;
                     }); // Consider initialize grid as int.
      TIMES.push_back(time);
      counter++;
      next_write_time += period;
    } else {
      std::cout << counter << " vs " << ceil(TIME / period) << std::endl;
      std::cout << "time: " << time << "\n";
      throw std::runtime_error("Counter excedded Data");
    }
  }
};

template <typename T> void tasep::Basic<T>::append_all() {
  append_details();
  append_trajectory();
}

template <typename T> void tasep::Basic<T>::simulation() {
  while (time < TIME) {
    iteration();
    (this->*analysis_function)();
    time += dt;
  }
}

template <typename T>
std::tuple<std::vector<uint8_t>, std::vector<uint16_t>, std::vector<T>>
tasep::Basic<T>::get_details() {
  std::cout << std::flush;
  return std::make_tuple(res_actions, res_sides, res_dts);
};

template <typename T>
std::tuple<std::vector<uint8_t>, std::vector<T>>
tasep::Basic<T>::get_trajectory() {
  std::cout << std::flush;
  return std::make_tuple(DATA, TIMES);
};




template class tasep::Basic<float>;
template class tasep::Basic<double>;