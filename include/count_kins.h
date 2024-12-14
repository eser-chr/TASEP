#pragma once

#include "timer.hpp"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric> // For partial sum
#include <random>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace tasep {

template <typename T> class CountKins {

public:
  CountKins(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
  void simulation();
  void printme();

  std::vector<uint16_t> KINS;
  std::vector<T> TIMES;

protected:
  int L;
  int ITERS;
  T kon; // _ _ _ + _ _   ->   _ _ _ + + _
  T koff;
  T kstep; // _ _ _ + _ _   ->   _ _ _ _ + _
  T q;
  T kq;
  uint16_t TOTAL_KINS = 0;
  std::vector<T> propensities;
  std::vector<T> sum_propensities;
  std::vector<bool> grid;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<T> dis;

  int _action, _side, _index, temp;
  T r1, r2, dt;
  T time = 0.0;
  int _iter = 0;

  enum ACTION { BIND = 0, UNBIND = 1, STEP = 2, DEACTIVATE = 3 };
  const int l_ghost = 1;
  const int r_ghost = 2;
  const int ghost = l_ghost + r_ghost;
  const int N_actions = 4;

private:
  void bind(int side);
  void unbind(int side);
  void step(int side);
  void deactivate(int side);

  void fix_boundaries();
  void fix_cumsum(int side);

  void iteration();
  void append_trajectory();
};

}; // namespace tasep


namespace fastTasep {

template <typename T> class CountKins {

public:
  CountKins(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
  void simulation();
  void printme();

  std::vector<uint16_t> KINS;
  std::vector<T> TIMES;

protected:
  int L; // 

  int ITERS;
  T kon; // _ _ _ + _ _   ->   _ _ _ + + _
  T koff;
  T kstep; // _ _ _ + _ _   ->   _ _ _ _ + _
  T q;
  T kq;
  std::vector<T> propensities;
  
  int ROWS, COLS;
  std::vector<T> sum_of_rows;
  std::vector<T> cumsum_of_rows;
  
  std::vector<bool> grid;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<T> dis;

  uint16_t TOTAL_KINS = 0;

  int _action, _side, _index, _INDEX, temp;
  int _iter = 0;
  T r1, r2, dt;
  T time = 0.0;

  enum ACTION { BIND = 0, UNBIND = 1, STEP = 2, DEACTIVATE = 3 };
  const int l_ghost = 1;
  int r_ghost = 2;
  const int ghost = l_ghost + r_ghost;
  const int N_actions = 4;

private:
  void set_index();
  void bind(int side);
  void unbind(int side);
  void step(int side);
  void deactivate(int side);

  void fix_boundaries();
  void fix_cumsum(int side);

  void iteration();
  void append_trajectory();
};

}; // end namespace fastTasep
