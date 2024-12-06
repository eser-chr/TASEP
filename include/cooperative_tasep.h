#ifndef COOP_TASEP
#define COOP_TASEP

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

template <typename T> class Basic {

public:
  Basic(int L, T TIME, T kon, T koff, T kstep, T q, T kq, bool trajectory,
        bool details, T period);
  virtual void simulation();

  std::tuple<std::vector<uint8_t>, std::vector<uint16_t>, std::vector<T>>
  get_details();
  std::tuple<std::vector<uint8_t>, std::vector<T>> get_trajectory();

protected:
  int L;
  T TIME;
  T kon; // _ _ _ + _ _   ->   _ _ _ + + _
  T koff;
  T kstep; // _ _ _ + _ _   ->   _ _ _ _ + _
  T q;
  T kq;
  std::vector<T> propensities;
  std::vector<T> sum_propensities;
  std::vector<bool> grid;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<T> dis;

  int _action, _side, _index, temp;
  T r1, r2, dt;
  T time = 0.0;

  enum ACTION { BIND = 0, UNBIND = 1, STEP = 2, DEACTIVATE = 3 };
  const int l_ghost = 1;
  const int r_ghost = 2;
  const int ghost = l_ghost + r_ghost;
  const int N_actions = 4;

  void (Basic::*analysis_function)();
  bool trajectory, details;
  T period;
  T next_write_time = 0.0;
  int counter = 0;
  std::vector<uint8_t> res_actions;
  std::vector<uint16_t> res_sides;
  std::vector<T> res_dts;
  std::vector<T> TIMES;
  std::vector<uint8_t> DATA;

private:
  virtual void bind(int side);
  virtual void unbind(int side);
  virtual void step(int side);
  virtual void deactivate(int side);

  virtual void fix_boundaries();
  virtual void fix_cumsum(int side);

  virtual void iteration();
  virtual void append_details();
  virtual void append_trajectory();
  virtual void append_all();
};

template <typename T> class BasicIteration {

public:
  BasicIteration(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
  void simulation();
  void printme();

  std::vector<uint8_t> DATA;
  std::vector<T> TIMES;
  std::vector<uint8_t> ACTION;
  std::vector<uint16_t> SIDE;

protected:
  int L;
  int ITERS;
  T kon; // _ _ _ + _ _   ->   _ _ _ + + _
  T koff;
  T kstep; // _ _ _ + _ _   ->   _ _ _ _ + _
  T q;
  T kq;
  std::vector<T> propensities;
  std::vector<T> sum_propensities;
  std::vector<uint8_t> grid;

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

#endif