#pragma once
#include "timer.hpp"
#include "bucket/bucket.h"
#include <iostream>
#include <random>
#include <stdexcept>
#include <memory>
#include <vector>

namespace fastTasep {

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
  int L; // 
  int ITERS;
  T kon; // _ _ _ + _ _   ->   _ _ _ + + _
  T koff;
  T kstep; // _ _ _ + _ _   ->   _ _ _ _ + _
  T q;
  T kq;
  std::vector<T> propensities;
  
  int ROWS, COLS;
  std::unique_ptr<bucket<T>> _manager; // manages the cumsum effectively
  std::vector<uint8_t> grid;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<T> dis;

  int _action, _side, _index, temp;
  int _iter = 0;
  T r1, r2, dt;
  T time = 0.0;

  enum ACTION { BIND = 0, UNBIND = 1, STEP = 2, DEACTIVATE = 3 };
  const int l_ghost = 1;
  int r_ghost = 2;
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

}; // end namespace fastTasep