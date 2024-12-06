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

template <typename T> class NeighborsBase {

public:
  NeighborsBase(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
  void simulation();
  void printme();

  //   std::vector<uint16_t> KINS;
  std::vector<int16_t> NEIGHBORS;
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

protected:
  virtual void bind(int side); // Only difference
  void unbind(int side);
  void step(int side);
  void deactivate(int side);

  void fix_boundaries();
  void fix_cumsum(int side);

  void iteration();
  void append_trajectory();
};

template <typename T> class Neighbors : public NeighborsBase<T> {
public:
  Neighbors(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);

private:
  int RNN, LNN, rnn, lnn;
  void bind(int side) override;
};

template <typename T> class NNeighbors : public NeighborsBase<T> {
public:
  NNeighbors(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);

private:
  int NN, rnn, lnn;
  void bind(int side) override;
};

}; // namespace tasep
