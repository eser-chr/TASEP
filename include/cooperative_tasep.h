#ifndef COOP_TASEP
#define COOP_TASEP

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

class Basic {
  using DATATYPE = double;

public:
  Basic(int L, DATATYPE T, DATATYPE kon, DATATYPE koff, DATATYPE kstep,
        DATATYPE q, DATATYPE kq, bool trajectory, bool details,
        DATATYPE period);
  virtual void simulation();

  std::tuple<std::vector<uint8_t>, std::vector<uint16_t>, std::vector<DATATYPE>>
  get_details();
  std::tuple<std::vector<uint8_t>, std::vector<DATATYPE>> get_trajectory();
  // std::tuple<std::vector<std::vector<uint8_t>>, std::vector<float>>
  // get_trajectory();

protected:
  int L;
  DATATYPE T;
  DATATYPE kon; // _ _ _ + _ _   ->   _ _ _ + + _
  DATATYPE koff;
  DATATYPE kstep; // _ _ _ + _ _   ->   _ _ _ _ + _
  DATATYPE q;
  DATATYPE kq;
  std::vector<DATATYPE> propensities;
  std::vector<DATATYPE> sum_propensities;
  std::vector<bool> grid;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<DATATYPE> dis;

  int _action, _side, _index, temp;
  DATATYPE r1, r2, dt;
  DATATYPE time = 0.0;

  enum ACTION { BIND = 0, UNBIND = 1, STEP = 2, DEACTIVATE = 3 };
  const int l_ghost = 1;
  const int r_ghost = 2;
  const int ghost = l_ghost + r_ghost;
  const int N_actions = 4;

  void (Basic::*analysis_function)();
  bool trajectory, details;
  DATATYPE period;
  DATATYPE next_write_time = 0.0;
  int counter = 0;
  std::vector<uint8_t> res_actions;
  std::vector<uint16_t> res_sides;
  std::vector<DATATYPE> res_dts;
  // std::vector<std::vector<uint8_t>> DATA;
  std::vector<DATATYPE> TIMES;
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

  // void simulation_details();
  // void simulation_trajectory();
  // void simulation_details_trajectory();
};

class BasicIteration {
  using DATATYPE = double;

public:
  BasicIteration(int L, int ITERS, DATATYPE kon, DATATYPE koff, DATATYPE kstep,
                 DATATYPE q, DATATYPE kq);
  void simulation();
  void printme();
  // std::tuple<std::vector<uint8_t>, std::vector<DATATYPE>> get_trajectory();
  // std::tuple<std::vector<uint8_t>, std::vector<DATATYPE>, std::vector<u_int8_t>,
  //            std::vector<uint16_t>>
  // get_trajectory();

  std::vector<uint8_t> DATA;
  std::vector<DATATYPE> TIMES;
  std::vector<uint8_t> ACTION;
  std::vector<uint16_t> SIDE;
protected:
  int L;
  int ITERS;
  DATATYPE kon; // _ _ _ + _ _   ->   _ _ _ + + _
  DATATYPE koff;
  DATATYPE kstep; // _ _ _ + _ _   ->   _ _ _ _ + _
  DATATYPE q;
  DATATYPE kq;
  std::vector<DATATYPE> propensities;
  std::vector<DATATYPE> sum_propensities;
  std::vector<uint8_t> grid;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<DATATYPE> dis;

  int _action, _side, _index, temp;
  DATATYPE r1, r2, dt;
  DATATYPE time = 0.0;
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

class CooperativeTasep {
public:
  CooperativeTasep(int L, int T, double kon, double koff, double kstep,
                   double q, double kq);

  virtual void bind(int side);
  virtual void unbind(int side);
  virtual void step(int side);
  virtual void deactivate(int side);

  virtual void fix_boundaries();
  virtual void fix_cumsum();

  virtual void iteration();
  virtual void append_to_data();

  virtual void simulation();

  std::tuple<std::vector<std::vector<int>>, std::vector<double>> get_results();

protected:
  int L;
  int T;
  double kon;
  double koff;
  double kstep;
  double q;
  double kq;
  std::vector<double> propensities;
  std::vector<double> sum_propensities;
  std::vector<bool> grid;
  std::vector<std::vector<int>> DATA;
  std::vector<double> TIMES;
  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<double> dis;
  double next_write_time = writing_period;
  enum ACTION { BIND = 0, UNBIND = 1, STEP = 2, DEACTIVATE = 3 };

  double time = 0;
  // int BOUND;
  // int UNBOUND;
  const int l_ghost = 1;
  const int r_ghost = 2;
  const int ghost = l_ghost + r_ghost;
  const int N_actions = 4;
  const double writing_period = 0.1;

  // Displacement
  const int __BIND = BIND + l_ghost;
  const int __UNBIND = UNBIND * (L + ghost) + l_ghost;
  const int __STEP = STEP * (L + ghost) + l_ghost;
  const int __DEACTIVATE = DEACTIVATE * (L + ghost) + l_ghost;

  div_t result;
  int _action, _side, _index;
  double r1, r2, dt;
};

class SpecificCooperativeTasep : CooperativeTasep {
public:
  SpecificCooperativeTasep(int L, int T, double kon, double koff, double kstep,
                           double q, double kq);
  void append_to_data() override;
  void bind(int side) override;
  void simulation() override;

  /*
  Returns
  DATA, ACTIVATION, nearest_neighbor, TIMES, res, sides,  dts

   */
  std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>,
             std::vector<std::vector<int>>, std::vector<double>,
             std::vector<int>, std::vector<int>, std::vector<double>>
  get_results();

private:
  std::vector<std::vector<int>> ACTIVATION;
  std::vector<std::vector<int>> nearest_neighbor;
  std::vector<int> res;
  std::vector<int> sides;
  std::vector<double> dts;
  int left_ptr, right_ptr;
};

#endif