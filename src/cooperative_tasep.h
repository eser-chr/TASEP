#ifndef COOP_TASEP
#define COOP_TASEP

#include <random>
#include <vector>
#include <numeric> // For partial sum
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <tuple>
#include <random>

class CooperativeTasep{
public:

CooperativeTasep(int L, int T, double kon, double koff, double kstep, double q, double kq);

virtual void bind(int side);
virtual void unbind(int side);
virtual void step(int side);
virtual void deactivate(int side);

virtual void fix_boundaries();
virtual void fix_cumsum();

virtual void iteration();
virtual void append_to_data();

virtual void simulation();

std::tuple< std::vector<std::vector<int>>, 
            std::vector<double> 
            > get_results();





protected:
    int L; // with 2 ghost layers
    // int N; //population
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
    enum ACTION{BIND=0, UNBIND=1, STEP=2, DEACTIVATE=3};

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
    const int __UNBIND = UNBIND*(L+ghost) +l_ghost;
    const int __STEP = STEP*(L+ghost) +l_ghost;
    const int __DEACTIVATE = DEACTIVATE*(L+ghost) +l_ghost;

    int _action, _side, _index;
    double r1, r2, dt;
};


class SpecificCooperativeTasep:CooperativeTasep{
public:
    SpecificCooperativeTasep(int L, int T, double kon, double koff, double kstep, double q, double kq);
    void append_to_data()override;
    void bind(int side)override;
    void simulation()override;

/*
Returns 
DATA, ACTIVATION, nearest_neighbor, TIMES, res, dts

 */
std::tuple< std::vector<std::vector<int>>, 
            std::vector<std::vector<int>>, 
            std::vector<std::vector<int>>, 
            std::vector<double>, 
            std::vector<int>,
            std::vector<double> 
            > get_results();

private:
    std::vector<std::vector<int>> ACTIVATION;
    std::vector<std::vector<int>> nearest_neighbor;
    std::vector<int> res;
    std::vector<double> dts;
    int left_ptr, right_ptr;
};

#endif