#pragma once

#include <memory>
#include <random>
#include <vector>

#include "bucket/bucket.h"
#include "debug_utils.hpp"

namespace fastTasep {
template <typename T>
class AbstractIteration {
   public:
    AbstractIteration(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
    virtual ~AbstractIteration() = default;
    void simulation();
    void printme();

   protected:
    int L;  //
    int ITERS;
    T kon;  // _ _ _ + _ _   ->   _ _ _ + + _
    T koff;
    T kstep;  // _ _ _ + _ _   ->   _ _ _ _ + _
    T q;
    T kq;
    std::vector<T> propensities;

    int ROWS, COLS;
    double Nactions_per_col;
    std::unique_ptr<bucket<T>> _manager;  // manages the cumsum effectively
    std::vector<uint8_t> grid;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<T> dis;

    int _action, _side, _index, temp;
    size_t _iter = 0;
    T r1, r2, dt;
    T time = 0.0;

    enum ACTION { BIND = 0, UNBIND = 1, STEP = 2, DEACTIVATE = 3 };
    const int l_ghost = 1;
    const int r_ghost = 2;
    const int ghost = l_ghost + r_ghost;
    const int N_actions = 4;
    virtual void bind(int side);
    virtual void unbind(int side);
    virtual void step(int side);
    virtual void deactivate(int side);
    using action_func = void (AbstractIteration<T>::*)(int);

    static const std::array<action_func, 4> actions;

   private:
    void fixBoundaries();
    void fixCumsum(int side);
    inline void executeAction(int action, int side) { (this->*actions[action])(side); }

    void iteration();
    virtual void append_trajectory() = 0;
};


template <typename T>
const std::array<typename AbstractIteration<T>::action_func, 4> AbstractIteration<T>::actions = {
    &AbstractIteration<T>::bind,
    &AbstractIteration<T>::unbind,
    &AbstractIteration<T>::step,
    &AbstractIteration<T>::deactivate
};


};  // namespace fastTasep