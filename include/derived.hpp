#pragma once

#include <memory>
#include <random>
#include <span>
#include <vector>

#include "abstract.hpp"
#include "bucket/bucket.h"
#include "debug_utils.hpp"
#include "memory_pool.hpp"

namespace fastTasep {

template <typename T>
class BasicIteration : public AbstractIteration<T> {
   public:
    BasicIteration(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
    // BasicIteration(int L, int ITERS, T kon, T koff, T kstep, T q, T kq, MemoryPool* pool);

    std::vector<uint8_t> DATA;
    // std::span<uint8_t> DATA;
    std::vector<T> TIMES;
    std::vector<uint8_t> ACTION;
    std::vector<uint16_t> SIDE;

   private:
    void append_trajectory();
};
template <typename T>
class Profile : public AbstractIteration<T> {
   public:
    Profile(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);

    std::vector<uint64_t> DATA;
    size_t iter_equil = 0;
    T t_equil = 0.0;
    T t_final = 0.0;

   private:
    bool is_equil_set = false;
    T T_equil();
    void append_trajectory();
};

template <typename T>
class CountKins : public AbstractIteration<T> {
   public:
    CountKins(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
    std::vector<uint16_t> KINS;
    std::vector<T> TIMES;

   private:
    uint16_t TOTAL_KINS = 0;
    void bind(int side) override;
    void unbind(int side) override;
    void step(int side) override;
    void append_trajectory();
};

template <typename T>
class Neighbors : public AbstractIteration<T> {
   public:
    Neighbors(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
    std::vector<int16_t> NEIGHBORS;
    std::vector<T> TIMES;
    std::vector<uint8_t> ACTION;
    std::vector<uint16_t> SIDE;

   private:
    int16_t RNN, LNN, rnn, lnn;
    void bind(int side) override;
    void append_trajectory();
};

template <typename T>
class NearestNeighbor : public AbstractIteration<T> {
   public:
    NearestNeighbor(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
    std::vector<int16_t> NEIGHBORS;
    std::vector<T> TIMES;
    std::vector<uint8_t> ACTION;
    std::vector<uint16_t> SIDE;

   private:
    int16_t NN, rnn, lnn;
    void bind(int side) override;
    void append_trajectory();
};

};  // end namespace fastTasep