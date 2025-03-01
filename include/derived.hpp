#pragma once

#include <memory>
#include <random>
#include <span>
#include <vector>

#include "abstract.hpp"
#include "bucket/bucket.h"
#include "debug_utils.hpp"

namespace fastTasep {

template <typename T>
class BasicIteration : public AbstractIteration<T> {
   public:
    BasicIteration(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
    std::vector<uint8_t> DATA;
    std::vector<T> TIMES;
    py::tuple export_python();

   private:
    void append_trajectory();
};

template <typename T>
class CountKins : public AbstractIteration<T> {
   public:
    CountKins(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
    std::vector<uint16_t> KINS;
    std::vector<T> TIMES;
    py::tuple export_python();

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
    std::vector<uint16_t> Left;
    std::vector<uint16_t> Right;
    py::tuple export_python();

   private:
    void bind(int side) override;
    void append_trajectory();
};

template <typename T>
class NearestNeighbor : public AbstractIteration<T> {
   public:
    NearestNeighbor(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);
    std::vector<uint16_t> NEIGHBORS;
    py::tuple export_python();

   private:
    void bind(int side) override;
    void append_trajectory();
};

};  // end namespace fastTasep

// template <typename T>
// class Profile : public AbstractIteration<T> {
//    public:
//     Profile(int L, int ITERS, T kon, T koff, T kstep, T q, T kq);

//     std::vector<uint64_t> DATA;
//     size_t iter_equil = 0;
//     T t_equil = 0.0;
//     T t_final = 0.0;
//         void export_python();

//    private:
//     bool is_equil_set = false;
//     T T_equil();
//     void append_trajectory();
// };