#pragma once
#include <random>
#include <pcg_random.hpp>

namespace fastTasep {

template <typename T>
class BaseRNG {
   public:
    virtual T random() = 0;
    virtual ~BaseRNG() = default;
};

template <typename T>
class MTRNG : public BaseRNG<T> {
   private:
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<T> dis;

   public:
    MTRNG() : gen(std::random_device{}()), dis(0.0, 1.0) {}
    inline T random() override { return dis(gen); }
};

template <typename T>
class PCGRNG : public BaseRNG<T> {
   private:
    pcg32 gen;
    std::uniform_real_distribution<T> dis;

   public:
    PCGRNG()
        : gen{pcg_extras::seed_seq_from<std::random_device>()},  // or any seeding strategy
          dis(static_cast<T>(0.0), static_cast<T>(1.0)) {}

    inline T random() override { return dis(gen); }
};
}  // namespace fastTasep