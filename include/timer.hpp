#pragma once
#include <chrono>

class MyTimer {
  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point end;
  //   std::chrono::milliseconds end;

public:
  MyTimer() : start(std::chrono::high_resolution_clock::now()) {}
  double get() {
    end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
  }
};