#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#define ENABLE_CHECKS

#ifdef ENABLE_CHECKS
#define ROW_CHECK(cond, msg)                                                   \
  if (!(cond)) throw std::runtime_error(msg);
#define VAL_CHECK(cond, msg)                                                   \
  if (!(cond)) throw std::runtime_error(msg);
#else
#define ROW_CHECK(cond, msg)
#define VAL_CHECK(cond, msg)
#endif

template <typename T> class bucket_ref {

private:
  std::vector<T> foo;

public:
  int _ROWS;
  int _COLS;
  int _size;
  const std::vector<T> &_vector;
  std::vector<T> _p_sums;
  std::vector<T> _p_cum_sums;
  using value_type = T;

//   bucket_ref() : _ROWS(0), _COLS(0), _vector(foo) {}

//   void reset(int ROWS, int COLS, const std::vector<T> &other) {
//     _ROWS = ROWS;
//     _COLS = COLS;
//     _vector = other;
//     _p_sums.resize(_ROWS);
//     _p_cum_sums.resize(_ROWS + 1);
//     update_sum();
//     update_cumsum();
//   }

  bucket_ref(int ROWS, int COLS, const std::vector<T> &other)
      : _ROWS(ROWS), _COLS(COLS), _vector(other) {
    _size = ROWS * _COLS;
    _p_sums.resize(_ROWS);
    _p_cum_sums.resize(_ROWS + 1);
    update_sum();
    update_cumsum();
  }

  struct sums {
    auto begin() { return _p_sums.begin(); }
    auto end() { return _p_sums.end(); }
    auto begin() const { return _p_sums.begin(); }
    auto end() const { return _p_sums.end(); }
    auto cbegin() const { return _p_sums.cbegin(); }
    auto cend() const { return _p_sums.cend(); }
  };
  class cumsums {
    auto begin() { return _p_cum_sums.begin(); }
    auto end() { return _p_cum_sums.end(); }
    auto begin() const { return _p_cum_sums.begin(); }
    auto end() const { return _p_cum_sums.end(); }
    auto cbegin() const { return _p_cum_sums.cbegin(); }
    auto cend() const { return _p_cum_sums.cend(); }
    void print() {
      for (const T i : _p_cum_sums) std::cout << i << ",";
      std::cout << std::endl;
    }
  };
  //----------------------------------------------
  void print() {
    for (const T i : _p_cum_sums) std::cout << i << ",";
    std::cout << std::endl;
  }

  void update_sum() {
    for (size_t row = 0; row < _ROWS; row++) update_sum_at_row(row);
  }

  void update_sum_at_row(int row) {
    ROW_CHECK(row >= 0 && row < _ROWS, "Row index out of range");
    _p_sums[row] = static_cast<T>(0);
    for (size_t i = _COLS * row; i < _COLS * (row + 1); i++) {
      _p_sums[row] += _vector[i];
    }
  }

  void update_cumsum() {
    _p_cum_sums[0] = static_cast<T>(0);

    for (size_t row = 0; row < _ROWS; row++) {
      _p_cum_sums[row + 1] = _p_cum_sums[row] + _p_sums[row];
    }
  }

  void update_cumsum_from_row(int row) {
    ROW_CHECK(row >= 0 && row < _ROWS, "Row index out of range");
    for (size_t l_row = row; l_row < _ROWS; l_row++) {
      _p_cum_sums[l_row + 1] = _p_cum_sums[l_row] + _p_sums[l_row];
    }
  }

  void update_single_row(int row) {
    ROW_CHECK(row >= 0 && row < _ROWS, "Row index out of range");
    T diff = _p_sums[row];
    update_sum_at_row(row);
    diff -= _p_sums[row];

    for (size_t l_row = row + 1; l_row < _ROWS + 1; l_row++) {
      _p_cum_sums[l_row] -= diff;
    }
  }

  void update_single_row_of_index(int index) {
    update_single_row(index / _COLS);
  }

  int find_upper_bound_in_cumsum(const T &val) const {
    auto it = std::upper_bound(_p_cum_sums.begin(), _p_cum_sums.end(), val);
    return static_cast<int>(std::distance(_p_cum_sums.begin(), it));
  }

  int find_upper_bound(const T &val) const {
    VAL_CHECK(
        val > 0,
        "In upper limit, the value passed is smaller than the first element")
    VAL_CHECK(val < _p_cum_sums.back(), "In upper limit, the value passed is "
                                        "bigger or equal to the last element")

    int row_index = std::distance(_p_cum_sums.begin(),
                                  std::upper_bound(_p_cum_sums.begin(),
                                                   _p_cum_sums.end(), val)) -
                    1;
    // if(row_index<0){
    //     return -1;
    // }else if(row_>=_ROWS){
    //     r = _ROWS -1;
    // }
    int index = row_index * _COLS;
    T temp = _p_cum_sums[row_index];
    // while (temp < val){
    //     temp+=_vector[index];
    //     index++;
    // }

    for (; index < (row_index + 1) * _COLS; index++) {
      temp += _vector[index];
      if (temp >= val) break;
    }
    if (index >= (row_index + 1) * _COLS) return -1;
    return index;
  }
};
