#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// using namespace std;
namespace py = pybind11;

template <typename T>
py::array_t<T> vector_to_numpy(std::vector<T> &&vec) {
  py::array_t<T> result({vec.size()}, {sizeof(T)});
  auto result_buffer = result.request();
  T *result_ptr = static_cast<T *>(result_buffer.ptr);
  std::move(vec.begin(), vec.end(), result_ptr);
  return result;
}

template <typename T>
py::array_t<T> vector_to_numpy(std::vector<T> &&vec, size_t rows,
                               size_t cols) {
  py::array_t<T> result({rows, cols});
  auto result_buffer = result.request();
  T *result_ptr = static_cast<T *>(result_buffer.ptr);
  std::move(vec.begin(), vec.end(), result_ptr);

  return result;
}

template <typename T, typename Allocator>
py::array_t<T> vector_to_numpy(std::vector<T, Allocator>&& vec,
                               size_t rows, size_t cols,
                               py::handle base) {
    T* data_ptr = vec.data();
    return py::array_t<T>(
        {rows, cols},                  // shape
        {cols * sizeof(T), sizeof(T)}, // strides (row-major)
        data_ptr,                      // pointer to data
        base                           // base object: the capsule
    );
}

template <typename T>
py::array_t<T> vector_to_numpy(T* data_ptr, 
                               size_t rows, size_t cols,
                               py::handle base) {
    // T* data_ptr = vec.data();
    return py::array_t<T>(
        {rows, cols},                  // shape
        {cols * sizeof(T), sizeof(T)}, // strides (row-major)
        data_ptr,                      // pointer to data
        base                           // base object: the capsule
    );
}


