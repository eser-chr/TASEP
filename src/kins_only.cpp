#include "derived.hpp"

template <typename T>
fastTasep::CountKins<T>::CountKins(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq) {
    KINS.resize(ITERS);
    TIMES.resize(ITERS);
}

template <typename T>
void fastTasep::CountKins<T>::append_trajectory() {
    KINS[this->_iter] = TOTAL_KINS;
    TIMES[this->_iter] = this->_time;
};

template<typename T>
void fastTasep::CountKins<T>::bind(int side){
  AbstractIteration<T>::bind(side);
  TOTAL_KINS++;
}

template<typename T>
void fastTasep::CountKins<T>::unbind(int side){
  AbstractIteration<T>::unbind(side);
  TOTAL_KINS--;
}

template<typename T>
void fastTasep::CountKins<T>::step(int side){
  AbstractIteration<T>::step(side);
  if(side == this->L)TOTAL_KINS--;
}

template <typename T>
py::tuple fastTasep::CountKins<T>::export_python() {
    return py::make_tuple(
        vector_to_numpy(std::move(KINS)), 
        vector_to_numpy(std::move(TIMES))
    );
}

template class fastTasep::CountKins<double>;
template class fastTasep::CountKins<float>;
