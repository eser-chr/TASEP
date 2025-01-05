#include "new.h"

template <typename T>
fastTasep::CountKins<T>::CountKins(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq) {
    KINS.resize(ITERS);
    TIMES.resize(ITERS);
}

template <typename T>
void fastTasep::CountKins<T>::append_trajectory() {
    KINS[this->_iter] = TOTAL_KINS;
    TIMES[this->_iter] = this->time;
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

template class fastTasep::CountKins<double>;
template class fastTasep::CountKins<float>;
