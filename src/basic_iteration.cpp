#include <cstring>

#include "derived.hpp"


template <typename T>
fastTasep::BasicIteration<T>::BasicIteration(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq) {

    DATA.resize(L * ITERS);
    TIMES.resize(ITERS);

}

// --------------------------------------------------------------------------

template <typename T>
void fastTasep::BasicIteration<T>::append_trajectory() {
    std::copy(this->grid.begin() + this->l_ghost, this->grid.end() - this->r_ghost,
              DATA.begin() + this->_iter * this->L);

    TIMES[this->_iter] = this->time;

};

template <typename T>
py::tuple fastTasep::BasicIteration<T>::export_python() {
    return py::make_tuple(
        vector_to_numpy(std::move(DATA), this->ITERS, this->L),
        vector_to_numpy(std::move(TIMES)));
}


template class fastTasep::BasicIteration<float>;
template class fastTasep::BasicIteration<double>;