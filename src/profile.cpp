#include "derived.hpp"

template <typename T>
fastTasep::Profile<T>::Profile(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq), DATA(L) {}

template <typename T>
void fastTasep::Profile<T>::append_trajectory() {
    if (this->_time >= this->Tequil) {
        for (size_t i = this->l_ghost; i < this->L + this->r_ghost; ++i) {
            DATA[i] += static_cast<double>(this->grid[i]) * this->dt;
        }
        ttotal += this->dt;
    }
}

template <typename T>
py::tuple fastTasep::Profile<T>::export_python() {
    for (auto &d : DATA) d /= ttotal;
    return py::make_tuple(vector_to_numpy(std::move(DATA)));
}

template class fastTasep::Profile<double>;
template class fastTasep::Profile<float>;