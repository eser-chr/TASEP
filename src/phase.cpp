#include "derived.hpp"

template <typename T>
fastTasep::Profile<T>::Profile(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq) {
    DATA.resize(L);
    t_equil = T_equil();
}

template <typename T>
void fastTasep::Profile<T>::append_trajectory() {
    if (this->time >= t_equil) {
        if(!is_equil_set) {
            iter_equil = this->_iter;
            is_equil_set = true;
        }
        for (int i = 0; i < this->L; i++) {
            DATA[i] += this->dt*this->grid[i];
        }
    }
    if(this->_iter == this->ITERS - 1) {
        t_final = this->time;
    }

}

template <typename T>
T fastTasep::Profile<T>::T_equil() {
    return 1.0 / (this->kon + this->koff);
}


template class fastTasep::Profile<double>;
template class fastTasep::Profile<float>;