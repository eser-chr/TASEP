#include <cstring>

#include "derived.hpp"


template <typename T>
fastTasep::BasicIteration<T>::BasicIteration(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq) {
    // DATA = std::span<uint8_t>(ptr, ITERS*L);

    // std::cout<<DATA.size()<<std::endl;
    // DATA.resize(ITERS);
    DATA.resize(L * ITERS);
    TIMES.resize(ITERS);
    ACTION.resize(ITERS);
    SIDE.resize(ITERS);
}

// --------------------------------------------------------------------------

template <typename T>
void fastTasep::BasicIteration<T>::append_trajectory() {
    std::copy(this->grid.begin() + this->l_ghost, this->grid.end() - this->r_ghost,
              DATA.begin() + this->_iter * this->L);

    TIMES[this->_iter] = this->time;
    ACTION[this->_iter] = static_cast<u_int8_t>(this->_action);
    SIDE[this->_iter] = static_cast<u_int16_t>(this->_side - this->l_ghost);
};

template class fastTasep::BasicIteration<float>;
template class fastTasep::BasicIteration<double>;