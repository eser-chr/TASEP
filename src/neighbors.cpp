#include "derived.hpp"

template <typename T>
fastTasep::Neighbors<T>::Neighbors(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq) {
    NEIGHBORS.reserve(2 * ITERS);
    TIMES.resize(ITERS);
    ACTION.resize(ITERS);
    SIDE.resize(ITERS);
}




template <typename T>
void fastTasep::Neighbors<T>::bind(int side) {
    AbstractIteration<T>::bind(side);
    RNN = 0;
    LNN = 0;
    lnn = (side - 1);
    rnn = (side + 1);

    while (lnn > this->l_ghost) {
        if (this->grid[lnn]) {
            LNN = lnn - side;
            break;
        }
        lnn--;
    }

    while (rnn < this->L + this->r_ghost) {
        if (this->grid[rnn]) {
            RNN = rnn - side;
            break;
        }
        rnn++;
    }
    if (LNN != 0 && RNN != 0) {
        ASSERT(LNN < 0, "LNN out of boundaries");
        ASSERT(RNN > 0, "RNN out of boundaries");
        this->NEIGHBORS.push_back(LNN);
        this->NEIGHBORS.push_back(RNN);
    }
}



template <typename T>
void fastTasep::Neighbors<T>::append_trajectory() {
    TIMES[this->_iter] = this->time;
    ACTION[this->_iter] = this->_action;
    SIDE[this->_iter] = this->_side;
};


template class fastTasep::Neighbors<double>;
template class fastTasep::Neighbors<float>;
