#include "derived.hpp"

template <typename T>
fastTasep::NearestNeighbor<T>::NearestNeighbor(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq) {
    NEIGHBORS.reserve(ITERS);
    TIMES.resize(ITERS);
}

template <typename T>
void fastTasep::NearestNeighbor<T>::bind(int side) {
    AbstractIteration<T>::bind(side);
    NN = 0;
    lnn = side - 1;
    rnn = side + 1;

    while (lnn > this->l_ghost && rnn < this->L + this->r_ghost) {
        if (this->grid[lnn]) {
            NN = lnn - side;
            break;
        };
        if (this->grid[rnn]) {
            NN = rnn - side;
            break;
        };

        lnn--;
        rnn++;
    }
    if (NN != 0) this->NEIGHBORS.push_back(NN);
}

template <typename T>
void fastTasep::NearestNeighbor<T>::append_trajectory() {
    TIMES[this->_iter] = this->time;
};

template <typename T>
py::tuple fastTasep::NearestNeighbor<T>::export_python() {
    return py::make_tuple(vector_to_numpy(std::move(NEIGHBORS)),
                          vector_to_numpy(std::move(TIMES)));
}

template class fastTasep::NearestNeighbor<double>;
template class fastTasep::NearestNeighbor<float>;
