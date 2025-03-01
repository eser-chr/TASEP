#include "derived.hpp"

template <typename T>
fastTasep::NearestNeighbor<T>::NearestNeighbor(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq), NEIGHBORS(L) {}

template <typename T>
void fastTasep::NearestNeighbor<T>::bind(int side) {
    AbstractIteration<T>::bind(side);
    bool isleft = false;
    bool isright = false;
    uint16_t lnn = side - 1;
    uint16_t rnn = side + 1;

    while (lnn > this->l_ghost && rnn < this->L + this->r_ghost) {
        if (this->grid[lnn]) {
            isleft = true;
            break;
        };
        if (this->grid[rnn]) {
            isright = true;
            break;
        };

        lnn--;
        rnn++;
    }
    if (isleft)
        NEIGHBORS[side - lnn]++;
    else if (isright)
        NEIGHBORS[rnn - side]++;
}

template <typename T>
void fastTasep::NearestNeighbor<T>::append_trajectory() {};

template <typename T>
py::tuple fastTasep::NearestNeighbor<T>::export_python() {
    return py::make_tuple(vector_to_numpy(std::move(NEIGHBORS)));
}

template class fastTasep::NearestNeighbor<double>;
template class fastTasep::NearestNeighbor<float>;
