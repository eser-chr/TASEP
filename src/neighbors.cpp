#include "derived.hpp"

template <typename T>
fastTasep::Neighbors<T>::Neighbors(int L, int ITERS, T kon, T koff, T kstep, T q, T kq)
    : AbstractIteration<T>(L, ITERS, kon, koff, kstep, q, kq), Left(L), Right(L) {}

template <typename T>
void fastTasep::Neighbors<T>::bind(int side) {
    AbstractIteration<T>::bind(side);

    bool isleft = false;
    bool isright = false;

    uint16_t lnn = (side - 1);
    uint16_t rnn = (side + 1);

    while (lnn > this->l_ghost) {
        if (this->grid[lnn]) {
            isleft = true;
            break;
        }
        lnn--;
    }

    while (rnn < this->L + this->r_ghost) {
        if (this->grid[rnn]) {
            isright = true;
            break;
        }
        rnn++;
    }
    if (isleft && isright) {
        ASSERT(LNN < 0, "LNN out of boundaries");
        ASSERT(RNN > 0, "RNN out of boundaries");
        Left[side - lnn]++;
        Right[rnn - side]++;
    }
}

template <typename T>
void fastTasep::Neighbors<T>::append_trajectory() {};

template <typename T>
py::tuple fastTasep::Neighbors<T>::export_python() {
    return py::make_tuple(vector_to_numpy(std::move(Left)), vector_to_numpy(std::move(Right)));
}

template class fastTasep::Neighbors<double>;
template class fastTasep::Neighbors<float>;
