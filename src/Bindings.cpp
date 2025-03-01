#include "derived.hpp"
#include "timer.hpp"
#include "utils.hpp"

#ifdef TIME_ME

template <typename T>
// void new_iter_sim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
py::tuple new_iter_sim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    MyTimer timer;
    // fastTasep::BasicIteration<T> sim(L, ITERS, kon, koff, kstep, q, kq);

    // std::size_t mem_bytes = ((L + 3) * sizeof(char) + sizeof(T)) * ITERS;
    // // mem_bytes+= ITERS* sizeof(T);
    // void *memBlock = std::malloc(mem_bytes);
    // // if (!memBlock) {
    // //     std::throw << "Failed to allocate memory pool.\n";
    // //     // return 1;
    // // }
    // uint8_t *DATA_ptr = static_cast<uint8_t*>(std::malloc(ITERS*L * sizeof(uint8_t)));
    // py::capsule pool_capsule(DATA_ptr, [](void *m) { std::free(m); });

    // MemoryPool pool(memBlock, mem_bytes);
    fastTasep::BasicIteration<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    timer.add_lap();
    sim.simulation();
    timer.add_lap();

    // auto array1 = vector_to_numpy(DATA_ptr, ITERS, L, pool_capsule);
    // auto array2 = vector_to_numpy(std::move(sim.TIMES));
    // auto array3 = vector_to_numpy(std::move(sim.ACTION));
    // auto array4 = vector_to_numpy(std::move(sim.SIDE));
    // pool_capsule); auto array3 = vector_to_numpy(std::move(sim.ACTION), /*rows*/ ITERS, /*cols*/
    // 1, pool_capsule); auto array4 = vector_to_numpy(std::move(sim.SIDE), /*rows*/ ITERS, /*cols*/
    // 1, pool_capsule);

    const auto &to_rtn = py::make_tuple(
        vector_to_numpy(std::move(sim.DATA), ITERS, L), vector_to_numpy(std::move(sim.TIMES)),
        vector_to_numpy(std::move(sim.ACTION)), vector_to_numpy(std::move(sim.SIDE)));
    timer.add_lap();
    timer.print_times();

    // std::free(memBlock);

    // return py::make_tuple(array1, array2, array3, array4);
    return to_rtn;
}

template <typename T>
py::tuple kins_time_sim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    MyTimer timer;
    fastTasep::CountKins<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    timer.add_lap();

    sim.simulation();
    timer.add_lap();

    const auto &to_rtn =
        py::make_tuple(vector_to_numpy(std::move(sim.KINS)), vector_to_numpy(std::move(sim.TIMES)));
    timer.add_lap();
    timer.print_times();
    return to_rtn;
}

template <typename T>
py::tuple neighbors_sim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    MyTimer timer;
    fastTasep::Neighbors<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    timer.add_lap();

    sim.simulation();
    timer.add_lap();

    const auto &to_rtn = py::make_tuple(
        vector_to_numpy(std::move(sim.NEIGHBORS)), vector_to_numpy(std::move(sim.TIMES)),
        vector_to_numpy(std::move(sim.ACTION)), vector_to_numpy(std::move(sim.SIDE)));
    timer.add_lap();
    timer.print_times();

    return to_rtn;
}

template <typename T>
py::tuple nneighbors_sim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    MyTimer timer;
    fastTasep::NearestNeighbor<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    timer.add_lap();

    sim.simulation();
    timer.add_lap();

    const auto &to_rtn = py::make_tuple(
        vector_to_numpy(std::move(sim.NEIGHBORS)), vector_to_numpy(std::move(sim.TIMES)),
        vector_to_numpy(std::move(sim.ACTION)), vector_to_numpy(std::move(sim.SIDE)));
    timer.add_lap();
    timer.print_times();
    return to_rtn;
}

//============================================
#else
//============================================

template <typename T>
py::tuple new_iter_sim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    fastTasep::BasicIteration<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    sim.simulation();
    return py::make_tuple(
        vector_to_numpy(std::move(sim.DATA), ITERS, L), vector_to_numpy(std::move(sim.TIMES)),
        vector_to_numpy(std::move(sim.ACTION)), vector_to_numpy(std::move(sim.SIDE)));
}

template <typename T>
py::tuple kins_time_sim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    fastTasep::CountKins<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    sim.simulation();
    return py::make_tuple(vector_to_numpy(std::move(sim.KINS)),
                          vector_to_numpy(std::move(sim.TIMES)));
}

template <typename T>
py::tuple neighbors_sim(int L, int ITERS, T kon, T koff, T kstep, T q, T kq) {
    fastTasep::Neighbors<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    sim.simulation();
    return py::make_tuple(
        vector_to_numpy(std::move(sim.NEIGHBORS)), vector_to_numpy(std::move(sim.TIMES)),
        vector_to_numpy(std::move(sim.ACTION)), vector_to_numpy(std::move(sim.SIDE)));
}

template <typename T>
py::tuple nneighbors_sim(int L, int ITERS, T kon, T koff, T kstep, T q, double kq) {
    fastTasep::NearestNeighbor<T> sim(L, ITERS, kon, koff, kstep, q, kq);
    sim.simulation();
    return py::make_tuple(
        vector_to_numpy(std::move(sim.NEIGHBORS)), vector_to_numpy(std::move(sim.TIMES)),
        vector_to_numpy(std::move(sim.ACTION)), vector_to_numpy(std::move(sim.SIDE)));
}

#endif

//=================================
//=================================

PYBIND11_MODULE(tasep, m) {
    m.def("D_kins_time", &kins_time_sim<double>, "Returns only the total number of kins");
    m.def("F_kins_time", &kins_time_sim<float>, "Returns only the total number of kins");
    m.def("D_neighbors", &neighbors_sim<double>, "Returns only the total number of kins");
    m.def("F_neighbors", &neighbors_sim<float>, "Returns only the total number of kins");
    m.def("D_nneighbors", &nneighbors_sim<double>, "Returns only the total number of kins");
    m.def("F_nneighbors", &nneighbors_sim<float>, "Returns only the total number of kins");
    m.def("D_new_itersim", &new_iter_sim<double>, "Returns only the total number of kins");
    m.def("F_new_itersim", &new_iter_sim<float>, "Returns only the total number of kins");
}
