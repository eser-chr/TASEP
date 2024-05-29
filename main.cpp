#include <iostream>
#include <vector>
#include <cmath>
#include <random>

typedef std::vector<std::vector<bool>> Matrix;
typedef std::vector<double> Vector;

std::tuple<Matrix, Vector, Vector, Vector> simulation(int L, double T, double kon, double koff, double kstep, double kq, double q) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    std::vector<double> propensities(4 * L, 0.0);
    for (int i = 0; i < L; ++i) {
        propensities[i] = kon;
    }

    double theo_dt = (1.0 / L * (q * kon + koff + kstep + kq)) * log(1.01);
    int theo_iter = int(T / theo_dt);
    Vector res(theo_iter, 0.0);
    Vector dts(theo_iter, 0.0);

    double period = 0.1;
    int blocks = int(T / period);

    int i = 0;
    double t = 0;
    int block = 0;
    double next_write_time = period;
    Matrix DATA(blocks, std::vector<bool>(L, false));
    Vector TIMES(blocks, 0.0);

    while (t < T) {
        double sum_propensities = 0.0;
        for (double prop : propensities) {
            sum_propensities += prop;
        }
        double r1 = dis(gen);
        double r2 = sum_propensities * dis(gen);

        double dt = (1.0 / sum_propensities) * log(1 / r1);
        int index = 0;
        double cumulative_sum = 0.0;
        for (int j = 0; j < propensities.size(); ++j) {
            cumulative_sum += propensities[j];
            if (cumulative_sum >= r2) {
                index = j;
                break;
            }
        }
        int action = index / L;
        int side = index % L;

        if (action == 0) {
            propensities[side] = 0;
            propensities[L + side] = koff;
            propensities[2 * L + side - 1] = 0;
            propensities[2 * L + side] = (1 - propensities[L + side + 1] / koff) * kstep;
            propensities[3 * L + side] = 0;
        } else if (action == 1) {
            propensities[side] = q * kon;
            propensities[L + side] = 0;
            propensities[2 * L + side] = 0;
            propensities[2 * L + side - 1] = kstep * (propensities[L + side - 1]) / koff;
            propensities[3 * L + side] = kq;
        } else if (action == 2) {
            propensities[side] = q * kon;
            propensities[side + 1] = 0;
            propensities[L + side] = 0;
            propensities[L + side + 1] = koff;
            propensities[2 * L + side - 1] = kstep * (propensities[L + side - 1]) / koff;
            propensities[2 * L + side] = 0;
            propensities[2 * L + side + 1] = (1 - propensities[L + side + 2] / koff) * kstep;
            propensities[3 * L + side] = kq;
            propensities[3 * L + side + 1] = 0;
        } else {
            propensities[side] = kon;
            propensities[3 * L + side] = 0;
        }

        propensities[0] = 0;
        propensities[L] = 0;
        propensities[2 * L] = 0;
        propensities[3 * L] = 0;
        propensities[L - 1] = 0;
        propensities[2 * L - 1] = koff;
        propensities[3 * L - 1] = 0;
        propensities[4 * L - 1] = 0;

        if (next_write_time < t) {
            for (int j = 0; j < L; ++j) {
                DATA[block][j] = (propensities[j] == 0);
            }
            TIMES[block] = t;
            block++;
            next_write_time += period;
        }

        if (i == theo_iter) {
            int new_elements_count = int(theo_iter * (T / t));
            res.insert(res.end(), new_elements_count, 0.0);
            dts.insert(dts.end(), new_elements_count, 0.0);
        }

        res[i] = action;
        dts[i] = dt;
        t += dt;
        i++;
    }

    return std::make_tuple(DATA, TIMES, res, dts);
}

