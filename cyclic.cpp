#include <iostream>
#include <vector>
#include <numeric>
#include <random>
#include <cmath>
#include <tuple>
#include <functional>

// using namespace std;

std::tuple< std::vector<std::vector<int>>, std::vector<double>, std::vector<int>, std::vector<double> > 
cyclic_simulation(int L, double T, double kon, double koff, double kstep, double kq, double q, int initial_density) {
    std::cout<<initial_density<<std::endl;
    std::vector<double> propensities(4 * L, 0.0);
    std::fill(propensities.begin(), propensities.begin()+L, kon);
    // propensities[0] = kon;
    // propensities[1] = kon * q;
    // propensities[L - 1] = koff;
    // propensities[2 * L - 1] = koff * kstep / (koff + kstep);
    // propensities[3 * L - 1] = kq;
    
    double theo_dt = (1.0 / (L * (q * kon + koff + kstep + kq))) * log(1.01);
    int theo_iter = static_cast<int>(T / theo_dt);
    
    std::vector<std::vector<int>> DATA;
    std::vector<double> TIMES;
    std::vector<int> res;
    std::vector<double> dts;
    
    double period = 0.1;
    int blocks = static_cast<int>(T / period);
    
    int i = 0;
    double t = 0;
    int block = 0;
    double next_write_time = period;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);


    // Populate the grid stochastically
    for(int side=0; side<L; side++){
        double r = dis(gen);
        if(100*r<initial_density){
            int front = (side + 1) % L;
            // int ffront = (side + 2) % L;
            int back = (side - 1 + L) % L; // Ensure back is non-negative
            propensities[side] = 0;
            propensities[L + side] = koff;
            propensities[2 * L + back] = 0;
            propensities[2 * L + side] = (1 - propensities[L + front] / koff) * kstep;
            propensities[3 * L + side] = 0;
        }
    }
    
    while (t < T) {
        std::vector<double> S(4 * L);
        partial_sum(propensities.begin(), propensities.end(), S.begin());
        double r1 = dis(gen);
        double r2 = S.back() * dis(gen);
        double dt = (1.0 / S.back()) * log(1 / r1);
        int index = std::distance(S.begin(), std::upper_bound(S.begin(), S.end(), r2));
       int action = index / L;
        int side = index % L;
        int front = (side + 1) % L;
        int ffront = (side + 2) % L;
        int back = (side - 1 + L) % L; // Ensure back is non-negative

        switch (action) {
            case 0:
                propensities[side] = 0;
                propensities[L + side] = koff;
                propensities[2 * L + back] = 0;
                propensities[2 * L + side] = (1 - propensities[L + front] / koff) * kstep;
                propensities[3 * L + side] = 0;
                break;
            case 1:
                propensities[side] = q * kon;
                propensities[L + side] = 0;
                propensities[2 * L + side] = 0;
                propensities[2 * L + back] = kstep * (propensities[L + back]) / koff;
                propensities[3 * L + side] = kq;
                break;
            case 2:
                propensities[side] = q * kon;
                propensities[front] = 0;
                propensities[L + side] = 0;
                propensities[L + front] = koff;
                propensities[2 * L + back] = kstep * (propensities[L + back]) / koff;
                propensities[2 * L + side] = 0;
                propensities[2 * L + front] = (1 - propensities[L + ffront] / koff) * kstep;
                propensities[3 * L + side] = kq;
                propensities[3 * L + front] = 0;
                break;
            default:
                propensities[side] = kon;
                propensities[3 * L + side] = 0;
        }

        
        // propensities[0] = 0;
        // propensities[L] = 0;
        // propensities[2 * L] = 0;
        // propensities[3 * L] = 0;
        // propensities[L - 1] = 0;
        // propensities[2 * L - 1] = koff;
        // propensities[3 * L - 1] = 0;
        // propensities[4 * L - 1] = 0;
        
        if (next_write_time < t) {
            std::vector<int> data;
            for (int i = 0; i < L; ++i)
                data.push_back(propensities[i] == 0 ? 1 : 0);
            DATA.push_back(data);
            TIMES.push_back(t);
            ++block;
            next_write_time += period;
        }
        
        if (i == theo_iter) {
            res.resize(static_cast<int>(theo_iter * (T / t)));
            dts.resize(static_cast<int>(theo_iter * (T / t)));
        }
        
        res.push_back(action);
        dts.push_back(dt);
        t += dt;
        ++i;
    }
    
    return std::make_tuple(DATA, TIMES, res, dts);
}



