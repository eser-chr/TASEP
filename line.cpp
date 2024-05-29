#include <iostream>
#include <vector>
#include <numeric>
#include <random>
#include <cmath>
#include <tuple>
#include <functional>
#include <chrono>

// using namespace std;

std::tuple< std::vector<std::vector<int>>, std::vector<double>, std::vector<int>, std::vector<double> > 
line_simulation(int L, double T, double kon, double koff, double kstep, double kq, double q) {

    std::vector<double> propensities(4 * L, 0.0);
    std::fill(propensities.begin()+1, propensities.begin()+L-1, kon);

    // propensities[L - 1] = koff;
    propensities[L - 1] = koff;
    
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
    
    while (t < T) {
        std::vector<double> S(4 * L);
        partial_sum(propensities.begin(), propensities.end(), S.begin());
        double r1 = dis(gen);
        double r2 = S.back() * dis(gen);
        double dt = (1.0 / S.back()) * log(1 / r1);
        int index = std::distance(S.begin(), std::upper_bound(S.begin(), S.end(), r2));
        int action = index / L;
        int side = index % L;
        
        switch(action){
            case 0:
                propensities[side] = 0;
                propensities[L + side] = koff;
                propensities[2 * L + side - 1] = 0;
                propensities[2 * L + side] = (1 - propensities[L + side + 1] / koff) * kstep;
                propensities[3 * L + side] = 0;
                break;
            case 1:
                propensities[side] = q * kon;
                propensities[L + side] = 0;
                propensities[2 * L + side] = 0;
                propensities[2 * L + side - 1] = kstep * (propensities[L + side - 1]) / koff;
                propensities[3 * L + side] = kq;
                break;
            case 2:
                propensities[side] = q * kon;
                propensities[side + 1] = 0;
                propensities[L + side] = 0;
                propensities[L + side + 1] = koff;
                propensities[2 * L + side - 1] = kstep * (propensities[L + side - 1]) / koff;
                propensities[2 * L + side] = 0;
                propensities[2 * L + side + 1] = (1 - propensities[L + side + 2] / koff) * kstep;
                propensities[3 * L + side] = kq;
                propensities[3 * L + side + 1] = 0;
                break;
            case 3:
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



struct side{
    double kon;
    double koff;
    double kstep;
    double kq;

    side(double kon) : kon(kon), koff(0), kstep(0), kq(0) {}
};

// void binding_action(std::vector<side>& propensities, int side, int koff, int kstep)

std::tuple< std::vector<std::vector<int>>, std::vector<double>, std::vector<int>, std::vector<double> > 
line_simulation_b(int L, double T, double kon, double koff, double kstep, double kq, double q) {

    // std::vector<double> propensities(4 * L, );
    std::vector<side> propensities(L, side(kon));


    // std::fill(propensities.begin()+1, propensities.begin()+L-1, kon);

    // propensities[L - 1] = koff;
    propensities[L - 1].koff = koff;
    
    // double theo_dt = (1.0 / (L * (q * kon + koff + kstep + kq))) * log(1.01);
    // int theo_iter = static_cast<int>(T / theo_dt);
    
    std::vector<std::vector<int>> DATA;
    std::vector<double> TIMES;
    std::vector<int> res;
    std::vector<double> dts;
    
    const double period = 0.1;
    int blocks = static_cast<int>(T / period);
    
    int i = 0;
    double t = 0;
    int block = 0;
    double next_write_time = period;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    while (t < T) {
        std::vector<double> S(4 * L);
        partial_sum(propensities.begin(), propensities.end(), S.begin());
        double r1 = dis(gen);
        double r2 = S.back() * dis(gen);
        double dt = (1.0 / S.back()) * log(1 / r1);
        int index = std::distance(S.begin(), std::upper_bound(S.begin(), S.end(), r2));
        int action = index / L;
        int side = index % L;
        
        switch(action){
            case 0: //Binding
                propensities[side].kon = 0;
                propensities[side].koff = koff;
                propensities[side - 1].kstep = 0;
                propensities[side].kstep = (1 - propensities[side + 1].koff / koff) * kstep;
                propensities[side].kq = 0;
                break;
            case 1:
                propensities[side].kon = q * kon;
                propensities[side].koff = 0;
                propensities[side].kstep = 0;
                propensities[side - 1].kstep = kstep * (propensities[side - 1].koff) / koff;
                propensities[side].kq = kq;
                break;
            case 2: // Stepping
                propensities[side].kon = q * kon;
                propensities[side + 1].kon = 0;
                propensities[side].koff = 0;
                propensities[side + 1].koff = koff;
                propensities[side - 1].kstep = kstep * (propensities[side - 1].koff) / koff;
                propensities[side].kstep = 0;
                propensities[side + 1].kstep = (1 - propensities[side + 2].koff / koff) * kstep;
                propensities[side].kq = kq;
                propensities[side + 1].kq = 0;
                break;
            case 3:
                propensities[side] = kon;
                propensities[3 * L + side] = 0;        
        }
        
        propensities[0].kon = 0;
        propensities[0].koff = 0;
        propensities[0].kstep = 0;
        propensities[0].kq = 0;

        propensities[L - 1].kon = 0;
        propensities[L - 1].koff = 0;
        propensities[L - 1].kstep = 0;
        propensities[L - 1].kq = 0;
        

        if (next_write_time < t) {
            std::vector<int> data;
            for (int i = 0; i < L; ++i)
                data.push_back(propensities[i].kon == 0 ? 1 : 0);
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


int main() {
    // Parameters
    int L = 1000;
    double T = 10.0;
    double kon = 0.1;
    double koff = 0.2;
    double kstep = 0.05;
    double kq = 0.01;
    double q = 0.5;

    // Timing for line_simulation
    auto start_time = std::chrono::high_resolution_clock::now();
    line_simulation(L, T, kon, koff, kstep, kq, q);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time taken by line_simulation: " << duration.count() << " milliseconds" << std::endl;

    // Timing for line_simulation_b
    start_time = std::chrono::high_resolution_clock::now();
    line_simulation_b(L, T, kon, koff, kstep, kq, q);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time taken by line_simulation_b: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}
