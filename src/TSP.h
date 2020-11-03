#ifndef KTH_AA_TSP_H
#define KTH_AA_TSP_H

#include <iostream>
#include <vector>
#include <chrono>

#include "utility.h"

using namespace std;
using namespace chrono;

class TSP {
public:
    static vector<int> travel(const vector<pair<double, double>> &cities);

    static vector<int> travel_naive(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> travel_nn(vector<pair<double, double>> cities, Grid<int> &distances);

    static vector<int> travel_cw(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> travel_cw_seq(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> local_2opt(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knearest,
                                  time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

    static vector<int> local_2opt_no_knn(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                                         time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

    static vector<int> local_3opt(vector<int> tour, Grid<int> &distances,
                                  time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

    static int tour_distance(const vector<pair<double, double>> &cities, vector<int> tour);

    static vector<pair<double, double>> create_n_cities(int n, int seed = 72);

};

#endif //KTH_AA_TSP_H
