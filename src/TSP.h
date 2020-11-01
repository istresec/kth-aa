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
    static vector<int> travel(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> travel_naive(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> travel_nn(vector<pair<double, double>> cities, Grid<int> &distances);

    static vector<int> travel_cw(const vector<pair<double, double>>& cities, Grid<int> &distances);

    static vector<int> travel_cw_seq(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> local_2opt(vector<int> tour, Grid<int> &distances, Grid<int> &knearest,
                                  time_point<system_clock, duration<long, ratio<1, 1000000000>>> *deadline);

    static int tour_distance(const vector<pair<double, double>> &cities, vector<int> tour);

    static vector<pair<double, double>> create_n_cities(int n);
};


#endif //KTH_AA_TSP_H
