#ifndef KTH_AA_TSP_H
#define KTH_AA_TSP_H

#include <iostream>
#include <vector>

using namespace std;

class TSP {
public:
    static vector<int> travel(vector<pair<double, double>> cities);

    static vector<int> travel_naive(vector<pair<double, double>> cities);

    static vector<int> travel_cw(vector<pair<double, double>> cities);

    static int tour_distance(vector<pair<double, double>> cities, vector<int> tour);

    static vector<pair<double, double>> create_n_cities(int n);
};


#endif //KTH_AA_TSP_H
