#include "TSP.h"

#include <iostream>
#include <vector>

using namespace std;

double squared_distance(pair<double, double> &t1, pair<double, double> &t2) {
    double dx = t1.first - t2.first;
    double dy = t1.second - t2.second;
    return dx * dx + dy * dy;
}

vector<int> TSP::travel(vector<pair<double, double>> cities) {
    return TSP::travel_naive(cities);
}

vector<int> TSP::travel_naive(vector<pair<double, double>> cities) {
    vector<int> tour = vector<int>(cities.size());
    vector<bool> used = vector<bool>(cities.size());
    tour[0] = 0;
    used[0] = true;
    for (int long i = 1; i <= cities.size() - 1; i++) {
        int best = -1;
        for (int j = 0; j <= cities.size() - 1; j++) {
            bool is_better = best == -1 or squared_distance(cities[tour[i - 1]], cities[j]) <
                                           squared_distance(cities[tour[i - 1]], cities[best]);
            if (not used[j] and is_better) {
                best = j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
    return tour;
}