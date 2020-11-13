#ifndef KTH_AA_GREEDY_H
#define KTH_AA_GREEDY_H

#include <vector>
#include <tuple>
#include <chrono>
#include "utility.h"

using namespace std;

template<class T, class U>
vector<T> travelGreedy(const vector<pair<double, double>> &cities, Grid<U> &distances) {
    vector<T> tour = vector<T>(cities.size());
    vector<bool> used = vector<bool>(cities.size());
    tour[0] = 0;
    used[0] = true;
    for (unsigned i = 1; i <= cities.size() - 1; i++) {
        int best = -1;
        for (unsigned j = 0; j <= cities.size() - 1; j++) {
            bool isBetter = best == -1 or distances[i - 1][j] < distances[i - 1][best];
            if (not used[j] and isBetter) {
                best = (int) j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
    return tour;
}

#endif //KTH_AA_GREEDY_H
