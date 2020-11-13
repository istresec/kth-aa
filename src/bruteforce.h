#ifndef KTH_AA_BRUTEFORCE_H
#define KTH_AA_BRUTEFORCE_H

#include <vector>
#include <tuple>
#include "utility.h"

using namespace std;

template<class T, class U>
vector<T> travelBruteforce(const vector<pair<double, double>> &cities, Grid<U> &distances) {
    auto cityIndices = vector<T>();
    cityIndices.reserve(cities.size());
    for (T i = 0; i < cities.size(); i++) cityIndices.emplace_back(i);

    auto best = cityIndices;
    U bestDistance = -1;
    do {
        auto currentDistance = tourDistance(distances, cityIndices);
        if (bestDistance == -1 or currentDistance < bestDistance) {
            bestDistance = currentDistance;
            best = cityIndices;
        }
    } while (next_permutation(cityIndices.begin() + 1, cityIndices.end()));

    return best;
}


#endif //KTH_AA_BRUTEFORCE_H
