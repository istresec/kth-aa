#ifndef KTH_AA_BRUTEFORCE_H
#define KTH_AA_BRUTEFORCE_H

#include <vector>
#include <tuple>
#include "utility.h"

using namespace std;

template<class T, class U>
vector<T> travel_bruteforce(const vector<pair<double, double>> &cities, Grid<U> &distances) {
    auto city_indices = vector<T>();
    city_indices.reserve(cities.size());
    for (T i = 0; i < cities.size(); i++) city_indices.emplace_back(i);

    auto best = city_indices;
    U best_distance = -1;
    do {
        auto current_distance = tour_distance(distances, city_indices);
        if (best_distance == -1 or current_distance < best_distance) {
            best_distance = current_distance;
            best = city_indices;
        }
    } while (next_permutation(city_indices.begin() + 1, city_indices.end()));

    return best;
}


#endif //KTH_AA_BRUTEFORCE_H
