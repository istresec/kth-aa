#ifndef KTH_AA_TSP_H
#define KTH_AA_TSP_H

#include "utility.h"
#include "christofides.h"
#include "chokolino.h"
#include "bruteforce.h"
#include "local2opt.h"
#include "local3opt_no_knn_sequential.h"

using namespace std;
using namespace chrono;

template<class T, class U>
vector<T> travel(const vector<pair<double, double>> &cities) {
    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(1950);
    auto distances = distance_matrix<U>(cities);
    vector<T> tour;
    if (cities.size() <= 2) {
        tour = vector<T>();
        for (T i = 0; i < cities.size(); i++) tour.emplace_back(i);
    } else if (cities.size() <= 13) {
        return travel_bruteforce<T, U>(cities, *distances);
    } else {
        tour = travel_christofides<T, U>(*distances);
        Grid<T> *knn = k_nearest_neighbors<T>(*distances, min((int) (cities.size() - 1), 300));
        tour = chokolino<T, U>(tour, *distances, *knn, &deadline);
        tour = local_2opt<T, U>(tour, *distances, *knn, &deadline);
        tour = local_3opt_no_knn_sequential<T, U>(tour, *distances, *knn, &deadline);
    }
    return tour;
}


#endif //KTH_AA_TSP_H
