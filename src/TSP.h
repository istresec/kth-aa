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
    auto distances = distanceMatrix<U>(cities);
    vector<T> tour;
    if (cities.size() <= 2) {
        tour = vector<T>();
        for (T i = 0; i < cities.size(); i++) tour.emplace_back(i);
    } else if (cities.size() <= 13) {
        return travelBruteforce<T, U>(cities, *distances);
    } else {
        tour = travelChristofides<T, U>(*distances);
        Grid<T> *knn = kNearestNeighbors<T>(*distances, min((int) (cities.size() - 1), 300));
        tour = chokolino<T, U>(tour, *distances, *knn, &deadline);
        tour = local2Opt<T, U>(tour, *distances, *knn, &deadline);
        tour = local3OptNoKnnSequential<T, U>(tour, *distances, *knn, &deadline);
    }
    return tour;
}


#endif //KTH_AA_TSP_H
