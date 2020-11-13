#ifndef KTH_AA_LOCAL2OPT_NO_KNN_H
#define KTH_AA_LOCAL2OPT_NO_KNN_H

#include <vector>
#include <chrono>
#include "utility.h"

using namespace std;
using namespace chrono;

template<class T, class U>
vector<T> local2OptNoKnn(vector<T> tour, Grid<U> &distances, Grid<T> &,
                         time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    int bestChange;
    do {
        bestChange = 0.;
        unsigned bestI, bestJ;
        int change;
        for (T i = 0; i < tour.size() - 2; i++) {
            for (T j = i + 2; j < tour.size(); j++) {
                change = (distances[tour[i]][tour[i + 1]] + distances[tour[j]][tour[(j + 1) % tour.size()]]) -
                         (distances[tour[i]][tour[j]] + distances[tour[i + 1]][tour[(j + 1) % tour.size()]]);
                if (change > bestChange) {
                    bestChange = change;
                    bestI = i;
                    bestJ = j;
                }
            }
        }
        if (bestChange > 0) {
            reverse(tour.begin() + bestI + 1, tour.begin() + bestJ + 1);
        }
    } while ((bestChange > 0) and
             (deadline == nullptr or (*deadline > steady_clock::now())));

    return tour;
}


#endif //KTH_AA_LOCAL2OPT_NO_KNN_H
