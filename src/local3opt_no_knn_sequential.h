#ifndef KTH_AA_LOCAL3OPT_NO_KNN_SEQUENTIAL_H
#define KTH_AA_LOCAL3OPT_NO_KNN_SEQUENTIAL_H

#include <chrono>
#include <vector>
#include "utility.h"

using namespace std;
using namespace chrono;

// 3-opt impl. based on pseudocode from https://en.wikipedia.org/wiki/3-opt
template<class T, class U>
vector<T> local_3opt_no_knn_sequential(vector<T> tour, Grid<U> &distances, Grid<T> &,
                                       time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    bool better;
    do {
        T n = tour.size();
        better = false;
        for (T i = 0; i < n; i++) {
            if (deadline != nullptr and not(*deadline > steady_clock::now()))
                break;
            for (T j = i + 2; j < n; j++) {
                for (T k = j + 2; k < n + (T) (i > 0); k++) {
                    if (reverse_segment_3opt_seq(&tour, i, j, k, distances, true) > 0)
                        better = true;
                }
            }
        }
    } while (better and (deadline == nullptr or (*deadline > steady_clock::now())));

    return tour;
}

#endif //KTH_AA_LOCAL3OPT_NO_KNN_SEQUENTIAL_H
