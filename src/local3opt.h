#ifndef KTH_AA_LOCAL3OPT_H
#define KTH_AA_LOCAL3OPT_H

#include <chrono>
#include <vector>
#include "utility.h"

using namespace std;
using namespace chrono;


template<class T, class U>
vector<T> local_3opt(vector<T> tour, Grid<U> &distances, Grid<T> &knn,
                     time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline,
                     bool reverse_seq, bool apply_seq) {
    bool better;
    int n = tour.size();
    do {
        better = false;
        int best_diff = 0, best_i, best_j, best_k;

        vector<T> city_to_tour_idx(n);
        for (T tour_idx = 0; tour_idx < n; tour_idx++) {
            city_to_tour_idx[tour[tour_idx]] = tour_idx;
        }

        for (T i = 0; i < n; i++) {
            if (deadline != nullptr and not(*deadline > steady_clock::now()))
                break;
            for (T j_knn = 0; j_knn < (int) knn.columns() - 1; j_knn++) {
                T j = city_to_tour_idx[knn[i][j_knn]];
                if (i >= j - 1) continue;

                for (T k_knn = j_knn + 1; k_knn < (int) knn.columns(); k_knn++) {
                    int k = city_to_tour_idx[knn[i][k_knn]];
                    if (j >= k - 1) continue;

                    for (T successor = 0; successor <= 1; successor++) {
                        T i_current = (i + successor) % n;
                        T j_current = (j + successor) % n;
                        T k_current = (k + successor) % n;
                        if (i_current >= j_current - 1 or j_current >= k_current - 1) continue;

                        int diff = reverse_seq
                                   ? reverse_segment_3opt_seq(&tour, i_current, j_current, k_current, distances,
                                                              apply_seq)
                                   : reverse_segment_3opt(&tour, i_current, j_current, k_current, distances, apply_seq);
                        if (diff > best_diff) {
                            if (not apply_seq) {
                                best_i = i_current, best_j = j_current, best_k = k_current;
                            }
                            better = true;
                        }
                    }
                }
            }
        }

        if (better) {
            if (not apply_seq) {
                reverse_seq
                ? reverse_segment_3opt_seq(&tour, best_i, best_j, best_k, distances, true)
                : reverse_segment_3opt(&tour, best_i, best_j, best_k, distances, true);
            }
        }

    } while (better and (deadline == nullptr or (*deadline > steady_clock::now())));

    return tour;
}

#endif //KTH_AA_LOCAL3OPT_H
