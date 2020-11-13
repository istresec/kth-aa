#ifndef KTH_AA_LOCAL2OPT_H
#define KTH_AA_LOCAL2OPT_H

#include <chrono>
#include <vector>
#include "utility.h"

using namespace std;
using namespace chrono;

template<class T, class U>
inline U local2OptChange(Grid<T> &tour, Grid<U> &distances, T i, T j) {
    return (distances[i][tour[i][1]] + distances[j][tour[j][1]]) -
           (distances[i][j] + distances[tour[i][1]][tour[j][1]]);
}

template<class T>
inline void local2OptUpdate(Grid<T> &tour, T i, T j) {
    T ii = tour[i][1], jj = tour[j][1];
    tour[i][1] = j;
    tour[jj][0] = ii;
    T next = tour[ii][1];
    while (next != j) {
        swap(tour[next][0], tour[next][1]);
        next = tour[next][0];
    }
    tour[j][1] = tour[j][0];
    tour[j][0] = i;
    tour[ii][0] = tour[ii][1];
    tour[ii][1] = jj;
}

template<class T, class U>
vector<T> local2Opt(vector<T> tourVector, Grid<U> &distances, Grid<T> &knn,
                    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    // TODO use a tree structure for storing tour so that reverse is fast
    bool sequential = true;

    T n = tourVector.size();
    auto tour = Grid<T>(n, 2);
    for (T i = 0; i < n; i++) {
        tour[tourVector[i]][0] = tourVector[i == 0 ? n - 1 : i - 1];
        tour[tourVector[i]][1] = tourVector[i == n - 1 ? 0 : i + 1];
    }

    U bestChange;
    auto shouldBeChecked = vector<bool>(n, true);
    do {
        bestChange = 0.;
        T best_i, best_j;
        for (uint16_t i = 0; i < n - 1; i++) {

            for (T j = 0; j < knn.columns(); j++) {
                T a = i, b = knn[i][j];
                if (i > b) continue;
                if (tour[i][0] == knn[i][j] or tour[i][1] == knn[i][j]) continue;

                U change1 = local2OptChange(tour, distances, a, b);
                U change2 = local2OptChange(tour, distances, tour[a][0], tour[b][0]);
                U change = change1 > change2 ? change1 : change2;
                T betterA = change1 > change2 ? a : tour[a][0];
                T betterB = change1 > change2 ? b : tour[b][0];

                if (sequential and change > 0) {
                    local2OptUpdate(tour, betterA, betterB);
                }
                if (change > bestChange) {
                    bestChange = change;
                    best_i = betterA;
                    best_j = betterB;
                }
            }
        }
        if (bestChange > 0) {
            if (not sequential) {
                local2OptUpdate(tour, best_i, best_j);
            }
        }
    } while ((bestChange > 0) & (deadline == nullptr or (steady_clock::now() < *deadline)));

    T current = 0;
    for (unsigned i = 0; i < n; i++) {
        tourVector[i] = current;
        current = tour[current][1];
    }
    return tourVector;
}

#endif //KTH_AA_LOCAL2OPT_H
