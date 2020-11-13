#ifndef KTH_AA_LOCAL3OPT_H
#define KTH_AA_LOCAL3OPT_H

#include <chrono>
#include <vector>
#include "utility.h"

using namespace std;
using namespace chrono;

template<class T, class U>
inline int reverseSegment3Opt(vector<T> *tour, int i, int j, int k, Grid<U> &distances, bool apply) {
    if (i < 0 or i >= j - 1 or j >= k - 1 or k > tour->size()) { // TODO remove
        throw logic_error("Inconsistent state! Me no likey.");
    }

    int a = (*tour)[(i - 1 + tour->size()) % tour->size()];
    int b = (*tour)[i];
    int c = (*tour)[j - 1];
    int d = (*tour)[j];
    int e = (*tour)[k - 1];
    int f = (*tour)[k % (*tour).size()];

    // original distance
    int d0 = distances[a][b] + distances[c][d] + distances[e][f];
    int d1 = distances[a][c] + distances[b][d] + distances[e][f];
    int d2 = distances[a][b] + distances[c][e] + distances[d][f];
    int d3 = distances[a][d] + distances[e][b] + distances[c][f];
    int d4 = distances[f][b] + distances[c][d] + distances[e][a];
    int d5 = distances[a][e] + distances[d][b] + distances[c][f];
    int d6 = distances[a][c] + distances[b][e] + distances[d][f];
    int d7 = distances[a][d] + distances[e][c] + distances[b][f];

    int best = min(min(min(min(min(min(d1, d2), d3), d4), d5), d6), d7);

    if (best >= d0) return 0;
    if (not apply) return d0 - best;

    if (best == d1) {
        reverse(tour->begin() + i, tour->begin() + j);
    } else if (best == d2) {
        reverse(tour->begin() + j, tour->begin() + k);
        reverse(tour->begin() + j, tour->begin() + k);
    } else if (best == d3) {
        vector<int> tempTour = vector<int>{};
        tempTour.insert(tempTour.end(), tour->begin() + j, tour->begin() + k);
        tempTour.insert(tempTour.end(), tour->begin() + i, tour->begin() + j);
        copy_n(tempTour.begin(), tempTour.size(), &(*tour)[i]);
    } else if (best == d4) {
        reverse(tour->begin() + i, tour->begin() + k);
    } else if (best == d5) {
        // get ac bd ef like in d1
        reverse(tour->begin() + i, tour->begin() + j);
        // get ae db cf
        reverse(tour->begin() + i, tour->begin() + k);
    } else if (best == d6) {
        // get ac bd ef like with d1
        reverse(tour->begin() + i, tour->begin() + j);
        // now reverse from d to e
        reverse(tour->begin() + j, tour->begin() + k);
    } else if (best == d7) {
        // get ab ce df like with d2
        reverse(tour->begin() + j, tour->begin() + k);
        // now reverse from d to b
        reverse(tour->begin() + i, tour->begin() + k);
    } else {
        throw runtime_error("inconsistent state");
    }


    return d0 - best;
}

template<class T, class U>
vector<T> local3Opt(vector<T> tour, Grid<U> &distances, Grid<T> &knn,
                    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline,
                    bool reverseSeq, bool applySeq) {
    bool better;
    int n = tour.size();
    do {
        better = false;
        int bestDiff = 0, bestI, bestJ, bestK;

        vector<T> cityToTourIdx(n);
        for (T tourIdx = 0; tourIdx < n; tourIdx++) {
            cityToTourIdx[tour[tourIdx]] = tourIdx;
        }

        for (T i = 0; i < n; i++) {
            if (deadline != nullptr and not(*deadline > steady_clock::now()))
                break;
            for (T jKnn = 0; jKnn < (int) knn.columns() - 1; jKnn++) {
                T j = cityToTourIdx[knn[i][jKnn]];
                if (i >= j - 1) continue;

                for (T kKnn = jKnn + 1; kKnn < (int) knn.columns(); kKnn++) {
                    int k = cityToTourIdx[knn[i][kKnn]];
                    if (j >= k - 1) continue;

                    for (T successor = 0; successor <= 1; successor++) {
                        T iCurrent = (i + successor) % n;
                        T jCurrent = (j + successor) % n;
                        T kCurrent = (k + successor) % n;
                        if (iCurrent >= jCurrent - 1 or jCurrent >= kCurrent - 1) continue;

                        int diff = reverseSeq
                                   ? reverseSegment3OptSeq(&tour, iCurrent, jCurrent, kCurrent, distances,
                                                           applySeq)
                                   : reverseSegment3Opt(&tour, iCurrent, jCurrent, kCurrent, distances, applySeq);
                        if (diff > bestDiff) {
                            if (not applySeq) {
                                bestI = iCurrent, bestJ = jCurrent, bestK = kCurrent;
                            }
                            better = true;
                        }
                    }
                }
            }
        }

        if (better) {
            if (not applySeq) {
                reverseSeq
                ? reverseSegment3OptSeq(&tour, bestI, bestJ, bestK, distances, true)
                : reverseSegment3Opt(&tour, bestI, bestJ, bestK, distances, true);
            }
        }

    } while (better and (deadline == nullptr or (*deadline > steady_clock::now())));

    return tour;
}

#endif //KTH_AA_LOCAL3OPT_H
