#ifndef KTH_AA_CHOKOLINO
#define KTH_AA_CHOKOLINO

#include <chrono>
#include <vector>
#include <tuple>
#include <stack>
#include <unordered_set>
#include "utility.h"

using namespace std;
using namespace chrono;

// insures that edge defined by two towns is unique for sets (the first in the pair is always smaller)
template<class T>
pair<T, T> createEdgePair(T t1, T t2) {
    if (t1 < t2) return pair(t1, t2);
    return pair(t2, t1);
}

template<class T, class U>
bool chooseY(vector<T> &tour, unordered_set<pair<T, T>, pair_hash> &tourEdges, vector<T> &cityToTourIdx,
             Grid<U> &distances, Grid<T> &knn, T tour1, T tour2I, U gain,
             unordered_set<pair<T, T>, pair_hash> &broken, unordered_set<pair<T, T>, pair_hash> &joined,
             time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

template<class T>
inline unordered_set<pair<T, T>, pair_hash> getTourEdges(const vector<T> &tour) {
    unordered_set<pair<T, T>, pair_hash> tourEdges;
    tourEdges.reserve(tour.size());
    for (T i = 0; i < tour.size(); i++) {
        tourEdges.emplace(createEdgePair(tour[i], tour[(i + 1) % tour.size()]));
    }
    return tourEdges;
}

template<class T>
inline bool generateTour(const vector<T> &tour,
                         const unordered_set<pair<T, T>, pair_hash> &added,
                         const unordered_set<pair<T, T>, pair_hash> &removed,
                         vector<T> &newTour) {
    T n = tour.size();
    auto cityToNeighbors = Grid<T>(n, 2); // city --> [pred, successor] initially
    for (T i = 0; i < n; i++) {
        cityToNeighbors[tour[i]][0] = tour[i == 0 ? n - 1 : i - 1];
        cityToNeighbors[tour[i]][1] = tour[i == n - 1 ? 0 : i + 1];
    }

    for (auto edgeToRemove: removed) {
        if (cityToNeighbors[edgeToRemove.first][0] == edgeToRemove.second) {
            cityToNeighbors[edgeToRemove.first][0] = n;
            cityToNeighbors[edgeToRemove.second][1] = n;
        } else {
            cityToNeighbors[edgeToRemove.first][1] = n;
            cityToNeighbors[edgeToRemove.second][0] = n;
        }
    }

    for (auto edgeToAdd: added) {
        if (cityToNeighbors[edgeToAdd.first][0] == n) {
            cityToNeighbors[edgeToAdd.first][0] = edgeToAdd.second;
        } else {
            cityToNeighbors[edgeToAdd.first][1] = edgeToAdd.second;
        }
        if (cityToNeighbors[edgeToAdd.second][0] == n) {
            cityToNeighbors[edgeToAdd.second][0] = edgeToAdd.first;
        } else {
            cityToNeighbors[edgeToAdd.second][1] = edgeToAdd.first;
        }
    }

    newTour.reserve(n);
    T startCity = 0;
    T prevCity;
    T nextCity = startCity;
    newTour.emplace_back(nextCity);
    for (T i = 0; i < n - 1; i++) {
        auto left = cityToNeighbors[nextCity][0];
        auto right = cityToNeighbors[nextCity][1];
        if (right == prevCity) {
            prevCity = nextCity;
            nextCity = left;
            newTour.emplace_back(nextCity);
        } else {
            prevCity = nextCity;
            nextCity = right;
            newTour.emplace_back(nextCity);
        }

        if (nextCity == startCity) {
            return false;
        }
    }

    if (cityToNeighbors[nextCity][0] != startCity and cityToNeighbors[nextCity][1] != startCity) {
        return false;
    }

    return true;
}

// chooses next edge to remove, part of the Lin-Kernighan algorithm implementation
template<class T, class U>
bool chooseX(vector<T> &tour, unordered_set<pair<T, T>, pair_hash> &tourEdges, vector<T> &cityToTourIdx,
             Grid<U> &distances, Grid<T> &knn, T tour1, T tourLast, U gain,
             unordered_set<pair<T, T>, pair_hash> &broken, unordered_set<pair<T, T>, pair_hash> &joined,
             time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {

    // neighbourhood of last edge in our tour (indices in tour)
    pair<T, T> neighbourhood = make_pair((T) (tourLast + 1) % tour.size(),
                                         (T) ((tourLast + tour.size() - 1) % tour.size()));

    // special case for x4 from Helsgaun
    if (broken.size() == 4) {
        // give choice priority to longer edge
        if (distances[tour[neighbourhood.first]][tourLast] > distances[tour[neighbourhood.first]][tourLast]) {
            swap(neighbourhood.first, neighbourhood.second);
        }
    }

    T city1 = tour[tour1];
    for (uint8_t nIdx = 0; nIdx <= 1; nIdx++) {
        T tour2I = nIdx == 0 ? neighbourhood.first : neighbourhood.second;
        T city2I = tour[tour2I];
        auto xi = createEdgePair(tour[tourLast], city2I);

        U gainI = gain + distances[tour[tourLast]][city2I];

        // verify that xi is not already in X or Y
        if (joined.find(xi) == joined.end() and broken.find(xi) == broken.end()) {
            auto removed = broken;
            removed.emplace(xi);

            auto yRelink = createEdgePair(city2I, city1);
            vector<T> newTour;
            bool isTour;
            U relink = gainI - distances[city2I][city1];
            if (relink > 0 and joined.find(yRelink) == joined.end() and broken.find(yRelink) == broken.end()) {
                auto added = joined;
                added.emplace(yRelink);
                isTour = generateTour(tour, added, removed, newTour);
            } else {
                isTour = false;
            }

            if (isTour) {
                tour = newTour;
                tourEdges = getTourEdges(tour);
                createCityToTourIdx(tour, cityToTourIdx);
                return true;
            } else {
                if (joined.size() > 3) continue;
                return chooseY(tour, tourEdges, cityToTourIdx, distances, knn, tour1, tour2I, gainI, removed,
                               joined, deadline);
            }
        }
    }

    return false;
}

// chooses next edge to add, part of the Lin-Kernighan algorithm implementation
template<class T, class U>
bool chooseY(vector<T> &tour, unordered_set<pair<T, T>, pair_hash> &tourEdges, vector<T> &cityToTourIdx,
             Grid<U> &distances, Grid<T> &knn, T tour1, T tour2I, U gain,
             unordered_set<pair<T, T>, pair_hash> &broken, unordered_set<pair<T, T>, pair_hash> &joined,
             time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {

    T city2I = tour[tour2I];

    // original LKH heuristic: check only 5 closest neighbours when doing y2, otherwise just 1
    // int top = broken.size() == 2 ? 5 : 1;
    // int top = (broken.size() == 2) ? 30 : ((broken.size() == 3) ? 5 : 5);
    int top = (broken.size() == 2) ? 20 : ((broken.size() == 3) ? 5 : ((broken.size() == 4) ? 3 : 1)); // tuning->kattis


    // try to find an edge to add from closest neighbours
    for (unsigned i = 0; i < knn.columns() and top >= 1; i++) {
        T city2IPlus1 = knn[city2I][i];
        auto yi = createEdgePair(city2I, city2IPlus1);

        // Three checks need to be done
        // 1. Gain still positive
        U gainI = gain - distances[city2I][city2IPlus1];
        if (gainI <= 0) {
            continue;
        }
        // 2. yi not already removed or added to tour
        if (joined.find(yi) != joined.end() or broken.find(yi) != broken.end()) {
            continue;
        }
        // 3. yi is not in tour_edges
        if (tourEdges.find(yi) != tourEdges.end()) {
            continue;
        }

        top--;
        auto added = joined;
        added.emplace(yi);
        if (chooseX(tour, tourEdges, cityToTourIdx, distances, knn, tour1, cityToTourIdx[city2IPlus1],
                    gainI, broken, added, deadline)) {
            return true;
        }
    }

    return false;
}

// main loop for the Lin-Kernighan algorithm
template<class T, class U>
bool linKernighanMainLoop(vector<T> &tour, unordered_set<pair<T, T>, pair_hash> &tourEdges,
                          vector<T> &cityToTourIdx, Grid<U> &distances, Grid<T> &knn,
                          time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    // steps refer to a basic outline of the LK algorithm found on page 12 of
    // "An Effective Implementation of the Lin-Kernighan Traveling Salesman Heuristic"
    // by Keld Helsgaun (Roskilde University, Denmark)
    // http://akira.ruc.dk/~keld/research/LKH/LKH-2.0/DOC/LKH_REPORT.pdf

    auto setX = unordered_set<pair<T, T>, pair_hash>(); // set of broken edges
    auto setY = unordered_set<pair<T, T>, pair_hash>(); // set of added edges

    bool anyImprovement = false;
    // step 2
    for (T tourI = 0; tourI < (T) tour.size(); tourI++) { // go through all possible towns for t1
        T city1 = tour[tourI];

        if (deadline != nullptr and (*deadline <= steady_clock::now())) {
            return false;
        }

        bool innerImprovement = false;
        // step 3
        for (uint8_t nIdx = 0; nIdx <= 1; nIdx++) { // get towns around t1 for t2
            T tourJ = (tourI + tour.size() + (nIdx == 0 ? -1 : 1)) % tour.size();
            if (innerImprovement) break;
            T city2 = tour[tourJ];

            // create edge x1 (edge between t1 and t2) and add it to set of edges X
            auto x1 = createEdgePair(city1, city2);
            setX = unordered_set<pair<T, T>, pair_hash>(); // or set_x.clear();
            setX.emplace(x1);

            // step 4
            // find t3 to create a second edge for 2-opt
            for (T neighborIdx = 0; neighborIdx < (T) knn.columns(); neighborIdx++) {
                // TODO consider using some better metric than knn
                //      (like sorting neighbors by 1-tree score or by possible relinking gain)
                auto city3 = knn[city2][neighborIdx];
                auto y1 = createEdgePair(city2, city3);
                U gain = distances[city1][city2] - distances[city2][city3];

                if (gain > 0) { // if the two edges would create a shorter tour go find more
                    // y1 is not in tour_edges, if it is continue
                    if (tourEdges.count(y1) > 0) {
                        continue;
                    }
                    // step 6
                    T tourK = cityToTourIdx[city3];
                    setY = unordered_set<pair<T, T>, pair_hash>(); // or set_y.clear()
                    setY.emplace(y1);
                    if (chooseX(tour, tourEdges, cityToTourIdx, distances, knn, tourI, tourK, gain, setX,
                                setY, deadline)) {
                        innerImprovement = true;
                        anyImprovement = true;
                        break;
                    }
                }
            }
        }
    }

    return anyImprovement;
}

// implementation based on outline from https://arthur.maheo.net/implementing-lin-kernighan-in-python/
// note: Chokolino is the codename for Lin Kernighan
template<class T, class U>
static vector<T> chokolino(vector<T> &tour, Grid<U> &distances, Grid<T> &knn,
                           time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    bool improved = true;
    unordered_set<pair<T, T>, pair_hash> tourEdges = getTourEdges(tour);
    vector<T> cityToTourIdx(tour.size());
    cityToTourIdx.reserve(tour.size());
    createCityToTourIdx(tour, cityToTourIdx);
    while ((deadline == nullptr or (*deadline > steady_clock::now())) and improved) {
        improved = linKernighanMainLoop(tour, tourEdges, cityToTourIdx, distances, knn, deadline);
    }
    return tour;
}

#endif // KTH_AA_CHOKOLINO