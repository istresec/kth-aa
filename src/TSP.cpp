#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <random>
#include <chrono>

#include "TSP.h"

using namespace std;
using namespace chrono;

vector<int> TSP::travel(const vector<pair<double, double>> &cities) {
    time_point<system_clock, duration<long, ratio<1, 1000000000>>> deadline = system_clock::now() +
                                                                              milliseconds(1950);
    auto distances = distance_matrix(cities);
    auto tour = travel_cw(cities, *distances);
    Grid<uint16_t> *knn = k_nearest_neighbors<uint16_t>(*distances, 200);
    tour = TSP::local_2opt(tour, *distances, *knn, &deadline);
//    Grid<uint16_t> *knn = nullptr;
//    tour = TSP::local_2opt_no_knn(tour, *distances, *knn, &deadline);
    return tour;
}

int TSP::tour_distance(const vector<pair<double, double>> &cities, vector<int> tour) {
    int distance = 0;
    auto tsp_table = distance_matrix(cities);
    for (unsigned i = 0; i < tour.size(); i++) {
        distance += (*tsp_table)[tour[i]][tour[(i + 1) % cities.size()]];
    }
    return distance;
}

vector<pair<double, double>> TSP::create_n_cities(int n, int seed) {
    double lower_bound = -10e6;
    double upper_bound = 10e6;
    vector<pair<double, double >> cities = vector<pair<double, double >>();
    cities.reserve(n);

    uniform_real_distribution<double> unif(lower_bound, upper_bound);
    default_random_engine generator(seed);
    for (int i = 0; i < n; i++) {
        cities.emplace_back(pair<double, double>(unif(generator), unif(generator)));
    }
    return cities;
}

inline bool compare_savings(tuple<int, int, int> t1, tuple<int, int, int> t2) {
    return get<2>(t1) > get<2>(t2);
}

vector<tuple<int, int, int>> savings(const vector<pair<double, double>> &cities, Grid<int> &distances) {
    vector<tuple<int, int, int >> savings = vector<tuple<int, int, int >>();
    for (unsigned i = 1; i < cities.size() - 1; i++) {
        for (unsigned j = i + 1; j < cities.size(); j++) {
            savings.emplace_back(make_tuple(i, j, distances[0][i] + distances[0][j] - distances[i][j]));
        }
    }
    sort(savings.begin(), savings.end(), compare_savings);
    return savings;
}

vector<int> TSP::travel_naive(const vector<pair<double, double>> &cities, Grid<int> &distances) {
    vector<int> tour = vector<int>(cities.size());
    vector<bool> used = vector<bool>(cities.size());
    tour[0] = 0;
    used[0] = true;
    for (unsigned i = 1; i <= cities.size() - 1; i++) {
        int best = -1;
        for (unsigned j = 0; j <= cities.size() - 1; j++) {
            bool is_better = best == -1 or distances[i - 1][j] < distances[i - 1][best];
            if (not used[j] and is_better) {
                best = (int) j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
    return tour;
}

// Nearest Neighbour Algorithm - pick random vertex and go the shortest distance you can
vector<int> TSP::travel_nn(vector<pair<double, double>> cities, Grid<int> &distances) {
    vector<int> tour = vector<int>();
    return tour;
}

// Clarke-Wright Savings Algorithm. City at index 0 used as hub.
vector<int> TSP::travel_cw(const vector<pair<double, double>> &cities, Grid<int> &distances) {
    // if only one city
    if (cities.size() == 1)
        return vector<int>{0};

    // get savings
    vector<tuple<int, int, int>> s = savings(cities, distances);

    // initialize tours
    vector<int> tours[cities.size()];
    for (int i = 1; i < (int) cities.size(); i++)
        tours[i] = vector<int>{i}; // instead of 0, i, 0 just use i

    // algorithm
    vector<int> temp_tour;
    vector<int> first, second;
    vector<int> empty_tour = vector<int>();
    int i, j;
    for (auto &it : s) {
        i = get<0>(it);
        j = get<1>(it);
        // if two distinct tours with endpoints i and j exist, if so combine them (O(1))
        if (not tours[i].empty() and not tours[j].empty() and tours[i].front() != tours[j].front()) {
            first = tours[i]; // remember tour with endpoint i
            second = tours[j]; // remember tour with endpoint j
            // remove tours with endpoints i and j while a new tour is constructed
            tours[first.front()] = tours[first.back()] = tours[second.front()] = tours[second.back()] = empty_tour;

            if (first.front() == i)
                reverse(first.begin(), first.end()); // reverse tour with endpoint i if it starts with i
            if (second.front() != j)
                reverse(second.begin(), second.end()); // reverse tour with j if it doesn't start with j

            // create new tour by joining the two
            first.insert(first.end(), second.begin(), second.end());

            // remember endpoints of the new tour in array of endpoints for quick access
            tours[first.front()] = first;
            tours[first.back()] = first;
        }
    }

    // create final tour
    vector<int> tour = vector<int>();
    for (i = 1; i < (int) cities.size(); i++) {
        if (not tours[i].empty())
            break;
    }
    // check where to place the hub to minimize tour distance (start/end?)
    if (distances[0][tours[i].front()] < distances[0][tours[i].back()]) {
        tour.emplace_back(0);
        tour.insert(tour.end(), tours[i].begin(), tours[i].end());
    } else {
        tour = tours[i];
        tour.emplace_back(0);
    }

    return tour;
}

// TODO: implement if needed
// Clarke-Wright Sequential Savings Algorithm. City at index 0 used as hub.
vector<int> TSP::travel_cw_seq(const vector<pair<double, double>> &cities, Grid<int> &distances) {
    // if only one city
    if (cities.size() == 1)
        return vector<int>{0};

    // get savings
    vector<tuple<int, int, int>> s = savings(cities, distances);

    // initialize tours
    vector<int> tours[cities.size()];
    for (int i = 1; i < (int) cities.size(); i++)
        tours[i] = vector<int>{i}; // instead of 0, i, 0 just use i

    // algorithm
    vector<int> temp_tour;
    vector<int> first, second;
    vector<int> empty_tour = vector<int>();
    int i, j;
    for (auto &it : s) {
        i = get<0>(it);
        j = get<1>(it);
        // if two distinct tours with endpoints i and j exist, if so combine them (O(1))
        if (not tours[i].empty() and not tours[j].empty() and tours[i].front() != tours[j].front()) {
            first = tours[i]; // remember tour with endpoint i
            second = tours[j]; // remember tour with endpoint j
            // remove tours with endpoints i and j while a new tour is constructed
            tours[first.front()] = tours[first.back()] = tours[second.front()] = tours[second.back()] = empty_tour;

            if (first.front() == i)
                reverse(first.begin(), first.end()); // reverse tour with endpoint i if it starts with i
            if (second.front() != j)
                reverse(second.begin(), second.end()); // reverse tour with j if it doesn't start with j

            // create new tour by joining the two
            first.insert(first.end(), second.begin(), second.end());

            // remember endpoints of the new tour in array of endpoints for quick access
            tours[first.front()] = first;
            tours[first.back()] = first;
        }
    }

    // create final tour
    vector<int> tour = vector<int>();
    for (i = 1; i < (int) cities.size(); i++) {
        if (not tours[i].empty())
            break;
    }
    // check where to place the hub to minimize tour distance (start/end?)
    if (distances[0][tours[i].front()] < distances[0][tours[i].back()]) {
        tour.emplace_back(0);
        tour.insert(tour.end(), tours[i].begin(), tours[i].end());
    } else {
        tour = tours[i];
        tour.emplace_back(0);
    }

    return tour;
}

template<class T, class U>
inline int local2OptChange(Grid<T> &tour, Grid<U> &distances, int i, int j) {
    return (distances[i][tour[i][1]] + distances[j][tour[j][1]]) -
           (distances[i][j] + distances[tour[i][1]][tour[j][1]]);
}

inline void local2OptUpdate(Grid<uint16_t> &tour, uint16_t i, uint16_t j) {
    uint16_t ii = tour[i][1], jj = tour[j][1];
    tour[i][1] = j;
    tour[jj][0] = ii;
    uint16_t next = tour[ii][1];
    while (next != j) {
        swap(tour[next][0], tour[next][1]);
        next = tour[next][0];
    }
    tour[j][1] = tour[j][0];
    tour[j][0] = i;
    tour[ii][0] = tour[ii][1];
    tour[ii][1] = jj;
}

vector<int> TSP::local_2opt(vector<int> tour_vector, Grid<int> &distances, Grid<uint16_t> &knn,
                            time_point<system_clock, duration<long, ratio<1, 1000000000>>> *deadline) {
// TODO use a tree structure for storing tour so that reverse is fast
    uint16_t n = tour_vector.size();
    auto tour = Grid<uint16_t>(n, 2);
    for (int i = 0; i < n; i++) {
        tour[tour_vector[i]][0] = tour_vector[i == 0 ? i - 1 + n : i - 1];
        tour[tour_vector[i]][1] = tour_vector[i == n - 1 ? i + 1 - n : i + 1];
    }

    int bestChange;
    auto shouldBeChecked = vector<bool>(n, true);
    do {
        bestChange = 0.;
        uint16_t best_i, best_j;
        for (uint16_t i = 0; i < n - 1; i++) {
//            if (not shouldBeChecked[i]) continue;

//            bool iCannotDoBetter = true;
            for (uint16_t j = 0; j < knn.columns(); j++) {
                uint16_t a = i, b = knn[i][j];
                if (i > b) continue;
                if (tour[i][0] == knn[i][j] or tour[i][1] == knn[i][j]) continue;

                int change1 = local2OptChange(tour, distances, a, b);
                int change2 = local2OptChange(tour, distances, tour[a][0], tour[b][0]);
                int change = change1 > change2 ? change1 : change2;
                uint16_t better_a = change1 > change2 ? a : tour[a][0];
                uint16_t better_b = change1 > change2 ? b : tour[b][0];

//                if (change > 0) { iCannotDoBetter = false; }
                if (change > bestChange) {
                    bestChange = change;
                    best_i = better_a;
                    best_j = better_b;
                }
            }
//            if (iCannotDoBetter) { shouldBeChecked[i] = false; }
        }
        if (bestChange > 0) {
//            for (uint16_t j = 0; j < knn.columns(); j++) {
//                shouldBeChecked[knn[best_i][j]] = true;
//                shouldBeChecked[knn[best_j][j]] = true;
//                shouldBeChecked[knn[tour[best_i][1]][j]] = true;
//                shouldBeChecked[knn[tour[best_j][1]][j]] = true;
//            }
            local2OptUpdate(tour, best_i, best_j);
        }
    } while ((bestChange > 0) & (deadline == nullptr or (system_clock::now() < *deadline)));

    uint16_t current = 0;
    for (unsigned i = 0; i < n; i++) {
        tour_vector[i] = current;
        current = tour[current][1];
    }
    return tour_vector;
}

vector<int> TSP::local_2opt_no_knn(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                                   time_point<system_clock, duration<long, ratio<1, 1000000000>>> *deadline) {
    int best_change;
    do {
        best_change = 0.;
        unsigned best_i, best_j;
        int change;
        for (unsigned i = 0; i < tour.size() - 2; i++) {
            for (unsigned j = i + 2; j < tour.size(); j++) {
                change = (distances[tour[i]][tour[i + 1]] + distances[tour[j]][tour[(j + 1) % tour.size()]]) -
                         (distances[tour[i]][tour[j]] + distances[tour[i + 1]][tour[(j + 1) % tour.size()]]);
                if (change > best_change) {
                    best_change = change;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        if (best_change > 0) {
            reverse(tour.begin() + best_i + 1, tour.begin() + best_j + 1);
        }
    } while ((best_change > 0) & (deadline == nullptr or (system_clock::now() < *deadline)));

    return tour;
}
