#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <random>
#include <chrono>

#include "TSP.h"
#include "utility.h"

using namespace std;
using namespace chrono;

int TSP::tour_distance(const vector<pair<double, double>> &cities, vector<int> tour) {
    int distance = 0;
    auto tsp_table = distance_matrix(cities);
    for (unsigned i = 0; i < tour.size() - 1; i++) {
        distance += (*tsp_table)[tour[i]][tour[i + 1]];
    }
    return distance;
}

vector<pair<double, double>> TSP::create_n_cities(int n, int seed) {
    double lower_bound = 0;
    double upper_bound = (double) n * 10;
    vector<pair<double, double >> cities = vector<pair<double, double >>();
    cities.reserve(n);

    uniform_real_distribution<double> unif(lower_bound, upper_bound);
    chrono::high_resolution_clock clock;
    default_random_engine generator;
    generator.seed(seed);
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

vector<int> TSP::travel(const vector<pair<double, double>> &cities) {
    auto distances = distance_matrix(cities);
    auto tour = TSP::travel_cw(cities, *distances);
    return tour;
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
    // get savings
    vector<tuple<int, int, int>> s = savings(cities, distances);

    // initialize tours
    vector<int> tours[cities.size()];
    for (int i = 1; i < cities.size(); i++)
        tours[i] = vector<int>{i}; // instead of 0, i, 0 just use i

    // algorithm
    vector<int> temp_tour;
    vector<int> first, second;
    int i, j;
    for (auto &it : s) {
        i = get<0>(it);
        j = get<1>(it);
        // if two distinct tours with endpoints i and j exist, if so combine them (O(1))
        if (not tours[i].empty() and not tours[j].empty() and tours[i].front() != tours[j].front()) {
            first = tours[i]; // remember tour with endpoint i
            second = tours[j]; // remember tour with endpoint j
            // remove tours with endpoints i and j while a new tour is constructed
            tours[first.front()] = tours[first.back()] = tours[second.front()] = tours[second.back()] = vector<int>();

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
    tour.emplace_back(0);
    for (i = 1; i < cities.size(); i++) {
        if (not tours[i].empty())
            break;
    }
    tour.insert(tour.end(), tours[i].begin(), tours[i].end());

    cout << "Returning Clarke Wright.\n";
    return tour;
}

// Clarke-Wright Sequential Savings Algorithm. City at index 0 used as hub.
vector<int> TSP::travel_cw_seq(const vector<pair<double, double>> &cities, Grid<int> &distances) {
    // get savings
    vector<tuple<int, int, int>> s = savings(cities, distances);

    // initialize tours
    vector<vector<int>> tours = vector<vector<int >>();
    for (int i = 1; i < cities.size(); i++) {
        tours.emplace_back(vector<int>{i}); // instead of 0, i, 0 just use i
    }

    // TODO: implement if needed

    return vector<int>();
}

vector<int> TSP::local_2opt(vector<int> tour, Grid<int> &distances, Grid<int> &knearest,
                            time_point<system_clock, duration<long, ratio<1, 1000000000>>> *deadline) {
    // TODO use a tree structure for storing tour so that reverse is fast
    int best_change;
    do {
        best_change = 0.;
        unsigned best_i, best_j;
        int change;
        for (unsigned i = 0; i < tour.size() - 2; i++) {
            for (unsigned j = i + 2; j < tour.size(); j++) {
                // TODO use knn to narrow the search. now one step is O(n^2), with knn it will be O(kn)
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
