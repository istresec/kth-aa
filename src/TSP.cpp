#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <random>

#include "TSP.h"
#include "utility.h"

using namespace std;

inline bool compare_savings(tuple<int, int, int> t1, tuple<int, int, int> t2) {
    return get<2>(t1) > get<2>(t2);
}

vector<tuple<int, int, int>> savings(const vector<pair<double, double>> &cities) {
    vector<tuple<int, int, int>> savings = vector<tuple<int, int, int>>();
    Grid<int> *tsp_table = distance_matrix(cities);
    for (unsigned i = 1; i < cities.size() - 1; i++) {
        for (unsigned j = i + 1; j < cities.size(); j++) {
            savings.emplace_back(make_tuple(i, j, (*tsp_table)[0][i] + (*tsp_table)[0][j] - (*tsp_table)[i][j]));
        }
    }
    sort(savings.begin(), savings.end(), compare_savings);
    return savings;
}

vector<int> TSP::travel(const vector<pair<double, double>> &cities) {
    return TSP::travel_naive(cities);
}

vector<int> TSP::travel_naive(vector<pair<double, double>> cities) {
    vector<int> tour = vector<int>(cities.size());
    vector<bool> used = vector<bool>(cities.size());
    tour[0] = 0;
    used[0] = true;
    for (unsigned i = 1; i <= cities.size() - 1; i++) {
        int best = -1;
        for (unsigned j = 0; j <= cities.size() - 1; j++) {
            bool is_better = best == -1 or squared_distance(cities[tour[i - 1]], cities[j]) <
                                           squared_distance(cities[tour[i - 1]], cities[best]);
            if (not used[j] and is_better) {
                best = (int) j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
    return tour;
}

// Clarke-Wright Savings Algorithms. City at index 0 used as hub.
vector<int> TSP::travel_cw(const vector<pair<double, double>> &cities) {
    // get savings
    vector<tuple<int, int, int>> s = savings(cities);

    // initialize tours
    vector<vector<int>> tours = vector<vector<int>>();
    for (int i = 1; i < cities.size(); i++)
        tours.emplace_back(vector<int>{i}); // instead of 0, i, 0 just use i

    // algorithm
    vector<int> indices;
    vector<int> temp_tour;
    int i, j, endpoint;
    for (auto &it : s) {
        indices = vector<int>();
        i = get<0>(it);
        j = get<1>(it);
        for (unsigned k = 0; k < tours.size(); k++) {
            endpoint = tours[k].back();
            if (endpoint == i or endpoint == j) {
                indices.emplace_back(k);
                if (indices.size() == 2) {
                    temp_tour = tours[indices[0]];
                    temp_tour.insert(temp_tour.end(), tours[indices[1]].begin(), tours[indices[1]].end());
                    tours[indices[0]] = temp_tour;
                    tours.erase(tours.begin() + indices[1]);
                    break;
                }
            }
        }
    }

    // create final tour
    vector<int> tour = vector<int>();
    tour.emplace_back(0);
    tour.insert(tour.end(), tours[0].begin(), tours[0].end());

    cout << "Returning Clarke Wright.\n";
    return tour;
}

int TSP::tour_distance(const vector<pair<double, double>> &cities, vector<int> tour) {
    int distance = 0;
    auto tsp_table = distance_matrix(cities);
    for (unsigned i = 0; i < tour.size() - 1; i++) {
        distance += (*tsp_table)[tour[i]][tour[i + 1]];
    }
    return distance;
}

vector<pair<double, double>> TSP::create_n_cities(int n) {
    double lower_bound = 0;
    double upper_bound = (double) n * 10;
    vector<pair<double, double>> cities = vector<pair<double, double>>();
    cities.reserve(n);

    uniform_real_distribution<double> unif(lower_bound, upper_bound);
    default_random_engine re;
    for (int i = 0; i < n; i++) {
        cities.emplace_back(pair<double, double>(unif(re), unif(re)));
    }
    return cities;
}
