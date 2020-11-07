#include <memory>
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <random>

#include "utility.h"

using namespace std;

vector<pair<double, double>> create_n_cities(int n, int seed) {
    double lower_bound = -1e6;
    double upper_bound = 1e6;
    vector<pair<double, double >> cities = vector<pair<double, double >>();
    cities.reserve(n);

    uniform_real_distribution<double> unif(lower_bound, upper_bound);
    default_random_engine generator(seed);
    for (int i = 0; i < n; i++) {
        cities.emplace_back(pair<double, double>(unif(generator), unif(generator)));
    }
    return cities;
}

int tour_distance(Grid<int> &distances, vector<int> tour) {
    int distance = 0;
    for (unsigned i = 0; i < tour.size(); i++) {
        distance += distances[tour[i]][tour[(i + 1) % distances.rows()]];
    }
    return distance;
}

void log_cities(const vector<pair<double, double>> &cities, const string &path, const string &name) {
    ofstream outfile(path, ios_base::app);
    outfile << "\ncities_name=" << name << "\n";
    outfile << "cities=[";
    for (unsigned i = 0; i < cities.size(); i++) {
        outfile << "(" << cities[i].first << "," << cities[i].second << "),";
    }
    outfile << "]\n";
    outfile.close();
}
