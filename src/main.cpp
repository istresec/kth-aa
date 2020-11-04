#include <iostream>
#include <chrono>
#include <array>

#include "TSP.h"
#include "VariadicTable.h"

using namespace std;
using namespace chrono;

bool DEBUG = false;
int DEADLINE = 1950;

void measure_preprocessing_time(vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt) {
    steady_clock::time_point start;
    duration<double> elapsed{};

    start = steady_clock::now();
    auto distances = distance_matrix(cities);
    elapsed = steady_clock::now() - start;
    vt.addRow("distance_matrix", -1, elapsed.count());

    start = steady_clock::now();
    k_nearest_neighbors<uint16_t>(*distances, 20);
    elapsed = steady_clock::now() - start;
    vt.addRow("knn k=20", -1, elapsed.count());

    start = steady_clock::now();
    k_nearest_neighbors<uint16_t>(*distances, 80);
    elapsed = steady_clock::now() - start;
    vt.addRow("knn k=80", -1, elapsed.count());

    start = steady_clock::now();
    k_nearest_neighbors<uint16_t>(*distances, 150);
    elapsed = steady_clock::now() - start;
    vt.addRow("knn k=150", -1, elapsed.count());

    start = steady_clock::now();
    k_nearest_neighbors<uint16_t>(*distances, 999);
    elapsed = steady_clock::now() - start;
    vt.addRow("knn k=999", -1, elapsed.count());
}

void demo_christofides(const string &alg_name, vector<pair<double, double>> &cities,
                       VariadicTable<string, int, double> &vt, bool use2opt = false, bool use3opt = false, int k = 80) {
    if (DEBUG) {
        cout << alg_name << ": ";
    }
    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(DEADLINE);
    auto start = steady_clock::now();
    auto distances = distance_matrix(cities);
    auto tour = TSP::travel_christofides<uint16_t, int>(*distances);

    Grid<uint16_t> *knn;
    if (use2opt or use3opt) {
        knn = k_nearest_neighbors<uint16_t>(*distances, k);
        if (use2opt) {
            tour = TSP::local_2opt(tour, *distances, *knn, &deadline);
        }
        if (use3opt) {
            tour = TSP::local_3opt_no_knn_sequential(tour, *distances, *knn, &deadline);
        }
    }

    duration<double> elapsed = steady_clock::now() - start;
    int distance = tour_distance(*distances, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }
    vt.addRow(alg_name, distance, elapsed.count());
}

void demo_cw_alg(const string &alg_name, vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt,
                 vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &, int), int hub = 0) {
    if (DEBUG) {
        cout << alg_name << ": ";
    }

    auto start = steady_clock::now();
    auto distances = distance_matrix(cities);
    auto tour = construction_alg(cities, *distances, hub);

    duration<double> elapsed = steady_clock::now() - start;
    int distance = tour_distance(*distances, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }
    vt.addRow(alg_name, distance, elapsed.count());
}

void demo_cw_opt(
        string alg_name,
        vector<pair<double, double>> &cities,
        VariadicTable<string, int, double> &vt,
        bool use_deadline,
        int k,
        vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &, int),
        vector<int> (*opt_alg)(vector<int>, Grid<int> &, Grid<uint16_t> &,
                               time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *),
        vector<int> (*opt_alg2)(vector<int>, Grid<int> &, Grid<uint16_t> &,
                                time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *) = nullptr) {
    alg_name = use_deadline ? alg_name : alg_name + " NO deadline";
    if (DEBUG) {
        cout << alg_name;
    }

    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(DEADLINE);
    auto start = steady_clock::now();
    auto distances = distance_matrix(cities);
    Grid<uint16_t> *knn = k_nearest_neighbors<uint16_t>(*distances, k);

    auto tour = construction_alg(cities, *distances, 0);
    tour = opt_alg(tour, *distances, *knn, use_deadline ? &deadline : nullptr);
    if (opt_alg2 != nullptr) {
        tour = opt_alg(tour, *distances, *knn, use_deadline ? &deadline : nullptr);
    }

    duration<double> elapsed = steady_clock::now() - start;
    int distance = tour_distance(*distances, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }
    vt.addRow(alg_name, distance, elapsed.count());
}

int main(int, char **) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    std::cout << "Let's travel tha woooooooooooooooooooorld!" << std::endl;

    int n = 1000;
    vector<int> tour;
    auto cities = create_n_cities(n, 720);
    VariadicTable<string, int, double> vt({"Algo", "Distance", "Time elapsed"});

    // logging
    if (DEBUG) {
        cout << "Generated " << n << " cities:\n";
        for (unsigned i = 0; i < cities.size(); i++) {
            cout << i << ' ' << get<0>(cities[i]) << ' ' << get<1>(cities[i]) << endl;
        }
        cout << "Distance matrix:" << endl;
        auto distances = distance_matrix(cities);
        distances->print();
        cout << "KNN matrix:" << endl;
        k_nearest_neighbors<uint16_t>(*distances, 10)->print();
    }

//    measure_preprocessing_time(cities, vt);
//    vt.print(cout);

//    demo_alg("Naive", cities, vt, TSP::travel_naive);
//    demo_cw_opt("Naive 2opt-20", cities, vt, true, 20, TSP::travel_naive, TSP::local_2opt);
//    demo_cw_opt("Naive 2opt-80", cities, vt, true, 80, TSP::travel_naive, TSP::local_2opt);
//    demo_cw_opt("Naive 3opt", cities, vt, true, 20, TSP::travel_naive, TSP::local_3opt);
//    vt.print(cout);

    demo_christofides("Christofides", cities, vt, false, false, 80);
//    demo_cw_alg("CW", cities, vt, TSP::travel_cw);
//    demo_cw_opt("CW 2opt-20", cities, vt, true, 20, TSP::travel_cw, TSP::local_2opt);
    demo_cw_opt("CW 2opt-80", cities, vt, true, 80, TSP::travel_cw, TSP::local_2opt);
//    demo_cw_opt("CW 2opt-150", cities, vt, true, 150, TSP::travel_cw, TSP::local_2opt);
//    demo_cw_opt("CW 2opt-999", cities, vt, false, 999, TSP::travel_cw, TSP::local_2opt);
//    demo_cw_opt("CW 2opt no knn", cities, vt, false, 20, TSP::travel_cw, TSP::local_2opt_no_knn);
//    demo_cw_opt("CW 3opt-20", cities, vt, true, 20, TSP::travel_cw, TSP::local_3opt);
//    demo_cw_opt("CW 3opt-80", cities, vt, true, 80, TSP::travel_cw, TSP::local_3opt);
//    demo_cw_opt("CW 3opt-150", cities, vt, true, 150, TSP::travel_cw, TSP::local_3opt);
//    demo_cw_opt("CW 3opt-200", cities, vt, true, 200, TSP::travel_cw, TSP::local_3opt);
//    demo_cw_opt("CW 3opt no knn", cities, vt, true, 20, TSP::travel_cw, TSP::local_3opt_no_knn);
//    demo_cw_opt("CW 3opt no knn sequential", cities, vt, true, 20, TSP::travel_cw, TSP::local_3opt_no_knn_sequential);
//    demo_cw_opt("CW 2opt-20 + 3opt-20", cities, vt, true, 20, TSP::travel_cw,
//                TSP::local_2opt, TSP::local_3opt);
//    demo_cw_opt("CW 2opt-20 + 3opt no knn sequential", cities, vt, true, 20, TSP::travel_cw,
//                TSP::local_2opt, TSP::local_3opt_no_knn_sequential);
    vt.print(cout);

    return 0;
}
