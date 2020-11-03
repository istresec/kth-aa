#include <iostream>
#include <chrono>
#include <array>

#include "TSP.h"
#include "VariadicTable.h"

using namespace std;
using namespace chrono;

bool DEBUG = false;

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

void demo_alg(const string &alg_name, vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt,
              vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &)) {
    if (DEBUG) {
        cout << alg_name << ": ";
    }

    auto start = steady_clock::now();
    auto distances = distance_matrix(cities);
    auto tour = construction_alg(cities, *distances);

    duration<double> elapsed = steady_clock::now() - start;
    int distance = TSP::tour_distance(cities, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }
    vt.addRow(alg_name, distance, elapsed.count());
}

void demo_2opt(
        const string &alg_name,
        vector<pair<double, double>> &cities,
        VariadicTable<string, int, double> &vt,
        bool use_deadline,
        int k,
        vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &),
        vector<int> (*opt_alg)(vector<int>, Grid<int> &, Grid<uint16_t> &,
                               time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *) = TSP::local_2opt) {
    if (DEBUG) {
        cout << alg_name;
    }

    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(1900);
    auto start = steady_clock::now();
    auto distances = distance_matrix(cities);
    Grid<uint16_t> *knn = k_nearest_neighbors<uint16_t>(*distances, k);

    auto tour = construction_alg(cities, *distances);
    tour = opt_alg(tour, *distances, *knn, use_deadline ? &deadline : nullptr);

    duration<double> elapsed = steady_clock::now() - start;
    int distance = TSP::tour_distance(cities, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }
    vt.addRow(alg_name, distance, elapsed.count());
}

void demo_3opt(
        const string &alg_name,
        vector<pair<double, double>> &cities,
        VariadicTable<string, int, double> &vt,
        bool use_deadline,
        vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &),
        vector<int> (*opt_alg)(vector<int>, Grid<int> &,
                               time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *) = TSP::local_3opt) {
    if (DEBUG) {
        cout << alg_name;
    }

    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + seconds(180);
    auto start = steady_clock::now();
    auto distances = distance_matrix(cities);

    auto tour = construction_alg(cities, *distances);
    tour = opt_alg(tour, *distances, use_deadline ? &deadline : nullptr);

    duration<double> elapsed = steady_clock::now() - start;
    int distance = TSP::tour_distance(cities, tour);
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
    auto cities = TSP::create_n_cities(n);
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

    measure_preprocessing_time(cities, vt);
    vt.print(cout);

//    demo_alg("Naive", cities, vt, TSP::travel_naive);
//    demo_2opt("Naive 2opt-20", cities, vt, true, 20, TSP::travel_naive);
//    demo_2opt("Naive 2opt-80", cities, vt, true, 80, TSP::travel_naive);
//    demo_3opt("Naive 3opt", cities, vt, true, TSP::travel_naive);
//    vt.print(cout);
    demo_alg("CW", cities, vt, TSP::travel_cw);
    demo_2opt("CW 2opt-20", cities, vt, false, 20, TSP::travel_cw);
    demo_2opt("CW 2opt-80", cities, vt, false, 80, TSP::travel_cw);
    demo_2opt("Naive 2opt no knn", cities, vt, false, 20, TSP::travel_cw, TSP::local_2opt_no_knn);
    demo_3opt("CW 3opt", cities, vt, false, TSP::travel_cw);
    vt.print(cout);

    // No deadlines (slow, therefore at end)
//    demo_2opt("Naive 2opt-20 NO deadline", cities, vt, false, 20, TSP::travel_naive);
//    vt.print(cout);
//    demo_2opt("Naive 2opt-80 NO deadline", cities, vt, false, 80, TSP::travel_naive);
//    vt.print(cout);
//    demo_2opt("Naive 2opt-999", cities, vt, false, 999, TSP::travel_naive);
//    vt.print(cout);
//    demo_2opt("CW 2opt-20 NO deadline", cities, vt, false, 20, TSP::travel_cw);
//    vt.print(cout);
//    demo_2opt("CW 2opt-80 NO deadline", cities, vt, false, 80, TSP::travel_cw);
//    vt.print(cout);
//    demo_2opt("CW 2opt-999", cities, vt, false, 999, TSP::travel_cw);
//    vt.print(cout);
//    demo_2opt("Naive 2opt no knn", cities, vt, false, 20, TSP::travel_naive, TSP::local_2opt_no_knn);
//    vt.print(cout);
//    demo_2opt("CW 2opt no knn", cities, vt, false, 20, TSP::travel_cw, TSP::local_2opt_no_knn);
//    vt.print(cout);

    return 0;
}
