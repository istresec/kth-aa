#include <iostream>
#include <chrono>
#include <array>

#include "TSP.h"
#include "utility.h"
#include "VariadicTable.h"

using namespace std;
using namespace chrono;

void measure_preprocessing_time(vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt) {
    system_clock::time_point start;
    duration<double> elapsed{};

    start = high_resolution_clock::now();
    auto distances = distance_matrix(cities);
    elapsed = high_resolution_clock::now() - start;
    vt.addRow("distance_matrix", -1, elapsed.count());

    start = high_resolution_clock::now();
    k_nearest_neighbors(cities, *distances, 20);
    elapsed = high_resolution_clock::now() - start;
    vt.addRow("knn k=20", -1, elapsed.count());

    start = high_resolution_clock::now();
    k_nearest_neighbors(cities, *distances, 80);
    elapsed = high_resolution_clock::now() - start;
    vt.addRow("knn k=80", -1, elapsed.count());

    start = high_resolution_clock::now();
    k_nearest_neighbors(cities, *distances, 150);
    elapsed = high_resolution_clock::now() - start;
    vt.addRow("knn k=150", -1, elapsed.count());

    start = high_resolution_clock::now();
    k_nearest_neighbors(cities, *distances, 999);
    elapsed = high_resolution_clock::now() - start;
    vt.addRow("knn k=999", -1, elapsed.count());
}

void demo_alg(const string &alg_name, vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt,
              vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &)) {
    cout << alg_name << ": ";

    auto start = high_resolution_clock::now();
    auto distances = distance_matrix(cities);
    auto tour = construction_alg(cities, *distances);

    duration<double> elapsed = high_resolution_clock::now() - start;
    int distance = TSP::tour_distance(cities, tour);
    for (int i: tour) {
        cout << i << " ";
    }
    cout << endl;
    vt.addRow(alg_name, distance, elapsed.count());
}

void demo_2opt(
        const string &alg_name,
        vector<pair<double, double>> &cities,
        VariadicTable<string, int, double> &vt,
        bool use_deadline,
        int k,
        vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &)) {
    cout << alg_name;

    time_point<system_clock, duration<long, ratio<1, 1000000000>>> deadline = system_clock::now() + milliseconds(1900);
    auto start = high_resolution_clock::now();
    auto distances = distance_matrix(cities);
    Grid<int> *knn = k_nearest_neighbors(cities, *distances, k);

    auto tour = construction_alg(cities, *distances);
    tour = TSP::local_2opt(tour, *distances, *knn, use_deadline ? &deadline : nullptr);

    duration<double> elapsed = high_resolution_clock::now() - start;
    int distance = TSP::tour_distance(cities, tour);
    for (int i: tour) {
        cout << i << " ";
    }
    cout << endl;
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
    {
        cout << "Generated " << n << " cities:\n";
        for (int i = 0; i < cities.size(); i++) {
            cout << i << ' ' << get<0>(cities[i]) << ' ' << get<1>(cities[i]) << endl;
        }
        cout << "Distance matrix:" << endl;
        auto distances = distance_matrix(cities);
        distances->print();
        cout << "KNN matrix:" << endl;
        k_nearest_neighbors(cities, *distances, 10)->print();
    }

    measure_preprocessing_time(cities, vt);
    vt.print(cout);

    demo_alg("Naive", cities, vt, TSP::travel_naive);
    demo_2opt("Naive 2opt-80", cities, vt, true, 999, TSP::travel_naive);
//    demo_2opt("Naive 2opt-80 NO deadline", cities, vt, false, 999, TSP::travel_naive);
    vt.print(cout);

    demo_alg("CW", cities, vt, TSP::travel_cw);
    demo_2opt("CW 2opt-80", cities, vt, true, 999, TSP::travel_cw);
    demo_2opt("CW 2opt-80 NO deadline", cities, vt, false, 999, TSP::travel_cw);
    vt.print(cout);

    return 0;
}
