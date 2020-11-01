#include <iostream>
#include <chrono>
#include <array>

#include "TSP.h"
#include "utility.h"
#include "VariadicTable.h"

using namespace std;
using namespace chrono;

[[maybe_unused]] int kattis_demo() {
    int n;
    vector<int> tour;
    vector<pair<double, double>> cities;

    cin >> n;
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        cities.emplace_back(x, y);
    }

    tour = TSP::travel(cities);
    for (int i = 0; i < n; i++) {
        cout << tour[i] << endl;
    }

    return 0;
}

void demo_alg(const string &alg_name, vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt,
              vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &)) {
    cout << alg_name << ": ";

    auto start = chrono::high_resolution_clock::now();
    auto distances = distance_matrix(cities);
    auto tour = construction_alg(cities, *distances);

    chrono::duration<double> elapsed = chrono::high_resolution_clock::now() - start;
    int distance = TSP::tour_distance(cities, tour);
    for (int i: tour) {
        cout << i << " ";
    }
    cout << endl;
    vt.addRow(alg_name, distance, elapsed.count());
}

void demo_2opt(const string &alg_name, vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt,
               bool use_deadline, vector<int> (*construction_alg)(const vector<pair<double, double>> &, Grid<int> &)) {
    cout << alg_name;

    time_point<system_clock, duration<long, ratio<1, 1000000000>>> deadline = system_clock::now() + milliseconds(1900);
    auto start = chrono::high_resolution_clock::now();
    auto distances = distance_matrix(cities);
    Grid<int> *knn = distances; // TODO Dummy assignment. knn not implemented

    auto tour = construction_alg(cities, *distances);
    tour = TSP::local_2opt(tour, *distances, *knn, use_deadline ? &deadline : nullptr);

    chrono::duration<double> elapsed = chrono::high_resolution_clock::now() - start;
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

    /* Kattis */
    // return kattis_demo();

    std::cout << "Let's travel tha woooooooooooooooooooorld!" << std::endl;

    int n = 1000;
    vector<int> tour;
    auto cities = TSP::create_n_cities(n);
    cout << "Generated " << n << " cities:\n";
    for (auto &city: cities) {
        cout << get<0>(city) << ' ' << get<1>(city) << endl;
    }
    VariadicTable<string, int, double> vt({"Algo", "Distance", "Time elapsed"});

    demo_alg("Naive", cities, vt, TSP::travel_naive);
    demo_2opt("Naive 2opt", cities, vt, true, TSP::travel_naive);
    demo_2opt("Naive 2opt NO deadline", cities, vt, false, TSP::travel_naive);
    vt.print(cout);

    demo_alg("CW", cities, vt, TSP::travel_cw);
    demo_2opt("CW 2opt", cities, vt, true, TSP::travel_cw);
    demo_2opt("CW 2opt NO deadline", cities, vt, false, TSP::travel_cw);
    vt.print(cout);

    return 0;
}

/*

input:
10
95.0129 61.5432
23.1139 79.1937
60.6843 92.1813
48.5982 73.8207
89.1299 17.6266
76.2097 40.5706
45.6468 93.5470
1.8504 91.6904
82.1407 41.0270
44.4703 89.3650

naive_output:
0
8
5
4
3
9
6
2
1
7

 */