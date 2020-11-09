#include <iostream>
#include <chrono>
#include <array>

#include "TSP.h"
//#include "TSP.cpp"
#include "utility.h"
#include "VariadicTable.h"

using namespace std;
using namespace chrono;

bool DEBUG = false;
int DEADLINE = 1950;
bool LOG = true;
string LOG_PATH = "/mnt/terra/xoding/kth-aa/logs/002.log";

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
                       VariadicTable<string, int, double> &vt, int k = 80, bool use2opt = false, bool use3opt = false,
                       bool use_lin_kernighan = false) {
    if (DEBUG) {
        cout << alg_name << ": ";
    }
    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(DEADLINE);
    auto start = steady_clock::now();
    auto distances = distance_matrix<double>(cities);
    auto tour = travel_christofides<uint16_t, double>(*distances);
    tour = vector<int>({0, 1, 3, 2}); // TODO remove

    if (LOG) {
        log_cities(cities, LOG_PATH, alg_name);
    }
    if (LOG) {
        auto distance = tour_distance<double>(*distances, tour);
        log_tour(tour, LOG_PATH, " 0 --> bare " + to_string(distance));
    }

    Grid<uint16_t> *knn;
    if (use2opt or use3opt or use_lin_kernighan) {
        knn = k_nearest_neighbors<uint16_t>(*distances, k);
        if (use2opt) {
            tour = local_2opt(tour, *distances, *knn, &deadline);
            if (LOG) {
                int distance = tour_distance(*distances, tour);
                log_tour(tour, LOG_PATH, " 1 --> 2opt " + to_string(distance));
            }
        }
        if (use3opt) {
            tour = local_3opt_no_knn_sequential(tour, *distances, *knn, &deadline);
//            tour = local_3opt(tour, *distances, *knn, &deadline);
            if (LOG) {
                int distance = tour_distance(*distances, tour);
                log_tour(tour, LOG_PATH, " 2 --> 3opt " + to_string(distance));
            }
        }
        if (use_lin_kernighan) {
//            tour = chokolino(tour, *distances, *knn, &deadline);
//            if (LOG) {
//                int distance = tour_distance(*distances, tour);
//                log_tour(tour, LOG_PATH, " 3 --> CHOKOLINO " + to_string(distance));
//            }
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

    if (LOG) {
        log_cities(cities, LOG_PATH, alg_name);
    }
    if (LOG) {
        log_tour(tour, LOG_PATH, "d=" + to_string(distance));
    }
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
    if (LOG) {
        log_cities(cities, LOG_PATH, alg_name);
    }
    if (LOG) {
        int distance = tour_distance(*distances, tour);
        log_tour(tour, LOG_PATH, "0 --> bare d=" + to_string(distance));
    }
    tour = opt_alg(tour, *distances, *knn, use_deadline ? &deadline : nullptr);
    if (LOG) {
        int distance = tour_distance(*distances, tour);
        log_tour(tour, LOG_PATH, "1 --> opt_1 d=" + to_string(distance));
    }
    if (opt_alg2 != nullptr) {
        tour = opt_alg(tour, *distances, *knn, use_deadline ? &deadline : nullptr);
        if (LOG) {
            int distance = tour_distance(*distances, tour);
            log_tour(tour, LOG_PATH, "2 --> opt_2 d=" + to_string(distance));
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

int main(int, char **) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    std::cout << "Let's travel tha woooooooooooooooooooorld!" << std::endl;

    int n = 1000;
    vector<int> tour;
//    auto cities = create_n_cities(n, 720);
    VariadicTable<string, int, double> vt({"Algo", "Distance", "Time elapsed"});

    cin >> n;
    vector<pair<double, double>> cities;
    for (int i = 0; i < n; i++) {
        double x, y;
        // scanf(" %lf %lf", &x, &y);
        cin >> x >> y;
        cities.emplace_back(x, y);
    }

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

//    demo_alg("Naive", cities, vt, travel_naive);
//    demo_cw_opt("Naive 2opt-20", cities, vt, true, 20, travel_naive, local_2opt);
//    demo_cw_opt("Naive 2opt-80", cities, vt, true, 80, travel_naive, local_2opt);
//    demo_cw_opt("Naive 3opt", cities, vt, true, 20, travel_naive, local_3opt);
//    vt.print(cout);

//    demo_christofides("Christofides", cities, vt, 150,false, false);
//    demo_christofides("Christofides 2opt-150", cities, vt, 150,true, false);
    demo_christofides("Christofides 2opt-150 3opt_no_knn_sequential", cities, vt, 150,true, true);
//    demo_christofides("Christofides 3opt-150_no_knn_sequential", cities, vt, 150,false, true);
    demo_christofides("Christofides chokolino-3", cities, vt, 3, false, false, true);
//    demo_cw_alg("CW", cities, vt, travel_cw);
//    demo_cw_opt("CW 2opt-20", cities, vt, true, 20, travel_cw, local_2opt);
//    demo_cw_opt("CW 2opt-80", cities, vt, true, 80, travel_cw, local_2opt);
//    demo_cw_opt("CW 2opt-150", cities, vt, true, 150, travel_cw, local_2opt);
//    demo_cw_opt("CW 2opt-999", cities, vt, false, 999, travel_cw, local_2opt);
//    demo_cw_opt("CW 2opt no knn", cities, vt, false, 20, travel_cw, local_2opt_no_knn);
//    demo_cw_opt("CW 3opt-20", cities, vt, true, 20, travel_cw, local_3opt);
//    demo_cw_opt("CW 3opt-80", cities, vt, true, 80, travel_cw, local_3opt);
//    demo_cw_opt("CW 3opt-150", cities, vt, true, 150, travel_cw, local_3opt);
//    demo_cw_opt("CW 3opt-200", cities, vt, true, 200, travel_cw, local_3opt);
//    demo_cw_opt("CW 3opt no knn", cities, vt, true, 20, travel_cw, local_3opt_no_knn);
//    demo_cw_opt("CW 3opt 150 knn sequential NO DEAD", cities, vt, false, 150, travel_cw, local_3opt_no_knn_sequential);
//    demo_cw_opt("CW 3opt 150 knn sequential", cities, vt, true, 150, travel_cw, local_3opt_no_knn_sequential);
//    demo_cw_opt("CW 2opt-150 + 3opt-150", cities, vt, true, 150, travel_cw,
//                local_2opt, local_3opt);
//    demo_cw_opt("CW 2opt-150 + 3opt no knn sequential", cities, vt, true, 150, travel_cw,
//                local_2opt, local_3opt_no_knn_sequential);
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