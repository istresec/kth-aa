#include <iostream>
#include <chrono>
#include <array>

#include "greedy.h"
#include "clarke_wright.h"
#include "christofides.h"
#include "local2opt.h"
//#include "local2opt_no_knn.h"
//#include "local3opt.h"
#include "local3opt_no_knn_sequential.h"
#include "chokolino.h"
#include "utility.h"
#include "VariadicTable.h"

using namespace std;
using namespace chrono;

bool DEBUG = true;
int DEADLINE = 1950;
bool LOG = false;
string LOG_PATH = "/mnt/terra/xoding/kth-aa/logs/002.log";

void measurePreprocessingTime(vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt) {
    steady_clock::time_point start;
    duration<double> elapsed{};

    start = steady_clock::now();
    auto distances = distanceMatrix(cities);
    elapsed = steady_clock::now() - start;
    vt.addRow("distanceMatrix", -1, elapsed.count());

    start = steady_clock::now();
    kNearestNeighbors<uint16_t>(*distances, 20);
    elapsed = steady_clock::now() - start;
    vt.addRow("knn k=20", -1, elapsed.count());

    start = steady_clock::now();
    kNearestNeighbors<uint16_t>(*distances, 80);
    elapsed = steady_clock::now() - start;
    vt.addRow("knn k=80", -1, elapsed.count());

    start = steady_clock::now();
    kNearestNeighbors<uint16_t>(*distances, 150);
    elapsed = steady_clock::now() - start;
    vt.addRow("knn k=150", -1, elapsed.count());

    start = steady_clock::now();
    kNearestNeighbors<uint16_t>(*distances, 999);
    elapsed = steady_clock::now() - start;
    vt.addRow("knn k=999", -1, elapsed.count());
}

template<class T=uint16_t, class U=int>
void
demoGreedy(const string &algName, vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt,
           vector<T> (*constructionAlg)(const vector<pair<double, double>> &, Grid<U> &), int k = 80,
           bool useLinKernighan = false, bool use2Opt = false, bool use3Opt = false,
           bool useDeadline = true) {
    if (DEBUG) {
        cout << algName << ": ";
    }
    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(DEADLINE);

    auto start = steady_clock::now();
    auto distances = distanceMatrix(cities);
    auto tour = constructionAlg(cities, *distances);

    if (LOG) {
        logCities(cities, LOG_PATH, algName);
    }
    if (LOG) {
        auto distance = tourDistance(*distances, tour);
        logTour(tour, LOG_PATH, " 0 --> bare " + to_string(distance));
    }

    Grid<uint16_t> *knn;
    if (use2Opt or use3Opt or useLinKernighan) {
        knn = kNearestNeighbors<uint16_t>(*distances, k);
        if (useLinKernighan) {
            tour = chokolino(tour, *distances, *knn, useDeadline ? &deadline : nullptr);
            if (LOG) {
                int distance = tourDistance(*distances, tour);
                logTour(tour, LOG_PATH, " 3 --> CHOKOLINO " + to_string(distance));
            }
        }
        if (use2Opt) {
            tour = local2Opt(tour, *distances, *knn, useDeadline ? &deadline : nullptr);
            if (LOG) {
                int distance = tourDistance(*distances, tour);
                logTour(tour, LOG_PATH, " 1 --> 2opt " + to_string(distance));
            }
        }
        if (use3Opt) {
            tour = local3OptNoKnnSequential(tour, *distances, *knn, useDeadline ? &deadline : nullptr);
            if (LOG) {
                int distance = tourDistance(*distances, tour);
                logTour(tour, LOG_PATH, " 2 --> 3opt " + to_string(distance));
            }
        }
    }

    duration<double> elapsed = steady_clock::now() - start;
    int distance = tourDistance(*distances, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }

    vt.addRow(algName, distance, elapsed.count());
}

void demoChristofides(const string &algName, vector<pair<double, double>> &cities,
                      VariadicTable<string, int, double> &vt, int k = 80, bool use2Opt = false, bool use3Opt = false,
                      bool useLinKernighan = false) {
    if (DEBUG) {
        cout << algName << ": ";
    }
    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(DEADLINE);
    auto start = steady_clock::now();
    auto distances = distanceMatrix<double>(cities);
    auto tour = travelChristofides<uint16_t, double>(*distances);

    if (LOG) {
        logCities(cities, LOG_PATH, algName);
    }
    if (LOG) {
        auto distance = tourDistance(*distances, tour);
        logTour(tour, LOG_PATH, " 0 --> bare " + to_string(distance));
    }

    Grid<uint16_t> *knn;
    if (use2Opt or use3Opt or useLinKernighan) {
        knn = kNearestNeighbors<uint16_t>(*distances, k);
        if (useLinKernighan) {
            tour = chokolino(tour, *distances, *knn, &deadline);
            if (LOG) {
                int distance = tourDistance(*distances, tour);
                logTour(tour, LOG_PATH, " 3 --> CHOKOLINO " + to_string(distance));
            }
        }
        if (use2Opt) {
            tour = local2Opt(tour, *distances, *knn, &deadline);
            if (LOG) {
                int distance = tourDistance(*distances, tour);
                logTour(tour, LOG_PATH, " 1 --> 2opt " + to_string(distance));
            }
        }
        if (use3Opt) {
            tour = local3OptNoKnnSequential(tour, *distances, *knn, &deadline);
//            tour = local3Opt(tour, *distances, *knn, &deadline);
            if (LOG) {
                int distance = tourDistance(*distances, tour);
                logTour(tour, LOG_PATH, " 2 --> 3opt " + to_string(distance));
            }
        }
    }

    duration<double> elapsed = steady_clock::now() - start;
    int distance = tourDistance(*distances, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }

    vt.addRow(algName, distance, elapsed.count());
}

void demoCwAlg(const string &algName, vector<pair<double, double>> &cities, VariadicTable<string, int, double> &vt,
               vector<int> (*constructionAlg)(const vector<pair<double, double>> &, Grid<int> &, int), int hub = 0) {
    if (DEBUG) {
        cout << algName << ": ";
    }

    auto start = steady_clock::now();
    auto distances = distanceMatrix(cities);
    auto tour = constructionAlg(cities, *distances, hub);

    duration<double> elapsed = steady_clock::now() - start;
    int distance = tourDistance(*distances, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }
    vt.addRow(algName, distance, elapsed.count());

    if (LOG) {
        logCities(cities, LOG_PATH, algName);
    }
    if (LOG) {
        logTour(tour, LOG_PATH, "d=" + to_string(distance));
    }
}

template<class T, class U>
void demoCwOpt(
        string algName,
        vector<pair<double, double>> &cities,
        VariadicTable<string, int, double> &vt,
        bool useDeadline,
        int k,
        vector<T> (*constructionAlg)(const vector<pair<double, double>> &, Grid<U> &, int),
        vector<T> (*optAlg)(vector<T>, Grid<U> &, Grid<T> &,
                            time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *) = nullptr,
        vector<T> (*optAlg2)(vector<T>, Grid<U> &, Grid<T> &,
                             time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *) = nullptr,
        vector<T> (*optAlgChokolino)(
                vector<T> &, Grid<U> &, Grid<T> &,
                time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *) = nullptr) {
    algName = useDeadline ? algName : algName + " NO deadline";
    if (DEBUG) {
        cout << algName;
    }

    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(DEADLINE);
    auto start = steady_clock::now();
    auto distances = distanceMatrix(cities);
    Grid<T> *knn = kNearestNeighbors<T>(*distances, k);

    auto tour = constructionAlg(cities, *distances, 0);
    if (LOG) {
        logCities(cities, LOG_PATH, algName);
    }
    if (LOG) {
        int distance = tourDistance(*distances, tour);
        logTour(tour, LOG_PATH, "0 --> bare d=" + to_string(distance));
    }
    if (optAlg != nullptr) {
        tour = optAlg(tour, *distances, *knn, useDeadline ? &deadline : nullptr);
        if (LOG) {
            int distance = tourDistance(*distances, tour);
            logTour(tour, LOG_PATH, "1 --> opt_1 d=" + to_string(distance));
        }
    }
    if (optAlg2 != nullptr) {
        tour = optAlg(tour, *distances, *knn, useDeadline ? &deadline : nullptr);
        if (LOG) {
            int distance = tourDistance(*distances, tour);
            logTour(tour, LOG_PATH, "2 --> opt_2 d=" + to_string(distance));
        }
    }
    if (optAlgChokolino != nullptr) {
        tour = optAlgChokolino(tour, *distances, *knn, useDeadline ? &deadline : nullptr);
        if (LOG) {
            int distance = tourDistance(*distances, tour);
            logTour(tour, LOG_PATH, "3 --> opt_2 d=" + to_string(distance));
        }
    }

    duration<double> elapsed = steady_clock::now() - start;
    int distance = tourDistance(*distances, tour);
    if (DEBUG) {
        for (int i: tour) {
            cout << i << " ";
        }
        cout << endl;
    }
    vt.addRow(algName, distance, elapsed.count());
}

int main(int, char **) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    std::cout << "Let's travel tha woooooooooooooooooooorld!" << std::endl;

    string citiesName = "Random -- ";
    int n = 1000;
    vector<int> tour;
    auto cities = createNCities(n, 721);
    VariadicTable<string, int, double> vt({"Algo", "Distance", "Time elapsed"});

//    getline(cin, cities_name);
//    cin >> n;
//    cities.clear();
//    for (int i = 0; i < n; i++) {
//        double x, y;
//        // scanf(" %lf %lf", &x, &y);
//        cin >> x >> y;
//        cities.emplace_back(x, y);
//    }

    // logging
    if (DEBUG) {
        cout << "Generated " << n << " cities:\n";
        for (unsigned i = 0; i < cities.size(); i++) {
            cout << i << ' ' << get<0>(cities[i]) << ' ' << get<1>(cities[i]) << endl;
        }
        cout << "Distance matrix:" << endl;
        auto distances = distanceMatrix(cities);
        distances->print();
        cout << "KNN matrix:" << endl;
        kNearestNeighbors<uint16_t>(*distances, 10)->print();
    }

    measurePreprocessingTime(cities, vt);
//    vt.print(cout);

    demoGreedy(citiesName + "Greedy", cities, vt, travelGreedy);
//    demoGreedy("Greedy 2opt-100", cities, vt, travelGreedy, 100, false, true, false);
//    demoGreedy("Greedy 2opt-80", cities, vt, travelGreedy, 80, false, true, false);
//    demoGreedy("Greedy 3opt-100", cities, vt, travelGreedy, 100, false, false, true);
//    demoGreedy("Greedy chokolino-100", cities, vt, travelGreedy, 100, true, false, false);
//    vt.print(cout);

    demoCwAlg(citiesName + "CW", cities, vt, travel_cw);
//    demoCwOpt<uint16_t, int>("CW 2opt no knn", cities, vt, true, 20, travel_cw, local2OptNoKnn);
//    demoCwOpt<uint16_t, int>("CW 2opt-100", cities, vt, true, 100, travel_cw, local2Opt);
//    demoCwOpt<uint16_t, int>("CW 3opt", cities, vt, true, 100, travel_cw, local3OptNoKnnSequential);
    demoCwOpt<uint16_t, int>("CW chokolino-100", cities, vt, true, 100, travel_cw, nullptr, nullptr, chokolino);
//    demoCwOpt<uint16_t, int>("CW 3opt", cities, vt, true, 150, travel_cw, local3OptNoKnnSequential);
//    vt.print(cout);

    demoChristofides(citiesName + "Christofides", cities, vt, 100, false, false);
    demoChristofides("Christofides 2opt-100", cities, vt, 100, true, false);
    demoChristofides("Christofides 2opt-100 3opt_no_knn_sequential", cities, vt, 100, true, true);
    demoChristofides("Christofides 3opt_no_knn_sequential", cities, vt, 100, false, true);
//    demoChristofides("Christofides 2opt-chokolino-100", cities, vt, 100, true, false, true);
//    vt.print(cout);

    demoChristofides(citiesName + "Christofides chokolino-100", cities, vt, 100, false, false, true);
//    cities = createNCities(n, 722);
//    demoChristofides("c2: Christofides chokolino-100", cities, vt, 100, false, false, true);
//    cities = createNCities(n, 723);
//    demoChristofides("c3: Christofides chokolino-100", cities, vt, 100, false, false, true);
//    cities = createNCities(n, 724);
//    demoChristofides("c4: Christofides chokolino-100", cities, vt, 100, false, false, true);

    demoChristofides(citiesName + "Christofides chokolino-2opt-3opt k=100", cities, vt, 100, true, true, true);
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