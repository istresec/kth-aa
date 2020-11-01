#include <iostream>
#include <chrono>
#include <array>

#include "TSP.h"
#include "utility.h"
#include "VariadicTable.h"

using namespace std;
using namespace chrono;

int main(int, char **) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    std::cout << "Let's travel tha woooooooooooooooooooorld!" << std::endl;

    int n;
    vector<int> tour;

//    cin >> n;
    n = 1000;
    vector<pair<double, double>> cities;
    VariadicTable<std::string, int, double> vt({"Algo", "Distance", "Time elapsed"});

    // input
//    for (int i = 0; i < n; i++) {
//        double x, y;
//        cin >> x >> y;
//        cities.emplace_back(x, y);
//    }

    // generate cities randomly
    cities = TSP::create_n_cities(n);
    cout << "Generated " << n << " cities:\n";
    for (auto &city: cities) {
        cout << get<0>(city) << ' ' << get<1>(city) << endl;
    }

    /* Naive tour */
    cout << "Naive/greedy tour: ";
    auto start = chrono::high_resolution_clock::now();
    auto distances = distance_matrix(cities);
    tour = TSP::travel_naive(cities, *distances);

    chrono::duration<double> elapsed = chrono::high_resolution_clock::now() - start;
    int distance = TSP::tour_distance(cities, tour);
    for (int i = 0; i < n; i++) { cout << tour[i] << " "; }
    cout << endl;
    cout << "Time elapsed: " << elapsed.count() << endl;
    cout << "Distance: " << distance << endl;
    vt.addRow("Naive", distance, elapsed.count());
    vt.print(cout);


    /* Naive with 2opt */
    cout << "Naive with 2opt: ";
    time_point deadline = system_clock::now() + milliseconds(2000);
    Grid<int> *knn = distances; // TODO Dummy assignment. knn not implemented
    tour = TSP::local_2opt(tour, *distances, *knn, &deadline);

    distance = TSP::tour_distance(cities, tour);
    elapsed = chrono::high_resolution_clock::now() - start;
    for (int i = 0; i < n; i++) { cout << tour[i] << " "; }
    cout << endl;
    cout << "Time elapsed: " << elapsed.count() << endl;
    cout << "Distance: " << distance << endl;
    vt.addRow("Naive 2opt", distance, elapsed.count());
    vt.print(cout);

    tour = TSP::local_2opt(tour, *distances, *knn, nullptr);
    distance = TSP::tour_distance(cities, tour);
    elapsed = chrono::high_resolution_clock::now() - start;
    vt.addRow("Naive 2opt no deadline", distance, elapsed.count());
    vt.print(cout);


    /* Clarke Wright */
    cout << "Clarke Wright tour:\n";
    start = chrono::high_resolution_clock::now();
    tour = TSP::travel_cw(cities);

    elapsed = chrono::high_resolution_clock::now() - start;
    distance = TSP::tour_distance(cities, tour);
    for (int i = 0; i < n; i++) { cout << tour[i] << " "; }
    cout << endl;
    cout << "Time elapsed: " << elapsed.count() << endl;
    cout << "Distance: " << distance << endl;
    vt.addRow("Clarke Wright", distance, elapsed.count());
    vt.print(cout);


    /* Clarke Wright with 2opt */
    cout << "Clarke Wright with 2opt: ";
    tour = TSP::local_2opt(tour, *distances, *knn, nullptr);

    elapsed = chrono::high_resolution_clock::now() - start;
    distance = TSP::tour_distance(cities, tour);
    for (int i = 0; i < n; i++) { cout << tour[i] << " "; }
    cout << endl;
    cout << "Time elapsed: " << elapsed.count() << endl;
    cout << "Distance: " << distance << endl;
    vt.addRow("Clarke Wright 2opt no deadline", distance, elapsed.count());
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