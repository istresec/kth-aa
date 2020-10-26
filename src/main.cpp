#include <iostream>
#include <chrono>
#include "TSP.h"
using namespace std;

int main(int, char **) { // (int argc, char **argv)
    std::cout << "Let's travel tha woooooooooooooooooooorld!" << std::endl;

    int n;
    vector<int> tour;

    cin >> n;
    vector<pair<double, double>> cities;

    // generate cities randomly
    cities = TSP::create_n_cities(n);
    cout << "Generated " << n << " cities:\n";
    for(auto & city: cities) {
        cout << get<0>(city) << ' ' << get<1>(city) << endl;
    }

    // input
//    for (int i = 0; i < n; i++) {
//        double x, y;
//        cin >> x >> y;
//        cities.emplace_back(x, y);
//    }

    cout << "Naive/greedy tour:\n";
    auto start = chrono::high_resolution_clock::now();
    tour = TSP::travel(cities);
    chrono::duration<double> elapsed = chrono::high_resolution_clock::now() - start;
    for (int i = 0; i < n; i++) {
        cout << tour[i] << endl;
    }
    cout << "Time elapsed: " << elapsed.count() << endl;
    cout << "Distance: " << TSP::tour_distance(cities, tour) << endl;

    cout << "Clarke Wright tour:\n";
    start = chrono::high_resolution_clock::now();
    tour = TSP::travel_cw(cities);
    elapsed = chrono::high_resolution_clock::now() - start;
    for (int i = 0; i < n; i++) {
        cout << tour[i] << endl;
    }
    cout << "Time elapsed: " << elapsed.count() << endl;
    cout << "Distance: " << TSP::tour_distance(cities, tour) << endl;

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