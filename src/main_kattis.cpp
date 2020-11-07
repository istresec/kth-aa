#include <iostream>
#include <array>

#include "TSP.h"
#include "utility.h"

using namespace std;

int main(int, char **) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    vector<int> tour;
    vector<pair<double, double>> cities;

    cin >> n;
    for (int i = 0; i < n; i++) {
        double x, y;
//        scanf(" %lf %lf", &x, &y);
        cin >> x >> y;
        cities.emplace_back(x, y);
    }

    tour = TSP::travel(cities);
    for (int i = 0; i < n; i++) {
        printf("%d\n", tour[i]);
    }

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