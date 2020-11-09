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

    tour = travel(cities);
    for (int i = 0; i < n; i++) {
        printf("%d\n", tour[i]);
    }

    return 0;
}