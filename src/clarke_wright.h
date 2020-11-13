#ifndef KTH_AA_CLARKE_WRIGHT_H
#define KTH_AA_CLARKE_WRIGHT_H

#include <vector>
#include <tuple>
#include "utility.h"

using namespace std;

template<class T>
inline bool compareSavings(tuple<int, int, T> t1, tuple<int, int, T> t2) {
    return get<2>(t1) > get<2>(t2);
}

template<class T>
vector<tuple<int, int, T>> savings(const vector<pair<double, double>> &cities, Grid<T> &distances, int hub) {
    vector<tuple<int, int, T >> savings = vector<tuple<int, int, T >>();
    for (T i = 0; i < cities.size() - 1; i++) {
        if (i == (T) hub) continue;
        for (unsigned j = i + 1; j < cities.size(); j++) {
            savings.emplace_back(make_tuple(i, j, distances[hub][i] + distances[hub][j] - distances[i][j]));
        }
    }
    sort(savings.begin(), savings.end(), compareSavings<T>);
    return savings;
}

// Clarke-Wright Savings Algorithm. City at index 0 used as hub.
template<class T, class U>
vector<T> travelCw(const vector<pair<double, double>> &cities, Grid<U> &distances, int hub) {
    // if only one city
    if (cities.size() == 1)
        return vector<T>{0};

    vector<tuple<int, int, U>> s = savings(cities, distances, hub);

    // initialize tours
    vector<T> tours[cities.size()];
    vector<T> emptyTour = vector<T>();
    for (T i = 0; i < (int) cities.size(); i++) {
        if (i == hub) {
            tours[i] = emptyTour;
        } else {
            tours[i] = vector<T>{i}; // instead of 0, i, 0 just use i
        }
    }

    // algorithm
    vector<T> tempTour;
    vector<T> first, second;
    int i, j;
    for (auto &it : s) {
        i = get<0>(it);
        j = get<1>(it);
        // if two distinct tours with endpoints i and j exist, if so combine them (O(1))
        if (not tours[i].empty() and not tours[j].empty() and tours[i].front() != tours[j].front()) {
            first = tours[i]; // remember tour with endpoint i
            second = tours[j]; // remember tour with endpoint j
            // remove tours with endpoints i and j while a new tour is constructed
            tours[first.front()] = tours[first.back()] = tours[second.front()] = tours[second.back()] = emptyTour;

            if (first.front() == i)
                reverse(first.begin(), first.end()); // reverse tour with endpoint i if it starts with i
            if (second.front() != j)
                reverse(second.begin(), second.end()); // reverse tour with j if it doesn't start with j

            // create new tour by joining the two
            first.insert(first.end(), second.begin(), second.end());

            // remember endpoints of the new tour in array of endpoints for quick access
            tours[first.front()] = first;
            tours[first.back()] = first;
        }
    }

    // create final tour
    vector<T> tour = vector<T>();
    for (i = 1; i < (int) cities.size(); i++) {
        if (not tours[i].empty())
            break;
    }

    tour.emplace_back(hub);
    tour.insert(tour.end(), tours[i].begin(), tours[i].end());

    return tour;
}

#endif //KTH_AA_CLARKE_WRIGHT_H
