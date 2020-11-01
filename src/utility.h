#ifndef KTH_AA_UTILITY_H
#define KTH_AA_UTILITY_H

#include <memory>
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <random>

using namespace std;

template<class T>
class Grid { // https://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
    size_t _rows;
    size_t _columns;
    T *data;

public:
    Grid(size_t rows, size_t columns)
            : _rows{rows},
              _columns{columns},
              data{new T[rows * columns]} {}

    ~Grid() {
        delete[] data;
    }

    [[nodiscard]] size_t rows() const { return _rows; }

    [[nodiscard]] size_t columns() const { return _columns; }

    T *operator[](size_t row) { return row * _columns + data; }

    T &operator()(size_t row, size_t column) {
        return data[row * _columns + column];
    }

    void print() {
        for (int r = 0; r < rows(); r++) {
            for (int c = 0; c < columns(); c++) {
                cout << (*this)[r][c] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
};

inline double squared_distance(const pair<double, double> &t1, const pair<double, double> &t2) {
    double dx = t1.first - t2.first;
    double dy = t1.second - t2.second;
    return dx * dx + dy * dy;
}

inline Grid<int> *distance_matrix(const vector<pair<double, double>> &cities) {
    auto matrix = new Grid<int>(cities.size(), cities.size());
    for (unsigned i = 0; i < cities.size() - 1; i++) {
        for (unsigned j = i + 1; j < cities.size(); j++) {
            (*matrix)[i][j] = (*matrix)[j][i] = round(sqrt(squared_distance(cities[i], cities[j])));
        }
    }
    return matrix;
}

inline Grid<int> *k_nearest_neighbors(const vector<pair<double, double>> &cities, Grid<int> &distances, int k) {
    auto knn = new Grid<int>(cities.size(), k);

    auto cmp = [](const pair<int, int> &a, const pair<int, int> &b) { return a.first > b.first; };
    auto heap = vector<pair<int, int>>();
    heap.reserve(cities.size());

    for (unsigned i = 0; i < cities.size(); i++) {
        heap.clear();
        for (unsigned j = 0; j < cities.size(); j++) {
            if (i != j) { heap.emplace_back(distances[i][j], j); }
        }
        make_heap(heap.begin(), heap.end(), cmp);

        for (int neighbor = 0; neighbor < k; neighbor++) {
            (*knn)[i][neighbor] = heap.front().second;
            pop_heap(heap.begin(), heap.end(), cmp);
            heap.pop_back();
        }
    }
    return knn;
}

#endif //KTH_AA_UTILITY_H
