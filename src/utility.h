#ifndef KTH_AA_UTILITY_H
#define KTH_AA_UTILITY_H

#include <memory>
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

struct pair_hash {
public:
    template<typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
};


template<class T>
class Grid { // https://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
    size_t rows_;
    size_t columns_;
    T *data;

public:
    Grid(size_t rows, size_t columns)
            : rows_{rows},
              columns_{columns},
              data{new T[rows * columns]} {}

    ~Grid() {
        delete[] data;
    }

    [[nodiscard]] size_t rows() const { return rows_; }

    [[nodiscard]] size_t columns() const { return columns_; }

    T *operator[](size_t row) { return row * columns_ + data; }

    T &operator()(size_t row, size_t column) {
        return data[row * columns_ + column];
    }

    void print() {
        for (size_t r = 0; r < rows(); r++) {
            for (size_t c = 0; c < columns(); c++) {
                cout << (*this)[r][c] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
};


class UnionFind {
public:
    explicit UnionFind(size_t n) : rank_(n, 0)
//    , size_(n, 1)
    {
        parent_.resize(n);
        while (n--) parent_[n] = n;
    }

//    int size(int x) {
//        return size_[x];
//    }

    int find(int x) {
        if (parent_[x] == x) return x;
        return parent_[x] = find(parent_[x]); // path compression
    }

    void merge(int x, int y) {
        x = find(x), y = find(y);
        if (rank_[x] > rank_[y]) {
            parent_[y] = x;
//            size_[x] = size_[x] + size_[y];
            return;
        }

        parent_[x] = y;
//        size_[y] = size_[x] + size_[y];
        if (rank_[x] == rank_[y]) {
            rank_[x]++;
        }
    }

    bool equiv(int x, int y) {
        return find(x) == find(y);
    }

private:
    vector<int> parent_;
    vector<int> rank_;
//    vector<int> size_;
};

inline double squaredDistance(const pair<double, double> &t1, const pair<double, double> &t2) {
    double dx = t1.first - t2.first;
    double dy = t1.second - t2.second;
    return dx * dx + dy * dy;
}

template<class T=int, class U=int>
U tourDistance(Grid<U> &distances, vector<T> tour) {
    U distance = 0;
    for (T i = 0; i < (T) tour.size(); i++) {
        distance += distances[tour[i]][tour[(i + 1) % distances.rows()]];
    }
    return distance;
}

vector<pair<double, double>> createNCities(int n, int seed) {
    double lowerBound = -1e6;
    double upperBound = 1e6;
    vector<pair<double, double >> cities = vector<pair<double, double >>();
    cities.reserve(n);

    uniform_real_distribution<double> unif(lowerBound, upperBound);
    default_random_engine generator(seed);
    for (int i = 0; i < n; i++) {
        cities.emplace_back(pair<double, double>(unif(generator), unif(generator)));
    }
    return cities;
}

template<class U=int>
inline Grid<U> *distanceMatrix(const vector<pair<double, double>> &cities) {
    auto n = cities.size();
    auto matrix = new Grid<U>(n, n);
    for (unsigned i = 0; i < n; i++) {
        (*matrix)[i][i] = 0;
    }
    if (is_integral<U>::value) {
        for (unsigned i = 0; i < n - 1; i++) {
            for (unsigned j = i + 1; j < n; j++) {
                (*matrix)[i][j] = (*matrix)[j][i] = round(sqrt(squaredDistance(cities[i], cities[j])));
            }
        }
    } else {
        for (unsigned i = 0; i < n - 1; i++) {
            for (unsigned j = i + 1; j < n; j++) {
                (*matrix)[i][j] = (*matrix)[j][i] = sqrt(squaredDistance(cities[i], cities[j]));
            }
        }
    }

    return matrix;
}

template<class T, class U>
inline Grid<T> *kNearestNeighbors(Grid<U> &distances, int k) {
    auto n = distances.rows();
    auto knn = new Grid<T>(n, k);

    auto cmp = [](const pair<U, T> &a, const pair<U, T> &b) { return a.first > b.first; };
    auto heap = vector<pair<U, T>>();
    heap.reserve(n);

    for (T i = 0; i < n; i++) {
        heap.clear();
        for (T j = 0; j < n; j++) {
            if (i != j) { heap.emplace_back(distances[i][j], j); }
        }
        make_heap(heap.begin(), heap.end(), cmp);

        for (T neighbor = 0; neighbor < k; neighbor++) {
            (*knn)[i][neighbor] = heap.front().second;
            pop_heap(heap.begin(), heap.end(), cmp);
            heap.pop_back();
        }
    }
    return knn;
}

// reverse tour segment if it makes the tour shorter given three vertices (each vertex is used for its left edge)
template<class T, class U>
inline int reverseSegment3OptSeq(vector<T> *tour, int i, int j, int k, Grid<U> &distances, bool) {
    int a = (*tour)[(i - 1 + tour->size()) % tour->size()];
    int b = (*tour)[i];
    int c = (*tour)[j - 1];
    int d = (*tour)[j];
    int e = (*tour)[k - 1];
    int f = (*tour)[k % (*tour).size()];

    // original distance
    int d0 = distances[a][b] + distances[c][d] + distances[e][f];


    int d1 = distances[a][c] + distances[b][d] + distances[e][f];
    if (d0 > d1) {
        reverse(tour->begin() + i, tour->begin() + j);
        return d0 - d1;
    }

    int d2 = distances[a][b] + distances[c][e] + distances[d][f];
    if (d0 > d2) {
        reverse(tour->begin() + j, tour->begin() + k);
        return d0 - d2;
    }

    int d4 = distances[f][b] + distances[c][d] + distances[e][a];
    if (d0 > d4) {
        reverse(tour->begin() + i, tour->begin() + k);
        return d0 - d2;
    }

    int d3 = distances[a][d] + distances[e][b] + distances[c][f];
    if (d0 > d3) {
        vector<int> tempTour = vector<int>{};
        tempTour.insert(tempTour.end(), tour->begin() + j, tour->begin() + k);
        tempTour.insert(tempTour.end(), tour->begin() + i, tour->begin() + j);
        copy_n(tempTour.begin(), tempTour.size(), &(*tour)[i]);
        return d0 - d3;
    }

    int d5 = distances[a][e] + distances[d][b] + distances[c][f];
    if (d0 > d5) {
        // get ac bd ef like in d1
        reverse(tour->begin() + i, tour->begin() + j);
        // get ae db cf
        reverse(tour->begin() + i, tour->begin() + k);
        return d0 - d5;
    }

    int d6 = distances[a][c] + distances[b][e] + distances[d][f];
    if (d0 > d6) {
        // get ac bd ef like with d1
        reverse(tour->begin() + i, tour->begin() + j);
        // now reverse from d to e
        reverse(tour->begin() + j, tour->begin() + k);
        return d0 - d6;
    }

    int d7 = distances[a][d] + distances[e][c] + distances[b][f];
    if (d0 > d7) {
        // get ab ce df like with d2
        reverse(tour->begin() + j, tour->begin() + k);
        // now reverse from d to b
        reverse(tour->begin() + i, tour->begin() + k);
        return d0 - d7;
    }

    return 0;
}

void logCities(const vector<pair<double, double>> &cities, const string &path, const string &name) {
    ofstream outfile(path, ios_base::app);
    outfile << "\ncities_name=" << name << "\n";
    outfile << "cities=[";
    for (const auto &city : cities) {
        outfile << "(" << city.first << "," << city.second << "),";
    }
    outfile << "]\n";
    outfile.close();
}

template<class T>
void logTour(vector<T> tour, const string &path, const string &name) {
    ofstream outfile(path, ios_base::app);
    outfile << "tour_name=" << name << "\n";
    outfile << "tour=[";
    for (auto c: tour) {
        outfile << c << ",";
    }
    outfile << "]\n";
    outfile.close();
}

template<class T>
inline void createCityToTourIdx(const vector<T> &tour, vector<T> &cityToTourIdx) {
    for (T tourIdx = 0; tourIdx < (T) tour.size(); tourIdx++) {
        cityToTourIdx[tour[tourIdx]] = tourIdx;
    }
}

#endif //KTH_AA_UTILITY_H