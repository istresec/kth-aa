#ifndef KTH_AA_CHOKOLINO
#define KTH_AA_CHOKOLINO

#include <chrono>
#include <vector>
#include <tuple>
#include <stack>
#include <unordered_set>
#include "utility.h"

using namespace std;
using namespace chrono;

// insures that edge defined by two towns is unique (the first in the pair is always smaller)
template<class T>
pair<T, T> create_edge_pair(T t1, T t2) {
    if (t1 < t2) return pair(t1, t2);
    return pair(t2, t1);
}

template<class T, class U>
bool choose_y(vector<T> &tour, unordered_set<pair<T, T>, pair_hash> &tour_edges, vector<T> &city_to_tour_idx,
              Grid<U> &distances, Grid<T> &knn, T tour_1, T tour_2i, U gain,
              unordered_set<pair<T, T>, pair_hash> &broken, unordered_set<pair<T, T>, pair_hash> &joined,
              time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

template<class T>
inline unordered_set<pair<T, T>, pair_hash> get_tour_edges(const vector<T> &tour) {
    unordered_set<pair<T, T>, pair_hash> tour_edges;
    tour_edges.reserve(tour.size());
    for (T i = 0; i < tour.size(); i++) {
        tour_edges.emplace(create_edge_pair(tour[i], tour[(i + 1) % tour.size()]));
    }
    return tour_edges;
}

template<class T>
inline bool generate_tour(const vector<T> &tour,
                          const unordered_set<pair<T, T>, pair_hash> &added,
                          const unordered_set<pair<T, T>, pair_hash> &removed,
                          vector<T> &new_tour) {
    T n = tour.size();
    auto city_to_neighbors = Grid<T>(n, 2); // city --> [pred, successor] initially
    for (T i = 0; i < n; i++) {
        city_to_neighbors[tour[i]][0] = tour[i == 0 ? n - 1 : i - 1];
        city_to_neighbors[tour[i]][1] = tour[i == n - 1 ? 0 : i + 1];
    }

    for (auto edge_to_remove: removed) {
        if (city_to_neighbors[edge_to_remove.first][0] == edge_to_remove.second) {
            city_to_neighbors[edge_to_remove.first][0] = n;
            city_to_neighbors[edge_to_remove.second][1] = n;
        } else {
            city_to_neighbors[edge_to_remove.first][1] = n;
            city_to_neighbors[edge_to_remove.second][0] = n;
        }
    }

    for (auto edge_to_add: added) {
        if (city_to_neighbors[edge_to_add.first][0] == n) {
            city_to_neighbors[edge_to_add.first][0] = edge_to_add.second;
        } else {
            city_to_neighbors[edge_to_add.first][1] = edge_to_add.second;
        }
        if (city_to_neighbors[edge_to_add.second][0] == n) {
            city_to_neighbors[edge_to_add.second][0] = edge_to_add.first;
        } else {
            city_to_neighbors[edge_to_add.second][1] = edge_to_add.first;
        }
    }

    new_tour.reserve(n);
    T start_city = 0;
    T prev_city;
    T next_city = start_city;
    new_tour.emplace_back(next_city);
    for (T i = 0; i < n - 1; i++) {
        auto left = city_to_neighbors[next_city][0];
        auto right = city_to_neighbors[next_city][1];
        if (right == prev_city) {
            prev_city = next_city;
            next_city = left;
            new_tour.emplace_back(next_city);
        } else {
            prev_city = next_city;
            next_city = right;
            new_tour.emplace_back(next_city);
        }

        if (next_city == start_city) {
            return false;
        }
    }

    if (city_to_neighbors[next_city][0] != start_city and city_to_neighbors[next_city][1] != start_city) {
        return false;
    }

    return true;
}

// chooses next edge to remove, part of the Lin-Kernighan algorithm implementation
// implementation based on https://arthur.maheo.net/implementing-lin-kernighan-in-python/
template<class T, class U>
bool choose_x(vector<T> &tour, unordered_set<pair<T, T>, pair_hash> &tour_edges, vector<T> &city_to_tour_idx,
              Grid<U> &distances, Grid<T> &knn, T tour_1, T tour_last, U gain,
              unordered_set<pair<T, T>, pair_hash> &broken, unordered_set<pair<T, T>, pair_hash> &joined,
              time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {

    // neighbourhood of last edge in our tour (indices in tour)
    pair<T, T> neighbourhood = make_pair((T) (tour_last + 1) % tour.size(),
                                         (T) ((tour_last + tour.size() - 1) % tour.size()));

    // special case for x4
    if (broken.size() == 4) {
        // give choice priority to longer edge
        if (distances[tour[neighbourhood.first]][tour_last] > distances[tour[neighbourhood.first]][tour_last]) {
            swap(neighbourhood.first, neighbourhood.second);
        }
    }

    T city_1 = tour[tour_1];
    for (uint8_t n_idx = 0; n_idx <= 1; n_idx++) {
        T tour_2i = n_idx == 0 ? neighbourhood.first : neighbourhood.second;
        T city_2i = tour[tour_2i];
        auto xi = create_edge_pair(tour[tour_last], city_2i);

        U gain_i = gain + distances[tour[tour_last]][city_2i];

        // verify that X and Y are disjoint and that xi is not already in there...
        if (joined.find(xi) == joined.end() and broken.find(xi) == broken.end()) {
            auto removed = broken;
            removed.emplace(xi);

            auto y_relink = create_edge_pair(city_2i, city_1);
            vector<T> new_tour;
            bool is_tour;
            U relink = gain_i - distances[city_2i][city_1];
            if (relink > 0 and joined.find(y_relink) == joined.end() and broken.find(y_relink) == broken.end()) {
                auto added = joined;
                added.emplace(y_relink);
                is_tour = generate_tour(tour, added, removed, new_tour);
            } else {
                is_tour = false;
            }

            if (is_tour) {
                tour = new_tour;
                tour_edges = get_tour_edges(tour);
                create_city_to_tour_idx(tour, city_to_tour_idx);
                return true;
            } else {
                if (joined.size() > 3) continue;
                return choose_y(tour, tour_edges, city_to_tour_idx, distances, knn, tour_1, tour_2i, gain_i, removed,
                                joined, deadline);
            }
        }
    }

    return false;
}

// chooses next edge to add, part of the Lin-Kernighan algorithm implementation
// implementation based on https://arthur.maheo.net/implementing-lin-kernighan-in-python/
template<class T, class U>
bool choose_y(vector<T> &tour, unordered_set<pair<T, T>, pair_hash> &tour_edges, vector<T> &city_to_tour_idx,
              Grid<U> &distances, Grid<T> &knn, T tour_1, T tour_2i, U gain,
              unordered_set<pair<T, T>, pair_hash> &broken, unordered_set<pair<T, T>, pair_hash> &joined,
              time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {

    T city_2i = tour[tour_2i];

    // original heuristic: check only 5 closest neighbours when doing y2, otherwise just 1
    // int top = broken.size() == 2 ? 5 : 1;
    // int top = (broken.size() == 2) ? 30 : ((broken.size() == 3) ? 5 : 5);
    int top = (broken.size() == 2) ? 20 : ((broken.size() == 3) ? 5 : ((broken.size() == 4) ? 3 : 1)); // tuning->kattis


    // try to find an edge to add from closest neighbours
    for (unsigned i = 0; i < knn.columns() and top >= 1; i++) {
        T city_2i_plus_1 = knn[city_2i][i];
        auto yi = create_edge_pair(city_2i, city_2i_plus_1);

        // Three checks need to be done
        // 1. Gain still positive
        U gain_i = gain - distances[city_2i][city_2i_plus_1];
        if (gain_i <= 0) {
            continue;
        }
        // 2. yi not already removed or added to tour
        if (joined.find(yi) != joined.end() or broken.find(yi) != broken.end()) {
            continue;
        }
        // 3. yi is not in tour_edges
        if (tour_edges.find(yi) != tour_edges.end()) {
            continue;
        }

        top--;
        auto added = joined;
        added.emplace(yi);
        if (choose_x(tour, tour_edges, city_to_tour_idx, distances, knn, tour_1, city_to_tour_idx[city_2i_plus_1],
                     gain_i, broken, added, deadline)) {
            return true;
        }
    }

    return false;
}

// main loop for the Lin-Kernighan algorithm
// implementation based on the outline of https://arthur.maheo.net/implementing-lin-kernighan-in-python/
template<class T, class U>
bool lin_kernighan_main_loop(vector<T> &tour, unordered_set<pair<T, T>, pair_hash> &tour_edges,
                             vector<T> &city_to_tour_idx, Grid<U> &distances, Grid<T> &knn,
                             time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    // steps refer to a basic outline of the LK algorithm found on page 12 of
    // "An Effective Implementation of the Lin-Kernighan Traveling Salesman Heuristic"
    // by Keld Helsgaun (Roskilde University, Denmark)
    // http://akira.ruc.dk/~keld/research/LKH/LKH-2.0/DOC/LKH_REPORT.pdf
    // TODO (they might be a description from the original work by Lin & Kernighan?)

    auto set_x = unordered_set<pair<T, T>, pair_hash>(); // set of broken edges
    auto set_y = unordered_set<pair<T, T>, pair_hash>(); // set of added edges

    bool any_improvement = false;
    // step 2
    //
    // (the following loops basically look for a viable 2-opt solution, but then go on a bit further...
    for (T tour_i = 0; tour_i < (T) tour.size(); tour_i++) { // go through all possible towns for t1
        T city_1 = tour[tour_i];

        if (deadline != nullptr and (*deadline <= steady_clock::now())) {
            return false;
        }

        bool inner_improvement = false;
        // step 3
        //
        for (uint8_t n_idx = 0; n_idx <= 1; n_idx++) { // get towns around t1 for t2
            T tour_j = (tour_i + tour.size() + (n_idx == 0 ? -1 : 1)) % tour.size();
            if (inner_improvement) break;
            T city_2 = tour[tour_j];

            // create edge x1 (edge between t1 and t2) and add it to set of edges X
            auto x1 = create_edge_pair(city_1, city_2);
            set_x = unordered_set<pair<T, T>, pair_hash>(); // or set_x.clear();
            set_x.emplace(x1);

            // step 4
            //
            // find t3 to create a second edge for 2-opt
            for (T neighbor_idx = 0; neighbor_idx < (T) knn.columns(); neighbor_idx++) {
                // TODO consider using some better metric than knn
                //      (like sorting neighbors by 1-tree score or by possible relinking gain)
                auto city_3 = knn[city_2][neighbor_idx];
                auto y1 = create_edge_pair(city_2, city_3);
                U gain = distances[city_1][city_2] - distances[city_2][city_3];

                if (gain > 0) { // if the two edges would create a shorter tour go find more
                    // y1 is not in tour_edges, if it is continue
                    if (tour_edges.count(y1) > 0) {
                        continue;
                    }
                    // step 6
                    //
                    T tour_k = city_to_tour_idx[city_3];
                    set_y = unordered_set<pair<T, T>, pair_hash>(); // or set_y.clear()
                    set_y.emplace(y1);
                    if (choose_x(tour, tour_edges, city_to_tour_idx, distances, knn, tour_i, tour_k, gain, set_x,
                                 set_y, deadline)) {
                        inner_improvement = true;
                        any_improvement = true;
                        break;
                    }
                }
            }
        }
    }

    return any_improvement;
}

// implementation based on https://arthur.maheo.net/implementing-lin-kernighan-in-python/
template<class T, class U>
static vector<T> chokolino(vector<T> &tour, Grid<U> &distances, Grid<T> &knn,
                           time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    bool improved = true;
    unordered_set<pair<T, T>, pair_hash> tour_edges = get_tour_edges(tour);
    vector<T> city_to_tour_idx(tour.size());
    city_to_tour_idx.reserve(tour.size());
    create_city_to_tour_idx(tour, city_to_tour_idx);
    while ((deadline == nullptr or (*deadline > steady_clock::now())) and improved) {
        improved = lin_kernighan_main_loop(tour, tour_edges, city_to_tour_idx, distances, knn, deadline);
    }
    return tour;
}

#endif // KTH_AA_CHOKOLINO