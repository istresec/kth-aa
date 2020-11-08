#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>

#include "TSP.h"
//#include <blossom5-v2.05.src/PerfectMatching.h>
//#include <blossom5-v2.05.src/GEOM/GeomPerfectMatching.h>
#include "blossom.h"


using namespace std;
using namespace chrono;

vector<int> TSP::travel(const vector<pair<double, double>> &cities) {
    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> deadline =
            steady_clock::now() + milliseconds(1980);
    auto distances = distance_matrix(cities);
    vector<int> tour;
    if (cities.size() <= 2) {
        tour = vector<int>();
        for (unsigned i = 0; i < cities.size(); i++) tour.emplace_back();
    } else if (cities.size() <= 13) {
        tour = travel_bruteforce(cities, *distances);
    } else {
//        tour = travel_cw(cities, *distances);
        tour = travel_christofides(*distances);
//    Grid<uint16_t> *knn = nullptr;
//    tour = TSP::local_2opt_no_knn(tour, *distances, *knn, &deadline);
        Grid<uint16_t> *knn = k_nearest_neighbors<uint16_t>(*distances, min((int) (cities.size() - 1), 50));
        tour = TSP::local_2opt(tour, *distances, *knn, &deadline);
        tour = TSP::local_3opt_no_knn_sequential(tour, *distances, *knn, &deadline);
    }
    return tour;
}

inline bool compare_savings(tuple<int, int, int> t1, tuple<int, int, int> t2) {
    return get<2>(t1) > get<2>(t2);
}

vector<int> TSP::travel_bruteforce(const vector<pair<double, double>> &cities, Grid<int> &distances) {
    auto city_indices = vector<int>();
    city_indices.reserve(cities.size());
    for (unsigned i = 0; i < cities.size(); i++) city_indices.emplace_back(i);

    auto best = city_indices;
    int best_distance = -1;
    do {
        int current_distance = tour_distance(distances, city_indices);
        if (best_distance == -1 or current_distance < best_distance) {
            best_distance = current_distance;
            best = city_indices;
        }
    } while (next_permutation(city_indices.begin() + 1, city_indices.end()));

    return best;
}

vector<tuple<int, int, int>> savings(const vector<pair<double, double>> &cities, Grid<int> &distances, int hub) {
    vector<tuple<int, int, int >> savings = vector<tuple<int, int, int >>();
    for (unsigned i = 0; i < cities.size() - 1; i++) {
        if (i == hub) continue;
        for (unsigned j = i + 1; j < cities.size(); j++) {
            savings.emplace_back(make_tuple(i, j, distances[hub][i] + distances[hub][j] - distances[i][j]));
        }
    }
    sort(savings.begin(), savings.end(), compare_savings);
    return savings;
}

vector<int> TSP::travel_naive(const vector<pair<double, double>> &cities, Grid<int> &distances) {
    vector<int> tour = vector<int>(cities.size());
    vector<bool> used = vector<bool>(cities.size());
    tour[0] = 0;
    used[0] = true;
    for (unsigned i = 1; i <= cities.size() - 1; i++) {
        int best = -1;
        for (unsigned j = 0; j <= cities.size() - 1; j++) {
            bool is_better = best == -1 or distances[i - 1][j] < distances[i - 1][best];
            if (not used[j] and is_better) {
                best = (int) j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
    return tour;
}

// Nearest Neighbour Algorithm - pick random vertex and go the shortest distance you can
vector<int> TSP::travel_nn(const vector<pair<double, double>> &cities, Grid<int> &distances) {
    vector<int> tour = vector<int>();
    return tour;
}

// Clarke-Wright Savings Algorithm. City at index 0 used as hub.
vector<int> TSP::travel_cw(const vector<pair<double, double>> &cities, Grid<int> &distances, int hub) {
    // if only one city
    if (cities.size() == 1)
        return vector<int>{0};

    vector<tuple<int, int, int>> s = savings(cities, distances, hub);

    // initialize tours
    vector<int> tours[cities.size()];
    vector<int> empty_tour = vector<int>();
    for (int i = 0; i < (int) cities.size(); i++) {
        if (i == hub) {
            tours[i] = empty_tour;
        } else {
            tours[i] = vector<int>{i}; // instead of 0, i, 0 just use i
        }
    }

    // algorithm
    vector<int> temp_tour;
    vector<int> first, second;
    int i, j;
    for (auto &it : s) {
        i = get<0>(it);
        j = get<1>(it);
        // if two distinct tours with endpoints i and j exist, if so combine them (O(1))
        if (not tours[i].empty() and not tours[j].empty() and tours[i].front() != tours[j].front()) {
            first = tours[i]; // remember tour with endpoint i
            second = tours[j]; // remember tour with endpoint j
            // remove tours with endpoints i and j while a new tour is constructed
            tours[first.front()] = tours[first.back()] = tours[second.front()] = tours[second.back()] = empty_tour;

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
    vector<int> tour = vector<int>();
    for (i = 1; i < (int) cities.size(); i++) {
        if (not tours[i].empty())
            break;
    }

    tour.emplace_back(hub);
    tour.insert(tour.end(), tours[i].begin(), tours[i].end());

    return tour;
}

template<class T, class U>
inline int local2OptChange(Grid<T> &tour, Grid<U> &distances, int i, int j) {
    return (distances[i][tour[i][1]] + distances[j][tour[j][1]]) -
           (distances[i][j] + distances[tour[i][1]][tour[j][1]]);
}

inline void local2OptUpdate(Grid<uint16_t> &tour, uint16_t i, uint16_t j) {
    uint16_t ii = tour[i][1], jj = tour[j][1];
    tour[i][1] = j;
    tour[jj][0] = ii;
    uint16_t next = tour[ii][1];
    while (next != j) {
        swap(tour[next][0], tour[next][1]);
        next = tour[next][0];
    }
    tour[j][1] = tour[j][0];
    tour[j][0] = i;
    tour[ii][0] = tour[ii][1];
    tour[ii][1] = jj;
}

vector<int> TSP::local_2opt(vector<int> tour_vector, Grid<int> &distances, Grid<uint16_t> &knn,
                            time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    // TODO use a tree structure for storing tour so that reverse is fast
    bool sequential = false;

    uint16_t n = tour_vector.size();
    auto tour = Grid<uint16_t>(n, 2);
    for (int i = 0; i < n; i++) {
        tour[tour_vector[i]][0] = tour_vector[i == 0 ? i - 1 + n : i - 1];
        tour[tour_vector[i]][1] = tour_vector[i == n - 1 ? i + 1 - n : i + 1];
    }

    int bestChange;
    auto shouldBeChecked = vector<bool>(n, true);
    do {
        bestChange = 0.;
        uint16_t best_i, best_j;
        for (uint16_t i = 0; i < n - 1; i++) {
//            if (not shouldBeChecked[i]) continue;

//            bool iCannotDoBetter = true;
            for (uint16_t j = 0; j < knn.columns(); j++) {
                uint16_t a = i, b = knn[i][j];
                if (i > b) continue;
                if (tour[i][0] == knn[i][j] or tour[i][1] == knn[i][j]) continue;

                int change1 = local2OptChange(tour, distances, a, b);
                int change2 = local2OptChange(tour, distances, tour[a][0], tour[b][0]);
                int change = change1 > change2 ? change1 : change2;
                uint16_t better_a = change1 > change2 ? a : tour[a][0];
                uint16_t better_b = change1 > change2 ? b : tour[b][0];

//                if (change > 0) { iCannotDoBetter = false; }
                if (sequential and change > 0) {
                    local2OptUpdate(tour, better_a, better_b);
                }
                if (change > bestChange) {
                    bestChange = change;
                    best_i = better_a;
                    best_j = better_b;
                }
            }
//            if (iCannotDoBetter) { shouldBeChecked[i] = false; }
        }
        if (bestChange > 0) {
//            for (uint16_t j = 0; j < knn.columns(); j++) {
//                shouldBeChecked[knn[best_i][j]] = true;
//                shouldBeChecked[knn[best_j][j]] = true;
//                shouldBeChecked[knn[tour[best_i][1]][j]] = true;
//                shouldBeChecked[knn[tour[best_j][1]][j]] = true;
//            }
            if (not sequential) {
                local2OptUpdate(tour, best_i, best_j);
            }
        }
    } while ((bestChange > 0) & (deadline == nullptr or (steady_clock::now() < *deadline)));

    uint16_t current = 0;
    for (unsigned i = 0; i < n; i++) {
        tour_vector[i] = current;
        current = tour[current][1];
    }
    return tour_vector;
}

vector<int> TSP::local_2opt_no_knn(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knearest,
                                   time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    int best_change;
    do {
        best_change = 0.;
        unsigned best_i, best_j;
        int change;
        for (unsigned i = 0; i < tour.size() - 2; i++) {
            for (unsigned j = i + 2; j < tour.size(); j++) {
                change = (distances[tour[i]][tour[i + 1]] + distances[tour[j]][tour[(j + 1) % tour.size()]]) -
                         (distances[tour[i]][tour[j]] + distances[tour[i + 1]][tour[(j + 1) % tour.size()]]);
                if (change > best_change) {
                    best_change = change;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        if (best_change > 0) {
            reverse(tour.begin() + best_i + 1, tour.begin() + best_j + 1);
        }
    } while ((best_change > 0) and
             (deadline == nullptr or (*deadline > steady_clock::now())));

    return tour;
}

// 3-opt impl. based on pseudocode from https://en.wikipedia.org/wiki/3-opt
vector<int> TSP::local_3opt_no_knn_sequential(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                                              time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    bool better;
    do {
        int n = tour.size();
        better = false;
        for (int i = 0; i < n; i++) {
            if (deadline != nullptr and not(*deadline > steady_clock::now()))
                break;
            for (int j = i + 2; j < n; j++) {
                for (int k = j + 2; k < n + (int) (i > 0); k++) {
                    if (reverse_segment_3opt_seq(&tour, i, j, k, distances, true) > 0)
                        better = true;
                }
            }
        }
    } while (better and (deadline == nullptr or (*deadline > steady_clock::now())));

    return tour;
}

vector<int> TSP::local_3opt_no_knn(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                                   time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    bool better;
    do {
        int n = tour.size();
        better = false;
        int best_diff = 0, best_i, best_j, best_k;
        for (int i = 0; i < n; i++) {
            for (int j = i + 2; j < n; j++) {
                for (int k = j + 2; k < n + (int) (i > 0); k++) {
                    if (deadline != nullptr and not(*deadline > steady_clock::now()))
                        break;
                    int diff = reverse_segment_3opt(&tour, i, j, k, distances, false);
                    if (diff > best_diff) {
                        best_i = i, best_j = j, best_k = k;
                        better = true;
                    }
                }
            }
        }
        if (better) {
            reverse_segment_3opt(&tour, best_i, best_j, best_k, distances, true);
        }
    } while (better and (deadline == nullptr or (*deadline > steady_clock::now())));

    return tour;
}

vector<int> TSP::local_3opt(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                            time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline) {
    bool reverse_seq = true;
    bool apply_seq = true;

    bool better;
    int n = tour.size();
    do {
        better = false;
        int best_diff = 0, best_i, best_j, best_k;

        vector<int> city_to_tour_idx(n);
        for (int tour_idx = 0; tour_idx < n; tour_idx++) {
            city_to_tour_idx[tour[tour_idx]] = tour_idx;
        }

        for (int i = 0; i < n; i++) {
            if (deadline != nullptr and not(*deadline > steady_clock::now()))
                break;
            for (int j_knn = 0; j_knn < knn.columns() - 1; j_knn++) {
                int j = city_to_tour_idx[knn[i][j_knn]];
                if (i >= j - 1) continue;

                for (int k_knn = j_knn + 1; k_knn < knn.columns(); k_knn++) {
                    int k = city_to_tour_idx[knn[i][k_knn]];
                    if (j >= k - 1) continue;

                    for (int successor = 0; successor <= 1; successor++) {
                        int i_current = (i + successor) % n;
                        int j_current = (j + successor) % n;
                        int k_current = (k + successor) % n;
                        if (i_current >= j_current - 1 or j_current >= k_current - 1) continue;

                        int diff = reverse_seq
                                   ? reverse_segment_3opt_seq(&tour, i_current, j_current, k_current, distances,
                                                              apply_seq)
                                   : reverse_segment_3opt(&tour, i_current, j_current, k_current, distances, apply_seq);
                        if (diff > best_diff) {
                            if (not apply_seq) {
                                best_i = i_current, best_j = j_current, best_k = k_current;
                            }
                            better = true;
                        }
                    }
                }
            }
        }

        if (better) {
            if (not apply_seq) {
                reverse_seq
                ? reverse_segment_3opt_seq(&tour, best_i, best_j, best_k, distances, true)
                : reverse_segment_3opt(&tour, best_i, best_j, best_k, distances, true);
            }
        }

    } while (better and (deadline == nullptr or (*deadline > steady_clock::now())));

    return tour;
}


template<class T, class U>
vector<int> TSP::travel_christofides(Grid<U> &distances) {
    int n = distances.rows();

    // *** Kruskal's algo ***

    // O(1)
    auto edges_heap = vector<tuple<T, T, U>>(); // i, j, distance
    edges_heap.reserve(n * (n - 1) / 2);

    // O(n^2)
    for (unsigned i = 0; i < n - 1; i++) {
        for (unsigned j = i + 1; j < n; j++) {
            edges_heap.emplace_back(make_tuple(i, j, distances[i][j]));
        }
    }

    // O(n^2)
    auto cmp = [](const tuple<T, T, U> &a, const tuple<T, T, U> &b) { return get<2>(a) > get<2>(b); };
    make_heap(edges_heap.begin(), edges_heap.end(), cmp);

    // O(n)
    UnionFind uf = UnionFind(n);

    // O(1)
    auto selected_edges = vector<tuple<T, T, U>>();
    selected_edges.reserve(n - 1); // there are exactly n-1 edges in a MSTree of a fully connected graph

    // O(1)
    auto city_mst_degree = vector<T>(n, 0);

    // O(TODO n^2*log(n^2))
    while (selected_edges.size() < n - 1) {
        // O(log(n**2))
        auto next_edge = edges_heap.front();
        pop_heap(edges_heap.begin(), edges_heap.end(), cmp);
        edges_heap.pop_back();

        if (!uf.equiv(get<0>(next_edge), get<1>(next_edge))) {
            // O(TODO)
            uf.merge(get<0>(next_edge), get<1>(next_edge));
            selected_edges.emplace_back(next_edge);
            city_mst_degree[get<0>(next_edge)]++;
            city_mst_degree[get<1>(next_edge)]++;
        }
    }

    // *** Minimum-weight perfect matching ***
    // O(n)
    auto unmatched = vector<int>();
    unmatched.reserve(n);
    for (int i = 0; i < n; i++) {
        if (city_mst_degree[i] % 2) unmatched.emplace_back(i);
    }
    if (unmatched.size() % 2) {
        throw logic_error("Inconsistent state: MST does not have an even number of odd degree nodes!");
    }

    bool use_pm = true;
    auto matched_edges = vector<tuple<T, T, U>>();
    if (use_pm) {
        // O(BLOSSOM_V) TODO
        struct PerfectMatching::Options options;
        options.verbose = false;
        struct GeomPerfectMatching::GPMOptions gpm_options;
        int *edges;
        int *weights;
        int pm_edges = unmatched.size() * (unmatched.size() - 1) / 2;
        int pm_nodes = unmatched.size();

        vector<PerfectMatching::EdgeId> edge_indices;
        auto *pm = new PerfectMatching(pm_nodes, pm_edges);
        for (unsigned i = 0; i < unmatched.size() - 1; i++) {
            for (unsigned j = i + 1; j < unmatched.size(); j++) {
                edge_indices.emplace_back(pm->AddEdge(i, j, distances[unmatched[i]][unmatched[j]]));
            }
        }

        pm->options = options;
        pm->Solve();
//        auto matched_edges = vector<tuple<T, T, U>>();
        matched_edges.reserve(unmatched.size() / 2);
        int idx = 0;
        for (unsigned i = 0; i < unmatched.size() - 1; i++) {
            for (unsigned j = i + 1; j < unmatched.size(); j++) {
                if (pm->GetSolution(edge_indices[idx++])) {
                    matched_edges.emplace_back(
                            tuple<T, T, U>(unmatched[i], unmatched[j], distances[unmatched[i]][unmatched[j]]));
                }
            }
        }
        if (matched_edges.size() != unmatched.size() / 2) {
            throw logic_error("Inconsistent state: Incorrect number of matched edges!");
        }
//        if (check_perfect_matching) {
//            int res = CheckPerfectMatchingOptimality(node_num, edge_num, edges, weights, pm);
//            printf("check optimality: res=%d (%s)\n", res, (res == 0) ? "ok" : ((res == 1) ? "error" : "fatal error"));
//        }
//        double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
//        printf("cost = %.1f\n", cost);
//        if (save_filename) SaveMatching(node_num, pm, save_filename);
        delete pm;


    } else {
        // *** Minimum-weight perfect matching ***

        // O(n)
        unordered_map<int, bool> has_been_matched; // store the cities for which the matching has to be found. maps to true if found
        has_been_matched.reserve(10 * n);
        for (int i = 0; i < n; i++) {
            if (city_mst_degree[i] % 2) has_been_matched[i] = false;
        }
        if (has_been_matched.size() % 2) {
            throw logic_error("Inconsistent state: MST does not have an even number of odd degree nodes!");
        }

        // TODO this is a greedy algorithm. Blossom might be much better (if it is not too slow).
        // O(n^2)
        auto matching_candidates = vector<tuple<T, T, U>>(); // i, j, distance
        matching_candidates.reserve(has_been_matched.size() * (has_been_matched.size() - 1) / 2);

        // O(n^2)
        for (auto city_a: has_been_matched) {
            for (auto city_b: has_been_matched) {
                if (city_a.first >= city_b.first) continue;
                matching_candidates.emplace_back(
                        make_tuple(city_a.first, city_b.first, distances[city_a.first][city_b.first]));
            }
        }

        // O(n^2)
        make_heap(matching_candidates.begin(), matching_candidates.end(), cmp);

        // O(1)
//        auto matched_edges = vector<tuple<T, T, U>>();
        matched_edges.reserve(has_been_matched.size() / 2);

        // O(n^2*log(n^2))
        unsigned matched = 0;
        while (has_been_matched.size() > matched) {
            // O(log(n^2))
            auto next_edge = matching_candidates.front();
            pop_heap(matching_candidates.begin(), matching_candidates.end(), cmp);
            matching_candidates.pop_back();

            // O(1)
            if (!has_been_matched[get<0>(next_edge)] and !has_been_matched[get<1>(next_edge)]) {
                has_been_matched[get<0>(next_edge)] = has_been_matched[get<1>(next_edge)] = true;
                matched += 2;
                matched_edges.emplace_back(next_edge);
            }
        }
    }

    // *** Find Euler ***

    // O(n)
    auto graph = vector<unordered_map<T, tuple<T, T, U>>>(n); // dont judge me :)
    for (auto um: graph) um.reserve(n);

    // O(n)
    for (auto e: selected_edges) {
        graph[get<0>(e)][get<1>(e)] = e;
        graph[get<1>(e)][get<0>(e)] = e;
    }

    // O(n)
    for (auto e: matched_edges) {
        graph[get<0>(e)][get<1>(e)] = e;
        graph[get<1>(e)][get<0>(e)] = e;
    }

    // O(1)
    auto eulerian_circuit = vector<T>();
    eulerian_circuit.reserve(10 * n);
    auto st = stack<T>();
    st.push(0);
    eulerian_circuit.emplace_back(0);

    // O(TODO)
    while (!st.empty()) {
        auto city = st.top();
        if (graph[city].empty()) {
            st.pop();
            eulerian_circuit.emplace_back(city);
            continue;
        } else {
            auto entry = graph[city].begin();
            auto next_city = entry->first;
            graph[city].erase(entry);
            graph[next_city].erase(city);
            st.push(next_city);
        }
    }

    // *** Find Hamiltonian ***

    // O(TODO)
    auto visited_cities = vector<bool>(n, false);
    auto hamiltonian = vector<int>();
    for (auto city: eulerian_circuit) {
        if (!visited_cities[city]) {
            visited_cities[city] = true;
            hamiltonian.emplace_back(city);
        }
    }

    return hamiltonian;
}

