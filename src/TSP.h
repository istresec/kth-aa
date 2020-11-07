#ifndef KTH_AA_TSP_H
#define KTH_AA_TSP_H

#include <iostream>
#include <vector>
#include <stack>
#include <chrono>
#include <blossom5-v2.05.src/PerfectMatching.h>
#include <blossom5-v2.05.src/GEOM/GeomPerfectMatching.h>

#include "utility.h"

using namespace std;
using namespace chrono;

class TSP {
public:
    static vector<int> travel(const vector<pair<double, double>> &cities);

    static vector<int> travel_naive(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> travel_nn(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> travel_bruteforce(const vector<pair<double, double>> &cities, Grid<int> &distances);

    static vector<int> travel_cw(const vector<pair<double, double>> &cities, Grid<int> &distances, int hub = 0);

    static vector<int> local_2opt(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knearest,
                                  time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

    static vector<int> local_2opt_no_knn(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                                         time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

    static vector<int> local_3opt(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                                  time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

    static vector<int> local_3opt_no_knn(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                                         time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

    static vector<int> local_3opt_no_knn_sequential(vector<int> tour, Grid<int> &distances, Grid<uint16_t> &knn,
                                                    time_point<steady_clock, duration<long long int, ratio<1, 1000000000>>> *deadline);

    template<class T=uint16_t, class U=int>
    static vector<int> travel_christofides(Grid<U> &distances);
};

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

    bool pm = true;
    if (pm) {
        // O(BLOSSOM_V) TODO
        struct PerfectMatching::Options options;
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
        auto matched_edges = vector<tuple<T, T, U>>();
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
        if (check_perfect_matching) {
            int res = CheckPerfectMatchingOptimality(node_num, edge_num, edges, weights, pm);
            printf("check optimality: res=%d (%s)\n", res, (res == 0) ? "ok" : ((res == 1) ? "error" : "fatal error"));
        }
        double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
        printf("cost = %.1f\n", cost);
        if (save_filename) SaveMatching(node_num, pm, save_filename);
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
        auto matched_edges = vector<tuple<T, T, U>>();
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

#endif //KTH_AA_TSP_H
