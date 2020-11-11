#ifndef KTH_AA_CHRISTOFIDES
#define KTH_AA_CHRISTOFIDES

#include <vector>
#include <tuple>
#include <stack>
#include "utility.h"
#include "blossom5_all_in_one_file.h"

using namespace std;

template<class T, class U>
vector<T> travel_christofides(Grid<U> &distances) {
    unsigned n = distances.rows();

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
    auto unmatched = vector<T>();
    unmatched.reserve(n);
    for (unsigned i = 0; i < n; i++) {
        if (city_mst_degree[i] % 2) unmatched.emplace_back(i);
    }
//    if (unmatched.size() % 2) {
//        throw logic_error("Inconsistent state: MST does not have an even number of odd degree nodes!");
//    }

    auto matched_edges = vector<tuple<T, T, U>>();
    // O(BLOSSOM_V) TODO
    struct PerfectMatching::Options options;
    options.verbose = false;
    struct GeomPerfectMatching::GPMOptions gpm_options;
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
    delete pm;

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
    auto hamiltonian = vector<T>();
    for (auto city: eulerian_circuit) {
        if (!visited_cities[city]) {
            visited_cities[city] = true;
            hamiltonian.emplace_back(city);
        }
    }

    return hamiltonian;
}

#endif //KTH_AA_CHRISTOFIDES