#ifndef KTH_AA_CHRISTOFIDES
#define KTH_AA_CHRISTOFIDES

template<class T, class U>
void kruskalsMst(Grid<U> &distances, unsigned int n, vector<tuple<T, T, U>> &selectedEdges,
                 vector<T> &cityMstDegree);

template<class T, class U>
void minimumWeightPerfectMatching(Grid<U> &distances, const vector<T> &unmatched,
                                  vector<tuple<T, T, U>> &matchedEdges);

template<class T, class U>
void findEulerianCircuit(unsigned int n, const vector<tuple<T, T, U>> &selectedEdges,
                         const vector<tuple<T, T, U>> &matchedEdges, vector<T> &eulerianCircuit);

#include <vector>
#include <tuple>
#include <stack>
#include <unordered_map>
#include "utility.h"
#include "blossom5_all_in_one_file.h"

using namespace std;

template<class T, class U>
vector<T> travelChristofides(Grid<U> &distances) {
    unsigned n = distances.rows();

    // *** Kruskal's algo ***
    auto selectedEdges = vector<tuple<T, T, U>>(); // O(1)
    auto cityMstDegree = vector<T>(n, 0); // O(1)
    kruskalsMst(distances, n, selectedEdges, cityMstDegree); // O(TODO BLOSSOM V)

    // *** Minimum-weight perfect matching ***
    // O(n)
    auto unmatched = vector<T>();
    unmatched.reserve(n);
    for (unsigned i = 0; i < n; i++) {
        if (cityMstDegree[i] % 2) unmatched.emplace_back(i);
    }

    // O(1)
    auto matchedEdges = vector<tuple<T, T, U>>();
    matchedEdges.reserve(unmatched.size() / 2);

    // O(TODO BLOSSOM V)
    minimumWeightPerfectMatching(distances, unmatched, matchedEdges);


    // *** Find Euler ***
    // O(1)
    auto eulerianCircuit = vector<T>();
    eulerianCircuit.reserve(10 * n);

    // O(n)
    findEulerianCircuit(n, selectedEdges, matchedEdges, eulerianCircuit);


    // *** Find Hamiltonian ***

    // O(TODO)
    auto visitedCities = vector<bool>(n, false);
    auto hamiltonian = vector<T>();
    for (auto city: eulerianCircuit) {
        if (!visitedCities[city]) {
            visitedCities[city] = true;
            hamiltonian.emplace_back(city);
        }
    }

    return hamiltonian;
}

template<class T, class U>
void kruskalsMst(Grid<U> &distances, unsigned int n, vector<tuple<T, T, U>> &selectedEdges,
                 vector<T> &cityMstDegree) {
    auto edgesHeap = vector<tuple<T, T, U>>(); // i, j, distance
    edgesHeap.reserve(n * (n - 1) / 2);

    // O(n^2)
    for (unsigned i = 0; i < n - 1; i++) {
        for (unsigned j = i + 1; j < n; j++) {
            edgesHeap.emplace_back(make_tuple(i, j, distances[i][j]));
        }
    }

    // O(n^2)
    auto cmp = [](const tuple<T, T, U> &a, const tuple<T, T, U> &b) { return get<2>(a) > get<2>(b); };
    make_heap(edgesHeap.begin(), edgesHeap.end(), cmp);

    // O(n)
    UnionFind uf = UnionFind(n);

    // O(1)
    selectedEdges.reserve(n - 1); // there are exactly n-1 edges in a MSTree of a fully connected graph

    // O(1)
    // O(n^2 * (log(n^2) + TODO))
    while (selectedEdges.size() < n - 1) {
        // O(log(n**2))
        auto nextEdge = edgesHeap.front();
        pop_heap(edgesHeap.begin(), edgesHeap.end(), cmp);
        edgesHeap.pop_back();

        if (!uf.equiv(get<0>(nextEdge), get<1>(nextEdge))) {
            // O(TODO)
            uf.merge(get<0>(nextEdge), get<1>(nextEdge));
            selectedEdges.emplace_back(nextEdge);
            cityMstDegree[get<0>(nextEdge)]++;
            cityMstDegree[get<1>(nextEdge)]++;
        }
    }
}

template<class T, class U>
void minimumWeightPerfectMatching(Grid<U> &distances, const vector<T> &unmatched,
                                  vector<tuple<T, T, U>> &matchedEdges) {

    // O(BLOSSOM_V) TODO
    struct PerfectMatching::Options options;
    options.verbose = false;
    int pmEdges = unmatched.size() * (unmatched.size() - 1) / 2;
    int pmNodes = unmatched.size();

    vector<PerfectMatching::EdgeId> edgeIndices;
    auto *pm = new PerfectMatching(pmNodes, pmEdges);
    for (unsigned i = 0; i < unmatched.size() - 1; i++) {
        for (unsigned j = i + 1; j < unmatched.size(); j++) {
            edgeIndices.emplace_back(pm->AddEdge(i, j, distances[unmatched[i]][unmatched[j]]));
        }
    }

    pm->options = options;
    pm->Solve();
    int idx = 0;
    for (unsigned i = 0; i < unmatched.size() - 1; i++) {
        for (unsigned j = i + 1; j < unmatched.size(); j++) {
            if (pm->GetSolution(edgeIndices[idx++])) {
                matchedEdges.emplace_back(
                        tuple<T, T, U>(unmatched[i], unmatched[j], distances[unmatched[i]][unmatched[j]]));
            }
        }
    }
    delete pm;
}

template<class T, class U>
void findEulerianCircuit(unsigned int n, const vector<tuple<T, T, U>> &selectedEdges,
                         const vector<tuple<T, T, U>> &matchedEdges, vector<T> &eulerianCircuit) {

    // O(n)
    auto graph = vector<unordered_map<T, tuple<T, T, U>>>(n); // dont judge me :)
    for (auto um: graph) um.reserve(n);

    // O(n)
    for (auto e: selectedEdges) {
        graph[get<0>(e)][get<1>(e)] = e;
        graph[get<1>(e)][get<0>(e)] = e;
    }

    // O(n)
    for (auto e: matchedEdges) {
        graph[get<0>(e)][get<1>(e)] = e;
        graph[get<1>(e)][get<0>(e)] = e;
    }

    auto st = stack<T>();
    st.push(0);
    eulerianCircuit.emplace_back(0);

    // O(TODO)
    while (!st.empty()) {
        auto city = st.top();
        if (graph[city].empty()) {
            st.pop();
            eulerianCircuit.emplace_back(city);
            continue;
        } else {
            auto entry = graph[city].begin();
            auto nextCity = entry->first;
            graph[city].erase(entry);
            graph[nextCity].erase(city);
            st.push(nextCity);
        }
    }
}

#endif //KTH_AA_CHRISTOFIDES