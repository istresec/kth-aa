#ifndef __BLOSSOM__
#define __BLOSSOM__
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef __MINCOST_H__
#define __MINCOST_H__
template<typename FlowType, typename CostType>
class MinCost {
public:
    typedef int NodeId;
    typedef int EdgeId;
    MinCost(int NodeNum, int edgeNumMax, void (*err_function)(const char *) = NULL);
    
    ~MinCost();
    void AddNodeExcess(NodeId i, FlowType excess);
    
    
    
    EdgeId AddEdge(NodeId i, NodeId j, FlowType cap, FlowType rev_cap, CostType cost);
    CostType Solve();
    
    FlowType GetRCap(EdgeId e);
    void SetRCap(EdgeId e, FlowType new_rcap);
    FlowType GetReverseRCap(EdgeId e);
    void SetReverseRCap(EdgeId e, FlowType new_rcap);
    void PushFlow(EdgeId e, FlowType delta);
    void UpdateCost(EdgeId e, FlowType cap_orig, CostType delta);
    CostType GetDual(NodeId i) { return nodes[i].pi; }
protected:
    
    struct Node;
    struct Arc;
    struct Node {
        Arc *firstNonsaturated;
        Arc *firstSaturated;
        Arc *parent;
        Node *next; 
        FlowType excess;
        CostType pi;
        int flag;
        union {
            int heap_ptr;
            Node *next_permanent;
        };
#ifdef MINCOST_DEBUG
        int			id;
#endif
    };
    struct Arc {
        Node *head;
        Arc *prev;
        Arc *next;
        Arc *sister;    
        FlowType r_cap;        
#ifdef MINCOST_DEBUG
        FlowType	cap_orig;
#endif
        CostType cost;
        CostType GetRCost() { return cost + head->pi - sister->head->pi; }
    };
    int nodeNum, edgeNum, edgeNumMax;
    Node *nodes;
    Arc *arcs;
    Node *firstActive;
    int counter;
    CostType cost;
    void (*error_function)(const char *);    
    
    
    
    struct PriorityQueue {
        PriorityQueue();
        ~PriorityQueue();
        void Reset();
        CostType GetKey(Node *i);
        void Add(Node *i, CostType key);
        void DecreaseKey(Node *i, CostType key);
        Node *RemoveMin(CostType &key);
    private:
        struct Item {
            Node *i;
            CostType key;
        } *array;
        int N, arraySize;
        void Swap(int k1, int k2);
    };
    PriorityQueue queue;
    
    void SetRCap(Arc *a, FlowType new_rcap);
    void PushFlow(Arc *a, FlowType delta);
    void Init();
    void DecreaseRCap(Arc *a, FlowType delta);
    void IncreaseRCap(Arc *a, FlowType delta);
    FlowType Augment(Node *start, Node *end);
    void Dijkstra(Node *start);
    void TestOptimality();
#ifdef MINCOST_DEBUG
    void TestCosts();
#endif
};
template<typename CostType>
class DualMinCost : private MinCost<int, CostType> {
public:
    typedef int NodeId;
    DualMinCost(int node_num, int constraint_num_max);
    ~DualMinCost();
    void AddUnaryTerm(NodeId i, int objective_coef);
    void SetLowerBound(NodeId, CostType cmin);
    void SetUpperBound(NodeId, CostType cmax);
    void AddConstraint(NodeId i, NodeId j, CostType cmax); 
    void Solve();
    CostType GetSolution(NodeId i);
private:
    NodeId source;
};
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::AddNodeExcess(NodeId _i, FlowType excess) {
    assert(_i >= 0 && _i < nodeNum);
    nodes[_i].excess += excess;
    if (nodes[_i].excess > 0 && !nodes[_i].next) {
        nodes[_i].next = firstActive;
        firstActive = &nodes[_i];
    }
}
template<typename FlowType, typename CostType>
inline typename MinCost<FlowType, CostType>::EdgeId MinCost<FlowType, CostType>::AddEdge(NodeId _i, NodeId _j, FlowType cap, FlowType rev_cap, CostType cost) {
    assert(_i >= 0 && _i < nodeNum);
    assert(_j >= 0 && _j < nodeNum);
    assert(_i != _j && edgeNum < edgeNumMax);
    assert(cap >= 0);
    assert(rev_cap >= 0);
    Arc *a = &arcs[2 * edgeNum];
    Arc *a_rev = a + 1;
    edgeNum++;
    Node *i = nodes + _i;
    Node *j = nodes + _j;
    a->sister = a_rev;
    a_rev->sister = a;
    if (cap > 0) {
        if (i->firstNonsaturated) i->firstNonsaturated->prev = a;
        a->next = i->firstNonsaturated;
        i->firstNonsaturated = a;
    } else {
        if (i->firstSaturated) i->firstSaturated->prev = a;
        a->next = i->firstSaturated;
        i->firstSaturated = a;
    }
    a->prev = NULL;
    if (rev_cap > 0) {
        if (j->firstNonsaturated) j->firstNonsaturated->prev = a_rev;
        a_rev->next = j->firstNonsaturated;
        j->firstNonsaturated = a_rev;
    } else {
        if (j->firstSaturated) j->firstSaturated->prev = a_rev;
        a_rev->next = j->firstSaturated;
        j->firstSaturated = a_rev;
    }
    a_rev->prev = NULL;
    a->head = j;
    a_rev->head = i;
    a->r_cap = cap;
    a_rev->r_cap = rev_cap;
    a->cost = cost;
    a_rev->cost = -cost;
#ifdef MINCOST_DEBUG
                                                                                                                            a->cap_orig = cap;
	a_rev->cap_orig = rev_cap;
#endif
    if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
    if (a_rev->r_cap > 0 && a_rev->GetRCost() < 0) PushFlow(a_rev, a_rev->r_cap);
    return edgeNum - 1;
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::DecreaseRCap(Arc *a, FlowType delta) {
    a->r_cap -= delta;
    if (a->r_cap == 0) {
        Node *i = a->sister->head;
        if (a->next) a->next->prev = a->prev;
        if (a->prev) a->prev->next = a->next;
        else i->firstNonsaturated = a->next;
        a->next = i->firstSaturated;
        if (a->next) a->next->prev = a;
        a->prev = NULL;
        i->firstSaturated = a;
    }
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::IncreaseRCap(Arc *a, FlowType delta) {
    if (a->r_cap == 0) {
        Node *i = a->sister->head;
        if (a->next) a->next->prev = a->prev;
        if (a->prev) a->prev->next = a->next;
        else i->firstSaturated = a->next;
        a->next = i->firstNonsaturated;
        if (a->next) a->next->prev = a;
        a->prev = NULL;
        i->firstNonsaturated = a;
    }
    a->r_cap += delta;
}
template<typename FlowType, typename CostType>
inline FlowType MinCost<FlowType, CostType>::GetRCap(EdgeId e) {
    Arc *a = &arcs[2 * e];
    return a->r_cap;
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::SetRCap(Arc *a, FlowType new_rcap) {
    assert(new_rcap >= 0);
#ifdef MINCOST_DEBUG
    a->cap_orig += new_rcap - a->r_cap;
#endif
    if (a->r_cap == 0) {
        Node *i = a->sister->head;
        if (a->next) a->next->prev = a->prev;
        if (a->prev) a->prev->next = a->next;
        else i->firstSaturated = a->next;
        a->next = i->firstNonsaturated;
        if (a->next) a->next->prev = a;
        a->prev = NULL;
        i->firstNonsaturated = a;
    }
    a->r_cap = new_rcap;
    if (a->r_cap == 0) {
        Node *i = a->sister->head;
        if (a->next) a->next->prev = a->prev;
        if (a->prev) a->prev->next = a->next;
        else i->firstNonsaturated = a->next;
        a->next = i->firstSaturated;
        if (a->next) a->next->prev = a;
        a->prev = NULL;
        i->firstSaturated = a;
    }
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::SetRCap(EdgeId e, FlowType new_rcap) {
    SetRCap(&arcs[2 * e], new_rcap);
}
template<typename FlowType, typename CostType>
inline FlowType MinCost<FlowType, CostType>::GetReverseRCap(EdgeId e) {
    Arc *a = &arcs[2 * e + 1];
    return a->r_cap;
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::SetReverseRCap(EdgeId e, FlowType new_rcap) {
    SetRCap(&arcs[2 * e + 1], new_rcap);
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::PushFlow(Arc *a, FlowType delta) {
    if (delta < 0) {
        a = a->sister;
        delta = -delta;
    }
    DecreaseRCap(a, delta);
    IncreaseRCap(a->sister, delta);
    a->head->excess += delta;
    a->sister->head->excess -= delta;
    cost += delta * a->cost;
    if (a->head->excess > 0 && !a->head->next) {
        a->head->next = firstActive;
        firstActive = a->head;
    }
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::PushFlow(EdgeId e, FlowType delta) {
    PushFlow(&arcs[2 * e], delta);
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::UpdateCost(EdgeId e, FlowType cap_orig, CostType delta) {
    Arc *a = &arcs[2 * e];
    cost += delta * (cap_orig - a->r_cap);
    a->cost += delta;
    a->sister->cost = -a->cost;
    if (a->GetRCost() > 0) a = a->sister;
    if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
}
template<typename FlowType, typename CostType>
inline MinCost<FlowType, CostType>::PriorityQueue::PriorityQueue() {
    N = 0;
    arraySize = 16;
    array = (Item *) malloc(arraySize * sizeof(Item));
}
template<typename FlowType, typename CostType>
inline MinCost<FlowType, CostType>::PriorityQueue::~PriorityQueue() {
    free(array);
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::PriorityQueue::Reset() {
    N = 0;
}
template<typename FlowType, typename CostType>
inline CostType MinCost<FlowType, CostType>::PriorityQueue::GetKey(Node *i) {
    return array[i->heap_ptr].key;
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::PriorityQueue::Swap(int k1, int k2) {
    Item *a = array + k1;
    Item *b = array + k2;
    a->i->heap_ptr = k2;
    b->i->heap_ptr = k1;
    Node *i = a->i;
    a->i = b->i;
    b->i = i;
    CostType key = a->key;
    a->key = b->key;
    b->key = key;
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::PriorityQueue::Add(Node *i, CostType key) {
    if (N == arraySize) {
        arraySize *= 2;
        array = (Item *) realloc(array, arraySize * sizeof(Item));
    }
    int k = i->heap_ptr = N++;
    array[k].i = i;
    array[k].key = key;
    while (k > 0) {
        int k2 = (k - 1) / 2;
        if (array[k2].key <= array[k].key) break;
        Swap(k, k2);
        k = k2;
    }
}
template<typename FlowType, typename CostType>
inline void MinCost<FlowType, CostType>::PriorityQueue::DecreaseKey(Node *i, CostType key) {
    int k = i->heap_ptr;
    array[k].key = key;
    while (k > 0) {
        int k2 = (k - 1) / 2;
        if (array[k2].key <= array[k].key) break;
        Swap(k, k2);
        k = k2;
    }
}
template<typename FlowType, typename CostType>
inline typename MinCost<FlowType, CostType>::Node *
MinCost<FlowType, CostType>::PriorityQueue::RemoveMin(CostType &key) {
    if (N == 0) return NULL;
    Swap(0, N - 1);
    N--;
    int k = 0;
    while (1) {
        int k1 = 2 * k + 1, k2 = k1 + 1;
        if (k1 >= N) break;
        int k_min = (k2 >= N || array[k1].key <= array[k2].key) ? k1 : k2;
        if (array[k].key <= array[k_min].key) break;
        Swap(k, k_min);
        k = k_min;
    }
    key = array[N].key;
    return array[N].i;
}
#endif
#ifndef HALSKDJDFHALSJASFDFASJGLA
#define HALSKDJDFHALSJASFDFASJGLA
#ifndef __BLOCK_H__
#define __BLOCK_H__
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
template<class Type>
class Block {
public:
    /* Constructor. Arguments are the block size and
       (optionally) the pointer to the function which
       will be called if allocation failed; the message
       passed to this function is "Not enough memory!" */
    Block(int size, void (*err_function)(const char *) = NULL) {
        first = last = NULL;
        block_size = size;
        error_function = err_function;
    }
    /* Destructor. Deallocates all items added so far */
    ~Block() {
        while (first) {
            block *next = first->next;
            delete[] ((char *) first);
            first = next;
        }
    }
    /* Allocates 'num' consecutive items; returns pointer
       to the first item. 'num' cannot be greater than the
       block size since items must fit in one block */
    Type *New(int num = 1) {
        Type *t;
        if (!last || last->current + num > last->last) {
            if (last && last->next) last = last->next;
            else {
                block *next = (block * )
                new char[sizeof(block) + (block_size - 1) * sizeof(Type)];
                if (!next) {
                    if (error_function) (*error_function)("Not enough memory!");
                    exit(1);
                }
                if (last) last->next = next;
                else first = next;
                last = next;
                last->current = &(last->data[0]);
                last->last = last->current + block_size;
                last->next = NULL;
            }
        }
        t = last->current;
        last->current += num;
        return t;
    }
    /* Returns the first item (or NULL, if no items were added) */
    Type *ScanFirst() {
        for (scan_current_block = first; scan_current_block; scan_current_block = scan_current_block->next) {
            scan_current_data = &(scan_current_block->data[0]);
            if (scan_current_data < scan_current_block->current) return scan_current_data++;
        }
        return NULL;
    }
    /* Returns the next item (or NULL, if all items have been read)
       Can be called only if previous ScanFirst() or ScanNext()
       call returned not NULL. */
    Type *ScanNext() {
        while (scan_current_data >= scan_current_block->current) {
            scan_current_block = scan_current_block->next;
            if (!scan_current_block) return NULL;
            scan_current_data = &(scan_current_block->data[0]);
        }
        return scan_current_data++;
    }
    /* Marks all elements as empty */
    void Reset() {
        block *b;
        if (!first) return;
        for (b = first;; b = b->next) {
            b->current = &(b->data[0]);
            if (b == last) break;
        }
        last = first;
    }
/***********************************************************************/
private:
    typedef struct block_st {
        Type *current, *last;
        struct block_st *next;
        Type data[1];
    } block;
    int block_size;
    block *first;
    block *last;
    block *scan_current_block;
    Type *scan_current_data;
    void (*error_function)(const char *);
};
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
template<class Type>
class DBlock {
public:
    /* Constructor. Arguments are the block size and
       (optionally) the pointer to the function which
       will be called if allocation failed; the message
       passed to this function is "Not enough memory!" */
    DBlock(int size, void (*err_function)(const char *) = NULL) {
        first = NULL;
        first_free = NULL;
        block_size = size;
        error_function = err_function;
    }
    /* Destructor. Deallocates all items added so far */
    ~DBlock() {
        while (first) {
            block *next = first->next;
            delete[] ((char *) first);
            first = next;
        }
    }
    /* Allocates one item */
    Type *New() {
        block_item *item;
        if (!first_free) {
            block *next = first;
            first = (block * )
            new char[sizeof(block) + (block_size - 1) * sizeof(block_item)];
            if (!first) {
                if (error_function) (*error_function)("Not enough memory!");
                exit(1);
            }
            first_free = &(first->data[0]);
            for (item = first_free; item < first_free + block_size - 1; item++)
                item->next_free = item + 1;
            item->next_free = NULL;
            first->next = next;
        }
        item = first_free;
        first_free = item->next_free;
        return (Type *) item;
    }
    /* Deletes an item allocated previously */
    void Delete(Type *t) {
        ((block_item *) t)->next_free = first_free;
        first_free = (block_item *) t;
    }
/***********************************************************************/
private:
    typedef union block_item_st {
        Type t;
        block_item_st *next_free;
    } block_item;
    typedef struct block_st {
        struct block_st *next;
        block_item data[1];
    } block;
    int block_size;
    block *first;
    block_item *first_free;
    void (*error_function)(const char *);
};
#endif
class PerfectMatching {
public:
#ifdef PERFECT_MATCHING_DOUBLE
                                                                                                                            typedef double REAL;
	#define PM_INFTY ((REAL)1e100)
#else
    typedef int REAL;
#define PM_INFTY (INT_MAX/2)
#endif
    typedef int NodeId;
    typedef int EdgeId;
    PerfectMatching(int nodeNum, int edgeNumMax);
    ~PerfectMatching();
    
    EdgeId AddEdge(NodeId i, NodeId j, REAL cost);
    
    
    
    void Solve(bool finish = true);
    
    
    int GetSolution(EdgeId e); 
    NodeId GetMatch(NodeId i); 
    
    
    
    
    void GetDualSolution(int *blossom_parents, REAL *twice_y);
    int GetBlossomNum();
    
    
    
    
    void StartUpdate();
    void FinishUpdate();
    
    REAL GetTwiceSum(
            NodeId i); 
    EdgeId AddNewEdge(NodeId i, NodeId j, REAL cost,
                      bool do_not_add_if_positive_slack = true); 
    void UpdateCost(EdgeId e, REAL delta_cost);
    
    
    
    
    struct Options {
        Options() : fractional_jumpstart(true),
                    dual_greedy_update_option(0),
                    dual_LP_threshold(0.00),
                    update_duals_before(false),
                    update_duals_after(false),
                    single_tree_threshold(1.00),
                    verbose(true) {}
        bool fractional_jumpstart; 
        int dual_greedy_update_option; 
        
        
        double dual_LP_threshold; 
        
        bool update_duals_before; 
        bool update_duals_after;  
        double single_tree_threshold; 
        bool verbose;
    } options;
    
    
    
    void Save(char *filename, int format = 0);
    
    
    
    
    
private:
    struct Node;
    struct Arc; 
    struct Edge; 
    struct Tree;
    struct TreeEdge;
    struct PQPointers;
    struct EdgeIterator;
    struct TreeEdgeIterator;
    struct LCATreeX;
    Node *nodes;
    Edge *edges;
    char *edges_orig;
    DBlock<Node> *blossoms;
    Tree *trees;
    DBlock<TreeEdge> *tree_edges;
    struct ExpandTmpItem {
        Node *i;
        Node *blossom_parent;
        Node *blossom_grandparent;
    };
    Block<ExpandTmpItem> *expand_tmp_list; 
    int node_num;
    int edge_num, edge_num_max;
    int tree_num, tree_num_max;
    Node *removed_first;
    int blossom_num;
    int removed_num;
    void *pq_buf;
    bool first_solve;
    
    struct Stat {
        int shrink_count;
        int expand_count;
        int grow_count;
        double shrink_time;
        double expand_time;
        double dual_time;
    } stat;
    
    void InitGreedy(bool allocate_trees = true);
    void InitGlobal(); 
    Node *FindBlossomRootInit(Edge *a0);
    void ShrinkInit(Edge *a0, Node *tree_root);
    void ExpandInit(Node *b);
    void AugmentBranchInit(Node *i0, Node *tree_root);
    void Finish(); 
    void ProcessNegativeEdge(Edge *a);
    void GetRealEndpoints(Edge *a, Node *&tail, Node *&head);
    Node *FindBlossomRoot(Edge *a0);
    void Shrink(Edge *a0);
    void Expand(Node *b);
    void Augment(Edge *a0);
    void AugmentBranch(Node *i0);
    void GrowNode(Node *i);
    void GrowTree(Node *root, bool new_subtree);
    bool ProcessEdge00(Edge *a, bool update_boundary_edge = true); 
    void ProcessSelfloop(Node *b, Edge *a);
    void AddTreeEdge(Tree *t0, Tree *t1);
    void ComputeEpsSingle(); 
    void ComputeEpsCC(); 
    void ComputeEpsSCC(); 
    void ComputeEpsGlobal(); 
    bool UpdateDuals();
    void FreeRemoved();
    void CommitEps();
    void ReallocateEdges();
    void PrintAll();
};
int CheckPerfectMatchingOptimality(int node_num, int edge_num, int *edges, int *weights, PerfectMatching *pm,
                                   PerfectMatching::REAL threshold = (PerfectMatching::REAL) (1e-10));
double ComputePerfectMatchingCost(int node_num, int edge_num, int *edges, int *weights, PerfectMatching *pm);
#endif
struct GPMKDTree;
class GeomPerfectMatching {
public:
    typedef int REAL; 
    typedef int PointId;
    
    GeomPerfectMatching(int pointNum, int DIM);
    ~GeomPerfectMatching();
    
    
    
    PointId AddPoint(REAL *coord);
    
    
    REAL SolveComplete();
    
    struct GPMOptions {
        GPMOptions() : init_Delaunay(true), init_KNN(0), init_greedy(true), iter_max(0) {}
        
        
        
        
        
        bool init_Delaunay;  
        int init_KNN;       
        bool init_greedy;    
        int iter_max;   
        
        
    };
    struct PerfectMatching::Options options;
    struct GPMOptions gpm_options;
    REAL Solve();
    
    void AddInitialEdge(PointId i, PointId j);
    
    
    
    
    PointId GetMatch(PointId p) { return matching[p]; }
    REAL Dist(REAL *coord1, REAL *coord2);
    REAL Dist(PointId p, PointId q);
    
    
    
    
    
private:
    friend struct GPMKDTree;
    struct Edge {
        PointId head[2];
        Edge *next[2];
    };
    struct Node {
        Edge *first[2];
        int is_marked;
    };
    Node *nodes;
    Block<Edge> *edges;
    REAL *coords; 
    REAL *sums; 
    PointId *matching; 
    int DIM;
    int node_num, node_num_max;
    int edge_num;
    double graph_update_time;
#define GPM_ROUND(X) (REAL)( ( ((REAL)1 / 2) == 0 ) ? ((X)+0.5) : (X) )
#define GPM_GET_NORM2(X, c) { int d; X = 0; for (d=0; d<DIM; d++) X += ((double)(c)[d])*(c)[d]; }
#define GPM_GET_NORM(x, c) { double X; GPM_GET_NORM2(X, c); X = sqrt(X); x = GPM_ROUND(X); }
#define GPM_GET_DIST2(X, c1, c2) { int d; X = 0; for (d=0; d<DIM; d++) X += ((double)((c1)[d]-(c2)[d]))*((c1)[d]-(c2)[d]); }
#define GPM_GET_DIST(x, c1, c2) { double X; GPM_GET_DIST2(X, c1, c2); X = sqrt(X); x = GPM_ROUND(X); }
#define GPM_CHECK_NORM(result, c, x){if (threshold <= 0) result = false; else {            double X;            GPM_GET_NORM2(X, c);            result = (4*X < (double)x*threshold);        }    }
#define GPM_CHECK_DIST(result, c1, c2, x)    {        if (threshold <= 0) result = false;        else        {            double X;            GPM_GET_DIST2(X, c1, c2);            result = (4*X < (double)x*threshold);        }    }
    double Norm2(REAL *coord) {
        double norm2;
        GPM_GET_NORM2(norm2, coord);
        return norm2;
    }
    REAL Norm(REAL *coord) {
        REAL norm;
        GPM_GET_NORM (norm, coord);
        return norm;
    }
    double Dist2(REAL *coord1, REAL *coord2) {
        double dist2;
        GPM_GET_DIST2(dist2, coord1, coord2);
        return dist2;
    }
    double Dist2(PointId p, PointId q) { return Dist2(coords + DIM * p, coords + DIM * q); }
    void CompleteInitialMatching(); 
    void InitKNN(int K);
    void InitDelaunay();
    REAL ComputeCost(PointId *matching);
};
#ifndef HFKSJHFKJHARBABDAKFAF
#define HFKSJHFKJHARBABDAKFAF
#define PQ_INTERLEAVED_MULTIPASS
template<typename REAL>
class PriorityQueue {
public:
    struct Item {
        REAL slack;
        Item *parentPQ;
        union {
            struct {
                Item *leftPQ;
                Item *rightPQ;
            };
            REAL y_saved; 
        };
    };
    static void *AllocateBuf();
    static void DeallocateBuf(void *buf);
    static void ResetItem(Item *i);
    static bool isReset(Item *i);
    
    void Reset();
    void Add(Item *i);
#define Remove(i, buf) _Remove(i)
    void _Remove(Item *i);
    void Decrease(Item *i_old, Item *i_new, void *buf);
    Item *GetMin();
    
    void Update(REAL delta);
    void Merge(PriorityQueue<REAL> &dest);
    
    
    Item *GetAndResetFirst();
    Item *GetAndResetNext();
    Item *GetFirst();
    Item *GetNext(Item *i);
    
private:
    struct Buf {
    };
    Item *rootPQ;
    void RemoveRoot();
};
template<typename REAL>
inline void *PriorityQueue<REAL>::AllocateBuf() {
    return NULL;
}
template<typename REAL>
inline void PriorityQueue<REAL>::DeallocateBuf(void *_buf) {
}
template<typename REAL>
inline void PriorityQueue<REAL>::ResetItem(Item *i) {
    i->parentPQ = NULL;
}
template<typename REAL>
inline bool PriorityQueue<REAL>::isReset(Item *i) {
    return (i->parentPQ == NULL);
}
template<typename REAL>
inline void PriorityQueue<REAL>::Reset() {
    rootPQ = NULL;
}
/*
template <typename REAL> inline void PriorityQueue<REAL>::RemoveRoot()
{
	Item* r = rootPQ;
	PriorityQueue<REAL> pq;
	pq.rootPQ = rootPQ;
	rootPQ = NULL;
	Item* i;
	for (i=pq.GetAndResetFirst(); i; i=pq.GetAndResetNext())
	{
		if (i != r) Add(i);
	}
	r->parentPQ = NULL;
}
*/
#define MERGE_PQ(i, j)    {        if (i->slack <= j->slack)        {            j->rightPQ = i->leftPQ;            if (j->rightPQ) j->rightPQ->parentPQ = j;            j->parentPQ = i;            i->leftPQ = j;        }        else        {            i->rightPQ = j->leftPQ;            if (i->rightPQ) i->rightPQ->parentPQ = i;            i->parentPQ = j;            j->leftPQ = i;            i = j;        }    }
template<typename REAL>
inline void PriorityQueue<REAL>::RemoveRoot() {
    Item *i = rootPQ->leftPQ;
    rootPQ->parentPQ = NULL;
    if (i) {
#ifdef PQ_MULTIPASS
                                                                                                                                while ( i->rightPQ )
		{
			Item** prev_ptr = &rootPQ;
			while ( 1 )
			{
				if (i->rightPQ)
				{
					Item* j = i->rightPQ;
					Item* next = j->rightPQ;
					MERGE_PQ(i, j);
					*prev_ptr = i;
					if (!next) { i->rightPQ = NULL; break; }
					prev_ptr = &i->rightPQ;
					i = next;
				}
				else
				{
					*prev_ptr = i;
					i->rightPQ = NULL;
					break;
				}
			}
			i = rootPQ;
		}
#endif
#ifdef PQ_INTERLEAVED_MULTIPASS
        while (i->rightPQ) {
            Item *prev = NULL;
            while (i) {
                Item *next;
                if (i->rightPQ) {
                    Item *j = i->rightPQ;
                    next = j->rightPQ;
                    MERGE_PQ(i, j);
                } else next = NULL;
                i->rightPQ = prev;
                prev = i;
                i = next;
            }
            i = prev;
        }
#endif
        i->parentPQ = i;
    }
    rootPQ = i;
}
template<typename REAL>
inline void PriorityQueue<REAL>::Add(Item *i) {
    if (!rootPQ) {
        rootPQ = i;
        i->parentPQ = i;
        i->leftPQ = i->rightPQ = NULL;
    } else if (i->slack <= rootPQ->slack) {
        rootPQ->parentPQ = i;
        i->leftPQ = rootPQ;
        i->rightPQ = NULL;
        rootPQ = i;
        i->parentPQ = i;
    } else {
        i->leftPQ = NULL;
        i->rightPQ = rootPQ->leftPQ;
        if (i->rightPQ) i->rightPQ->parentPQ = i;
        rootPQ->leftPQ = i;
        i->parentPQ = rootPQ;
    }
}
template<typename REAL>
inline void PriorityQueue<REAL>::_Remove(Item *i) {
    Item *p = i->parentPQ;
    if (p == i) RemoveRoot();
    else {
        if (i->rightPQ) i->rightPQ->parentPQ = p;
        if (p->leftPQ == i) p->leftPQ = i->rightPQ;
        else p->rightPQ = i->rightPQ;
        if (i->leftPQ) {
            i->parentPQ = i;
            i->rightPQ = NULL;
            PriorityQueue<REAL> pq;
            pq.rootPQ = i;
            pq.RemoveRoot();
            pq.Merge(*this);
        } else i->parentPQ = NULL;
    }
}
template<typename REAL>
inline void PriorityQueue<REAL>::Decrease(Item *i_old, Item *i_new, void *_buf) {
    if (i_old->parentPQ == i_old) {
        if (i_old != i_new) {
            rootPQ = i_new;
            i_new->parentPQ = i_new;
            i_new->leftPQ = i_old->leftPQ;
            i_new->rightPQ = NULL;
            if (i_new->leftPQ) i_new->leftPQ->parentPQ = i_new;
            i_old->parentPQ = NULL;
        }
    } else {
        Remove(i_old, _buf);
        Add(i_new);
    }
}
template<typename REAL>
inline typename PriorityQueue<REAL>::Item *PriorityQueue<REAL>::GetMin() {
    return rootPQ;
}
template<typename REAL>
inline void PriorityQueue<REAL>::Merge(PriorityQueue<REAL> &dest) {
    if (!rootPQ) return;
    if (!dest.rootPQ) dest.rootPQ = rootPQ;
    else {
        if (rootPQ->slack < dest.rootPQ->slack) {
            Item *j = rootPQ;
            rootPQ = dest.rootPQ;
            dest.rootPQ = j;
        }
        rootPQ->rightPQ = dest.rootPQ->leftPQ;
        if (rootPQ->rightPQ) rootPQ->rightPQ->parentPQ = rootPQ;
        rootPQ->parentPQ = dest.rootPQ;
        dest.rootPQ->leftPQ = rootPQ;
    }
    rootPQ = NULL;
}
template<typename REAL>
inline void PriorityQueue<REAL>::Update(REAL delta) {
    if (!rootPQ) return;
    Item *i = rootPQ;
    while (i->leftPQ) i = i->leftPQ;
    while (1) {
        i->slack += delta;
        if (i->rightPQ) {
            i = i->rightPQ;
            while (i->leftPQ) i = i->leftPQ;
        } else {
            while (1) {
                Item *j = i;
                i = i->parentPQ;
                if (i == j) return;
                if (i->leftPQ == j) break;
            }
        }
    }
}
struct PerfectMatching::Edge : PriorityQueue<REAL>::Item {
    Node *head[2];
    Node *head0[2];
    Edge *next[2];
    Edge *prev[2];
};
template<typename REAL>
inline typename PriorityQueue<REAL>::Item *PriorityQueue<REAL>::GetAndResetFirst() {
    if (!rootPQ) return NULL;
    return GetAndResetNext();
}
template<typename REAL>
inline typename PriorityQueue<REAL>::Item *PriorityQueue<REAL>::GetAndResetNext() {
    if (!rootPQ) return NULL;
    Item *result = rootPQ;
    result->parentPQ = NULL;
    Item *i = rootPQ->leftPQ;
    if (!i) rootPQ = result->rightPQ;
    else {
        rootPQ = i;
        while (i->rightPQ) i = i->rightPQ;
        i->rightPQ = result->rightPQ;
    }
    return result;
}
template<typename REAL>
inline typename PriorityQueue<REAL>::Item *PriorityQueue<REAL>::GetFirst() {
    if (!rootPQ) return NULL;
    Item *i = rootPQ;
    while (i->leftPQ) i = i->leftPQ;
    return i;
}
template<typename REAL>
inline typename PriorityQueue<REAL>::Item *PriorityQueue<REAL>::GetNext(Item *i) {
    if (i->rightPQ) {
        i = i->rightPQ;
        while (i->leftPQ) i = i->leftPQ;
        return i;
    }
    while (1) {
        Item *j = i;
        i = i->parentPQ;
        if (i == j) return NULL;
        if (i->leftPQ == j) return i;
    }
}
#endif
#ifndef NJAKSDTHASKJERAXJGFBZJDLAGZ
#define NJAKSDTHASKJERAXJGFBZJDLAGZ
#if defined (PM_TIMER_MSVC) || defined (PM_TIMER_CLOCK_GETTIME) || defined (PM_TIMER_GETRUSAGE) || defined (PM_TIMER_EXTERNAL) || defined (PM_TIMER_NONE)
#else
#ifdef _MSC_VER
#define PM_TIMER_MSVC
#elif defined(__APPLE_CC__)
#define PM_TIMER_GETRUSAGE
#else
#define PM_TIMER_CLOCK_GETTIME
#endif
#endif
#ifdef PM_TIMER_MSVC
                                                                                                                        #include <windows.h>
	inline double get_time()
	{
		LARGE_INTEGER t, frequency;
		QueryPerformanceCounter(&t);
		QueryPerformanceFrequency(&frequency);
		return (double)t.QuadPart/(double)frequency.QuadPart;
	}
#endif
#ifdef PM_TIMER_CLOCK_GETTIME
#include <time.h>
inline double get_time() {
    struct timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    return (double) t.tv_nsec * 1.00E-9 + (double) t.tv_sec;
}
#endif
#ifdef PM_TIMER_GETRUSAGE
                                                                                                                        #include <sys/resource.h>
	inline double get_time()
	{
		struct rusage t;
		getrusage (RUSAGE_SELF, &t);
		return (double)t.ru_utime.tv_usec*1.00E-6 + (double)t.ru_utime.tv_sec;
	}
#endif
#ifdef PM_TIMER_EXTERNAL
extern double get_time();
#endif
#ifdef PM_TIMER_NONE
inline double get_time() { return 0; }
#endif
#endif
#ifndef NJKASKJTASJNBAJSNRBAJS
#define NJKASKJTASJNBAJSNRBAJS
#ifndef GNAKDLATHJSTHAJSRNAKSJDA
#define GNAKDLATHJSTHAJSRNAKSJDA
class LCATree {
public:
    typedef void *NodeId; 
    typedef int PreorderId;
    LCATree(int node_num_max);
    ~LCATree();
    
    
    PreorderId Add(NodeId i, NodeId i_parent);
    PreorderId AddRoot(NodeId i); 
    PreorderId GetLCA(PreorderId i, PreorderId j);
    
    
    
    
    void GetPenultimateNodes(PreorderId &i, PreorderId &j);
private:
    int n, n_max, K, k_max;
    int **array;
    NodeId *buf0;
    int *buf1;
    NodeId *parent_current;
    int *child_current;
    int *parents;
    int _GetLCA(int i, int j); 
    int GetLCADirect(int i, int j);
};
inline LCATree::LCATree(int _n_max) : n_max(_n_max), array(NULL) {
#ifdef LCA_BLOCKS
                                                                                                                            K = -2;
	n = n_max;
	while (n > 0) { K ++; n /= 2; }
	if (K < 1) K = 1;
#else
    K = 1;
    n = 0;
#endif
    parents = new int[n_max];
    buf0 = new NodeId[n_max];
    buf1 = new int[n_max];
    parent_current = buf0;
    child_current = buf1;
}
inline LCATree::~LCATree() {
    int k;
    delete[] parents;
    if (buf0) delete[] buf0;
    if (buf1) delete[] buf1;
    if (array) {
        for (k = 1; k <= k_max; k++) delete[] array[k];
        delete[] array;
    }
}
inline LCATree::PreorderId LCATree::Add(NodeId i, NodeId i_parent) {
    assert(n < n_max);
    if (n == 0) {
        *parent_current = i;
        *(++parent_current) = i_parent;
        parents[0] = -1;
    } else {
        if (i == *parent_current) {
            int c = *child_current--;
            while (1) {
                int c_next = parents[c];
                parents[c] = n;
                if (c_next < 0) break;
                c = c_next;
            }
            parent_current--;
        }
        if (i_parent == *parent_current) parents[n] = *child_current;
        else {
            *(++parent_current) = i_parent;
            parents[n] = -1;
            child_current++;
        }
    }
    *child_current = n;
    return n++;
}
inline LCATree::PreorderId LCATree::AddRoot(NodeId i) {
    assert(n < n_max);
    if (n > 0) {
        if (i != *parent_current || parent_current != buf0 + 1) {
            printf("Error in LCATree construction: wrong sequence of calls!\n");
            exit(1);
        }
        int c = *child_current--;
        while (1) {
            int c_next = parents[c];
            parents[c] = n;
            if (c_next < 0) break;
            c = c_next;
        }
        child_current++;
    }
    parents[n++] = -1;
    delete[] buf0;
    buf0 = NULL;
    delete[] buf1;
    buf1 = NULL;
    
    int b, k = 1, block_num = (n - 1) / K + 1;
    if (block_num < 3) return n - 1;
    int d = (block_num - 1) / 4;
    while (d) {
        k++;
        d >>= 1;
    }
    k_max = k;
    array = new int *[k_max + 1];
    array[0] = parents;
    for (k = 1, d = 2; k <= k_max; k++, d *= 2) {
        array[k] = new int[block_num - d];
        if (k == 1) {
            for (b = 0; b < block_num - d; b++) array[1][b] = GetLCADirect((b + 1) * K - 1, (b + d) * K);
        } else {
            for (b = 0; b < block_num - d; b++) {
                int i = array[k - 1][b];
                int j = array[k - 1][b + d / 2];
                if (i < j) i = j;
                j = array[1][b + d / 2 - 1];
                array[k][b] = (i > j) ? i : j;
            }
        }
    }
    return n - 1;
}
inline int LCATree::GetLCADirect(int i, int j) {
    while (i < j) i = parents[i];
    return i;
}
inline int LCATree::_GetLCA(int i, int j) {
#ifdef LCA_BLOCKS
                                                                                                                            int bi = i/K, bj = j/K;
	if (bi == bj) return GetLCADirect(i, j);
	int i_last = (bi+1)*K-1, j_first = bj*K;
	i = GetLCADirect(i, i_last);
	j = GetLCADirect(j_first, j);
	if (i < j) i = j;
	
	if (j_first - i_last == 1) j = parents[i_last];
	else
	{
		int k = 1, d = (bj-bi)/4;
		while (d) { k ++; d >>= 1; }
		int diff = 1<<k;
		
		j = (array[k][bi] > array[k][bj-diff]) ? array[k][bi] : array[k][bj-diff];
	}
	return (i > j) ? i : j;
#else
    if (j == i) return i;
    int k = 0, d = (j - i) / 2;
    while (d) {
        k++;
        d >>= 1;
    }
    int diff = 1 << k;
    
    return (array[k][i] > array[k][j - diff]) ? array[k][i] : array[k][j - diff];
#endif
}
inline LCATree::PreorderId LCATree::GetLCA(PreorderId i, PreorderId j) {
    if (i > j) {
        PreorderId k = i;
        i = j;
        j = k;
    }
    return _GetLCA(i, j);
}
inline void LCATree::GetPenultimateNodes(PreorderId &_i, PreorderId &_j) {
    int i, j, d, swap;
    if (_i < _j) {
        i = _i;
        j = _j;
        swap = 0;
    } else {
        i = _j;
        j = _i;
        swap = 1;
    }
    int r = _GetLCA(i, j);
    assert(i != r && j != r);
    while (parents[i] != r) {
        int i0 = parents[i];
        d = (j - i0) / 2;
        while ((i = _GetLCA(i0, i0 + d)) == r) d /= 2;
    }
    while (parents[j] != r) {
        int j0 = parents[j];
        d = (r - j0) / 2;
        while ((j = _GetLCA(j0, j0 + d)) == r) d /= 2;
    }
    if (swap == 0) {
        _i = i;
        _j = j;
    } else {
        _j = i;
        _i = j;
    }
}
#endif
struct Neighbors {
    typedef GeomPerfectMatching::REAL REAL;
    typedef GeomPerfectMatching::PointId PointId;
    Neighbors();
    ~Neighbors();
    int GetNum() { return num; }
    void Init(int K, PointId *array);
    void Add(PointId p, double dist);
    double GetMax();
private:
    void Swap(int k1, int k2);
    PointId *array;
    double *dist_array;
    int num, K, K_max;
};
struct GPMKDTree {
    typedef GeomPerfectMatching::REAL REAL;
    typedef GeomPerfectMatching::PointId PointId;
    GPMKDTree(int D, int point_num, REAL *coords, GeomPerfectMatching *GPM);
    ~GPMKDTree();
    
    void AddPerfectMatching(PointId *rev_mapping);
    void ComputeKNN(PointId p, int K, PointId *neighbors);
    
    void AddNegativeEdges(PointId p, PerfectMatching *pm);
    
public:
#define CHILDREN_MAX 2
    struct Node {
        Node *parent;
        int d; 
#define IS_LEAF(i) ((i)->d < 0)
        union {
            struct 
            {
                REAL coord;
                Node *first_child; 
            };
            struct 
            {
                PointId points[CHILDREN_MAX]; 
            };
        };
        int order;
    } *nodes;
    
    Node **rev_mapping;
    int D, DIM, point_num, node_num;
    REAL sum_max;
    REAL *traversing_buf;
    GeomPerfectMatching *GPM;
    Neighbors neighbors;
};
#endif
inline GeomPerfectMatching::REAL GeomPerfectMatching::Dist(REAL *coord1, REAL *coord2) {
    REAL dist;
    GPM_GET_DIST (dist, coord1, coord2);
    return dist;
}
inline GeomPerfectMatching::REAL GeomPerfectMatching::Dist(PointId p, PointId q) {
    return Dist(coords + DIM * p, coords + DIM * q);
}
void GeomPerfectMatching::CompleteInitialMatching() {
    if (options.verbose) printf("adding edges to make sure that a perfect matching exists...");
    PointId p, q;
    Edge *e;
    double len, len_min;
    int unmatched_num = 0, edge_num0 = edge_num;
    
    for (p = 0; p < node_num; p++) {
        if (nodes[p].is_marked) continue;
        q = -1;
        for (e = nodes[p].first[0]; e; e = e->next[0]) {
            if (nodes[e->head[0]].is_marked) continue;
            len = Dist2(p, e->head[0]);
            if (q < 0 || len_min > len) {
                q = e->head[0];
                len_min = len;
            }
        }
        if (q >= 0) {
            nodes[p].is_marked = nodes[q].is_marked = 1;
        } else unmatched_num++;
    }
    if (unmatched_num == 0) {
        for (p = 0; p < node_num; p++) nodes[p].is_marked = 0;
        return;
    }
    
    REAL *unmatched_coords = new REAL[unmatched_num * DIM];
    int *rev_mapping = new int[unmatched_num];
    unmatched_num = 0;
    for (p = 0; p < node_num; p++) {
        if (nodes[p].is_marked) nodes[p].is_marked = 0;
        else {
            memcpy(unmatched_coords + unmatched_num * DIM, coords + p * DIM, DIM * sizeof(REAL));
            rev_mapping[unmatched_num++] = p;
        }
    }
    GPMKDTree *kd_tree = new GPMKDTree(DIM, unmatched_num, unmatched_coords, this);
    kd_tree->AddPerfectMatching(rev_mapping);
    delete kd_tree;
    delete[] unmatched_coords;
    delete[] rev_mapping;
    if (options.verbose) printf("done (%d edges)\n", edge_num - edge_num0);
}
void GeomPerfectMatching::InitKNN(int K) {
    if (node_num != node_num_max) {
        printf("InitKNN() cannot be called before all points have been added!\n");
        exit(1);
    }
    if (options.verbose) printf("adding K nearest neighbors (K=%d)\n", K);
    int dir, k;
    PointId p;
    Edge *e;
    if (K > node_num - 1) K = node_num - 1;
    GPMKDTree *kd_tree = new GPMKDTree(DIM, node_num, coords, this);
    PointId *neighbors = new PointId[K];
    for (p = 0; p < node_num; p++) {
        for (dir = 0; dir < 2; dir++)
            for (e = nodes[p].first[dir]; e; e = e->next[dir]) {
                nodes[e->head[dir]].is_marked = 1;
            }
        kd_tree->ComputeKNN(p, K, neighbors);
        for (k = 0; k < K; k++) {
            if (nodes[neighbors[k]].is_marked) continue;
            AddInitialEdge(p, neighbors[k]);
            nodes[neighbors[k]].is_marked = 1;
        }
        for (dir = 0; dir < 2; dir++)
            for (e = nodes[p].first[dir]; e; e = e->next[dir]) {
                nodes[e->head[dir]].is_marked = 0;
            }
    }
    delete kd_tree;
    delete[] neighbors;
}
#ifdef DELAUNAY_TRIANGLE
                                                                                                                        #ifdef _MSC_VER
#pragma warning(disable: 4311)
#pragma warning(disable: 4312)
#endif
extern "C" {
#define ANSI_DECLARATORS
#define TRILIBRARY
#define NO_TIMER
#define main NO_MAIN_FUNCTION
#include "../triangle/triangle.c"
}
void GeomPerfectMatching::InitDelaunay()
{
	if (node_num < 16) return;
	if (options.verbose) printf("adding edges in Delaunay triangulation\n");
	int k;
	struct triangulateio in, out, vorout;
	in.numberofpoints = node_num;
	in.numberofpointattributes = 0;
	in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
	for (k=0; k<2*node_num; k++) in.pointlist[k] = coords[k];
	in.pointattributelist = NULL;
	in.pointmarkerlist = NULL;
	in.numberofsegments = 0;
	in.numberofholes = 0;
	in.numberofregions = 0;
	in.regionlist = 0;
	out.pointlist = (REAL *) NULL;
	out.pointattributelist = (REAL *) NULL;
	out.pointmarkerlist = (int *) NULL;
	out.trianglelist = (int *) NULL;
	out.triangleattributelist = (REAL *) NULL;
	out.neighborlist = (int *) NULL;
	out.segmentlist = (int *) NULL;
	out.segmentmarkerlist = (int *) NULL;
	out.edgelist = (int *) NULL;
	out.edgemarkerlist = (int *) NULL;
	vorout.pointlist = (REAL *) NULL;
	vorout.pointattributelist = (REAL *) NULL;
	vorout.edgelist = (int *) NULL;
	vorout.normlist = (REAL *) NULL;
	triangulate("pczAevn", &in, &out, &vorout);
	free(in.pointlist);
	free(out.pointlist);
	free(out.pointmarkerlist);
	free(out.trianglelist);
	free(out.neighborlist);
	free(out.segmentlist);
	free(out.segmentmarkerlist);
	free(out.edgemarkerlist);
	free(vorout.pointlist);
	free(vorout.pointattributelist);
	free(vorout.edgelist);
	free(vorout.normlist);
	for (k=0; k<out.numberofedges; k++) AddInitialEdge(out.edgelist[2*k], out.edgelist[2*k+1]);
	free(out.edgelist);
}
#else
void GeomPerfectMatching::InitDelaunay() {
    printf("You need to download the 'Triangle' software from \n\twww.cs.cmu.edu/~quake/triangle.html ,\nextract it to the directory GeomPerfectMatching and define DELAUNARY_TRIANG in GeomPerfectMatching.h\n");
    exit(1);
}
#endif
struct Node {
    PerfectMatching::REAL sum; 
    Node *match;
    Node *parent;
    Node *child;
    Node *sibling;
    int lca_preorder;
};
int CheckPerfectMatchingOptimality(int node_num, int edge_num, int *edges, int *weights, PerfectMatching *pm,
                                   PerfectMatching::REAL threshold) {
    int _i, _j, _e;
    Node *i;
    Node *j;
    int blossom_num = pm->GetBlossomNum();
    int *blossom_parents = new int[node_num + blossom_num];
    PerfectMatching::REAL *twice_y = new PerfectMatching::REAL[node_num + blossom_num];
    PerfectMatching::REAL y_blossom_min = 0;
    PerfectMatching::REAL slack_min = 0;
    PerfectMatching::REAL active_slack_max = 0;
    
    pm->GetDualSolution(blossom_parents, twice_y);
    Node *nodes = new Node[node_num + blossom_num + 1];
    memset(nodes, 0, (node_num + blossom_num + 1) * sizeof(Node));
    Node *ROOT = nodes + node_num + blossom_num;
    for (_i = 0, i = nodes; _i < node_num + blossom_num; _i++, i++) {
        i->sum = twice_y[_i];
        if (_i >= node_num && y_blossom_min > i->sum) y_blossom_min = i->sum;
        if (blossom_parents[_i] >= 0) {
            if (blossom_parents[_i] < node_num || blossom_parents[_i] >= node_num + blossom_num) {
                delete[] nodes;
                delete[] blossom_parents;
                delete[] twice_y;
                return 2;
            }
            i->parent = nodes + blossom_parents[_i];
            i->sibling = i->parent->child;
            i->parent->child = i;
        }
    }
    delete[] blossom_parents;
    delete[] twice_y;
    for (i = nodes; i < nodes + node_num + blossom_num; i++) {
        if (!i->parent) {
            i->parent = ROOT;
            i->sibling = ROOT->child;
            ROOT->child = i;
        }
    }
    LCATree *lca_tree = new LCATree(node_num + blossom_num + 1);
    Node **rev_mapping = new Node *[node_num + blossom_num];
    i = ROOT;
    while (1) {
        if (i->child) {
            if (i < nodes + node_num) {
                delete[] nodes;
                delete lca_tree;
                delete[] rev_mapping;
                return 2;
            }
            i->child->sum += i->sum;
            i = i->child;
        } else {
            if (i >= nodes + node_num) {
                delete[] nodes;
                delete lca_tree;
                delete[] rev_mapping;
                return 2;
            }
            while (1) {
                i->lca_preorder = lca_tree->Add(i, i->parent);
                rev_mapping[i->lca_preorder] = i;
                if (i->sibling) break;
                i = i->parent;
                if (i == ROOT) {
                    i->lca_preorder = lca_tree->AddRoot(i);
                    break;
                }
            }
            if (i == ROOT) break;
            i = i->sibling;
            i->sum += i->parent->sum;
        }
    }
    int matched_num = 0;
    for (_e = 0; _e < edge_num; _e++) {
        _i = edges[2 * _e];
        _j = edges[2 * _e + 1];
        if (_i < 0 || _j < 0 || _i >= node_num || _j >= node_num || _i == _j) {
            delete[] nodes;
            delete lca_tree;
            delete[] rev_mapping;
            return 2;
        }
        int lca_i = nodes[_i].lca_preorder;
        int lca_j = nodes[_j].lca_preorder;
        lca_tree->GetPenultimateNodes(lca_i, lca_j);
        i = rev_mapping[lca_i];
        j = rev_mapping[lca_j];
        PerfectMatching::REAL twice_slack =
                2 * weights[_e] - (nodes[_i].sum - i->parent->sum) - (nodes[_j].sum - j->parent->sum);
        if (slack_min > twice_slack) slack_min = twice_slack;
        if (pm->GetSolution(_e)) {
            if (pm->GetMatch(_i) != _j || pm->GetMatch(_j) != _i || i->match || j->match) {
                delete[] nodes;
                delete lca_tree;
                delete[] rev_mapping;
                return 2;
            }
            i->match = j;
            j->match = i;
            if (active_slack_max < twice_slack) active_slack_max = twice_slack;
            matched_num += 2;
        }
    }
    delete[] nodes;
    delete lca_tree;
    delete[] rev_mapping;
    if (matched_num != node_num) return 2;
    if (y_blossom_min < -threshold || slack_min < -threshold || active_slack_max > threshold) {
        printf("ERROR in CheckPerfectMatchingOptimality():\n");
        if (((PerfectMatching::REAL) 1 / 2) == 0)
            printf("\ty_blossom_min=%d\n\tslack_min=%d\n\tactive_slack_max=%d\n", (int) y_blossom_min, (int) slack_min,
                   (int) active_slack_max);
        else
            printf("\ty_blossom_min=%.15f\n\tslack_min=%.15f\n\tactive_slack_max=%.15f\n", (double) y_blossom_min,
                   (double) slack_min, (double) active_slack_max);
        return 1;
    }
    return 0;
}
double ComputePerfectMatchingCost(int node_num, int edge_num, int *edges, int *weights, PerfectMatching *pm) {
    int i;
    int j;
    int e;
    double cost = 0;
    int *nodes = new int[node_num];
    memset(nodes, 0, node_num * sizeof(int));
    for (e = 0; e < edge_num; e++) {
        if (pm->GetSolution(e)) {
            i = edges[2 * e];
            j = edges[2 * e + 1];
            nodes[i]++;
            nodes[j]++;
            cost += weights[e];
        }
    }
    for (i = 0; i < node_num; i++) {
        if (nodes[i] != 1) {
            printf("ComputeCost(): degree = %d!\n", nodes[i]);
            exit(1);
        }
    }
    delete[] nodes;
    return cost;
}
#ifndef ASKHAKJSNTJAKSNBBAVASRA
#define ASKHAKJSNTJAKSNBBAVASRA
#ifdef _MSC_VER
                                                                                                                        #pragma warning(disable: 4311)
#pragma warning(disable: 4312)
#endif
#define LCA_REPAIRS
#define IS_INT ( ((REAL)1 / 2) == 0 )
#define COST_FACTOR 2
#define PM_THRESHOLD ((REAL)1e-12)
struct PerfectMatching::Node {
    unsigned int is_outer: 1; 
    unsigned int flag: 2; 
    unsigned int is_tree_root: 1;
    unsigned int is_processed: 1;
    unsigned int is_blossom: 1;
    unsigned int is_marked: 1;
    unsigned int is_removed: 1;
    Edge *first[2];
    union {
        Arc *match; 
        Node *blossom_grandparent;
    };
    REAL y;
    union {
        struct 
        {
            Arc *blossom_sibling;
            Node *blossom_parent;
            union {
                Edge *blossom_selfloops;
                Node *blossom_ptr; 
#ifdef LCA_REPAIRS
                int lca_preorder; 
#endif
            };
            REAL blossom_eps; 
            
        };
        struct 
        {
            union {
                struct 
                {
                    Node *first_tree_child;
                    Node *tree_sibling_prev; 
                    Node *tree_sibling_next;
                };
                Arc *tree_parent; 
            };
            union {
                Tree *tree;
                Edge *best_edge;  
#ifdef LCA_REPAIRS
                int lca_size; 
                LCATreeX *lca;  
#endif
            };
        };
    };
};
typedef unsigned long POINTER_TYPE;
extern char dummy_array[2 * (sizeof(void *) == sizeof(POINTER_TYPE)) - 1];
#define ARC_TO_EDGE_PTR(a)       ( (Edge*) ( ((POINTER_TYPE)(a)) & (~1)      ) )
#define ARC_TO_EDGE_DIR(a)       ( (int)   ( ((POINTER_TYPE)(a)) & 1         ) )
#define EDGE_DIR_TO_ARC(a, dir)  ( (Arc*)  ( (char*)(a) + (dir)) )
#define ARC_REV(a) ( (Arc*) ( ((POINTER_TYPE)(a)) ^ 1 ) )
#define ARC_TAIL(a)  (ARC_TO_EDGE_PTR(a)->head [1-ARC_TO_EDGE_DIR(a)])
#define ARC_TAIL0(a) (ARC_TO_EDGE_PTR(a)->head0[1-ARC_TO_EDGE_DIR(a)])
#define ARC_HEAD(a)  (ARC_TO_EDGE_PTR(a)->head [ARC_TO_EDGE_DIR(a)])
#define ARC_HEAD0(a) (ARC_TO_EDGE_PTR(a)->head0[ARC_TO_EDGE_DIR(a)])
struct PerfectMatching::PQPointers {
    PriorityQueue<REAL> pq00; 
    union {
        PriorityQueue<REAL> pq01[2]; 
        struct 
        {
            PriorityQueue<REAL> pq0; 
            PriorityQueue<REAL> pq_blossoms;
        };
    };
};
struct PerfectMatching::Tree : PQPointers {
    REAL eps;
    TreeEdge *first[2];
    Node *root;
    PQPointers *pq_current;
    int dir_current;
    
    
    REAL eps_delta;
    Tree *next;
    union {
        int id;
        TreeEdge *dfs_parent;
    };
};
struct PerfectMatching::TreeEdge : PQPointers {
    Tree *head[2];
    TreeEdge *next[2];
};
#define GET_PENULTIMATE_BLOSSOM(j)    {        Node* jtmp1 = j;        while ( 1 )        {            if (!j->blossom_grandparent->is_outer) j = j->blossom_grandparent;            else if (j->blossom_grandparent != j->blossom_parent) j->blossom_grandparent = j->blossom_parent;            else break;        }        Node* jtmp2;        for ( ; jtmp1!=j; jtmp1=jtmp2)        {            jtmp2 = jtmp1->blossom_grandparent;            jtmp1->blossom_grandparent = j;        }    }
#define GET_PENULTIMATE_BLOSSOM2(j)    {        Node* jtmp1 = j;        Node* jtmp_prev = NULL;        while ( 1 )        {            if (!j->blossom_grandparent->is_outer) { jtmp_prev = j; j = j->blossom_grandparent; }            else if (j->blossom_grandparent != j->blossom_parent) j->blossom_grandparent = j->blossom_parent;            else break;        }        if (jtmp_prev)        {            Node* jtmp2;            for ( ; jtmp1!=jtmp_prev; jtmp1=jtmp2)            {                jtmp2 = jtmp1->blossom_grandparent;                jtmp1->blossom_grandparent = jtmp_prev;            }        }    }
struct PerfectMatching::EdgeIterator {
    Edge *a_last;
    int start_flag;
};
#define FOR_ALL_EDGES(i, a, dir, I)    for ( dir = (i->first[0]) ? 0 : 1, I.a_last = a = i->first[dir], I.start_flag = (a) ? 0 : 1;        a != I.a_last || (I.start_flag ++ == 0) || (dir ++ == 0 && (I.a_last = a = i->first[1]));        a = a->next[dir] )
#define CONTINUE_FOR_ALL_EDGES(i, a, dir, I)    for ( a = a->next[dir];        a != I.a_last || (I.start_flag ++ == 0) || (dir ++ == 0 && (I.a_last = a = i->first[1]));        a = a->next[dir] )
#define REMOVE_EDGE(i, a, dir)    {        if ((a)->prev[dir]==(a)) (i)->first[dir] = NULL;        else        {            (a)->prev[dir]->next[dir] = (a)->next[dir];            (a)->next[dir]->prev[dir] = (a)->prev[dir];            (i)->first[dir] = (a)->next[dir];        }    }
#define ADD_EDGE(i, a, dir)    {        if ((i)->first[dir])        {            (a)->prev[dir] = (i)->first[dir]->prev[dir];            (a)->next[dir] = (i)->first[dir];            (i)->first[dir]->prev[dir]->next[dir] = (a);            (i)->first[dir]->prev[dir] = (a);        }        else (i)->first[dir] = (a)->prev[dir] = (a)->next[dir] = (a);        (a)->head[1-(dir)] = (i);    }
#define MOVE_EDGE(i_old, i_new, a, dir)    {        REMOVE_EDGE(i_old, a, dir);        ADD_EDGE(i_new, a, dir);    }
#define GET_OUTER_HEAD(a, dir, j)    {        j = (a)->head[dir];        if (!j->is_outer)        {            Node* j_orig = j;            GET_PENULTIMATE_BLOSSOM(j);            j = j->blossom_parent;            int dir_rev = 1 - (dir);            MOVE_EDGE(j_orig, j, a, dir_rev);        }    }
#define GET_TREE_PARENT(child, parent)    {        Arc* a = (child)->tree_parent;        Edge* e = ARC_TO_EDGE_PTR(a);        int dir = ARC_TO_EDGE_DIR(a);        GET_OUTER_HEAD(e, dir, parent);    }
struct PerfectMatching::TreeEdgeIterator {
    TreeEdge **e_ptr;
};
#define FOR_ALL_TREE_EDGES(t, e, dir)    for ( dir = (t->first[0]) ? 0 : 1, e = t->first[dir];          e || (dir ++ == 0 && (e = t->first[1]));          e = e->next[dir] )
#define FOR_ALL_TREE_EDGES_X(t, e, dir, T)    for ( dir = (t->first[0]) ? 0 : 1, T.e_ptr = &t->first[dir], e = *T.e_ptr;          e || (dir ++ == 0 && (e = *(T.e_ptr = &t->first[1])));          e = *T.e_ptr )    if (e->head[dir] == NULL) { *T.e_ptr = e->next[dir]; tree_edges->Delete(e); }    else if ((T.e_ptr = &e->next[dir]))
#define MOVE_NODE_IN_TREE(i)    {        if ((i)->first_tree_child) (i) = (i)->first_tree_child;        else        {            while (!(i)->is_tree_root && !(i)->tree_sibling_next) { (i) = ARC_HEAD((i)->match); GET_TREE_PARENT(i, i); }            if ((i)->is_tree_root) break;            (i) = (i)->tree_sibling_next;        }    }
#define ADD_TREE_CHILD(i, j)    {        (j)->flag = 0;        (j)->tree = (i)->tree;        (j)->first_tree_child = NULL;        (j)->tree_sibling_next = (i)->first_tree_child;        if ((i)->first_tree_child)        {            (j)->tree_sibling_prev = (i)->first_tree_child->tree_sibling_prev;            (i)->first_tree_child->tree_sibling_prev = j;        }        else        {            (j)->tree_sibling_prev = j;        }        (i)->first_tree_child = j;    }
#define REMOVE_FROM_TREE(i)    {        if ((i)->tree_sibling_next) (i)->tree_sibling_next->tree_sibling_prev = (i)->tree_sibling_prev;        else        {            Node* i_NEXT = ARC_HEAD((i)->match); i_NEXT = ARC_HEAD(i_NEXT->tree_parent); i_NEXT = i_NEXT->first_tree_child;            i_NEXT->tree_sibling_prev = (i)->tree_sibling_prev;        }        if ((i)->tree_sibling_prev->tree_sibling_next) (i)->tree_sibling_prev->tree_sibling_next = (i)->tree_sibling_next;        else        {            Node* i_PARENT = ARC_HEAD((i)->match); i_PARENT = ARC_HEAD(i_PARENT->tree_parent);            i_PARENT->first_tree_child = (i)->tree_sibling_next;        }    }
#endif
inline void PerfectMatching::ProcessSelfloop(Node *b, Edge *a) {
    int dir;
    Node *j;
    Node *prev[2];
    for (dir = 0; dir < 2; dir++) {
        j = a->head[dir];
        GET_PENULTIMATE_BLOSSOM(j);
        prev[dir] = j;
    }
    if (prev[0] != prev[1]) {
        ADD_EDGE(prev[0], a, 1);
        ADD_EDGE(prev[1], a, 0);
        a->slack -= 2 * prev[0]->blossom_eps;
    } else {
        a->next[0] = prev[0]->blossom_selfloops;
        prev[0]->blossom_selfloops = a;
    }
}
void PerfectMatching::Expand(Node *b) {
    assert(b->is_blossom);
    assert(b->is_outer);
    assert(b->flag == 1);
    double start_time = get_time();
    Node *i;
    Node *j;
    Node *k;
    Edge *a;
    EdgeIterator I;
    int dir;
    ExpandTmpItem *tmp_item;
    Tree *t = b->tree;
    REAL eps = t->eps;
    Edge *a_augment = NULL;
    GET_TREE_PARENT(b, i);
    a = ARC_TO_EDGE_PTR(b->tree_parent);
    dir = ARC_TO_EDGE_DIR(b->tree_parent);
    j = a->head0[1 - dir];
    GET_PENULTIMATE_BLOSSOM(j);
    MOVE_EDGE(b, j, a, dir);
    a = ARC_TO_EDGE_PTR(b->match);
    dir = ARC_TO_EDGE_DIR(b->match);
    k = a->head0[1 - dir];
    GET_PENULTIMATE_BLOSSOM(k);
    MOVE_EDGE(b, k, a, dir);
    i = ARC_HEAD(k->blossom_sibling);
    while (1) {
        tmp_item = expand_tmp_list->New();
        tmp_item->i = i;
        tmp_item->blossom_parent = i->blossom_parent;
        tmp_item->blossom_grandparent = i->blossom_grandparent;
        i->flag = 2;
        
        i->is_outer = 1;
        while ((a = i->blossom_selfloops)) {
            i->blossom_selfloops = a->next[0];
            ProcessSelfloop(i, a);
        }
        i->is_outer = 0;
        if (i == k) break;
        i->match = i->blossom_sibling;
        j = ARC_HEAD(i->match);
        tmp_item = expand_tmp_list->New();
        tmp_item->i = j;
        tmp_item->blossom_parent = j->blossom_parent;
        tmp_item->blossom_grandparent = j->blossom_grandparent;
        j->flag = 2;
        
        j->is_outer = 1;
        while ((a = j->blossom_selfloops)) {
            j->blossom_selfloops = a->next[0];
            ProcessSelfloop(j, a);
        }
        j->is_outer = 0;
        j->match = ARC_REV(i->match);
        i = ARC_HEAD(j->blossom_sibling);
    }
    k->match = b->match;
    i = ARC_TAIL(b->tree_parent);
    Arc *aa = i->blossom_sibling;
    i->flag = 1;
    i->tree = b->tree;
    i->y += b->tree->eps;
    i->tree_parent = b->tree_parent;
    if (i != k) {
        Node **i_ptr;
        if (i->match == aa) {
            i = ARC_HEAD(i->match);
            i_ptr = &j;
            while (1) {
                aa = i->blossom_sibling;
                i->flag = 0;
                i->tree = b->tree;
                i->y -= t->eps;
                *i_ptr = i;
                i_ptr = &i->first_tree_child;
                i->tree_sibling_prev = i;
                i->tree_sibling_next = NULL;
                i = ARC_HEAD(aa);
                i->flag = 1;
                i->tree = b->tree;
                i->y += t->eps;
                i->tree_parent = ARC_REV(aa);
                if (i == k) break;
                i = ARC_HEAD(i->match);
            }
            *i_ptr = ARC_HEAD(k->match);
        } else {
            i = k;
            j = ARC_HEAD(k->match);
            do {
                i->tree_parent = i->blossom_sibling;
                i->flag = 1;
                i->tree = b->tree;
                i->y += b->tree->eps;
                i = ARC_HEAD(i->tree_parent);
                i->flag = 0;
                i->tree = b->tree;
                i->y -= b->tree->eps;
                i->first_tree_child = j;
                j = i;
                i->tree_sibling_prev = i;
                i->tree_sibling_next = NULL;
                i = ARC_HEAD(i->match);
            } while (i->flag != 1);
        }
        i = ARC_HEAD(k->match);
        j->tree_sibling_prev = i->tree_sibling_prev;
        j->tree_sibling_next = i->tree_sibling_next;
        if (i->tree_sibling_prev->tree_sibling_next) i->tree_sibling_prev->tree_sibling_next = j;
        else
            ARC_HEAD(b->tree_parent)->first_tree_child = j;
        if (i->tree_sibling_next) i->tree_sibling_next->tree_sibling_prev = j;
        else
            ARC_HEAD(b->tree_parent)->first_tree_child->tree_sibling_prev = j;
        i->tree_sibling_prev = i;
        i->tree_sibling_next = NULL;
    }
    
    i = k;
    while (1) {
        
        if (i->is_blossom) {
            a = ARC_TO_EDGE_PTR(i->match);
            REAL tmp = a->slack;
            a->slack = i->y;
            i->y = tmp;
            t->pq_blossoms.Add(a);
        }
        FOR_ALL_EDGES(i, a, dir, I) {
            j = a->head[dir];
            if (j->flag != 0) a->slack -= eps;
        }
        i->is_processed = 1;
        if (i->tree_parent == b->tree_parent) break;
        i = ARC_HEAD(i->tree_parent);
        
        FOR_ALL_EDGES(i, a, dir, I) {
            j = a->head[dir];
            if (j->flag == 2) {
                a->slack += eps;
                t->pq0.Add(a);
            } else if (j->flag == 0 && i < j) {
                a->slack += 2 * eps;
                t->pq00.Add(a);
            }
        }
        i->is_processed = 1;
        i = ARC_HEAD(i->match);
    }
    
    for (tmp_item = expand_tmp_list->ScanFirst(); tmp_item; tmp_item = expand_tmp_list->ScanNext()) {
        i = tmp_item->i;
        j = tmp_item->blossom_parent;
        tmp_item->blossom_parent = i->blossom_parent;
        i->blossom_parent = j;
        j = tmp_item->blossom_grandparent;
        tmp_item->blossom_grandparent = i->blossom_grandparent;
        i->blossom_grandparent = j;
    }
    for (dir = 0; dir < 2; dir++) {
        if (!b->first[dir]) continue;
        b->first[dir]->prev[dir]->next[dir] = NULL;
        Edge *a_next;
        for (a = b->first[dir]; a; a = a_next) {
            a_next = a->next[dir];
            i = a->head0[1 - dir];
            GET_PENULTIMATE_BLOSSOM2(i);
            ADD_EDGE(i, a, dir);
            GET_OUTER_HEAD(a, dir, j);
            if (i->flag == 1) continue;
            if (j->flag == 0 && j->tree != t) j->tree->pq_current->pq01[1 - j->tree->dir_current].Remove(a, pq_buf);
            if (i->flag == 2) {
                a->slack += eps;
                if (j->flag == 0) j->tree->pq0.Add(a);
            } else {
                a->slack += 2 * eps;
                if (j->flag == 2) t->pq0.Add(a);
                else if (j->flag == 0) {
                    if (j->tree != t) {
                        if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                        if (a->slack <= j->tree->eps + eps) a_augment = a;
                    }
                    j->tree->pq_current->pq00.Add(a);
                } else if (j->tree != t) {
                    if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                    j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
                }
            }
        }
    }
    for (tmp_item = expand_tmp_list->ScanFirst(); tmp_item; tmp_item = expand_tmp_list->ScanNext()) {
        i = tmp_item->i;
        i->blossom_parent = tmp_item->blossom_parent;
        i->blossom_grandparent = tmp_item->blossom_grandparent;
        i->is_outer = 1;
    }
    expand_tmp_list->Reset();
    b->is_removed = 1;
    b->tree_sibling_next = removed_first;
    removed_first = b;
    removed_num++;
    if (4 * removed_num > node_num) FreeRemoved();
    blossom_num--;
    stat.expand_count++;
    stat.expand_time += get_time() - start_time;
    if (a_augment) Augment(a_augment);
}
void PerfectMatching::FreeRemoved() {
    Node *i0;
    Node *i;
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        for (i = i0; !i->is_outer && !i->is_marked; i = i->blossom_parent) {
            i->is_marked = 1;
            if (i->blossom_grandparent->is_removed) i->blossom_grandparent = i->blossom_parent;
        }
    }
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        for (i = i0; !i->is_outer && i->is_marked; i = i->blossom_parent) {
            i->is_marked = 0;
        }
    }
    while ((i = removed_first)) {
        removed_first = i->tree_sibling_next;
        blossoms->Delete(i);
        removed_num--;
    }
    assert(removed_num == 0);
}
struct PerfectMatching::LCATreeX : LCATree {
    LCATreeX(int size) : LCATree(size) { rev_mapping = new Node *[size]; }
    ~LCATreeX() { delete[] rev_mapping; }
    Node **rev_mapping;
};
void PerfectMatching::StartUpdate() {
    Node *i0;
    Node *i;
    Node *j;
    Node *b;
    while ((i = removed_first)) {
        removed_first = i->tree_sibling_next;
        blossoms->Delete(i);
        removed_num--;
    }
    Edge *a;
    Edge *selfloop_first = NULL;
    Edge *selfloop_last = NULL;
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        i0->is_processed = 0;
        if (i0->is_outer) continue;
        i0->is_tree_root = 0;
        i0->blossom_ptr = NULL;
        i = i0;
        while (1) {
            j = i->blossom_parent;
            j->is_processed = 0;
            if (j->is_outer) {
                j->first_tree_child = i;
                break;
            }
            if (j->is_marked) break;
            if ((a = j->blossom_selfloops)) {
                if (selfloop_last) selfloop_last->next[1] = a;
                else selfloop_first = a;
                selfloop_last = a;
                a->next[1] = NULL;
            }
            j->blossom_ptr = i;
            i = j;
        }
        b = (i->blossom_parent->is_outer) ? i->blossom_parent : i->blossom_parent->blossom_grandparent;
#ifdef LCA_REPAIRS
        if (!b->is_marked) {
            b->lca_size = 1;
            b->is_marked = 1;
        }
#endif
        while (1) {
#ifdef LCA_REPAIRS
            b->lca_size++;
#endif
            ARC_TO_EDGE_PTR(i->blossom_sibling)->y_saved = i->y;
            i->y += i->blossom_parent->y;
            if (!i->is_blossom) break;
            i->is_marked = 1;
            j = i;
            i = i->blossom_ptr;
            j->blossom_grandparent = b;
        }
        i->blossom_grandparent = b;
    }
#ifdef LCA_REPAIRS
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        if (i0->is_outer) continue;
        b = i0->blossom_grandparent;
        if (!b->is_marked) continue;
        b->is_marked = 0;
        LCATreeX *lca = new LCATreeX(b->lca_size);
        b->blossom_ptr = b->first_tree_child;
        i = b;
        while (1) {
            if (i->blossom_ptr) i = i->blossom_ptr;
            else {
                while (1) {
                    if (i->is_outer) break;
                    i->lca_preorder = lca->Add(i, i->blossom_parent);
                    lca->rev_mapping[i->lca_preorder] = i;
                    i = ARC_HEAD(i->blossom_sibling);
                    if (i != i->blossom_parent->blossom_ptr) break;
                    i = i->blossom_parent;
                }
                if (i->is_outer) {
                    lca->AddRoot(i);
                    break;
                }
            }
        }
        b->lca = lca;
    }
#endif
    while ((a = selfloop_first)) {
        selfloop_first = a->next[1];
        do {
            Edge *a_next = a->next[0];
#ifdef LCA_REPAIRS
            int _i = a->head0[1]->lca_preorder;
            int _j = a->head0[0]->lca_preorder;
            Node *b = a->head0[1]->blossom_grandparent;
            b->lca->GetPenultimateNodes(_i, _j);
            i = b->lca->rev_mapping[_i];
            j = b->lca->rev_mapping[_j];
#else
            GetRealEndpoints(a, i, j);
#endif
            ADD_EDGE(i, a, 0);
            ADD_EDGE(j, a, 1);
            a->slack -= 2 * i->blossom_eps;
            a = a_next;
        } while (a);
    }
    /*
	for (i0=nodes; i0<nodes+node_num; i0++)
	{
		if (i0->is_outer) continue;
		b = i0->blossom_grandparent;
		if (b->lca)
		{
			delete b->lca;
			b->lca = NULL;
		}
	}
	*/
    nodes[node_num].first_tree_child = NULL;
}
void PerfectMatching::FinishUpdate() {
    Node *i0;
    Node *i;
    Node *j;
    Edge *a;
    EdgeIterator I;
    int dir;
    Tree *t;
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        if (i0->is_outer) continue;
#ifdef LCA_REPAIRS
        if (i0->blossom_grandparent->lca) {
            delete i0->blossom_grandparent->lca;
            i0->blossom_grandparent->lca = NULL;
        }
#endif
        
        if (!i0->blossom_grandparent->is_removed) {
            i = i0;
            do {
                i->y = ARC_TO_EDGE_PTR(i->blossom_sibling)->y_saved;
                i->is_marked = 0;
                i->blossom_selfloops = NULL;
                i = i->blossom_parent;
            } while (i->is_marked);
            continue;
        }
        
        i = i0->blossom_parent;
        while (1) {
            if (i->is_removed && !i->is_outer) break;
            REAL y_parent = (i->is_outer) ? 0 : i->blossom_parent->y;
            for (dir = 0; dir < 2; dir++) {
                if (!i->first[dir]) continue;
                i->first[dir]->prev[dir]->next[dir] = NULL;
                Edge *a_next;
                for (a = i->first[dir]; a; a = a_next) {
                    a_next = a->next[dir];
                    j = a->head0[1 - dir];
                    ADD_EDGE(j, a, dir);
                    a->slack += j->blossom_parent->y - y_parent;
                }
                i->first[dir] = NULL;
            }
            if (i->is_removed) break;
            j = i->blossom_parent;
            i->is_removed = 1;
            i->tree_sibling_next = removed_first;
            removed_first = i;
            i = j;
        }
        i0->y = ARC_TO_EDGE_PTR(i0->blossom_sibling)->y_saved;
        i0->is_outer = 1;
        i0->flag = 2;
        i0->is_tree_root = 1;
    }
    Node *blossom_list = nodes[node_num].first_tree_child;
    for (i = nodes; i < nodes + node_num; i++) {
        if (!i->is_tree_root) continue;
        i->first_tree_child = nodes[node_num].first_tree_child;
        nodes[node_num].first_tree_child = i;
        REAL slack_min = PM_INFTY;
        FOR_ALL_EDGES(i, a, dir, I) {
            if (slack_min > a->slack) slack_min = a->slack;
        }
        i->y += slack_min;
        FOR_ALL_EDGES(i, a, dir, I) a->slack -= slack_min;
    }
    tree_num = 0;
    for (i = nodes[node_num].first_tree_child; i != blossom_list; i = i->first_tree_child) {
        tree_num++;
        if (!i->is_tree_root) continue;
        FOR_ALL_EDGES(i, a, dir, I) {
            j = a->head[dir];
            if (a->slack <= 0 && j->is_tree_root) {
                i->is_tree_root = j->is_tree_root = 0;
                i->match = EDGE_DIR_TO_ARC(a, dir);
                j->match = EDGE_DIR_TO_ARC(a, 1 - dir);
                tree_num -= 2;
                break;
            }
        }
    }
    for (; i; i = i->first_tree_child) {
        if (i->is_removed) {
            i->is_tree_root = 0;
            continue;
        }
        tree_num++;
    }
    if (tree_num > tree_num_max) {
        if (trees) free(trees);
        tree_num_max = tree_num;
        trees = (Tree *) malloc(tree_num_max * sizeof(Tree));
    }
    t = trees;
    Node *last_root = &nodes[node_num];
    Node *i_next;
    for (i = nodes; i; i = i_next) {
        if (!i->is_blossom) i_next = (i < nodes + node_num) ? (i + 1) : blossom_list;
        else i_next = i->first_tree_child;
        if (!i->is_tree_root) continue;
        i->flag = 0;
        i->first_tree_child = NULL;
        i->tree_sibling_prev = last_root;
        last_root->tree_sibling_next = i;
        last_root = i;
        i->tree = t;
        t->root = i;
        t->eps = 0;
        t->first[0] = t->first[1] = NULL;
        t->pq_current = NULL;
        t->pq00.Reset();
        t->pq0.Reset();
        t->pq_blossoms.Reset();
        t++;
    }
    assert(t == trees + tree_num);
    last_root->tree_sibling_next = NULL;
    while ((i = removed_first)) {
        removed_first = i->tree_sibling_next;
        blossoms->Delete(i);
        blossom_num--;
    }
}
PerfectMatching::REAL PerfectMatching::GetTwiceSum(NodeId i) {
    assert(i >= 0 && i < node_num);
    return nodes[i].y;
}
inline void PerfectMatching::ProcessNegativeEdge(Edge *a) {
    int dir;
    Node *i;
    for (dir = 0; dir < 2; dir++) {
        i = a->head0[dir];
        if (i->is_outer) {
            if (!i->is_tree_root) {
                i->is_tree_root = 1;
                i = ARC_HEAD(i->match);
                assert(!i->is_tree_root && i->is_outer);
                i->is_tree_root = 1;
                if (i->is_blossom) {
                    i->first_tree_child = nodes[node_num].first_tree_child;
                    nodes[node_num].first_tree_child = i;
                }
            }
            return;
        }
        if (i->blossom_grandparent->is_removed) return;
    }
    Node *b = i->blossom_grandparent;
    assert(b->is_outer);
    if (!b->is_tree_root) {
        b->is_tree_root = 1;
        i = ARC_HEAD(b->match);
        assert(!i->is_tree_root && i->is_outer);
        i->is_tree_root = 1;
        if (i->is_blossom) {
            i->first_tree_child = nodes[node_num].first_tree_child;
            nodes[node_num].first_tree_child = i;
        }
    }
    b->is_removed = 1;
    b->tree_sibling_next = removed_first;
    removed_first = b;
}
PerfectMatching::EdgeId PerfectMatching::AddNewEdge(NodeId _i, NodeId _j, REAL cost, bool do_not_add_if_positive_slack) {
    assert(_i >= 0 && _i < node_num && _j >= 0 && _j < node_num && _i != _j);
    if (edge_num >= edge_num_max) ReallocateEdges();
    Node *i = nodes + _i;
    Node *j = nodes + _j;
    Edge *a = edges + edge_num;
    a->slack = cost * COST_FACTOR;
    a->head0[0] = j;
    a->head0[1] = i;
    Node *bi = (i->is_outer) ? i : i->blossom_grandparent;
    Node *bj = (j->is_outer) ? j : j->blossom_grandparent;
    if (bi == bj) {
#ifdef LCA_REPAIRS
        int _i = i->lca_preorder;
        int _j = j->lca_preorder;
        bi->lca->GetPenultimateNodes(_i, _j);
        i = bi->lca->rev_mapping[_i];
        j = bi->lca->rev_mapping[_j];
#else
        GetRealEndpoints(a, i, j);
#endif
        a->slack += i->blossom_parent->y + j->blossom_parent->y;
    } else {
        i = bi;
        j = bj;
    }
    a->slack -= a->head0[0]->y + a->head0[1]->y;
    if (a->slack >= 0 && do_not_add_if_positive_slack) return -1;
    ADD_EDGE(i, a, 0);
    ADD_EDGE(j, a, 1);
    PriorityQueue<REAL>::ResetItem(a);
    if (a->slack < 0) {
        ProcessNegativeEdge(a);
    }
    return edge_num++;
}
void PerfectMatching::UpdateCost(EdgeId e, REAL delta_cost) {
    assert(e >= 0 && e < edge_num);
    Edge *a = edges + e;
    a->slack += delta_cost * COST_FACTOR;
    if (a->slack == 0) return;
    if (a->slack > 0) {
        Node *i = a->head[1];
        Node *j = a->head[0];
        if (i->is_outer) {
            if (ARC_TO_EDGE_PTR(i->match) != a && ARC_TO_EDGE_PTR(j->match) != a) return;
        } else {
            if (ARC_TO_EDGE_PTR(i->blossom_sibling) != a && ARC_TO_EDGE_PTR(j->blossom_sibling) != a) return;
        }
    }
    ProcessNegativeEdge(a);
}
PerfectMatching::PerfectMatching(int nodeNum, int edgeNumMax)
        : node_num(nodeNum),
          edge_num(0),
          edge_num_max(edgeNumMax),
          trees(NULL),
          tree_num_max(0),
          removed_first(NULL),
          blossom_num(0),
          removed_num(0),
          first_solve(true) {
    if (node_num & 1) {
        printf("# of nodes is odd: perfect matching cannot exist\n");
        exit(1);
    }
    nodes = (Node *) malloc((node_num + 1) * sizeof(Node));
    edges_orig = (char *) malloc(edge_num_max * sizeof(Edge) + 1);
    edges = (Edge *) ((((POINTER_TYPE) edges_orig) & 1) ? (edges_orig + 1) : edges_orig);
    memset(nodes, 0, (node_num + 1) * sizeof(Node));
    blossoms = new DBlock<Node>(256);
    tree_edges = new DBlock<TreeEdge>(256);
    expand_tmp_list = new Block<ExpandTmpItem>(256);
    pq_buf = PriorityQueue<REAL>::AllocateBuf();
}
void PerfectMatching::Save(char *filename, int format) {
    if (!first_solve) {
        printf("Save() cannot be called after Solve()!\n");
        exit(1);
    }
    int e;
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        printf("Can't open %s\n", filename);
        exit(1);
    }
    if (format == 0) {
        fprintf(fp, "p edge %d %d\n", node_num, edge_num);
        for (e = 0; e < edge_num; e++) {
            fprintf(fp, "e %d %d %d\n", 1 + (int) (edges[e].head0[1] - nodes), 1 + (int) (edges[e].head0[0] - nodes),
                    (int) edges[e].slack / COST_FACTOR);
        }
    } else {
        fprintf(fp, "%d %d\n", node_num, edge_num);
        for (e = 0; e < edge_num; e++) {
            fprintf(fp, "%d %d %d\n", (int) (edges[e].head0[1] - nodes), (int) (edges[e].head0[0] - nodes),
                    (int) edges[e].slack / COST_FACTOR);
        }
    }
    fclose(fp);
}
PerfectMatching::~PerfectMatching() {
    free(nodes);
    free(edges_orig);
    delete blossoms;
    delete tree_edges;
    delete expand_tmp_list;
    if (trees) free(trees);
    PriorityQueue<REAL>::DeallocateBuf(pq_buf);
}
PerfectMatching::EdgeId PerfectMatching::AddEdge(NodeId _i, NodeId _j, REAL cost) {
    if (_i < 0 || _i >= node_num || _j < 0 || _j > node_num || _i == _j) {
        printf("wrong node id's! (%d,%d)\n", _i, _j);
        exit(1);
    }
    if (edge_num >= edge_num_max) ReallocateEdges();
    Node *i = nodes + _i;
    Node *j = nodes + _j;
    Edge *a = edges + edge_num;
    ADD_EDGE(i, a, 0);
    ADD_EDGE(j, a, 1);
    a->head0[0] = j;
    a->head0[1] = i;
    a->slack = cost * COST_FACTOR;
    PriorityQueue<REAL>::ResetItem(a);
    return edge_num++;
}
int PerfectMatching::GetSolution(EdgeId e) {
    assert(e >= 0 && e < edge_num);
    Edge *a = edges + e;
    return (a->head0[1]->match == EDGE_DIR_TO_ARC(a, 0)) ? 1 : 0;
}
PerfectMatching::NodeId PerfectMatching::GetMatch(NodeId i) {
    assert(i >= 0 && i < node_num);
    return (int) (ARC_HEAD0(nodes[i].match) - nodes);
}
void PerfectMatching::GetRealEndpoints(Edge *a, Node *&tail, Node *&head) {
    Node *i;
    Node *j;
    int delta = 0;
    for (i = a->head0[1]; !i->is_outer; i = i->blossom_parent, delta--) {}
    for (j = a->head0[0]; !j->is_outer; j = j->blossom_parent, delta++) {}
    if (i == j) {
        i = a->head0[1];
        j = a->head0[0];
        while (delta < 0) {
            i = i->blossom_parent;
            delta++;
        }
        while (delta > 0) {
            j = j->blossom_parent;
            delta--;
        }
        while (i->blossom_parent != j->blossom_parent) {
            i = i->blossom_parent;
            j = j->blossom_parent;
        }
    }
    tail = i;
    head = j;
    assert((i->is_outer && j->is_outer) || (i->blossom_parent == j->blossom_parent && !i->is_outer && !j->is_outer));
}
void PerfectMatching::ReallocateEdges() {
    printf("Warning: reallocating edges. Increasing edge_num_max in the constructor may improve memory efficiency!\n");
    edge_num_max = edge_num_max * 3 / 2 + 16;
    char *edges_orig_old = edges_orig;
    Edge *edges_old = edges;
    edges_orig = (char *) realloc(edges_orig_old, edge_num_max * sizeof(Edge) + 1);
    edges = (Edge *) ((((POINTER_TYPE) edges_orig_old) & 1) ? (edges_orig + 1) : edges_orig);
    if (((POINTER_TYPE) edges) & 1) {
        char *edges_orig_old2 = edges_orig;
        Edge *edges_old2 = edges;
        edges_orig = (char *) malloc(edge_num_max * sizeof(Edge) + 1);
        edges = (Edge *) ((((POINTER_TYPE) edges_orig_old) & 1) ? (edges_orig + 1) : edges_orig);
        memcpy(edges, edges_old2, edge_num * sizeof(Edge));
        free(edges_orig_old2);
    }
#define UPDATE_EDGE_PTR(ptr) ptr = (Edge*)((char*)(ptr) + ((char*)edges - (char*)edges_old))
#define UPDATE_ARC_PTR(ptr) ptr = (Arc*)((char*)(ptr) + ((char*)edges - (char*)edges_old))
    Node *i;
    Edge *a;
    for (a = edges; a < edges + edge_num; a++) {
        if (a->next[0]) UPDATE_EDGE_PTR(a->next[0]);
        if (a->next[1]) UPDATE_EDGE_PTR(a->next[1]);
        if (a->prev[0]) UPDATE_EDGE_PTR(a->prev[0]);
        if (a->prev[1]) UPDATE_EDGE_PTR(a->prev[1]);
    }
    if (first_solve) {
        for (i = nodes; i < nodes + node_num; i++) {
            if (i->first[0]) UPDATE_EDGE_PTR(i->first[0]);
            if (i->first[1]) UPDATE_EDGE_PTR(i->first[1]);
        }
    } else {
        Node *i0;
        for (i0 = nodes; i0 < nodes + node_num; i0++) {
            i = i0;
            while (1) {
                if (i->is_outer) {
                    UPDATE_ARC_PTR(i->match);
                    if (i->first[0]) UPDATE_EDGE_PTR(i->first[0]);
                    if (i->first[1]) UPDATE_EDGE_PTR(i->first[1]);
                    break;
                }
                UPDATE_ARC_PTR(i->blossom_sibling);
                if (i->first[0]) UPDATE_EDGE_PTR(i->first[0]);
                if (i->first[1]) UPDATE_EDGE_PTR(i->first[1]);
                i = i->blossom_parent;
                if (i->is_outer) {
                    if (i->is_marked) break;
                    i->is_marked = 1;
                } else {
                    if (!i->is_marked) break;
                    i->is_marked = 0;
                }
            }
        }
        for (i0 = nodes; i0 < nodes + node_num; i0++) {
            i = i0;
            while (1) {
                if (i->is_outer) break;
                i = i->blossom_parent;
                if (i->is_outer) {
                    if (!i->is_marked) break;
                    i->is_marked = 0;
                } else {
                    if (i->is_marked) break;
                    i->is_marked = 1;
                }
            }
        }
    }
}
int PerfectMatching::GetBlossomNum() {
    return blossom_num;
}
void PerfectMatching::GetDualSolution(int *blossom_parents, REAL *twice_y) {
    int _i0, id = node_num;
    int *child_ptr;
    Node *i0;
    Node *i;
    int *tmp_array = new int[blossom_num];
    int *tmp_array_ptr = tmp_array;
    for (_i0 = 0, i0 = nodes; _i0 < node_num; _i0++, i0++) {
        twice_y[_i0] = i0->y;
        if (i0->is_outer) {
            blossom_parents[_i0] = -1;
            continue;
        }
        child_ptr = &blossom_parents[_i0];
        i = i0->blossom_parent;
        while (1) {
            if (i->is_marked) {
                *child_ptr = i->lca_preorder;
                break;
            }
            i->is_marked = 1;
            *tmp_array_ptr++ = i->lca_preorder;
            *child_ptr = i->lca_preorder = id++;
            child_ptr = &blossom_parents[i->lca_preorder];
            twice_y[i->lca_preorder] = i->y;
            if (i->is_outer) {
                *child_ptr = -1;
                break;
            }
            i = i->blossom_parent;
        }
    }
    assert(id == node_num + blossom_num && tmp_array_ptr == tmp_array + blossom_num);
    tmp_array_ptr = tmp_array;
    for (_i0 = 0, i0 = nodes; _i0 < node_num; _i0++, i0++) {
        if (i0->is_outer) continue;
        i = i0->blossom_parent;
        while (1) {
            if (!i->is_marked) break;
            i->is_marked = 0;
            i->lca_preorder = *tmp_array_ptr++;
            if (i->is_outer) break;
            i = i->blossom_parent;
        }
    }
    delete[] tmp_array;
}
void PerfectMatching::Finish() {
#define IS_VALID_MATCH(i) ((Edge*)(i->match) >= edges && (Edge*)(i->match) < edges + edge_num)
    Node *i0;
    Node *i;
    Node *j;
    Node *k;
    Node *b;
    Node *b_prev;
    Node *b_prev_prev;
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        if (IS_VALID_MATCH(i0)) continue;
        b_prev = NULL;
        b = i0;
        do {
            b->blossom_grandparent = b_prev;
            b_prev = b;
            b = b->blossom_parent;
        } while (!IS_VALID_MATCH(b));
        b_prev_prev = b_prev->blossom_grandparent;
        while (1) {
            for (k = ARC_TAIL0(b->match); k->blossom_parent != b; k = k->blossom_parent) {}
            k->match = b->match;
            i = ARC_HEAD(k->blossom_sibling);
            while (i != k) {
                i->match = i->blossom_sibling;
                j = ARC_HEAD(i->match);
                j->match = ARC_REV(i->match);
                i = ARC_HEAD(j->blossom_sibling);
            }
            b = b_prev;
            if (!b->is_blossom) break;
            b_prev = b_prev_prev;
            b_prev_prev = b_prev->blossom_grandparent;
        }
    }
}
void PerfectMatching::AddTreeEdge(Tree *t0, Tree *t1) {
    TreeEdge *e = tree_edges->New();
    e->head[0] = t1;
    e->head[1] = t0;
    e->next[0] = t0->first[0];
    t0->first[0] = e;
    e->next[1] = t1->first[1];
    t1->first[1] = e;
    e->pq00.Reset();
    e->pq01[0].Reset();
    e->pq01[1].Reset();
    t1->pq_current = e;
    t1->dir_current = 0;
}
bool PerfectMatching::ProcessEdge00(Edge *a, bool update_boundary_edge) {
    int dir;
    Node *j;
    Node *prev[2];
    Node *last[2];
    for (dir = 0; dir < 2; dir++) {
        if (a->head[dir]->is_outer) {
            prev[dir] = NULL;
            last[dir] = a->head[dir];
        } else {
            j = a->head[dir];
            GET_PENULTIMATE_BLOSSOM(j);
            prev[dir] = j;
            last[dir] = prev[dir]->blossom_parent;
            
        }
    }
    if (last[0] != last[1]) {
        for (dir = 0; dir < 2; dir++) {
            j = a->head[dir];
            if (j != last[dir]) {
                int dir_rev = 1 - dir;
                MOVE_EDGE(j, last[dir], a, dir_rev);
            }
        }
        if (update_boundary_edge) a->slack -= 2 * a->head[0]->tree->eps;
        return true;
    }
    if (prev[0] != prev[1]) {
        for (dir = 0; dir < 2; dir++) {
            j = a->head[dir];
            if (j != prev[dir]) {
                int dir_rev = 1 - dir;
                MOVE_EDGE(j, prev[dir], a, dir_rev);
            }
        }
        a->slack -= 2 * prev[0]->blossom_eps;
        return false;
    }
    for (dir = 0; dir < 2; dir++) {
        j = a->head[1 - dir];
        REMOVE_EDGE(j, a, dir);
    }
    a->next[0] = prev[0]->blossom_selfloops;
    prev[0]->blossom_selfloops = a;
    return false;
}
inline void PerfectMatching::AugmentBranch(Node *i0) {
    int dir;
    Tree *t = i0->tree;
    Node *r = t->root;
    Node *tree_root_prev = r->tree_sibling_prev;
    Node *i;
    Node *j;
    Edge *a;
    EdgeIterator I;
    Arc *aa;
    REAL eps = t->eps;
    PriorityQueue<REAL>::Item *q;
    TreeEdge *e;
    TreeEdgeIterator T;
    Tree *t2;
    t = r->tree;
    t->pq_current = t;
    FOR_ALL_TREE_EDGES_X(t, e, dir, T) {
            t2 = e->head[dir];
            e->head[1 - dir] = NULL; 
            t2->pq_current = e;
            t2->dir_current = dir;
        }
    i = r->first_tree_child;
    if (i)
        while (1) {
            Node *i0 = i;
            i = ARC_HEAD(i->match);
            if (i->is_processed) {
                if (i->is_blossom) {
                    a = ARC_TO_EDGE_PTR(i->match);
                    REAL tmp = a->slack;
                    a->slack = i->y;
                    i->y = tmp;
                    PriorityQueue<REAL>::ResetItem(a);
                }
                FOR_ALL_EDGES(i, a, dir, I) {
                    GET_OUTER_HEAD(a, dir, j);
                    if (j->flag == 0 && j->is_processed) {
                        if (j->tree != t) {
                            a->slack += eps;
                            if (PriorityQueue<REAL>::isReset(a)) j->tree->pq0.Add(a);
                        }
                    } else a->slack += eps;
                }
            }
            i = i0;
            MOVE_NODE_IN_TREE(i);
        }
    
    FOR_ALL_TREE_EDGES(t, e, dir) {
        t2 = e->head[dir];
        t2->pq_current = NULL;
        e->pq01[1 - dir].Merge(t2->pq0);
        for (q = e->pq00.GetFirst(); q; q = e->pq00.GetNext(q)) {
            q->slack -= eps;
            int dir2;
            for (dir2 = 0; dir2 < 2; dir2++) GET_OUTER_HEAD((Edge *) q, dir2, j);
        }
        e->pq00.Merge(t2->pq0);
        for (q = e->pq01[dir].GetAndResetFirst(); q; q = e->pq01[dir].GetAndResetNext()) {
            q->slack -= eps;
            int dir2;
            for (dir2 = 0; dir2 < 2; dir2++) GET_OUTER_HEAD((Edge *) q, dir2, j);
        }
    }
    for (q = t->pq0.GetAndResetFirst(); q; q = t->pq0.GetAndResetNext()) {
        q->slack -= eps;
        int dir2;
        for (dir2 = 0; dir2 < 2; dir2++) GET_OUTER_HEAD((Edge *) q, dir2, j);
    }
    for (q = t->pq00.GetAndResetFirst(); q; q = t->pq00.GetAndResetNext()) {
        ProcessEdge00((Edge *) q);
    }
    
    r->flag = 2;
    r->is_processed = 0;
    i = r->first_tree_child;
    r->y += eps;
    if (i)
        while (1) {
            j = ARC_HEAD(i->match);
            j->flag = 2;
            i->flag = 2;
            j->is_processed = 0;
            i->is_processed = 0;
            j->y -= eps;
            i->y += eps;
            MOVE_NODE_IN_TREE(i);
        }
    
    i = i0;
    if (!i0->is_tree_root) {
        j = ARC_HEAD(i0->match);
        GET_TREE_PARENT(j, i);
        j->match = aa = j->tree_parent;
        while (!i->is_tree_root) {
            j = ARC_HEAD(i->match);
            i->match = ARC_REV(aa);
            GET_TREE_PARENT(j, i);
            j->match = aa = j->tree_parent;
        }
        i->match = ARC_REV(aa);
    }
    r->is_tree_root = 0;
    tree_root_prev->tree_sibling_next = r->tree_sibling_next;
    if (r->tree_sibling_next) r->tree_sibling_next->tree_sibling_prev = tree_root_prev;
    tree_num--;
}
void PerfectMatching::Augment(Edge *a) {
    Node *j;
    int dir;
    for (dir = 0; dir < 2; dir++) {
        GET_OUTER_HEAD(a, dir, j);
        AugmentBranch(j);
        j->match = EDGE_DIR_TO_ARC(a, 1 - dir);
    }
    if (options.verbose) {
        int k = 1;
        while (k < tree_num) k *= 2;
        if (k == tree_num || tree_num <= 8 || (tree_num <= 64 && (tree_num % 8) == 0)) {
            printf("%d.", tree_num);
            fflush(stdout);
        }
    }
}
inline void PerfectMatching::GrowNode(Node *i) {
    
    
    Edge *a;
    EdgeIterator I;
    int dir;
    Node *j;
    Tree *t = i->tree;
    REAL eps = t->eps;
    Edge *a_augment = NULL;
    FOR_ALL_EDGES(i, a, dir, I) {
        GET_OUTER_HEAD(a, dir, j);
        if (j->flag == 2) {
            a->slack += eps;
            if (a->slack > 0) {
                t->pq0.Add(a);
            } else {
                j->flag = 1;
                j->tree = i->tree;
                j->tree_parent = EDGE_DIR_TO_ARC(a, 1 - dir);
                j->y += eps;
                j = ARC_HEAD(j->match);
                j->y -= eps;
                ADD_TREE_CHILD(i, j);
            }
        } else {
            if (j->flag == 0 && j->is_processed) {
                if (!PriorityQueue<REAL>::isReset(a)) j->tree->pq0.Remove(a, pq_buf);
                if (a->slack <= j->tree->eps && j->tree != t) a_augment = a;
                a->slack += eps;
                if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                j->tree->pq_current->pq00.Add(a);
            } else {
                a->slack += eps;
                if (j->flag == 1 && j->tree != t) {
                    if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                    j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
                }
            }
        }
    }
    
    i->is_processed = 1;
    if (!i->is_tree_root) {
        j = ARC_HEAD(i->match);
        
        j->is_processed = 1;
        if (j->is_blossom) {
            a = ARC_TO_EDGE_PTR(i->match);
            REAL tmp = a->slack;
            a->slack = j->y;
            j->y = tmp;
            t->pq_blossoms.Add(a);
        }
    }
    if (a_augment) Augment(a_augment);
    stat.grow_count++;
}
void PerfectMatching::GrowTree(Node *r, bool new_subtree) {
    
    Node *i = r;
    Node *j;
    Node *stop = r->tree_sibling_next;
    if (new_subtree && r->first_tree_child) stop = r->first_tree_child;
    Edge *a;
    EdgeIterator I;
    int dir;
    Tree *t = r->tree;
    REAL eps = t->eps;
    int tree_num0 = tree_num;
    while (1) {
        if (!i->is_tree_root) {
            
            i = ARC_HEAD(i->match);
            FOR_ALL_EDGES(i, a, dir, I) {
                GET_OUTER_HEAD(a, dir, j);
                if (j->flag == 2) a->slack -= eps;
                else {
                    if (j->flag == 0 && j->is_processed) {
                        if (!PriorityQueue<REAL>::isReset(a)) j->tree->pq0.Remove(a, pq_buf);
                        a->slack -= eps;
                        if (j->tree != t) {
                            if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                            j->tree->pq_current->pq01[1 - j->tree->dir_current].Add(a);
                        }
                    } else a->slack -= eps;
                }
            }
            i = ARC_HEAD(i->match);
        }
        
        GrowNode(i);
        if (tree_num != tree_num0) break;
        if (i->first_tree_child) i = i->first_tree_child;
        else {
            while (i != r && !i->tree_sibling_next) {
                i = ARC_HEAD(i->match);
                GET_TREE_PARENT(i, i);
            }
            i = i->tree_sibling_next;
        }
        if (i == stop) break;
    }
}
void PerfectMatching::Solve(bool finish) {
    Node *i;
    Node *j;
    Node *r;
    Node *r2;
    Node *r3 = NULL; 
    PriorityQueue<REAL>::Item *q;
    Edge *a;
    Tree *t;
    Tree *t2;
    TreeEdge *e;
    TreeEdgeIterator T;
    int dir;
    REAL eps;
    double start_time = get_time();
    if (IS_INT) {
        if (options.dual_greedy_update_option == 2) {
            printf("Fixed eps approach can only be used with floating point REAL!\n");
            printf("Change REAL to double in PerfectMatching.h and recompile\n");
            exit(1);
        }
        if (options.dual_LP_threshold > 0) {
            printf("LP approach can only be used with floating point REAL!\n");
            printf("Change REAL to double in PerfectMatching.h and recompile\n");
            exit(1);
        }
    }
    if (options.verbose) {
        printf("perfect matching with %d nodes and %d edges\n", node_num, edge_num);
        fflush(stdout);
    }
    if (first_solve) {
        if (options.verbose) {
            printf("    starting init...");
            fflush(stdout);
        }
        if (options.fractional_jumpstart) InitGlobal();
        else InitGreedy();
        if (options.verbose) printf("done [%.3f secs]. ", get_time() - start_time);
        first_solve = false;
    } else if (options.verbose) printf("    solving updated problem. ");
    if (options.verbose) {
        printf("%d trees\n    .", tree_num);
        fflush(stdout);
    }
    memset(&stat, 0, sizeof(Stat));
    
    
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        
        t = r->tree;
        
        EdgeIterator I;
        FOR_ALL_EDGES(r, a, dir, I) {
            j = a->head[dir];
            if (j->flag == 2) t->pq0.Add(a);
            else if (j->is_processed) {
                
                if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                j->tree->pq_current->pq00.Add(a);
            }
        }
        r->is_processed = 1;
        FOR_ALL_TREE_EDGES(t, e, dir) e->head[dir]->pq_current = NULL;
    }
    
    
    
    while (1) {
        int tree_num0 = tree_num;
        Stat stat0 = stat;
        REAL delta = 0;
        for (r = nodes[node_num].tree_sibling_next; r;) {
            r2 = r->tree_sibling_next;
            if (r2) r3 = r2->tree_sibling_next;
            t = r->tree;
            int tree_num1 = tree_num;
            
            
            
            t->pq_current = t;
            if (options.update_duals_before) {
                eps = PM_INFTY;
                Edge *a_augment = NULL;
                REAL eps_augment = PM_INFTY;
                if ((q = t->pq0.GetMin())) eps = q->slack;
                if ((q = t->pq_blossoms.GetMin()) && eps > q->slack) eps = q->slack;
                while ((q = t->pq00.GetMin())) {
                    if (ProcessEdge00((Edge *) q, false)) break;
                    t->pq00.Remove(q, pq_buf);
                }
                if (q && 2 * eps > q->slack) eps = q->slack / 2;
                FOR_ALL_TREE_EDGES_X(t, e, dir, T) {
                        t2 = e->head[dir];
                        t2->pq_current = e;
                        t2->dir_current = dir;
                        if ((q = e->pq00.GetMin()) && (!a_augment || eps_augment > q->slack - t2->eps)) {
                            a_augment = (Edge *) q;
                            eps_augment = q->slack - t2->eps;
                        }
                        if ((q = e->pq01[dir].GetMin()) && eps > q->slack + t2->eps) eps = q->slack + t2->eps;
                    }
                if (eps > eps_augment) eps = eps_augment;
                if (eps > t->eps) {
                    delta += eps - t->eps;
                    t->eps = eps;
                }
                if (a_augment && eps_augment <= t->eps) Augment(a_augment);
            } else {
                FOR_ALL_TREE_EDGES_X(t, e, dir, T) {
                        t2 = e->head[dir];
                        t2->pq_current = e;
                        t2->dir_current = dir;
                        if ((q = e->pq00.GetMin()) && (q->slack - t->eps <= t2->eps)) {
                            Augment((Edge *) q);
                            break;
                        }
                    }
            }
            
            
            
            eps = t->eps;
            REAL twice_eps = 2 * eps;
            while (tree_num1 == tree_num) {
                if ((q = t->pq0.GetMin()) && q->slack <= t->eps) {
                    a = (Edge *) q;
                    dir = (a->head[1]->flag == 2 && a->head[1]->is_outer) ? 1 : 0;
                    GET_OUTER_HEAD(a, 1 - dir, i);
                    j = a->head[dir];
                    
                    j->flag = 1;
                    j->tree = i->tree;
                    j->tree_parent = EDGE_DIR_TO_ARC(a, 1 - dir);
                    j->y += eps;
                    j = ARC_HEAD(j->match);
                    j->y -= eps;
                    ADD_TREE_CHILD(i, j);
                    GrowTree(j, true);
                } else if ((q = t->pq00.GetMin()) && q->slack <= twice_eps) {
                    t->pq00.Remove(q, pq_buf);
                    a = (Edge *) q;
                    if (ProcessEdge00(a)) Shrink(a);
                } else if ((q = t->pq_blossoms.GetMin()) && q->slack <= eps) {
                    t->pq_blossoms.Remove(q, pq_buf);
                    a = (Edge *) q;
                    j = (a->head[0]->flag == 1) ? a->head[0] : a->head[1];
                    REAL tmp = a->slack;
                    a->slack = j->y;
                    j->y = tmp;
                    Expand(j);
                } else break;
            }
            
            
            
            if (tree_num1 == tree_num) {
                t->pq_current = NULL;
                if (options.update_duals_after) {
                    eps = PM_INFTY;
                    Edge *a_augment = NULL;
                    REAL eps_augment = PM_INFTY;
                    if ((q = t->pq0.GetMin())) eps = q->slack;
                    if ((q = t->pq_blossoms.GetMin()) && eps > q->slack) eps = q->slack;
                    while ((q = t->pq00.GetMin())) {
                        if (ProcessEdge00((Edge *) q, false)) break;
                        t->pq00.Remove(q, pq_buf);
                    }
                    if (q && 2 * eps > q->slack) eps = q->slack / 2;
                    FOR_ALL_TREE_EDGES(t, e, dir) {
                        t2 = e->head[dir];
                        e->head[dir]->pq_current = NULL;
                        if ((q = e->pq00.GetMin()) && (!a_augment || eps_augment > q->slack - t2->eps)) {
                            a_augment = (Edge *) q;
                            eps_augment = q->slack - t2->eps;
                        }
                        if ((q = e->pq01[dir].GetMin()) && eps > q->slack + t2->eps) eps = q->slack + t2->eps;
                    }
                    if (eps > eps_augment) eps = eps_augment;
                    bool progress = false;
                    if (eps > t->eps) {
                        delta += eps - t->eps;
                        t->eps = eps;
                        progress = true;
                    }
                    if (a_augment && eps_augment <= t->eps) Augment(a_augment);
                    else if (progress && tree_num >= options.single_tree_threshold * node_num) {
                        
                        r = t->root;
                        continue;
                    }
                } else {
                    FOR_ALL_TREE_EDGES(t, e, dir) e->head[dir]->pq_current = NULL;
                }
            }
            
            
            
            r = r2;
            if (r && !r->is_tree_root) r = r3;
        }
        if (tree_num == 0) break;
        if (tree_num == tree_num0)
            
            
            
        {
            if (!UpdateDuals()) {
                if (!IS_INT && delta <= PM_THRESHOLD) 
                {
                    
                    int dual_greedy_update_option = options.dual_greedy_update_option;
                    options.dual_greedy_update_option = 2;
                    UpdateDuals();
                    options.dual_greedy_update_option = dual_greedy_update_option;
                }
            }
        }
    }
    if (finish) Finish();
    if (options.verbose) {
        printf("\ndone [%.3f secs]. %d grows, %d expands, %d shrinks\n", get_time() - start_time, stat.grow_count,
               stat.expand_count, stat.shrink_count);
        printf("    expands: [%.3f secs], shrinks: [%.3f secs], dual updates: [%.3f secs]\n", stat.expand_time,
               stat.shrink_time, stat.dual_time);
        fflush(stdout);
    }
}
PerfectMatching::Node *PerfectMatching::FindBlossomRoot(Edge *a0) {
    Node *i;
    Node *j;
    Node *_i[2];
    Node *r;
    int branch;
    _i[0] = ARC_HEAD(a0);
    _i[1] = ARC_TAIL(a0);
    branch = 0;
    while (1) {
        if (_i[branch]->is_marked) {
            r = _i[branch];
            j = _i[1 - branch];
            break;
        }
        _i[branch]->is_marked = 1;
        if (_i[branch]->is_tree_root) {
            j = _i[branch];
            i = _i[1 - branch];
            while (!i->is_marked) {
                i->is_marked = 1;
                i = ARC_HEAD(i->match);
                GET_TREE_PARENT(i, i);
            }
            r = i;
            break;
        }
        i = ARC_HEAD(_i[branch]->match);
        GET_TREE_PARENT(i, _i[branch]);
        branch = 1 - branch;
    }
    i = r;
    while (i != j) {
        i = ARC_HEAD(i->match);
        i = ARC_HEAD(i->tree_parent);
        i->is_marked = 0;
    }
    
    i = ARC_HEAD(a0);
    while (i != r) {
        i->is_marked = 0;
        i->is_outer = 0;
        i = ARC_HEAD(i->match);
        i->is_outer = 0;
        i = ARC_HEAD(i->tree_parent);
    }
    i = ARC_TAIL(a0);
    while (i != r) {
        i->is_marked = 0;
        i->is_outer = 0;
        i = ARC_HEAD(i->match);
        i->is_outer = 0;
        i = ARC_HEAD(i->tree_parent);
    }
    r->is_marked = 0;
    r->is_outer = 0;
    return r;
}
void PerfectMatching::Shrink(Edge *a0) {
    
    
    double start_time = get_time();
    int branch, dir;
    Node *r;
    Node *i;
    Node *j;
    Edge *a;
    Edge **a_inner_ptr;
    Arc *a_prev;
    Node *b = blossoms->New();
    Edge *a_augment = NULL;
    Edge *b_match;
    b->first[0] = b->first[1] = NULL;
    
    r = FindBlossomRoot(a0);
    Tree *t = r->tree;
    REAL eps = t->eps;
    b->first_tree_child = NULL;
    i = ARC_HEAD(a0);
    branch = 0;
    while (1) {
        if (i == r && branch) break;
        i->is_marked = 1;
        if (i == r) {
            branch = 1;
            i = ARC_TAIL(a0);
            continue;
        }
        
        REMOVE_FROM_TREE(i);
        
        if (i->first_tree_child) {
            j = i->first_tree_child;
            if (!b->first_tree_child) b->first_tree_child = j;
            else {
                Node *j_last = j->tree_sibling_prev;
                j->tree_sibling_prev = b->first_tree_child->tree_sibling_prev;
                b->first_tree_child->tree_sibling_prev->tree_sibling_next = j;
                b->first_tree_child->tree_sibling_prev = j_last;
            }
        }
        
        i = ARC_HEAD(i->match);
        i->is_marked = 1;
        if (i->is_blossom) {
            a = ARC_TO_EDGE_PTR(i->match);
            t->pq_blossoms.Remove(a, pq_buf);
            REAL tmp = a->slack;
            a->slack = i->y;
            i->y = tmp;
        }
        i = ARC_HEAD(i->tree_parent);
    }
    
    if (i->first_tree_child) {
        j = i->first_tree_child;
        if (!b->first_tree_child) b->first_tree_child = j;
        else {
            Node *j_last = j->tree_sibling_prev;
            j->tree_sibling_prev = b->first_tree_child->tree_sibling_prev;
            b->first_tree_child->tree_sibling_prev->tree_sibling_next = j;
            b->first_tree_child->tree_sibling_prev = j_last;
        }
    }
    
    b->is_removed = 0;
    b->is_outer = 1;
    b->flag = 0;
    b->is_blossom = 1;
    b->is_tree_root = r->is_tree_root;
    b->is_processed = 1;
    b->tree = t;
    b->y = -eps;
    b->is_marked = 0;
    
    b->tree_sibling_prev = r->tree_sibling_prev;
    b->tree_sibling_next = r->tree_sibling_next;
    Node *b_parent = NULL;
    if (!b->is_tree_root) {
        b_parent = ARC_HEAD(r->match);
        GET_TREE_PARENT(b_parent, b_parent);
    }
    if (b->tree_sibling_prev->tree_sibling_next) b->tree_sibling_prev->tree_sibling_next = b;
    else b_parent->first_tree_child = b;
    if (b->tree_sibling_next) b->tree_sibling_next->tree_sibling_prev = b;
    else if (b_parent) b_parent->first_tree_child->tree_sibling_prev = b;
    if (b->is_tree_root) {
        b->tree->root = b;
        b_match = NULL;
    } else {
        b->match = r->match;
        b_match = ARC_TO_EDGE_PTR(b->match);
    }
    REAL b_match_slack = 0; 
    if (b_match && ARC_HEAD(b->match)->is_blossom) {
        b_match_slack = b_match->slack;
        b_match->slack = ARC_HEAD(b->match)->y;
    }
    
    branch = 0;
    a_prev = EDGE_DIR_TO_ARC(a0, 0);
    i = ARC_HEAD(a_prev);
    while (1) {
        
        if (i->flag == 0) i->y += eps;
        else i->y -= eps;
        i->is_processed = 0;
        if (i->flag == 1) {
            Edge *a_prev;
            for (dir = 0; dir < 2; dir++)
                if (i->first[dir]) {
                    for (a_inner_ptr = &i->first[dir], a = *a_inner_ptr, a_prev = a->prev[dir], a_prev->next[dir] = NULL; a; a = *a_inner_ptr) {
                        Node *j0 = a->head[dir];
                        for (j = j0; !j->is_outer && !j->is_marked; j = j->blossom_parent) {}
                        if (j != j0) { /*assert(j->flag == 0);*/ int dir_rev = 1 - dir;
                            MOVE_EDGE(j0, j, a, dir_rev);
                        }
                        if (j->is_marked) 
                        {
                            a_inner_ptr = &a->next[dir];
                            a->prev[dir] = a_prev;
                            a_prev = a;
                            if (j->flag == 1) a->slack += eps;
                        } else 
                        {
                            *a_inner_ptr = a->next[dir];
                            ADD_EDGE(b, a, dir);
                            if (j->flag == 0 && j->tree != t) {
                                j->tree->pq_current->pq01[1 - j->tree->dir_current].Remove(a, pq_buf);
                                if (a->slack + eps <= j->tree->eps) a_augment = a;
                            }
                            a->slack += 2 * eps;
                            if (j->flag == 2) t->pq0.Add(a);
                            else if (j->flag == 0) {
                                if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                                j->tree->pq_current->pq00.Add(a);
                            } else if (j->tree != t) {
                                if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                                j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
                            }
                        }
                    }
                    if (i->first[dir]) {
                        a_prev->next[dir] = i->first[dir];
                        i->first[dir]->prev[dir] = a_prev;
                    }
                }
        }
        Arc *a_next = (i->flag == 0) ? i->match : i->tree_parent;
        i->blossom_parent = b;
        i->match = NULL;
        i->blossom_grandparent = b;
        i->blossom_selfloops = NULL;
        if (branch == 0) {
            i->blossom_sibling = a_next;
            if (i == r) {
                branch = 1;
                a_prev = ARC_REV(a0);
                i = ARC_HEAD(a_prev);
                if (i == r) break;
            } else {
                a_prev = i->blossom_sibling;
                i = ARC_HEAD(a_prev);
            }
        } else {
            i->blossom_sibling = ARC_REV(a_prev);
            a_prev = a_next;
            i = ARC_HEAD(a_prev);
            if (i == r) break;
        }
    }
    i->blossom_sibling = ARC_REV(a_prev);
    r->is_tree_root = 0;
    for (i = ARC_HEAD(r->blossom_sibling);; i = ARC_HEAD(i->blossom_sibling)) {
        i->is_marked = 0;
        i->blossom_eps = eps;
        if (i == r) break;
    }
    if (b_match) {
        if (ARC_HEAD(b->match)->is_blossom) {
            b_match->slack = b_match_slack;
        }
        dir = ARC_TO_EDGE_DIR(b->match);
        
        MOVE_EDGE(r, b, b_match, dir);
    }
    stat.shrink_count++;
    blossom_num++;
    stat.shrink_time += get_time() - start_time;
    if (a_augment) Augment(a_augment);
}
void PerfectMatching::InitGreedy(bool allocate_trees) {
    Node *i;
    int dir;
    Edge *a;
    EdgeIterator I;
    Tree *t = NULL;
    Node *last_root = &nodes[node_num];
    REAL slack_min;
    for (i = nodes; i < nodes + node_num; i++) i->y = PM_INFTY;
    for (a = edges; a < edges + edge_num; a++) {
        if (a->head[0]->y > a->slack) a->head[0]->y = a->slack;
        if (a->head[1]->y > a->slack) a->head[1]->y = a->slack;
    }
    for (a = edges; a < edges + edge_num; a++) {
        i = a->head[0];
        if (!i->is_outer) {
            i->is_outer = 1;
            i->y /= 2;
        }
        a->slack -= i->y;
        i = a->head[1];
        if (!i->is_outer) {
            i->is_outer = 1;
            i->y /= 2;
        }
        a->slack -= i->y;
    }
    tree_num = node_num;
    for (i = nodes; i < nodes + node_num; i++) {
        if (i->flag == 2) continue;
        slack_min = PM_INFTY;
        FOR_ALL_EDGES(i, a, dir, I) if (slack_min > a->slack) slack_min = a->slack;
        i->y += slack_min;
        FOR_ALL_EDGES(i, a, dir, I) {
            if (a->slack <= slack_min && i->flag == 0 && a->head[dir]->flag == 0) {
                i->flag = 2;
                a->head[dir]->flag = 2;
                i->match = EDGE_DIR_TO_ARC(a, dir);
                a->head[dir]->match = EDGE_DIR_TO_ARC(a, 1 - dir);
                tree_num -= 2;
            }
            a->slack -= slack_min;
        }
    }
    if (allocate_trees) {
        if (tree_num > tree_num_max) {
            if (trees) free(trees);
            tree_num_max = tree_num;
            trees = (Tree *) malloc(tree_num_max * sizeof(Tree));
        }
        t = trees;
    }
    for (i = nodes; i < nodes + node_num; i++) {
        if (i->flag != 0) continue;
        i->is_tree_root = 1;
        i->first_tree_child = NULL;
        i->tree_sibling_prev = last_root;
        last_root->tree_sibling_next = i;
        last_root = i;
        if (allocate_trees) {
            i->tree = t;
            t->root = i;
            t->eps = 0;
            t->first[0] = t->first[1] = NULL;
            t->pq_current = NULL;
            t->pq00.Reset();
            t->pq0.Reset();
            t->pq_blossoms.Reset();
            t++;
        }
    }
    last_root->tree_sibling_next = NULL;
}
PerfectMatching::Node *PerfectMatching::FindBlossomRootInit(Edge *a0) {
    Node *i;
    Node *j;
    Node *_i[2];
    Node *r;
    int branch;
    _i[0] = ARC_HEAD(a0);
    _i[1] = ARC_TAIL(a0);
    branch = 0;
    while (1) {
        if (!_i[branch]->is_outer) {
            r = _i[branch];
            j = _i[1 - branch];
            break;
        }
        _i[branch]->is_outer = 0;
        if (_i[branch]->is_tree_root) {
            j = _i[branch];
            i = _i[1 - branch];
            while (i->is_outer) {
                i->is_outer = 0;
                i = ARC_HEAD(i->match);
                i->is_outer = 0;
                i = ARC_HEAD(i->tree_parent);
            }
            r = i;
            break;
        }
        i = ARC_HEAD(_i[branch]->match);
        i->is_outer = 0;
        _i[branch] = ARC_HEAD(i->tree_parent);
        branch = 1 - branch;
    }
    i = r;
    while (i != j) {
        i = ARC_HEAD(i->match);
        i->is_outer = 1;
        i = ARC_HEAD(i->tree_parent);
        i->is_outer = 1;
    }
    return r;
}
void PerfectMatching::ShrinkInit(Edge *a0, Node *tree_root) {
    int branch, flag;
    Node *i;
    Node *j;
    Node *r;
    Arc *a_prev;
    Arc *aa;
    tree_root->flag = 2;
    i = tree_root->first_tree_child;
    if (i)
        while (1) {
            ARC_HEAD(i->match)->flag = 2;
            i->flag = 2;
            MOVE_NODE_IN_TREE(i);
        }
    r = FindBlossomRootInit(a0);
    if (!r->is_tree_root) {
        j = ARC_HEAD(r->match);
        j->match = aa = j->tree_parent;
        i = ARC_HEAD(aa);
        while (!i->is_tree_root) {
            j = ARC_HEAD(i->match);
            i->match = ARC_REV(aa);
            j->match = aa = j->tree_parent;
            i = ARC_HEAD(aa);
        }
        i->match = ARC_REV(aa);
    }
    tree_root->is_tree_root = 0;
    branch = 0;
    flag = 0;
    a_prev = EDGE_DIR_TO_ARC(a0, 0);
    i = ARC_HEAD(a_prev);
    while (1) {
        Arc *a_next = (flag == 0) ? i->match : i->tree_parent;
        flag = 1 - flag;
        i->flag = 0;
        i->match = NULL;
        if (branch == 0) {
            i->blossom_sibling = a_next;
            if (i == r) {
                branch = 1;
                flag = 0;
                a_prev = ARC_REV(a0);
                i = ARC_HEAD(a_prev);
                if (i == r) break;
            } else {
                a_prev = i->blossom_sibling;
                i = ARC_HEAD(a_prev);
            }
        } else {
            i->blossom_sibling = ARC_REV(a_prev);
            a_prev = a_next;
            i = ARC_HEAD(a_prev);
            if (i == r) break;
        }
    }
    i->blossom_sibling = ARC_REV(a_prev);
}
void PerfectMatching::ExpandInit(Node *k) {
    Node *i = ARC_HEAD(k->blossom_sibling);
    Node *j;
    while (1) {
        i->flag = 2;
        i->is_outer = 1;
        if (i == k) break;
        i->match = i->blossom_sibling;
        j = ARC_HEAD(i->match);
        j->flag = 2;
        j->is_outer = 1;
        j->match = ARC_REV(i->match);
        i = ARC_HEAD(j->blossom_sibling);
    }
}
void PerfectMatching::AugmentBranchInit(Node *i0, Node *r) {
    Node *tree_root_prev = r->tree_sibling_prev;
    Node *i;
    Node *j;
    Arc *aa;
    r->flag = 2;
    i = r->first_tree_child;
    if (i)
        while (1) {
            ARC_HEAD(i->match)->flag = 2;
            i->flag = 2;
            MOVE_NODE_IN_TREE(i);
        }
    i = i0;
    if (!i0->is_tree_root) {
        j = ARC_HEAD(i0->match);
        j->match = aa = j->tree_parent;
        i = ARC_HEAD(aa);
        while (!i->is_tree_root) {
            j = ARC_HEAD(i->match);
            i->match = ARC_REV(aa);
            j->match = aa = j->tree_parent;
            i = ARC_HEAD(aa);
        }
        i->match = ARC_REV(aa);
    }
    r->is_tree_root = 0;
    tree_root_prev->tree_sibling_next = r->tree_sibling_next;
    if (r->tree_sibling_next) r->tree_sibling_next->tree_sibling_prev = tree_root_prev;
    tree_num--;
}
void PerfectMatching::InitGlobal() {
    Node *i;
    Node *j;
    Node *r;
    Node *r2;
    Node *r3 = NULL; 
    Edge *a;
    EdgeIterator I;
    int dir;
    Tree TREE;
    enum {
        NONE, AUGMENT, SHRINK
    } flag;
    InitGreedy();
    for (i = nodes; i < nodes + node_num; i++) i->best_edge = NULL;
    PriorityQueue<REAL> pq;
    for (r = nodes[node_num].tree_sibling_next; r;) {
        r2 = r->tree_sibling_next;
        if (r2) r3 = r2->tree_sibling_next;
        i = r;
        pq.Reset();
        r->tree = &TREE;
        REAL eps = 0;
        Arc *critical_arc = NULL;
        REAL critical_eps = PM_INFTY;
        flag = NONE;
        Node *branch_root = i;
        while (1) {
            i->is_processed = 1;
            i->y -= eps;
            if (!i->is_tree_root) ARC_HEAD(i->match)->y += eps;
            FOR_ALL_EDGES(i, a, dir, I) {
                a->slack += eps;
                j = a->head[dir];
                if (j->tree == &TREE) {
                    
                    if (j->flag == 0) {
                        REAL slack = a->slack;
                        if (!j->is_processed) slack += eps;
                        if (2 * critical_eps > slack || critical_arc == NULL) {
                            flag = SHRINK;
                            critical_eps = slack / 2;
                            critical_arc = EDGE_DIR_TO_ARC(a, dir);
                            if (critical_eps <= eps) break;
                            
                        }
                    }
                } else if (j->flag == 0) {
                    
                    if (critical_eps >= a->slack || critical_arc == NULL) {
                        flag = AUGMENT;
                        critical_eps = a->slack;
                        critical_arc = EDGE_DIR_TO_ARC(a, dir);
                        if (critical_eps <= eps) break;
                        
                    }
                } else {
                    
                    if (a->slack > eps) {
                        if (a->slack < critical_eps) {
                            if (j->best_edge == NULL) {
                                j->best_edge = a;
                                pq.Add(a);
                            } else {
                                if (a->slack < j->best_edge->slack) {
                                    pq.Decrease(j->best_edge, a, pq_buf);
                                    j->best_edge = a;
                                }
                            }
                        }
                    } else {
                        assert(j->flag == 2 && !j->is_blossom && !ARC_HEAD(j->match)->is_blossom);
                        if (j->best_edge) pq.Remove(j->best_edge, pq_buf);
                        j->flag = 1;
                        j->tree = i->tree;
                        j->tree_parent = EDGE_DIR_TO_ARC(a, 1 - dir);
                        j = ARC_HEAD(j->match);
                        if (j->best_edge) pq.Remove(j->best_edge, pq_buf);
                        ADD_TREE_CHILD(i, j);
                    }
                }
            }
            if (dir < 2 && a) {
                Edge *atmp = a;
                int dirtmp = dir;
                CONTINUE_FOR_ALL_EDGES(i, atmp, dirtmp, I) atmp->slack += eps;
                break;
            }
            
            if (i->first_tree_child) i = i->first_tree_child;
            else {
                while (i != branch_root && !i->tree_sibling_next) {
                    i = ARC_HEAD(i->match);
                    i = ARC_HEAD(i->tree_parent);
                }
                if (i == branch_root) {
                    PriorityQueue<REAL>::Item *q = pq.GetMin();
                    if (q == NULL || q->slack >= critical_eps) {
                        eps = critical_eps;
                        break;
                    }
                    pq.Remove(q, pq_buf);
                    a = (Edge *) q;
                    dir = (a->head[0]->flag == 2) ? 0 : 1;
                    j = a->head[0];
                    Arc *aa = EDGE_DIR_TO_ARC(a, dir);
                    eps = a->slack;
                    assert(eps < critical_eps);
                    
                    i = ARC_TAIL(aa);
                    j = ARC_HEAD(aa);
                    assert(j->flag == 2 && !j->is_blossom && !ARC_HEAD(j->match)->is_blossom);
                    j->flag = 1;
                    j->tree = i->tree;
                    j->tree_parent = ARC_REV(aa);
                    j = ARC_HEAD(j->match);
                    if (j->best_edge) pq.Remove(j->best_edge, pq_buf);
                    ADD_TREE_CHILD(i, j);
                    i = branch_root = j;
                    continue;
                }
                i = i->tree_sibling_next;
            }
        }
        
        i = r;
        while (1) {
            if (i->is_processed) {
                i->y += eps;
                if (!i->is_tree_root) {
                    j = ARC_HEAD(i->match);
                    j->y -= eps;
                    REAL delta = eps - ARC_TO_EDGE_PTR(i->match)->slack;
                    FOR_ALL_EDGES(j, a, dir, I) a->slack += delta;
                    j->best_edge = NULL;
                }
                FOR_ALL_EDGES(i, a, dir, I) {
                    if (!PriorityQueue<REAL>::isReset(a)) {
                        assert(a->head[dir]->flag == 2 && a->head[dir]->best_edge == a);
                        a->head[dir]->best_edge = NULL;
                        PriorityQueue<REAL>::ResetItem(a);
                    }
                    a->slack -= eps;
                }
                i->is_processed = 0;
            } else {
                if (!i->is_tree_root) ARC_HEAD(i->match)->best_edge = NULL;
            }
            i->best_edge = NULL;
            MOVE_NODE_IN_TREE(i);
        }
        i = ARC_TAIL(critical_arc);
        j = ARC_HEAD(critical_arc);
        if (flag == SHRINK) {
            
            ShrinkInit(ARC_TO_EDGE_PTR(critical_arc), r);
        } else {
            
            AugmentBranchInit(i, r);
            if (j->is_outer) {
                AugmentBranchInit(j, j);
            } else {
                ExpandInit(j);
                tree_num--;
            }
            i->match = critical_arc;
            j->match = ARC_REV(critical_arc);
        }
        r = r2;
        if (r && !r->is_tree_root) r = r3;
    }
    if (tree_num > tree_num_max) {
        if (trees) free(trees);
        tree_num_max = tree_num;
        trees = (Tree *) malloc(tree_num_max * sizeof(Tree));
    }
    Tree *t = trees;
    for (r = nodes; r < nodes + node_num; r++) {
        if (!r->is_outer) {
            ExpandInit(r);
            r->is_tree_root = 1;
            r->flag = 0;
            r->first_tree_child = NULL;
            if (t == trees) {
                nodes[node_num].tree_sibling_next = r;
                r->tree_sibling_prev = &nodes[node_num];
            } else {
                (t - 1)->root->tree_sibling_next = r;
                r->tree_sibling_prev = (t - 1)->root;
            }
            r->tree = t;
            t->root = r;
            t->eps = 0;
            t->first[0] = t->first[1] = NULL;
            t->pq_current = NULL;
            t->pq00.Reset();
            t->pq0.Reset();
            t->pq_blossoms.Reset();
            t++;
        }
    }
    assert(t == trees + tree_num);
    if (t == trees) nodes[node_num].tree_sibling_next = NULL;
    else (t - 1)->root->tree_sibling_next = NULL;
}
void PerfectMatching::ComputeEpsGlobal() {
    Node *r;
    PriorityQueue<REAL>::Item *q;
    Tree *t;
    Tree *t2;
    TreeEdge *e;
    int i, j, k, N = 0, E = 0;
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        t->id = N;
        N += 2;
        for (k = 0; k < 2; k++)
            for (e = t->first[k]; e; e = e->next[k]) E += 6;
    }
    DualMinCost<REAL> *m = new DualMinCost<REAL>(N, E);
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        i = t->id;
        m->AddUnaryTerm(i, -1);
        m->SetLowerBound(i, 0);
        m->AddUnaryTerm(i + 1, 1);
        m->SetUpperBound(i + 1, 0);
        if (t->eps_delta < PM_INFTY) {
            m->SetUpperBound(i, t->eps_delta);
            m->SetLowerBound(i + 1, -t->eps_delta);
        }
        for (e = t->first[0]; e; e = e->next[0]) {
            t2 = e->head[0];
            if (t2 == NULL) continue;
            j = e->head[0]->id;
            if ((q = e->pq01[0].GetMin())) {
                m->AddConstraint(j, i, q->slack - t->eps + t2->eps);
                m->AddConstraint(i + 1, j + 1, q->slack - t->eps + t2->eps);
            }
            if ((q = e->pq01[1].GetMin())) {
                m->AddConstraint(i, j, q->slack - t2->eps + t->eps);
                m->AddConstraint(j + 1, i + 1, q->slack - t2->eps + t->eps);
            }
            if ((q = e->pq00.GetMin())) {
                m->AddConstraint(i + 1, j, q->slack - t->eps - t2->eps);
                m->AddConstraint(j + 1, i, q->slack - t->eps - t2->eps);
            }
        }
    }
    m->Solve();
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        i = t->id;
        t->eps_delta = (m->GetSolution(i) - m->GetSolution(i + 1)) / 2;
    }
    delete m;
}
void PerfectMatching::ComputeEpsSingle() {
    Node *r;
    PriorityQueue<REAL>::Item *q;
    Tree *t;
    Tree *t2;
    TreeEdge *e;
    REAL eps = PM_INFTY;
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        if (eps > t->eps_delta) eps = t->eps_delta;
        for (e = t->first[0]; e; e = e->next[0]) {
            t2 = e->head[0];
            if ((q = e->pq00.GetMin()) && 2 * eps > q->slack - t->eps - t2->eps) {
                eps = (q->slack - t->eps - t2->eps) / 2;
            }
        }
    }
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        r->tree->eps_delta = eps;
    }
}
void PerfectMatching::ComputeEpsCC() {
    Node *r;
    PriorityQueue<REAL>::Item *q;
    Tree *t0;
    Tree *t;
    Tree *t2;
    Tree *t_next;
    TreeEdge *e;
    REAL eps, eps2;
    Tree *queue_last;
    int dir;
    Tree *FIXED_TREE = trees - 1;
    int component_num = 0;
    TreeEdge **e_ptr;
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t0 = r->tree;
        t0->next = NULL;
    }
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t0 = r->tree;
        if (t0->next) continue;
        eps = t0->eps_delta;
        t0->next = queue_last = t = t0;
        while (1) {
            for (dir = 0; dir < 2; dir++)
                for (e_ptr = &t->first[dir], e = *e_ptr; e; e = *e_ptr) {
                    t2 = e->head[dir];
                    if (t2 == NULL) {
                        *e_ptr = e->next[dir];
                        tree_edges->Delete(e);
                        continue;
                    }
                    e_ptr = &e->next[dir];
                    REAL eps00 = ((q = e->pq00.GetMin())) ? (q->slack - t->eps - t2->eps) : PM_INFTY;
                    if (t2->next && t2->next != FIXED_TREE) {
                        if (2 * eps > eps00) eps = eps00 / 2;
                        continue;
                    }
                    REAL eps01[2];
                    eps01[dir] = ((q = e->pq01[dir].GetMin())) ? (q->slack - t->eps + t2->eps) : PM_INFTY;
                    eps01[1 - dir] = ((q = e->pq01[1 - dir].GetMin())) ? (q->slack - t2->eps + t->eps) : PM_INFTY;
                    if (t2->next == FIXED_TREE) eps2 = t2->eps_delta;
                    else if (eps01[0] > 0 && eps01[1] > 0) eps2 = 0;
                    else {
                        queue_last->next = t2;
                        queue_last = t2;
                        t2->next = t2;
                        if (eps > eps00) eps = eps00;
                        if (eps > t2->eps_delta) eps = t2->eps_delta;
                        continue;
                    }
                    if (eps > eps00 - eps2) eps = eps00 - eps2;
                    if (eps > eps2 + eps01[dir]) eps = eps2 + eps01[dir];
                }
            if (t->next == t) break;
            t = t->next;
        }
        for (t = t0;; t = t_next) {
            t->eps_delta = eps;
            t_next = t->next;
            t->next = FIXED_TREE;
            if (t_next == t) break;
        }
        component_num++;
    }
    
}
void PerfectMatching::ComputeEpsSCC() {
    PriorityQueue<REAL>::Item *q;
    Node *r;
    Tree *t0;
    Tree *t;
    Tree *t2;
    TreeEdge *e;
    TreeEdge **e_ptr;
    REAL eps;
    int c, dir;
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t0 = r->tree;
        t0->dfs_parent = NULL;
        for (dir = 0; dir < 2; dir++)
            for (e_ptr = &t0->first[dir], e = *e_ptr; e; e = *e_ptr) {
                t2 = e->head[dir];
                if (t2 == NULL) {
                    *e_ptr = e->next[dir];
                    tree_edges->Delete(e);
                    continue;
                }
                e_ptr = &e->next[dir];
            }
    }
    Tree *stack = NULL;
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t0 = r->tree;
        if (t0->dfs_parent) continue;
        t = t0;
        e = (t->first[0]) ? t->first[0] : t->first[1];
        t->dfs_parent = (TreeEdge *) trees;
        while (1) {
            if (e == NULL) {
                t->next = stack;
                stack = t;
                if (t == t0) break;
                e = t->dfs_parent;
                if (t == e->head[0]) {
                    t = e->head[1];
                    e = (e->next[0]) ? e->next[0] : t->first[1];
                } else {
                    t = e->head[0];
                    e = e->next[1];
                }
                continue;
            }
            if (e->head[1] == t) {
                if (e->head[0]->dfs_parent || !(q = e->pq01[0].GetMin()) || q->slack - t->eps + e->head[0]->eps > 0) {
                    e = (e->next[0]) ? e->next[0] : t->first[1];
                    continue;
                }
                t = e->head[0];
            } else {
                if (e->head[1]->dfs_parent || !(q = e->pq01[1].GetMin()) || q->slack - t->eps + e->head[1]->eps > 0) {
                    e = e->next[1];
                    continue;
                }
                t = e->head[1];
            }
            t->dfs_parent = e;
            e = (t->first[0]) ? t->first[0] : t->first[1];
        }
    }
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) r->tree->dfs_parent = NULL;
    int component_num = 0;
    while (stack) {
        t0 = stack;
        stack = t0->next;
        if (t0->dfs_parent) continue;
        t = t0;
        e = (t->first[0]) ? t->first[0] : t->first[1];
        t->dfs_parent = (TreeEdge *) trees;
        while (1) {
            if (e == NULL) {
                e = t->dfs_parent;
                t->dfs_parent = (TreeEdge *) ((char *) trees + component_num);
                if (t == t0) break;
                if (t == e->head[0]) {
                    t = e->head[1];
                    e = (e->next[0]) ? e->next[0] : t->first[1];
                } else {
                    t = e->head[0];
                    e = e->next[1];
                }
                continue;
            }
            if (e->head[1] == t) {
                if (e->head[0]->dfs_parent || !(q = e->pq01[1].GetMin()) || q->slack - e->head[0]->eps + t->eps > 0) {
                    e = (e->next[0]) ? e->next[0] : t->first[1];
                    continue;
                }
                t = e->head[0];
            } else {
                if (e->head[1]->dfs_parent || !(q = e->pq01[0].GetMin()) || q->slack - e->head[1]->eps + t->eps > 0) {
                    e = e->next[1];
                    continue;
                }
                t = e->head[1];
            }
            t->dfs_parent = e;
            e = (t->first[0]) ? t->first[0] : t->first[1];
        }
        component_num++;
    }
    Tree **array = new Tree *[component_num];
    memset(array, 0, component_num * sizeof(Tree *));
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        t->id = (int) ((char *) t->dfs_parent - (char *) trees);
        t->next = array[t->id];
        array[t->id] = t;
    }
    for (c = component_num - 1; c >= 0; c--) {
        eps = PM_INFTY;
        for (t = array[c]; t; t = t->next) {
            if (eps > t->eps_delta) eps = t->eps_delta;
            FOR_ALL_TREE_EDGES(t, e, dir) {
                t2 = e->head[dir];
                REAL eps00 = (q = e->pq00.GetMin()) ? (q->slack - t->eps - t2->eps) : PM_INFTY;
                REAL eps01[2];
                eps01[dir] = ((q = e->pq01[dir].GetMin())) ? (q->slack - t->eps + t2->eps) : PM_INFTY;
                eps01[1 - dir] = ((q = e->pq01[1 - dir].GetMin())) ? (q->slack - t2->eps + t->eps) : PM_INFTY;
                if (t2->id < c) {
                    if (eps > eps01[dir]) eps = eps01[dir];
                    if (eps > eps00) eps = eps00;
                } else if (t2->id == c) {
                    if (2 * eps > eps00) eps = eps00 / 2;
                } else {
                    if (eps > eps01[dir] + t2->eps_delta) eps = eps01[dir] + t2->eps_delta;
                    if (eps > eps00 - t2->eps_delta) eps = eps00 - t2->eps_delta;
                }
            }
        }
        for (t = array[c]; t; t = t->next) t->eps_delta = eps;
    }
    delete[] array;
    
}
void PerfectMatching::CommitEps() {
    printf("CommitEps()\n");
    Node *i;
    Node *j;
    Node *r;
    int dir;
    Edge *a;
    EdgeIterator I;
    Tree *t;
    TreeEdge *e;
    TreeEdge **e_ptr;
    REAL eps, eps2;
    PriorityQueue<REAL>::Item *q;
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        eps = t->eps;
        i = r;
        while (1) {
            i->y += eps;
            if (!i->is_tree_root) {
                Node *i0 = i;
                i = ARC_HEAD(i0->match);
                if (i->is_blossom) ARC_TO_EDGE_PTR(i0->match)->slack -= eps;
                else i->y -= eps;
                FOR_ALL_EDGES(i, a, dir, I) {
                    GET_OUTER_HEAD(a, dir, j);
                    a->slack += eps;
                    if (j->flag == 0) a->slack -= j->tree->eps;
                }
                i = i0;
            }
            MOVE_NODE_IN_TREE(i);
        }
        t->pq0.Update(-eps);
        PriorityQueue<REAL> pq00 = t->pq00;
        t->pq00.Reset();
        for (q = pq00.GetAndResetFirst(); q; q = pq00.GetAndResetNext()) {
            a = (Edge *) q;
            if (ProcessEdge00(a)) t->pq00.Add(a);
        }
        for (e_ptr = &t->first[0], e = *e_ptr; e; e = *e_ptr) {
            if (e->head[0] == NULL) {
                *e_ptr = e->next[0];
                tree_edges->Delete(e);
                continue;
            }
            e_ptr = &e->next[0];
            eps2 = e->head[0]->eps;
            e->pq00.Update(-eps - eps2);
        }
    }
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) r->tree->eps = 0;
}
bool PerfectMatching::UpdateDuals() {
    Node *r;
    double start_time = get_time();
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        Tree *t = r->tree;
        PriorityQueue<REAL>::Item *q;
        REAL eps = PM_INFTY;
        if ((q = t->pq0.GetMin())) eps = q->slack;
        if ((q = t->pq_blossoms.GetMin()) && eps > q->slack) eps = q->slack;
        while ((q = t->pq00.GetMin())) {
            if (ProcessEdge00((Edge *) q, false)) break;
            t->pq00.Remove(q, pq_buf);
        }
        if (q && 2 * eps > q->slack) eps = q->slack / 2;
        t->eps_delta = eps - t->eps;
    }
    if (tree_num >= options.dual_LP_threshold * node_num) {
        if (options.dual_greedy_update_option == 0) ComputeEpsCC();
        else if (options.dual_greedy_update_option == 1) ComputeEpsSCC();
        else ComputeEpsSingle();
    } else ComputeEpsGlobal();
    REAL delta = 0;
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        if (r->tree->eps_delta > 0) {
            delta += r->tree->eps_delta;
            r->tree->eps += r->tree->eps_delta;
        }
    }
    stat.dual_time += get_time() - start_time;
    return (delta > PM_THRESHOLD);
}
template<typename FlowType, typename CostType>
MinCost<FlowType, CostType>::MinCost(int _nodeNum, int _edgeNumMax, void (*err_function)(const char *))
        : nodeNum(_nodeNum),
          edgeNum(0),
          edgeNumMax(_edgeNumMax),
          counter(0),
          cost(0),
          error_function(err_function) {
    nodes = (Node *) malloc(nodeNum * sizeof(Node));
    arcs = (Arc *) malloc(2 * edgeNumMax * sizeof(Arc));
    if (!nodes || !arcs) {
        if (error_function) (*error_function)("Not enough memory!");
        exit(1);
    }
    memset(nodes, 0, nodeNum * sizeof(Node));
    memset(arcs, 0, 2 * edgeNumMax * sizeof(Arc));
    firstActive = &nodes[nodeNum];
#ifdef MINCOST_DEBUG
    for (int i=0; i<nodeNum; i++) nodes[i].id = i;
#endif
}
template<typename FlowType, typename CostType>
MinCost<FlowType, CostType>::~MinCost() {
    free(nodes);
    free(arcs);
}
template<typename FlowType, typename CostType>
void MinCost<FlowType, CostType>::Init() {
    Node *i;
    Arc *a;
    for (a = arcs; a < arcs + 2 * edgeNum; a++) {
        if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
    }
    Node **lastActivePtr = &firstActive;
    for (i = nodes; i < nodes + nodeNum; i++) {
        if (i->excess > 0) {
            *lastActivePtr = i;
            lastActivePtr = &i->next;
        } else i->next = NULL;
    }
    *lastActivePtr = &nodes[nodeNum];
}
template<typename FlowType, typename CostType>
FlowType MinCost<FlowType, CostType>::Augment(Node *start, Node *end) {
    FlowType delta = (start->excess < -end->excess) ? start->excess : -end->excess;
    Arc *a;
    for (a = end->parent; a; a = a->sister->head->parent) {
        if (delta > a->r_cap) delta = a->r_cap;
    }
    assert(delta > 0);
    end->excess += delta;
    for (a = end->parent; a; a = a->head->parent) {
        DecreaseRCap(a, delta);
        a = a->sister;
        IncreaseRCap(a, delta);
    }
    start->excess -= delta;
    return delta;
}
template<typename FlowType, typename CostType>
void MinCost<FlowType, CostType>::Dijkstra(Node *start) {
    assert(start->excess > 0);
    Node *i;
    Node *j;
    Arc *a;
    CostType d;
    Node *permanentNodes;
    int FLAG0 = ++counter; 
    int FLAG1 = ++counter; 
    start->parent = NULL;
    start->flag = FLAG1;
    queue.Reset();
    queue.Add(start, 0);
    permanentNodes = NULL;
    while ((i = queue.RemoveMin(d))) {
        if (i->excess < 0) {
            FlowType delta = Augment(start, i);
            cost += delta * (d - i->pi + start->pi);
            for (i = permanentNodes; i; i = i->next_permanent) i->pi += d;
            break;
        }
        i->pi -= d;
        i->flag = FLAG0;
        i->next_permanent = permanentNodes;
        permanentNodes = i;
        for (a = i->firstNonsaturated; a; a = a->next) {
            j = a->head;
            if (j->flag == FLAG0) continue;
            d = a->GetRCost();
            if (j->flag == FLAG1) {
                if (d >= queue.GetKey(j)) continue;
                queue.DecreaseKey(j, d);
            } else {
                queue.Add(j, d);
                j->flag = FLAG1;
            }
            j->parent = a;
        }
    }
}
template<typename FlowType, typename CostType>
CostType MinCost<FlowType, CostType>::Solve() {
    Node *i;
    
    while (1) {
        i = firstActive;
        if (i == &nodes[nodeNum]) break;
        firstActive = i->next;
        i->next = NULL;
        if (i->excess > 0) {
            Dijkstra(i);
            if (i->excess > 0 && !i->next) {
                i->next = firstActive;
                firstActive = i;
            }
        }
    }
#ifdef MINCOST_DEBUG
                                                                                                                            TestOptimality();
	TestCosts();
#endif
    return cost;
}
template<typename FlowType, typename CostType>
void MinCost<FlowType, CostType>::TestOptimality() {
    Node *i;
    Arc *a;
    for (i = nodes; i < nodes + nodeNum; i++) {
        if (i->excess != 0) {
            assert(0);
        }
        for (a = i->firstSaturated; a; a = a->next) {
            if (a->r_cap != 0) {
                assert(0);
            }
        }
        for (a = i->firstNonsaturated; a; a = a->next) {
            CostType c = a->GetRCost();
            if (a->r_cap <= 0 || a->GetRCost() < -1e-5) {
                assert(0);
            }
        }
    }
}
#ifdef MINCOST_DEBUG
                                                                                                                        template <typename FlowType, typename CostType>
	void MinCost<FlowType, CostType>::TestCosts()
{
	Arc* a;
	CostType _cost = 0;
	for (a=arcs; a<arcs+2*edgeNum; a+=2)
	{
		assert(a->r_cap + a->sister->r_cap == a->cap_orig + a->sister->cap_orig);
		_cost += a->cost*(a->cap_orig - a->r_cap);
	}
	CostType delta = cost - _cost;
	if (delta < 0) delta = -delta;
	if (delta >= 1e-5)
	{
		assert(0);
	}
}
#endif
#define FLOW_INFTY ((int)0x00fffffff)
template<typename CostType>
DualMinCost<CostType>::DualMinCost(int _nodeNum, int _edgeNumMax)
        : MinCost<int, CostType>(_nodeNum + 1, _edgeNumMax + 2 * _nodeNum) {
    source = _nodeNum;
}
template<typename CostType>
DualMinCost<CostType>::~DualMinCost() {
}
template<typename CostType>
void DualMinCost<CostType>::AddUnaryTerm(NodeId i, int objective_coef) {
    MinCost<int, CostType>::AddNodeExcess(i, objective_coef);
    MinCost<int, CostType>::AddNodeExcess(source, -objective_coef);
}
template<typename CostType>
void DualMinCost<CostType>::SetLowerBound(NodeId i, CostType cmin) {
    DualMinCost<CostType>::AddEdge(i, source, FLOW_INFTY, 0, -cmin);
}
template<typename CostType>
void DualMinCost<CostType>::SetUpperBound(NodeId i, CostType cmax) {
    DualMinCost<CostType>::AddEdge(source, i, FLOW_INFTY, 0, cmax);
}
template<typename CostType>
void DualMinCost<CostType>::AddConstraint(NodeId i, NodeId j, CostType cmax) {
    DualMinCost<CostType>::AddEdge(i, j, FLOW_INFTY, 0, cmax);
}
template<typename CostType>
void DualMinCost<CostType>::Solve() {
    MinCost<int, CostType>::Solve();
}
template<typename CostType>
CostType DualMinCost<CostType>::GetSolution(NodeId i) {
    return MinCost<int, CostType>::nodes[source].pi - MinCost<int, CostType>::nodes[i].pi;
}
#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif
template class MinCost<int, int>;
template class MinCost<int, double>;
template class DualMinCost<int>;
template class DualMinCost<double>;
GeomPerfectMatching::GeomPerfectMatching(int nodeNum, int _DIM)
        : DIM(_DIM),
          node_num(0),
          node_num_max(nodeNum),
          edge_num(0) {
    if (node_num_max < 1) {
        printf("too few nodes\n");
        exit(1);
    }
    if (node_num_max & 1) {
        printf("# of points is odd: perfect matching cannot exist\n");
        exit(1);
    }
    nodes = (Node *) malloc(node_num_max * sizeof(Node));
    memset(nodes, 0, node_num_max * sizeof(Node));
    edges = new Block<Edge>(512);
    coords = (REAL *) malloc((DIM + 1) * node_num_max * sizeof(REAL));
    sums = coords + DIM * node_num_max;
    matching = (int *) malloc(node_num_max * sizeof(int));
    int i;
    for (i = 0; i < node_num_max; i++) matching[i] = i;
}
GeomPerfectMatching::~GeomPerfectMatching() {
    free(nodes);
    delete edges;
    free(coords);
    free(matching);
}
GeomPerfectMatching::PointId GeomPerfectMatching::AddPoint(REAL *coord) {
    if (node_num >= node_num_max) {
        printf("Error: you are trying to add too many points!\n");
        exit(1);
    }
    memcpy(coords + DIM * node_num, coord, DIM * sizeof(REAL));
    return node_num++;
}
void GeomPerfectMatching::AddInitialEdge(PointId _i, PointId _j) {
    assert(_i >= 0 && _i < node_num_max && _j >= 0 && _j < node_num_max && _i != _j);
    if (_j < _i) {
        int _k = _i;
        _i = _j;
        _j = _k;
    }
    Node *i = nodes + _i;
    Node *j = nodes + _j;
    Edge *e = edges->New();
    edge_num++;
    e->head[1] = _i;
    e->head[0] = _j;
    e->next[0] = i->first[0];
    e->next[1] = j->first[1];
    i->first[0] = e;
    j->first[1] = e;
}
GeomPerfectMatching::REAL GeomPerfectMatching::ComputeCost(PointId *matching) {
    if (node_num != node_num_max) {
        printf("ComputeCost() cannot be called before all points have been added!\n");
        exit(1);
    }
    REAL cost = 0;
    int i;
    for (i = 0; i < node_num; i++) {
        if (matching[i] == i || matching[i] < 0 || matching[i] >= node_num || matching[matching[i]] != i) {
            printf("ComputeCost(): not a valid matching!\n");
            exit(1);
        }
        if (matching[i] > i) {
            cost += Dist(i, matching[i]);
        }
    }
    return cost;
}
template<typename Type>
inline void sort(Type *array, int array_skip, int *mapping, int N) {
    
    int i;
    int *mappingSrc = mapping;
    int *mappingDst = mapping + N;
    int *pSrc1;
    int *pSrc2;
    int *pSrc1End;
    int *pSrc2End;
    int *pDst;
    for (i = 0; i < (N & (~1)); i += 2) {
        if (array[array_skip * i] < array[array_skip * (i + 1)]) {
            mappingSrc[i] = i;
            mappingSrc[i + 1] = i + 1;
        } else {
            mappingSrc[i] = i + 1;
            mappingSrc[i + 1] = i;
        }
    }
    if (i != N) mappingSrc[i] = i;
    int step;
    for (step = 2; step < N; step *= 2) {
        pSrc2End = mappingSrc;
        pDst = mappingDst;
        while (1) {
            pSrc1 = pSrc2End;
            pSrc1End = pSrc1 + step;
            if (pSrc1End >= mappingSrc + N) {
                memcpy(pDst, pSrc1, (int) ((char *) (mappingSrc + N) - (char *) pSrc1));
                break;
            }
            pSrc2 = pSrc1End;
            pSrc2End = pSrc2 + step;
            if (pSrc2End > mappingSrc + N) pSrc2End = mappingSrc + N;
            while (1) {
                if (array[(array_skip) * (*pSrc1)] < array[array_skip * (*pSrc2)]) {
                    *pDst++ = *pSrc1++;
                    if (pSrc1 == pSrc1End) {
                        memcpy(pDst, pSrc2, (int) ((char *) pSrc2End - (char *) pSrc2));
                        pDst = (int *) ((char *) pDst + (int) ((char *) pSrc2End - (char *) pSrc2));
                        break;
                    }
                } else {
                    *pDst++ = *pSrc2++;
                    if (pSrc2 == pSrc2End) {
                        memcpy(pDst, pSrc1, (int) ((char *) pSrc1End - (char *) pSrc1));
                        pDst = (int *) ((char *) pDst + (int) ((char *) pSrc1End - (char *) pSrc1));
                        break;
                    }
                }
            }
        }
        pDst = mappingDst;
        mappingDst = mappingSrc;
        mappingSrc = pDst;
    }
    if (mappingSrc != mapping) memcpy(mapping, mappingSrc, N * sizeof(int));
}
#define NEIGHBOR_PARENT(k) (((k)-1)>>1)
#define NEIGHBOR_FIRST_CHILD(k) (((k)<<1)+1)
Neighbors::Neighbors() {
    K_max = 0;
    dist_array = NULL;
}
Neighbors::~Neighbors() {
    if (dist_array) delete[] dist_array;
}
void Neighbors::Init(int _K, PointId *_array) {
    K = _K;
    array = _array;
    num = 0;
    if (K > K_max) {
        if (dist_array) delete[] dist_array;
        K_max = K;
        dist_array = new double[K_max];
    }
}
inline void Neighbors::Swap(int k1, int k2) {
    PointId p = array[k1];
    array[k1] = array[k2];
    array[k2] = p;
    double d = dist_array[k1];
    dist_array[k1] = dist_array[k2];
    dist_array[k2] = d;
}
inline void Neighbors::Add(PointId p, double dist) {
    int k;
    if (num < K) {
        k = num++;
        array[k] = p;
        dist_array[k] = dist;
        while (k > 0) {
            int k_parent = NEIGHBOR_PARENT(k);
            if (dist_array[k] <= dist_array[k_parent]) break;
            Swap(k, k_parent);
            k = k_parent;
        }
    } else {
        if (dist_array[0] <= dist) return;
        array[0] = p;
        dist_array[0] = dist;
        k = 0;
        while (1) {
            int k_child = NEIGHBOR_FIRST_CHILD(k);
            if (k_child >= K) break;
            if (k_child + 1 < K && dist_array[k_child + 1] > dist_array[k_child]) k_child++;
            if (dist_array[k] >= dist_array[k_child]) break;
            Swap(k, k_child);
            k = k_child;
        }
    }
    
}
inline double Neighbors::GetMax() {
    assert(num > 0);
    return dist_array[0];
}
GPMKDTree::GPMKDTree(int _D, int _point_num, REAL *coords, GeomPerfectMatching *_GPM)
        : D(_D), DIM(_GPM->DIM), point_num(_point_num), GPM(_GPM) {
    Node *i;
    Node *j;
    int d, d0, k;
    int *mapping = new int[(D + 2) * point_num];
    int *buf = mapping + D * point_num;
    int *marking = buf + point_num;
    memset(marking, 0, point_num * sizeof(int));
    int *ptr = mapping;
    int node_num_max = 4 * point_num / 3 + 2;
    nodes = (Node *) malloc(node_num_max * sizeof(Node));
    rev_mapping = (Node **) malloc(point_num * sizeof(Node *));
    memset(rev_mapping, 0, point_num * sizeof(Node *));
    REAL **coords_array = new REAL *[D];
    int *skip_array = new int[D];
    for (d = 0; d < DIM; d++) {
        coords_array[d] = coords + d;
        skip_array[d] = DIM;
    }
    if (d < D) {
        coords_array[d] = GPM->sums;
        skip_array[d] = 1;
    }
    for (d = 0; d < D; d++) {
        sort<REAL>(coords_array[d], skip_array[d], ptr, point_num);
        if (d == DIM) sum_max = GPM->sums[ptr[point_num - 1]];
        ptr += point_num;
    }
    nodes[0].parent = NULL;
    nodes[0].order = 0;
    nodes[0].d = point_num;
    nodes[0].first_child = NULL;
    node_num = 1;
    Node *first_unprocessed = &nodes[0];
    while ((i = first_unprocessed)) {
        first_unprocessed = i->first_child;
        int start = i->order;
        int num0 = i->d, num;
        if ((DIM == D && num0 <= 2) || (DIM < D && num0 <= CHILDREN_MAX)) {
            
            i->d = -num0;
            for (k = 0; k < num0; k++) {
                i->points[k] = mapping[start + k];
                rev_mapping[mapping[start + k]] = i;
            }
            continue;
        }
        
        if (node_num + 2 > node_num_max) {
            node_num_max = 3 * node_num_max / 2 + 16;
            Node *nodes_old = nodes;
            nodes = (Node *) realloc(nodes, node_num_max * sizeof(Node));
#define UPDATE_NODE_PTR(ptr) ptr = (Node*)((char*)ptr + ((char*)nodes-(char*)nodes_old))
            UPDATE_NODE_PTR(i);
            if (first_unprocessed) UPDATE_NODE_PTR(first_unprocessed);
            for (k = 0; k < node_num; k++) {
                if (nodes[k].parent) UPDATE_NODE_PTR(nodes[k].parent);
                if (nodes[k].d >= 0 && nodes[k].first_child) UPDATE_NODE_PTR(nodes[k].first_child);
            }
            for (k = 0; k < point_num; k++) if (rev_mapping[k]) UPDATE_NODE_PTR(rev_mapping[k]);
        }
        
        
        
        
        
        const int FRACTION = 20;
        int num_min = 1;
        if (num_min < num0 / FRACTION) num_min = num0 / FRACTION;
        int num_max = num0 - 1;
        if (num_max > (FRACTION - 1) * num0 / FRACTION) num_max = (FRACTION - 1) * num0 / FRACTION;
        int num_max_DIM = num0 - 1;
        if (D > DIM && (!i->parent || !i->parent->parent || !i->parent->parent->parent)) d0 = D - 1;
        else {
            d0 = -1;
            REAL diff_max = 0;
            for (d = 0; d < D; d++) {
                ptr = mapping + d * point_num;
                REAL c_min, c_max;
                c_min = coords_array[d][ptr[start + num_min] * skip_array[d]];
                c_max = coords_array[d][ptr[start + ((d < DIM) ? num_max : num_max_DIM)] * skip_array[d]];
                REAL diff = c_max - c_min;
                
                if (d0 < 0 || diff_max < diff) {
                    d0 = d;
                    diff_max = diff;
                }
            }
        }
        ptr = mapping + d0 * point_num;
        REAL c_min, c_max;
        if (d < DIM) {
            c_min = coords_array[d0][ptr[start + num_min] * skip_array[d0]];
            c_max = coords_array[d0][ptr[start + num_max] * skip_array[d0]];
        } else {
            c_min = coords_array[d0][ptr[start + 1] * skip_array[d0]];
            c_max = coords_array[d0][ptr[start + num0 - 1] *
                                     skip_array[d0]]; 
        }
        REAL split = (c_min + c_max) / 2;
        for (num = num_min; num < num_max; num++) {
            if (coords_array[d0][ptr[start + num] * skip_array[d0]] > split) break;
        }
        
        
        ptr = mapping + d0 * point_num;
        for (k = start; k < start + num; k++) marking[ptr[k]] = 1;
        for (d = 0; d < D; d++) {
            if (d == d0) continue;
            ptr = mapping + d * point_num;
            int *bufa = buf;
            int *bufb = buf + num;
            for (k = start; k < start + num0; k++) {
                if (marking[ptr[k]]) *bufa++ = ptr[k];
                else *bufb++ = ptr[k];
            }
            memcpy(ptr + start, buf, num0 * sizeof(int));
        }
        ptr = mapping + d0 * point_num;
        for (k = start; k < start + num; k++) marking[ptr[k]] = 0;
        i->d = d0;
        PointId p = ptr[start + num - ((num >= num0 - num) ? 1 : 0)];
        i->coord = coords_array[d0][p * skip_array[d0]];
        i->first_child = j = &nodes[node_num];
        node_num += 2;
        j->parent = (j + 1)->parent = i;
        j->order = start;
        j->d = num;
        (j + 1)->order = start + num;
        (j + 1)->d = num0 - num;
        (j + 1)->first_child = first_unprocessed;
        j->first_child = j + 1;
        first_unprocessed = j;
    }
    delete[] coords_array;
    delete[] skip_array;
    delete[] mapping;
    
    
    int depth = 0, depth_max = 0;
    
    i = &nodes[0];
    k = 0;
    if (D > DIM) {
        
        while (1) {
            if (!IS_LEAF(i)) {
                i = i->first_child;
                depth++;
                if (depth_max < depth) depth_max = depth;
            } else {
                while (1) {
                    i->order = k++;
                    if (!i->parent) break;
                    if (i->parent->first_child == i) {
                        i++;
                        break;
                    }
                    i = i->parent;
                    depth--;
                }
                if (!i->parent) break;
            }
        }
    } else {
        
        while (1) {
            if (!IS_LEAF(i)) {
                i = i->first_child;
                depth++;
                if (depth_max < depth) depth_max = depth;
            } else {
                i->order = k++;
                while (1) {
                    if (!i->parent) break;
                    if (i->parent->first_child == i) {
                        i->parent->order = k++;
                        i++;
                        break;
                    }
                    i = i->parent;
                    depth--;
                }
                if (!i->parent) break;
            }
        }
    }
    traversing_buf = (REAL *) malloc((D + depth_max + 2) * sizeof(REAL));
}
GPMKDTree::~GPMKDTree() {
    free(nodes);
    free(rev_mapping);
    free(traversing_buf);
}
void GPMKDTree::AddPerfectMatching(PointId *rev_mapping) {
    Node *i;
    int k;
    PointId p, q = -1;
    i = &nodes[0];
    do {
        if (IS_LEAF(i)) {
            for (k = 0; k < -i->d; k++) {
                p = i->points[k];
                if (q < 0) q = p;
                else {
                    GPM->AddInitialEdge(rev_mapping[p], rev_mapping[q]);
                    q = -1;
                }
            }
        } else {
            i = i->first_child;
            continue;
        }
        while (i->parent) {
            if (i->parent->first_child == i) {
                i++;
                break;
            }
            i = i->parent;
        }
    } while (i->parent);
}
#define MOVE_DOWN_LEFT(i)    {        *stack ++ = current_diff[i->d];        if (current_diff[i->d] <= 0 && (diff = i->coord - coord0[i->d]) < 0) current_diff[i->d] = diff;        i = i->first_child;    }
#define MOVE_DOWN_RIGHT(i)    {        *stack ++ = current_diff[i->d];        if (current_diff[i->d] >= 0 && (diff = i->coord - coord0[i->d]) > 0) current_diff[i->d] = diff;        i = i->first_child+1;    }
#define MOVE_LEFT(i)    {        int d_prev = i->parent->d;        current_diff[d_prev] = stack[-1];        if (current_diff[d_prev] <= 0 && (diff = i->parent->coord - coord0[d_prev]) < 0) current_diff[d_prev] = diff;        i --;    }
#define MOVE_RIGHT(i)    {        int d_prev = i->parent->d;        current_diff[d_prev] = stack[-1];        if (current_diff[d_prev] >= 0 && (diff = i->parent->coord - coord0[d_prev]) > 0) current_diff[d_prev] = diff;        i ++;    }
#define MOVE_UP(i)    {        i = i->parent;        current_diff[i->d] = *(-- stack);    }
#define MOVE_DOWN_LEFT_X(i)    {        *stack ++ = current_diff[i->d];        if (i->d < D-1) { if (current_diff[i->d] <= 0 && (diff = i->coord - coord0[i->d]) < 0) current_diff[i->d] = diff; }        else current_diff[i->d] = i->coord;        i = i->first_child;    }
#define MOVE_DOWN_RIGHT_X(i)    {        *stack ++ = current_diff[i->d];        if (i->d < D-1) { if (current_diff[i->d] >= 0 && (diff = i->coord - coord0[i->d]) > 0) current_diff[i->d] = diff; }        i = i->first_child+1;    }
#define MOVE_LEFT_X(i)    {        int d_prev = i->parent->d;        current_diff[d_prev] = stack[-1];        if (d_prev < D-1) { if (current_diff[d_prev] <= 0 && (diff = i->parent->coord - coord0[d_prev]) < 0) current_diff[d_prev] = diff; }        else current_diff[d_prev] = i->parent->coord;        i --;    }
#define MOVE_RIGHT_X(i)    {        int d_prev = i->parent->d;        current_diff[d_prev] = stack[-1];        if (d_prev < D-1) { if (current_diff[d_prev] >= 0 && (diff = i->parent->coord - coord0[d_prev]) > 0) current_diff[d_prev] = diff; }        i ++;    }
#define MOVE_UP_X(i)    {        i = i->parent;        current_diff[i->d] = *(-- stack);    }
/*
void GPMKDTree::ComputeKNN(PointId p0, int K, PointId* neighbors_array)
{
	int p;
	REAL* coord0 = GPM->coords + p0*DIM;
	neighbors.Init(K, neighbors_array);
	for (p=0; p<point_num; p++)
	if (p != p0)
	{
		neighbors.Add(p, GPM->Dist(coord0, GPM->coords+p*DIM));
	}
}
*/
void GPMKDTree::ComputeKNN(PointId p0, int K, PointId *neighbors_array) {
    int neighbor_num = 0;
    Node *i = rev_mapping[p0];
    int k, order0 = i->order;
    REAL diff;
    REAL *coords = GPM->coords;
    REAL *coord0 = GPM->coords + p0 * DIM;
    REAL *current_diff = traversing_buf;
    REAL *stack = traversing_buf + D;
    neighbors.Init(K, neighbors_array);
    for (k = 0; k < D; k++) current_diff[k] = 0;
    i = &nodes[0];
    do {
        if (IS_LEAF(i)) {
            for (k = 0; k < -i->d; k++) {
                PointId p = i->points[k];
                if (p == p0) continue;
                double dist2;
                REAL *coord2 = coords + p * DIM;
                GPM_GET_DIST2(dist2, coord0, coord2);
                neighbors.Add(p, dist2);
            }
        } else {
            if (neighbors.GetNum() < K || GPM->Norm2(current_diff) < neighbors.GetMax()) {
                if (i->order > order0) {
                    MOVE_DOWN_LEFT(i);
                } else {
                    MOVE_DOWN_RIGHT(i);
                }
                continue;
            }
        }
        while (i->parent) {
            if (i->parent->order > order0) {
                if (i->parent->first_child == i) {
                    MOVE_RIGHT(i);
                    break;
                }
            } else {
                if (i->parent->first_child != i) {
                    MOVE_LEFT(i);
                    break;
                }
            }
            MOVE_UP(i);
        }
    } while (i->parent);
}
/*
void GPMKDTree::AddNegativeEdges(PointId p, PerfectMatching* pm)
{
	PointId q;
	for (q=p+1; q<node_num; q++)
	{
		if (GPM->nodes[q].is_marked) continue;
		REAL len = GPM->Dist(p, q);
		if (2*len - GPM->sums[p] - GPM->sums[q] < 0)
		{
			if (pm->AddNewEdge(p, q, len, true)>=0)
			{
				GPM->AddInitialEdge(p, q);
			}
		}
	}
}
*/
void GPMKDTree::AddNegativeEdges(PointId p0, PerfectMatching *pm) {
    Node *i = rev_mapping[p0];
    int k, order0 = i->order;
    bool check;
    REAL diff;
    REAL *coords = GPM->coords;
    REAL *coord0 = coords + p0 * DIM;
    REAL *sums = GPM->sums;
    REAL sum0 = sums[p0];
    REAL *current_diff = traversing_buf;
    REAL *stack = traversing_buf + D;
    for (k = 0; i->points[k] != p0; k++) {}
    for (k++; k < -i->d; k++) {
        PointId p = i->points[k];
        REAL len = GPM->Dist(coord0, GPM->coords + p * DIM);
        if (2 * len - GPM->sums[p] < GPM->sums[p0]) {
            double start_time = get_time();
            if (pm->AddNewEdge(p0, p, len, true) >= 0) GPM->AddInitialEdge(p0, p);
            GPM->graph_update_time += get_time() - start_time;
        }
    }
    for (k = 0; k < DIM; k++) current_diff[k] = 0;
    current_diff[k] = sum_max;
    i = &nodes[0];
    do {
        if (i->order > order0) {
            if (IS_LEAF(i)) {
                for (k = 0; k < -i->d; k++) {
                    PointId p = i->points[k];
                    if (!GPM->nodes[p].is_marked) {
                        
                        
                        
                        
                        {
                            REAL len = GPM->Dist(coord0, GPM->coords + p * DIM);
                            if (2 * len - GPM->sums[p] < GPM->sums[p0]) {
                                double start_time = get_time();
                                if (pm->AddNewEdge(p0, p, len, true) >= 0) GPM->AddInitialEdge(p0, p);
                                GPM->graph_update_time += get_time() - start_time;
                            }
                        }
                    }
                }
            } else {
                REAL threshold = current_diff[D - 1] + sum0;
                GPM_CHECK_NORM(check, current_diff, threshold);
                if (check) {
                    MOVE_DOWN_LEFT_X(i);
                    continue;
                }
            }
        }
        while (i->parent) {
            if (i->parent->first_child == i) {
                MOVE_RIGHT_X(i);
                break;
            }
            MOVE_UP_X(i);
        }
    } while (i->parent);
}
GeomPerfectMatching::REAL GeomPerfectMatching::SolveComplete() {
    if (node_num != node_num_max) {
        printf("ComputeCost() cannot be called before all points have been added!\n");
        exit(1);
    }
    PointId p, q;
    int e = 0, E = node_num * (node_num - 1) / 2;
    PerfectMatching *pm = new PerfectMatching(node_num, E);
    for (p = 0; p < node_num; p++) {
        for (q = p + 1; q < node_num; q++) {
            pm->AddEdge(p, q, Dist(p, q));
        }
    }
    pm->options = options;
    pm->Solve();
    for (p = 0; p < node_num; p++) {
        for (q = p + 1; q < node_num; q++) {
            if (pm->GetSolution(e++)) {
                matching[p] = q;
                matching[q] = p;
            }
        }
    }
    delete pm;
    return ComputeCost(matching);
}
GeomPerfectMatching::REAL GeomPerfectMatching::Solve() {
    double start_time = get_time();
    double perfect_matching_time = 0;
    double negative_edges_time = 0;
    if (options.verbose) {
        printf("starting geometric matching with %d points\n", node_num);
        fflush(stdout);
    }
    PointId p, q;
    Edge *e;
    int _e;
    int iter;
    bool success = false;
    PerfectMatching *pm = NULL;
    GPMKDTree *kd_tree;
    double init_matching_time = get_time();
    if (gpm_options.init_Delaunay) InitDelaunay();
    if (gpm_options.init_KNN > 0) InitKNN(gpm_options.init_KNN);
    if (gpm_options.init_greedy) CompleteInitialMatching();
    init_matching_time = get_time() - init_matching_time;
    graph_update_time = 0;
    int iter_max = gpm_options.iter_max;
    for (iter = 0; iter_max <= 0 || iter < iter_max; iter++) {
        if (pm) {
            double negative_edges_start_time = get_time();
            int edge_num0 = edge_num;
            pm->StartUpdate();
            for (p = 0; p < node_num; p++) {
                PerfectMatching::REAL s = pm->GetTwiceSum(p);
                if (((REAL) 1 / 2) == 0 && ((PerfectMatching::REAL) 1 / 2) != 0) sums[p] = (REAL) ceil((double) s);
                else sums[p] = (REAL) s;
            }
            if (options.verbose) {
                printf("building kd_tree...");
                fflush(stdout);
            }
            {
                kd_tree = new GPMKDTree(DIM + 1, node_num, coords, this);
            }
            if (options.verbose) {
                printf(" done. Now adding negative edges:\n    ");
                fflush(stdout);
            }
            for (p = 0; p < node_num; p++) {
                if (options.verbose && (p % (node_num / 72) == 0)) {
                    printf("+");
                    fflush(stdout);
                }
                for (e = nodes[p].first[0]; e; e = e->next[0]) nodes[e->head[0]].is_marked = 1;
                for (e = nodes[p].first[1]; e; e = e->next[1]) nodes[e->head[1]].is_marked = 1;
                kd_tree->AddNegativeEdges(p, pm);
                for (e = nodes[p].first[0]; e; e = e->next[0]) nodes[e->head[0]].is_marked = 0;
                for (e = nodes[p].first[1]; e; e = e->next[1]) nodes[e->head[1]].is_marked = 0;
            }
            delete kd_tree;
            
            if (0) 
            {
                delete pm;
                pm = NULL;
            } else {
                pm->FinishUpdate();
                if (edge_num0 == edge_num) success = true;
            }
            if (options.verbose) {
                printf("\ndone (%d edges added)\n", edge_num - edge_num0);
                fflush(stdout);
            }
            negative_edges_time += get_time() - negative_edges_start_time;
        }
        if (!pm) {
            int E = 5 * node_num;
            if (E < 5 * edge_num / 4) E = 5 * edge_num / 4;
            pm = new PerfectMatching(node_num, E);
            for (e = edges->ScanFirst(); e; e = edges->ScanNext()) {
                p = e->head[1];
                q = e->head[0];
                pm->AddEdge(p, q, Dist(p, q));
            }
        }
        if (options.verbose) printf("iter %d: ", iter + 1);
        pm->options = options;
        double perfect_matching_start = get_time();
        pm->Solve();
        perfect_matching_time += get_time() - perfect_matching_start;
        if (success) break;
    }
    for (_e = 0, e = edges->ScanFirst(); e; _e++, e = edges->ScanNext()) {
        if (pm->GetSolution(_e)) {
            p = e->head[1];
            q = e->head[0];
            matching[p] = q;
            matching[q] = p;
        }
    }
    delete pm;
    REAL cost = ComputeCost(matching);
    if (options.verbose) {
        printf("geometric matching finished [%.3f secs]. cost=%.1f \n", get_time() - start_time, (double) cost);
        printf("    selecting initial edges: [%.3f secs], perfect matching: [%.3f secs]\n", init_matching_time,
               perfect_matching_time);
        printf("    pricing: [%.3f secs] including graph updates: [%.3f secs]\n", negative_edges_time,
               graph_update_time);
        fflush(stdout);
    }
    return cost;
}
#endif
