/*
 * Fragments_rw.h — Random Walks variant of FragmentCollection
 *
 * Identical to Fragments_rw.h except:
*    - next_step() is replaced by do_random_walk()
 *   - do_random_walk picks ONE random binary edge per step (true random walk)
 *   - do_random_walk tracks visited_edges to halt if backtracking revisits an edge.
 */

#ifndef SPECTRUM_FRAGMENTS_V3_H
#define SPECTRUM_FRAGMENTS_V3_H

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <unordered_set>
#include <string>
#include <cctype>
#include <list>
#include <algorithm>
#include <random>
#include "../HyperGraphs/HyperGraph.h"
#include "fragment_utils_dispatch.h"
#include "../Utils/map_utils.h"
#include "../Utils/vector_utils.h"
#include <omp.h>

using namespace std;
using NodeId = size_t;
using NodeName = string;
using EdgeId = size_t;
using PredicateId = size_t;
using Predicate = string;
using Location = vector<EdgeId>;
using Fragment = vector<PredicateId>;

class FragmentCollection{
public:
    unordered_map<Fragment, vector<Location>, VectorHash> fragment_count;
    unordered_map<Fragment, vector<Location>, VectorHash> fragment_cyclic_count;
private:
    unordered_map<NodeId, pair<vector<EdgeId>,vector<EdgeId>>> unary_binary_edges;
    unordered_map<EdgeId, PredicateId> edge_to_predicate;
    unordered_set<PredicateId> unary_predicates;

public:
    HyperGraph hg;
    FragmentCollection(const HyperGraph& hg, size_t number_of_random_walks, size_t max_depth, bool force_cyclic);
    FragmentCollection();
    void set_initial_fragments(NodeId node, vector<EdgeId> unaries, const vector<EdgeId>& binaries);
    void storing_fragments(const map<Fragment, Location>& fragments);
    void storing_cyclic_fragments(const map<Fragment, Location>& fragments);
    void get_fragments(size_t number_of_random_walks, size_t max_depth, bool force_cyclic);
    unordered_map<Fragment, vector<Location>, VectorHash>& get_fragment_count();
    unordered_map<Fragment, vector<Location>, VectorHash>& get_fragment_cyclic_count();
    void set_fragment_count(unordered_map<Fragment, vector<Location>, VectorHash> fragment_count);
    void set_fragment_cyclic_count(unordered_map<Fragment, vector<Location>, VectorHash> fragment_cyclic_count);
    string get_predicate_from_id(PredicateId predicate);
    PredicateId get_predicate_id(Predicate predicate);
    
    // Random walk step — replaces next_step for the RW variant
    void do_random_walk(std::mt19937& gen,
                       NodeId node,
                       NodeId starting_node,
                       int depth,
                       const map<vector<PredicateId>, Location> &fragments,
                       size_t max_depth,
                       std::unordered_set<EdgeId> visited_edges);

    bool is_unary(PredicateId predicate);
    bool is_unary_edge(EdgeId edge);
    // Random-walk fragments always alternate unary/binary slots, even on binary-only data
    bool unary_layout() const { return true; }
    vector<NodeId> get_edge(EdgeId edge);
    pair<vector<EdgeId>,vector<EdgeId>> get_unary_binary_edges(NodeId node);
};
#endif //SPECTRUM_FRAGMENTS_V3_H
