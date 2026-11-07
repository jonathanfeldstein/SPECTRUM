#ifndef SPECTRUM_FRAGMENTS_H
#define SPECTRUM_FRAGMENTS_H

#include <iostream>
#include <map>
#include <vector>
#include <boost/container/small_vector.hpp>
#include <set>
#include <unordered_set>
#include <string>
#include <cctype>
#include <list>
#include <algorithm>
#include <random>
#include "../HyperGraphs/HyperGraph.h"
#include "fragment_utils_dispatch.h"
#include "../Utils/vector_utils.h"
#include "../Utils/profiling_utils.h"

#include <omp.h>

using namespace std;
using NodeId = size_t;
using NodeName = string;
using EdgeId = size_t;
using PredicateId = size_t;
using Predicate = string;
using Fragment = std::vector<PredicateId>;
using Location = std::vector<EdgeId>;

using SmallFragment =
    boost::container::small_vector<PredicateId,8>;

using SmallLocation =
    boost::container::small_vector<EdgeId,8>;

class FragmentCollection{
public:
    unordered_map<Fragment, vector<Location>, VectorHash> fragment_count;
    unordered_map<Fragment, vector<Location>, VectorHash> fragment_cyclic_count;
private:
    bool has_unaries_;
    unordered_map<NodeId, pair<vector<EdgeId>,vector<EdgeId>>> unary_binary_edges;
    unordered_map<NodeId, vector<EdgeId>> unary_edges;
    unordered_map<NodeId, vector<EdgeId>> binary_edges;
    unordered_map<NodeId,unordered_map<NodeId, vector<EdgeId>>> shared_edges;

    unordered_map<NodeId,unordered_set<NodeId>> connected_nodes;
    
    unordered_map<EdgeId, PredicateId> edge_to_predicate;
    unordered_set<PredicateId> unary_predicates;

public:
    HyperGraph hg;
    FragmentCollection(const HyperGraph& hg, size_t number_of_random_walks, size_t max_depth, bool force_cyclic, bool has_unaries);
    FragmentCollection();

    // Fragment storage (vector and single-pair overloads)
    void storing_fragments(
        const vector<pair<SmallFragment, SmallLocation>>& fragments,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
        RecursiveProfiling& profile);

    void storing_fragments(
        const pair<SmallFragment, SmallLocation>& fragment,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
        RecursiveProfiling& profile);

    void storing_cyclic_fragments(
        const vector<pair<SmallFragment, SmallLocation>>& fragments,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count,
        RecursiveProfiling& profile);

    void storing_cyclic_fragments(
        const pair<SmallFragment, SmallLocation>& fragment,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count,
        RecursiveProfiling& profile);

    void canonicalize_fragments();

    // Fragment discovery
    void get_fragments(size_t number_of_random_walks, size_t max_depth, bool force_cyclic);

    void set_initial_fragments(
        NodeId node,
        const vector<EdgeId>& binaries,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count,
        RecursiveProfiling& recursive_profile);

    // Accessors
    unordered_map<Fragment, vector<Location>, VectorHash>& get_fragment_count();
    unordered_map<Fragment, vector<Location>, VectorHash>& get_fragment_cyclic_count();
    string get_predicate_from_id(PredicateId predicate);
    PredicateId get_predicate_id(Predicate predicate);

    // Graph helpers
    vector<EdgeId> find_cyclic_edges(
        NodeId node,
        NodeId starting_node,
        const unordered_set<NodeId>& connected_nodes,
        const boost::container::small_vector<EdgeId,3>& previous_edges,
        RecursiveProfiling& profile
    );

    vector<EdgeId> get_edges_for_next_step(
        NodeId node,
        const boost::container::small_vector<EdgeId,3>& previous_edges,
        const vector<EdgeId>& forbidden_edges,
        size_t number_of_random_walks,
        RecursiveProfiling& profile
    );

    // Edge selection helper (shared between next_step variants)
    vector<EdgeId> select_binaries_for_step(
        size_t number_of_random_walks,
        NodeId node,
        NodeId starting_node,
        int depth,
        boost::container::small_vector<EdgeId,3>& previous_edges,
        const unordered_set<NodeId>& connected_nodes,
        bool force_cyclic,
        RecursiveProfiling& profile
    );

    // Recursive walk steps (separate due to different fragment types)
    void next_step_no_unaries(
        size_t number_of_random_walks,
        NodeId node,
        NodeId starting_node,
        int depth,
        pair<SmallFragment, SmallLocation>& fragment,
        boost::container::small_vector<EdgeId,3>& previous_edges,
        const unordered_set<NodeId>& connected_nodes,
        size_t max_depth,
        bool force_cyclic,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count,
        RecursiveProfiling& profile);

    void next_step_with_unaries(
        size_t number_of_random_walks,
        NodeId node,
        NodeId starting_node,
        int depth,
        const vector<pair<SmallFragment, SmallLocation>>& fragments,
        boost::container::small_vector<EdgeId,3>& previous_edges,
        const unordered_set<NodeId>& connected_nodes,
        size_t max_depth,
        bool force_cyclic,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
        unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count,
        RecursiveProfiling& profile);
        
    bool is_unary(PredicateId predicate);
    bool is_unary_edge(EdgeId edge);
    // Mined fragments have unary slots only when the data has unary predicates
    bool unary_layout() const { return has_unaries_; }
    const vector<NodeId>& get_edge(EdgeId edge);
    pair<vector<EdgeId>,vector<EdgeId>> get_unary_binary_edges(NodeId node);
};
#endif //SPECTRUM_FRAGMENTS_H
