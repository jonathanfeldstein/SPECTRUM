#ifndef SPECTRUM_HYPERGRAPH_H
#define SPECTRUM_HYPERGRAPH_H

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
#include "../Utils/map_utils.h"
#include <fstream>
#include "../Utils/read_files_utils.h"
#include "Relation.h"
#include "../Utils/exceptions.h"
#include <stdexcept>
#include <cerrno>
#include <cstring>
#include "../Utils/set_utils.h"
#include "../Utils/vector_utils.h"

using namespace std;
using NodeId = size_t;
using EdgeId = size_t;
using PredicateId = size_t;
using Predicate = string;
using NodeName = string;


class HyperGraph{
private:
    // MEMBERS
    unordered_map<EdgeId, vector<NodeId>> edges; // edge_id : list(node_id)
    // node_id : node_type
    unordered_map<NodeId, vector<EdgeId>> memberships; // node_id : list(edge_id)
    // set(node_type)
    unordered_map<EdgeId, PredicateId> predicates;
    unordered_set<PredicateId> unary_predicates;
    unordered_map<PredicateId, Predicate> predicate_id_to_name;
    unordered_map<Predicate, PredicateId> predicate_name_to_id;
    unordered_map<NodeId,unordered_map<NodeId, vector<EdgeId>>> shared_edges;

    unordered_map<NodeId,unordered_set<NodeId>> connected_nodes;

protected:
    unordered_map<EdgeId, double> edge_weights; // edge_id : weight
    unordered_map<NodeId, NodeName> node_ids_names; // node_id : node_name
    unordered_map<NodeName, NodeId> node_names_ids;

public:
    HyperGraph();
    HyperGraph(string const& db_file_path, string const& info_file_path, bool safe);

    ~HyperGraph();

    void add_edge(EdgeId edge_id, PredicateId predicate, vector<NodeId> node_ids, double weight);
    unordered_map<EdgeId, vector<NodeId>>& get_edges();
    const vector<NodeId>& get_edge(EdgeId edge_id);
    unordered_set<NodeId> get_nodes();
    string_view get_predicate(EdgeId edge_id);
    unordered_map<EdgeId, PredicateId> get_predicates();
    unordered_map<PredicateId , Predicate> get_predicate_id_to_name();
    unordered_map<Predicate, PredicateId> get_predicate_name_to_id();
    unordered_map<NodeId, vector<EdgeId>> get_memberships();
    PredicateId get_predicate_id(EdgeId edge);
    PredicateId get_predicate_id(Predicate predicate);
    const vector<EdgeId>& get_memberships(NodeId node_id);
    size_t number_of_nodes(); //Number of nodes
    size_t number_of_edges(); //Number of edges
    double get_edge_weight(EdgeId edge_id);
    unordered_map<EdgeId, double> get_edge_weights();
    void print();
    bool is_unary(PredicateId predicate);
    bool is_unary_edge(EdgeId edge);
    bool has_unaries() const;
    string get_predicate_from_id(PredicateId predicate);
    pair<vector<size_t>, vector<size_t>> get_unary_binary_edges(NodeId node);
    const NodeId& get_next_node(NodeId start, EdgeId edge);
    const unordered_map<NodeId, vector<EdgeId>>& get_shared_edges(NodeId node);
    const unordered_set<NodeId>& get_connected_nodes(NodeId node);
    bool has_shared_edges(NodeId node1, NodeId node2);

};

#endif //SPECTRUM_HYPERGRAPH_H
