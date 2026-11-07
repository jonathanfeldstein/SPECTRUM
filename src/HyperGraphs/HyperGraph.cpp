#include "HyperGraph.h"

using namespace std;

std::size_t getUniqueID(const std::string& str) {
    std::hash<std::string> hash_fn;
    return hash_fn(str);
}

HyperGraph::HyperGraph() = default;

HyperGraph::HyperGraph(string const& db_file_path, string const& info_file_path, bool safe) {
    if(!file_exists(db_file_path)){
        throw FileNotFoundException(db_file_path);
    } else if(!file_exists(info_file_path)){
        throw FileNotFoundException(info_file_path);
    }else{
        fstream info_file;
        info_file.open(info_file_path, ios::in);
        if(info_file.is_open()) {
            string line;
            PredicateId predicate_id{1};  // Start predicate id from 1
            // initialize the predicate to id and id to predicate bimap

            while (getline(info_file, line)) {
                if ((line[0] == '/' && line[1] == '/') || line.empty()) { // Escaping commented lines and empty lines
                    continue;
                }
                // line is formatted like predicate_name(arg1, arg2, ...)
                // strip the line to extract the predicate name
                Relation relation= parse_line_info(line, safe);
                this->predicate_id_to_name[predicate_id] = relation.predicate;
                this->predicate_name_to_id[relation.predicate] = predicate_id;
                if(relation.arguments.size() == 1){
                    this->unary_predicates.insert(predicate_id);
                }
                predicate_id++;
            }
            info_file.close();
        }
        fstream db_file;
        db_file.open(db_file_path, ios::in);
        if(db_file.is_open()){
            string line;
            EdgeId edge_id{0};
            while(getline(db_file, line)){
                if((line[0] == '/' && line[1] == '/') || line.empty()){ // Escaping commented lines and empty lines
                    continue;
                }
                GroundRelation relation = parse_line_db(line, safe);
                vector<NodeId> node_ids_in_edge;
                for(auto &argument: relation.arguments){
                    if(!this->node_names_ids.count(argument)){
                        NodeId node_id = node_names_ids.size();
                        this->node_names_ids[argument] = node_id;
                        this->node_ids_names[node_id] = argument;
                    }
                    node_ids_in_edge.push_back(this->node_names_ids[argument]);
                }
                add_edge(edge_id, predicate_name_to_id[relation.predicate], node_ids_in_edge, relation.weight);
                edge_id++;
            }
            db_file.close();

            this->shared_edges.reserve(this->memberships.size());
            this->connected_nodes.reserve(this->memberships.size());
            for(auto const& [edge_id, edge] :this->edges){
                for(size_t i = 0; i < edge.size(); i++){
                    for(size_t j = 0; j < edge.size(); j++){
                        if(edge[i] != edge[j]){
                            this->shared_edges[edge[i]][edge[j]].emplace_back(edge_id);
                            this->connected_nodes[edge[i]].insert(edge[j]);
                        }
                    }
                }
            }
        }else{
            throw FileNotOpenedException(db_file_path);
        }
        if(this->number_of_nodes() < 3){
            throw HyperGraphSizeException();
        }
    }
}


HyperGraph::~HyperGraph() = default;

void HyperGraph::add_edge(EdgeId edge_id, PredicateId predicate, vector<NodeId> node_ids, double weight) {
    this->edges.insert({edge_id, node_ids});
    this->edge_weights[edge_id] = weight;
    this->predicates[edge_id] = predicate;
    for(size_t i = 0; i < node_ids.size(); i++){
        this->memberships[node_ids[i]].emplace_back(edge_id);
    }
}

const unordered_map<NodeId, vector<EdgeId>>& HyperGraph::get_shared_edges(NodeId node){
    auto it = shared_edges.find(node);
    return it->second;
}

bool HyperGraph::has_shared_edges(NodeId node1, NodeId node2)
{
    
    if(this->shared_edges[node1].find(node2) != this->shared_edges[node1].end()){
        return true;
    }
    return false;
}

const unordered_set<NodeId>& HyperGraph::get_connected_nodes(NodeId node){
    return this->connected_nodes.at(node);
}

unordered_map<EdgeId, vector<NodeId>>& HyperGraph::get_edges() {
    return this->edges;
}
const vector<NodeId>& HyperGraph::get_edge(EdgeId edge_id) {
    return this->edges.at(edge_id);
}

unordered_set<NodeId> HyperGraph::get_nodes() {
    return get_keys(this->memberships);
}

unordered_map<EdgeId, PredicateId> HyperGraph::get_predicates(){
    return this->predicates;
}
PredicateId HyperGraph::get_predicate_id(EdgeId edge) {
    return this->predicates[edge];
}

PredicateId HyperGraph::get_predicate_id(Predicate predicate) {
    return this->predicate_name_to_id[predicate];
}

string_view HyperGraph::get_predicate(EdgeId edge_id) {
    return this->predicate_id_to_name[this->predicates[edge_id]];
}


unordered_map<NodeId, vector<EdgeId>> HyperGraph::get_memberships() {
    return this->memberships;
}
const vector<EdgeId>& HyperGraph::get_memberships(NodeId node_id){
    return this->memberships[node_id];
}

size_t HyperGraph::number_of_nodes() {
    return get_nodes().size();
}

size_t HyperGraph::number_of_edges() {
    return edges.size();
}

double HyperGraph::get_edge_weight(EdgeId edge_id) {
    return this->edge_weights[edge_id];
}

void HyperGraph::print() {
    cout << endl;
    cout << "Edges (edge id | node ids)\n";
    for(const auto& edge: this->edges){
        cout << edge.first << " | ";
        for (auto node : edge.second) {
            cout << node << " ";
        }
        cout << "\n";
    }
    cout << endl;
    cout << "Predicates (edge id | predicate)\n";
    for(const auto& predicate: this->predicates){
        cout << predicate.first << " | " << predicate.second << "\n";
    }
    cout << endl;
}


pair<vector<EdgeId>, vector<EdgeId>> HyperGraph::get_unary_binary_edges(NodeId node) {
    set<EdgeId> unaries;
    set<EdgeId> binaries;
    const vector<EdgeId>& all_edges = get_memberships(node);
    for(const auto& edge:all_edges){
        if(is_unary(predicates[edge])){
            unaries.insert(edge);
        }else{
            binaries.insert(edge);
        }
    }
    vector<EdgeId> unary_edges(unaries.begin(), unaries.end());
    vector<EdgeId> binary_edges(binaries.begin(), binaries.end());
    return {unary_edges, binary_edges};
}

const NodeId& HyperGraph::get_next_node(NodeId start, EdgeId edge) {
    const vector<NodeId>& members = get_edge(edge);
    if(start == members[0]){
        return members[1];
    }else{
        return members[0];
    }
}

bool HyperGraph::is_unary(PredicateId predicate) {
    return has(unary_predicates, predicate);
}

string HyperGraph::get_predicate_from_id(PredicateId predicate) {
    return this->predicate_id_to_name[predicate];
}

unordered_map<EdgeId, double> HyperGraph::get_edge_weights() {
    return this->edge_weights;
}

unordered_map<PredicateId, Predicate> HyperGraph::get_predicate_id_to_name() {
    return this->predicate_id_to_name;
}

unordered_map<Predicate, PredicateId> HyperGraph::get_predicate_name_to_id() {
    return this->predicate_name_to_id;
}

bool HyperGraph::is_unary_edge(EdgeId edge) {
    return is_unary(predicates[edge]);
}

bool HyperGraph::has_unaries() const {
    return !this->unary_predicates.empty();
}
