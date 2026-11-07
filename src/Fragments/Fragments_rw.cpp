/*
 * Fragments_rw.cpp — Random Walks variant of FragmentCollection
 *
 * Key differences from Fragments.cpp:
 *   1. get_fragments(): Performs N independent random walks per node
 *      (each picking ONE random binary edge) instead of fanning out
 *      across all edges with budget splitting.
 *   2. do_random_walk(): Replaces next_step(). No previous_edges
 *      filtering (backtracking allowed). No force_cyclic bias.
 *      Cyclic fragments are still detected naturally when a walk
 *      returns to its starting node at depth 3-4.
 *   3. All other methods (storing_fragments, storing_cyclic_fragments,
 *      set_initial_fragments, etc.) are identical to the original.
 */

#include "Fragments_rw.h"

#include <utility>

// ============================================================
// UNCHANGED METHODS (identical to Fragments.cpp)
// ============================================================

void FragmentCollection::storing_fragments(const map<Fragment, Location>& fragments){
    unordered_map<Fragment, vector<Location>, VectorHash> local_fragment_count;

    for (const auto &fragment: fragments) {
        Fragment reversed_fragment = fragment.first;
        std::reverse(reversed_fragment.begin(), reversed_fragment.end());
        if (fragment.first < reversed_fragment) {
            local_fragment_count[fragment.first].push_back(fragment.second);
        } else if(fragment.first > reversed_fragment){
            vector<NodeId> reversed_location = fragment.second;
            std::reverse(reversed_location.begin(), reversed_location.end());
            local_fragment_count[reversed_fragment].push_back(reversed_location);
        }else{
            vector<NodeId> reversed_location = fragment.second;
            std::reverse(reversed_location.begin(), reversed_location.end());
            if(reversed_location<fragment.second){
                local_fragment_count[reversed_fragment].push_back(reversed_location);
            }else{
                local_fragment_count[fragment.first].push_back(fragment.second);
            }
        }
    }

#pragma omp critical
    {
        for (const auto &entry : local_fragment_count) {
            this->fragment_count[entry.first].insert(this->fragment_count[entry.first].end(),entry.second.begin(), entry.second.end());
        }
    }
}

void FragmentCollection::storing_cyclic_fragments(const map<Fragment, Location>& fragments){
    unordered_map<Fragment, vector<Location>, VectorHash> local_fragment_cyclic_count;

    for (const auto &fragment: fragments) {
        Fragment new_fragment = fragment.first;
        Fragment new_location = fragment.second;
        tie(new_fragment, new_location) = unary_utils::order_cyclic_fragment_and_location(new_fragment, new_location);
        vector<pair<Fragment, Location>> sub_patterns = unary_utils::find_bodies_of_cyclic_fragment(new_fragment, new_location);
        for(auto sub_pattern:sub_patterns){
            map<Fragment, Location> local_sub_pattern{};
            local_sub_pattern.insert({sub_pattern.first, sub_pattern.second});
            storing_fragments(local_sub_pattern);
        }
        local_fragment_cyclic_count[new_fragment].push_back(new_location);
    }

#pragma omp critical
    {
        for (const auto &entry : local_fragment_cyclic_count) {
            this->fragment_cyclic_count[entry.first].insert(this->fragment_cyclic_count[entry.first].end(),entry.second.begin(), entry.second.end());
        }
    }
}

void FragmentCollection::set_initial_fragments(NodeId node, vector<EdgeId> unaries, const vector<EdgeId> &binaries) {
    map<Fragment, Location> unary_fragments{};
    for(int i=0; i<unaries.size(); i++){
        for(int j = i+1; j<unaries.size(); j++){
            Fragment fragment{};
            fragment.push_back(edge_to_predicate[unaries[i]]);
            fragment.push_back(edge_to_predicate[unaries[j]]);
            Location unary_location{unaries[i], unaries[j]};
            unary_fragments.insert({fragment, unary_location});
        }
    }
    storing_fragments(unary_fragments);
    map<NodeId,vector<EdgeId>> node_edge_memberships{};
    for(auto binary:binaries){
        NodeId next = hg.get_next_node(node, binary);
        node_edge_memberships[next].push_back(binary);
    }
    for(auto members:node_edge_memberships){
        if(members.second.size()>1){
            map<Fragment, Location> binary_fragments{};
            for(int i=0; i<members.second.size(); i++){
                for(int j = i+1; j<members.second.size(); j++){
                    Fragment fragment{};
                    fragment.push_back(edge_to_predicate[members.second[i]]);
                    fragment.push_back(edge_to_predicate[members.second[j]]);
                    Location binary_location{members.second[i], members.second[j]};
                    binary_fragments.insert({fragment, binary_location});
                }
            }
            storing_fragments(binary_fragments);
        }
    }
}

FragmentCollection::FragmentCollection(const HyperGraph& hg, size_t number_of_random_walks, size_t max_depth, bool force_cyclic){
    this->hg = hg;
    this->edge_to_predicate = this->hg.get_predicates();
    for(auto node:this->hg.get_nodes()){
        vector<EdgeId> unaries;
        vector<EdgeId> binaries;
        unordered_set<EdgeId> previous_edges{};
        tie(unaries, binaries) = this->hg.get_unary_binary_edges(node);
        this->unary_binary_edges[node] = make_pair(unaries, binaries);
    }
    get_fragments(number_of_random_walks, max_depth, force_cyclic);
}

unordered_map<Fragment, vector<Location>, VectorHash>& FragmentCollection::get_fragment_count() {
    vector<pair<const Fragment, vector<Location>>*> items;
    for(auto& entry : this->fragment_count) {
        items.push_back(&entry);
    }
#pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < items.size(); ++i) {
        set<Location> unique_locations{items[i]->second.begin(), items[i]->second.end()};
        items[i]->second = vector<Location>(unique_locations.begin(), unique_locations.end());
    }
    return this->fragment_count;
}

bool FragmentCollection::is_unary(PredicateId predicate) {
    return this->hg.is_unary(predicate);
}

unordered_map<Fragment, vector<Location>, VectorHash>& FragmentCollection::get_fragment_cyclic_count() {
    vector<pair<const Fragment, vector<Location>>*> items;
    for(auto& entry : this->fragment_cyclic_count) {
        items.push_back(&entry);
    }
#pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < items.size(); ++i) {
        set<Location> unique_locations{items[i]->second.begin(), items[i]->second.end()};
        items[i]->second = vector<Location>(unique_locations.begin(), unique_locations.end());
    }
    return this->fragment_cyclic_count;
}

string FragmentCollection::get_predicate_from_id(PredicateId predicate) {
    return this->hg.get_predicate_from_id(predicate);
}

vector<NodeId> FragmentCollection::get_edge(EdgeId edge) {
    return this->hg.get_edge(edge);
}

PredicateId FragmentCollection::get_predicate_id(Predicate predicate) {
    return this->hg.get_predicate_id(std::move(predicate));
}

bool FragmentCollection::is_unary_edge(EdgeId edge) {
    return is_unary(this->hg.get_predicate_id(edge));
}

pair<vector<EdgeId>, vector<EdgeId>> FragmentCollection::get_unary_binary_edges(NodeId node) {
    return this->unary_binary_edges[node];
}

void FragmentCollection::set_fragment_count(unordered_map<Fragment, vector<Location>, VectorHash> new_fragment_count) {
    this->fragment_count = new_fragment_count;
}

void FragmentCollection::set_fragment_cyclic_count(unordered_map<Fragment, vector<Location>, VectorHash> new_fragment_cyclic_count) {
    this->fragment_cyclic_count = new_fragment_cyclic_count;
}

FragmentCollection::FragmentCollection() =default;


// ============================================================
// MODIFIED METHODS — True Random Walks (with Early Termination)
// ============================================================

void FragmentCollection::get_fragments(size_t number_of_random_walks, size_t max_depth, bool force_cyclic){
    // force_cyclic is ignored in the RW variant (always treated as false)
    unordered_set<NodeId> nodes_set = hg.get_nodes();
    vector<NodeId> nodes(nodes_set.begin(), nodes_set.end());
    nodes.reserve(nodes_set.size());

#pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < nodes.size(); i++) {
        NodeId node = nodes[i];
        vector<EdgeId> unaries = unary_binary_edges[node].first;
        vector<EdgeId> binaries = unary_binary_edges[node].second;

        // Single Unary predicates (same as original)
        map<Fragment, Location> initial_unary_fragment{};
        for(auto unary : unaries){
            Fragment fragment{edge_to_predicate[unary]};
            Location unary_location{unary};
            initial_unary_fragment.insert({fragment, unary_location});
        }
        storing_fragments(initial_unary_fragment);
        set_initial_fragments(node, unaries, binaries);

        if(binaries.empty()) continue;

        // Add the null/empty fragment prefix
        Fragment null_fragment{0};
        Location null_location{};
        initial_unary_fragment.insert({null_fragment, null_location});

        // Thread-local RNG seeded uniquely per node
        std::mt19937 gen(std::random_device{}() ^ (i * 2654435761u));
        std::uniform_int_distribution<size_t> dist(0, binaries.size() - 1);

        // Perform number_of_random_walks INDEPENDENT walks from this node
        for(size_t walk = 0; walk < (size_t)number_of_random_walks; walk++){
            // Pick ONE random binary edge to start the walk
            EdgeId chosen_binary = binaries[dist(gen)];
            NodeId next_node = hg.get_next_node(node, chosen_binary);

            // Extend initial fragments with the chosen binary edge
            map<Fragment, Location> walk_fragments{};
            for(const auto& fragment : initial_unary_fragment){
                Fragment new_fragment = fragment.first;
                new_fragment.push_back(edge_to_predicate[chosen_binary]);
                Location new_location = fragment.second;
                new_location.push_back(chosen_binary);
                walk_fragments.insert({new_fragment, new_location});
            }

            // Initialize visited edges
            std::unordered_set<EdgeId> visited_edges;
            for(auto unary: unaries) {
                visited_edges.insert(unary);
            }
            visited_edges.insert(chosen_binary);

            // Continue the random walk from the next node
            do_random_walk(gen,
                          next_node,
                          node,
                          1,
                          walk_fragments,
                          max_depth,
                          visited_edges);
        }
    }
}

void FragmentCollection::do_random_walk(std::mt19937& gen,
                                        NodeId node,
                                        NodeId starting_node,
                                        int depth,
                                        const map<vector<PredicateId>, Location> &fragments,
                                        size_t max_depth,
                                        std::unordered_set<EdgeId> visited_edges){

    vector<EdgeId> node_unaries = unary_binary_edges[node].first;
    vector<EdgeId> node_binaries = unary_binary_edges[node].second;

    // Check for early termination if any unary is revisited
    bool corrupted = false;
    for(auto unary : node_unaries) {
        if (visited_edges.find(unary) != visited_edges.end()) {
            corrupted = true;
            break;
        }
    }

    if (corrupted) {
        return; // Early termination due to visited edge
    }

    // Add unaries of current node to visited set
    for(auto unary : node_unaries) {
        visited_edges.insert(unary);
    }

    // Natural cyclic fragment detection (no forced bias)
    if((depth == 3 || depth == 4) && starting_node == node){
        map<vector<PredicateId>,Location> cyclic_fragments{};
        for(const auto& fragment:fragments){
            Fragment new_fragment = fragment.first;
            Location new_location = fragment.second;
            cyclic_fragments.insert({new_fragment, new_location});
        }
        storing_cyclic_fragments(cyclic_fragments);
    }

    if(depth <= 3 && depth < (int)max_depth){
        map<vector<PredicateId>, Location> new_fragments{};

        // Extend fragments with null separator (no-unary option)
        for(const auto& fragment:fragments){
            vector<PredicateId> new_fragment = fragment.first;
            new_fragment.push_back(0);
            new_fragments.insert({new_fragment,fragment.second});
        }
        // Extend with unary predicates at the current node
        for(auto unary:node_unaries){
            for(const auto& fragment:fragments){
                vector<PredicateId> new_fragment = fragment.first;
                new_fragment.push_back(edge_to_predicate[unary]);
                Location new_location = fragment.second;
                // Prevent duplicate edge IDs within a single fragment's location
                if(!has(fragment.second,unary)){
                    new_location.push_back(unary);
                    new_fragments.insert({new_fragment,new_location});
                }
            }
        }

        // Store intermediate fragments
        storing_fragments(new_fragments);

        // Pick ONE random binary edge for the next step
        if(!node_binaries.empty()){
            std::uniform_int_distribution<size_t> dist(0, node_binaries.size() - 1);
            EdgeId chosen_binary = node_binaries[dist(gen)];
            NodeId next_node = hg.get_next_node(node, chosen_binary);

            // Corrupted binary edge check
            if (visited_edges.find(chosen_binary) != visited_edges.end()) {
                return; // corrupted, early termination
            }

            map<vector<PredicateId>,Location> extended_fragments{};
            for(const auto& fragment:new_fragments){
                vector<PredicateId> new_fragment = fragment.first;
                new_fragment.push_back(edge_to_predicate[chosen_binary]);
                Location new_location = fragment.second;
                new_location.push_back(chosen_binary);
                extended_fragments.insert({new_fragment,new_location});
            }

            visited_edges.insert(chosen_binary);

            do_random_walk(gen,
                          next_node,
                          starting_node,
                          depth+1,
                          extended_fragments,
                          max_depth,
                          visited_edges);
        }
    }
}
