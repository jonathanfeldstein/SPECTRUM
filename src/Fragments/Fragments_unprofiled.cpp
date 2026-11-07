#include "Fragments.h"

#include <utility>

#include "../Utils/profiling_utils.h"


void FragmentCollection::storing_fragments(
    const vector<pair<SmallFragment, SmallLocation>>& fragments,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count)
{

    for (const auto& fragment : fragments) {

        Fragment materialized_fragment(
            fragment.first.begin(),
            fragment.first.end()
        );

        Location materialized_location(
            fragment.second.begin(),
            fragment.second.end()
        );

        local_fragment_count[
            std::move(materialized_fragment)
        ].push_back(
            std::move(materialized_location)
        );
    }

}

void FragmentCollection::storing_cyclic_fragments(
    const vector<pair<SmallFragment, SmallLocation>>& fragments,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count)
{    
    for (const auto& fragment : fragments) {

        Fragment materialized_fragment(
            fragment.first.begin(),
            fragment.first.end()
        );

        Location materialized_location(
            fragment.second.begin(),
            fragment.second.end()
        );

        local_fragment_cyclic_count[
            std::move(materialized_fragment)
        ].push_back(
            std::move(materialized_location)
        );
    }

}

void FragmentCollection::canonicalize_fragments()
{
    const int num_threads =
        omp_get_max_threads();

    //
    // Build vector of cyclic entries
    //
    std::vector<std::pair<Fragment,std::vector<Location>*>> cyclic_entries;

    cyclic_entries.reserve(fragment_cyclic_count.size());

    size_t cyclic_location_count = 0;

    for (auto& entry : fragment_cyclic_count) {
        cyclic_entries.emplace_back(entry.first, &entry.second);
    }

    unordered_map<Fragment, vector<Location>, VectorHash> canonical_cyclic_count;
    canonical_cyclic_count.reserve(this->fragment_cyclic_count.size()*5);
    std::unordered_map<Fragment,vector<Location>, VectorHash> canonical_fragment_count;
    canonical_fragment_count.reserve(this->fragment_count.size()*2);

    vector<unordered_map<Fragment,vector<Location>,VectorHash>> local_cyclic_maps(num_threads);

    for (auto& map : local_cyclic_maps)
    {
        map.reserve(cyclic_entries.size() / num_threads);
    }

    vector<unordered_map<Fragment,vector<Location>,VectorHash>> local_fragment_maps(num_threads);
    for (auto& map : local_fragment_maps)
    {
        map.reserve(this->fragment_cyclic_count.size()*5);
    }

    Timer cyclic_timer("Cyclic Canonicalisation");

#pragma omp parallel
{
    const int tid = omp_get_thread_num();

    auto& local_cyclic_map = local_cyclic_maps[tid];

    auto& local_map = local_fragment_maps[tid];
 
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < cyclic_entries.size(); ++i) {

        const Fragment& fragment = cyclic_entries[i].first;

        auto ordered_indices_possibilities = find_canonical_orders(fragment);


        auto& locations = *cyclic_entries[i].second;

        for (auto& location : locations) {
            Fragment canonical_fragment = fragment;

            Location canonical_location = location;

            std::tie(canonical_fragment, canonical_location) = reorder_cyclic_fragment(fragment,location,ordered_indices_possibilities);

            
            auto sub_patterns = find_bodies_of_cyclic_fragment(fragment, location);

            for (auto& sub_pattern : sub_patterns) {
                local_map[sub_pattern.first].emplace_back(sub_pattern.second);
            }

            local_cyclic_map[canonical_fragment].push_back(std::move(canonical_location));
        }

    }

}

    parallel_merge_maps(
    local_cyclic_maps,
    canonical_cyclic_count,
    num_threads
    );

    parallel_merge_maps(
        local_fragment_maps,
        fragment_count,
        num_threads
    );
    
    cyclic_timer.Stop();

    size_t location_count = 0;

    vector<pair<Fragment,vector<Location>*>> fragment_entries;
    fragment_entries.reserve(fragment_count.size());
    for (auto& entry : fragment_count) {
        fragment_entries.emplace_back(entry.first,&entry.second);
    }

    vector<unordered_map<Fragment,vector<Location>,VectorHash>> local_normal_fragment_maps(num_threads);
    for (auto& map : local_normal_fragment_maps)
    {
        map.reserve(this->fragment_count.size()*2);
    }

    Timer normal_fragment("Normal Fragments Canonicalisation");

#pragma omp parallel
{
        unordered_map<Fragment, vector<Location>, VectorHash> local_map;
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < fragment_entries.size(); ++i)
    {

        const Fragment& fragment = fragment_entries[i].first;

        auto& locations = *fragment_entries[i].second;

        if (locations.size() < 2) {
            continue;
        }

        Fragment reversed_fragment = fragment;

        reverse(reversed_fragment.begin(), reversed_fragment.end());

        if (reversed_fragment < fragment) {
            auto& target = local_map[reversed_fragment];
            target.reserve(target.size() + locations.size());

            for (auto& location : locations)
            {
                reverse(location.begin(),location.end());
            }
            target.insert(
                target.end(),
                std::make_move_iterator(locations.begin()),
                std::make_move_iterator(locations.end())
            );
        } else if (fragment < reversed_fragment) {
            auto& target = local_map[fragment];

            target.insert(
                target.end(),
                std::make_move_iterator(locations.begin()),
                std::make_move_iterator(locations.end())
            );
        }else {
            auto& target = local_map[fragment];

            target.reserve(target.size() + locations.size());

            for (auto& location : locations) {
                Location reversed_location = location;

                std::reverse(
                    reversed_location.begin(),
                    reversed_location.end()
                );

                if (reversed_location < location)
                {
                    target.push_back(
                        std::move(reversed_location)
                    );
                }
                else
                {
                    target.push_back(
                        std::move(location)
                    );
                }
            }
        }
    }

}

    parallel_merge_maps(
    local_normal_fragment_maps,
    canonical_fragment_count,
    num_threads
    );

    normal_fragment.Stop();


    fragment_count = std::move(canonical_fragment_count);

    fragment_cyclic_count = std::move(canonical_cyclic_count);
}
void FragmentCollection::get_fragments(
    size_t number_of_random_walks,
    size_t max_depth,
    bool force_cyclic)
{
    unordered_set<NodeId> nodes_set = hg.get_nodes();
    vector<NodeId> nodes(nodes_set.begin(), nodes_set.end());
    nodes.reserve(nodes_set.size());

    const int num_threads = omp_get_max_threads();

    vector<unordered_map<Fragment,vector<Location>,VectorHash>> local_fragment_counts(num_threads);
    vector<unordered_map<Fragment,vector<Location>,VectorHash>> local_fragment_cyclic_counts(num_threads);

#pragma omp parallel
{
    const int tid = omp_get_thread_num();

    auto& local_fragment_count =
        local_fragment_counts[tid];

    auto& local_fragment_cyclic_count =
        local_fragment_cyclic_counts[tid];

    local_fragment_count.reserve(50000);
    local_fragment_cyclic_count.reserve(50000);

#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < nodes.size(); i++) {



        NodeId node =nodes[i];

        const vector<EdgeId>& unaries = unary_edges.at(node);
        const vector<EdgeId>& binaries = binary_edges.at(node);

        boost::container::small_vector<EdgeId,3> previous_edges;

        vector<pair<SmallFragment, SmallLocation>> initial_unary_fragment;
        initial_unary_fragment.reserve(unaries.size() * unaries.size());

        for (auto unary : unaries) {
            previous_edges.push_back(unary);

            SmallFragment fragment{
                edge_to_predicate[unary]
            };

            SmallLocation unary_location{
                unary
            };
            initial_unary_fragment.emplace_back(fragment,unary_location);
        }

        storing_fragments(
            initial_unary_fragment,
            local_fragment_count
        );

        set_initial_fragments(
            node,
            unaries,
            binaries,
            local_fragment_count
        );


        SmallFragment null_fragment{0};
        SmallLocation null_location{};
        initial_unary_fragment.emplace_back(null_fragment, null_location);

        const auto& connected_nodes = hg.get_connected_nodes(node);


        for (auto binary: binaries) {

            const NodeId& next_node = hg.get_next_node(node,binary);

            vector<pair<SmallFragment, SmallLocation>> final_fragments;
            final_fragments.reserve(initial_unary_fragment.size());

            for (const auto& fragment: initial_unary_fragment) {

                SmallFragment new_fragment;

                new_fragment.insert(
                    new_fragment.end(),
                    fragment.first.begin(),
                    fragment.first.end()
                );
                new_fragment.push_back(hg.get_predicate_id(binary));

                SmallLocation new_location;


                new_location.insert(
                    new_location.end(),
                    fragment.second.begin(),
                    fragment.second.end()
                );
                new_location.push_back(binary);

                final_fragments.emplace_back(std::move(new_fragment),std::move(new_location));

            }

            size_t left_walks = 1 + ((number_of_random_walks - 1)/binaries.size());

            previous_edges.push_back(binary);


            next_step(
                left_walks,
                next_node,
                node,
                1,
                final_fragments,
                previous_edges,
                connected_nodes,
                max_depth,
                force_cyclic,
                local_fragment_count,
                local_fragment_cyclic_count
            );

            previous_edges.pop_back();

        }

    }

}

    parallel_merge_maps(
    local_fragment_counts,
    this->fragment_count,
    num_threads
    );

    parallel_merge_maps(
        local_fragment_cyclic_counts,
        this->fragment_cyclic_count,
        num_threads
    );
    
    canonicalize_fragments();

}

void FragmentCollection::set_initial_fragments(NodeId node, 
                                               const vector<EdgeId> &unaries, 
                                               const vector<EdgeId> &binaries,
                                               unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count 
                                               ) 
{
    size_t unary_pairs = (unaries.size() * (unaries.size() - 1)) / 2;

    vector<pair<SmallFragment, SmallLocation>> unary_fragments;
    unary_fragments.reserve(unary_pairs);

    for(int i=0; i<unaries.size(); i++){
        for(int j = i+1; j<unaries.size(); j++){
            SmallFragment fragment{};
            fragment.push_back(edge_to_predicate[unaries[i]]);
            fragment.push_back(edge_to_predicate[unaries[j]]);
            SmallLocation unary_location{
                unaries[i],
                unaries[j]
            };
            unary_fragments.emplace_back(fragment, unary_location);
        }
    }
    storing_fragments(unary_fragments, local_fragment_count);
    // Find all neighbors of node and save all binary edges that connect neighbor to node
    unordered_map<NodeId, vector<EdgeId>> node_edge_memberships;
    node_edge_memberships.reserve(binaries.size());

    for(auto binary:binaries){
        const NodeId& next = hg.get_next_node(node, binary);
        node_edge_memberships[next].push_back(binary);
    }
    //for each neighbor, if it is connected by more than 1 binary edge, create all fragments of pairs of binary edges connecting to the node
    for(auto members:node_edge_memberships){
        if(members.second.size()>1){
            
            size_t pair_count = (members.second.size() * (members.second.size()-1)) / 2;
            vector<pair<SmallFragment, SmallLocation>> binary_fragments{};
            binary_fragments.reserve(pair_count);

            for(int i=0; i<members.second.size(); i++){
                for(int j = i+1; j<members.second.size(); j++){

                    SmallFragment fragment{};
                    fragment.push_back(edge_to_predicate[members.second[i]]);
                    fragment.push_back(edge_to_predicate[members.second[j]]);

                    SmallLocation binary_location{
                        members.second[i],
                        members.second[j]
                    };
                    binary_fragments.emplace_back(fragment, binary_location);
                }
            }
            storing_fragments(binary_fragments, local_fragment_count);
        }
    }
}
void FragmentCollection::next_step(
    size_t number_of_random_walks,
    NodeId node,
    NodeId starting_node,
    int depth,
    // const vector<pair<Fragment, Location>>& fragments,
    const vector<pair<SmallFragment, SmallLocation>>& fragments,
    boost::container::small_vector<EdgeId,3>& previous_edges,
    const unordered_set<NodeId>& connected_nodes,
    size_t max_depth,
    bool force_cyclic,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count)
{

    if ((depth >= 3) && starting_node == node) {

        storing_cyclic_fragments(
            fragments,
            local_fragment_count,
            local_fragment_cyclic_count
        );

    }

    if (depth < max_depth) {


        vector<pair<SmallFragment, SmallLocation>> new_fragments;
        new_fragments.reserve(fragments.size() * 2);


        const vector<EdgeId>& unaries = unary_edges.at(node);

        vector<EdgeId> new_unaries;
        new_unaries.reserve(unaries.size());
        vector<EdgeId> new_binaries;
        new_binaries.reserve(binary_edges.at(node).size());



        for (auto u : unaries) {
            if (!has(previous_edges,u)) {
                new_unaries.push_back(u);
            }
        }

        vector<EdgeId> cyclic_binaries;
        cyclic_binaries.reserve(new_binaries.size());
        vector<EdgeId> non_cyclic_binaries;
        non_cyclic_binaries.reserve(new_binaries.size());
        vector<EdgeId> final_binaries;
        final_binaries.reserve(number_of_random_walks);

        bool store_fragments = true;

        if (force_cyclic && depth >= 2) {
            

            if(connected_nodes.find(node) != connected_nodes.end()){

                const unordered_map<NodeId, vector<EdgeId>>& shared_edges = hg.get_shared_edges(node);

                for (auto& b : shared_edges.at(starting_node)) {
                    if (!has(previous_edges,b)) {
                        cyclic_binaries.push_back(b);
                    }
                }
            } 


            if (cyclic_binaries.size() >= number_of_random_walks) {

                final_binaries = get_random_elements(
                    cyclic_binaries,
                    number_of_random_walks
                );

            
            } else {

                size_t needed =
                number_of_random_walks - cyclic_binaries.size();

                const auto& edges =
                    binary_edges.at(node);

                size_t forbidden_estimate =
                    previous_edges.size()
                    + cyclic_binaries.size();

                size_t available_estimate =
                    edges.size() > forbidden_estimate
                        ? edges.size() - forbidden_estimate
                        : 0;

                if (needed <= 4 && available_estimate > needed * 3) {
                     auto non_cyclic_final_binaries = direct_random_non_cyclic_selection(
                        edges,
                        previous_edges,
                        cyclic_binaries,
                        needed
                    );
                    final_binaries = std::move(cyclic_binaries);
                    final_binaries.insert(final_binaries.end(), 
                                        std::make_move_iterator(non_cyclic_final_binaries.begin()),
                                        std::make_move_iterator(non_cyclic_final_binaries.end())
                                        );

                } else {
                    vector<EdgeId> previous_edges_cpy(previous_edges.begin(),previous_edges.end());
                    previous_edges_cpy.insert(previous_edges_cpy.end(),cyclic_binaries.begin(),cyclic_binaries.end()); 
                
                    non_cyclic_binaries = edges;
                    non_cyclic_binaries.erase(
                        std::remove_if(
                            non_cyclic_binaries.begin(),
                            non_cyclic_binaries.end(),
                            [&](EdgeId e)
                            {
                                return has(previous_edges_cpy,e);
                            }
                        ),
                        non_cyclic_binaries.end()
                    );

                    final_binaries = std::move(cyclic_binaries);

                    if (!non_cyclic_binaries.empty()) {

                        if(non_cyclic_binaries.size() > (number_of_random_walks - cyclic_binaries.size())) {


                            auto non_cyclic_final_binaries = get_random_elements(non_cyclic_binaries,
                                                                                number_of_random_walks - cyclic_binaries.size());
                            


                            final_binaries.insert(final_binaries.end(), 
                                                std::make_move_iterator(non_cyclic_final_binaries.begin()),
                                                std::make_move_iterator(non_cyclic_final_binaries.end())
                                                );
                        }else {
                            final_binaries.insert(final_binaries.end(), 
                            std::make_move_iterator(non_cyclic_binaries.begin()),
                            std::make_move_iterator(non_cyclic_binaries.end())
                            );
                        }

                        

                    }
                }
            }

        } else {

            const vector<EdgeId>& binaries = binary_edges.at(node);


            new_binaries = binaries;
            new_binaries.erase(
                std::remove_if(
                    new_binaries.begin(),
                    new_binaries.end(),
                    [&](EdgeId e)
                    {
                        return has(previous_edges,e);
                    }
                ),
                new_binaries.end()
            );

            if (new_binaries.size() > number_of_random_walks) {

                final_binaries =
                    get_random_elements(
                        new_binaries,
                        number_of_random_walks
                    );



            }else {

                final_binaries = std::move(new_binaries);

            }
        }


        for (const auto& fragment : fragments) {

            SmallFragment new_fragment;

            new_fragment.insert(
                new_fragment.end(),
                fragment.first.begin(),
                fragment.first.end()
            );
            new_fragment.push_back(0);
            new_fragments.emplace_back(std::move(new_fragment),fragment.second);

        }

        for (auto unary : new_unaries) {

            previous_edges.push_back(unary);



            for (const auto& fragment : fragments) {

                if (!has(fragment.second,unary)) {
                    SmallFragment new_fragment;

                    new_fragment.insert(
                        new_fragment.end(),
                        fragment.first.begin(),
                        fragment.first.end()
                    );
                    new_fragment.push_back(edge_to_predicate[unary]);

                    SmallLocation new_location;

                    new_location.insert(
                        new_location.end(),
                        fragment.second.begin(),
                        fragment.second.end()
                    );
                    new_location.push_back(unary);
                    new_fragments.emplace_back(std::move(new_fragment),std::move(new_location));

                }
            }
        }


        if (store_fragments) {
            storing_fragments(new_fragments, local_fragment_count);
        }


        for (auto binary : final_binaries) {

            const NodeId& next_node = hg.get_next_node(node,binary);

            vector<pair<SmallFragment,SmallLocation>> final_fragments;
            final_fragments.reserve(new_fragments.size());


            for (const auto& fragment : new_fragments) {

                SmallFragment new_fragment;

                new_fragment.insert(
                    new_fragment.end(),
                    fragment.first.begin(),
                    fragment.first.end()
                );
                new_fragment.push_back(edge_to_predicate[binary]);

                SmallLocation new_location;

                new_location.insert(
                    new_location.end(),
                    fragment.second.begin(),
                    fragment.second.end()
                );
                new_location.push_back(binary);

                final_fragments.emplace_back(std::move(new_fragment), std::move(new_location));
                
            }


            previous_edges.push_back(binary);


            next_step(
                1 +((number_of_random_walks - 1)/final_binaries.size()),
                next_node,
                starting_node,
                depth + 1,
                final_fragments,
                previous_edges,
                connected_nodes,
                max_depth,
                force_cyclic,
                local_fragment_count,
                local_fragment_cyclic_count
            );

            previous_edges.pop_back();
        }
    }
}
FragmentCollection::FragmentCollection(const HyperGraph& hg, size_t number_of_random_walks, size_t max_depth, bool force_cyclic){
    this->hg = hg;
    this->edge_to_predicate = this->hg.get_predicates();
    for(auto node:this->hg.get_nodes()){
        vector<EdgeId> unaries;
        vector<EdgeId> binaries;
        tie(unaries, binaries) = this->hg.get_unary_binary_edges(node);
        this->unary_edges[node] = unaries;
        this->binary_edges[node] = binaries;
        this->unary_binary_edges[node] = make_pair(unaries, binaries);
    }
    get_fragments(number_of_random_walks, max_depth, force_cyclic);
}

unordered_map<Fragment, vector<Location>, VectorHash>&
FragmentCollection::get_fragment_count()
{
    vector<vector<Location>*> items;
    items.reserve(fragment_count.size());

    size_t location_count = 0;

    for (auto& entry : fragment_count) {
        items.push_back(&entry.second);
        location_count += entry.second.size();
    }

    std::cout << "Before deduplication location count:" << location_count << "\n";


#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < items.size(); ++i)
    {
        auto& vec = *items[i];

        if (vec.size() <= 1) {
            continue;
        }

        std::sort(vec.begin(), vec.end());

        vec.erase(
            std::unique(vec.begin(), vec.end()),
            vec.end()
        );
    }

    size_t unique_location_count = 0;

    for (auto& entry : fragment_count) {
        unique_location_count += entry.second.size();
    }

    std::cout << "After cyclic location count:" << unique_location_count << "\n";

    return fragment_count;
}

bool FragmentCollection::is_unary(PredicateId predicate) {
    return this->hg.is_unary(predicate);
}

unordered_map<Fragment, vector<Location>, VectorHash>&
FragmentCollection::get_fragment_cyclic_count()
{
    vector<vector<Location>*> items;
    items.reserve(fragment_cyclic_count.size());

    size_t location_count = 0;

    for (auto& entry : fragment_cyclic_count) {
        items.push_back(&entry.second);
        location_count += entry.second.size();
    }

    std::cout << "Before deduplication cyclic location count:" << location_count << "\n";


#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < items.size(); ++i)
    {
        auto& vec = *items[i];

        if (vec.size() <= 1) {
            continue;
        }else if (vec.size() >2){
            std::sort(vec.begin(), vec.end());
        }

        vec.erase(
            std::unique(vec.begin(), vec.end()),
            vec.end()
        );
    }

    size_t unique_location_count = 0;

    for (auto& entry : fragment_cyclic_count) {
        unique_location_count += entry.second.size();
    }

    std::cout << "After deduplication cyclic location count:" << unique_location_count << "\n";

    return fragment_cyclic_count;
}

string FragmentCollection::get_predicate_from_id(PredicateId predicate) {
    return this->hg.get_predicate_from_id(predicate);
}

const vector<NodeId>& FragmentCollection::get_edge(EdgeId edge) {
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

FragmentCollection::FragmentCollection() =default;
