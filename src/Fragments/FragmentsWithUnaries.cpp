#include "Fragments.h"

void FragmentCollection::get_fragments_with_unaries(
    size_t number_of_random_walks,
    size_t max_depth,
    bool force_cyclic)
{
    unordered_set<NodeId> nodes_set = hg.get_nodes();
    vector<NodeId> nodes(nodes_set.begin(), nodes_set.end());
    nodes.reserve(nodes_set.size());

    const int num_threads = omp_get_max_threads();

    vector<FragmentProfiling> thread_profiles(num_threads);
    vector<RecursiveProfiling> recursive_profiles(num_threads);
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

        auto& profile = thread_profiles[tid];
        auto& recursive_profile = recursive_profiles[tid];

        START_TIMER(setup_start)

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
            local_fragment_count,
            recursive_profile
        );

        set_initial_fragments_with_unaries(
            node,
            unaries,
            binaries,
            local_fragment_count,
            recursive_profile
        );

        SmallFragment null_fragment{0};
        SmallLocation null_location{};
        initial_unary_fragment.emplace_back(null_fragment, null_location);

        STOP_TIMER(
            profile.initial_time,
            setup_start
        )

        DEBUG_ONLY(profile.nodes_processed++;)

        const auto& connected_nodes = hg.get_connected_nodes(node);

        START_TIMER(recursion_start)

        for (auto binary: binaries) {

            START_TIMER(construction_start)

            const NodeId& next_node = hg.get_next_node(node,binary);

            for (auto& fragment: initial_unary_fragment) {

                fragment.first.push_back(hg.get_predicate_id(binary));

                fragment.second.push_back(binary);

            }

            size_t left_walks = 1 + ((number_of_random_walks - 1)/binaries.size());
            

            STOP_TIMER(
                recursive_profile.fragment_construction_time,
                construction_start
            )

            previous_edges.push_back(binary);

            next_step_with_unaries(
                left_walks,
                next_node,
                node,
                1,
                initial_unary_fragment,
                previous_edges,
                connected_nodes,
                max_depth,
                force_cyclic,
                local_fragment_count,
                local_fragment_cyclic_count,
                recursive_profile
            );

            for (auto& fragment: initial_unary_fragment) {
                fragment.first.pop_back();
                fragment.second.pop_back();
            }

            previous_edges.pop_back();

            DEBUG_ONLY(profile.recursive_calls++;)
        }

        STOP_TIMER(
            profile.recursive_time,
            recursion_start
        )
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


    double canonical_start =
        omp_get_wtime();
    
    canonicalize_fragments();

    double canonical_time =
        omp_get_wtime() -
        canonical_start;

    DEBUG_ONLY(print_profiling(thread_profiles, recursive_profiles, canonical_time);)

}

void FragmentCollection::set_initial_fragments_with_unaries(NodeId node, 
                                               const vector<EdgeId> &unaries, 
                                               const vector<EdgeId> &binaries,
                                               unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count, 
                                               RecursiveProfiling& recursive_profile) 
{
    size_t unary_pairs = (unaries.size() * (unaries.size() - 1)) / 2;

    vector<pair<SmallFragment, SmallLocation>> unary_fragments;
    unary_fragments.reserve(unary_pairs);

    for(int i=0; i<unaries.size(); i++){
        for(int j = i+1; j<unaries.size(); j++){
            SmallFragment fragment{};
            fragment.push_back(edge_to_predicate[unaries[i]]);
            fragment.push_back(edge_to_predicate[unaries[j]]);
            SmallLocation unary_location{unaries[i], unaries[j]};
            unary_fragments.emplace_back(fragment, unary_location);
        }
    }
    storing_fragments(unary_fragments, local_fragment_count, recursive_profile);
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
            storing_fragments(binary_fragments, local_fragment_count, recursive_profile);
        }
    }
}

void FragmentCollection::next_step_with_unaries(
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
    RecursiveProfiling& profile)
{

    if ((depth >= 3) && starting_node == node) {

        storing_cyclic_fragments(
            fragments,
            local_fragment_count,
            local_fragment_cyclic_count,
            profile
        );

    }

    if (depth < max_depth) {

        START_TIMER(graph_start)

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

        STOP_TIMER(
            profile.graph_lookup_time,
            graph_start
        )

        vector<EdgeId> cyclic_binaries;
        cyclic_binaries.reserve(new_binaries.size());

        vector<EdgeId> non_cyclic_binaries;
        non_cyclic_binaries.reserve(new_binaries.size());

        vector<EdgeId> final_binaries;
        final_binaries.reserve(number_of_random_walks);

        bool store_fragments = true;
 
        if (force_cyclic && depth >= 2) {

            auto cyclic_binaries =
                find_cyclic_edges(
                    node,
                    starting_node,
                    connected_nodes,
                    previous_edges,
                    profile
                );

            if (cyclic_binaries.size() >= number_of_random_walks) {
                
                START_TIMER(random_time)

                final_binaries =
                    get_random_elements(
                        cyclic_binaries,
                        number_of_random_walks
                    );

                STOP_TIMER(profile.random_selection_time, random_time)

            } else {

                auto non_cyclic_binaries =
                    get_edges_for_next_step(
                        node,
                        previous_edges,
                        cyclic_binaries,
                        number_of_random_walks
                            - cyclic_binaries.size(),
                        profile
                    );

                final_binaries =
                    std::move(cyclic_binaries);

                final_binaries.insert(
                    final_binaries.end(),
                    std::make_move_iterator(
                        non_cyclic_binaries.begin()),
                    std::make_move_iterator(
                        non_cyclic_binaries.end())
                );
            }

        } else {

            static const vector<EdgeId> empty;

            final_binaries =
                get_edges_for_next_step(
                    node,
                    previous_edges,
                    empty,
                    number_of_random_walks,
                    profile
                );

        }

        START_TIMER(construction_start)

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

        STOP_TIMER(profile.fragment_construction_time, construction_start)

        if (store_fragments) {
            storing_fragments(new_fragments, local_fragment_count, profile);
        }


        for (auto binary : final_binaries) {

            START_TIMER(construction_start_2)

            const NodeId& next_node = hg.get_next_node(node,binary);

            for (auto& fragment: new_fragments) {

                fragment.first.push_back(hg.get_predicate_id(binary));

                fragment.second.push_back(binary);

            }
            STOP_TIMER(profile.fragment_construction_time, construction_start_2)

            previous_edges.push_back(binary);


            next_step_with_unaries(
                1 +((number_of_random_walks - 1)/final_binaries.size()),
                next_node,
                starting_node,
                depth + 1,
                new_fragments,
                previous_edges,
                connected_nodes,
                max_depth,
                force_cyclic,
                local_fragment_count,
                local_fragment_cyclic_count,
                profile
            );

            for (auto& fragment: new_fragments) {

                fragment.first.pop_back();

                fragment.second.pop_back();

            }

            previous_edges.pop_back();
        }
    }
}

