#include "Fragments.h"

#include <utility>

// ============================================================
// Constructor
// ============================================================

FragmentCollection::FragmentCollection(
    const HyperGraph& hg,
    size_t number_of_random_walks,
    size_t max_depth,
    bool force_cyclic,
    bool has_unaries)
    : has_unaries_(has_unaries)
{
    this->hg = hg;
    this->edge_to_predicate = this->hg.get_predicates();
    for(auto node : this->hg.get_nodes()){
        vector<EdgeId> unaries;
        vector<EdgeId> binaries;
        tie(unaries, binaries) = this->hg.get_unary_binary_edges(node);
        this->unary_edges[node] = unaries;
        this->binary_edges[node] = binaries;
        this->unary_binary_edges[node] = make_pair(unaries, binaries);
    }
    get_fragments(number_of_random_walks, max_depth, force_cyclic);
}

FragmentCollection::FragmentCollection() = default;

// ============================================================
// Fragment Storage
// ============================================================

void FragmentCollection::storing_fragments(
    const vector<pair<SmallFragment, SmallLocation>>& fragments,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
    RecursiveProfiling& profile)
{
    START_TIMER(start)

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

    STOP_TIMER(
        profile.storing_fragments_time,
        start
    )

    DEBUG_ONLY(profile.fragments_stored +=
        fragments.size();)
}

void FragmentCollection::storing_fragments(
    const pair<SmallFragment, SmallLocation>& fragment,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
    RecursiveProfiling& profile)
{
    START_TIMER(start)
    Fragment materialized_fragment(fragment.first.begin(), fragment.first.end());
    Location materialized_location(fragment.second.begin(), fragment.second.end());
    local_fragment_count[std::move(materialized_fragment)].push_back(std::move(materialized_location));
    STOP_TIMER(profile.storing_fragments_time, start)
    DEBUG_ONLY(profile.fragments_stored += 1;)
}

void FragmentCollection::storing_cyclic_fragments(
    const vector<pair<SmallFragment, SmallLocation>>& fragments,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count,
    RecursiveProfiling& profile)
{
    START_TIMER(start)

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

    STOP_TIMER(
        profile.storing_cyclic_time,
        start
    )

    DEBUG_ONLY(profile.cyclic_fragments_stored +=
        fragments.size();)
}

void FragmentCollection::storing_cyclic_fragments(
    const pair<SmallFragment, SmallLocation>& fragment,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count,
    RecursiveProfiling& profile)
{
    START_TIMER(start)
    Fragment materialized_fragment(fragment.first.begin(), fragment.first.end());
    Location materialized_location(fragment.second.begin(), fragment.second.end());
    local_fragment_cyclic_count[std::move(materialized_fragment)].push_back(std::move(materialized_location));
    STOP_TIMER(profile.storing_cyclic_time, start)
    DEBUG_ONLY(profile.cyclic_fragments_stored += 1;)
}

// ============================================================
// Canonicalization
// ============================================================

void FragmentCollection::canonicalize_fragments()
{
    const int num_threads =
        omp_get_max_threads();

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

        auto ordered_indices_possibilities = has_unaries_ ? unary_utils::find_canonical_orders(fragment)
                                                           : plain_utils::find_canonical_orders(fragment);


        auto& locations = *cyclic_entries[i].second;

        for (auto& location : locations) {
            Fragment canonical_fragment = fragment;

            Location canonical_location = location;

            std::tie(canonical_fragment, canonical_location) = has_unaries_
                ? unary_utils::reorder_cyclic_fragment(fragment,location,ordered_indices_possibilities)
                : plain_utils::reorder_cyclic_fragment(fragment,location,ordered_indices_possibilities);

            
            auto sub_patterns = has_unaries_ ? unary_utils::find_bodies_of_cyclic_fragment(fragment, location)
                                              : plain_utils::find_bodies_of_cyclic_fragment(fragment, location);

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
    const int tid = omp_get_thread_num();

    auto& local_map = local_normal_fragment_maps[tid];

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

    fragment_count.swap(canonical_fragment_count);

    fragment_cyclic_count.swap(canonical_cyclic_count);

}

// ============================================================
// Initial Fragment Setup (merged: handles both unary and non-unary)
// ============================================================

void FragmentCollection::set_initial_fragments(
    NodeId node,
    const vector<EdgeId>& binaries,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_count,
    unordered_map<Fragment, vector<Location>, VectorHash>& local_fragment_cyclic_count,
    RecursiveProfiling& recursive_profile)
{
    // Unary pair fragments (only when has_unaries_)
    if (has_unaries_) {
        const vector<EdgeId>& unaries = unary_edges.at(node);
        size_t unary_pairs = (unaries.size() * (unaries.size() - 1)) / 2;
        vector<pair<SmallFragment, SmallLocation>> unary_fragments;
        unary_fragments.reserve(unary_pairs);
        for (size_t i = 0; i < unaries.size(); i++) {
            for (size_t j = i + 1; j < unaries.size(); j++) {
                SmallFragment fragment{};
                // canonical order (by predicate id), so each pair of predicates has one pattern
                EdgeId first = unaries[i], second = unaries[j];
                if (edge_to_predicate[second] < edge_to_predicate[first]) std::swap(first, second);
                fragment.push_back(edge_to_predicate[first]);
                fragment.push_back(edge_to_predicate[second]);
                SmallLocation unary_location{first, second};
                unary_fragments.emplace_back(fragment, unary_location);
            }
        }
        storing_fragments(unary_fragments, local_fragment_count, recursive_profile);
    }

    // Binary pair fragments (shared between both paths)
    // Find all neighbors of node and save all binary edges that connect neighbor to node
    unordered_map<NodeId, vector<EdgeId>> node_edge_memberships;
    node_edge_memberships.reserve(binaries.size());

    for (auto binary : binaries) {
        const NodeId& next = hg.get_next_node(node, binary);
        node_edge_memberships[next].push_back(binary);
    }
    //for each neighbor, if it is connected by more than 1 binary edge, create all fragments of pairs of binary edges connecting to the node
    for (auto& members : node_edge_memberships) {
        if (members.second.size() > 1) {
            size_t pair_count = (members.second.size() * (members.second.size() - 1)) / 2;
            vector<pair<SmallFragment, SmallLocation>> binary_fragments{};
            binary_fragments.reserve(pair_count);

            for (size_t i = 0; i < members.second.size(); i++) {
                for (size_t j = i + 1; j < members.second.size(); j++) {
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
            storing_cyclic_fragments(binary_fragments, local_fragment_count, local_fragment_cyclic_count, recursive_profile);
        }
    }
}

// ============================================================
// Fragment Discovery (merged: dispatches based on has_unaries_)
// ============================================================

void FragmentCollection::get_fragments(
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

        NodeId node = nodes[i];

        const vector<EdgeId>& binaries = binary_edges.at(node);

        boost::container::small_vector<EdgeId,3> previous_edges;

        set_initial_fragments(
            node,
            binaries,
            local_fragment_count,
            local_fragment_cyclic_count,
            recursive_profile
        );

        if (has_unaries_) {
            // ---- With-unaries path ----
            // Build initial unary fragment vector (one entry per unary edge)
            const vector<EdgeId>& unaries = unary_edges.at(node);

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
                initial_unary_fragment.emplace_back(fragment, unary_location);
            }

            storing_fragments(
                initial_unary_fragment,
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

            for (auto binary : binaries) {

                START_TIMER(construction_start)

                const NodeId& next_node = hg.get_next_node(node, binary);

                for (auto& fragment : initial_unary_fragment) {
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

                for (auto& fragment : initial_unary_fragment) {
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
        } else {
            // ---- No-unaries path ----
            pair<SmallFragment, SmallLocation> final_fragment;

            STOP_TIMER(
                profile.initial_time,
                setup_start
            )

            DEBUG_ONLY(profile.nodes_processed++;)

            const auto& connected_nodes = hg.get_connected_nodes(node);

            START_TIMER(recursion_start)

            for (auto binary : binaries) {

                START_TIMER(construction_start)

                const NodeId& next_node = hg.get_next_node(node, binary);

                final_fragment.first.push_back(hg.get_predicate_id(binary));
                final_fragment.second.push_back(binary);

                size_t left_walks = 1 + ((number_of_random_walks - 1)/binaries.size());
                

                STOP_TIMER(
                    recursive_profile.fragment_construction_time,
                    construction_start
                )

                previous_edges.push_back(binary);

                next_step_no_unaries(
                    left_walks,
                    next_node,
                    node,
                    1,
                    final_fragment,
                    previous_edges,
                    connected_nodes,
                    max_depth,
                    force_cyclic,
                    local_fragment_count,
                    local_fragment_cyclic_count,
                    recursive_profile
                );

                final_fragment.first.pop_back();
                final_fragment.second.pop_back();

                previous_edges.pop_back();

                DEBUG_ONLY(profile.recursive_calls++;)
            }

            STOP_TIMER(
                profile.recursive_time,
                recursion_start
            )
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


    double canonical_start =
        omp_get_wtime();
    
    canonicalize_fragments();

    double canonical_time =
        omp_get_wtime() -
        canonical_start;

    DEBUG_ONLY(print_profiling(thread_profiles, recursive_profiles, canonical_time);)

}

// ============================================================
// Edge Selection Helper (shared between next_step variants)
// ============================================================

vector<EdgeId> FragmentCollection::select_binaries_for_step(
    size_t number_of_random_walks,
    NodeId node,
    NodeId starting_node,
    int depth,
    boost::container::small_vector<EdgeId,3>& previous_edges,
    const unordered_set<NodeId>& connected_nodes,
    bool force_cyclic,
    RecursiveProfiling& profile)
{
    vector<EdgeId> final_binaries;
    final_binaries.reserve(number_of_random_walks);

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

    return final_binaries;
}

// ============================================================
// Recursive Walk: No Unaries (single fragment)
// ============================================================

void FragmentCollection::next_step_no_unaries(
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
    RecursiveProfiling& profile)
{

    if ((depth >= 3) && starting_node == node) {

        storing_cyclic_fragments(
            fragment,
            local_fragment_count,
            local_fragment_cyclic_count,
            profile
        );

    }

    if (depth < max_depth) {

        vector<EdgeId> final_binaries = select_binaries_for_step(
            number_of_random_walks, node, starting_node, depth,
            previous_edges, connected_nodes, force_cyclic, profile
        );

        storing_fragments(fragment, local_fragment_count, profile);

        for (auto binary : final_binaries) {

            START_TIMER(construction_start_2)

            const NodeId& next_node = hg.get_next_node(node, binary);

            fragment.first.push_back(hg.get_predicate_id(binary));
            fragment.second.push_back(binary);

            STOP_TIMER(profile.fragment_construction_time, construction_start_2)

            previous_edges.push_back(binary);

            next_step_no_unaries(
                1 + ((number_of_random_walks - 1)/final_binaries.size()),
                next_node,
                starting_node,
                depth + 1,
                fragment,
                previous_edges,
                connected_nodes,
                max_depth,
                force_cyclic,
                local_fragment_count,
                local_fragment_cyclic_count,
                profile
            );

            fragment.first.pop_back();
            fragment.second.pop_back();

            previous_edges.pop_back();
        }
    }
}

// ============================================================
// Recursive Walk: With Unaries (vector of fragments)
// ============================================================

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

        for (auto u : unaries) {
            if (!has(previous_edges, u)) {
                new_unaries.push_back(u);
            }
        }

        STOP_TIMER(
            profile.graph_lookup_time,
            graph_start
        )

        vector<EdgeId> final_binaries = select_binaries_for_step(
            number_of_random_walks, node, starting_node, depth,
            previous_edges, connected_nodes, force_cyclic, profile
        );

        START_TIMER(construction_start)

        // Extend each existing fragment with a null (no-unary) separator
        for (const auto& fragment : fragments) {

            SmallFragment new_fragment;

            new_fragment.insert(
                new_fragment.end(),
                fragment.first.begin(),
                fragment.first.end()
            );
            new_fragment.push_back(0);
            new_fragments.emplace_back(std::move(new_fragment), fragment.second);

        }

        // Extend each existing fragment with each new unary edge
        const size_t previous_edges_mark = previous_edges.size();
        for (auto unary : new_unaries) {

            previous_edges.push_back(unary);

            for (const auto& fragment : fragments) {

                if (!has(fragment.second, unary)) {
                    
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
                    new_fragments.emplace_back(std::move(new_fragment), std::move(new_location));

                }
            }
        }

        STOP_TIMER(profile.fragment_construction_time, construction_start)
        previous_edges.resize(previous_edges_mark); // unaries of this node are only excluded on this path

        storing_fragments(new_fragments, local_fragment_count, profile);

        for (auto binary : final_binaries) {

            START_TIMER(construction_start_2)

            const NodeId& next_node = hg.get_next_node(node, binary);

            for (auto& fragment : new_fragments) {

                fragment.first.push_back(hg.get_predicate_id(binary));

                fragment.second.push_back(binary);

            }
            STOP_TIMER(profile.fragment_construction_time, construction_start_2)

            previous_edges.push_back(binary);

            next_step_with_unaries(
                1 + ((number_of_random_walks - 1)/final_binaries.size()),
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

            for (auto& fragment : new_fragments) {

                fragment.first.pop_back();

                fragment.second.pop_back();

            }

            previous_edges.pop_back();
        }
    }
}

// ============================================================
// Graph Helpers
// ============================================================

vector<EdgeId>
FragmentCollection::find_cyclic_edges(
    NodeId node,
    NodeId starting_node,
    const unordered_set<NodeId>& connected_nodes,
    const boost::container::small_vector<EdgeId,3>& previous_edges,
    RecursiveProfiling& profile)
{
    START_TIMER(start)

    vector<EdgeId> cyclic_binaries;

    if (connected_nodes.find(node) != connected_nodes.end()) {

        const auto& shared_edges =
            hg.get_shared_edges(node);

        auto it =
            shared_edges.find(starting_node);

        if (it != shared_edges.end()) {

            cyclic_binaries.reserve(
                it->second.size()
            );

            for (auto edge : it->second) {

                if (!has(previous_edges, edge)) {
                    cyclic_binaries.push_back(edge);
                }
            }
        }
        
    }

    STOP_TIMER(profile.find_cyclic_edges_1, start)

    return cyclic_binaries;
}

vector<EdgeId>
FragmentCollection::get_edges_for_next_step(
    NodeId node,
    const boost::container::small_vector<EdgeId,3>& previous_edges,
    const vector<EdgeId>& forbidden_edges,
    size_t needed_edges,
    RecursiveProfiling& profile)
{
    START_TIMER(start)

    const auto& edges =
        binary_edges.at(node);

    size_t forbidden_estimate =
        previous_edges.size()
        + forbidden_edges.size();

    size_t available_estimate =
        edges.size() > forbidden_estimate
            ? edges.size() - forbidden_estimate
            : 0;

    if (needed_edges <= 4 &&
        available_estimate > needed_edges * 3)
    {
        auto result =
            direct_random_non_cyclic_selection(
                edges,
                previous_edges,
                forbidden_edges,
                needed_edges
            );

        STOP_TIMER(profile.find_cyclic_edges_2, start)

        return result;
    }

    vector<EdgeId> forbidden(
        previous_edges.begin(),
        previous_edges.end()
    );

    forbidden.insert(
        forbidden.end(),
        forbidden_edges.begin(),
        forbidden_edges.end()
    );

    vector<EdgeId> candidates = edges;

    candidates.erase(
        std::remove_if(
            candidates.begin(),
            candidates.end(),
            [&](EdgeId edge)
            {
                return has(forbidden, edge);
            }),
        candidates.end()
    );

    if (candidates.size() > needed_edges) {

        START_TIMER(random_time)

        candidates =
            get_random_elements(
                candidates,
                needed_edges
            );

        STOP_TIMER(profile.random_selection_time, random_time)

    }

    STOP_TIMER(profile.find_cyclic_edges_3, start)

    return candidates;
}

// ============================================================
// Accessors
// ============================================================

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
