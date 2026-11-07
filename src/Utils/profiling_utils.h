#ifndef SPECTRUM_PROFILING_UTILS_H
#define SPECTRUM_PROFILING_UTILS_H

#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <omp.h>

#ifdef ENABLE_DEBUG

constexpr bool debug = true;

#define DEBUG_ONLY(...) __VA_ARGS__

#define START_TIMER(name) \
    double name = omp_get_wtime();

#define STOP_TIMER(target,name) \
    target += omp_get_wtime() - name;

#else

constexpr bool debug = false;

#define DEBUG_ONLY(...)

#define START_TIMER(name)

#define STOP_TIMER(target,name)

#endif

using namespace std;

class Timer{
public:
    Timer(string function_name){
        function = function_name;
        startTime = std::chrono::high_resolution_clock::now();
    }

    ~Timer(){
    }

    void Stop(){
        auto endTime = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(startTime).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
        auto duration = (double)(end - start)/1000000.0;
        cout << function
             << " (Wall Time) : "
             << duration
             << std::setprecision(3)
             << " s"
             << endl;
    }

private:
    string function;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
};

struct FragmentProfiling {

    double initial_time = 0.0;
    double recursive_time = 0.0;

    size_t nodes_processed = 0;
    size_t recursive_calls = 0;
};

struct RecursiveProfiling {

    double storing_fragments_time = 0.0;
    double storing_cyclic_time = 0.0;
    double graph_lookup_time = 0.0;
    double fragment_construction_time = 0.0;
    double subpattern_create_time = 0.0;
    double subpattern_insert_time = 0.0;
    double edge_filter_time = 0.0;
    double random_selection_time = 0.0;
    double find_cyclic_edges = 0.0;
    double find_cyclic_edges_1 = 0.0;
    double find_cyclic_edges_2 = 0.0;
    double find_cyclic_edges_3 = 0.0;

    size_t fragments_stored = 0;
    size_t cyclic_fragments_stored = 0;
};

inline void print_profiling(
    const vector<FragmentProfiling>& thread_profiles,
    const vector<RecursiveProfiling>& recursive_profiles,
    double canonical_time)
{
#ifdef ENABLE_DEBUG

    FragmentProfiling total{};

    RecursiveProfiling total_recursive{};

    for (const auto& p : thread_profiles)
    {
        total.initial_time += p.initial_time;
        total.recursive_time += p.recursive_time;
        total.nodes_processed += p.nodes_processed;
        total.recursive_calls += p.recursive_calls;
    }

    for (const auto& p : recursive_profiles)
    {
        total_recursive.graph_lookup_time +=
            p.graph_lookup_time;

        total_recursive.edge_filter_time +=
            p.edge_filter_time;

        total_recursive.storing_cyclic_time +=
            p.storing_cyclic_time;

        total_recursive.storing_fragments_time +=
            p.storing_fragments_time;

        total_recursive.random_selection_time +=
            p.random_selection_time;

        total_recursive.find_cyclic_edges +=
            p.find_cyclic_edges;

        total_recursive.fragment_construction_time +=
            p.fragment_construction_time;

        total_recursive.fragments_stored +=
            p.fragments_stored;

        total_recursive.cyclic_fragments_stored +=
            p.cyclic_fragments_stored;

        total_recursive.find_cyclic_edges_1 +=
            p.find_cyclic_edges_1;

        total_recursive.find_cyclic_edges_2 +=
            p.find_cyclic_edges_2;

        total_recursive.find_cyclic_edges_3 +=
            p.find_cyclic_edges_3;
    }

    std::cout
        << "\n=== Fragment Profiling ===\n";

    std::cout
        << "Nodes processed: "
        << total.nodes_processed
        << "\n";

    std::cout
        << "Recursive calls: "
        << total.recursive_calls
        << "\n";

    std::cout
        << "Initial setup time: "
        << total.initial_time
        << " s\n";

    std::cout
        << "Recursive time (CPU time): "
        << total.recursive_time
        << " s\n";

    std::cout
        << "Canonicalization time (Wall time): "
        << canonical_time
        << " s\n";

    std::cout
        << "=========================\n";

    std::cout
        << "\n=== Recursive Breakdown ===\n";

    std::cout
        << "storing_fragments time: "
        << total_recursive.storing_fragments_time
        << "\n";

    std::cout
        << "storing cyclic fragments time: "
        << total_recursive.storing_cyclic_time
        << "\n";

    std::cout
        << "random selection time (CPU time): "
        << total_recursive.random_selection_time
        << "\n";

    std::cout
        << "finding edges time (CPU time): "
        << total_recursive.edge_filter_time
        << "\n";

    std::cout
        << "finding cyclic edges time (CPU time): "
        << total_recursive.find_cyclic_edges
        << "\n";

    std::cout
        << "finding cyclic edges part 1 (CPU time): "
        << total_recursive.find_cyclic_edges_1
        << "\n";

    std::cout
        << "finding cyclic edges part 2 (CPU time): "
        << total_recursive.find_cyclic_edges_2
        << "\n";

    std::cout
        << "finding cyclic edges part 3 (CPU time): "
        << total_recursive.find_cyclic_edges_3
        << "\n";

    std::cout
        << "graph lookups (CPU time): "
        << total_recursive.graph_lookup_time
        << " s\n";

    std::cout
        << "fragment construction (CPU time): "
        << total_recursive.fragment_construction_time
        << " s\n";

    std::cout
        << "fragments stored: "
        << total_recursive.fragments_stored
        << "\n";

    std::cout
        << "cyclic fragments stored: "
        << total_recursive.cyclic_fragments_stored
        << "\n";

    std::cout
        << "============================\n";

#else

    (void)thread_profiles;
    (void)recursive_profiles;
    (void)canonical_time;

#endif
}

#endif // SPECTRUM_PROFILING_UTILS_H