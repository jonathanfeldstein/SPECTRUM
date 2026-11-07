#ifndef SPECTRUM_MAP_UTILS_H
#define SPECTRUM_MAP_UTILS_H

#include <map>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include "exceptions.h"
#include <functional>
#include <boost/container/small_vector.hpp>

using namespace std;

struct VectorHash
{
    template<typename Container>
    size_t operator()(
        const Container& v) const
    {
        size_t h = 0;

        for (const auto& x : v)
        {
            h ^=
                std::hash<
                    typename Container::value_type
                >{}(x)
                +
                0x9e3779b9
                +
                (h << 6)
                +
                (h >> 2);
        }

        return h;
    }
};
template<typename K, typename V, typename Hash>
void parallel_merge_maps(
    std::vector<std::unordered_map<K,std::vector<V>,Hash>>& local_maps,
    std::unordered_map<K,std::vector<V>,Hash>& global_map,
    int num_threads)
{
    //------------------------------------------------------------------
    // key -> threads containing key
    //------------------------------------------------------------------

    std::unordered_map<K,std::vector<int>,Hash> key_threads;

    size_t estimated_keys = 0;

    for (int t = 0; t < num_threads; ++t) {
        estimated_keys += local_maps[t].size();
    }

    key_threads.reserve(estimated_keys);

    for (int t = 0; t < num_threads; ++t) {

        for (auto& entry : local_maps[t]) {

            auto& thread_list =
                key_threads[entry.first];

            thread_list.push_back(t);
        }
    }

    //------------------------------------------------------------------
    // materialise key list
    //------------------------------------------------------------------

    std::vector<K> keys;

    keys.reserve(key_threads.size());

    global_map.reserve(keys.size());


    for (auto& entry : key_threads) {

        keys.push_back(entry.first);

        global_map.try_emplace(entry.first);
    }


    //------------------------------------------------------------------
    // parallel merge
    //------------------------------------------------------------------

#pragma omp parallel for schedule(static)

    for (size_t i = 0; i < keys.size(); ++i) {

        const K& key = keys[i];

        auto& target =
            global_map.find(key)->second;

        const auto& threads =
            key_threads.find(key)->second;

        //--------------------------------------------------------------
        // total size
        //--------------------------------------------------------------

        size_t total_size = 0;

        for (int tid : threads) {

            auto it =
                local_maps[tid].find(key);

            total_size += it->second.size();
        }

        target.reserve(total_size);

        //--------------------------------------------------------------
        // move
        //--------------------------------------------------------------

        for (int tid : threads) {

            auto it =
                local_maps[tid].find(key);

            target.insert(
                target.end(),
                std::make_move_iterator(
                    it->second.begin()
                ),
                std::make_move_iterator(
                    it->second.end()
                )
            );
        }
    }
}


struct PairHash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

template<typename TK, typename TV>
static set<TK> get_keys(map<TK, TV> const& input_map) {
    set<TK> retval;
    for (auto const& element : input_map) {
        retval.insert(element.first);
    }
    return retval;
}
template<typename TK, typename TV>
static unordered_set<TK> get_keys(unordered_map<TK, TV> const& input_map) {
    unordered_set<TK> retval;
    for (auto const& element : input_map) {
        retval.insert(element.first);
    }
    return retval;
}
template<typename TK, typename TV>
static std::vector<TV> get_values(std::map<TK, TV> const& input_map) {
    vector<TV> retval;
    for (auto const& element : input_map) {
        retval.push_back(element.second);
    }
    return retval;
}

template<typename TK, typename TV>
static TV get_value(std::map<TK, TV> const& input_map, TK const&key) {
    try {
        TV value = input_map.at(key);
        return value;
    }
    catch (const std::out_of_range&) {
        throw KeyNotFoundException(key);
    }
}

bool cmp(const pair<string, int> &a,
         const pair<string, int> &b);

vector<pair<string, int> > sort(map<string, int>& M);

vector<pair<string, double> > normalize_and_sort_counts(map<string, int>& M, size_t normalization_constant);
#endif //SPECTRUM_MAP_UTILS_H
