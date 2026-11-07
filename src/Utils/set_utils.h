#ifndef SPECTRUM_SET_UTILS_H
#define SPECTRUM_SET_UTILS_H

#include <set>
#include <algorithm>

using namespace std;

template<typename T>
bool has(const set<T>& checked_set, const T& element){
    return checked_set.find(element) != checked_set.end();
}
template<typename T>
bool has(const unordered_set<T>& checked_set, const T& element){
    return checked_set.find(element) != checked_set.end();
}
#endif //SPECTRUM_SET_UTILS_H
