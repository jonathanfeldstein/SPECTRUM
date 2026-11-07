#ifndef SPECTRUM_RELATION_H
#define SPECTRUM_RELATION_H

#include <vector>
#include <string>
#include <set>

using Predicate = std::string;

struct GroundRelation{
    Predicate predicate;
    std::vector<std::string> arguments;
    double weight;
};

struct Relation{
    Predicate predicate;
    std::vector<std::string> arguments;
};
#endif //SPECTRUM_RELATION_H
