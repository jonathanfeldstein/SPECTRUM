#ifndef SPECTRUM_CANDIDATERULES_H
#define SPECTRUM_CANDIDATERULES_H

#ifdef USE_RANDOM_WALKS_V3
#include "../Fragments/Fragments_v3.h"
#elif defined(USE_RANDOM_WALKS)
#include "../Fragments/Fragments_rw.h"
#else
#include "../Fragments/Fragments.h"
#endif
#include "../Utils/profiling_utils.h"
#include <sstream>
#include <iomanip>

using PredicateId = std::size_t ;
using Fragment = std::vector<PredicateId>;
using Location = std::vector<EdgeId>;


struct Rule{
    Fragment head;
    Fragment body;
    Fragment rule; //Reordered to HEAD<-BODY
    vector<Location> locations;
    unordered_map<NodeId, size_t> covered_locations_count; // The nodes that it predicts (i.e. different heads)
    unordered_set<NodeId> covered_locations; // The maximum number of nodes the rule could have predicted
    double utility;
    double fidelity;
    unordered_map<EdgeId, double> local_structure_penalty;
};


class CandidateRuleCollection{
private:
    bool has_unaries_;
    FragmentCollection fragment_collection;
    unordered_map<PredicateId, vector<Rule>> candidate_rules; // KEY Rule in Head <- Body shape
    HyperGraph hg;
    unordered_map<PredicateId, Predicate> predicate_id_to_name;
    unordered_map<Predicate, PredicateId> predicate_name_to_id;
    unordered_map<EdgeId, vector<NodeId>> relations;

    void calculate_prior_probabilities(const vector<string>& target_predicates, unordered_set<PredicateId>& target_predicates_ids, unordered_map<PredicateId, double>& prior_probability);
    void process_normal_fragments(const vector<Fragment>& constrained_fragments, unordered_map<Fragment, vector<Location>, VectorHash>& fragment_count, const unordered_set<PredicateId>& target_predicates_ids, const unordered_map<PredicateId, double>& prior_probability, double min_fidelity);
    void process_cyclic_fragments(const vector<Fragment>& constrained_cyclic_fragments, unordered_map<Fragment, vector<Location>, VectorHash>& fragment_count, unordered_map<Fragment, vector<Location>, VectorHash>& fragment_cyclic_count, const unordered_set<PredicateId>& target_predicates_ids, const unordered_map<PredicateId, double>& prior_probability, double min_fidelity);
    void calculate_initial_utilities();
public:
    CandidateRuleCollection(FragmentCollection&& fragment_collection, double min_fidelity, const vector<string>& target_predicates, bool has_unaries);
    vector<Fragment> get_constrained_fragments(const unordered_set<PredicateId>& target_predicates);
    vector<Fragment> get_constrained_cyclic_fragments();
    void calculate_coverage();
    vector<Rule> get_final_rules(size_t M, const vector<string>& target_predicates, size_t max_number_of_rules_per_predicate);
    static pair<double,size_t>  calculate_utility(const vector<Rule>& rules);
    static double calculate_complexity(const vector<Rule>& rules);
    void print_rules(const string& logic,
                     const string& output_file_name,
                     size_t M,
                     const vector<string>& target_predicates,
                    size_t max_number_of_rules_per_predicate);
    string print_rule(Rule rule, const string& logic);
    void calculate_structure_penalties();
    static double compute_local_utility(EdgeId head, const vector<Rule>& rules);
    static std::map<std::pair<Fragment, Fragment>, size_t> calculate_symmetry_factors(const Fragment& fragment);
};
std::string doubleToStringWithTwoDecimals(double value);
#endif //SPECTRUM_CANDIDATERULES_H
