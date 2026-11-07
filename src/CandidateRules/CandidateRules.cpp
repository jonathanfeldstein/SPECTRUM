#include "CandidateRules.h"

#include <utility>

bool compare_rule_utility(const Rule& rule1, const Rule& rule2){return rule1.utility>rule2.utility;}
bool compare_rule_fidelity(const Rule& rule1, const Rule& rule2){return rule1.fidelity>rule2.fidelity;}
std::string doubleToStringWithTwoDecimals(double value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << value;
    return oss.str();
}
vector<Fragment> CandidateRuleCollection::get_constrained_fragments(
    const unordered_set<PredicateId>& target_predicates
) {
    vector<const pair<const Fragment, vector<Location>>*> items;
    items.reserve(fragment_collection.fragment_count.size());

    for (auto& entry : fragment_collection.fragment_count) {
        items.push_back(&entry);
    }

    const int num_threads = omp_get_max_threads();

    vector<vector<Fragment>> local_fragments(num_threads);

    size_t estimated_per_thread =
        (items.size() / num_threads) + 1;

    for (auto& vec : local_fragments) {
        vec.reserve(estimated_per_thread);
    }

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < items.size(); ++i) {

        const int tid = omp_get_thread_num();

        const Fragment& fragment = items[i]->first;
        const vector<Location>& locations = items[i]->second;

        const bool constrained =
            (
                (
                    fragment.front() != 0 &&
                    fragment.back()  != 0 &&
                    fragment.size() > 2
                )
                ||
                ( 
                    fragment.size() == 2 &&
                    fragment.front() == 0 &&
                    fragment.back()  == 0
                )
                ||
                ( // two unary atoms on one node: unary -> unary rule
                    has_unaries_ &&
                    fragment.size() == 2 &&
                    fragment.front() != 0 &&
                    fragment.back()  != 0
                )
            )
            &&
            locations.size() > 5
            &&
            (
                has(target_predicates, fragment.front()) ||
                has(target_predicates, fragment.back())
            );

        if (constrained) {
            local_fragments[tid].push_back(fragment);
        }
    }

    unordered_set<Fragment, VectorHash> unique_fragments;

    size_t total_estimate = 0;

    for (const auto& vec : local_fragments) {
        total_estimate += vec.size();
    }

    unique_fragments.reserve(total_estimate);

    for (auto& vec : local_fragments) {
        for (auto& frag : vec) {
            unique_fragments.insert(std::move(frag));
        }
    }

    vector<Fragment> unique_framents_vector;
    unique_framents_vector.reserve(unique_fragments.size());

    for (auto& frag : unique_fragments) {
        unique_framents_vector.push_back(std::move(frag));
    }

    return unique_framents_vector;
}

vector<Fragment> CandidateRuleCollection::get_constrained_cyclic_fragments() {

    vector<const pair<const Fragment, vector<Location>>*> items;

    items.reserve(fragment_collection.fragment_cyclic_count.size());

    for (auto& entry : fragment_collection.fragment_cyclic_count) {
        items.push_back(&entry);
    }

    const int num_threads = omp_get_max_threads();

    vector<vector<Fragment>> local_fragments(num_threads);

    size_t estimated_per_thread =
        (items.size() / num_threads) + 1;

    for (auto& vec : local_fragments) {
        vec.reserve(estimated_per_thread);
    }

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < items.size(); ++i) {

        const int tid = omp_get_thread_num();

        const Fragment& fragment = items[i]->first;
        const vector<Location>& locations = items[i]->second;

        if (locations.size() >= 10) {
            local_fragments[tid].push_back(fragment);
        }
    }

    unordered_set<Fragment, VectorHash> unique_fragments;

    size_t total_estimate = 0;

    for (const auto& vec : local_fragments) {
        total_estimate += vec.size();
    }

    unique_fragments.reserve(total_estimate);

    for (auto& vec : local_fragments) {
        for (auto& frag : vec) {
            unique_fragments.insert(std::move(frag));
        }
    }

    vector<Fragment> constrained_fragments_vector;

    constrained_fragments_vector.reserve(unique_fragments.size());

    for (auto& frag : unique_fragments) {
        constrained_fragments_vector.push_back(std::move(frag));
    }

    return constrained_fragments_vector;
}


void CandidateRuleCollection::calculate_prior_probabilities(const vector<string>& target_predicates, unordered_set<PredicateId>& target_predicates_ids, unordered_map<PredicateId, double>& prior_probability) {
    unordered_map<PredicateId, double> target_counts{};
    double total_count{0};
    for(const auto& predicate:target_predicates){
        target_predicates_ids.insert(fragment_collection.get_predicate_id(predicate));
    }
    for(auto edge:fragment_collection.hg.get_predicates()){
        if(has(target_predicates_ids,edge.second)){
            target_counts[edge.second]++;
            total_count++;
        }
    }
    size_t not_empty_nodes{0};
    if (has_unaries_) {
        if(fragment_collection.is_unary(*target_predicates_ids.begin())){
            for(auto node:fragment_collection.hg.get_nodes()){
                const auto& members = fragment_collection.hg.get_memberships(node);
                for(auto member:members){
                    if(fragment_collection.is_unary_edge(member)){
                        not_empty_nodes++;
                        break;
                    }
                }
            }
        }
    }
    total_count += (double)(fragment_collection.hg.number_of_nodes() - not_empty_nodes);
    for(auto predicate:target_counts){
        prior_probability[predicate.first]=total_count/predicate.second;
    }
}

void CandidateRuleCollection::process_normal_fragments(const vector<Fragment>& constrained_fragments, unordered_map<Fragment, vector<Location>, VectorHash>& fragment_count, const unordered_set<PredicateId>& target_predicates_ids, const unordered_map<PredicateId, double>& prior_probability, double min_fidelity) {
    int num_threads = omp_get_max_threads();
    vector<unordered_map<PredicateId, vector<Rule>>> local_candidate_rules(num_threads);
    
#pragma omp parallel for schedule(dynamic)
    for(int i =0; i<constrained_fragments.size(); i++){
        int tid = omp_get_thread_num();
        Fragment fragment = constrained_fragments[i];
        auto TT_count = (double)fragment_count[fragment].size();
        //Option 1 HEAD <- BODY
        // IF HEAD IN TARGETS
        size_t symmetry_factor = 1;
        Fragment reversed_fragment = fragment;
        reverse(reversed_fragment.begin(),reversed_fragment.end());
        if(has(target_predicates_ids, fragment.front())){
            Fragment head_1{fragment.front()};
            Fragment body_1(fragment.begin() + 1, fragment.end());
            if(fragment == reversed_fragment){
                symmetry_factor = 2;
            }
            if (has_unaries_) {
                body_1.insert(body_1.begin(),0); //Replace removed unary edge by NULL Predicate
            }
            body_1 = fragment_collection.unary_layout() ? unary_utils::order_fragment(body_1) : plain_utils::order_fragment(body_1); // Sort body to canonical shape
            // Degenerate case of Unary->Unary or Binary-> Binary Rule
            if(fragment.size()==2){
                if(has_unaries_ && fragment_collection.is_unary(fragment.front())){
                    body_1 = Fragment{fragment.back()};
                }else{
                    head_1 = Fragment{0,fragment.front(),0};
                    body_1 = Fragment{0,fragment.back(), 0};
                }
            }

            double TF_count_1 = (double)fragment_count[body_1].size() - TT_count;
            double fidelity_1 = (TT_count*symmetry_factor*prior_probability.at(head_1[0]))/(TF_count_1 + TT_count);
            if(fidelity_1 > min_fidelity) {
                Rule rule_front{};
                rule_front.body = body_1;
                rule_front.head = head_1;
                rule_front.rule = fragment;
                rule_front.fidelity = fidelity_1;
                rule_front.locations = fragment_count[fragment];
                local_candidate_rules[tid][rule_front.rule[0]].push_back(std::move(rule_front));
            }
        }

        if(has(target_predicates_ids, fragment.back())){
            if(reversed_fragment != fragment){
                //Option 2 BODY -> HEAD
                Fragment rule_2 = fragment;
                reverse(rule_2.begin(),rule_2.end()); // Enforce HEAD <- BODY shape, as body might get reversed for utility calculation
                Fragment head_2{rule_2.front()};
                Fragment body_2(rule_2.begin() + 1, rule_2.end());
                if (has_unaries_) {
                    body_2.insert(body_2.begin(),0); //Replace removed unary edge by NULL Predicate
                }
                body_2 = fragment_collection.unary_layout() ? unary_utils::order_fragment(body_2) : plain_utils::order_fragment(body_2); // Sort body to canonical shape

                // Degenerate case of Unary->Unary or Binary-> Binary Rule
                if(rule_2.size()==2){
                    if(has_unaries_ && fragment_collection.is_unary(rule_2.front())){
                        body_2 = Fragment{rule_2.back()};
                    }else{
                        head_2 = Fragment{0,rule_2.front(),0};
                        body_2 = Fragment{0,rule_2.back(), 0};
                    }
                }
                double TF_count_2 = (double)fragment_count[body_2].size() - TT_count;
                double fidelity_2 = (TT_count*prior_probability.at(head_2[0]))/(TF_count_2 + TT_count);
                if(fidelity_2 > min_fidelity){
                    Rule rule_back{};
                    rule_back.body=body_2;
                    rule_back.head=head_2;
                    rule_back.rule=rule_2;
                    rule_back.fidelity=fidelity_2;
                    // Reorder locations to HEAD <- BODY shape
                    for(const auto& location: fragment_count[fragment]){
                        Location reversed_location = location;
                        reverse(reversed_location.begin(),reversed_location.end());
                        rule_back.locations.push_back(reversed_location);
                    }
                    local_candidate_rules[tid][rule_back.rule[0]].push_back(std::move(rule_back));
                }
            }
        }
    }
    for (int t = 0; t < num_threads; ++t) {
        for (auto& [pred, rules] : local_candidate_rules[t]) {
            auto& dest = candidate_rules[pred];
            dest.insert(dest.end(), rules.begin(), rules.end());
        }
    }
}

void CandidateRuleCollection::process_cyclic_fragments(const vector<Fragment>& constrained_cyclic_fragments, unordered_map<Fragment, vector<Location>, VectorHash>& fragment_count, unordered_map<Fragment, vector<Location>, VectorHash>& fragment_cyclic_count, const unordered_set<PredicateId>& target_predicates_ids, const unordered_map<PredicateId, double>& prior_probability, double min_fidelity) {
    int num_threads = omp_get_max_threads();
    vector<unordered_map<PredicateId, vector<Rule>>> local_cyclic_candidate_rules(num_threads);
#pragma omp parallel for schedule(dynamic)
    for(int j = 0; j < constrained_cyclic_fragments.size(); j++) {
        int tid = omp_get_thread_num();
        set<Fragment> tested_fragments{};
        Fragment fragment = constrained_cyclic_fragments[j];
        auto TT_count = (double)fragment_cyclic_count[fragment].size();
        auto symmetry_map = calculate_symmetry_factors(fragment);
        if(!has_unaries_ || fragment.size() % 2 == 0){
            for(int i=0; i<fragment.size(); i++) {
                if (fragment[i] != 0 && has(target_predicates_ids, fragment[i])) {
                    Fragment head_cyclic{fragment[i]};
                    Fragment body_cyclic;
                    body_cyclic.reserve(fragment.size()-1);
                    for(int j=1; j<(fragment.size()); j++){
                        body_cyclic.push_back(fragment[(i+j) % fragment.size()]);
                    }
                    Fragment reversed_body = body_cyclic;
                    std::reverse(reversed_body.begin(), reversed_body.end());
                    bool reverse_location = false;
                    if (body_cyclic > reversed_body) { // Find canonical order
                        body_cyclic = reversed_body;
                        reverse_location = true;
                    }
                    double symmetry_factor = (double)symmetry_map[{head_cyclic, body_cyclic}];
                    Fragment rule;
                    rule.reserve(fragment.size());
                    rule.push_back(head_cyclic[0]);
                    rule.insert(rule.end(), body_cyclic.begin(), body_cyclic.end());
                    if (has_unaries_ && i % 2 == 0) { // i.e. the head is a unary predicate
                        body_cyclic.push_back(0);
                        body_cyclic.insert(body_cyclic.begin(), 0);
                    }
                    double TF_count_cyclic = (double)fragment_count[body_cyclic].size() - TT_count;
                    double fidelity_cyclic = (TT_count *symmetry_factor*prior_probability.at(head_cyclic[0]))/ (TF_count_cyclic + TT_count);
                    if (fidelity_cyclic > min_fidelity && !has(tested_fragments, rule)) {
                        tested_fragments.insert(rule);
                        Rule rule_cyclic{};
                        rule_cyclic.body = body_cyclic;
                        rule_cyclic.head = head_cyclic;
                        rule_cyclic.rule = rule;
                        rule_cyclic.fidelity = fidelity_cyclic;
                        for(const auto& location:fragment_cyclic_count[fragment]){
                            Location ordered_location{};
                            Fragment ordered_fragment{};
                            if(reverse_location){
                                vector<size_t> order;
                                order.reserve(fragment.size());
                                for(int k = 0; k<fragment.size();k++){
                                    order.push_back((size_t)(i+(fragment.size()-k))%fragment.size());
                                }
                                tie(ordered_fragment, ordered_location) = fragment_collection.unary_layout() ? unary_utils::reorder_cyclic_fragment(fragment,location, {order}) : plain_utils::reorder_cyclic_fragment(fragment,location, {order});

                            }else{
                                vector<size_t> order;
                                order.reserve(fragment.size());
                                for(int k = 0; k<fragment.size();k++){
                                    order.push_back((size_t)(i+k)%fragment.size());
                                }
                                tie(ordered_fragment, ordered_location) = fragment_collection.unary_layout() ? unary_utils::reorder_cyclic_fragment(fragment,location, {order}) : plain_utils::reorder_cyclic_fragment(fragment,location, {order});
                            }
                            rule_cyclic.locations.push_back(ordered_location);
                        }
                        local_cyclic_candidate_rules[tid][rule_cyclic.rule[0]].push_back(std::move(rule_cyclic));
                    }
                }
            }
        }
    }
    for (int t = 0; t < num_threads; ++t) {
        for (auto& [pred, rules] : local_cyclic_candidate_rules[t]) {
            auto& dest = candidate_rules[pred];
            dest.insert(dest.end(), rules.begin(), rules.end());
        }
    }
}

void CandidateRuleCollection::calculate_initial_utilities() {
    Timer utility_calculation("Utility Calculation");
    for (auto& entry : candidate_rules)
    {
        auto& rules = entry.second;

        size_t heads = 0;

    #pragma omp parallel for reduction(+:heads) schedule(dynamic)
        for (size_t i = 0; i < rules.size(); ++i)
        {
            auto utility_pair =
                calculate_utility({rules[i]});

            rules[i].utility = utility_pair.first;

            heads += utility_pair.second;
        }
    }
    utility_calculation.Stop();
}

CandidateRuleCollection::CandidateRuleCollection(FragmentCollection&& fragments, 
                                                 double min_fidelity, 
                                                 const vector<string>& target_predicates,
                                                 bool has_unaries) : has_unaries_(has_unaries) {
    fragment_collection = std::move(fragments);
    
    unordered_set<PredicateId> target_predicates_ids;
    unordered_map<PredicateId, double> prior_probability;
    calculate_prior_probabilities(target_predicates, target_predicates_ids, prior_probability);
    
    Timer get_frags("Getting Fragments");
    unordered_map<Fragment, vector<Location>, VectorHash>& fragment_count = fragment_collection.get_fragment_count();
    unordered_map<Fragment, vector<Location>, VectorHash>& fragment_cyclic_count = fragment_collection.get_fragment_cyclic_count();
    cout << fragment_count.size()+fragment_cyclic_count.size() << " Fragments found" << endl;

    vector<Fragment> constrained_fragments = get_constrained_fragments(target_predicates_ids);
    vector<Fragment> constrained_cyclic_fragments = get_constrained_cyclic_fragments();
    cout << constrained_fragments.size() + constrained_cyclic_fragments.size() << " Constrained Fragments found" << endl;
    get_frags.Stop();
    
    process_normal_fragments(constrained_fragments, fragment_count, target_predicates_ids, prior_probability, min_fidelity);
    process_cyclic_fragments(constrained_cyclic_fragments, fragment_count, fragment_cyclic_count, target_predicates_ids, prior_probability, min_fidelity);
    
    Timer coverage_calculation("Coverage Calculation");
    size_t total_number_of_rules{0};
    for(auto & [predicate, rules]:candidate_rules){
        total_number_of_rules += rules.size();
    }
    cout << "Total number of rules " << total_number_of_rules << endl;
    calculate_coverage();
    coverage_calculation.Stop();
    
    calculate_initial_utilities();
}


void CandidateRuleCollection::calculate_coverage()
{
    for (auto& entry : this->candidate_rules)
    {
        auto& rules = entry.second;

#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < rules.size(); ++i)
        {
            Rule& rule = rules[i];

            rule.covered_locations.clear();
            rule.covered_locations_count.clear();

            for (const auto& location : rule.locations)
            {
                ++rule.covered_locations_count[location[0]];
                rule.covered_locations.insert(location[0]);
            }
        }
    }
}

vector<Rule> CandidateRuleCollection::get_final_rules(size_t M, const vector<string>& target_predicates, size_t max_number_of_rules_per_predicate) {
    map<PredicateId, vector<Rule>> final_rules{};
    map<PredicateId, vector<Rule>> left_candidate_rules{};
    vector<PredicateId> predicates_of_interest{};
    vector<PredicateId> informative_predicates{};
    unordered_map<PredicateId, unordered_set<PredicateId >> remaining_predicates{};
    for(auto predicate:this->fragment_collection.hg.get_predicate_name_to_id()){
        if(has(target_predicates, predicate.first)){
            predicates_of_interest.push_back(predicate.second);
        }else{
            informative_predicates.push_back(predicate.second);
        }
    }
    for(auto predicate:predicates_of_interest){
        remaining_predicates[predicate] = unordered_set<PredicateId>(informative_predicates.begin(), informative_predicates.end());
    }
    vector<Rule> final_order;
    map<PredicateId, pair<double, Rule>> next_best_candidates{}; // Predicate: {utility_increase, {candidate_rules}}
    // Initialise the final rules with at least one rule per predicate
    for(const auto& predicate:predicates_of_interest){
        for(const auto& rule:candidate_rules[predicate]){
            left_candidate_rules[predicate].push_back(rule);
        }
    }
    for(auto predicate:predicates_of_interest){
        if(!left_candidate_rules[predicate].empty()) {
            if (final_rules[predicate].empty()) {
                sort(left_candidate_rules[predicate].begin(), left_candidate_rules[predicate].end(),
                     compare_rule_fidelity);
                int number_of_rules = min((size_t)max_number_of_rules_per_predicate,
                                          left_candidate_rules[predicate].size());
                left_candidate_rules[predicate] = vector<Rule>(left_candidate_rules[predicate].begin(),
                                                               left_candidate_rules[predicate].begin() + number_of_rules);
                final_order.insert(final_order.end(),left_candidate_rules[predicate].begin(),left_candidate_rules[predicate].end());
                continue;
                double max_utility = left_candidate_rules[predicate][0].utility;
                Rule best_rule = left_candidate_rules[predicate][0];
                int rule_complexity = 0;
                int best_rule_location = 0;
                for (auto atom: best_rule.rule) {
                    if (atom != 0) {
                        rule_complexity++;
                    }
                }
                int rule_location=0;
                for (const auto &rule: left_candidate_rules[predicate]) {
                    if (rule.utility == max_utility) {
                        int new_rule_complexity = 0;
                        for (auto atom: rule.rule) {
                            if (atom != 0) {
                                new_rule_complexity++;
                            }
                        }
                        if (new_rule_complexity < rule_complexity) {
                            best_rule = rule;
                            rule_complexity = new_rule_complexity;
                            best_rule_location = rule_location;
                        }
                        rule_location++;
                    } else {
                        break;
                    }
                }
                final_rules[predicate] = vector<Rule>{best_rule};
                final_order.push_back(best_rule);
                left_candidate_rules[predicate].erase(left_candidate_rules[predicate].begin() + best_rule_location);
                for(auto remaining_predicate:best_rule.body){
                    remaining_predicates[predicate].erase(remaining_predicate);
                }
                vector<PredicateId> remaining_predicates_vector(remaining_predicates[predicate].begin(), remaining_predicates[predicate].end());
            }
            auto min_utility_pair = calculate_utility(final_rules[predicate]);
            double min_utility = min_utility_pair.first;
            double max_utility = min_utility;
            Rule most_promising_rule{};
            for (size_t i=0; i<left_candidate_rules[predicate].size(); i++) {
                const auto &rule = left_candidate_rules[predicate][i];
                vector<Rule> rule_set = final_rules[predicate];
                rule_set.push_back(rule);
                auto utility_pair = calculate_utility(rule_set);
                double utility = utility_pair.first;
                if (utility > max_utility) {
                    max_utility = utility;
                    most_promising_rule = rule;
                } else if (utility == max_utility) {
                    int new_rule_complexity = 0;
                    for (auto atom: rule.rule) {
                        if (atom != 0) {
                            new_rule_complexity++;
                        }
                    }
                    int rule_complexity = 0;
                    for (auto atom: most_promising_rule.rule) {
                        if (atom != 0) {
                            rule_complexity++;
                        }
                    }
                    if (new_rule_complexity < rule_complexity) {
                        most_promising_rule = rule;
                    }
                }
            }
            next_best_candidates[predicate] = {max_utility - min_utility, most_promising_rule};
        }
    }
    size_t N = 0;
    for(const auto& rule_sets:final_rules){
        N += rule_sets.second.size();
    }
    N=M;
    while(N < M){
        // Find the maximum utility increase
        double utility_increase{0};
        PredicateId maximal_increasing_predicate{};
        Rule maximal_increasing_rule{};
        for(const auto& candidate:next_best_candidates){
            if(candidate.second.first > utility_increase){
                maximal_increasing_predicate = candidate.first;
                maximal_increasing_rule = candidate.second.second;
                utility_increase = candidate.second.first;
            }
        }
        if(utility_increase > 0){
            final_rules[maximal_increasing_predicate].push_back(maximal_increasing_rule);
            final_order.push_back(maximal_increasing_rule);
            int j = 0;
            for(const auto& rule:left_candidate_rules[maximal_increasing_predicate]){
                if(rule.rule == maximal_increasing_rule.rule){
                    left_candidate_rules[maximal_increasing_predicate].erase(left_candidate_rules[maximal_increasing_predicate].begin()+j);
                    break;
                }
                j++;
            }
            auto min_utility_pair = calculate_utility(final_rules[maximal_increasing_predicate]);
            double min_utility = min_utility_pair.first;
            double max_utility = min_utility;
            Rule most_promising_rule{};
#pragma omp parallel for schedule(dynamic)
            for(int i=0; i<left_candidate_rules[maximal_increasing_predicate].size(); i++){
                vector<Rule> rule_set = final_rules[maximal_increasing_predicate];
                rule_set.push_back(left_candidate_rules[maximal_increasing_predicate][i]);
                auto utility_pair = calculate_utility(rule_set);
                double utility = utility_pair.first;
                if(utility > max_utility){
                    max_utility = utility;
                    most_promising_rule = left_candidate_rules[maximal_increasing_predicate][i];
                } else if (utility == max_utility) {
                    int new_rule_complexity = 0;
                    for (auto atom: left_candidate_rules[maximal_increasing_predicate][i].rule) {
                        if (atom != 0) {
                            new_rule_complexity++;
                        }
                    }
                    int rule_complexity = 0;
                    for (auto atom: most_promising_rule.rule) {
                        if (atom != 0) {
                            rule_complexity++;
                        }
                    }
                    if (new_rule_complexity < rule_complexity) {
                        most_promising_rule = left_candidate_rules[maximal_increasing_predicate][i];
                    }
                }
            }
            next_best_candidates[maximal_increasing_predicate]={max_utility-min_utility, most_promising_rule};
            N++;
        }else{
            break;
        }

    }
    return final_order;
}

pair<double, size_t>
CandidateRuleCollection::calculate_utility(
    const vector<Rule>& rules
) {
    double global_precision = 0.0;

    for (const auto& rule : rules) {
        global_precision += rule.fidelity;
    }

    const double global_complexity =
        calculate_complexity(rules);

    unordered_map<EdgeId, size_t> local_recall_degrees;

    size_t estimated_size = 0;

    for (const auto& rule : rules) {
        estimated_size += rule.covered_locations.size();
    }

    local_recall_degrees.reserve(estimated_size);

    // Aggregate coverage counts once
    for (const auto& rule : rules) {

        for (const auto& [head, count]
             : rule.covered_locations_count)
        {
            local_recall_degrees[head] += count;
        }
    }

    double global_recall = 0.0;

    for (const auto& [head, local_recall_degree]
         : local_recall_degrees)
    {
        global_recall +=
            log(1.0 + (double)local_recall_degree);
    }

    const double global_utility =
        global_precision *
        global_recall *
        global_complexity;

    return {
        global_utility,
        local_recall_degrees.size()
    };
}

void CandidateRuleCollection::print_rules(const string& logic,
                                          const string& output_file_name,
                                          size_t M,
                                          const vector<string>& target_predicates,
                                        size_t max_number_of_rules_per_predicate) {

    string model_filename = output_file_name;
    if(logic == "psl") {
        model_filename += ".json";
    }else{
        model_filename += "." + logic;
    }
    ofstream file;
    file.open(model_filename);
    
    // Add stats file output
    ofstream stats_file;
    stats_file.open(output_file_name + "_stats.csv");
    stats_file << "rule,utility,fidelity,recall\n";

    vector<Rule> rules = get_final_rules(M, target_predicates, max_number_of_rules_per_predicate);
    map<PredicateId,vector<Rule>> utility_normalisation_map{};
    for(auto rule:rules){
        utility_normalisation_map[rule.head[0]].push_back(rule);
    }
    for(auto& [predicate, normalised_rules]:utility_normalisation_map){
        double utility_sum = 0;
        for(const auto& rule:normalised_rules){
            utility_sum += rule.utility;
        }
        for(auto& rule:normalised_rules){
            rule.utility = rule.utility/utility_sum;
        }
    }
    for(auto& rule:rules){
        for(auto&  normalised_rule:utility_normalisation_map[rule.head[0]]){
            if(normalised_rule.rule == rule.rule){
                rule.utility = normalised_rule.utility;
            }
        }
    }
    if(file.is_open()){
        if(logic == "psl"){
            file << "{\n";
            file << "    \"rules\": [\n";
        }
        for(int i =0; i< rules.size(); i++){
            string printable_rule = print_rule(rules[i], logic);
            if(logic == "psl") {
                file << "        " << printable_rule;
                if (i != rules.size() - 1) {
                    file << ",\n";
                } else {
                    file << "\n";
                }
            }else{
                file << printable_rule;
                if(i != rules.size()-1){
                    file << "\n";
                }
            }
            // Strip out the `"fidelity: ` prefix for CSV if it's psl
            string csv_rule = printable_rule;
            if (logic == "psl") {
                size_t colon_pos = printable_rule.find(": ");
                if (colon_pos != string::npos) {
                    csv_rule = printable_rule.substr(colon_pos + 2);
                    // Also strip the trailing `^2"`
                    size_t hat_pos = csv_rule.rfind(" ^2\"");
                    if (hat_pos != string::npos) {
                        csv_rule = csv_rule.substr(0, hat_pos);
                    }
                }
            }
            // Recall as in calculate_utility: sum over covered head facts of ln(1 + groundings entailing it)
            double recall = 0.0;
            for (const auto& [head, count] : rules[i].covered_locations_count) {
                recall += log(1.0 + (double)count);
            }
            stats_file << "\"" << csv_rule << "\"," << rules[i].utility << "," << rules[i].fidelity << "," << recall << "\n";
        }
        if(logic == "psl"){
            file << "    ]\n";
            file << "}\n";
        }
    }else{
        throw FileNotOpenedException(model_filename);
    }
    stats_file.close();
}

double CandidateRuleCollection::compute_local_utility(
    EdgeId head,
    const vector<Rule>& rules
) {

    double total_recall_degree = 0.0;
    double total_local_utility = 0.0;
    double normalization_factor = 0.0;

    for (const auto& rule : rules) {

        auto covered_it =
            rule.covered_locations.find(head);

        if (covered_it == rule.covered_locations.end()) {
            continue;
        }

        auto recall_it =
            rule.covered_locations_count.find(head);

        if (recall_it ==
            rule.covered_locations_count.end()) {
            continue;
        }

        auto structure_penalty_it =
            rule.local_structure_penalty.find(head);

        if (structure_penalty_it ==
            rule.local_structure_penalty.end()) {
            continue;
        }

        const double recall_degree =
            static_cast<double>(recall_it->second);

        total_local_utility +=
            rule.fidelity * structure_penalty_it->second;

        total_recall_degree += recall_degree;

        normalization_factor += 1.0;
    }

    if (normalization_factor == 0.0) {
        return 0.0;
    }

    const double global_coverage_mass =
        log(1.0 + total_recall_degree);

    return
        (global_coverage_mass * total_local_utility)
        / normalization_factor;
}

string CandidateRuleCollection::print_rule(Rule rule, const string& logic) {
    map<NodeId, int> nodes_to_var{};
    vector<EdgeId> location = *rule.locations.begin();
    int var_counter = 0;
    for(int i = (int)location.size()-1; i>-1; i--){
        const vector<NodeId>& nodes = this->fragment_collection.get_edge(location[i]);
        for(auto node:nodes){
            if (nodes_to_var.find(node) == nodes_to_var.end()) {
                nodes_to_var[node] = var_counter;
                var_counter++;
            }
        }
    }

    std::vector<std::string> vars(var_counter);
    for (int i = 0; i < var_counter; ++i) {
        if (logic == "psl") {
            if (i < 26) {
                vars[i] = std::string(1, 'A' + i);
            } else {
                vars[i] = "V" + std::to_string(i);
            }
        } else {
            vars[i] = "v" + std::to_string(i);
        }
    }
    
    string printable_rule{};
    int k = (int)location.size()-1;
    if(logic == "psl"){
        printable_rule += "\"" + doubleToStringWithTwoDecimals(rule.fidelity) + ": ";
    }
    for(int i = (int)rule.rule.size() - 1; i > 0; i--) {
        if (rule.rule[i] != 0) {
            if(logic != "psl"){
                printable_rule += "!";
            }
            std::string predicate = this->fragment_collection.get_predicate_from_id(rule.rule[i]);
            size_t pos = predicate.find_last_not_of("0123456789");
            std::string clean_predicate = predicate.substr(0, pos + 1);
            std::string constant = predicate.substr(pos + 1);
            printable_rule += clean_predicate;
            printable_rule += "(";
            if (has_unaries_ && this->fragment_collection.is_unary(rule.rule[i])) {
                printable_rule += vars[nodes_to_var[this->fragment_collection.get_edge(location[k])[0]]];
            } else {
                printable_rule += vars[nodes_to_var[this->fragment_collection.get_edge(location[k])[0]]];
                printable_rule += ",";
                printable_rule += vars[nodes_to_var[this->fragment_collection.get_edge(location[k])[1]]];
            }
            if(!constant.empty()){
                printable_rule += ",'" + constant + "'";
            }
            printable_rule += ")";
            k--;
            if (i > 1) {
                if (has_unaries_) {
                    for(int j=i-1; j>0; j--){
                        if(rule.rule[j] != 0){
                            if(logic == "psl")
                                printable_rule += " & ";
                            else
                                printable_rule += " v ";
                            break;
                        }
                    }
                } else {
                    if(logic == "psl")
                        printable_rule += " & ";
                    else
                        printable_rule += " v ";
                }
            }
        }
    }
    if(logic == "psl")
        printable_rule += " >> ";
    else
        printable_rule += " v ";
    std::string predicate = this->fragment_collection.get_predicate_from_id(rule.rule[0]);
    size_t pos = predicate.find_last_not_of("0123456789");
    std::string clean_predicate = predicate.substr(0, pos + 1);
    std::string constant = predicate.substr(pos + 1);
    printable_rule += clean_predicate;
    printable_rule += "(";
    if (has_unaries_ && this->fragment_collection.is_unary(rule.rule[0])) {
        NodeId node = nodes_to_var[this->fragment_collection.get_edge(location[k])[0]];
        printable_rule += vars[node];
    } else {
        printable_rule += vars[nodes_to_var[this->fragment_collection.get_edge(location[k])[0]]];
        printable_rule += ",";
        printable_rule += vars[nodes_to_var[this->fragment_collection.get_edge(location[k])[1]]];
    }
    if(!constant.empty()){
        printable_rule += ",'" + constant + "'";
    }
    printable_rule += ")";
    if(logic == "psl"){
        printable_rule += " ^2\"";
    }
    return printable_rule;
}

void CandidateRuleCollection::calculate_structure_penalties()
{
    for (auto& entry : this->candidate_rules)
    {
        auto& rules = entry.second;

#pragma omp parallel for schedule(static)
        for (size_t j = 0; j < rules.size(); ++j)
        {
            Rule& rule = rules[j];

            unordered_map<EdgeId, vector<double>>
                local_structure_penalty;

            local_structure_penalty.reserve(
                rule.locations.size()
            );

            if (rule.rule.size() == 2) {

                for (const auto& location : rule.locations) {

                    double structure_penalty;

                    const bool unary_edge =
                        this->fragment_collection
                        .is_unary_edge(location[0]);

                    const vector<NodeId>& nodes =
                        this->fragment_collection
                        .hg.get_edge(location[0]);

                    if (has_unaries_ && unary_edge) {

                        auto [unaries, binaries] =
                            this->fragment_collection
                            .get_unary_binary_edges(nodes[0]);

                        structure_penalty =
                            2.0 / (double)unaries.size();

                    } else {

                        auto [unaries_1, binaries_1] =
                            this->fragment_collection
                            .get_unary_binary_edges(nodes[0]);

                        auto [unaries_2, binaries_2] =
                            this->fragment_collection
                            .get_unary_binary_edges(nodes[1]);

                        structure_penalty =
                            1.0 /
                            (
                                (double)binaries_2.size()
                                +
                                1.0 /
                                (double)binaries_1.size()
                            );
                    }

                    local_structure_penalty
                        [location[0]]
                        .push_back(structure_penalty);
                }

            } else {

                for (const auto& location : rule.locations) {

                    double structure_penalty = 1.0;

                    if (has_unaries_) {
                        if (
                            rule.rule.size() == 6 ||
                            rule.rule.size() == 8
                        ) {
                            structure_penalty = 2.0;
                        }
                    } else {
                        if (
                            rule.rule.size() == 3 ||
                            rule.rule.size() == 4
                        ) {
                            structure_penalty = 2.0;
                        }
                    }

                    NodeId starting_node =
                        this->fragment_collection
                        .get_edge(location[0])[0];

                    NodeId current_node =
                        starting_node;

                    bool last_edge_binary = false;

                    for (
                        size_t i = 0;
                        i < location.size() - 1;
                        ++i
                    ) {

                        if (
                            i == 0 &&
                            (!has_unaries_ || !this->fragment_collection
                            .is_unary_edge(location[0]))
                        ) {

                            last_edge_binary = true;

                            const auto& nodes_1 =
                                this->fragment_collection
                                .hg.get_edge(location[0]);

                            const auto& nodes_2 =
                                this->fragment_collection
                                .hg.get_edge(location[1]);

                            current_node =
                                (
                                    nodes_1[1] == nodes_2[0] ||
                                    nodes_1[1] == nodes_2[1]
                                )
                                ? nodes_1[1]
                                : nodes_1[0];
                        }

                        auto [unaries, binaries] =
                            this->fragment_collection
                            .get_unary_binary_edges(
                                current_node
                            );

                        const size_t number_of_binaries =
                            binaries.size();

                        if (number_of_binaries == 0) {
                            continue;
                        }

                        const size_t number_of_unaries =
                            unaries.size();

                        const bool unary_edge =
                            this->fragment_collection
                            .is_unary_edge(location[i]);

                        if (has_unaries_ && unary_edge) {

                            structure_penalty *=
                                2.0 /
                                (double)(
                                    number_of_unaries + 1
                                );

                        } else {

                            if (
                                current_node ==
                                starting_node
                            ) {

                                structure_penalty *=
                                    1.0 /
                                    (double)(
                                        number_of_binaries
                                    );

                            } else {

                                structure_penalty *=
                                    2.0 /
                                    (double)(
                                        number_of_binaries
                                    );
                            }

                            if (has_unaries_ && last_edge_binary) {

                                structure_penalty *=
                                    1.0 /
                                    (double)(
                                        number_of_unaries + 1
                                    );
                            }

                            last_edge_binary = true;

                            current_node =
                                this->fragment_collection
                                .hg.get_next_node(
                                    current_node,
                                    location[i]
                                );
                        }
                    }

                    local_structure_penalty
                        [location[0]]
                        .push_back(structure_penalty);
                }
            }

            rule.local_structure_penalty.clear();

            for (
                const auto& [edge, penalties]
                : local_structure_penalty
            ) {

                const double average =
                    accumulate(
                        penalties.begin(),
                        penalties.end(),
                        0.0
                    ) /
                    (double)penalties.size();

                rule.local_structure_penalty[edge] =
                    average;
            }
        }
    }
}


double CandidateRuleCollection::calculate_complexity(
    const vector<Rule>& rules
) {

    double complexity = 1.0;

    for (const auto& rule : rules) {

        complexity *= exp(
            -static_cast<double>(
                rule.locations[0].size()
            )
        );
    }

    return pow(complexity, 1.0 / rules.size());
}

map<pair<Fragment, Fragment>, size_t> CandidateRuleCollection::calculate_symmetry_factors(const Fragment& fragment) {
    map<pair<Fragment, Fragment>, size_t> symmetry_map;
    
    for (size_t i = 0; i < fragment.size(); ++i) {
        if(fragment[i] == 0){
            continue;
        }
        Fragment head{fragment[i]};
        Fragment body;
        body.reserve(fragment.size() - 1);
        for (size_t k = 1; k < fragment.size(); ++k) {
            body.push_back(fragment[(i + k) % fragment.size()]);
        }
        
        Fragment body_rev = body;
        std::reverse(body_rev.begin(), body_rev.end());
        
        Fragment canonical_body = (body > body_rev) ? body_rev : body;
        
        symmetry_map[{head, canonical_body}]++;
    }
    
    return symmetry_map;
}
