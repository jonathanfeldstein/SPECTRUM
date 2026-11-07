#include "Graph.hpp"
#include "RuleParser.hpp"
#include "Evaluator.hpp"
#include <iostream>
#include <iomanip>
#include <string>

// Evaluates SPECTRUM's PSL rules (rules/psl/train.json) on a dataset's test split using the
// NCRL protocol. Results are identical to NCRL's kg_completion.py run on the same rules.
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <dataset_dir> <rules_json_path> [top_k=200]\n";
        return 1;
    }

    std::string dataset_dir = argv[1];
    std::string rules_json_path = argv[2];
    int top_k = argc > 3 ? std::stoi(argv[3]) : 200;

    Graph graph;
    graph.load_dataset(dataset_dir);

    auto rules = RuleParser::parse(rules_json_path, top_k, graph.train_relations);

    EvalResult res = Evaluator::evaluate(graph, rules);

    std::cout << "Rules evaluated: " << res.num_rules << " (top " << top_k << " per head relation; skipped "
              << RuleParser::skipped_non_chain << " non-chain, " << RuleParser::skipped_relation
              << " with unsupported relations)\n";
    std::cout << "Scored queries: " << res.scored_queries << " / " << res.total_queries
              << " (queries whose relation has no rule are not scored)\n";
    std::cout << std::fixed << std::setprecision(4) << "\n";

    std::cout << "MRR:     " << res.mrr << "\n";
    for(int k=1; k<=10; ++k) {
        std::cout << "Hit@" << k << ":   " << res.hits[k-1] << "\n";
    }

    return 0;
}
