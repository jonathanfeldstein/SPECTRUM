#pragma once
#include "Graph.hpp"
#include "RuleParser.hpp"
#include <vector>
#include <omp.h>
#include <algorithm>
#include <iostream>

struct EvalResult {
    double mrr = 0.0;
    std::vector<double> hits = std::vector<double>(10, 0.0);  // hits[k-1] = Hit@k
    int num_rules = 0;
    int scored_queries = 0;
    int total_queries = 0;
};

class Evaluator {
public:
    static EvalResult evaluate(Graph& graph, const std::vector<Rule>& rules) {
        std::unordered_map<int, std::vector<Rule>> rules_by_head;
        for (const auto& r : rules) {
            if (graph.rel2idx.count(r.head_rel)) {
                rules_by_head[graph.rel2idx[r.head_rel]].push_back(r);
            }
        }

        double sum_mrr = 0;
        double sum_hits[10] = {0};
        int evaluated_count = 0;

        int num_entities = graph.ent2idx.size();
        int total_queries = graph.eval_queries.size();

        #pragma omp parallel for reduction(+:sum_mrr, sum_hits, evaluated_count)
        for (int i = 0; i < total_queries; ++i) {
            const auto& query = graph.eval_queries[i];
            int q_h = query.h;
            int q_r = query.r;
            int q_t = query.t;

            if (rules_by_head.find(q_r) == rules_by_head.end()) continue;

            std::vector<double> scores(num_entities, 0.0);

            for (const auto& rule : rules_by_head[q_r]) {
                // Count body paths first, then weight them: NCRL computes
                // conf * (M_1 * ... * M_n) and sums rules in order, and the
                // exact same floating-point operations are needed for ties to agree.
                std::unordered_map<int, double> current_reach;
                current_reach[q_h] = 1.0;

                for (const std::string& b_rel_str : rule.body_rels) {
                    if (graph.rel2idx.find(b_rel_str) == graph.rel2idx.end()) {
                        current_reach.clear();
                        break;
                    }
                    int b_rel = graph.rel2idx[b_rel_str];
                    std::unordered_map<int, double> next_reach;
                    for (auto const& pair : current_reach) {
                        int node = pair.first;
                        double val = pair.second;
                        if (b_rel < graph.adj.size() && node < graph.adj[b_rel].size()) {
                            for (int neighbor : graph.adj[b_rel][node]) {
                                next_reach[neighbor] += val;
                            }
                        }
                    }
                    current_reach = std::move(next_reach);
                }

                for (auto const& pair : current_reach) {
                    scores[pair.first] += pair.second * rule.confidence;
                }
            }

            int S_better = 0;
            double target_score = scores[q_t];

            for (int e = 0; e < num_entities; ++e) {
                if (e == q_t) continue;
                if (graph.gt.count({q_h, q_r, e})) continue; // filter
                
                if (scores[e] > target_score) S_better++;
            }

            // Rank = 1 + number of (filtered) candidates scoring strictly higher, as in NCRL
            double rank = S_better + 1.0;
            sum_mrr += 1.0 / rank;
            for(int k=1; k<=10; ++k) {
                if (rank <= k) sum_hits[k-1] += 1.0;
            }

            evaluated_count++;
        }

        EvalResult res;
        res.num_rules = rules.size();
        res.scored_queries = evaluated_count;
        res.total_queries = total_queries;
        if (evaluated_count > 0) {
            res.mrr = sum_mrr / evaluated_count;
            for(int k=0; k<10; ++k) {
                res.hits[k] = sum_hits[k] / evaluated_count;
            }
        }
        return res;
    }
};
