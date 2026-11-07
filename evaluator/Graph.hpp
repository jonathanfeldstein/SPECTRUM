#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

struct Triple {
    int h;
    int r;
    int t;
    bool operator==(const Triple& other) const {
        return h == other.h && r == other.r && t == other.t;
    }
};

namespace std {
    template<> struct hash<Triple> {
        size_t operator()(const Triple& t) const {
            return std::hash<int>()(t.h) ^ (std::hash<int>()(t.r) << 1) ^ (std::hash<int>()(t.t) << 2);
        }
    };
}

// Loads a dataset folder (entities.txt, relations.txt, facts.txt, train.txt, valid.txt, test.txt;
// one tab-separated "a<TAB>r<TAB>b" triple per line, i.e. r(a, b) in train.db) following the
// NCRL evaluation protocol (Cheng et al., ICLR 2023):
//  - a line "a r b" is the query (h=b, r, t=a): the evaluator predicts a given r and b
//  - rule bodies are matched over the path graph facts + inverse(facts) + train + valid
//  - edges are deduplicated
//  - the filter set is facts + inverse(facts) + train + valid + test
class Graph {
public:
    std::unordered_map<std::string, int> ent2idx;
    std::unordered_map<std::string, int> rel2idx;
    std::vector<std::string> idx2rel;

    // relations occurring in train.txt; rules using any other relation are not evaluated
    std::unordered_set<std::string> train_relations;

    // evaluation ground truth filter set
    std::unordered_set<Triple> gt;

    // adjacency list: rel -> head -> sorted unique tails
    std::vector<std::vector<std::vector<int>>> adj;

    std::vector<Triple> eval_queries;

    void load_dataset(const std::string& data_dir) {
        std::string root = data_dir + "/";

        for (const auto& e : read_lines(root + "entities.txt")) {
            ent2idx[e] = ent2idx.size();
        }

        std::vector<std::string> base_rels;
        for (const auto& r : read_lines(root + "relations.txt")) {
            if (rel2idx.count(r)) continue;
            rel2idx[r] = rel2idx.size();
            base_rels.push_back(r);
        }
        idx2rel = base_rels;
        for (const auto& r : base_rels) {
            rel2idx["inv_" + r] = rel2idx.size();
            idx2rel.push_back("inv_" + r);
        }
        adj.resize(rel2idx.size());

        auto facts = load_triples(root + "facts.txt");
        auto train = load_triples(root + "train.txt");
        auto valid = load_triples(root + "valid.txt");
        auto test = load_triples(root + "test.txt");

        std::vector<Triple> inv_facts;
        for (const auto& x : facts) inv_facts.push_back({x.t, inverse(x.r), x.h});

        for (const auto& x : train) train_relations.insert(idx2rel[x.r]);

        for (const auto* part : {&facts, &inv_facts, &train, &valid}) {
            for (const auto& x : *part) {
                if (adj[x.r].size() <= (size_t)x.h) adj[x.r].resize(x.h + 1);
                adj[x.r][x.h].push_back(x.t);
            }
        }
        for (auto& by_head : adj) {
            for (auto& tails : by_head) {
                std::sort(tails.begin(), tails.end());
                tails.erase(std::unique(tails.begin(), tails.end()), tails.end());
            }
        }

        for (const auto* part : {&facts, &inv_facts, &train, &valid, &test}) {
            for (const auto& x : *part) gt.insert(x);
        }

        eval_queries = test;

        std::cout << "Loaded " << ent2idx.size() << " entities and " << base_rels.size() << " relations.\n";
        std::cout << "Path graph: facts " << facts.size() << " (+ inverses), train " << train.size()
                  << ", valid " << valid.size() << "\n";
        std::cout << "Test queries: " << eval_queries.size() << "\n";
    }

private:
    static std::string strip(const std::string& s) {
        size_t first = s.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        size_t last = s.find_last_not_of(" \t\r\n");
        return s.substr(first, last - first + 1);
    }

    static std::vector<std::string> read_lines(const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Cannot open " + path);
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in, line)) {
            line = strip(line);
            if (!line.empty()) lines.push_back(line);
        }
        return lines;
    }

    int inverse(int r) const {
        int n_base = rel2idx.size() / 2;
        return r < n_base ? r + n_base : r - n_base;
    }

    std::vector<Triple> load_triples(const std::string& path) {
        std::vector<Triple> res;
        for (const auto& line : read_lines(path)) {
            std::vector<std::string> parts;
            std::stringstream ss(line);
            std::string p;
            while (std::getline(ss, p, '\t')) parts.push_back(p);
            if (parts.size() < 3) throw std::runtime_error("Malformed line in " + path + ": " + line);
            auto e1 = ent2idx.find(parts[0]);
            auto r = rel2idx.find(parts[1]);
            auto e2 = ent2idx.find(parts[2]);
            if (e1 == ent2idx.end() || r == rel2idx.end() || e2 == ent2idx.end()) {
                throw std::runtime_error("Unknown entity or relation in " + path + ": " + line);
            }
            res.push_back({e2->second, r->second, e1->second});
        }
        return res;
    }
};
