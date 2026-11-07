#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#include <cstdio>
#include <algorithm>
#include <sstream>
#include <stdexcept>

struct Rule {
    double confidence;
    std::string head_rel;
    std::vector<std::string> body_rels;
};

class RuleParser {
private:

    struct Atom {
        std::string pred;
        std::string arg1;
        std::string arg2;
    };

    static std::string trim(const std::string& s)
    {
        size_t first = s.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            return "";

        size_t last = s.find_last_not_of(" \t\r\n");

        return s.substr(first, last - first + 1);
    }

    static Atom parse_atom(const std::string& atom)
    {
        std::string s = trim(atom);

        size_t l = s.find('(');
        size_t c = s.find(',', l);
        size_t r = s.find(')', c);

        if (l == std::string::npos ||
            c == std::string::npos ||
            r == std::string::npos)
        {
            throw std::runtime_error(
                "Cannot parse atom: " + s);
        }

        Atom a;

        a.pred = trim(s.substr(0, l));
        a.arg1 = trim(s.substr(l + 1, c - l - 1));
        a.arg2 = trim(s.substr(c + 1, r - c - 1));

        return a;
    }

    // NCRL's translator writes the weight as "%.3f" and kg_completion.py
    // reads back only its first 5 characters (float(conf[0:5])).
    static double ncrl_confidence(double weight)
    {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.3f", weight);
        return std::stod(std::string(buf).substr(0, 5));
    }

    static std::vector<std::pair<std::string,std::string>>
    canonical_patterns(
        const std::string& X,
        const std::string& Y,
        size_t body_len)
    {
        std::vector<std::pair<std::string,std::string>> patterns;

        //
        // Length 1:
        //
        //    body(Y,X)
        //

        if (body_len == 1)
        {
            patterns.push_back({Y, X});
            return patterns;
        }

        //
        // Create
        //
        // __Z1 ... __Zn
        //

        std::vector<std::string> z;

        for (size_t i = 1; i < body_len; i++)
        {
            z.push_back("__Z" + std::to_string(i));
        }

        //
        // Reverse so nearest to Y comes first
        //

        std::reverse(z.begin(), z.end());

        //
        // (__Z_last, Y)
        //

        patterns.push_back({z[0], Y});

        //
        // (__Zi,__Zi-1)
        //

        for (size_t i = 0; i + 1 < z.size(); i++)
        {
            patterns.push_back(
            {
                z[i + 1],
                z[i]
            });
        }

        //
        // (X,__Z1)
        //

        patterns.push_back(
        {
            X,
            z.back()
        });

        return patterns;
    }

public:

    // rules dropped during parsing: not a simple chain / relation not evaluable
    static inline int skipped_non_chain = 0;
    static inline int skipped_relation = 0;

    static std::string unary_to_binary(std::string s)
    {
        std::regex re(
            "([a-zA-Z0-9_&]+)\\(([A-Z])\\)(?![,\\)])");

        return std::regex_replace(
            s,
            re,
            "$1($2,$2)");
    }

    static std::string reconstruct_predicates(std::string s)
    {
        std::regex re1(
            "([^\\(]+)\\(([A-Z]),([A-Z]),'(.*?)'\\)");

        s = std::regex_replace(
            s,
            re1,
            "$1$4($2,$3)");

        std::regex re2(
            "([^\\(]+)\\(([A-Z]),'(.*?)'\\)");

        s = std::regex_replace(
            s,
            re2,
            "$1$3($2)");

        return s;
    }

    static std::vector<Rule> parse(
        const std::string& json_path,
        int top_k,
        const std::unordered_set<std::string>& train_relations){
    std::vector<Rule> parsed_rules;

    std::ifstream in(json_path);

    if (!in)
        return parsed_rules;

    std::string json_str(
        (std::istreambuf_iterator<char>(in)),
        std::istreambuf_iterator<char>());

    size_t pos = 0;

    while ((pos = json_str.find("\"", pos)) != std::string::npos)
    {
        size_t end_pos = json_str.find("\"", pos + 1);

        if (end_pos == std::string::npos)
            break;

        std::string r =
            json_str.substr(
                pos + 1,
                end_pos - pos - 1);

        pos = end_pos + 1;

        if (r.find(" >> ") == std::string::npos)
            continue;

        if (r.size() > 2 &&
            r.substr(r.size() - 2) == "^2")
        {
            r = r.substr(0, r.size() - 2);
        }

        size_t colon_pos = r.find(": ");

        if (colon_pos == std::string::npos)
            continue;

        double confidence =
            ncrl_confidence(
                std::stod(r.substr(0, colon_pos)));

        std::string body_head =
            r.substr(colon_pos + 2);

        size_t arrow_pos =
            body_head.find(" >> ");

        if (arrow_pos == std::string::npos)
            continue;

        std::string body_str =
            body_head.substr(0, arrow_pos);

        std::string head_str =
            body_head.substr(arrow_pos + 4);

        body_str =
            reconstruct_predicates(
                unary_to_binary(body_str));

        head_str =
            reconstruct_predicates(
                unary_to_binary(head_str));

        Atom head = parse_atom(head_str);

        std::vector<Atom> body_atoms;

        std::stringstream ss(body_str);

        std::string atom;

        // Atoms are separated by " & "; relation names may themselves contain '&'
        // (UMLS: interacts&with), so a bare '&' must not split.
        size_t start = 0;
        size_t sep;
        while ((sep = body_str.find(" & ", start)) != std::string::npos)
        {
            body_atoms.push_back(
                parse_atom(body_str.substr(start, sep - start)));
            start = sep + 3;
        }
        body_atoms.push_back(
            parse_atom(body_str.substr(start)));

        // Rules NCRL's translator (exp/run_ncrl_eval.py) never passes on:
        // relations absent from train.txt, and relation names containing '&'.
        bool keep = true;
        std::vector<std::string> preds{head.pred};
        for (const auto& a : body_atoms)
            preds.push_back(a.pred);
        for (const auto& p : preds)
        {
            if (!train_relations.count(p) ||
                p.find('&') != std::string::npos)
            {
                keep = false;
            }
        }
        if (!keep)
        {
            skipped_relation++;
            continue;
        }

        //
        // Collect variables
        //

        std::vector<std::string> vars;

        auto add_var =
        [&](const std::string& v)
        {
            if (std::find(
                    vars.begin(),
                    vars.end(),
                    v) == vars.end())
            {
                vars.push_back(v);
            }
        };

        add_var(head.arg1);
        add_var(head.arg2);

        for (const auto& a : body_atoms)
        {
            add_var(a.arg1);
            add_var(a.arg2);
        }

        std::vector<std::string> internal;

        for (const auto& v : vars)
        {
            if (v != head.arg1 &&
                v != head.arg2)
            {
                internal.push_back(v);
            }
        }

        //
        // Rule must be a simple chain
        //

        if (internal.size() != body_atoms.size() - 1)
        {
            skipped_non_chain++;

            continue;
        }

        //
        // Build __Zi mapping
        //

        std::vector<std::string> z;

        for (size_t i = 1;
             i < body_atoms.size();
             i++)
        {
            z.push_back(
                "__Z" +
                std::to_string(i));
        }

        std::reverse(
            z.begin(),
            z.end());

        std::unordered_map<
            std::string,
            std::string> mapping;

        mapping[head.arg1] =
            head.arg1;

        mapping[head.arg2] =
            head.arg2;

        for (size_t i = 0;
             i < internal.size();
             i++)
        {
            mapping[z[i]] =
                internal[i];
        }

        auto patterns =
            canonical_patterns(
                head.arg1,
                head.arg2,
                body_atoms.size());

        std::vector<Atom> remaining =
            body_atoms;

        Rule rule;

        rule.confidence =
            confidence;

        rule.head_rel =
            head.pred;

        bool success = true;

        //
        // Match each expected edge
        //

        for (const auto& expected :
             patterns)
        {
            std::string from =
                mapping[expected.first];

            std::string to =
                mapping[expected.second];

            bool found = false;

            for (size_t i = 0;
                 i < remaining.size();
                 i++)
            {
                const Atom& a =
                    remaining[i];

                //
                // Forward edge
                //

                if (a.arg1 == from &&
                    a.arg2 == to)
                {
                    rule.body_rels.push_back(
                        a.pred);

                    remaining.erase(
                        remaining.begin() + i);

                    found = true;
                    break;
                }

                //
                // Reverse edge
                //

                if (a.arg1 == to &&
                    a.arg2 == from)
                {
                    rule.body_rels.push_back(
                        "inv_" + a.pred);

                    remaining.erase(
                        remaining.begin() + i);

                    found = true;
                    break;
                }
            }

            if (!found)
            {
                success = false;

                skipped_non_chain++;

                break;
            }
        }

        if (success)
        {
            parsed_rules.push_back(rule);
        }
    }

    if (top_k > 0)
    {
        std::unordered_map<
            std::string,
            int> head_counts;

        std::vector<Rule> filtered;

        for (const auto& r : parsed_rules)
        {
            if (head_counts[r.head_rel] < top_k)
            {
                filtered.push_back(r);
                head_counts[r.head_rel]++;
            }
        }

        return filtered;
    }

    return parsed_rules;
}
};