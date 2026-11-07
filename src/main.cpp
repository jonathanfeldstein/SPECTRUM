#include "Utils/profiling_utils.h"
#include "CandidateRules/CandidateRules.h"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"


using namespace Eigen;

using namespace std;

int main(int argc, char** argv)
{
    boost::property_tree::ptree config;
    boost::property_tree::read_json(argv[1], config);
    int M = stoi(argv[2]);
    double epsilon = stod(argv[3]);
    int rule_length = stoi(argv[4]);
    size_t top_rules = stoi(argv[5]);
    double min_fidelity = argc > 6 ? stod(argv[6]) : 1.0;
    auto logic = config.get<string>("logic");
    auto data_folder = config.get<string>("data_folder");
    auto info_file = config.get<string>("info_file");
    vector<string> data_files;
    for (const auto& item : config.get_child("data_files")) {
        data_files.push_back(item.second.get_value<string>());
    }
    vector<string> target_predicates{};
    for (const auto& item : config.get_child("targets")) {
        target_predicates.push_back(item.second.get_value<string>());
    }
    auto output_folder = config.get<string>("output_folder");
    auto debug = config.get<bool>("debug");
    auto verbose = config.get<bool>("verbose");
    auto force_cyclic = config.get<bool>("force_cyclic");

    // Load the entire database as a single hypergraph
    for(const auto& data_file : data_files){
        if(debug) {
            cout << "DEBUG MODE - Running on a single file: " << data_file << "\n";
        }else{
            cout << "Running on file: " << data_file << "\n";
        }
        Timer time_optics("SPECTRUM");
        Timer data_loading("Data Loading");
        string db_file_path = data_folder + data_file + ".db";
        HyperGraph complete_hg(db_file_path, info_file, false);
        data_loading.Stop();
        Timer fragment_timer("Fragment Creation");
        int max_number_of_random_walks = ceil((M*log(M))/(pow(epsilon,2)*(int)complete_hg.number_of_nodes()));
        if(verbose){
            cout << "Max number of random walks: " << max_number_of_random_walks << "\n";
        }
        bool has_unaries = complete_hg.has_unaries();
#ifdef USE_RANDOM_WALKS
        FragmentCollection fragments(complete_hg, max_number_of_random_walks, rule_length, force_cyclic);
#else
        FragmentCollection fragments(complete_hg, max_number_of_random_walks, rule_length, force_cyclic, has_unaries);
#endif
        if(verbose){
            fragment_timer.Stop();
        }
        Timer candidate_time("Candidate Creation");
        CandidateRuleCollection candidates(std::move(fragments), min_fidelity, target_predicates, has_unaries);
        if(verbose){
            candidate_time.Stop();
        }
        Timer candidate_optimization("Optimizing Candidates");
        string output_file_name = output_folder +  logic + "/" + data_file;
        candidates.print_rules(logic,
                               output_file_name,
                               M,
                               target_predicates,
                                top_rules);
        if(verbose){
            candidate_optimization.Stop();
        }
        time_optics.Stop();
        if(debug){
            break;
        }
    }

    return 0;
}
