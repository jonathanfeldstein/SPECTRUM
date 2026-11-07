#include "exceptions.h"

ArgumentsException::ArgumentsException(const std::string &message) noexcept {
    this->message_ = message + " Usage: (string)<path_to_data> (string)<path_to_info_file>"
                               " (string)<path_to_output> (double)<epsilon> (double)<alpha>"
                               " [Optional](bool)<multiprocessing> [Optional](bool)<verbose>"
                               " [Optional](bool)<safe>"
                               " [Optional](double)<max_lambda2> [Optional](int)<min_cluster_size>";
}

FileNotFoundException::FileNotFoundException(const std::string &filename) noexcept : filename(filename){
    this->message_ = "File not found: " + filename;
}

KeyNotFoundException::KeyNotFoundException(const std::string &key) noexcept : key(key){
    this->message_ = "Key not found: " + key;
}

FileNotOpenedException::FileNotOpenedException(const std::string &filename) noexcept {
    this->message_ = "File did not open correctly: " + filename;
}

DatabaseSyntaxException::DatabaseSyntaxException(const std::string &line) noexcept {
    this->message_ = "Database Syntax incorrect! Line provided was:\n\n"
                     + line + "\n\n" +
                     "Expected: predicateName(argument1, argument2, ...):weight \nWhere:\n"
                     "predicateName is an alpha-numeric string, , only special character supported is '_'. \n"
                     "Followed by a pair of parentheses, containing the arguments (constants, e.g. Alice, ToyStory1), \n"
                     "which are alpha-numeric strings without special characters, separated by commas. \n"
                     "Make sure no trailing comma exists in the parentheses.\n"
                     "This is followed by an optional colon (:). This is only needed if a weight is provided, otherwise it is ignored. \n"
                     "If you provide a weight, you should provide a weight between 0.0 and 1.0, both 0 and 1 are also accepted.\n"
                     "Careful: white spaces will be ignored! ";
}

InfoSyntaxException::InfoSyntaxException(const std::string &line) noexcept {
    this->message_ = "Info File Syntax incorrect! Line provided was:\n\n"
                     + line + "\n\n" +
                     "Expected: predicateName(argumentType1, argumentType2, ...) \nWhere:\n"
                     "predicateName is an alpha-numeric string, only special character supported is '_'. \n"
                     "Followed by a pair of parentheses, containing the argument types (e.g. Person, Movie), \n"
                     "which are alpha-numeric strings without special characters, separated by commas. \n"
                     "Make sure no trailing comma exists in the parentheses.\n"
                     "Careful: white spaces will be ignored! ";
}

