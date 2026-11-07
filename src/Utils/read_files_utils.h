#ifndef SPECTRUM_READ_FILES_UTILS_H
#define SPECTRUM_READ_FILES_UTILS_H

#include <string>
#include <algorithm>
#include <sys/stat.h>
#include <regex>
#include "exceptions.h"

#include "../HyperGraphs/Relation.h"

GroundRelation parse_line_db(std::string line, bool safe);
Relation parse_line_info(std::string line, bool safe);
bool file_exists(const std::string& path_to_file);

#endif //SPECTRUM_READ_FILES_UTILS_H
