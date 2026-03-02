#ifndef PARSE_LITERALS_H
#define PARSE_LITERALS_H

#include <vector>
#include <string>

int literals_find_datatype(int n, const std::vector<std::string> &literals);
int literal_get_datatype(const std::string &literal);
int literal_to_int(const std::string &literal);
double literal_to_double(const std::string &literal);

#endif
