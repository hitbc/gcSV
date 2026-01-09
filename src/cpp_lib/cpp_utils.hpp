/*
 * cpp_utils.hpp
 *
 *  Created on: 2021-10-27
 *      Author: fenghe
 */

#ifndef CPP_LIB_CPP_UTILS_HPP_
#define CPP_LIB_CPP_UTILS_HPP_

#include "string"
#include "vector"

void split_string(std::vector<std::string> &item_value, const char * split_line, const char *split_str);
void split_string_SKIP_blank(std::vector<std::string> &item_value, const char * split_line, const char *split_str);

bool endswith(const std::string&str, const std::string&suffix, int start, int end);
bool startswith(const std::string&str, const std::string&suffix, int start, int end);

void load_string_list_from_file(const char * string_fn, std::vector<std::string> &v);
void load_string_list_from_file_MAX_line(const char * string_fn, std::vector<std::string> &v, int64_t max_line_number);
void load_int_list_from_file(char * int_fn, std::vector<int> &v);

//dump to file, return the true load data size(in byte)
void vector_dump_bin(const char * path_name, const char * fn, void * data, uint64_t data_size);
//load from file, return the true load data size(in byte)
uint64_t vector_load_bin(const char * path_name, const char * fn, void ** data);

#endif /* CPP_LIB_CPP_UTILS_HPP_ */
