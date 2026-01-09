/*
 * cpp_utils.cpp
 *
 *  Created on: 2021-10-27
 *      Author: fenghe
 */

#include "../clib/utils.h"
#include "cpp_utils.hpp"
#include <cstring>
#include <iostream>
#include <fstream>

void split_string(std::vector<std::string> &item_value, const char * split_line, const char *split_str){
	item_value.clear();
	std::string s(split_line);
	std::string separator(split_str);
	int t = 1, pre = 0;
	for(int item_idx = 1; item_idx < 10000000; item_idx++){
        t = s.find(separator, pre);
        if (t != -1)
        	item_value.push_back(s.substr(pre, t - pre));
        else
            break;
        pre = t + separator.size();
	}
    if (pre < s.size())
    	item_value.push_back(s.substr(pre));
}

void split_string_SKIP_blank(std::vector<std::string> &item_value, const char * split_line, const char *split_str){
	item_value.clear();
	std::string s(split_line);
	std::string separator(split_str);
	int t = 1, pre = 0;
	for(int item_idx = 1; item_idx < 10000000; item_idx++){
        t = s.find(separator, pre);
        if (t != -1){
        	if(t - pre > 0)
        		item_value.push_back(s.substr(pre, t - pre));
        }
        else
            break;
        pre = t + separator.size();
	}
    if (pre < s.size())
    	item_value.push_back(s.substr(pre));
}

#define ADJUST_INDICES(start, end, len)     \
    if (end > len)                          \
        end = len;                          \
    else if (end < 0) {                     \
        end += len;                         \
        if (end < 0)                        \
        end = 0;                            \
    }                                       \
    if (start < 0) {                        \
        start += len;                       \
        if (start < 0)                      \
        start = 0;                          \
    }

int _string_tailmatch(const std::string&self, const std::string&substr, int start, int end, int direction)
{
	int selflen = (int)self.size();
	int slen = (int)substr.size();
	const char* str = self.c_str();
	const char* sub = substr.c_str();
	ADJUST_INDICES(start, end, selflen);
	if (direction < 0){
		if (start + slen>selflen)				return 0;
	}
	else{
		if (end - start<slen || start>selflen)	return 0;
		if (end - slen > start)					start = end - slen;
	}
	if (end - start >= slen)					return !memcmp(str + start, sub, slen);
	return 0;
}

bool endswith(const std::string&str, const std::string&suffix, int start, int end)
{
	int result = _string_tailmatch(str, suffix, start, end, +1);
	return static_cast<bool>(result);
}

bool startswith(const std::string&str, const std::string&suffix, int start, int end)
{
	int result = _string_tailmatch(str, suffix, start, end, -1);
	return static_cast<bool>(result);
}

#define MAX_LINE_LENGTH 10000000
void load_string_list_from_file_MAX_line(const char * string_fn, std::vector<std::string> &v, int64_t max_line_number){
	//get file names
	char *temp = new char[MAX_LINE_LENGTH];//10M
	fexist_check(string_fn);
	std::ifstream name_list_File(string_fn);
	int64_t line_number = 0;
	while(true){
		name_list_File.getline(temp, MAX_LINE_LENGTH);
		if(name_list_File.eof())	break;
		if(max_line_number == line_number++) break;
		v.emplace_back(temp);
	}
	name_list_File.close();
}

void load_string_list_from_file(const char * string_fn, std::vector<std::string> &v){
	load_string_list_from_file_MAX_line(string_fn, v, MAX_int64_t);
}

void load_int_list_from_file(char * int_fn, std::vector<int> &v){
	//get file names
	char *temp = new char[MAX_LINE_LENGTH];//10M
	fexist_check(int_fn);
	std::ifstream name_list_File(int_fn);
	while(true){
		name_list_File.getline(temp, MAX_LINE_LENGTH);
		if(*temp == 0)	break;
		v.emplace_back(atoi(temp));
	}
	name_list_File.close();
}

void vector_dump_bin(const char * path_name, const char * fn, void * data, uint64_t data_size){
	std::string full_name(path_name);
	full_name += "/";
	full_name += fn;
	FILE *f_dump = xopen (full_name.c_str(), "wb");
	fwrite(&data_size, 1, 8, f_dump);
	fwrite(data, 1, data_size, f_dump);
	fclose(f_dump);
}

//load from file, return the true load data size(in byte)
uint64_t vector_load_bin(const char * path_name, const char * fn, void ** data){
	std::string full_name(path_name);
	full_name += "/";
	full_name += fn;
	FILE *f_dump = xopen (full_name.c_str(), "rb");

	uint64_t data_size;
	xread(&data_size, 1, 8, f_dump);//equal to fread ....
	(*data) = (void* ) xmalloc (data_size);
	xread(*data, 1, data_size, f_dump);

	fclose(f_dump);
	return data_size;
}

