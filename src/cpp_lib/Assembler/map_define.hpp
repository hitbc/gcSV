/*
 * map_define.hpp
 *
 *  Created on: 2021年8月9日
 *      Author: fenghe
 */

#ifndef CPP_LIB_ASSEMBLER_MAP_DEFINE_HPP_
#define CPP_LIB_ASSEMBLER_MAP_DEFINE_HPP_

#include <set>
#include <unordered_map>
#include <string>

typedef int32_t pos_t;

typedef std::unordered_map<std::string, unsigned> str_uint_map_t;
// maps kmers to supporting reads
typedef std::unordered_map<std::string, std::set<unsigned>> str_set_uint_map_t;
typedef std::unordered_map<std::string, std::pair<unsigned, unsigned>> str_pair_uint_map_t;

#endif /* CPP_LIB_ASSEMBLER_MAP_DEFINE_HPP_ */
