/*
 * get_option_cpp.hpp
 *
 *  Created on: 2021-4-12
 *      Author: fenghe
 */

#ifndef GET_OPTION_CPP_HPP_
#define GET_OPTION_CPP_HPP_

#include <string>
#include<string.h>
#include <vector>
#include <map>
#include <getopt.h>
extern "C"{
#include "../clib/desc.h"
#include "../clib/utils.h"
}

struct ARG_TYPE{
	static const int BLANK = 0;
	static const int INT = 1;
	static const int STR  = 2;
	static const int DOUBLE = 3;
};

struct option_item{

	// for blank type
	option_item(const char * long_option_, char short_option_, const char * help_msg_){
		is_required_argument = false;
		new_item(long_option_, short_option_, help_msg_);
	}
	// for int type
	option_item(const char * long_option_, char short_option_, const char * help_msg_, bool has_default_, int default_value_){
		is_required_argument = true; arg_type = ARG_TYPE::INT;	has_default_arg = has_default_;
		new_item(long_option_, short_option_,  help_msg_);
		if(has_default_arg)	default_value_int = default_value_;
	}
	//for type string
	option_item(const char * long_option_, char short_option_, const char * help_msg_, bool has_default_, const char * default_value_){
		is_required_argument = true; arg_type = ARG_TYPE::STR;	has_default_arg = has_default_;
		new_item(long_option_, short_option_, help_msg_);
		if(has_default_arg)		default_value_str.append(default_value_);
	}
	//for double
	option_item(const char * long_option_, char short_option_, const char * help_msg_, bool has_default_, double default_value_){
		is_required_argument = true; arg_type = ARG_TYPE::DOUBLE;	has_default_arg = has_default_;
		new_item(long_option_, short_option_, help_msg_);
		if(has_default_arg)	default_value_double = default_value_;
	}

	// adding long/short argument and is_required_argument_
	void new_item(const char * long_option_, char short_option_, const char * help_msg_){
		long_option.clear();
		long_option.append(long_option_);

		short_option = short_option_;

		help_msg.clear();
		other_help_msg.clear();
		for(uint32_t i = 0; i < 5; i++)
			help_msg.append(" ");

		help_msg.append(help_msg_);
		//help_msg.append("\n");
	}

	void add_help_msg(const char * help_msg_){
		other_help_msg.append("");
		for(uint32_t i = 0; i < 5 ; i++)
			other_help_msg.append(" ");
		other_help_msg.append("");
		other_help_msg.append(help_msg_);
		other_help_msg.append("\n");
	}

	std::string long_option;
	char short_option;

	//argument
	bool is_required_argument = false;
	int arg_type = ARG_TYPE::BLANK;

	//default argument
	bool has_default_arg = true;
	int default_value_int = 0;
	std::string default_value_str;
	double default_value_double = 0;

	//argument pointer
	bool has_arg_pointer = false;
	const void * arg_p = NULL;

	void set_arg_pointer(void * p){
		if(has_arg_pointer){ show_option(stderr); xassert(0, "Conflict reset argument pointer"); }
		has_arg_pointer = true;
		arg_p = p;
	}

	void set_default_option_value(){
		if(has_arg_pointer){
			if(arg_type == ARG_TYPE::BLANK){ bool * c_arg = (bool *)arg_p;	c_arg[0] = false;	}
			else if(has_default_arg){
					 if(arg_type == ARG_TYPE::INT){ int * c_arg = (int *)arg_p;	c_arg[0] = default_value_int;	}
				else if(arg_type == ARG_TYPE::STR){ const char**c_arg = (const char **)arg_p; *c_arg = ( char*)xmalloc(default_value_str.size() + 1);	strcpy((char *)(*c_arg), default_value_str.c_str());}
				else if(arg_type == ARG_TYPE::DOUBLE){ double * c_arg = (double *)arg_p;	c_arg[0] = default_value_double;	}
			}
		}
	}

	bool set_option_value(char * option_char){
		if(has_arg_pointer){
				 if(arg_type == ARG_TYPE::BLANK){ bool * c_arg = (bool *)arg_p;	c_arg[0] = true;	}
			else if(arg_type == ARG_TYPE::INT){ int * c_arg = (int *)arg_p;	c_arg[0] = atoi(option_char);	}
			else if(arg_type == ARG_TYPE::STR){ const char**c_arg = (const char **)arg_p; *c_arg = option_char;	}
			else if(arg_type == ARG_TYPE::DOUBLE){ double * c_arg = (double *)arg_p; c_arg[0] = atof(option_char);	}
			return true;
		}
		return false;
	}

	//help message
	std::string help_msg;
	std::string other_help_msg;

	void show_option(FILE * output){
		fprintf(output, "%c %s %s\n", short_option, long_option.c_str(), help_msg.c_str());
	}

	void show_current_value(FILE * output){
		fprintf(output, "\t\t\t%s", long_option.c_str());
		if(has_arg_pointer){
			if(arg_type == ARG_TYPE::BLANK){ fprintf(output, "=%s", (((bool *)arg_p)[0])==true?"true":"false");	}
			else if(arg_type == ARG_TYPE::INT){ fprintf(output, "=%d", ((int *)arg_p)[0]); }
			else if(arg_type == ARG_TYPE::STR){ fprintf(output, "=%s", ((const char **)arg_p)[0]); }
			else if(arg_type == ARG_TYPE::DOUBLE){ fprintf(output, "=%f", ((double *)arg_p)[0]); }
		}
		fprintf(output, "\n");
	}
};

struct options_list{

	std::vector<option_item> l;
	std::string short_option_str;
	std::vector<option> long_option;
	std::string title_str;
	std::map<char, int> conflict_map;

	bool already_parse = false;

	void show_command(FILE * output, int argc, char *argv[]){
		fprintf(output, "command=");
		for(int i = 0; i < argc; i++){
			fprintf(output, "%s ",  argv[i]);
		}
		fprintf(output, "\n");
	}

	void add_title_string(const char * title_str_){
		title_str.append(title_str_);
	}

	//adding an option; overload functions
	void add_option(const char * long_option_, char short_option_, const char * help_msg_)
	{l.emplace_back(long_option_, short_option_, help_msg_); already_parse = false;}// for simple option
	void add_option(const char * long_option_, char short_option_, const char * help_msg_, bool has_default_, int    default_value_)
	{l.emplace_back(long_option_, short_option_, help_msg_, has_default_, default_value_);already_parse = false;}// for int type option
	void add_option(const char * long_option_, char short_option_, const char * help_msg_, bool has_default_, const char*  default_value_)
	{l.emplace_back(long_option_, short_option_, help_msg_, has_default_, default_value_);already_parse = false;}// for string type option
	void add_option(const char * long_option_, char short_option_, const char * help_msg_, bool has_default_, double default_value_)
	{l.emplace_back(long_option_, short_option_, help_msg_, has_default_, default_value_);already_parse = false;}// for int type option

	void set_arg_pointer_back(void * p){
		l.back().set_arg_pointer(p);
	}

	void add_help_msg_back(const char * help_msg_){
		l.back().add_help_msg(help_msg_);
	}

	//used when all options are simple registration option
	int default_option_handler(int argc, char *argv[]){
		char c;
		const char *short_option = ""; option *long_option = NULL;
		get_short_long_option_string(short_option, long_option);
		while((c = getopt_long(argc + 1 , argv - 1, short_option, long_option, NULL)) != -1)
		{
			if(false == set_option(c, optarg)){
				switch(c){
				default: fprintf(stderr, "[main:] Unknown option: %c\n", c); return output_usage(); break;
				}
			}
			if(c == 'h') return 1;//return for help
		}
		return 0;
	}

	void parse(){
		add_option("help", 			'h', "show this message");
		short_option_str.clear();
		long_option.clear();
		conflict_map.clear();
		for(uint32_t i = 0; i < l.size(); i++){
			//conflict check
			auto conflict_map_it = conflict_map.find(l[i].short_option);
			if(conflict_map_it!=conflict_map.end()){
				fprintf(stderr, "FATAL ERROR: Short option conflict:\n");
				l[conflict_map_it->second].show_option(stderr);
				fprintf(stderr, "with:\n");
				l[i].show_option(stderr);
				xassert(0, "");
			}
			else
				conflict_map[l[i].short_option] = i;
			//add to short_option_str
			short_option_str += l[i].short_option;
			if(l[i].is_required_argument)
				short_option_str += ':';
			//add to long option
			long_option.emplace_back();
			auto &c_l = long_option.back();
			c_l.name = l[i].long_option.c_str();
			c_l.has_arg = (l[i].is_required_argument)?required_argument:no_argument;
			c_l.flag = NULL;
			c_l.val = (int)l[i].short_option;
			//set default value
			l[i].set_default_option_value();
		}
	}

	void get_short_long_option_string(const char * &short_op, option * &long_op){
		xassert(already_parse == false, "Don`t add new option after option parsing"); parse();
		short_op = short_option_str.c_str(); long_op = &(long_option[0]);
		already_parse = true;
	}

	bool set_option(int short_option, char * option_char){
		auto conflict_map_it = conflict_map.find(short_option);
		if(conflict_map_it == conflict_map.end())
			return false;
		int idx = conflict_map_it->second;
		if(idx == (int)l.size() - 1){//help
			output_usage();
			return true;
		}
		return l[idx].set_option_value(option_char);
	}

	int output_usage(){
		fprintf(stderr, "\n");
		fprintf(stderr, "Program:   %s\n", PACKAGE_NAME);
		fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION);
		fprintf(stderr, "Contact:   %s\n\n", CONTACT);
		fprintf(stderr, "%s", title_str.c_str());
		//output options
		fprintf(stderr, "\n");
		for(uint32_t i = 0; i < l.size(); i++){
			fprintf(stderr, "\t\t\t\t-%c --%s", l[i].short_option, l[i].long_option.c_str());
			if(l[i].long_option.size() < 15){
				uint32_t long_opt_size_left = 15 - l[i].long_option.size();
				for(uint32_t i = 0; i < long_opt_size_left ; i++)
					fprintf(stderr, " ");
			}
			if(l[i].is_required_argument)
			{
				switch(l[i].arg_type){
					case ARG_TYPE::BLANK :fprintf(stderr, "\t[bool]"); break;
					case ARG_TYPE::INT :fprintf(stderr, "\t[INT]"); break;
					case ARG_TYPE::STR :fprintf(stderr, "\t[STR]"); break;
					case ARG_TYPE::DOUBLE :fprintf(stderr, "\t[DOUBLE]"); break;
					default :fprintf(stderr, "\t\t"); break;
				}
			}else
				fprintf(stderr, "\t     ");
			fprintf(stderr, "%s ", l[i].help_msg.c_str());
			if(l[i].has_default_arg){
				switch(l[i].arg_type){
					case ARG_TYPE::BLANK :fprintf(stderr, "\t"); break;
					case ARG_TYPE::INT :fprintf(stderr, "[%d]", l[i].default_value_int); break;
					case ARG_TYPE::STR :fprintf(stderr, "[%s]", l[i].default_value_str.c_str()); break;
					case ARG_TYPE::DOUBLE :fprintf(stderr, "[%f]", l[i].default_value_double); break;
					default :fprintf(stderr, "\t"); break;
				}
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "%s", l[i].other_help_msg.c_str());
		}
		return 1;
	}

	void show_c_value(FILE * output){
		fprintf(output, "\n\t\tCurrent parameters:\n");
		for(auto & i : l){
			i.show_current_value(output);
		}
	}
};

struct SECOND_LEVEL_FUNC{

	SECOND_LEVEL_FUNC(const char * command_, const char * help_msg_, int (*run_method_) (int argc, char *argv[])){
		command.clear();
		command.append(command_);
		help_msg.clear();

		if(command.size() < 20){
			for(uint32_t i = 0; i < 20 - command.size(); i++)
				help_msg.append(" ");
		}

		help_msg.append(help_msg_);
		help_msg.append("\n");
		run_method = run_method_;
	}

	void add_help_msg(const char * help_msg_){
		for(uint32_t i = 0; i < 28 ; i++)
			help_msg.append(" ");
		help_msg.append("- ");
		help_msg.append(help_msg_);
		help_msg.append("\n");
	}

	std::string command;
	std::string help_msg;
	int (*run_method) (int argc, char *argv[]);
};

struct COMMAND_HANDLER{
	std::vector<SECOND_LEVEL_FUNC> f_l;//function list
	std::string main_command_name;

	void add_function(const char * command_, const char * help_msg_, int (*run_method_) (int argc, char *argv[])){
		f_l.emplace_back(command_, help_msg_, run_method_);
	}
	void add_help_msg_back(const char * help_msg_){
		f_l.back().add_help_msg(help_msg_);
	}

	void set_main_command(const char * main_command_name_){
		main_command_name.append(main_command_name_);
	}

	int usage()
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "Program:	%s\n", PACKAGE_NAME);
		fprintf(stderr, "Version:	%s\n", PACKAGE_VERSION);
		fprintf(stderr, "Contact:	%s\n\n", CONTACT);

		fprintf(stderr, "Usage:		%s %s <command> [options]\n\n", PACKAGE_NAME, main_command_name.c_str());
		fprintf(stderr, "Command list: \n");
		for(auto & f: f_l)
			fprintf(stderr, "        %s%s", f.command.c_str(), f.help_msg.c_str());
		fprintf(stderr, "        --help	            show this message\n");

		return 0;
	}

	int run(int argc, char *argv[]){
		if (argc < 2)	return usage();
		bool successful_run = false;
		int return_int = 0;
		for(auto & f: f_l){ if(strcmp(argv[1], f.command.c_str()) == 0)	{ return_int = f.run_method(argc - 1, argv + 1); successful_run = true;}}
		if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0){
			return usage();
		}
		if(!successful_run){
			fprintf(stderr, "[Waring!!!] wrong command: '%s'\n", argv[1]);
			return usage();
		}
		if(return_int != 0){
			fprintf(stderr, "[Waring] Run error command: '%s'\n", argv[1]);
			return 1;
		}
		fprintf(stderr, "Successful run command: '%s'\n", argv[0]);
		return 0;

	}

};


#endif /* GET_OPTION_CPP_HPP_ */
