#ifndef RIPSER_CONFIG
#define RIPSER_CONFIG

#include <algorithm>
#include <limits>
#include <fstream>
#include <vector>


struct ripser_config {
	std::string input_path;
	std::string input_type;
	std::string output_path;
	int dim_max;
	float ratio;
	float config_threshold;
	bool use_enclosing_threshold;
	bool use_union_find;
	bool print_progress;
	std::vector<std::pair<int, int>> absolute_subcomplex;
	std::vector<std::pair<int, int>> config_relative_subcomplex;
	std::vector<std::pair<int, int>> relative_subcomplex;

	ripser_config()
		: input_path(""),
		  input_type("lower_distance_matrix"),
		  output_path(""),
		  dim_max(1),
		  ratio(1.0),
		  config_threshold(-1.0),
		  use_enclosing_threshold(false),
		  use_union_find(false),
		  print_progress(false),
		  absolute_subcomplex(),
		  config_relative_subcomplex(),
		  relative_subcomplex()
  { }
};


/* **************************************************************************
 * Config Parsing
 * *************************************************************************/

std::pair<int, int> parse_interval(std::string tok) {
	size_t dash_pos = tok.find("-");
	if(dash_pos == std::string::npos) {
		return std::make_pair(std::stoi(tok), std::stoi(tok));
	} else {
		int start = -1;
		int end = -1;
		std::string string_start = tok.substr(0, dash_pos);
		if(!string_start.empty()) {
			start = std::stoi(string_start);
		}
		if(dash_pos + 1 < tok.length()) {
			std::string string_end = tok.substr(dash_pos + 1, tok.length());
			if(!string_end.empty()) {
				end = std::stoi(string_end);
			}
		}
		return std::make_pair(start, end);
	}
}

int map_relative_to_absolute(ripser_config& config, int x) {
	int abs_cummu = 0;
	for(auto& abs_intv : config.absolute_subcomplex) {
		abs_cummu += abs_intv.second - abs_intv.first + 1;
		if(abs_intv.first <= x && abs_intv.second >= x) {
			return abs_cummu - (abs_intv.second - x + 1);
		}
	}
	return x;
}

ripser_config read_config(char* configpath) {
	std::ifstream file_stream(configpath);
	if(file_stream.fail()) {
		std::cerr << "error: couldn't open config file " << configpath << std::endl;
		exit(-1);
	}
	ripser_config config;
	std::string line;
	while(getline(file_stream, line)) {
		line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
		if(line.empty() || line[0] == '#' || line[0] == ';' || line[0] == '[') {
		  continue;
		}
		size_t delim_pos = line.find("=");
		std::string name = line.substr(0, delim_pos);
		std::string string_value = line.substr(delim_pos + 1);
		if(name == "input_path") {
			config.input_path = string_value;
		}
		if(name == "output_path") {
			config.output_path = string_value;
		}
		if(name == "input_type") {
			config.input_type = string_value;
		}
		if(name == "dim_max") {
			config.dim_max = std::stoi(string_value);
		}
		if(name == "ratio") {
			config.ratio = std::stod(string_value);
		}
		if(name == "use_enclosing_threshold") {
			config.use_enclosing_threshold = (string_value == "true") || (string_value == "1");
		}
		if(name == "use_union_find") {
			config.use_union_find = (string_value == "true" || (string_value == "1"));
		}
		if(name == "threshold") {
			config.config_threshold = std::stod(string_value);
		}
		if(name == "print_progress") {
			config.print_progress = (string_value == "true") || (string_value == "1");
		}
		if(name == "absolute_subcomplex") {
			if(!string_value.empty()) {
				size_t comma_pos = 0;
				std::string tok;
				while((comma_pos = string_value.find(",")) != std::string::npos) {
					tok = string_value.substr(0, comma_pos);
					config.absolute_subcomplex.push_back(parse_interval(tok));
					string_value.erase(0, comma_pos + 1);
				}
				config.absolute_subcomplex.push_back(parse_interval(string_value));
			}
		}
		if(name == "relative_subcomplex") {
			if(!string_value.empty()) {
				size_t comma_pos = 0;
				std::string tok;
				while((comma_pos = string_value.find(",")) != std::string::npos) {
					tok = string_value.substr(0, comma_pos);
					config.config_relative_subcomplex.push_back(parse_interval(tok));
					string_value.erase(0, comma_pos + 1);
				}
				config.config_relative_subcomplex.push_back(parse_interval(string_value));
			}
		}
	}
	// Translate the relative part to account for the absolute intervals
	for(auto& rel_int : config.config_relative_subcomplex) {
		int a = map_relative_to_absolute(config, rel_int.first);
		int b = map_relative_to_absolute(config, rel_int.second);
		config.relative_subcomplex.push_back(std::make_pair(a, b));
	}
	return config;
}

#endif
