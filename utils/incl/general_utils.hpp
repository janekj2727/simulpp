/*
* Header file for general_utils.cpp
* Author: J. Janek
*/
#ifndef GEN_UTILS_HEAD
#define GEN_UTILS_HEAD

#include <vector>

// based on https://en.wikipedia.org/wiki/ANSI_escape_code
// enumerates for terminal colors
enum bg_cols { bg_black = 40, bg_red, bg_green, bg_yellow, bg_blue, bg_magenta, bg_cyan, bg_white, bg_gray = 100, bg_bright_red, bg_bright_green, bg_bright_yellow, bg_bright_blue, bg_bright_magenta, bg_bright_cyan, bg_bright_white };
enum fg_cols { fg_black = 30, fg_red, fg_green, fg_yellow, fg_blue, fg_magenta, fg_cyan, fg_white, fg_gray = 90, fg_bright_red, fg_bright_green, fg_bright_yellow, fg_bright_blue, fg_bright_magenta, fg_bright_cyan, fg_bright_white };

// terminal defines
#define COL_ESCAPE_IN "\033["
#define COL_ESCAPE_OUT 'm'
#define TEXTBF "\033[1m"
#define COL_RESET "\033[0m"

void print_progress_bar(int percentage);
void print_warning(int error_level, std::string str1, std::string str2 = "", std::string str3 = "", std::string str4 = "");
int basic_statistics(std::ofstream &fprt, char *cpa_name, std::vector<std::string> qnames);

#endif