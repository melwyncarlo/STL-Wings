#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#ifndef MISC_H
#define MISC_H

using namespace std;


void space();
bool isNum(string str);
string charAndUnderscore(string str);
string repeated_string(string str, int repeats);
string decimal_point(string num_str, int points);
string remove_spaces(string str);
double arithmetic_operation(double val1, double val2, string op);


/* ---------------------------------------------------------------------------------- */


void write_to_file(string directory, string file_intro, string file, string extension, string data);
bool equation_valid(string input_str);


#endif



