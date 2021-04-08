#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "MISC.H"

using namespace std;


void space()
{
	cout << endl << endl;
}


bool isNum(string str)
{
	int i = 0;
	int i2 = str.length();
	bool dot_found = false;
	for (i=0; i<i2; i++)
	{
		if ((str[i] == '-'))
		{
			if (i != 0)
			{
				return false;
			}
		}
		else if ((str[i] == '.'))
		{
			if (dot_found)
			{
				return false;
			}
			else
			{
				dot_found = true;
			}
		}
		else if ((!isdigit(str[i])))
		{
			return false;
		}
	}
	return true;
}


string charAndUnderscore(string str)
{
	string temp = "";
	for (int i=0; i<(str.length()); i++)
	{
		if ((isalpha(str[i])))
		{
			temp += str.substr(i, 1);
		}
		else if ((isdigit(str[i])))
		{
			if ((str[i] != '-') && (str[i] != '.'))
			{
				temp += str.substr(i, 1);
			}
			else
			{
				temp += "_";
			}
		}
		else
		{
			temp += "_";
		}
	}
	return temp;
}


string repeated_string(string str, int repeats)
{
	string new_str = "";
	if (repeats != 0)
	{
		for (int i=1; i<=repeats; i++)
		{
			new_str += str;
		}
	}
	return new_str;
}


string decimal_point(string num_str, int points)
{
	string new_str;
	if (isNum(num_str))
	{
		string temp_str;
		string find_char = ".";
		size_t found = num_str.find(find_char);
		if (found != string::npos)
		{
			temp_str = num_str.substr(found + 1);
		}
		else
		{
			temp_str = "";
		}
		if (points == 0)
		{
		    new_str = num_str;
		}
		else if (temp_str.length() >= points)
		{
		    if (temp_str.length() != points)
		    {
    			if (stoi(temp_str.substr(points, 1)) >= 5)
    			{
    			    string temp_str2 = num_str.substr(0, found);
    			    string temp_str3 = temp_str.substr(0, points);
    			    int temp4 = stoi(temp_str3) + 1;
    			    int temp5;
    			    if (to_string(temp4).length() != temp_str3.length())
    			    {
    			        temp5 = stoi(temp_str2) + 1;
    			        temp_str2 = to_string(temp5);
    			        temp_str3 = to_string(temp4).substr(1, to_string(temp4).length()-1);
    			    }
    			    else
    			    {
    			        temp_str3 = to_string(temp4).substr(0);
    			    }
    			    new_str = temp_str2 + "." + temp_str3;
    			}
    			else
    			{
    			    new_str = num_str.substr(0, found) + "." + temp_str.substr(0, points);
    			}
		    }
		    else
		    {
		        new_str = num_str.substr(0, found) + "." + temp_str.substr(0, points);
		    }
		}
		else
		{
			new_str = num_str.substr(0, found) + "." + temp_str + repeated_string("0", points - temp_str.length());
		}
	}
	else
	{
		new_str = "0";
	}
	
	return new_str;
}


string remove_spaces(string str)
{
	string temp = "";
	for (int i=0; i<(str.length()); i++)
	{
		if (str[i] != ' ')
		{
			temp += str[i];
		}
	}

	return temp;
}


double arithmetic_operation(double val1, double val2, string op)
{
	if (op.compare("+") == 0)
	{
		return (val1 + val2);
	}
	else if (op.compare("-") == 0)
	{
		return (val1 - val2);
	}
	else if (op.compare("*") == 0)
	{
		return (val1 * val2);
	}
	else if (op.compare("/") == 0)
	{
		return (val1 / val2);
	}
	else if (op.compare("^") == 0)
	{
		return (pow(val1, val2));
	}
	else
	{
		return 0;
	}
}



/* ---------------------------------------------------------------------------------- */


void write_to_file(string directory, string file_intro, string file, string extension, string data)
{
	if (file_intro.compare("") != 0)
	{
		file_intro = file_intro + "_";
	}
	ofstream filestream;
	filestream.open (directory + "/" + file_intro + file + "." + extension);
	filestream << data;
	filestream.close();
}


bool equation_valid(string input_str)
{
	string str = remove_spaces(input_str);
	int log_open = 0;
	int nest_open = 0;
	bool nest_opened_recently = false;
	bool keyword = false;
	bool operator_used_recently = false;
	int temp_var_01 = 0;
	int temp_var_02 = 0;
	for (int j=0; j<(str.length()); j++)
	{
		if (str[j] == '(')
		{
			temp_var_01++;
		}
		else if (str[j] == ')')
		{
			temp_var_02++;
		}
	}
	if ((str[0] != '(') || (str[str.length() - 1] != ')'))
	{
		return false;
	}
	else if (temp_var_01 != temp_var_02)
	{
		return false;
	}
	else if ((str[0] == '(') && (str[1] == ')'))
	{
	    return false;
	}
	else
	{
		nest_opened_recently = true;
		for (int i=1; i<(str.length() - 1); i++)
		{
			string temp_str(1, str[i]);
			if (str[i] == '(')
			{
				if (!operator_used_recently)
				{
					if (!nest_opened_recently)
					{
						return false;
					}
				}
				nest_open++;
				nest_opened_recently = true;
				operator_used_recently = false;
				keyword = false;
			}
			else if (str[i] == ')')
			{
				nest_open--;
				if (nest_open < 0)
				{
					return false;
				}
				else if (operator_used_recently)
				{
					return false;
				}
				nest_opened_recently = false;
				operator_used_recently = false;
				keyword = false;
			}
			else if (str[i] == '>')
			{
				if ((log_open == 0) || (str[i+1] == '('))
				{
					return false;
				}
				else if ((str[i+1] != '+') || (str[i+1] != '-') || (str[i+1] != '*') || (str[i+1] != '/') || (str[i+1] != '^'))
				{
					return false;
				}
				else
				{
					log_open--;
				}
				nest_opened_recently = false;
			}
			else
			{
				if ((str[i] == '+') || (str[i] == '-') || (str[i] == '*') || (str[i] == '/') || (str[i] == '^'))
				{
					operator_used_recently = true;
				}
				else if ((str[i] == 'e') || (str[i] == 'x'))
				{
					if ((str[i-1] == '(') || (operator_used_recently))
					{
						if ((str[i+1] != ')') && (str[i+1] == '+') && (str[i+1] == '-') && (str[i+1] == '*') && (str[i+1] == '/') && (str[i+1] == '^'))
						{
							return false;
						}
					}
					else
					{
						return false;
					}
					operator_used_recently = false;
				}
				else if (str[i] == 'p')
				{
					if (str[i+1] != 'i')
					{
						return false;
					}
					i += 1;
					operator_used_recently = false;
				}
				else if (str[i] == 'l')
				{
					if ((str[i+1] == 'n'))
					{
						if ((str[i+2] == '<'))
						{
							if ((str[i+3] == '+') || (str[i+3] == '-') || (str[i+3] == '*') || (str[i+3] == '/') || (str[i+3] == '^') || (str[i+3] == ')'))
							{
								return false;
							}
							else
							{
								log_open++;
							    return true;
							}
						}
						else
						{
							return false;
						}
					}
					else if ((str[i+1] == 'o'))
					{
						if ((str[i+2] == 'g'))
						{
							if ((str[i+3] == '<'))
							{
								if ((str[i+4] == '+') || (str[i+4] == '-') || (str[i+4] == '*') || (str[i+4] == '/') || (str[i+4] == '^') || (str[i+4] == ')'))
								{
									return false;
								}
								else
								{
									log_open++;
									return true;
								}
							}
							else
							{
								return false;
							}
						}
						else
						{
							return false;
						}
					}
					else
					{
						return false;
					}
				}
				else if (isNum(temp_str))
				{
					
					operator_used_recently = false;
				}
				else
				{
					return false;
				}
				
				nest_opened_recently = false;
			}
		}
		return true;
	}
}



