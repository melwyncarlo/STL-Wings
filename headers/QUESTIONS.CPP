#include <iostream>
#include <string>
#include <stdlib.h>
#include "MISC.H"
#include "QUESTIONS.H"

using namespace std;

/*
Span Sections
Wing Twist (Wing Twist Position 0.5l or 0.25l, etc.)
Wing Thickness
(Create an airfoil for each section and twist, (thickness and profile))

For options, allow '0' as default
*/



void func_project_name(string &input_1)
{
	cout << "Note 1. The name must contain only alphanumeric characters." << endl;
	cout << "Note 2. The first character must be a letter." << endl;
	repeat:
	cout << "Enter a name for this Project: ";
	cin >> input_1;
	if (!isalpha(input_1[0]))
	{
		cout << endl << "Incorrect Project Name. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		for (int i=0; i<(input_1.length()); i++)
		{
			if ((!isdigit(input_1[i])) && (!isalpha(input_1[i])))
			{
				cout << endl << "Incorrect Project Name. Please try again." << endl << endl;
				goto repeat;
			}
		}
	}
}


void func_airfoil_name(string &input_1)
{
	cout << "Enter Airfoil Name: ";
	cin >> input_1;
	return;
}


void func_file_name(string &input_1)
{
	cout << "Enter File Name to save (without extension): ";
	cin >> input_1;
	return;
}


void func_nodes_n(int &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter the minimum number of nodes (points / 0.5 * Segments) required around an airfoil (in a general sense): ";
	cin >> temp_var;
	if ((temp_var < 50.0) || (temp_var > 10000.0) || ((temp_var - ((int)temp_var)) != 0.0) || (!isNum(to_string(temp_var))))
	{
		cout << endl << "Minimum nodes must be an integer between 50 and 10,000. ";
		cout << "Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
	return;
}


void func_le_nodes(double &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter a factor of the minimum number of nodes required concentrated around the Leading Edge of an airfoil: ";
	cin >> temp_var;
	if ((temp_var < 0.5) || (temp_var > 50.0) || (!isNum(to_string(temp_var))))
	{
		cout << endl << "Factor must be and integer between 0.5 and 50. ";
		cout << "Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
}


void func_resolution(int &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter number of pixels ('x' value) for an 'x' by 'x' image: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var > 10000.0) || (temp_var < 1.0))
	{
		cout << endl << "Number of Pixels (px) must be between 1px and 10,000px. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
}


void func_thickness(double &input_1, int input_2) // Input-2 is Airfoil Type
{
	double temp_var;
	repeat:
	cout << "Enter maximum thickness as percentage of chord: ";
	cin >> temp_var;
	if (input_2 == 4)
	{
		if (((temp_var - ((int)temp_var)) > 0.0) || (temp_var > 50.0) || (temp_var < 0.00000001) || (!isNum(to_string(temp_var))))
		{
			cout << endl << "Maximum Thickness (\%) must be an integer between (> 0) and 50 for NACA Airfoils. ";
			cout << "Please try again." << endl << endl;
			goto repeat;
		}
		else
		{
			input_1 = (int)temp_var;
		}
	}
	else
	{
		if ((!isNum(to_string(temp_var))) || (temp_var > 50.0) || (temp_var < 0.00001))
		{
			cout << endl << "Maximum Thickness (\%) must be between 0 and 50. Please try again." << endl << endl;
			goto repeat;
		}
		else
		{
			input_1 = temp_var;
		}
	}
}


void func_camber_height(double &input_1, double input_2) // Input-2 is Thickness Percentage
{
	double temp_var;
	repeat:
	cout << "Enter maximum camber as percentage of chord: ";
	cin >> temp_var;
	if (((temp_var + input_2) > 50.0) || (temp_var < 0.0) || (!isNum(to_string(temp_var))))
	{
		cout << endl << "Maximum Camber (\%) exceeds the limits (0 to ~30). Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
}


void func_camber_height_naca(double &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter maximum camber as percentage of chord: ";
	cin >> temp_var;
	if (((temp_var - ((int)temp_var)) > 0.0) || (temp_var > 9.0) || (temp_var < 0) || (!isNum(to_string(temp_var))))
	{
		cout << endl << "Maximum Camber (\%) must be between 0 and 9 for NACA airfoils. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
}


void func_camber_point(double &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter maximum camber distance from leading edge as percentage of chord: ";
	cin >> temp_var;
	if ((temp_var > 99) || (temp_var < 0.0) || (!isNum(to_string(temp_var))))
	{
		cout << endl << "Maximum Camber Distance (\%) must be between 0 and 99. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
}


void func_camber_point_naca(double &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter maximum camber distance from leading edge as percentage of chord: ";
	cin >> temp_var;
	if (((temp_var - ((int)temp_var)) > 0.0) || ((((int)temp_var) % 10) != 0) || (temp_var < 0.0) || (temp_var > 90.0) || (!isNum(to_string(temp_var))))
	{
		cout << endl << "Maximum Camber Distance (\%) must be and integer between 0 and 90 for NACA airfoils. ";
		cout << "Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
}


void func_chord(double &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter wing chord length in millimetres (mm): ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < 0) || (temp_var > 10000))
	{
		cout << endl << "The wing chord must be a number between 0 mm and 10,000 mm. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
}


void func_span(double &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter wing span length in millimetres (mm): ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < 0) || (temp_var > 50000))
	{
		cout << endl << "The wing chord must be a number between 0 mm and 50,000 mm. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
}


void func_sections_n(int &input_1, int input_2) // Input-2 is Maximum Number of Sections
{
	double temp_var;
	repeat:
	cout << "Enter Number of Sections: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < 0) || (temp_var > input_2))
	{
		cout << endl << "The number of sections must be a number between 1 and 10. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
}


/*
	Input-2	: Point of Starting Section Edge
	Input-3	: Wing Span
*/
void func_section_edge(double &input_1, double input_2, double input_3)
{
	double temp_var;
	repeat:
	cout << "Distance (mm) from point of starting section edge to point of ending section edge: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < 0.000001) || (temp_var > (input_3 - input_2)))
	{
		cout << endl << "The distance value must be within the limits of the wing span. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
}


void func_section_nodes(int &input_1, int input_2, int input_3)
{
	double temp_var;
	repeat:
	cout << "Enter minimum number of Section Nodes (points / 0.5 * Segments) required for Section-" + to_string(input_2) + ": ";
	cin >> temp_var;
	if ((temp_var < 10.0) || (temp_var > 10000.0) || ((temp_var - ((int)temp_var)) != 0.0) || (!isNum(to_string(temp_var))))
	{
		cout << endl << "Minimum Section Nodes must be an integer between 10 and 10,000. ";
		cout << "Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
	return;
}


void func_taper_ratio(double &input_1)
{
	double temp_var;
	cout << "For a straight wing, taper ratio = 1.0" << endl;
	repeat:
	cout << "Enter wing taper ratio: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < 0))
	{
		cout << endl << "The wing taper ratio must be a number greater than 0. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
	return;
}


void func_twist_extent(double &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter wing twist angle (degrees): ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < -90) || (temp_var > 90))
	{
		cout << endl << "The wing twist angle must be a number between -90 and +90. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
	return;
}


void func_twist_point_X(double &input_1)
{
	double temp_var;
	repeat:
	cout << "Enter X-position of twist from leading edge as percentage (\%) of chord: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < 0) || (temp_var > 100))
	{
		cout << endl << "The X-position of wing twist must be a percentage value 0\% and 100\%. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
	return;
}


void func_twist_point_Y(double &input_1)
{
	double temp_var;
	cout << "Option '10' implies a Y-Position value at half-thickness point along the Y-axis at the Twist X-Position Point." << endl;
	repeat:
	cout << "Enter Y-position of twist of values between '-1 to 1' (including '0') or '10': ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < -1) || ((temp_var > 1) && (temp_var != 10)))
	{
		cout << endl << "The Y-position of wing twist must be a value within the set limits. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
	return;
}


void func_span_division_edge(double &input_1)
{
	double temp_var;
	repeat:
	cout << "E.g. 1\% and 5\% division implies dividing the span length into 100 and 20 divisions respectively." << endl;
	cout << "Enter the Resolution of Span Function as Percentage Length (\%) of the divisions of the Span Length: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var < 0.00000001) || (temp_var > 100))
	{
		cout << endl << "The Resolution Percentage value must be real numbers greater than zero (> 0) and less than or equal to hundred.";
		cout << "(<= 100). Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = temp_var;
	}
}


/* ---------------------------------------------------------------------------------- */


void func_airfoil_type(int &input_1)
{
	double temp_var;
	cout << "Types of Airfoils - " << endl;
	cout << "1. Flat Plate" << endl;
	cout << "2. Parabolically Cambered Plate" << endl;
	cout << "3. Quadratically Cambered Plate" << endl;
	cout << "4. Cambered Four Digit NACA Airfoil" << endl;
	repeat:
	cout << "Enter type (number) of Airfoil: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var <= 0.0) || (temp_var > 4.0) || ((temp_var - ((int)temp_var)) != 0.0))
	{
		cout << endl << endl << "Incorrect Airfoil Type. Please try again." << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
	return;
}


void func_stl_choice(int &input_1)
{
	double temp_var;
	cout << "1. Yes" << endl;
	cout << "2. No" << endl;
	repeat:
	cout << "Do you wish to create an STL File? Enter '1' or '2': ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || ((temp_var - ((int)temp_var)) > 0.0)  || (temp_var < 1) || (temp_var > 2))
	{
		cout << endl << endl << "The option must be an integer value of '1' or '2'. Please try again." << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
	return;
}


void func_wing_type(int &input_1)
{
	double temp_var;
	cout << "1. Straight / Simply Tapered" << endl;
	cout << "2. Mathematical Function based" << endl;
	cout << "3. Tapered at Sections (Maximum of 10 Sections)" << endl;
	cout << "4. Different Airfoil Profile at Sections (Maximum of 10 Sections)" << endl;
	cout << "5. Combination of Options 1, 2 and 3 at Sections (Maximum of 10 Sections)" << endl;
	repeat:
	cout << "Enter Type (Number) of Wing Span: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var <= 0.0) || (temp_var > 5.0) || ((temp_var - ((int)temp_var)) != 0.0))
	{
		cout << endl << "Incorrect Span Type. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
	return;
}


void func_section_style(int &input_1, int input_2) // Input-2 is Section Number
{
	double temp_var;
	cout << "1. Using a Mathematical Function" << endl;
	cout << "2. Using Taper Ratio (Taper Ratio w.r.t. the Chord Length of Previous Section.";
	cout << "   For a Straight Wing Section, Taper Ratio = 1.0.)" << endl;
	cout << "3. Using an Airfoil Profile and a Chord Value" << endl;
	repeat:
	cout << "Which of the above process (Number) would you like to follow for Section-" + to_string(input_2) + ": ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var <= 0.0) || (temp_var > 3.0) || ((temp_var - ((int)temp_var)) != 0.0))
	{
		cout << endl << "Incorrect Section Process Type. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
	return;
}


void func_airfoil_dev_choice(int &input_1)
{
	double temp_var;
	cout << "1. Default Wing Root Airfoil Profile" << endl;
	cout << "2. New Airfoil Profile" << endl;
	repeat:
	cout << "Enter the type (number) of Airfoil Profile Required for this section: ";
	cin >> temp_var;
	if ((!isNum(to_string(temp_var))) || (temp_var <= 0.0) || (temp_var > 2.0) || ((temp_var - ((int)temp_var)) != 0.0))
	{
		cout << endl << "Incorrect Airfoil Profile Type. Please try again." << endl << endl;
		goto repeat;
	}
	else
	{
		input_1 = (int)temp_var;
	}
	return;
}


/* ---------------------------------------------------------------------------------- */


void func_wing_equation_intro()
{
	cout << "Type 'e' for Exponential" << endl;
	cout << "Type 'ln<f(x)>' for Natural Logarithm" << endl;
	cout << "Type 'log<f(x)>' for Common (Base-10) Logarithm" << endl;
	cout << "All functions must be in lowercase only" << endl;
	cout << "All functions must be encased in brackets and seperated for summation" << endl;
	cout << "For e.g., 'a + (b * (c - d))' translates to '(a) + (b*c) - (b*d)'" << endl;
	cout << "Boundary Condition - 1: The origin (x=0) is taken at point 0.5c (50\% of chord)" << endl;
	cout << "Boundary Condition - 2: The initial amplitudes are taken at point 0 and point c" << endl;
	cout << "Boundary Condition - 3: The wing span function, say, f(x), must always have a real value";
}


void func_wing_equation(string &input_1, int input_2) // Input-2 is Section Number
{
	repeat:
	if (input_2 == 0)
	{
		cout << "Enter wing span function: ";
	}
	else
	{
		cout << "Enter wing span function for Section - " + to_string(input_2) + ": ";
	}
	cin >> input_1;
	if (!equation_valid(input_1))
	{
		cout << endl << "The equation is invalid. Please follow the format and try again." << endl << endl;
		goto repeat;
	}
}



