#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include "headers/QUESTIONS.H"
#include "headers/MISC.H"
#include "headers/AIRFOILS.H"
#include "headers/SVGAERO.H"
#include "headers/STLAERO.H"

using namespace std;


const int MAX_SECTIONS = 10;

/*
	Number of Wing Span Function Types:
	1 = LHS
	2 = RHS
*/
const int M = 2;


/*
	g++ -o test airfoil_plot_V6.0.cpp headers/QUESTIONS.CPP headers/MISC.CPP headers/AIRFOILS.CPP headers/SVGAERO.CPP headers/STLAERO.CPP
*/


double** airfoil_type_01(double input01, int input02)
{
	cout << endl << "--------------------------------------------------------------------------------------------------";
	cout << endl << endl << "Processing . . . . . " << endl << endl;
	return (flatPlate(input01, input02));
}


double** airfoil_type_02(double &input01, double input02, int input03)
{
	cout << endl;
	func_camber_height(input01, input02);
	cout << endl << "--------------------------------------------------------------------------------------------------";
	cout << endl << endl << "Processing . . . . . " << endl << endl;
	return (parabolicallyCamberedPlate(input01, input02, input03));
}


double** airfoil_type_03(double &input01, double &input02, double input03, int input04)
{
	cout << endl;
	func_camber_height(input01, input03);
	cout << endl;
	func_camber_point(input02);
	cout << endl << "--------------------------------------------------------------------------------------------------";
	cout << endl << endl << "Processing . . . . . " << endl << endl;
	return (quadraticallyCamberedPlate(input01, input02, input03, input04));
}


double** airfoil_type_04(double &input01, double &input02, double input03, int input04, double input05, string &input06, string &input07)
{
	func_camber_height_naca(input01);
	cout << endl;
	func_camber_point_naca(input02);
	cout << endl << "--------------------------------------------------------------------------------------------------";
	cout << endl << endl << "Processing . . . . . " << endl << endl;
	input06 = "NACA-" + nacaCode(input01, input02, input03);
	input07 = "NACA_" + nacaCode(input01, input02, input03);
	return (nacaAirfoil(input01, input02, input03, input04, input05));
}


/*
	Input-1	: Airfoil Type
	Input-2	: Camber Height Percentage
	Input-3	: Camber Position Percentage
	Input-4	: Thickness Percentage
	Input-5	: Total Number of Nodes on Airfoil
	Input-6	: Factor representing Percentage of Total Nodes around Leading Edge
	Input-7	: Airfoil Name
	Input-8	: File Name
*/
double** airfoil_types(int input1, double &input2, double &input3, double input4, int input5, double input6, string &input7, string &input8)
{
	if (input1 == 1)
	{
		return airfoil_type_01(input4, input5);
	}
	else if (input1 == 2)
	{
		return airfoil_type_02(input2, input4, input5);
	}
	else if (input1 == 3)
	{
		return airfoil_type_03(input2, input3, input4, input5);
	}
	else
	{
		return airfoil_type_04(input2, input3, input4, input5, input6, input7, input8);
	}
}


/*
	Input-1	: Taper Ratio
	Input-2	: Twist Angle
	Input-3	: Twist Position X
	Input-4	: Twist Position Y
*/
void span_type_01(double &input1, double &input2, double &input3, double &input4)
{
	func_taper_ratio(input1);
	space();
	func_twist_extent(input2);
	space();
	func_twist_point_X(input3);
	space();
	func_twist_point_Y(input4);
	space();

	input4 = 100 * input4;
}


/*
	Input-1	: Equation String LHS
	Input-2	: Equation String RHS
	Input-3	: Section Number
	Input-4	: Twist Angle
	Input-5	: Twist Position X
	Input-6	: Twist Position Y
*/
void span_type_02(string &input1, string &input2, int input3, double &input4, double &input5, double &input6)
{
	func_wing_equation_intro();
	space();
	cout << "[ For the LHS Equation ]";
	func_wing_equation(input1, input3);
	space();
	cout << "[ For the RHS Equation ]";
	func_wing_equation(input2, input3);
	space();
	func_twist_extent(input4);
	space();
	func_twist_point_X(input5);
	space();
	func_twist_point_Y(input6);
	space();

	input4 = 100 * input6;
}


/*
	Input-1	: Number of Sections
	Input-2	: Section Length
	Input-3	: Taper Ratio
	Input-4	: Twist Angle
	Input-5	: Twist Position X
	Input-6	: Twist Position Y
	Input-7	: Wing Span
*/
void span_type_03_04_05(int &input1, double (&input2)[MAX_SECTIONS + 1], double (&input3)[MAX_SECTIONS], double &input4, double &input5, double &input6, double input7, int mode)
{
	func_sections_n(input1, MAX_SECTIONS);
	space();

	if (mode == 0)
	{
		for (int i=0; i<input1; i++)
		{
			func_taper_ratio(input3[i]);
		}
		space();
	}
	
	input2[0] = 0;
	double j_sum = 0;
	for (int j=1; j<input1; j++)
	{
		func_section_edge(input2[j], j_sum, input7);
		j_sum += input2[j];
	}
	input2[input1] = input7 - j_sum;
	cout << "For the Final Section: Distance (mm) from point of starting section edge to point of ending section edge is ";
	cout << input7 - j_sum << " mm." << endl;
	space();
	func_twist_extent(input4);
	space();
	func_twist_point_X(input5);
	space();
	func_twist_point_Y(input6);
	space();

	input6 = 100 * input6;
}


string equation_solver(string input_str, string x)
{
	/*string str = remove_spaces(input_str);
	string nests[50], operator_sign_local, operator_sign_global;
	int nest_n = 0;
	double  = 0;
	double accumulator = 0;
	double total = 0;
	bool nest_open = false;
	if (equation_valid(str))
	{
		for (int i=0; i<(str.length); i++)
		{
			if (str[i] == '(')
			{
				nest_open = "true";
			}
			else if (str[i] == ')')
			{
				nest_open = "false";
				nest_n++;
				if (nest_n == 1)
				{
					total = accumulator;
				}
				else
				{
					total = arithmetic_operation(total, accumulator, operator_sign_global);
				}
			}
			else
			{
				if (nest_open)
				{
					if ((str[i] == '+') || (str[i] == '-') || (str[i] == '*') || (str[i] == '/') || (str[i] == '^'))
					{
						operator_sign_local = str[i];
					}
					else if (str[i] == 'x')
					{
						accumulator = arithmetic_operation(accumulator, x, operator_sign_local);
					}
					else if (str[i] == 'e')
					{
						
					}
					else if (str[i] == 'l')
					{
						if (str[i+1] == 'n')
						{
							int sub_func_length = 1;
							
							i += (sub_func_length - 1);
						}
						else if (str[i+1] == 'o')
						{
							int sub_func_length = 1;
							
							i += (sub_func_length - 1);
						}
					}
					else
					{
						
					}
					accumulator = ;
				}
				else
				{
					operator_sign_global = str.substr(i, 1);
				}
			}
		}
	}
	else
	{
		return "NaN";
	}*/
	return input_str;
}


double find_y_pos(double **pts, double x_pos)
{
	double y1 = 0;
	double y2 = 0;
	double temp_var_01 = -100;
	bool one_found = false;
	bool direction_turned = false;
	for (int i=1; i<(pts[0][0]); i++)
	{
		double temp_var_02 = abs(pts[i][0] - x_pos);
		if (!one_found)
		{
			if (temp_var_02 >= temp_var_01)
			{
				y1 = pts[i][1];
				one_found = true;
			}
			else
			{
				temp_var_01 = temp_var_02;
			}
		}
		else
		{
			if (!direction_turned)
			{
				if (temp_var_02 < temp_var_01)
				{
					direction_turned = true;
				}
				temp_var_01 = temp_var_02;
			}
			else
			{
				if (temp_var_02 >= temp_var_01)
				{
					y2 = pts[i][1];
					break;
				}
				else
				{
					temp_var_01 = temp_var_02;
				}
			}
		}
	}
	return ((y1 + y2) / 2);
}


void message(string str, bool hr = true)
{
	space();
	if (hr)
	{
		cout << "--------------------------------------------------------------------------------------------------";
		space();
	}
	if (str.compare("") != 0)
	{
		cout << str;
		space();
	}
}



int main ()
{
	string project_name = "";
	string airfoil_name = "";
	string airfoil_name_temp = "";
	string file_name = "";
	string file_name_temp = "";
	string main_data = "";
	string main_data_temp = "";
	string svg_data_01 = "";
	string svg_data_02 = "";
	string svg_data_03 = "";
	string svg_data_01_temp = "";
	string svg_data_02_temp = "";
	string svg_data_03_temp = "";
	string stl_data = "";
	int airfoil_type, airfoil_type_temp;
	double grid_type;
	int resolution;
	int nodes;
	double factor, factor_temp;
	double t_percent, t_percent_temp, m_percent, m_percent_temp, p_percent, p_percent_temp;
	cout << endl << endl << "Welcome to Plotting Simple Airfoils" << endl << endl;

	func_project_name(project_name);
	space();
	
	func_airfoil_type(airfoil_type);
	if (airfoil_type != 4)
	{
		func_airfoil_name(airfoil_name);
		func_file_name(file_name);
	}
	cout << endl << "--------------------------------------------------------------------------------------------------";
	space();
	func_nodes_n(nodes);
	cout << endl;
	if (airfoil_type == 4)
	{
		func_le_nodes(factor);
		cout << endl;
	}
	func_resolution(resolution);
	space();
	func_thickness(t_percent, airfoil_type);
	space();
	
	double** points, **points_temp;
	double** section_airfoil_points[MAX_SECTIONS + 1];

	points = airfoil_types(airfoil_type, m_percent, p_percent, t_percent, nodes, factor, airfoil_name, file_name);
	// OR
	// section_airfoil_points[0] = airfoil_types(airfoil_type, m_percent, p_percent, t_percent, nodes, factor, airfoil_name, file_name);

	main_data = airfoil_name + "\n";
	main_data += data_points_to_string(points, points[0][0]);
	svg_data_01 = create_svg_airfoil(points, points[0][0], airfoil_name, t_percent, m_percent, p_percent, resolution, 1);
	svg_data_02 = create_svg_airfoil(points, points[0][0], airfoil_name, t_percent, m_percent, p_percent, resolution, 2);
	svg_data_03 = create_svg_airfoil(points, points[0][0], airfoil_name, t_percent, m_percent, p_percent, resolution, 3);
	
	cout << "Nearly Done . . . . . " << endl << endl;
	int stl_opt, span_type, sections_n;
	double wing_chord, wing_chord_1, wing_chord_2, wing_span, twist, twist_pos_x, twist_pos_y;
	double taper_ratio[MAX_SECTIONS], sections_l[MAX_SECTIONS + 1], section_l_sum;	
	double divisions, chord_left, chord_right;
	double division_l_local_percent, section_l_local_percent, local_twist;
	string wing_equations[M];
	func_stl_choice(stl_opt);

	if (stl_opt == 1)
	{
		section_airfoil_points[0] = points;

		cout << endl << "--------------------------------------------------------------------------------------------------";
		space();
		func_chord(wing_chord);
		space();
		func_wing_type(span_type);
		space();
		func_span(wing_span);
		space();

		stl_data = start_stl(project_name);
		stl_data += profile_stl(points, 0, false, t_percent, t_percent, wing_chord, wing_chord);

		if (span_type == 1)
		{
			span_type_01(taper_ratio[0], twist, twist_pos_x, twist_pos_y);
			if (twist_pos_y == 1000)
			{
				twist_pos_y = 100 * find_y_pos(points, twist_pos_x / 100);
			}

			wing_chord_1 = wing_chord;
			wing_chord_2 = wing_chord * taper_ratio[0];
			section_airfoil_points[1] = twisted_airfoil_points(points, twist, twist_pos_x / 100, twist_pos_y / 100);
			//stl_data = create_stl(points, double t_percent, wing_chord_1, wing_chord_2, wing_span, project_name);
			
			message("Processing again . . . . . ");

			stl_data += profile_stl(section_airfoil_points[1], wing_span, true, t_percent, t_percent, wing_chord, wing_chord_2);
			stl_data += span_stl(section_airfoil_points[0], section_airfoil_points[1], 0, wing_span, wing_chord, wing_chord_1, 
					wing_chord_2, t_percent);
			
		}
		else if (span_type == 2)
		{
			span_type_02(wing_equations[0], wing_equations[1], 0, twist, twist_pos_x, twist_pos_y);
			if (twist_pos_y == 1000)
			{
				twist_pos_y = 100 * find_y_pos(points, twist_pos_x / 100);
			}

			func_span_division_edge(divisions);
			
			double divisions_l = wing_span / divisions;

			message("Processing . . . . . ");

			for (int i=0; i<divisions; i++)
			{
				if (i == 0)
				{
					wing_chord_1 = wing_chord;
				}
				else
				{
					wing_chord_1 = wing_chord_2;
				}

				chord_left = stod(equation_solver(wing_equations[0], to_string(divisions_l * (i + 1))));
				chord_right = stod(equation_solver(wing_equations[1], to_string(divisions_l * (i + 1))));
				double wing_chord_2 = chord_left + chord_right;

				// Division Length as Percentage of the Wing Span
				division_l_local_percent = ((divisions_l * (i + 1)) / wing_span) * 100;
				local_twist = twist * (division_l_local_percent / 100);
				section_airfoil_points[1] = twisted_airfoil_points(points, local_twist, twist_pos_x / 100, twist_pos_y / 100);

				stl_data += profile_stl(section_airfoil_points[1], divisions_l * (i + 1), true, t_percent, t_percent, 
						wing_chord, wing_chord_2);
				stl_data += span_stl(section_airfoil_points[0], section_airfoil_points[1], divisions_l * i, divisions_l * (i + 1), 
						wing_chord, wing_chord_1, wing_chord_2, t_percent);
			}
		}
		else if (span_type == 3)
		{
			span_type_03_04_05(sections_n, sections_l, taper_ratio, twist, twist_pos_x, twist_pos_y, wing_span, 0);
			if (twist_pos_y == 1000)
			{
				twist_pos_y = 100 * find_y_pos(points, twist_pos_x / 100);
			}
		
			message("Processing . . . . . ");
	
			for (int i=0; i<sections_n; i++)
			{
				if (i == 0)
				{
					wing_chord_1 = wing_chord;
				}
				else
				{
					wing_chord_1 = wing_chord_2;
				}
				wing_chord_2 = wing_chord * taper_ratio[i];

				section_l_sum += sections_l[i + 1];
				section_l_local_percent = (section_l_sum / wing_span) * 100;
				local_twist = twist * (section_l_local_percent / 100);
				section_airfoil_points[i + 1] = twisted_airfoil_points(points, local_twist, twist_pos_x / 100, twist_pos_y / 100);
				
				stl_data += profile_stl(section_airfoil_points[i + 1], section_l_sum, true, t_percent, t_percent, wing_chord, 
						wing_chord_2);
				stl_data += span_stl(section_airfoil_points[i], section_airfoil_points[i + 1], section_l_sum - sections_l[i + 1], 
						section_l_sum, wing_chord, wing_chord_1, wing_chord_2, t_percent);
			}
		}
		else if (span_type == 4)
		{
			span_type_03_04_05(sections_n, sections_l, taper_ratio, twist, twist_pos_x, twist_pos_y, wing_span, 1);
			if (twist_pos_y == 1000)
			{
				twist_pos_y = 100 * find_y_pos(points, twist_pos_x / 100);
			}

			message("Processing . . . . . ");

			section_l_sum = 0;
			for (int i=0; i<sections_n; i++)
			{
				cout << "[ PERTAINING TO SECTION-" << i + 1 << " ]";
				space();
				func_airfoil_type(airfoil_type_temp);
				space();
				if (airfoil_type_temp == 4)
				{
					func_le_nodes(factor_temp);
					space();
				}
				else
				{
					func_airfoil_name(airfoil_name_temp);
					space();
					func_file_name(file_name_temp);
					space();
				}
				func_thickness(t_percent_temp, airfoil_type_temp);
				space();

				if (i == 0)
				{
					wing_chord_1 = wing_chord;
				}
				else
				{
					wing_chord_1 = wing_chord_2;
				}
				
				points_temp = airfoil_types(airfoil_type_temp, m_percent_temp, p_percent_temp, t_percent_temp, 
						nodes, factor_temp, airfoil_name_temp, file_name_temp);

				func_chord(wing_chord_2);

				section_l_sum += sections_l[i + 1];
				section_l_local_percent = (section_l_sum / wing_span) * 100;
				local_twist = twist * (section_l_local_percent / 100);
				section_airfoil_points[i + 1] = twisted_airfoil_points(points_temp, local_twist, twist_pos_x / 100, twist_pos_y / 100);

				message("Processing Supplementary Files . . . . . ");

				main_data_temp = airfoil_name_temp + "\n";
				main_data_temp += data_points_to_string(points_temp, points[0][0]);
				svg_data_01_temp = create_svg_airfoil(points_temp, points[0][0], airfoil_name_temp, t_percent_temp, m_percent_temp, 
				p_percent_temp, resolution, 1);
				svg_data_02_temp = create_svg_airfoil(points_temp, points[0][0], airfoil_name_temp, t_percent_temp, m_percent_temp, 
				p_percent_temp, resolution, 2);
				svg_data_03_temp = create_svg_airfoil(points_temp, points[0][0], airfoil_name_temp, t_percent_temp, m_percent_temp, 
				p_percent_temp, resolution, 3);
				write_to_file("storage/airfoil_data", project_name, file_name_temp + "_SECTION_" + to_string(i + 1), 
						"dat", main_data_temp);
				write_to_file("storage/svg_data", project_name, file_name_temp + "_SECTION_" + to_string(i + 1) + 
						"_Simple_Grid", "svg", svg_data_01_temp);
				write_to_file("storage/svg_data", project_name, file_name_temp + "_SECTION_" + to_string(i + 1) + 
						"_Fine_Grid", "svg", svg_data_02_temp);
				write_to_file("storage/svg_data", project_name, file_name_temp + "_SECTION_" + to_string(i + 1) + 
						"_Finer_Grid", "svg", svg_data_03_temp);

				stl_data += profile_stl(section_airfoil_points[i + 1], section_l_sum, true, t_percent, t_percent_temp, wing_chord, 
						wing_chord_2);
				stl_data += span_stl(section_airfoil_points[i], section_airfoil_points[i + 1], section_l_sum - sections_l[i + 1], 
						section_l_sum, wing_chord, wing_chord_1, wing_chord_2, t_percent);
			}
		}
		else
		{
			// Taper Ratio is NOT used in this case.
			span_type_03_04_05(sections_n, sections_l, taper_ratio, twist, twist_pos_x, twist_pos_y, wing_span, 1);
			if (twist_pos_y == 1000)
			{
				twist_pos_y = 100 * find_y_pos(points, twist_pos_x / 100);
			}

			section_airfoil_points[0] = points;
			section_l_sum = 0;
			for (int i=0; i<sections_n; i++)
			{
				section_l_sum += sections_l[i + 1];	
				section_l_local_percent = (section_l_sum / wing_span) * 100;
				local_twist = twist * (section_l_local_percent / 100);
				double local_twist_ex = local_twist;

				if (i == 0)
				{
					wing_chord_1 = wing_chord;
				}
				else
				{
					wing_chord_1 = wing_chord_2;
				}
				int section_type;
				double local_taper_ratio;
				message("");
				cout << "[ PERTAINING TO SECTION-" << i + 1 << " ]";
				space();
				func_section_style(section_type, i + 1);
				space();

				int profile_type;
				func_airfoil_dev_choice(profile_type);
				space();
				if (profile_type == 1)
				{
					points_temp = points;
					t_percent_temp = t_percent;
				}
				else
				{
					func_airfoil_type(airfoil_type_temp);
					space();
					if (airfoil_type_temp == 4)
					{
						func_le_nodes(factor_temp);
						space();
					}
					else
					{
						func_airfoil_name(airfoil_name_temp);
						space();
						func_file_name(file_name_temp);
						space();
					}
					func_thickness(t_percent_temp, airfoil_type_temp);
					space();

					message("Processing Supplementary Files . . . . . ");
					message(" ", false);
					
					points_temp = airfoil_types(airfoil_type_temp, m_percent_temp, p_percent_temp, t_percent_temp, 
							nodes, factor_temp, airfoil_name_temp, file_name_temp);

					main_data_temp = airfoil_name_temp + "\n";
					main_data_temp += data_points_to_string(points_temp, points[0][0]);
					svg_data_01_temp = create_svg_airfoil(points_temp, points[0][0], airfoil_name_temp, t_percent_temp, 
								m_percent_temp, p_percent_temp, resolution, 1);
					svg_data_02_temp = create_svg_airfoil(points_temp, points[0][0], airfoil_name_temp, t_percent_temp, 
								m_percent_temp, p_percent_temp, resolution, 2);
					svg_data_03_temp = create_svg_airfoil(points_temp, points[0][0], airfoil_name_temp, t_percent_temp, 
								m_percent_temp, p_percent_temp, resolution, 3);
					write_to_file("storage/airfoil_data", project_name, file_name_temp + "_SECTION_" + to_string(i + 1), 
							"dat", main_data_temp);
					write_to_file("storage/svg_data", project_name, file_name_temp + "_SECTION_" + to_string(i + 1) + 
							"_Simple_Grid", "svg", svg_data_01_temp);
					write_to_file("storage/svg_data", project_name, file_name_temp + "_SECTION_" + to_string(i + 1) + 
							"_Fine_Grid", "svg", svg_data_02_temp);
					write_to_file("storage/svg_data", project_name, file_name_temp + "_SECTION_" + to_string(i + 1) + 
							"_Finer_Grid", "svg", svg_data_03_temp);
				}
				
				if (section_type == 1) // Mathematical Function
				{
					func_wing_equation_intro();
					space();
					cout << "[ For the LHS Equation ]";
					func_wing_equation(wing_equations[0], i + 1);
					space();
					cout << "[ For the RHS Equation ]";
					func_wing_equation(wing_equations[1], i + 1);
					space();
					func_span_division_edge(divisions);
					space();
					
					double divisions_l = wing_span / divisions;

					message("Processing . . . . . ", false);
  
					double local_twist_temp;
					double** section_airfoil_points_temp_ex, **section_airfoil_points_temp;
					section_airfoil_points_temp_ex = section_airfoil_points[i];
					for (int j=0; j<divisions; j++)
					{
						wing_chord_1 = wing_chord_2;

						chord_left = stod(equation_solver(wing_equations[0], to_string(divisions_l * (j + 1))));
						chord_right = stod(equation_solver(wing_equations[1], to_string(divisions_l * (j + 1))));
						wing_chord_2 = chord_left + chord_right;

						// Division Length as Percentage of the Wing Span
						division_l_local_percent = ((divisions_l * (j + 1)) / wing_span) * 100;
						local_twist_temp = local_twist_ex + ((local_twist - local_twist_ex) 
									* (division_l_local_percent / 100));
						section_airfoil_points_temp = twisted_airfoil_points(points_temp, local_twist_temp, twist_pos_x / 100, 
									twist_pos_y / 100);

						stl_data += profile_stl(section_airfoil_points_temp, divisions_l * (j + 1), true, t_percent, 
								t_percent_temp, wing_chord, wing_chord_2);
						stl_data += span_stl(section_airfoil_points_temp_ex, section_airfoil_points_temp, divisions_l * j, 
								divisions_l * (j + 1), wing_chord, wing_chord_1, wing_chord_2, t_percent);

						section_airfoil_points_temp_ex = section_airfoil_points_temp;
						if (j == (divisions - 1))
						{
							section_airfoil_points[i + 1] = section_airfoil_points_temp;
						}
					}
				}
				else if (section_type == 2) //  Taper Ratio
				{
					func_taper_ratio(local_taper_ratio);

					message("Processing . . . . . ", false);

					wing_chord_2 = wing_chord_1 * local_taper_ratio;

					section_airfoil_points[i + 1] = twisted_airfoil_points(points_temp, local_twist, twist_pos_x / 100, 
									twist_pos_y / 100);
					
					stl_data += profile_stl(section_airfoil_points[i + 1], section_l_sum, true, t_percent, t_percent_temp, 
							wing_chord, wing_chord_2);
					stl_data += span_stl(section_airfoil_points[i], section_airfoil_points[i + 1], 
							section_l_sum - sections_l[i + 1], section_l_sum, wing_chord, wing_chord_1, wing_chord_2, 
							t_percent);
					
				}
				else // if (section_type == 3) // Airfoil Profile
				{
					func_chord(wing_chord_2);

					message("Processing . . . . . ", false);

					section_airfoil_points[i + 1] = twisted_airfoil_points(points_temp, local_twist, twist_pos_x / 100, 
									twist_pos_y / 100);
					
					stl_data += profile_stl(section_airfoil_points[i + 1], section_l_sum, true, t_percent, t_percent_temp, 
							wing_chord, wing_chord_2);
					stl_data += span_stl(section_airfoil_points[i], section_airfoil_points[i + 1], 
							section_l_sum - sections_l[i + 1], section_l_sum, wing_chord, wing_chord_1, wing_chord_2, 
							t_percent);
				}

				if (i == (sections_n - 1))
				{
					message("");
				}
			}
		}

		stl_data += end_stl(project_name);
		
		message("Processing again . . . . . ", false);

		write_to_file("storage/airfoil_data", project_name, file_name, "dat", main_data);
		write_to_file("storage/svg_data", project_name, file_name + "_Simple_Grid", "svg", svg_data_01);
		write_to_file("storage/svg_data", project_name, file_name + "_Fine_Grid", "svg", svg_data_02);
		write_to_file("storage/svg_data", project_name, file_name + "_Finer_Grid", "svg", svg_data_03);
		write_to_file("storage/stl_data", "", project_name, "stl", stl_data);
		write_to_file("storage/stl_data", "", project_name, "txt", stl_data);
		
		cout << "Done!" << endl << endl;
		cout << "The following main files can be found in their respective directories (along with additional section-based files):" << endl;
		cout << "1. Airfoil Points Data : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/airfoil_data/" + project_name + "_" + file_name + ".dat"; 
		cout << endl;
		cout << "2. SVG Data (HTML) : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/svg_data/" + project_name + "_" + file_name;
		cout << "_Simple_Grid.svg" << endl;
		cout << "3. SVG Data : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/svg_data/" + project_name + "_" + file_name + "_Fine_Grid.svg";
		cout << endl;
		cout << "4. SVG Data : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/svg_data/" + project_name + "_" + file_name + "_Finer_Grid.svg";
		cout << endl;
		cout << "5. STL Data : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/stl_data/" + project_name + ".stl" << endl;
		cout << "6. STL Data (Text) : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/stl_data/" + project_name + ".txt" << endl;
		cout << endl << "--------------------------------------------------------------------------------------------------" << endl << endl;

	}
	else
	{

		message("Processing again . . . . . ");

		write_to_file("storage/airfoil_data", project_name, file_name, "dat", main_data);
		write_to_file("storage/svg_data", project_name, file_name + "_Simple_Grid", "svg", svg_data_01);
		write_to_file("storage/svg_data", project_name, file_name + "_Fine_Grid", "svg", svg_data_02);
		write_to_file("storage/svg_data", project_name, file_name + "_Finer_Grid", "svg", svg_data_03);

		cout << "Done!" << endl << endl;
		cout << "The following main files can be found in their respective directories (along with additional section-based files):" << endl;
		cout << "1. Airfoil Points Data : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/airfoil_data/" + project_name + "_" + file_name;
		cout << ".dat" << endl;
		cout << "2. SVG Data (HTML) : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/svg_data/" + project_name + "_" + file_name;
		cout << "_Simple_Grid.svg" << endl;
		cout << "3. SVG Data : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/svg_data/" + project_name + "_" + file_name + "_Fine_Grid.svg";
		cout << endl;
		cout << "4. SVG Data : [airfoil_plot_V3.0.cpp DIRECTORY]/storage/svg_data/" + project_name + "_" + file_name + "_Finer_Grid.svg";
		cout << endl;
		cout << endl << "--------------------------------------------------------------------------------------------------" << endl << endl;

	}

	return 0;
}
