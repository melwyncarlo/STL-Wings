#include <iostream>
#include <string>
#include <stdlib.h>
#include "MISC.H"

#ifndef QUESTIONS_H
#define QUESTIONS_H

using namespace std;


void func_project_name(string &input_1);
void func_airfoil_name(string &input_1);
void func_file_name(string &input_1);
void func_nodes_n(int &input_1);
void func_le_nodes(double &input_1);
void func_resolution(int &input_1);
void func_thickness(double &input_1, int input_2); // Input-2 is Airfoil Type
void func_camber_height(double &input_1, double input_2); // Input-2 is Thickness Percentage
void func_camber_height_naca(double &input_1);
void func_camber_point(double &input_1);
void func_camber_point_naca(double &input_1);
void func_chord(double &input_1);
void func_span(double &input_1);
void func_sections_n(int &input_1, int input_2);
void func_section_edge(double &input_1, double input_2, double input_3); // Input-2 is Point of Starting Section Edge, and Input-3 is Wing Span
void func_section_nodes(int &input_1, int input_2, int input_3);
void func_taper_ratio(double &input_1);
void func_twist_extent(double &input_1);
void func_twist_point_X(double &input_1);
void func_twist_point_Y(double &input_1);
void func_span_division_edge(double &input_1);


/* ---------------------------------------------------------------------------------- */


void func_airfoil_type(int &input_1);
void func_stl_choice(int &input_1);
void func_wing_type(int &input_1);
void func_section_style(int &input_1, int input_2); // Input-2 is Section Number
void func_airfoil_dev_choice(int &input_1);


/* ---------------------------------------------------------------------------------- */


void func_wing_equation_intro();
void func_wing_equation(string &input_1, int input_2); // Input-2 is Section Number


#endif


