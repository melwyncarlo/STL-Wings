#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "MISC.H"

#ifndef STLAERO_H
#define STLAERO_H

using namespace std;


#define PI 3.14159265
const int ARRAY_NUM = 19;
const int MAX_DEC_PTS = 6; // Maximum Decimal Points


string stl_type_value(double val);
string stl_facet_code(double f1, double f2, double f3, double vx1, double vy1, double vz1, double vx2, double vy2, double vz2, double vx3, double vy3, double vz3);
int centroid_array_num(double x_val, double normalised_x[], int size);
double normalise_angle(double angle);
double** twisted_airfoil_points(double **points, double twist_angle, double twist_pos_x, double twist_pos_y);
double absolute_angle(double input_angle, double dx, double dz);
int quadrant_mode(double dx, double dz);
int dir_val(double dx, double dz, int mode);
int dir_val_angle(double angle, int mode);
double cartesian_angle(double dx, double dy);


/* ---------------------------------------------------------------------------------- */


double roundn(double val1, double val2);
int ccw_z_sub_func(double val1z, double val2z, int opt1, int opt2, int mode, bool lhs);
int ccw_yz_sub_func(double val1y, double val2y, double val1z, double val2z, int opt1, int opt2, int mode, bool lhs);
int ccw_xyz_sub_func(double val1x, double val2x, double val1y, double val2y, double val1z, double val2z, int opt1, int opt2, int mode, bool lhs);
void ccw_func(int &output1, int &output2, int &output3, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz);
void unit_vector_func(double &output1, double &output2, double &output3, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz);
void ccw_vertices(double &ax, double &ay, double &az, double &bx, double &by, double &bz, double &cx, double &cy, double &cz);
void fill_vector(double (&output)[9], double val1, double val2, double val3, double val4, double val5, double val6, double val7, double val8, double val9);


/* ---------------------------------------------------------------------------------- */



string start_stl(string name);
string end_stl(string name);
void stl_gaps(double &input1, double &input2, double input3, double input4, double input5);
void stl_preprocess(double **points, double (&input1)[ARRAY_NUM], double (&input2)[ARRAY_NUM], double (&input3)[ARRAY_NUM], double input4, double input5, double input6, double &input7, double &input8);
string profile_stl(double **pts, double y_pos, bool y_dir_positive, double t_main_percent, double t_by_c_percent, double main_chord, double relative_chord);
string span_stl(double **pts1, double **pts2, double y1, double y2, double input3, double input4, double input5, double input6);


#endif



