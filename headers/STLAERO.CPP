#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "MISC.H"
#include "STLAERO.H"


using namespace std;


#define PI 3.14159265


string stl_type_value(double val)
{
	string str = "";
	int number_of_digits = to_string(abs((int)val)).length();
	int division_factor = pow(10, number_of_digits - 1);
	double val_single_precision = val / division_factor;
	int e_factor = (int)(log10(division_factor));
	string e_value = repeated_string("0", 3 - to_string(e_factor).length()) + to_string(e_factor);
	str = decimal_point(to_string(val_single_precision), 10) + "e+" + e_value;

	return str;
}


string stl_facet_code(double f1, double f2, double f3, double vx1, double vy1, double vz1, double vx2, double vy2, double vz2, double vx3, double vy3, double vz3)
{
	string str = "";
	str += "     facet normal " + stl_type_value(f1) + " " + stl_type_value(f2) + " " + stl_type_value(f3) + " \n";
	str += "          outer loop \n";
	str += "               vertex " + stl_type_value(vx1) + " " + stl_type_value(vy1) + " " + stl_type_value(vz1) + " \n";
	str += "               vertex " + stl_type_value(vx2) + " " + stl_type_value(vy2) + " " + stl_type_value(vz2) + " \n";
	str += "               vertex " + stl_type_value(vx3) + " " + stl_type_value(vy3) + " " + stl_type_value(vz3) + " \n";
	str += "          end loop \n";
	str += "     end facet \n\n";

	return str;
}


int centroid_array_num(double x_val, double normalised_x[], int size)
{
	int modified_x_val = (int)(x_val * 100);
	if (x_val >= normalised_x[1])
	{
		return 0;
	}
	else if (x_val < normalised_x[size-1])
	{
		return (size-1);
	}
	else
	{
		for (int i=2; i<size; i++)
		{
			int modified_norm_x = (int)(normalised_x[i] * 100);
			if (modified_x_val >= modified_norm_x)
			{
				return (i - 1);
			}
		}
	}
}


double normalise_angle(double angle)
{
	double new_angle;
	if ((angle > 360)) // && (angle <= 450))
	{
		new_angle = angle - 360;
	}
	else
	{
		new_angle = angle;
	}

	if (new_angle <= 90)
	{
		return (new_angle);
	}
	else if (new_angle <= 180)
	{
		return (180 - new_angle);
	}
	else if (new_angle <= 270)
	{
		return (new_angle - 180);
	}	
	else // if (new_angle <= 360)
	{
		return (360 - new_angle);
	}
}


double absolute_angle(double input_angle, double dx, double dz)
{
	if ((dx >= 0) && (dz >= 0))
	{
		return (input_angle);
	}
	else if ((dx < 0) && (dz >= 0))
	{
		return (180 - input_angle);
	}
	else if ((dx < 0) && (dz < 0))
	{
		return (180 + input_angle);
	}
	else // if ((dx >= 0) && (dz < 0))
	{
		return (360 - input_angle);
	}
}


/*
   	Absolute Twist Angle is Positive (+ve) for CW Twist and Negative (-ve) for CCW Twist
	Twist Postion-Y	should be 0 as Default, for a normalised airfoil
*/
double** twisted_airfoil_points(double **points, double absolute_twist_angle, double twist_pos_x, double twist_pos_y)
{
	if (absolute_twist_angle != 0)
	{

		double **array = new double*[(int)points[0][0]];
		for (int n=0; n<=(int)(points[0][0]); n++)
		{
			array[n] = new double[2];
		}

		array[0][0] = points[0][0];

		for (int i=1; i<=(points[0][0]); i++)
		{
			/*double abs_angle;
			double dx = points[i][0] - twist_pos_x;
			double dy = points[i][1] - twist_pos_y;
			double r_val = pow(pow(dx, 2) + pow(dy, 2), 0.5);
			double angle = atan(abs(dy) / abs(dx)) * 180 / PI;
			double dx2, dy2, twist_angle, normalised_twist_angle;

			abs_angle = absolute_angle(angle, dx, dy);

			twist_angle = abs_angle - absolute_twist_angle;

			if (twist_angle < 0)
			{
				twist_angle = 360 - (abs(twist_angle));
			}

			normalised_twist_angle = normalise_angle(twist_angle);
			
			dx2 = abs(r_val * cos(normalised_twist_angle * PI / 180));
			dy2 = abs(r_val * sin(normalised_twist_angle * PI / 180));
			
			if (twist_angle <= 90)
			{
			    dx2 = twist_pos_x + dx2;
			    dy2 = twist_pos_y + dy2;
			}
			else if (twist_angle <= 180)
			{
			    dx2 = twist_pos_x - dx2;
			    dy2 = twist_pos_y + dy2;
			}
			else if (twist_angle <= 270)
			{
			    dx2 = twist_pos_x - dx2;
			    dy2 = twist_pos_y - dy2;
			}
			else
			{
			    dx2 = twist_pos_x + dx2;
			    dy2 = twist_pos_y - dy2;
			}
			
			array[i][0] = dx2;
			array[i][1] = dy2;*/

			array[i][0] = twist_pos_x + ((points[i][0] * cos(-absolute_twist_angle * PI / 180)) 
						- (points[i][1] * sin(-absolute_twist_angle * PI / 180)));
			array[i][1] = twist_pos_y + ((points[i][1] * cos(-absolute_twist_angle * PI / 180)) 
						+ (points[i][0] * sin(-absolute_twist_angle * PI / 180)));
		}

		return array;

	}
	else
	{
		return points;
	}

}


int quadrant_mode(double dx, double dz)
{
	if ((dx >= 0) && (dz < 0))
	{
		return 3;
	}
	else if ((dx < 0) && (dz >= 0))
	{
		return 1;
	}
	else if ((dx < 0) && (dz < 0))
	{
		return 2;
	}
	else // if ((dx >= 0) && (dz >= 0))
	{
		return 4;
	}
}


/*
	Mode-1	: X-Axis
	Mode-2	: Z-Axis
*/
int dir_val(double dx, double dz, int mode)
{
	int x, z;
	if ((dx >= 0) && (dz < 0))
	{
		x = -1;
		z = -1;
	}
	else if ((dx < 0) && (dz >= 0))
	{
		x = 1;
		z = 1;
	}
	else if ((dx < 0) && (dz < 0))
	{
		x = -1;
		z = 1;
	}
	else // if ((dx >= 0) && (dz >= 0))
	{
		x = 1;
		z = -1;
	}

	if (mode == 2)
	{
		return z;
	}
	else // if (mode == 1)
	{
		return x;
	}
}


int dir_val_angle(double angle, int mode)
{
	int x, z;
	double new_angle;
	if ((angle > 360)) // && (angle <= 450))
	{
		new_angle = angle - 360;
	}
	else
	{
		new_angle = angle;
	}

	if (new_angle <= 90)
	{
		x = 1;
		z = 1;
	}
	else if (new_angle <= 180)
	{
		x = -1;
		z = 1;
	}
	else if (new_angle <= 270)
	{
		x = -1;
		z = -1;
	}
	else // if (new_angle <= 360)
	{
		x = 1;
		z = -1;
	}

	if (mode == 2)
	{
		return z;
	}
	else // if (mode == 1)
	{
		return x;
	}
}


double cartesian_angle(double x1, double y1, double x2, double y2)
{
	if ((stod(decimal_point(to_string(x1), MAX_DEC_PTS)) == 0) && (stod(decimal_point(to_string(x2), MAX_DEC_PTS)) == 0))
	{
		return 90;
	}
	else if ((x1 == x2) && (y1 == y2))
	{
		return 0;
	}
	else
	{
		double angle =  (atan(stod(decimal_point(to_string(abs(y2 - y1)), MAX_DEC_PTS)) / stod(decimal_point(to_string(abs(x2 - x1)), 
				MAX_DEC_PTS))) * 180 / PI);
		return (absolute_angle(angle, x2 -  x1, y2 - y1));
	}
}



/* ---------------------------------------------------------------------------------- */


double roundn(double val1, double val2)
{
	return (((int)(rint(val1 * pow(10, val2)))) / pow(10, val2));
}


int ccw_z_sub_func(double val1z, double val2z, int opt1, int opt2, int mode, bool lhs)
{
	int output_a, output_b;
	if (opt2 < opt1)
	{
		double temp = opt1;
		opt1 = opt2;
		opt2 = temp;
	}

	if (val1z > val2z)
	{
		if (lhs)
		{
			output_a = opt1;
			output_b = opt2;
		}
		else
		{
			output_a = opt2;
			output_b = opt1;
		}
	}
	else if (val2z > val1z)
	{
		if (lhs)
		{
			output_a = opt2;
			output_b = opt1;
		}
		else
		{
			output_a = opt1;
			output_b = opt2;
		}
	}
	else // if (val1z == val2z)
	{
		output_a = opt1;
		output_b = opt1;
	}

	if (mode == 1)
	{
		return output_a;
	}
	else // if (mode == 2)
	{
		return output_b;
	}
}


int ccw_yz_sub_func(double val1y, double val2y, double val1z, double val2z, int opt1, int opt2, int mode, bool lhs)
{
	int output_a, output_b;
	if (opt2 < opt1)
	{
		double temp = opt1;
		opt1 = opt2;
		opt2 = temp;
	}

	if (val1y > val2y)
	{
		if (lhs)
		{
			output_a = opt1;
			output_b = opt2;
		}
		else
		{
			output_a = opt2;
			output_b = opt1;
		}
	}
	else if (val2y > val1y)
	{
		if (lhs)
		{
			output_a = opt2;
			output_b = opt1;
		}
		else
		{
			output_a = opt1;
			output_b = opt2;
		}
	}
	else // if (val1y == val2y)
	{
		output_a = ccw_z_sub_func(val1z, val2z, opt1, opt2, 1, lhs);
		output_b = ccw_z_sub_func(val1z, val2z, opt1, opt2, 2, lhs);
	}

	if (mode == 1)
	{
		return output_a;
	}
	else // if (mode == 2)
	{
		return output_b;
	}
}


int ccw_xyz_sub_func(double val1x, double val2x, double val1y, double val2y, double val1z, double val2z, int opt1, int opt2, int mode, bool lhs)
{
	int output_a, output_b;
	if (opt2 < opt1)
	{
		double temp = opt1;
		opt1 = opt2;
		opt2 = temp;
	}

	if (val1x < val2x)
	{
		output_a = opt1;
		output_b = opt2;
	}
	else if (val2x < val1x)
	{
		output_a = opt2;
		output_b = opt1;
	}
	else // if (val1x == val2x)
	{
		output_a = ccw_yz_sub_func(val1y, val2y, val1z, val2z, opt1, opt2, 1, lhs);
		output_b = ccw_yz_sub_func(val1y, val2y, val1z, val2z, opt1, opt2, 2, lhs);
	}

	if (mode == 1)
	{
		return output_a;
	}
	else // if (mode == 2)
	{
		return output_b;
	}
}


void ccw_func(int &output1, int &output2, int &output3, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
	double determinant;
	if ((bx == 0) && (cx == 0)) // OR ((ax == 0) && (bx == 0) && (cx == 0))
	{
	    
		determinant = (ax * ((by * cz) - (bz * cy))) - (ay * ((bx * cz) - (bz * cx))) + (az * ((bx * cy) - (by * cx)));
	}
	else
	{
		double arbitrary_val = 100;
		determinant = (ax * ((by * arbitrary_val) - (arbitrary_val * cy))) 
		            - (ay * ((bx * arbitrary_val) - (arbitrary_val * cx))) 
		            + (arbitrary_val * ((bx * cy) - (by * cx)));
	}

	output1 = 1;
	output2 = 2;
	output3 = 3;

	/*if (determinant >= 0)
	{cout << "OPT-1" << endl;
		output1 = 1;
		output2 = 2;
		output3 = 3;
	}
	else // if (determinant < 0)
	{cout << "OPT-2" << endl;
		output1 = 1;
		output2 = 3;
		output3 = 2;
	}*/
}


void unit_vector_func(double &output1, double &output2, double &output3, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, int qmode)
{
	int i, j, k, sign;
	
	double v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, nx, ny, nz;

	if (qmode == 2) // ((qmode == 3) || (qmode == 4))
	{
		sign = -1;
	}
	else
	{
		sign = 1;
	}

	v1x = ax;
	v1y = ay;
	v1z = az;
	v2x = bx;
	v2y = by;
	v2z = bz;
	v3x = cx;
	v3y = cy;
	v3z = cz;
	
	ccw_vertices(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z);

	/*cout << "FUNCTION" << ":   " << v1x << ", " << v1y << ", " << v1z << ";   ";
	cout << v2x << ", " << v2y << ", " << v2z << ";   ";
	cout << v3x << ", " << v3y << ", " << v3z << ";   ";*/

	double abx = v2x - v1x;
	double aby = v2y - v1y;
	double abz = v2z - v1z;

	double acx = v3x - v1x;
	double acy = v3y - v1y;
	double acz = v3z - v1z;

	double ux = (aby * acz) - (abz * acy);
	double uy = (abx * acz) - (abz * acx);
	double uz = (abx * acy) - (aby * acx);
	double u = pow(pow(ux, 2) + pow(uy, 2) + pow(uz, 2), 0.5);

	nx = sign * (ux / u);
	ny = sign * (uy / u);
	nz = sign * (uz / u);

	if (abs(nx) == 0)
	{
		nx = abs(nx);
	}
	if (abs(ny) == 0)
	{
		ny = abs(ny);
	}
	if (abs(nz) == 0)
	{
		nz = abs(nz);
	}

	cout << nx << ", " << ny << ", " << nz << endl;

	output1 = nx;
	output2 = ny;
	output3 = nz;
}


void ccw_vertices(double &ax, double &ay, double &az, double &bx, double &by, double &bz, double &cx, double &cy, double &cz)
{
	int i, j, k;
	double vectors[3][3];

	vectors[0][0] = ax;
	vectors[0][1] = ay;
	vectors[0][2] = az;
	vectors[1][0] = bx;
	vectors[1][1] = by;
	vectors[1][2] = bz;
	vectors[2][0] = cx;
	vectors[2][1] = cy;
	vectors[2][2] = cz;

	ccw_func(i, j, k, vectors[0][0], vectors[0][1], vectors[0][2], vectors[1][0], vectors[1][1], 
				vectors[1][2], vectors[2][0], vectors[2][1], vectors[2][2]);

	ax = vectors[i - 1][0];
	ay = vectors[i - 1][1];
	az = vectors[i - 1][2];
	
	bx = vectors[j - 1][0];
	by = vectors[j - 1][1];
	bz = vectors[j - 1][2];

	cx = vectors[k - 1][0];
	cy = vectors[k - 1][1];
	cz = vectors[k - 1][2];
}


void fill_vector(double (&output)[9], double val1, double val2, double val3, double val4, double val5, double val6, double val7, double val8, double val9)
{
	output[0] = val1;
	output[1] = val2;
	output[2] = val3;
	output[3] = val4;
	output[4] = val5;
	output[5] = val6;
	output[6] = val7;
	output[7] = val8;
	output[8] = val9;
}


/* ---------------------------------------------------------------------------------- */



string start_stl(string name)
{
	string solidName = charAndUnderscore(name);
	return ("solid " + solidName + "\n\n");
}


string end_stl(string name)
{
	return ("\n\nendsolid " + name + "\n");
}


/*
	Input-1	: X-Gap	[REF]
	Input-2	: Z-Gap [REF]
	Input-3	: Main Chord Length
	Input-4	: Relative Chord Length
	Input-5	: Thickness by Percentage of Main Chord
*/
void stl_gaps(double &input1, double &input2, double input3, double input4, double input5)
{
	double t = input5 / 100.;
	int relative_chord_factor;
	if ((int)(input3 / input4) == 1)
	{
		relative_chord_factor = 0;
	}
	else
	{
		relative_chord_factor = 1;
	}

	input1 = (0.2 * input3) + ((input3 * relative_chord_factor) / 2) - ((input4 * relative_chord_factor) / 2);
	input2 = (0.1 + (t / 2)) * input3;
}


/*
	Input-1	: Centroid Normalised X		[REF]
	Input-2	: Centroid X Vertex		[REF]
	Input-3	: Centroid Z Vertex		[REF]
	Input-4	: Thickness by Percentage of Main Chord
	Input-5	: Main Chord Length
	Input-6	: Relative Chord Length
	Input-7	: X-Gap
	Input-8	: Z-Gap
*/
void stl_preprocess(double **points, double (&input1)[ARRAY_NUM], double (&input2)[ARRAY_NUM], double (&input3)[ARRAY_NUM], double input4, double input5, double input6, double &input7, double &input8)
{
	stl_gaps(input7, input8, input5, input6, input4);

	double x_max, x_min, twisted_x_gap;
	double temp_x = 10;
	x_max = points[1][0];
	for (int i=1; i<((int)((points[0][0]) / 1.9)); i++)
	{
		if (points[i][0] > temp_x)
		{
			x_min = temp_x;
			break;
		}
		else
		{
			temp_x = points[i][0];
		}
	}

	for (int k=0; k<ARRAY_NUM; k++)
	{
		input2[k] = input7 + ((0.95 - (0.05 * k)) * (x_max - x_min) * input6);
		input1[k] = (0.95 - (0.05 * k)) * (x_max - x_min);

		double sol1, sol2;
		double temp_dx_old = 1;
		for (int j=1; j<((int)((points[0][0]) / 1.8)); j++)
		{
			double temp_dx_new = abs(points[j][0] - (0.95 - (0.05 * k)) * (x_max - x_min));
			int temp_array_num = (int)(points[0][0] - (j - 1));
			if (temp_dx_new <= temp_dx_old)
			{
				temp_dx_old = temp_dx_new;
			}
			else
			{
				sol1 = points[j][1];
				break;
			}
		}
		temp_dx_old = 1;
		for (int j2=1; j2<((int)((points[0][0]) / 1.8)); j2++)
		{
			double temp_dx_new = abs(points[(int)(points[0][0] - (j2 - 1))][0] - (0.95 - (0.05 * k)) * (x_max - x_min));
			if (temp_dx_new <= temp_dx_old)
			{
				temp_dx_old = temp_dx_new;
			}
			else
			{
				sol2 = points[(int)(points[0][0] - (j2 - 1))][1];
				break;
			}
		}
		input3[k] = input8 + (((sol1 + sol2) / 2) * input6);
	}
}


/*string profile_stl(double **pts, double y_pos, bool y_dir_positive, double t_main_percent, double t_by_c_percent, double main_chord, double relative_chord)
{
	string code = "";
	double centroid_normalised_x[ARRAY_NUM], centroid_vertex_x[ARRAY_NUM], centroid_vertex_z[ARRAY_NUM], x_gap, z_gap;
	stl_preprocess(pts,centroid_normalised_x, centroid_vertex_x, centroid_vertex_z, t_main_percent, main_chord, relative_chord, x_gap, z_gap);

	int dir = 1;
	if (!y_dir_positive)
	{
		dir = -1;
	}

	int centroid_num, centroid_num_prev;
	for (int i=1; i<((int)(pts[0][0])); i++)
	{
		if (i == 1)
		{
			int temp_array_num = (int)(pts[0][0]);
			code += stl_facet_code(0, dir * 1, 0, centroid_vertex_x[0], y_pos, centroid_vertex_z[0], 
				x_gap + (pts[temp_array_num][0] * relative_chord), y_pos, 
				z_gap + (pts[temp_array_num][1] * relative_chord), x_gap + (pts[1][0] * relative_chord), 
				y_pos, z_gap + (pts[1][1] * relative_chord));

			centroid_num_prev = centroid_array_num(pts[i+1][0], centroid_normalised_x, ARRAY_NUM);
		}


		centroid_num = centroid_array_num(pts[i+1][0], centroid_normalised_x, ARRAY_NUM);

		if (centroid_num != centroid_num_prev)
		{
			code += stl_facet_code(0, dir * 1, 0, centroid_vertex_x[centroid_num], y_pos, 
				centroid_vertex_z[centroid_num], centroid_vertex_x[centroid_num_prev], y_pos, 
				centroid_vertex_z[centroid_num_prev], x_gap + (pts[i][0] * relative_chord), y_pos, 
				z_gap + (pts[i][1] * relative_chord));
		}
		centroid_num_prev = centroid_num;


		code += stl_facet_code(0, dir * 1, 0, centroid_vertex_x[centroid_num], y_pos, 
			centroid_vertex_z[centroid_num], x_gap + (pts[i][0] * relative_chord), y_pos, 
			z_gap + (pts[i][1] * relative_chord), x_gap + (pts[i+1][0] * relative_chord), y_pos, 
			z_gap + (pts[i+1][1] * relative_chord));

	}

	return (code + "\n\n");
}*/


string profile_stl(double **pts, double y_pos, bool y_dir_positive, double t_main_percent, double t_by_c_percent, double main_chord, double relative_chord)
{

	string code = "";
	double centroid_normalised_x[ARRAY_NUM], centroid_vertex_x[ARRAY_NUM], centroid_vertex_z[ARRAY_NUM], x_gap, z_gap;
	stl_preprocess(pts,centroid_normalised_x, centroid_vertex_x, centroid_vertex_z, t_main_percent, main_chord, relative_chord, x_gap, z_gap);


	int dir = 1;
	if (!y_dir_positive)
	{
		dir = -1;
	}

	int centroid_num, centroid_num_prev;
	for (int i=1; i<((int)(pts[0][0])); i++)
	{
		if (i == 1)
		{
			int temp_array_num = (int)(pts[0][0]);
			if (!y_dir_positive)
			{
				code += stl_facet_code(0, -1, 0, centroid_vertex_x[0], y_pos, centroid_vertex_z[0], 
					x_gap + (pts[temp_array_num][0] * relative_chord), y_pos, 
					z_gap + (pts[temp_array_num][1] * relative_chord), x_gap + (pts[1][0] * relative_chord), 
					y_pos, z_gap + (pts[1][1] * relative_chord));
			}
			else
			{
				code += stl_facet_code(0, 1, 0, centroid_vertex_x[0], y_pos, centroid_vertex_z[0], 
					x_gap + (pts[1][0] * relative_chord), y_pos,  
					z_gap + (pts[1][1] * relative_chord), x_gap + (pts[temp_array_num][0] * relative_chord), 
					y_pos, z_gap + (pts[temp_array_num][1] * relative_chord));
			}

			centroid_num_prev = centroid_array_num(pts[i+1][0], centroid_normalised_x, ARRAY_NUM);
		}


		centroid_num = centroid_array_num(pts[i+1][0], centroid_normalised_x, ARRAY_NUM);

		if (centroid_num != centroid_num_prev)
		{
			if (!y_dir_positive)
			{
				code += stl_facet_code(0, -1, 0, centroid_vertex_x[centroid_num], y_pos, 
					centroid_vertex_z[centroid_num], centroid_vertex_x[centroid_num_prev], y_pos, 
					centroid_vertex_z[centroid_num_prev], x_gap + (pts[i][0] * relative_chord), y_pos, 
					z_gap + (pts[i][1] * relative_chord));
			}
			else
			{
				code += stl_facet_code(0, 1, 0, centroid_vertex_x[centroid_num], y_pos, 
					centroid_vertex_z[centroid_num], x_gap + (pts[i][0] * relative_chord), y_pos, 
					z_gap + (pts[i][1] * relative_chord), centroid_vertex_x[centroid_num_prev], y_pos, 
					centroid_vertex_z[centroid_num_prev]);
			}
		}
		centroid_num_prev = centroid_num;


		if (!y_dir_positive)
		{
			code += stl_facet_code(0, -1, 0, centroid_vertex_x[centroid_num], y_pos, 
				centroid_vertex_z[centroid_num], x_gap + (pts[i][0] * relative_chord), y_pos, 
				z_gap + (pts[i][1] * relative_chord), x_gap + (pts[i+1][0] * relative_chord), y_pos, 
				z_gap + (pts[i+1][1] * relative_chord));
		}
		else
		{
			code += stl_facet_code(0, 1, 0, centroid_vertex_x[centroid_num], y_pos, 
				centroid_vertex_z[centroid_num], x_gap + (pts[i+1][0] * relative_chord), y_pos, 
				z_gap + (pts[i+1][1] * relative_chord), x_gap + (pts[i][0] * relative_chord), y_pos, 
				z_gap + (pts[i][1] * relative_chord));
		}

	}

	return (code + "\n\n");
}


/*
	Input-1a: X-Gap of Previous Airfoil
	Input-1b: Z-Gap of Previous Airfoil
	Input-2a: X-Gap of Current Airfoil
	Input-2b: Z-Gap of Current Airfoils
	Input-3	: Main Chord Length
	Input-4	: Relative Chord Length of Previous Airfoil
	Input-5	: Relative Chord Length of Current Airfoil
	Input-6	: Thickness by Percentage of Root Airfoil
*/
string span_stl(double **pts1, double **pts2, double y1, double y2, double input3, double input4, double input5, double input6)
{
	string code = "";
	double input1a, input1b, input2a, input2b;
	stl_gaps(input1a, input1b, input3, input4, input6);
	stl_gaps(input2a, input2b, input3, input5, input6);
	
	
	int centroid_num, centroid_num_prev;
	for (int i=1; i<((int)(pts1[0][0])); i++)
	{
		// double theta = atan(abs(pts[i+1][1] - pts[i][1]) / abs(pts[i+1][0] - pts[i][0])) * 180 / PI;
		int dirx, diry;
		double theta_normalised;
		double dx1 = pts1[i+1][0] - pts1[i][0];
		double dy1 = pts1[i+1][1] - pts1[i][1];
		double dx2 = pts2[i+1][0] - pts2[i][0];
		double dy2 = pts2[i+1][1] - pts2[i][1];
		double r1 = pow(pow(dx1, 2) + pow(dy1, 2), 0.5);
		double r2 = pow(pow(dx2, 2) + pow(dy2, 2), 0.5);
		double r = (r1 + r2) / 2;

		// double theta1 = absolute_angle(atan(dy1 / dx1) * 180 / PI, dx1, dy1);
		double theta1 = cartesian_angle(pts1[i][0], pts1[i][1], pts1[i+1][0], pts1[i+1][1]);
		/*double theta1 = atan(abs(dy1) / abs(dx1)) * 180 / PI;
		if (dy1 < 0)
		{
			theta1 = -theta1;
		}*/

		// double theta2 = absolute_angle(atan(dy2 / dx2) * 180 / PI, dx2, dy2);
		double theta2 = cartesian_angle(pts2[i][0], pts2[i][1], pts2[i+1][0], pts2[i+1][1]);
		/*double theta2 = atan(abs(dy2) / abs(dx2)) * 180 / PI;
		if (dy2 < 0)
		{
			theta2 = -theta2;
		}*/

		double theta = (theta1 + theta2) / 2;
		/*double theta = (theta1 + theta2) / 2;
		if (theta < 0)
		{
			theta = 360 - abs(theta);
		}*/

		theta_normalised = normalise_angle(theta);

		//cout << i << ".  " << theta1 << " - " << theta2 << " - " << theta << " - ";
		//cout << theta_normalised << ";   ";
		//cout << dx1 << " - " << dy1 << ";   " << pts2[i+1][0] << " - " << pts2[i][0] << " - ";
		//cout << pts2[i+1][1] << " - " << pts2[i][1] << endl;

		double unit_vector_x, unit_vector_y, unit_vector_z;
		/*double prev_
		double prev_i_x*/

		/*double unit_vector_x = stod(decimal_point(to_string( dir_val(dir_val_angle(theta, 1) * r * cos(theta_normalised * PI / 180), 
					dir_val_angle(theta, 2) * r * sin(theta_normalised * PI / 180), 1) 
					* (pow(2,0.5) * cos(theta_normalised * PI / 180)) ), MAX_DEC_PTS));
		double X_VAR = abs(((pts1[i+1][0] + pts2[i+1][0]) / 2) - ((pts1[i][0] + pts2[i][0]) / 2));
		double Y_VAR = abs(y2 - y1);
		double Z_VAR = abs(((pts1[i+1][1] + pts2[i+1][1]) / 2) - ((pts1[i][1] + pts2[i][1]) / 2));
		double V_VAR = pow(pow(X_VAR, 2) + pow(Y_VAR, 2) + pow(Z_VAR, 2), 0.5);
		double unit_vector_x2 = X_VAR / V_VAR;
		
		double unit_vector_z = stod(decimal_point(to_string( dir_val(dir_val_angle(theta, 1) * r * cos(theta_normalised * PI / 180), 
					dir_val_angle(theta, 2) * r * sin(theta_normalised * PI / 180), 2) 
					* (pow(2,0.5) * sin(theta_normalised * PI / 180)) ), MAX_DEC_PTS));*/

		//int quad_mode = quadrant_mode(dir_val_angle(theta, 1) * r * cos(theta_normalised * PI / 180), 
							//dir_val_angle(theta, 2) * r * sin(theta_normalised * PI / 180));
		int quad_mode;
		if ((((pts1[i+1][0] + pts2[i+1][0]) / 2) - ((pts1[i][0] + pts2[i][0]) / 2)) < 0)
		{
			quad_mode = 1;
		}
		else
		{
			quad_mode = 2;
		}
		cout << "Q=" << quad_mode << " - " << input1a << ", " << input1b << ", " << input2a << ", " << input2b << endl;

		double vector[9];

		if (i == 1)
		{
			unit_vector_func(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i][0] * input4), y1, 
					input1b + (pts1[i][1] * input4), input2a + (pts2[i][0] * input5), y2, 
					input2b + (pts2[i][1] * input5), input1a + (pts1[i+1][0] * input4), y1, 
					input1b + (pts1[i+1][1] * input4), quad_mode);

			fill_vector(vector, input1a + (pts1[i][0] * input4), y1, input1b + (pts1[i][1] * input4), input2a + (pts2[i][0] * input5), 
				y2, input2b + (pts2[i][1] * input5), input1a + (pts1[i+1][0] * input4), y1, input1b + (pts1[i+1][1] * input4));

			//ccw_vertices(vector[0], vector[1], vector[2], vector[3], vector[4], vector[5], vector[6], vector[7], vector[8]);

			code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, vector[0], vector[1], vector[2], vector[3], vector[4], 
				vector[5], vector[6], vector[7], vector[8]);

			/*code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i][0] * input4), y1, 
				input1b + (pts1[i][1] * input4), input2a + (pts2[i][0] * input5), y2, 
				input2b + (pts2[i][1] * input5), input1a + (pts1[i+1][0] * input4), y1, 
				input1b + (pts1[i+1][1] * input4));*/

		}
		else
		{
			unit_vector_func(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i][0] * input4), y1, 
					input1b + (pts1[i][1] * input4), input2a + (pts2[i-1][0] * input5), y2, 
					input2b + (pts2[i-1][1] * input5), input2a + (pts2[i][0] * input5), y2, 
					input2b + (pts2[i][1] * input5), 1);

			fill_vector(vector, input1a + (pts1[i][0] * input4), y1, input1b + (pts1[i][1] * input4), input2a + (pts2[i-1][0] * input5), 
				y2, input2b + (pts2[i-1][1] * input5), input2a + (pts2[i][0] * input5), y2, input2b + (pts2[i][1] * input5));

			//ccw_vertices(vector[0], vector[1], vector[2], vector[3], vector[4], vector[5], vector[6], vector[7], vector[8]);

			code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, vector[0], vector[1], vector[2], vector[3], vector[4], 
				vector[5], vector[6], vector[7], vector[8]);

			/*code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i][0] * input4), y1, 
				input1b + (pts1[i][1] * input4), input2a + (pts2[i-1][0] * input5), y2, 
				input2b + (pts2[i-1][1] * input5), input2a + (pts2[i][0] * input5), y2, 
				input2b + (pts2[i][1] * input5));*/

			//cout << input1a + (pts1[i][0] * input4) << ", " << y1 << ", " << input1b + (pts1[i][1] * input4) << "; ";
			//cout << input2a + (pts2[i-1][0] * input5) << ", " << y2 << ", " << input2b + (pts2[i-1][1] * input5) << "; ";
			//cout << input2a + (pts2[i][0] * input5) << ", " << y2 << ", " << input2b + (pts2[i][1] * input5) << "; ";
			cout << theta_normalised << ":   " << vector[0] << ", " << vector[1] << ", " << vector[2] << ";   ";
			cout << vector[3] << ", " << vector[4] << ", " << vector[5] << ";   ";
			cout << vector[6] << ", " << vector[7] << ", " << vector[8] << ";   ";
			cout << unit_vector_x << ", " << unit_vector_y << ", " << unit_vector_z << endl;

			unit_vector_func(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i][0] * input4), y1, 
					input1b + (pts1[i][1] * input4), input2a + (pts2[i][0] * input5), y2, 
					input2b + (pts2[i][1] * input5), input1a + (pts1[i+1][0] * input4), y1, 
					input1b + (pts1[i+1][1] * input4), 1);

			fill_vector(vector, input1a + (pts1[i][0] * input4), y1, input1b + (pts1[i][1] * input4), input2a + (pts2[i][0] * input5), 
				y2, input2b + (pts2[i][1] * input5), input1a + (pts1[i+1][0] * input4), y1, input1b + (pts1[i+1][1] * input4));

			//ccw_vertices(vector[0], vector[1], vector[2], vector[3], vector[4], vector[5], vector[6], vector[7], vector[8]);

			code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, vector[0], vector[1], vector[2], vector[3], vector[4], 
				vector[5], vector[6], vector[7], vector[8]);

			/*code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i][0] * input4), y1, 
				input1b + (pts1[i][1] * input4), input2a + (pts2[i][0] * input5), y2, 
				input2b + (pts2[i][1] * input5), input1a + (pts1[i+1][0] * input4), y1, 
				input1b + (pts1[i+1][1] * input4));*/

			//cout << input1a + (pts1[i][0] * input4) << ", " << y1 << ", " << input1b + (pts1[i][1] * input4) << "; ";
			//cout << input2a + (pts2[i][0] * input5) << ", " << y2 << ", " << input2b + (pts2[i][1] * input5) << "; ";
			//cout << input1a + (pts1[i+1][0] * input4) << ", " << y1 << ", " << input1b + (pts1[i+1][1] * input4) << "; ";
			cout << theta_normalised << ":   " << vector[0] << ", " << vector[1] << ", " << vector[2] << ";   ";
			cout << vector[3] << ", " << vector[4] << ", " << vector[5] << ";   ";
			cout << vector[6] << ", " << vector[7] << ", " << vector[8] << ";   ";
			cout << unit_vector_x << ", " << unit_vector_y << ", " << unit_vector_z << endl;

		}

		if (i == ((int)(pts1[0][0] - 1)))
		{
			unit_vector_func(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i+1][0] * input4), y1, 
					input1b + (pts1[i+1][1] * input4), input2a + (pts2[i][0] * input5), y2, 
					input2b + (pts2[i][1] * input5), input2a + (pts2[i+1][0] * input5), y2, 
					input2b + (pts2[i+1][1] * input5), quad_mode);

			code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i+1][0] * input4), y1, 
				input1b + (pts1[i+1][1] * input4), input2a + (pts2[i][0] * input5), y2, 
				input2b + (pts2[i][1] * input5), input2a + (pts2[i+1][0] * input5), y2, 
				input2b + (pts2[i+1][1] * input5));

			cout << unit_vector_x << ", " << unit_vector_y << ", " << unit_vector_z << endl;

			unit_vector_func(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i+1][0] * input4), y1, 
					input1b + (pts1[i+1][1] * input4), input2a + (pts2[1][0] * input5), y2, 
					input2b + (pts2[1][1] * input5), input2a + (pts2[i+1][0] * input5), y2, 
					input2b + (pts2[i+1][1] * input5), quad_mode);

			code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i+1][0] * input4), y1, 
				input1b + (pts1[i+1][1] * input4), input2a + (pts2[1][0] * input5), y2, 
				input2b + (pts2[1][1] * input5), input2a + (pts2[i+1][0] * input5), y2, 
				input2b + (pts2[i+1][1] * input5));

			cout << unit_vector_x << ", " << unit_vector_y << ", " << unit_vector_z << endl;

			unit_vector_func(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i+1][0] * input4), y1, 
					input1b + (pts1[i+1][1] * input4), input1a + (pts1[1][0] * input4), y1, 
					input1b + (pts1[1][1] * input4), input2a + (pts2[1][0] * input5), y2, 
					input2b + (pts2[1][1] * input5), quad_mode);

			code += stl_facet_code(unit_vector_x, unit_vector_y, unit_vector_z, input1a + (pts1[i+1][0] * input4), y1, 
				input1b + (pts1[i+1][1] * input4), input1a + (pts1[1][0] * input4), y1, 
				input1b + (pts1[1][1] * input4), input2a + (pts2[1][0] * input5), y2, 
				input2b + (pts2[1][1] * input5));

			cout << unit_vector_x << ", " << unit_vector_y << ", " << unit_vector_z << endl;
		}

		//cout << theta << " - " << theta_normalised << " --- " << unit_vector_x << " - " << unit_vector_y << " - " << unit_vector_z << endl;
	}

	return (code + "\n\n");
}





