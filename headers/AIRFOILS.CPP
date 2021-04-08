#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include "MISC.H"
#include "AIRFOILS.H"

using namespace std;

#define PI 3.14159265


string nacaCode(double m_by_c_percent, double p_by_c_percent, double t_by_c_percent)
{
	string temp_var;
	if (to_string((int)(t_by_c_percent)).length() < 2)
	{
		temp_var = "0" + to_string((int)t_by_c_percent);
	}
	else
	{
		temp_var = to_string((int)t_by_c_percent);
	}
	return to_string((int)m_by_c_percent) + to_string((int)(p_by_c_percent / 10)) + temp_var;
}


double mu(double m, double x)
{
	return (m * (1 - pow(((2 * (x - 0.5)) / 1),2)));
}


double mu_dash(double m, double x, int px)
{
	if (px == 0)
	{
		return ((4 * m * (x + 0.5)) / pow(1,2));
	}
	else
	{
		return ((-4 * m * (x - 0.5)) / pow(1,2));
	}
}


double mu_left(double m, double p, double x)
{
	if ((m == 0) || (p == 0))
	{
		return 0;
	}
	else
	{
		return (m / (pow(p,2))) * ((2 * p * (x / 1)) - (pow((x / 1),2)));
	}
}


double mu_right(double m, double p, double x)
{
	if ((m == 0) || (p == 0))
	{
		return 0;
	}
	else
	{
		return (m / (pow((1 - p),2))) * ((1 - (2 * p)) + (2 * p * (x / 1)) - (pow((x / 1),2)));
	}
}


double mu_dash_left(double m, double p, double x)
{
	if ((m == 0) || (p == 0))
	{
		return 0;
	}
	else
	{
		return ((2 * m * x) / (p * 1));
	}
}


double mu_dash_right(double m, double p, double x)
{
	if ((m == 0) || (p == 0))
	{
		return 0;
	}
	else
	{
		return (((2 * m) * ((p * 1) - 1) * (x - 1)) / (pow(1,2) * pow((1 - p),2)));
	}
}


double naca_mu_dash_left(double m, double p, double x)
{
	if ((m == 0) || (p == 0))
	{
		return 0;
	}
	else
	{
		return (((2 * m) / (pow(p,2))) * (p - (x / 1)));
	}
}


double naca_mu_dash_right(double m, double p, double x)
{
	if ((m == 0) || (p == 0))
	{
		return 0;
	}
	else
	{
		return (((2 * m) / (pow((1 - p),2))) * (p - (x / 1)));
	}
}


double naca_symmetric_airfoil_function(double t, double x)
{
	return (5 * t * ((0.2969 * pow(x,0.5)) - (0.1260 * x) - (0.3516 * pow(x,2)) + (0.2843 * pow(x,3)) - (0.1015 * pow(x,4))));
}


/*
	Input-1	: Theta Variable	[REF]
	Input-2	: PYC Variable		[REF]
	Input-3	: Maximum Camber Height Decimal Value
	Input-4	: Maximum Camber Position Decimal Value
	Input-5	: X-Position
*/
void naca_mu_functions(double &input1, double &input2, double input3, double input4, double input5)
{
	if (input5 >= input4)
	{
		input1 = atan(naca_mu_dash_right(input3, input4, input5));
		input2 = mu_right(input3, input4, input5);
	}
	else
	{
		input1 = atan(naca_mu_dash_left(input3, input4, input5));
		input2 = mu_left(input3, input4, input5);
	}
}


/* ---------------------------------------------------------------------------------- */


double** flatPlate(double t_by_c_percent, int nodes)
{
	int n_nodes = (int)(1.5 * nodes);
	double **array = new double*[n_nodes];
	for (int k=0; k<n_nodes; k++)
	{
		array[k] = new double[N];
	}
	double t = t_by_c_percent / 100.;
	int remainder_nodes = (int)(nodes - 7);
	int half_nodes = (int)((remainder_nodes / 2) + 1);
	int segments = (int)(half_nodes + 1);
	double segment_width = 1. / segments;
	int points_array_size;

	array[0][2] = 2;
	array[0][3] = 2;
	array[1][2] = 1;
	array[1][3] = 0;
	array[2][2] = 0;
	array[2][3] = 0;
	array[1][0] = 1;
	array[1][1] = 0;
	array[2][0] = 1;
	array[2][1] = t / 2;
	points_array_size = 2;
	
	for (int i=1; i<=half_nodes; i++)
	{
		array[i+2][0] = 1 - (i * segment_width);
		array[i+2][1] = t / 2;
		points_array_size++;
	}
	
	array[points_array_size+1][0] = 0;
	array[points_array_size+1][1] = t / 2;
	array[points_array_size+2][0] = 0;
	array[points_array_size+2][1] = 0;
	array[points_array_size+3][0] = 0;
	array[points_array_size+3][1] = -t / 2;
	points_array_size += 3;

	int temp_var = points_array_size;

	for (int j=1; j<=half_nodes; j++)
	{
		array[j+temp_var][0] = 0 + (j * segment_width);
		array[j+temp_var][1] = -t / 2;
		points_array_size++;
	}

	array[points_array_size+1][0] = 1;
	array[points_array_size+1][1] = -t / 2;
	array[points_array_size+2][0] = 1;
	array[points_array_size+2][1] = 0;
	points_array_size += 2;
	array[0][0] = points_array_size;
	array[0][1] = points_array_size;
	
	return array;
}


double** parabolicallyCamberedPlate(double m_by_c_percent, double t_by_c_percent, int nodes)
{
	int n_nodes = (int)(1.5 * nodes);
	double **array = new double*[n_nodes];
	for (int n=0; n<n_nodes; n++)
	{
		array[n] = new double[N];
	}
	double m = m_by_c_percent / 100.;
	double t = t_by_c_percent / 100.;
	int remainder_nodes = nodes - 2;
	int half_nodes = (int)((remainder_nodes / 2) + 1);
	int quarter_nodes = (int)((half_nodes / 2) + 1);
	int segments = quarter_nodes;
	double segment_width = 0.5 / segments;
	int points_array_size = 0;
	int line_array_size = 0;
	double prx, pry, plx, ply, theta, gamma, dx, dy;
	double px_i, py_i, px_f, py_f;

	int i;
	for (i=0; i<quarter_nodes; i++)
	{
		prx = 1 - (i * segment_width);
		plx = 1 - ((i + 1) * segment_width);
		pry = mu(m, prx);
		ply = mu(m, plx);
		if (i == 0)
		{
			px_f = prx;
			py_f = pry;
			theta = atan(abs(mu_dash(m, prx, 1) - mu_dash(m, plx, 1)) / abs(prx - plx)) * 180 / PI;
		}
		else
		{
			theta = atan(abs(pry - ply) / abs(prx - plx)) * 180 / PI;
		}
		gamma = 90. - theta;
		dx = (t / 2.) * cos(gamma * PI / 180);
		dy = (t / 2.) * sin(gamma * PI / 180);
		array[i+1][0] = prx + dx;
		array[i+1][1] = pry + dy;
		array[i+1][2] = prx;
		array[i+1][3] = pry;
		array[i+1][4] = theta;
		points_array_size++;
		line_array_size++;
	}

	array[points_array_size+1][0] = 0.5;
	array[points_array_size+1][1] = m + (t / 2);
	array[line_array_size+1][2] = 0.5;
	array[line_array_size+1][3] = m;
	array[points_array_size+1][4] = 0;
	points_array_size++;
	line_array_size++;
	int temp_var_05 = points_array_size;
	int temp_var_06 = line_array_size;
	int temp_var_07 = 1;

	for (int j=(quarter_nodes-1); j>0; j--)
	{
		prx = 0 + (j * segment_width);
		plx = 0 + ((j - 1) * segment_width);
		pry = mu(m, prx);
		ply = mu(m, plx);
		theta = atan(abs(pry - ply) / abs(prx - plx)) * 180 / PI;
		gamma = 90 - theta;
		dx = (t / 2) * cos(gamma * PI / 180);
		dy = (t / 2) * sin(gamma * PI / 180);
		array[temp_var_05+temp_var_07][0] = prx - dx;
		array[temp_var_05+temp_var_07][1] = pry + dy;
		array[temp_var_06+temp_var_07][2] = prx;
		array[temp_var_06+temp_var_07][3] = pry;
		array[temp_var_05+temp_var_07][4] = theta;
		if (j == 1)
		{
			px_i = plx;
			py_i = ply;
			theta = atan(abs(mu_dash(m, prx, 0) - mu_dash(m, plx, 0)) / abs(prx - plx)) * 180 / PI;
			gamma = 90 - theta;
			dx = (t / 2) * cos(gamma * PI / 180);
			dy = (t / 2) * sin(gamma * PI / 180);
			array[temp_var_05+temp_var_07+1][0] = px_i - dx;
			array[temp_var_05+temp_var_07+1][1] = py_i + dy;
			array[temp_var_05+temp_var_07+1][4] = theta;
			points_array_size++;
		}
		temp_var_07++;
		points_array_size++;
		line_array_size++;
	}

	array[line_array_size+1][2] = 0;
	array[line_array_size+1][3] = 0;
	line_array_size++;
	array[0][2] = line_array_size;
	array[0][3] = line_array_size;
	temp_var_05 = points_array_size;
	temp_var_07 = 1;

	for (int k=0; k<quarter_nodes; k++)
	{
		prx = 0 + ((k + 1) * segment_width);
		plx = 0 + (k * segment_width);
		pry = mu(m, prx);
		ply = mu(m, plx);
		if (i == 0)
		{
			px_f = prx;
			py_f = pry;
			theta = atan(abs(mu_dash(m, prx, 0) - mu_dash(m, plx, 0)) / abs(prx - plx)) * 180 / PI;
		}
		else
		{
			theta = atan(abs(pry - ply) / abs(prx - plx)) * 180 / PI;
		}
		gamma = 90 - theta;
		dx = (t / 2) * cos(gamma * PI / 180);
		dy = (t / 2) * sin(gamma * PI / 180);
		array[temp_var_05+temp_var_07][0] = plx + dx;
		array[temp_var_05+temp_var_07][1] = ply - dy;
		array[temp_var_05+temp_var_07][4] = theta;
		temp_var_07++;
		points_array_size++;
	}

	array[points_array_size+1][0] = 0.5;
	array[points_array_size+1][1] = m - (t / 2);
	array[points_array_size+1][4] = theta;
	points_array_size++;
	temp_var_05 = points_array_size;
	temp_var_07 = 1;

	for (int l=(quarter_nodes-1); l>0; l--)
	{
		prx = 1 - ((l - 1) * segment_width);
		plx = 1 - (l * segment_width);
		pry = mu(m, prx);
		ply = mu(m, plx);
		theta = atan(abs(pry - ply) / abs(prx - plx)) * 180 / PI;
		gamma = 90 - theta;
		dx = (t / 2) * cos(gamma * PI / 180);
		dy = (t / 2) * sin(gamma * PI / 180);
		array[temp_var_05+temp_var_07][0] = plx - dx;
		array[temp_var_05+temp_var_07][1] = ply - dy;
		array[temp_var_05+temp_var_07][4] = theta;
		if (l == 1)
		{
			theta = atan(abs(mu_dash(m, prx, 1) - mu_dash(m, plx, 1)) / abs(prx - plx)) * 180 / PI;
			gamma = 90 - theta;
			dx = (t / 2) * cos(gamma * PI / 180);
			dy = (t / 2) * sin(gamma * PI / 180);
			array[temp_var_05+temp_var_07+1][0] = px_f - dx;
			array[temp_var_05+temp_var_07+1][1] = py_f - dy;
			array[temp_var_05+temp_var_07+1][4] = theta;
			points_array_size++;
		}
		temp_var_07++;
		points_array_size++;
	}

	array[0][0] = points_array_size;
	array[0][1] = points_array_size;

	return array;
}


double** quadraticallyCamberedPlate(double m_by_c_percent, double p_by_c_percent, double t_by_c_percent, int nodes)
{
	int n_nodes = (int)(1.5 * nodes);
	double **array = new double*[n_nodes];
	for (int n=0; n<n_nodes; n++)
	{
		array[n] = new double[N];
	}
	double m = m_by_c_percent / 100.;
	double p = p_by_c_percent / 100.;
	double t = t_by_c_percent / 100.;
	int remainder_nodes = nodes - 2;
	int half_nodes = (int)((remainder_nodes / 2) + 1);
	int quarter_nodes_left = (int)((half_nodes * p) + 1);
	int quarter_nodes_right = (int)((half_nodes * (1 - p)) + 1);
	int segments_left = quarter_nodes_left;
	int segments_right = quarter_nodes_right;
	double segment_width_left = p / segments_left;
	double segment_width_right = (1 - p) / segments_right;
	int points_array_size = 0;
	int line_array_size = 0;
	double temp_var_01, temp_var_02, temp_var_03, temp_var_04;
	double prx, pry, plx, ply, theta, gamma, dx, dy;
	double px_i, py_i, px_f, py_f;

	int i;
	for (i=0; i<quarter_nodes_right; i++)
	{
		prx = 1 - (i * segment_width_right);
		plx = 1 - ((i + 1) * segment_width_right);
		pry = mu_right(m, p, prx);
		ply = mu_right(m, p, plx);
		if (i == 0)
		{
			px_f = prx;
			py_f = pry;
			theta = atan(abs(mu_dash_right(m, p, prx) - mu_dash_right(m, p, plx)) / abs(prx - plx)) * 180 / PI;
		}
		else
		{
			theta = atan(abs(pry - ply) / abs(prx - plx)) * 180 / PI;
		}
		gamma = 90. - theta;
		dx = (t / 2.) * cos(gamma * PI / 180);
		dy = (t / 2.) * sin(gamma * PI / 180);
		array[i+1][0] = prx + dx;
		array[i+1][1] = pry + dy;
		array[i+1][2] = prx;
		array[i+1][3] = pry;
		array[i+1][4] = theta;
		points_array_size++;
		line_array_size++;
	}

	array[points_array_size+1][0] = p;
	array[points_array_size+1][1] = m + (t / 2);
	array[line_array_size+1][2] = p;
	array[line_array_size+1][3] = m;
	array[points_array_size+1][4] = theta;
	
	points_array_size++;
	line_array_size++;
	int temp_var_05 = points_array_size;
	int temp_var_06 = line_array_size;
	int temp_var_07 = 1;

	for (int j=(quarter_nodes_left-1); j>0; j--)
	{
		prx = 0 + (j * segment_width_left);
		plx = 0 + ((j - 1) * segment_width_left);
		pry = mu_left(m, p, prx);
		ply = mu_left(m, p, plx);
		theta = atan(abs(pry - ply) / abs(prx - plx)) * 180 / PI;
		gamma = 90 - theta;
		dx = (t / 2) * cos(gamma * PI / 180);
		dy = (t / 2) * sin(gamma * PI / 180);
		array[temp_var_05+temp_var_07][0] = prx - dx;
		array[temp_var_05+temp_var_07][1] = pry + dy;
		array[temp_var_06+temp_var_07][2] = prx;
		array[temp_var_06+temp_var_07][3] = pry;
		array[temp_var_05+temp_var_07][4] = theta;
		if (j == 1)
		{
			px_i = plx;
			py_i = ply;
			theta = atan(abs(mu_dash_left(m, p, prx) - mu_dash_left(m, p, plx)) / abs(prx - plx)) * 180 / PI;
			gamma = 90 - theta;
			dx = (t / 2) * cos(gamma * PI / 180);
			dy = (t / 2) * sin(gamma * PI / 180);
			array[temp_var_05+temp_var_07+1][0] = px_i - dx;
			array[temp_var_05+temp_var_07+1][1] = py_i + dy;
			array[temp_var_05+temp_var_07+1][4] = theta;
			points_array_size++;
		}
		temp_var_07++;
		points_array_size++;
		line_array_size++;
	}

	array[line_array_size+1][2] = 0;
	array[line_array_size+1][3] = 0;
	line_array_size++;
	array[0][2] = line_array_size;
	array[0][3] = line_array_size;
	temp_var_05 = points_array_size;
	temp_var_07 = 1;

	for (int k=0; k<quarter_nodes_left; k++)
	{
		prx = 0 + ((k + 1) * segment_width_left);
		plx = 0 + (k * segment_width_left);
		pry = mu_left(m, p, prx);
		ply = mu_left(m, p, plx);
		if (i == 0)
		{
			px_f = prx;
			py_f = pry;
			theta = atan(abs(mu_dash_left(m, p, prx) - mu_dash_left(m, p, plx)) / abs(prx - plx)) * 180 / PI;
		}
		else
		{
			theta = atan(abs(pry - ply) / abs(prx - plx)) * 180 / PI;
		}
		gamma = 90 - theta;
		dx = (t / 2) * cos(gamma * PI / 180);
		dy = (t / 2) * sin(gamma * PI / 180);
		array[temp_var_05+temp_var_07][0] = plx + dx;
		array[temp_var_05+temp_var_07][1] = ply - dy;
		array[temp_var_05+temp_var_07][4] = theta;
		temp_var_07++;
		points_array_size++;
	}

	array[points_array_size+1][0] = p;
	array[points_array_size+1][1] = m - (t / 2);
	array[points_array_size+1][4] = theta;
	points_array_size++;
	temp_var_05 = points_array_size;
	temp_var_07 = 1;

	for (int l=(quarter_nodes_right-1); l>0; l--)
	{
		prx = 1 - ((l - 1) * segment_width_right);
		plx = 1 - (l * segment_width_right);
		pry = mu_right(m, p, prx);
		ply = mu_right(m, p, plx);
		theta = atan(abs(pry - ply) / abs(prx - plx)) * 180 / PI;
		gamma = 90 - theta;
		dx = (t / 2) * cos(gamma * PI / 180);
		dy = (t / 2) * sin(gamma * PI / 180);
		array[temp_var_05+temp_var_07][0] = plx - dx;
		array[temp_var_05+temp_var_07][1] = ply - dy;
		array[temp_var_05+temp_var_07][4] = theta;
		if (l == 1)
		{
			theta = atan(abs(mu_dash_right(m, p, prx) - mu_dash_right(m, p, plx)) / abs(prx - plx)) * 180 / PI;
			gamma = 90 - theta;
			dx = (t / 2) * cos(gamma * PI / 180);
			dy = (t / 2) * sin(gamma * PI / 180);
			array[temp_var_05+temp_var_07+1][0] = px_f - dx;
			array[temp_var_05+temp_var_07+1][1] = py_f - dy;
			array[temp_var_05+temp_var_07+1][4] = theta;
			points_array_size++;
		}
		temp_var_07++;
		points_array_size++;
	}

	array[0][0] = points_array_size;
	array[0][1] = points_array_size;

	return array;
}


double** nacaAirfoil(double m_by_c_percent, double p_by_c_percent, double t_by_c_percent, int min_nodes, double leading_edge_node_factor)
{
	int n_nodes = (int)(min_nodes + (min_nodes * leading_edge_node_factor) + 100);
	double **array = new double*[n_nodes];
	for (int n=0; n<n_nodes; n++)
	{
		array[n] = new double[N];
	}
	double m = m_by_c_percent / 100.;
	double p = p_by_c_percent / 100.;
	double t = t_by_c_percent / 100.;
	if ((m == 0) || (p == 0))
	{
	    m = 0;
	    p = 0;
	}

	int remainder_nodes = min_nodes;
	int half_nodes = (int)(remainder_nodes / 2);
	int quarter_nodes_right = (int)(half_nodes / 2);
	int quarter_nodes_left_right = (int)(quarter_nodes_right / 2);
	int quarter_nodes_left_left = (int)((leading_edge_node_factor * min_nodes) / 2);
	int segments_left_right = quarter_nodes_left_right - 1;
	int segments_left_left = quarter_nodes_left_left - 1;
	int segments_right = quarter_nodes_right - 1;
	double segment_width_left_right = 0.25 / segments_left_right;
	double segment_width_left_left = 0.25 / segments_left_left;
	double segment_width_right = 0.5 / segments_right;
	int points_array_size = 0;
	int line_array_size = 0;
	double px, theta, pyc, pyt;

	for (int i=0; i<quarter_nodes_right; i++)
	{
		px = 1 - (i * segment_width_right);
		naca_mu_functions(theta, pyc, m, p, px);
		pyt = naca_symmetric_airfoil_function(t, px);
		array[i+1][0] = px - (pyt * sin(theta));
		array[i+1][1] = pyc + (pyt * cos(theta));
		array[i+1][2] = px;
		array[i+1][3] = pyc;
		array[i+1][4] = theta;
		points_array_size++;
		line_array_size++;
	}

	int temp_var_01 = points_array_size;
	int temp_var_02 = line_array_size;
	int temp_var_03 = 1;

	for (int j1=1; j1<quarter_nodes_left_right; j1++)
	{
		px = 0.5 - (j1 * segment_width_left_right);
		naca_mu_functions(theta, pyc, m, p, px);
		pyt = naca_symmetric_airfoil_function(t, px);
		array[temp_var_01+temp_var_03][0] = px - (pyt * sin(theta));
		array[temp_var_01+temp_var_03][1] = pyc + (pyt * cos(theta));
		array[temp_var_02+temp_var_03][2] = px;
		array[temp_var_02+temp_var_03][3] = pyc;
		array[temp_var_01+temp_var_03][4] = theta;
		temp_var_03++;
		points_array_size++;
		line_array_size++;
	}

	temp_var_01 = points_array_size;
	temp_var_02 = line_array_size;
	temp_var_03 = 1;

	for (int j2=1; j2<quarter_nodes_left_left; j2++)
	{
		px = 0.25 - (j2 * segment_width_left_left);
		naca_mu_functions(theta, pyc, m, p, px);
		pyt = naca_symmetric_airfoil_function(t, px);
		array[temp_var_01+temp_var_03][0] = px - (pyt * sin(theta));
		array[temp_var_01+temp_var_03][1] = pyc + (pyt * cos(theta));
		array[temp_var_02+temp_var_03][2] = px;
		array[temp_var_02+temp_var_03][3] = pyc;
		array[temp_var_01+temp_var_03][4] = theta;
		temp_var_03++;
		points_array_size++;
		line_array_size++;
	}

	array[0][2] = line_array_size;
	array[0][3] = line_array_size;
	temp_var_01 = points_array_size;
	temp_var_03 = 1;

	for (int k1=0; k1<quarter_nodes_left_left; k1++)
	{
		px = 0 + (k1 * segment_width_left_left);
		naca_mu_functions(theta, pyc, m, p, px);
		pyt = naca_symmetric_airfoil_function(t, px);
		array[temp_var_01+temp_var_03][0] = px + (pyt * sin(theta));
		array[temp_var_01+temp_var_03][1] = pyc - (pyt * cos(theta));
		array[temp_var_01+temp_var_03][4] = theta;
		temp_var_03++;
		points_array_size++;
	}

	temp_var_01 = points_array_size;
	temp_var_03 = 1;

	for (int k2=1; k2<quarter_nodes_left_right; k2++)
	{
		px = 0.25 + (k2 * segment_width_left_right);
		naca_mu_functions(theta, pyc, m, p, px);
		pyt = naca_symmetric_airfoil_function(t, px);
		array[temp_var_01+temp_var_03][0] = px + (pyt * sin(theta));
		array[temp_var_01+temp_var_03][1] = pyc - (pyt * cos(theta));
		array[temp_var_01+temp_var_03][4] = theta;
		temp_var_03++;
		points_array_size++;
	}

	temp_var_01 = points_array_size;
	temp_var_03 = 1;

	for (int l=1; l<quarter_nodes_right; l++)
	{
		px = 0.5 + (l * segment_width_right);
		naca_mu_functions(theta, pyc, m, p, px);
		pyt = naca_symmetric_airfoil_function(t, px);
		array[temp_var_01+temp_var_03][0] = px + (pyt * sin(theta));
		array[temp_var_01+temp_var_03][1] = pyc - (pyt * cos(theta));
		array[temp_var_01+temp_var_03][4] = theta;
		temp_var_03++;
		points_array_size++;
	}

	array[0][0] = points_array_size;
	array[0][1] = points_array_size;

	return array;
}




