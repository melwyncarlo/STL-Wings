#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include "MISC.H"

#ifndef AIRFOILS_H
#define AIRFOILS_H

using namespace std;

#define PI 3.14159265
const int N = 5;


string nacaCode(double m_by_c_percent, double p_by_c_percent, double t_by_c_percent);
double mu(double m, double x);
double mu_dash(double m, double x, int px);
double mu_left(double m, double p, double x);
double mu_right(double m, double p, double x);
double mu_dash_left(double m, double p, double x);
double mu_dash_right(double m, double p, double x);
double naca_mu_dash_left(double m, double p, double x);
double naca_mu_dash_right(double m, double p, double x);
double naca_symmetric_airfoil_function(double t, double x);
void naca_mu_functions(double &input1, double &input2, double input3, double input4, double input5);


/* ---------------------------------------------------------------------------------- */


double** flatPlate(double t_by_c_percent, int nodes);
double** parabolicallyCamberedPlate(double m_by_c_percent, double t_by_c_percent, int nodes);
double** quadraticallyCamberedPlate(double m_by_c_percent, double p_by_c_percent, double t_by_c_percent, int nodes);
double** nacaAirfoil(double m_by_c_percent, double p_by_c_percent, double t_by_c_percent, int min_nodes, double leading_edge_node_factor);


#endif



