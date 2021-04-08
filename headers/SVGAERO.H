#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include "MISC.H"

#ifndef SVGAERO_H
#define SVGAERO_H

using namespace std;


string data_points_to_string(double **points, int size);
string svgPolyLine_code(double param1, double param2, double param3, double param4, double param5);
string svgText_code(double x, double y, bool bold, string info);


/* ---------------------------------------------------------------------------------- */


string create_svg_airfoil(double **points, int size, string name, double t, double m, double p, double pixel, int mode);


#endif



