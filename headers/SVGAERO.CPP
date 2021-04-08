#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include "MISC.H"
#include "SVGAERO.H"

using namespace std;


string data_points_to_string(double **points, int size)
{
	string data = "";
	int i;
	for (i=1; i<=size; i++)
	{
		data += decimal_point(to_string(points[i][0]), 5) + " " + decimal_point(to_string(points[i][1]), 5) + "\n";
	}
	return data;
}


/*
	Param-1	: X1-Pos
	Param-2	: Y1-Pos
	Param-3	: X2-Pos
	Param-4	: Y2-Pos
	Param-5	: Stroke Width
*/
string svgPolyLine_code(double param1, double param2, double param3, double param4, double param5)
{
	string temp = "";
	temp += "<polyline points='" + to_string(param1) + "," + to_string(param2) + " ";
	temp += to_string(param3) + "," + to_string(param4) + "' style='fill:none;stroke:black;stroke-width:";
	temp += to_string(param5) + "' /> \n";
	return temp;
}


string svgText_code(double x, double y, bool bold, string info)
{
	string temp = "<text x='" + to_string(x) + "' y='" + to_string(y) + "' style='fill:black;' ";
	if (bold)
	{
		temp += "font-weight='bold'>";
	}
	else
	{
		temp += ">";
	}
	temp += info + " </text>";
	return temp;
}


/* ---------------------------------------------------------------------------------- */


string create_svg_airfoil(double **points, int size, string name, double t, double m, double p, double pixel, int mode)
{
	string temp_var_01, temp_var_02;
	if (m == 0)
	{
		temp_var_01 = "n/a";
	}
	else
	{
		temp_var_01 = decimal_point(to_string(m/100.),2) + " c";
	}
	if (p == 0)
	{
		temp_var_02 = "n/a";
	}
	else
	{
		temp_var_02 = decimal_point(to_string(p/100.),2) + " c";
	}
	string code = "";

	string code_start = "<?xml version='1.0' encoding='UTF-8'?> \n \n <!DOCTYPE svg PUBLIC ";
	code_start += "'-//W3C//DTD SVG 1.1//EN' 'http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd'> \n \n";
	code_start += "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='1100' height='1225'> \n \n";

	string type = "Simple Grid";
	if ((mode == 2) || (mode == 3))
	{
		type = "Fine Grid";
		for (int k=0; k<=20; k++)
		{
			if (mode == 3)
			{
				type = "Finer Grid";
				if (k != 0)
				{
					for (int l=1; l<4; l++)
					{
						code += svgPolyLine_code(50 + (k * (((pixel / 20) / 4) * l)), 0, 
							50 + (k * (((pixel / 20) / 4) * l)), pixel + 100, 0.25);

						code += svgPolyLine_code(25, 50 + (k * (((pixel / 20) / 4) * l)), pixel + 100,
							50 + (k * (((pixel / 20) / 4) * l)), 0.25);
					}
				}
			}
			if (k == 10)
			{
				continue;
			}
			else
			{
				code += svgPolyLine_code(50 + ((pixel / 20) * k), 0, 50 + ((pixel / 20) * k), pixel + 100, 
					(1 * abs((k % 2) - 1)) + (0.5 * (k % 2)));

				code += svgPolyLine_code(2, 50 + ((pixel / 20) * k), pixel + 100, 50 + ((pixel / 20) * k), 
					(0.5 * abs((k % 2) - 1)) + (1 * (k % 2)));
			}
		}		
	}

	code += "<polygon points='";
	int i;
	for (i=1; i<=size; i++)
	{
		double i2 = (50 + (points[i][0] * pixel));
		double i3 = (550 - (points[i][1] * pixel));
		code += to_string(i2) + "," + to_string(i3) + " ";
	}
	code += "' style='fill:silver;stroke:dimgrey;stroke-width:1' />";
	
	code += "<polygon points='0,0 " + to_string(pixel + 100) + ",0 " + to_string(pixel + 100) + "," + to_string(pixel + 225) + 
		" 0," + to_string(pixel + 225) + " 0,0' ";
	code += "style='fill:none;stroke:black;stroke-width:1' />";

	code += svgPolyLine_code(10,550,1100,550,1);
	code += svgPolyLine_code(25,0,25,1100,1);
	code += svgPolyLine_code(550,0,550,1100,1);
	code += svgPolyLine_code(400,1100,1100,1100,1);

	for (int m=0; m<=20; m++)
	{
		double w1, w2;
		w1 = 1;
		w2 = 1;
		if (m == 0)
		{
			w1 = 0;
			w2 = 2;
		}
		else if (m == 10)
		{
			w1 = 3;
			w2 = 3;
		}
		else if (m == 20)
		{
			w1 = 0;
			w2 = 2;
		}
		code += svgPolyLine_code(50 + ((pixel / 20) * m), (50 + (pixel / 2)) - 15, 50 + ((pixel / 20) * m), 
			(50 + (pixel / 2)) + 15, w1);

		code += svgPolyLine_code(10, 50 + ((pixel / 20) * m), 40, 50 + ((pixel / 20) * m), w2);
	}

	code += "<polygon points='0," + to_string(pixel + 100) + " 550," + to_string(pixel + 100) + " 550," + to_string(pixel + 225) + 
		" 0," + to_string(pixel + 225) + " 0,";
	code += to_string(pixel + 100) + "' style='fill:#FFFF99;stroke:black;stroke-width:1' />";

	code += svgText_code(50,(pixel + 100 + 30),true, "Airfoil Name");
	code += svgText_code(50,(pixel + 100 + 55),true, "Thickness");
	code += svgText_code(50,(pixel + 100 + 80),true, "Max. Camber");
	code += svgText_code(50,(pixel + 100 + 105),true, "Max. Camber Point");
	code += svgText_code(220,(pixel + 100 + 30),true, ":");
	code += svgText_code(220,(pixel + 100 + 55),true, ":");
	code += svgText_code(220,(pixel + 100 + 80),true, ":");
	code += svgText_code(220,(pixel + 100 + 105),true, ":");
	code += svgText_code(240,(pixel + 100 + 30),false, name);
	code += svgText_code(240,(pixel + 100 + 55),false, decimal_point(to_string(t/100.),2) + " c");
	code += svgText_code(240,(pixel + 100 + 80),false, temp_var_01);
	code += svgText_code(240,(pixel + 100 + 105),false, temp_var_02);

	code += svgText_code(650,(pixel + 100 + 30),true, "Type");
	code += svgText_code(650,(pixel + 100 + 55),true, "Resolution");
	code += svgText_code(650,(pixel + 100 + 80),true, "Total No. of Nodes");
	code += svgText_code(820,(pixel + 100 + 30),true, ":");
	code += svgText_code(820,(pixel + 100 + 55),true, ":");
	code += svgText_code(820,(pixel + 100 + 80),true, ":");
	code += svgText_code(840,(pixel + 100 + 30),false, type);
	code += svgText_code(840,(pixel + 100 + 55),false, to_string((int)pixel) + "px by " + to_string((int)pixel) + "px");
	code += svgText_code(840,(pixel + 100 + 80),false, to_string(size - 1));
	

	code += "<polyline points='";
	int j;
	int j2 = (int)(points[0][3]);
	for (j=1; j<=j2; j++)
	{
		double j3 = (50 + (points[j][2] * pixel));
		double j4 = (550 - (points[j][3] * pixel));
		code += to_string(j3) + "," + to_string(j4) + " ";
	}
	code += "' style='fill:none;stroke:dimgrey;stroke-width:1' />";

	code += "Sorry, your browser does not support inline SVG.";
	code += "\n \n </svg>";

	return (code_start + code);
}



