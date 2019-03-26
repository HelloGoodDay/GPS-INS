
#include "global_variables.h"

double degree2arc(double degree)
{
	return (degree*PI/180);
}
double arc2degree(double arc)
{
	double degree = arc/PI*180;
	if(degree > 180)
		degree -= 360;
	if(degree < -180)
		degree += 360;
	return degree;
}