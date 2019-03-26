#ifndef _global_variables_H
#define _global_variables_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.1415926535898
#define EARTH_V  7.292115E-5
#define EARTH_a  6378137.
#define EARTH_b  6356752.3141
//椭球扁率
#define EARTH_f  1 - EARTH_b / EARTH_a  
//椭球第一偏心率
#define EARTH_e (sqrt(EARTH_a*EARTH_a - EARTH_b*EARTH_b)/EARTH_a)
#define GM 3.986005e+14

#define FILEMAX 256

//transform degree to arc
double degree2arc(double degree);
//degree [-180, 180]
double arc2degree(double arc);

#endif