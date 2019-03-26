#ifndef _global_variables_H
#define _global_variables_H

#define DATASUM 2000000
#define SOD     86400

#include "matrix.cc"
#include "iomanip"
#include "cmatrix"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include "coordinate_system_difine.h"
#include "coordinate_transformation.h"
using namespace std;

// ins parameters
#define PI 3.1415926535898
#define LIGHT_V  2.99792458E+08
#define EARTH_V  7.2921151467E-5
#define EARTH_a  6378137.
#define EARTH_b  6356752.3142
//#define EARTH_f  1 - EARTH_b / EARTH_a  
#define EARTH_f  1.0/298.257223563
#define EARTH_e (sqrt(EARTH_a*EARTH_a - EARTH_b*EARTH_b)/EARTH_a)
#define GM 3.986004418e+14
#define F_WGS84 1/298.257224

// GNSS parameters
#define EARTH_R 6378137.
#define GPS_NUM 32
#define GPS_L1 1.57542E+09
#define GPS_L2 1.2276E+09
#define LAMDA_Lw (LIGHT_V/(GPS_L1 - GPS_L2))
#define ELEV_MASK 20.0

#define _priorSigma 0.6
#define BAD_DCB_FLAG 70.0
#define Residual_Flag 5.0
#define BAD_P4_FLAG 150.0
#define ARC_LENGTH 40

// matrix control parameter
extern double MQk;
extern double MPk;
extern double MRk;

// extern parameters
extern double INS_INTERVAL;
extern double INTERVALS;
typedef techsoft::matrix<double>  dMatrix;
extern int para_sum;
extern dMatrix F;
extern dMatrix H;
extern dMatrix P0;
extern dMatrix DX0;
extern dMatrix G;
extern dMatrix Q;
extern dMatrix R;
extern dMatrix X;
extern dMatrix Z;

//transform degree to arc 
double degree2arc(double degree);
//degree [-180, 180]
double arc2degree(double arc);
//degree [0, 360]
double arc2degree2(double arc);

// matrix inv
void lubksb(const double *A, int n, const int *indx, double *b);
int matinv(double *A, int n);
int matinv(double *A, int n);
int ludcmp(double *A, int n, int *indx, double *d);
void matcpy(double *A, const double *B, int n, int m);
double *mat(int n, int m);
int *imat(int n, int m);
void matinv(dMatrix &A, int n);

#endif