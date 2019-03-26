#pragma once
#ifndef VEC_H
#define VEC_H

#include "global_variables.h"

void multiply(double *vec1, double a, double *vsum);
void cross(double *vec1, double *vec2, double *vsum);
void add(double *vec1, double *vec2, double *vsum);
void subtract(double *vec1, double *vec2, double *vsum);
void division(double *vec1, double a, double *vsum);

#endif