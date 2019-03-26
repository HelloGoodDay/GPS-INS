#pragma once
#ifndef VEC_H
#define VEC_H

#include "global_variables.h"

class Vec
{
private:
	int len;
	double vec[4];

public:
	Vec(int len0 = 3)
	{
		if(len > 4)
		{
			cout<<"length is beyond limit!"<<endl;
		}
		len = len0;
		for (int i = 0;i<len;i++)
			vec[i] = 0.0;
	}

	int length();
	void display();
	double &operator[](int n);
	
	// if vector in the left has values, need to delete it before operator
	//Vec &operator = (const Vec &a);
	// addition
	Vec operator+ (const Vec &a);
	// subtraction
	Vec operator- (const Vec a);
	// cross
	Vec operator* (const Vec a);
	// multiply, attention: pre-multiply
	Vec operator* (double a);
	// division
	Vec operator/ (double a);
};

// multiply, friend function 
Vec operator*(double a, Vec &b);

#endif