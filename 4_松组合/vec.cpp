#include "vec.h"

int Vec::length()
{
	return len;
}

void Vec::display()
{
	for (int i=0; i<len; i++)
		printf("%12.4lf", vec[i]);
	printf("\n");
}

double &Vec::operator[](int n)
{
	if(n < 0 && n >= len)
	{
		printf("Index is greater than upper limit!\n");
		return vec[0];
	}
	return vec[n]; 
}
/*
Vec &Vec::operator = (const Vec &a)
{
	if (&a != this)
	{
		if(len != a.len)
		{
			len = a.len;
		}
		for (int i=0; i< len; i++)
			vec[i] = a.vec[i];
	}
	return *this;
}*/
// addition
Vec Vec::operator+ (const Vec &a)
{
	Vec c = Vec(a.len);
	if(a.len != this->len)
	{
		printf("vector length doesn't match!\n");
		return c;
	}

	for(int i=0; i<a.len; i++)
		c.vec[i] = a.vec[i] + this->vec[i];
	return c;
}
// subtraction
Vec Vec::operator- (const Vec a)
{
	Vec c = Vec(a.len);
	if(a.len != this->len)
	{
		printf("vector length doesn't match!\n");
		return c;
	}

	for(int i=0; i<a.len; i++)
		c.vec[i] = this->vec[i] - a.vec[i];
	return c;
}
// cross
Vec Vec::operator* (const Vec a)
{
	if(this->len != a.len)
	{
		printf("vector length doesn't match!\n");
		return *this;
	}

	Vec tmp = Vec(this->len);
	tmp[0] = this->vec[1]*a.vec[2] - a.vec[1]*this->vec[2];
	tmp[1] = this->vec[2]*a.vec[0] - a.vec[2]*this->vec[0];
	tmp[2] = this->vec[0]*a.vec[1] - a.vec[0]*this->vec[1];
	return tmp;
}
// multiply
Vec Vec::operator* (const double a)
{
	Vec tmp = Vec(this->len);
	for(int i=0; i<len; i++)
		tmp.vec[i] = this->vec[i] * a;
	return tmp;
}

Vec operator*(double a, Vec &b)
{
	return b * a;
}

Vec Vec::operator/ (double a)
{
	Vec tmp = Vec(this->len);
	for(int i=0; i<len; i++)
		tmp.vec[i] = this->vec[i] / a;
	return tmp;
}
