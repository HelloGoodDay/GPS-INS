#include "vec.h"

void multiply(double *vec1, double a, double *vsum)
{
	int i;
	if(vsum == NULL || vec1 == NULL)
	{
		printf("vector is empty, please check\n");
		return;
	}
	for(i = 0; i<3; i++)
	{
		*vsum = *vec1 * a;
		vsum++; vec1++;
	}
}

void cross(double *vec1, double *vec2, double *vsum)
{
	if(vsum == NULL || vec1 == NULL || vec2 == NULL)
	{
		printf("vector is empty, please check\n");
		return;
	}
	vsum[0] = vec1[1]*vec2[2] - vec2[1]*vec1[2];
	vsum[1] = vec1[2]*vec2[0] - vec2[2]*vec1[0];
	vsum[2] = vec1[0]*vec2[1] - vec2[0]*vec1[1];
}

void add(double *vec1, double *vec2, double *vsum)
{
	int i = 0;
	if(vsum == NULL || vec1 == NULL || vec2 == NULL)
	{
		printf("vector is empty, please check\n");
		return;
	}
	for (i = 0;i<3; i++)
		vsum[i] = vec1[i] + vec2[i];
}

void subtract(double *vec1, double *vec2, double *vsum)
{
	int i;
	if(vsum == NULL || vec1 == NULL || vec2 == NULL)
	{
		printf("vector is empty, please check\n");
		return;
	}
	for (i = 0;i<3; i++)
		vsum[i] = vec1[i] - vec2[i];
}

void division(double *vec1, double a, double *vsum)
{
	int i;
	if(vsum== NULL || vec1 == NULL)
	{
		printf("vector is empty, please check\n");
		return;
	}
	for(i = 0; i<3; i++)
	{
		*vsum = *vec1 / a;
		vsum++; vec1++;
	}
}