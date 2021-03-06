#ifndef INS_H
#define INS_H

/*
b frame: ǰ����
n frame: NEU
*/

#include "vec.h"
#include "global_variables.h"

typedef struct init_state_tag
{
	//time
	double t;
	//latitude, longitude(rad), height(m)
	double phi, lamda, h;
	//velocity in north, east and down
	double v_n, v_e, v_u;
	//posture, unit -- rad, ǰ����
	double roll, pitch, yaw;
	init_state_tag()
	{
		t = 91620.0;
		phi = degree2arc(23.1373950708);
		lamda = degree2arc(113.3713651222);
		h = 2.175;
		v_n = 0.0;
		v_e = 0.0;
		v_u = 0.0;
		roll = degree2arc(0.0107951084511778);
		pitch = degree2arc(-2.14251290749072);
		yaw = degree2arc(-75.7498049314083);
	}
}InitState;

typedef struct State_tag
{
	int nepoch;
	double *t;
	// roll, pitch, yaw
	double *att; // attitude, arc
	// NEU
	double *vel; // velocity, m/s
	// BLH
	double *pos; // position, arc, arc, m
	State_tag()
	{
		nepoch = 0;
		t = new double[DATASUM]; 
		att = new double[DATASUM * 3];
		vel = new double[DATASUM * 3];
		pos = new double[DATASUM * 3];
	}
	~State_tag()
	{
		delete []t;
		delete []att;
		delete []vel;
		delete []pos;
	}
}InsState;

typedef struct EarthPara_tag
{
	double RM;
	double RN;
	Vec wie;
	Vec win;
	Vec gcc;
	EarthPara_tag()
	{
		RM = RN = 0.0;
	}
}EarthPara;

typedef struct InsData_tag
{
	int ncount;
	double *t; //time
	double *dg; // angle increment
	double *da; // acceleration increment
	InsData_tag()
	{
		ncount = 0;
		t = new double[DATASUM];
		dg = new double[DATASUM * 3];
		da = new double[DATASUM * 3];
	}
	~InsData_tag()
	{
		delete []t;
		delete []dg;
		delete []da;
	}
}InsData; // read increment in each interval

class INS
{
private:

public:
	// get [I - (O.5phi)X]
	static void halfrv2m(Vec rv, double mat[3][3]);
	// [I - (rv)X]
	static void rv2m(Vec rv, double mat[3][3]);
	// transform q into attitude
	static void q2att(Vec q, Vec &att);
	// rotating vector to q
	static void rv2q(Vec rv, Vec &q);
	// transform attitude into q
	static void att2q(Vec att, Vec &q);
	// transform q into matrix
	static void q2mat(Vec q, double mat[3][3]);
	// transform attitude into DCM
	static void att2m(Vec att, double mat[3][3]);
	// transform DCB into attitude
	static void m2att(double mat[3][3], Vec &att);
	// anti-symmetric matrix
	static void asym(Vec rv, double a[3][3]);
	// multiply for q
	static Vec q_multiply(Vec p, Vec q);
	// normalize q
	static void q_norm(Vec &q);
	
	// get gravity
	static double get_gravity(double B, double h);
	// get C_bb
	static Vec get_Cbb(double mat[3][3], double Cbn[3][3], Vec V_bb);

	// get earth parameter
	static void earth_parameter(double phi, double h, Vec vel, EarthPara *eth);
	// update algorithm in one epoch
	static int ins_update(Vec dam, Vec dam_1, Vec dgm, Vec dgm_1, Vec att, Vec vel, Vec pos, Vec &att2, Vec &vel2, Vec &pos2, double interval);
};
#endif
