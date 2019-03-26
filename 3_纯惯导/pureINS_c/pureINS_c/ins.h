#ifndef INS_H
#define INS_H

#define DATASUM 1000000
#define INTERVAL 0.01
#include "vec.h"
#include "global_variables.h"

typedef struct init_state_tag
{
	//time
	double t;
	//latitude, longitude(rad), height(m)
	double phi, lamda, h;
	//velocity in north, east and up
	double v_n, v_e, v_u;
	//posture, unit -- rad
	double roll, pitch, yaw;
}InitState;

typedef struct State_tag
{
	int nepoch;
	double t;
	double att[3]; // attitude, arc
	// roll, pitch, yaw
	double vel[3]; // velocity, m/s
	double pos[3]; // position, arc, arc, m
}InsState;

typedef struct EarthPara_tag
{
	double RM;
	double RN;
	double win[3];
	double gcc[3];
}EarthPara;

typedef struct InsData_tag
{
	int flag;
	double t; //time
	double dg[3]; // angle increment
	// 陀螺仪,rad/s，注意是增量
	double da[3]; // acceleration increment
	// 加速度计,m/s/s，注意是增量
}InsData;


// init state
void init_state(InitState * state0);
// get gravity
double get_gravity(double B, double h);
// get [I - (O.5phi)X]
void rv2m(double rv[3], double mat[3][3]);
// get C_bb
void get_Cbb(double mat[3][3], double Cbn[3][3],double V_bb[3], double Cbb[3]);
// transform attitude into DCM
void att2m(double rv[3], double mat[3][3]);
// rotating vector to q
void rv2q(double rv[3], double q[4]);
// transform attitude into q
void att2q(double att[3], double q[4]);
// transform q into attitude
void q2att(double q[3], double att[3]);
// multiply for q
void q_multiply(double p[4], double q[4], double qq[4]);
// normalize q
void q_norm(double q[4]);
// get earth parameter
void earth_parameter(double phi, double h, double vel[3], EarthPara *eth);

// binary -> read file, success - 1; fail - 0
int read_file_2(char* filename, InsData *ins_data);
// txt -> read file, success - 1; fail - 0
int read_file(char* filename, InsData *ins_data);
// update velocity
int update_vel(InsData *ins0, InsData *ins1, InsState *stt0, EarthPara *eth, InsState *stt1);
// update position
int update_pos(InsData *ins0, InsData *ins1, InsState *stt0, EarthPara *eth, InsState *stt1);
// update attitude
int update_att(InsData *ins0, InsData *ins1, InsState *stt0, EarthPara *eth, InsState *stt1);
// update algorithm in one epoch
int ins_update(InsData *ins0, InsData *ins1, InsState *stt0, InsState *stt1);
// pure ins algorithm
int ins_algorithm(InitState state0, InsData *ins_data, InsState *ins_state);

#endif
