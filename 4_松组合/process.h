#ifndef _PREPROCESS_H
#define _PREPROCESS_H

#include "global_variables.h"
#include "ins.h"

typedef struct gnss_data_tag
{
	int nepoch;
	float *gpstime;
	double *lat; //deg
	double *lon; //deg
	double *h;   //height
	double *Nstd;
	double *Estd;
	double *Ustd;
	gnss_data_tag()
	{
		nepoch = 0;
		gpstime = new float[DATASUM];
		lat = new double[DATASUM];
		lon = new double[DATASUM];
		h   = new double[DATASUM];
		Nstd = new double[DATASUM];
		Estd = new double[DATASUM];
		Ustd = new double[DATASUM];
	}
	~gnss_data_tag()
	{
		delete []lat;
		delete []lon;
		delete []h;
		delete []Nstd;
		delete []Estd;
		delete []Ustd;
	}
}GnssData;

typedef struct parameters_tag
{
	Vec antenna_arm;
	double start_time;
	double BLH[3];
	double BLH_std[3];
	double Vn[3];
	double Vn_std[3];
	// Roll, pitch, heading
	double att[3]; 
	double att_std[3];
	double ARW;
	double VRW;
	double G_std;
	double A_std;
	double G_k;
	double A_k;
	double rel_t; // relative time
}Paras;

typedef struct tag_info
{
	double dt;
	Vec att; // roll, pitch, yaw
	Vec vel; // vn, ve, vd
	Vec pos; // N, E, D (m,m,m)
	Vec dgm;
	Vec dam;
	tag_info()
	{
		dt = 0.0;
		for(int i = 0;i<3;i++)
		{
			att[i] = 0.0;
			vel[i] = 0.0;
			pos[i] = 0.0;
			dgm[i] = 0.0;
			dam[i] = 0.0;
		}
	}
}Info;

typedef struct tag_corrimu
{
	Vec dgym;
	Vec dacc;
	Vec Kgym;
	Vec Kacc;
	tag_corrimu()
	{
		for(int i = 0;i<3; i++)
		{
			dgym[i] = 0.0;
			dacc[i] = 0.0;
			Kgym[i] = 0.0;
			Kacc[i] = 0.0;
		}
	}
}CorrImu;

// read ins file
int read_imufile(string filename, InsData *imu_data);
// read gnss navigation result
int read_gnssfile(string filename, GnssData *gnss_data);
// initial statement
void init_parameters(Paras *init_state);
// initial imu state
void init_ins_state(Paras *init_state, InitState *ins_state0);
// output imu result
void output_ins_result(InsState *imu_state);
// norm little att difference
double norm_arc(double deg1);

// initial matrix
void init_matrix(Paras *init_paras);
// update matrix
void update_matrix(Paras *init_paras);
// get IMU data
int get_imu_data(Info *pre_info, Info *cur_info, InsData *imu_data, InitState *state0, int iepoch, double dt);
// get GNSS data
int get_gps_data(Vec &gps_pos, GnssData *gnss_result, double t);
// correct attitude
void correct_att(Vec datt, Vec &att);
// correct arm level error
void correct_arm(Vec &gps_pos, Paras *init_paras, Info *pre_info);
// correct gym and acc
void correct_imu(Info * info, CorrImu * corr);
// update observation matrix
void update_HR_matrix(Vec att, Paras *init_para, GnssData *gnss_result, int iepoch);
// update parameter matrix
void update_transform_matrix(double rel_t, Info *pre_info);
// update Z matrix
void update_Z(Vec &ins_pos, Paras *init_paras, Info *cur_info);
// unscented Kalman filter
void kalman_filter(char mode);
// feedback
void feedback(CorrImu *corr, Info *cur_info);

// integrate GPS and IMU
int integration(GnssData *gnss_result, InsData *imue_data, Paras *init_paras, InsState *integrate_result);
// GPS/IMU process
int algorithm();
// compare result
void compare_result(InsState *real_state, InsState *com_result);


#endif