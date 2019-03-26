#include "ins.h"

void init_state(InitState * state0)
{	
	state0->t = 91620.0;
	state0->phi = degree2arc(23.1373950708);
	state0->lamda = degree2arc(113.3713651222);
	state0->h = 2.175;
	state0->v_n = 0.0;
	state0->v_e = 0.0;
	state0->v_u = 0.0;
	state0->roll = degree2arc(0.0107951084511778);
	state0->pitch = degree2arc(-2.14251290749072);
	state0->yaw = degree2arc(-75.7498049314083);
}

int read_file_2(char* filename, InsData *ins_data)
{
	FILE *fp;
	int ncount = 0;
	double buf[7];
	// 以"r"模式打开时，读到0x1A会意外终止
	if((fp = fopen(filename, "rb")) == NULL)
	{
		printf("can't open file %s\n", filename);
		return 0;
	}
	
	while(fread(&buf, sizeof(double), 7, fp))
	{
		(ins_data + ncount)->t = buf[0];
		(ins_data + ncount)->dg[0] = buf[1];
		(ins_data + ncount)->dg[1] = buf[2];
		(ins_data + ncount)->dg[2] = buf[3];
		(ins_data + ncount)->da[0] = buf[4];
		(ins_data + ncount)->da[1] = buf[5];
		(ins_data + ncount)->da[2] = buf[6];
		(ins_data + ncount)->flag = 1;
		ncount++;
	}
	fclose(fp);
	printf("read file successfully, the epoch number is %d\n", ncount);
	return ncount;
}

int read_file(char* filename, InsData *ins_data)
{
	FILE *fp;
	int ncount = 0;
	char line[256];
	double buf[7];
	// 以"r"模式打开时，读到0x1A会意外终止
	if((fp = fopen(filename, "r")) == NULL)
	{
		printf("can't open file %s\n", filename);
		return 0;
	}

	while(fgets(line, 256, fp))
	{
		sscanf(line, "%lf%lf%lf%lf%lf%lf%lf", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5], &buf[6]);
		(ins_data + ncount)->t = buf[0];
		(ins_data + ncount)->dg[0] = buf[1] * INTERVAL;
		(ins_data + ncount)->dg[1] = buf[2] * INTERVAL;
		(ins_data + ncount)->dg[2] = buf[3] * INTERVAL;
		(ins_data + ncount)->da[0] = buf[4] * INTERVAL;
		(ins_data + ncount)->da[1] = buf[5] * INTERVAL;
		(ins_data + ncount)->da[2] = buf[6] * INTERVAL;
		(ins_data + ncount)->flag = 1;
		ncount++;
	}
	fclose(fp);
	printf("read file successfully, the epoch number is %d\n", ncount);
	return ncount;
}

double get_gravity(double B, double h)
{
	double a = EARTH_a;
	double b = EARTH_b;
	double f = EARTH_f;
	double m = EARTH_V*EARTH_V*a*a*b/GM;
	double ga = 9.7803267715;
	double gb = 9.8321863685;
	double gl = (a*ga*cos(B)*cos(B) + b*gb*sin(B)*sin(B))/
		sqrt(a*a*cos(B)*cos(B) + b*b*sin(B)*sin(B));
	double g = gl*(1 - 2*h/a*(1+f+m-2*f*sin(B)*sin(B)) + 3*h*h/a/a);
	return g;
}

void rv2m(double rv[3], double mat[3][3])
{
	double cross_m[3][3];
	double rv2[3];
	division(rv, 2.0, rv2);
	cross_m[0][0] = 0;      cross_m[0][1] = -rv2[2]; cross_m[0][2] = rv2[1];
	cross_m[1][0] = rv2[2]; cross_m[1][1] = 0;       cross_m[1][2] = -rv2[0];
	cross_m[2][0] = -rv2[1];cross_m[2][1] = rv2[0];  cross_m[2][2] = 0;
	//get mat
	mat[0][0] = 1 - cross_m[0][0];
	mat[0][1] = 0 - cross_m[0][1]; 
	mat[0][2] = 0 - cross_m[0][2]; 
	mat[1][0] = 0 - cross_m[1][0]; 
	mat[1][1] = 1 - cross_m[1][1]; 
	mat[1][2] = 0 - cross_m[1][2]; 
	mat[2][0] = 0 - cross_m[2][0]; 
	mat[2][1] = 0 - cross_m[2][1]; 
	mat[2][2] = 1 - cross_m[2][2]; 
}

void get_Cbb(double mat[3][3], double Cbn[3][3],double V_bb[3], double Cbb[3])
{
	int i, j, k;
	// mat * Cbn
	double tmp[3][3];
	for(i = 0; i<3; i++)
	{
		for (j = 0;j<3; j++)
		{
			tmp[i][j] = 0.0;
			for (k = 0;k<3;k++)
				tmp[i][j] += (mat[i][k] * Cbn[k][j]);
		}
	}
	// tmp * V_bb
	Cbb[0] = tmp[0][0]*V_bb[0] + tmp[0][1]*V_bb[1] + tmp[0][2]*V_bb[2];
	Cbb[1] = tmp[1][0]*V_bb[0] + tmp[1][1]*V_bb[1] + tmp[1][2]*V_bb[2];
	Cbb[2] = tmp[2][0]*V_bb[0] + tmp[2][1]*V_bb[1] + tmp[2][2]*V_bb[2];
}

void att2m(double rv[3], double mat[3][3])
{
	double v1 = rv[0];
	double v2 = rv[1];
	double v3 = rv[2];
	mat[0][0] = cos(v2) * cos(v3);
	mat[0][1] = -cos(v1)*sin(v3) + sin(v1)*sin(v2)*cos(v3);
	mat[0][2] =  sin(v1)*sin(v3) + cos(v1)*sin(v2)*cos(v3);
	mat[1][0] = cos(v2)*sin(v3);
	mat[1][1] =  cos(v1)*cos(v3) + sin(v1)*sin(v2)*sin(v3);
	mat[1][2] = -sin(v1)*cos(v3) + cos(v1)*sin(v2)*sin(v3);
	mat[2][0] = -sin(v2);
	mat[2][1] = sin(v1)*cos(v2);
	mat[2][2] = cos(v1)*cos(v2);
}

void rv2q(double rv[3], double q[4])
{
	double rtm = sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2]);
	double q0, s;
	int i;
	if(rtm < 1e-8)  //mode is small, taylor expansion, for low velocity update
	{
		q0 = 1 - rtm*rtm*(1/8.0 - rtm*rtm/384.0);
		s = 1/2.0 - rtm*rtm/48;
	}
	else
	{
		q0 = cos(rtm / 2);
		s = sin(rtm/2) / rtm;
	}
	q[0] = q0; 
	for (i=1;i<4;i++)
		q[i] = rv[i - 1] * s;
}

void att2q(double att[3], double q[4])
{
	double ang1 = att[0];
	double ang2 = att[1];
	double ang3 = att[2];
	q[0] = cos(ang1/2)*cos(ang2/2)*cos(ang3/2) + sin(ang1/2)*sin(ang2/2)*sin(ang3/2);
	q[1] = sin(ang1/2)*cos(ang2/2)*cos(ang3/2) - cos(ang1/2)*sin(ang2/2)*sin(ang3/2);
	q[2] = cos(ang1/2)*sin(ang2/2)*cos(ang3/2) + sin(ang1/2)*cos(ang2/2)*sin(ang3/2);
	q[3] = cos(ang1/2)*cos(ang2/2)*sin(ang3/2) - sin(ang1/2)*sin(ang2/2)*cos(ang3/2);
}

void q2att(double q[3], double att[3])
{
	// q2DCM
	double DCM[3][3];
	DCM[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
	DCM[0][1] = 2*(q[1]*q[2] - q[0]*q[3]);
	DCM[0][2] = 2*(q[1]*q[3] + q[0]*q[2]);
	DCM[1][0] = 2*(q[1]*q[2] + q[0]*q[3]);
	DCM[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
	DCM[1][2] = 2*(q[2]*q[3] - q[0]*q[1]);
	DCM[2][0] = 2*(q[1]*q[3] - q[0]*q[2]);
	DCM[2][1] = 2*(q[2]*q[3] + q[0]*q[1]);
	DCM[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
	//DCM2att
	att[0] = atan2(DCM[2][1], DCM[2][2]);
	att[1] = atan(-DCM[2][0] / sqrt(DCM[2][1]*DCM[2][1] + DCM[2][2]*DCM[2][2]));
	att[2] = atan2(DCM[1][0], DCM[0][0]);
}

void q_multiply(double p[4], double q[4], double qq[4])
{
	qq[0] = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3];
	qq[1] = p[1]*q[0] + p[0]*q[1] - p[3]*q[2] + p[2]*q[3];
	qq[2] = p[2]*q[0] + p[3]*q[1] + p[0]*q[2] - p[1]*q[3];
	qq[3] = p[3]*q[0] - p[2]*q[1] + p[1]*q[2] + p[0]*q[3];
}

void q_norm(double q[4])
{
	int i;
	double qs = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
	for (i = 0;i<4;i ++)
		q[i] = q[i] / qs;
}

void earth_parameter(double phi, double h, double vel[3], EarthPara *eth)
{
	double e2 = EARTH_e * EARTH_e;
	double RN = EARTH_a / sqrt(1 - e2 * sin(phi) * sin(phi));
	double RM = RN * (1-e2) / (1 - e2 * sin(phi) * sin(phi)); 
	double wie[3], wen[3], win[3];
	double gl_v[3], gcc[3], v_tmp[3];
	int i;

	wie[0] = EARTH_V*cos(phi); wie[1] = 0; wie[2] = -EARTH_V*sin(phi);
	wen[0] = vel[1]/(RN+h); wen[1] = -vel[0]/(RM+h); wen[2] = -vel[1]/(RN+h)*tan(phi);
	add(wie, wen, win);
	
	gl_v[0] = 0; gl_v[1] = 0;gl_v[2] = get_gravity(phi, h);
	add(win, wie, gcc);
	cross(gcc, vel, v_tmp);
	subtract(gl_v, v_tmp, gcc);

	eth->RM = RM;
	eth->RN = RN;
	for(i = 0; i<3; i++)
	{
		eth->win[i] = win[i];
		eth->gcc[i] = gcc[i];
	}
}

int ins_algorithm(InitState state0, InsData *ins_data, InsState *ins_state)
{
	int i, iepoch, nepoch;
	int first_flag = 1;
	InsData ins_data0, ins_data1;
	InsState ins_state0, ins_state1;
	double check_dt, check_dg;
	
	// check
	check_dt = (ins_data+1)->t - ins_data->t;
	check_dg = ins_data->da[2] / INTERVAL;
	if ( fabs(check_dt - INTERVAL) > 1e-4)
	{
		printf("\ninterval error!\n\n");
		return 0;
	}
	if ( fabs(check_dg + 9.8) > 0.5)
	{
		printf("\ndata error, please check file read function\n\n");
		return 0;
	}
	// initial state
	ins_state->att[0] = state0.roll;
	ins_state->att[1] = state0.pitch;
	ins_state->att[2] = state0.yaw;
	ins_state->vel[0] = state0.v_n;
	ins_state->vel[1] = state0.v_e;
	ins_state->vel[2] = state0.v_u;
	ins_state->pos[0] = state0.phi;
	ins_state->pos[1] = state0.lamda;
	ins_state->pos[2] = state0.h;
	ins_state->t = state0.t;

	iepoch = nepoch = 1;
	memset(&ins_state1, 0, sizeof(InsState));
	
	//--------------------------- precess --------------------------------------
	while((ins_data+nepoch)->flag == 1)
	{
		if ((ins_data+nepoch)->t <= state0.t)
		{
			nepoch++;
			continue;
		}
		
		if(first_flag)
		{
			for (i = 0; i<3; i++)
			{
				ins_data0.da[i] = 0;
				ins_data0.dg[i] = 0;
				ins_data1.da[i] = (ins_data + nepoch)->da[i];
				ins_data1.dg[i] = (ins_data + nepoch)->dg[i];
				ins_state0.att[i] = ins_state->att[i];
				ins_state0.vel[i] = ins_state->vel[i];
				ins_state0.pos[i] = ins_state->pos[i];
			}
			first_flag = 0;
		}
		else
		{
			for (i = 0; i<3; i++)
			{
				ins_data0.da[i] = ins_data1.da[i];
				ins_data0.dg[i] = ins_data1.dg[i];
				ins_data1.da[i] = (ins_data + nepoch)->da[i];
				ins_data1.dg[i] = (ins_data + nepoch)->dg[i];
				ins_state0.att[i] = ins_state1.att[i];
				ins_state0.vel[i] = ins_state1.vel[i];
				ins_state0.pos[i] = ins_state1.pos[i];
			}
		}
		// ins algorithm
		ins_update(&ins_data0, &ins_data1, &ins_state0, &ins_state1);
		// update
		for(i=0; i<3;i++)
		{
			(ins_state+iepoch)->att[i] = ins_state1.att[i];
			(ins_state+iepoch)->vel[i] = ins_state1.vel[i];
			(ins_state+iepoch)->pos[i] = ins_state1.pos[i];
		}
		(ins_state+iepoch)->t = (ins_data+nepoch)->t; 
		(ins_state+iepoch)->nepoch = iepoch;
		nepoch++;
		iepoch++;

	}
	return 1;
}

int update_vel(InsData *ins0, InsData *ins1, InsState *stt0, EarthPara *eth, InsState *stt1)
{
	double v_gcor[3], kesy_v[3], v_bb[3], v_f[3];
	double mat[3][3], C_bn[3][3];
	double temp_v[3];
	double vel[3];
	int i;
	double interval = INTERVAL;
	// gravity compensation
	// -> v_gcor = eth.gcc * interval;
	multiply(eth->gcc, interval, v_gcor);
	// cone and sculling error compensation
	multiply(eth->win, interval, kesy_v);
	rv2m(kesy_v, mat);
	att2m(stt0->att, C_bn);
	// get v_bb
	// -> v_bb = dam + dgm * dam / 2.0 + (dgm_1*dam + dam_1*dgm) / 12.0;
	cross(ins0->dg, ins1->da, temp_v);
	cross(ins0->da, ins1->dg, v_bb);
	add(temp_v, v_bb, v_bb);
	division(v_bb, 12.0, v_bb);
	cross(ins1->dg, ins1->da, temp_v);
	division(temp_v, 2.0, temp_v);
	add(ins1->da, temp_v, temp_v);
	add(temp_v, v_bb, v_bb);
	// get Cbb
	// -> v_f = get_Cbb(mat, C_bn, v_bb);
	get_Cbb(mat, C_bn, v_bb, v_f);
	// update 
	// ->vel2 = vel + v_f + v_gcor;
	add(stt0->vel, v_f, vel);
	add(vel, v_gcor, vel);
	for(i = 0;i<3;i++)
		stt1->vel[i] = vel[i];
	return 1;
}

int update_pos(InsData *ins0, InsData *ins1, InsState *stt0, EarthPara *eth, InsState *stt1)
{
	double vel_ave[3];
	double h_ave, phi_ave, dh[3];
	// height
	add(stt0->vel, stt1->vel, vel_ave);
	division(vel_ave, 2, vel_ave);
	stt1->pos[2] = stt0->pos[2] - vel_ave[2] * INTERVAL;
	// B
	h_ave = (stt0->pos[2] + stt1->pos[2])/ 2.0;
	stt1->pos[0] = stt0->pos[0] + vel_ave[0] * INTERVAL / (eth->RM + h_ave);
	// L
	phi_ave = (stt0->pos[0] + stt1->pos[0]) / 2.0;
	earth_parameter(phi_ave, h_ave, vel_ave, eth);
	stt1->pos[1] = stt0->pos[1] + vel_ave[1] * INTERVAL/ ((eth->RN + h_ave)*cos(phi_ave));
	return 1;
}

int update_att(InsData *ins0, InsData *ins1, InsState *stt0, EarthPara *eth, InsState *stt1)
{
	double q_bn[4], q_bb[4], q_nn[4], q_temp[4];
	double phi_bb[3];
	double kesy[3];
	// q in the first epoch
	att2q(stt0->att, q_bn);
	//get qbb
	// -> phi_bb = dgm + dgm_1 * dgm / 12.0;
	cross(ins0->dg, ins1->dg, phi_bb);
	division(phi_bb, 12.0, phi_bb);
	add(ins1->dg, phi_bb, phi_bb);
	rv2q(phi_bb, q_bb);
	//get qnn
	multiply(eth->win, -1 * INTERVAL, kesy);
	rv2q(kesy, q_nn);
	//update
	q_multiply(q_nn, q_bn, q_temp);
	q_multiply(q_temp, q_bb, q_bn);
	q_norm(q_bn);
	q2att(q_bn, stt1->att);
	return 1;
}

int ins_update(InsData *ins0, InsData *ins1, InsState *stt0, InsState *stt1)
{
	EarthPara eth;
	double interval = INTERVAL;

	earth_parameter(stt0->pos[0], stt0->pos[2], stt0->vel, &eth);
	// velocity update---------------------------------------------
	update_vel(ins0, ins1, stt0, &eth, stt1);
	// position update---------------------------------------------
	update_pos(ins0, ins1, stt0, &eth, stt1);
	// posture update---------------------------------------------
	update_att(ins0, ins1, stt0, &eth, stt1);
	return 1;
}

