#include "ins.h"

int INS::read_file(string filename, InsData *ins_data)
{
	FILE *fp;
	// 以"r"模式打开时，读到0x1A会意外终止
	if((fp = fopen(filename.c_str(), "rb")) == NULL)
	{
		cout<<"can't open file "<<filename.c_str()<<endl;
		return 0;
	}

	int ncount = 0;
	double buf[7];
	while(!feof (fp))
	{
		try
		{
			fread(&buf, sizeof(double), 7, fp);
		}
		catch(...)
		{
			cout<<"file format error!"<<endl;
			return 0;
		}
		ins_data->t[ncount] = buf[0];
		ins_data->dg[ncount*3 + 0] = buf[1];
		ins_data->dg[ncount*3 + 1] = buf[2];
		ins_data->dg[ncount*3 + 2] = buf[3];
		ins_data->da[ncount*3 + 0] = buf[4];
		ins_data->da[ncount*3 + 1] = buf[5];
		ins_data->da[ncount*3 + 2] = buf[6];
		ncount++;
	}
	ins_data->ncount = ncount;
	fclose(fp);
	cout<<"read file successfully, the epoch number is "<<ncount<<endl;
	return 1;
}

double INS::get_gravity(double B, double h)
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
	
	/*double a0 = 9.780325;
	double a1 = 0.00530240;
	double a2 = 0.00000582;
	double gl = a0 * (1+ a1*sin(B)*sin(B) - a2*sin(2*B)*sin(2*B));
	double gh = gl - 3.08e-6 * h;
	return gh;*/
}

void INS::rv2m(Vec rv, double mat[3][3])
{
	double cross_m[3][3];
	rv = rv/2;
	cross_m[0][0] = 0;     cross_m[0][1] = -rv[2]; cross_m[0][2] = rv[1];
	cross_m[1][0] = rv[2]; cross_m[1][1] = 0;      cross_m[1][2] = -rv[0];
	cross_m[2][0] = -rv[1];cross_m[2][1] = rv[0];  cross_m[2][2] = 0;
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

Vec INS::get_Cbb(double mat[3][3], double Cbn[3][3], Vec V_bb)
{
	// mat * Cbn
	double tmp[3][3];
	for(int i = 0; i<3; i++)
	{
		for (int j = 0;j<3; j++)
		{
			tmp[i][j] = 0.0;
			for (int k = 0;k<3;k++)
				tmp[i][j] += (mat[i][k] * Cbn[k][j]);
		}
	}
	// tmp * V_bb
	Vec Cbb(3);
	Cbb[0] = tmp[0][0]*V_bb[0] + tmp[0][1]*V_bb[1] + tmp[0][2]*V_bb[2];
	Cbb[1] = tmp[1][0]*V_bb[0] + tmp[1][1]*V_bb[1] + tmp[1][2]*V_bb[2];
	Cbb[2] = tmp[2][0]*V_bb[0] + tmp[2][1]*V_bb[1] + tmp[2][2]*V_bb[2];
	return Cbb;
}

void INS::att2m(Vec rv, double mat[3][3])
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

void INS::rv2q(Vec rv, Vec &q)
{
	double rtm = sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2]);
	double q0, s;
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
	for (int i=1;i<4;i++)
		q[i] = rv[i - 1] * s;
}

void INS::att2q(Vec att, Vec &q)
{
	double ang1 = att[0];
	double ang2 = att[1];
	double ang3 = att[2];
	q[0] = cos(ang1/2)*cos(ang2/2)*cos(ang3/2) + sin(ang1/2)*sin(ang2/2)*sin(ang3/2);
	q[1] = sin(ang1/2)*cos(ang2/2)*cos(ang3/2) - cos(ang1/2)*sin(ang2/2)*sin(ang3/2);
	q[2] = cos(ang1/2)*sin(ang2/2)*cos(ang3/2) + sin(ang1/2)*cos(ang2/2)*sin(ang3/2);
	q[3] = cos(ang1/2)*cos(ang2/2)*sin(ang3/2) - sin(ang1/2)*sin(ang2/2)*cos(ang3/2);
}

void INS::q2att(Vec q, Vec &att)
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

Vec INS::q_multiply(Vec p, Vec q)
{
	Vec qq(4);
	qq[0] = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3];
	qq[1] = p[1]*q[0] + p[0]*q[1] - p[3]*q[2] + p[2]*q[3];
	qq[2] = p[2]*q[0] + p[3]*q[1] + p[0]*q[2] - p[1]*q[3];
	qq[3] = p[3]*q[0] - p[2]*q[1] + p[1]*q[2] + p[0]*q[3];
	return qq;
}

void INS::q_norm(Vec &q)
{
	double qs = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
	for (int i = 0;i<4;i ++)
		q[i] = q[i] / qs;
}

void INS::earth_parameter(double phi, double h, Vec vel, EarthPara *eth)
{
	double e2 = EARTH_e * EARTH_e;
	double RN = EARTH_a / sqrt(1 - e2 * sin(phi) * sin(phi));
	double RM = RN * (1-e2) / (1 - e2 * sin(phi) * sin(phi)); 
	Vec wie, wen, win;
	wie[0] = EARTH_V*cos(phi); wie[1] = 0; wie[2] = -EARTH_V*sin(phi);
	wen[0] = vel[1]/(RN+h); wen[1] = -vel[0]/(RM+h); wen[2] = -vel[1]/(RN+h)*tan(phi);
	win = wie + wen;
	
	Vec gl_v;
	gl_v[0] = 0; gl_v[1] = 0;gl_v[2] = get_gravity(phi, h);
	Vec gcc = gl_v - (wie+ win) * vel;

	eth->RM = RM;
	eth->RN = RN;
	eth->win = win;
	eth->gcc = gcc;
}

int INS::ins_algorithm(InitState state0, InsData *ins_data, InsState *ins_state)
{
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
	ins_state->t[0] = state0.t;
	
	//--------------------------- precess --------------------------------------
	bool first_flag = true;
	Vec dam, dam_1, dgm, dgm_1;
	Vec att, vel, pos;
	Vec att2, vel2, pos2;
	int nepoch = 1;
	for (int iepoch = 0;iepoch < ins_data->ncount;iepoch++)
	{
		if (ins_data->t[iepoch] <= state0.t)
			continue;
		
		if(first_flag)
		{
			dam_1 = dam;
			dgm_1 = dgm;
			for (int i = 0; i<3; i++)
			{
				dgm[i] = ins_data->dg[iepoch*3 + i];
				dam[i] = ins_data->da[iepoch*3 + i];
				att[i] = ins_state->att[i];
				vel[i] = ins_state->vel[i];
				pos[i] = ins_state->pos[i];
			}
			first_flag = false;
		}
		else
		{
			dam_1 = dam;
			dgm_1 = dgm;
			for (int i = 0; i<3; i++)
			{
				dgm[i] = ins_data->dg[iepoch*3 + i];
				dam[i] = ins_data->da[iepoch*3 + i];
				att[i] = att2[i];
				vel[i] = vel2[i];
				pos[i] = pos2[i];
			}
		}
		// ins algorithm
		ins_update(dam, dam_1, dgm, dgm_1, att, vel, pos, att2, vel2, pos2);
		// update
		for(int i=0; i<3;i++)
		{
			ins_state->att[nepoch*3 + i] = att2[i];
			ins_state->vel[nepoch*3 + i] = vel2[i];
			ins_state->pos[nepoch*3 + i] = pos2[i];
		}
		ins_state->t[nepoch] = ins_data->t[iepoch]; 
		nepoch++;
	}
	ins_state->nepoch = nepoch;
	return 1;
}

int INS::ins_update(Vec dam, Vec dam_1, Vec dgm, Vec dgm_1, Vec att, Vec vel, Vec pos, Vec &att2, Vec &vel2, Vec &pos2)
{
	EarthPara eth;
	double interval = INTERVAL;
	// velocity update---------------------------------------------
	// gravity compensation
	earth_parameter(pos[0], pos[2], vel, &eth);
	Vec v_gcor = eth.gcc * interval;
	// cone and sculling error compensation
	Vec kesy_v = eth.win * interval;
	double mat[3][3];
	rv2m(kesy_v, mat);
	double C_bn[3][3];
	att2m(att, C_bn);
	Vec v_bb = dam + dgm * dam / 2.0 + (dgm_1*dam + dam_1*dgm) / 12.0;
	Vec v_f = get_Cbb(mat, C_bn, v_bb);
	// update
	vel2 = vel + v_f + v_gcor;

	// position update---------------------------------------------
	Vec vel_ave = (vel + vel2) / 2.0;
	pos2[2] = pos[2] - vel_ave[2] * interval;
	double h_ave = (pos[2] + pos2[2])/ 2.0;
	pos2[0] = pos[0] + vel_ave[0] * interval / (eth.RM + h_ave);
	double phi_ave = (pos[0] + pos2[0]) / 2.0;
	earth_parameter(phi_ave, h_ave, vel_ave, &eth);
	pos2[1] = pos[1] + vel_ave[1] * interval / ((eth.RN + h_ave)*cos(phi_ave));

	// posture update---------------------------------------------
	// q in the first epoch
	Vec q_bn, q_bb(4), q_nn(4);
	att2q(att, q_bn);
	//get qbb
	Vec phi_bb = dgm + dgm_1 * dgm / 12.0;
	rv2q(phi_bb, q_bb);
	//get qnn
	Vec kesy = -1 * eth.win * interval;
	rv2q(kesy, q_nn);
	//update
	q_bn = q_multiply((q_multiply(q_nn, q_bn)), q_bb);
	q_norm(q_bn);
	Vec tatt2(4);
	q2att(q_bn, tatt2);
	att2 = tatt2;

	return 1;
}

