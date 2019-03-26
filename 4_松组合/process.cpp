#include "process.h"
#include "ins.h"

int read_imufile(string filename, InsData *imu_data)
{
	FILE *fp;
	if((fp = fopen(filename.c_str(), "rb")) == NULL)
	{
		cout<<"can't open imu file"<<endl;
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
		imu_data->t[ncount] = buf[0];
		for (int i = 0;i<3;i++)
		{
			imu_data->dg[ncount*3 + i] = buf[i + 1]; // angle increment
			imu_data->da[ncount*3 + i] = buf[i + 4]; // acceleration increment
		}
		ncount++;
	}
	imu_data->ncount = ncount;
	cout<<"read imu file successfully, the epoch number is "<<ncount<<endl;
	fclose(fp);
	return 1;
}

int read_gnssfile(string filename, GnssData *gnss_data)
{
	INTERVALS = 1.0;
	FILE *fp;
	if((fp = fopen(filename.c_str(), "r")) == NULL)
	{
		cout<<"can't open gnss result file"<<endl;
		return 0;
	}

	char line[256];
	double tmp[7];
	int nepoch = 0;
	while(fgets(line, 255, fp))
	{
		if(strstr(line, "sec"))
			break;
	}
	while(fgets(line, 255, fp))
	{
		sscanf(line, "%lf%lf%lf%lf%lf%lf%lf", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6]);
		gnss_data->gpstime[nepoch] = tmp[0];
		gnss_data->lat[nepoch] = degree2arc(tmp[1]);
		gnss_data->lon[nepoch] = degree2arc(tmp[2]);
		gnss_data->Nstd[nepoch] = tmp[4];
		gnss_data->Estd[nepoch] = tmp[5];
		gnss_data->Ustd[nepoch] = tmp[6];
		gnss_data->h[nepoch] = tmp[3];
		nepoch++;
	}
	gnss_data->nepoch = nepoch;
	cout<<"read gnss result successfullly, the number of epoch is "<<nepoch<<endl;
	fclose(fp);
	return 1;
}

int read_result(string filename, InsState *real_state)
{
	FILE *fp;
	if ((fp = fopen(filename.c_str(), "r")) == NULL)
	{
		cout<<"can't read result file"<<endl;
		return 0;
	}

	char line[255+1];
	while(fgets(line, 255, fp))
	{
		if(strstr(line, " (sec) "))
			break;
	}
	int iepoch = 0; double last_gpstime = 0.0;
	while(fgets(line, 255, fp))
	{
		double gpstime, att[3], vel[3], pos[3];
		sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n", &gpstime, &pos[0], &pos[1], &pos[2], 
			&vel[0], &vel[1], &vel[2], &att[0], &att[1], &att[2]);
		//if(fabs(gpstime - last_gpstime) < INS_INTERVAL*0.8)
		//	continue;
		real_state->t[iepoch] = gpstime;
		real_state->att[iepoch*3 + 0] = degree2arc(att[0]);
		real_state->att[iepoch*3 + 1] = degree2arc(att[1]);
		real_state->att[iepoch*3 + 2] = degree2arc(att[2]);
		real_state->vel[iepoch*3 + 0] = vel[0];
		real_state->vel[iepoch*3 + 1] = vel[1];
		real_state->vel[iepoch*3 + 2] = vel[2];
		real_state->pos[iepoch*3 + 0] = degree2arc(pos[0]);
		real_state->pos[iepoch*3 + 1] = degree2arc(pos[1]);
		real_state->pos[iepoch*3 + 2] = pos[2];
		iepoch++;
		last_gpstime = gpstime;
	}
	real_state->nepoch = iepoch;
	cout<<"read result completely"<<endl;
	fclose(fp);
	return 1;
}

void init_parameters(Paras *init_state)
{
	init_state->antenna_arm[0] = 0.194;
	init_state->antenna_arm[1] = 0.278;
	init_state->antenna_arm[2] = -0.9;
	init_state->start_time = 442336;
	init_state->BLH[0] = degree2arc(30.5631420296);
	init_state->BLH[1] = degree2arc(114.4697206604);
	init_state->BLH[2] = 14.238;
	init_state->BLH_std[0] = 0.003;
	init_state->BLH_std[1] = 0.004;
	init_state->BLH_std[2] = 0.007;
	init_state->Vn[0] = 0.001;
	init_state->Vn[1] = -0.000;
	init_state->Vn[2] = -0.001;
	init_state->Vn_std[0] = 0.001;
	init_state->Vn_std[1] = 0.001;
	init_state->Vn_std[2] = 0.001;
	init_state->att[0] = degree2arc(0.95529927);
	init_state->att[1] = degree2arc(-0.44491903);
	init_state->att[2] = degree2arc(102.81089380);
	init_state->att_std[0] = degree2arc(0.001);
	init_state->att_std[1] = degree2arc(0.001);
	init_state->att_std[2] = degree2arc(0.038);
	init_state->ARW = degree2arc(0.0022/60);
	init_state->VRW = 0.00075/60;
	init_state->G_std = degree2arc(0.005/3600);
	init_state->A_std = 25E-5;
	init_state->G_k = 10E-6;
	init_state->A_k = 10E-6;
	init_state->rel_t = 1000;
}

void init_ins_state(Paras *init_state, InitState *ins_state0)
{
	ins_state0->phi = init_state->BLH[0];
	ins_state0->lamda = init_state->BLH[1];
	ins_state0->h = init_state->BLH[2];
	ins_state0->t = init_state->start_time;
	ins_state0->v_n = init_state->Vn[0];
	ins_state0->v_e = init_state->Vn[1];
	ins_state0->v_u = init_state->Vn[2];
	ins_state0->roll  = init_state->att[0];
	ins_state0->pitch = init_state->att[1];
	ins_state0->yaw   = init_state->att[2];
}

void output_ins_result(InsState *imu_state)
{
	FILE *fp;
	string filename = "../../data/ins_result.txt";
	if((fp = fopen(filename.c_str(), "w")) == NULL)
	{
		cout<<"can't open gnss result file"<<endl;
		return;
	}

	int nepoch = 0;
	nepoch = 1000;
	for (int iepoch = 0;iepoch < nepoch; iepoch++)
	{
		fprintf(fp, "%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf\n",
			imu_state->att[iepoch*3], imu_state->att[iepoch*3+1], imu_state->att[iepoch*3+2],
			imu_state->vel[iepoch*3], imu_state->vel[iepoch*3+1], imu_state->vel[iepoch*3+2],
			imu_state->pos[iepoch*3], imu_state->pos[iepoch*3+1], imu_state->pos[iepoch*3+2]);
	}
	fclose(fp);
}

void init_matrix(Paras *init_paras)
{
	F = dMatrix(para_sum, para_sum);
	H = dMatrix(3, para_sum);
	P0 = dMatrix(para_sum, para_sum);
	DX0 = dMatrix(para_sum, para_sum);
	G = dMatrix(para_sum, para_sum);
	Q = dMatrix(para_sum, para_sum);
	R = dMatrix(3, 3);
	X = dMatrix(para_sum, 1);
	Z = dMatrix(3, 1);
	// P
	DX0(0,0) = init_paras->att_std[0] * init_paras->att_std[0]; 
	DX0(1,1) = init_paras->att_std[1] * init_paras->att_std[1]; 
	DX0(2,2) = init_paras->att_std[2] * init_paras->att_std[2];
	DX0(3,3) = init_paras->Vn_std[0] * init_paras->Vn_std[0];
	DX0(4,4) = init_paras->Vn_std[1] * init_paras->Vn_std[1];
	DX0(5,5) = init_paras->Vn_std[2] * init_paras->Vn_std[2];
	DX0(6,6) = init_paras->BLH_std[0] * init_paras->BLH_std[0];
	DX0(7,7) = init_paras->BLH_std[1] * init_paras->BLH_std[1];
	DX0(8,8) = init_paras->BLH_std[2] * init_paras->BLH_std[2];
	DX0( 9, 9) = init_paras->G_std * init_paras->G_std;
	DX0(10,10) = init_paras->G_std * init_paras->G_std;
	DX0(11,11) = init_paras->G_std * init_paras->G_std;
	DX0(12,12) = init_paras->A_std * init_paras->A_std;
	DX0(13,13) = init_paras->A_std * init_paras->A_std;
	DX0(14,14) = init_paras->A_std * init_paras->A_std;
	DX0(15,15) = init_paras->G_k * init_paras->G_k;
	DX0(16,16) = init_paras->G_k * init_paras->G_k;
	DX0(17,17) = init_paras->G_k * init_paras->G_k;
	DX0(18,18) = init_paras->A_k * init_paras->A_k;
	DX0(19,19) = init_paras->A_k * init_paras->A_k;
	DX0(20,20) = init_paras->A_k * init_paras->A_k;
	for(int i = 0;i<para_sum;i++)
		for(int j = 0;j<para_sum;j++)
			P0(i,j) = MPk * DX0(i, j);
	// Q
	Q(0,0) = init_paras->ARW * init_paras->ARW;
	Q(1,1) = init_paras->ARW * init_paras->ARW;
	Q(2,2) = init_paras->ARW * init_paras->ARW;
	Q(3,3) = init_paras->VRW * init_paras->VRW;
	Q(4,4) = init_paras->VRW * init_paras->VRW;
	Q(5,5) = init_paras->VRW * init_paras->VRW;
	Q(9, 9 ) = 2 * init_paras->G_std * init_paras->G_std / init_paras->rel_t;
	Q(10,10) = 2 * init_paras->G_std * init_paras->G_std / init_paras->rel_t;
	Q(11,11) = 2 * init_paras->G_std * init_paras->G_std / init_paras->rel_t;
	Q(12,12) = 2 * init_paras->A_std * init_paras->A_std / init_paras->rel_t;
	Q(13,13) = 2 * init_paras->A_std * init_paras->A_std / init_paras->rel_t;
	Q(14,14) = 2 * init_paras->A_std * init_paras->A_std / init_paras->rel_t;
	Q(15,15) = 2 * init_paras->G_k * init_paras->G_k / init_paras->rel_t;
	Q(16,16) = 2 * init_paras->G_k * init_paras->G_k / init_paras->rel_t;
	Q(17,17) = 2 * init_paras->G_k * init_paras->G_k / init_paras->rel_t;
	Q(18,18) = 2 * init_paras->A_k * init_paras->A_k / init_paras->rel_t;
	Q(19,19) = 2 * init_paras->A_k * init_paras->A_k / init_paras->rel_t;
	Q(20,20) = 2 * init_paras->A_k * init_paras->A_k / init_paras->rel_t;
	for(int i = 0;i<para_sum;i++)
		for(int j = 0;j<para_sum;j++)
			Q(i,j) = MQk * Q(i, j);
	// R
	R(0,0) = 0.0085 * 0.0085;
	R(1,1) = 0.01 * 0.01;
	R(2,2) = 0.016 * 0.016;
	// H
	H(0,6) = 1;
	H(1,7) = 1;
	H(2,8) = 1;
}

void update_matrix(Paras *init_paras)
{
	for(int i = 0;i<para_sum;i++)
		X(i, 0) = 0.0;
}

double norm_arc(double deg1)
{
	if(deg1 > 180)
		deg1 -= 360;
	if(deg1 < -180)
		deg1 += 360;
	return deg1;
}

void update_transform_matrix(double rel_t, Info *pre_info)
{
	// parameters
	double Cbn[3][3];
	INS::att2m(pre_info->att, Cbn);
	EarthPara eth;
	double phi = pre_info->pos[0];
	double h = pre_info->pos[2];
	INS::earth_parameter(phi, h, pre_info->vel, &eth);
	double RM = eth.RM + h;
	double RN = eth.RN + h;
	double tl = tan(phi);
	double secl = 1.0 / cos(phi);
	double ve_rn = pre_info->vel[1]/RN;
	double ve_rn2 = pre_info->vel[1]/RN/RN;
	double vn_rm = pre_info->vel[0]/RM;
	double vn_rm2 = pre_info->vel[0]/RM/RM;
	double omega_sl = EARTH_V * sin(phi);
	double omega_cl = EARTH_V * cos(phi);
	double vn = pre_info->vel[0], ve = pre_info->vel[1], vd = pre_info->vel[2];
	
	// calculate matrix
	double Maa[3][3], Mav[3][3], Map[3][3];
	double Mva[3][3], Mvv[3][3], Mvp[3][3];
	double Mpa[3][3], Mpv[3][3], Mpp[3][3];
	// Maa
	INS::asym(eth.win, Maa);
	for(int i = 0;i<3;i++)
		for(int j = 0;j<3;j++)
			Maa[i][j] = -Maa[i][j];
	// Mav
	Mav[0][0] = 0.0;
	Mav[0][1] = 1.0/RN;
	Mav[0][2] = 0.0;
	Mav[1][0] = -1.0/RM;
	Mav[1][1] = 0.0;
	Mav[1][2] = 0.0;
	Mav[2][0] = 0.0;
	Mav[2][1] = -tl/RN;
	Mav[2][2] = 0.0;
	// Map
	Map[0][0] = -omega_sl / RM;
	Map[0][1] = 0.0;
	Map[0][2] = ve_rn2;
	Map[1][0] = 0.0;
	Map[1][1] = 0.0;
	Map[1][2] = -vn_rm2;
	Map[2][0] = -omega_cl/RM - ve_rn*secl*secl/RM;
	Map[2][1] = 0.0;
	Map[2][2] = -ve_rn2 * tl;
	// Mva
	Vec fb = pre_info->dam / pre_info->dt;
	Vec fn;
	fn[0] = Cbn[0][0]*fb[0] + Cbn[0][1]*fb[1] + Cbn[0][2]*fb[2];
	fn[1] = Cbn[1][0]*fb[0] + Cbn[1][1]*fb[1] + Cbn[1][2]*fb[2];
	fn[2] = Cbn[2][0]*fb[0] + Cbn[2][1]*fb[1] + Cbn[2][2]*fb[2];
	INS::asym(fn, Mva);

	// Mvv
	Mvv[0][0] = vd/RM;
	Mvv[0][1] = -2* (omega_sl + ve*tl/RN);
	Mvv[0][2] = vn/RM;
	Mvv[1][0] = 2*omega_sl + ve*tl/RN;
	Mvv[1][1] = (vn*tl + vd) / RN;
	Mvv[1][2] = 2*omega_cl + ve/RN;
	Mvv[2][0] = -2*vn_rm;
	Mvv[2][1] = -2* (omega_cl + ve_rn);
	Mvv[2][2] = 0.0;
	// Mvp
	Mvp[0][0] = -2*ve*omega_cl/RM - ve*ve*secl*secl/RN/RM;
	Mvp[0][1] = 0.0;
	Mvp[0][2] = -ve_rn2* ve*tl + vn_rm2*vd;
	Mvp[1][0] = 2*(-vd*omega_sl + vn*omega_cl)/RM + vn*ve*secl*secl/RN/RM;
	Mvp[1][1] = 0.0;
	Mvp[1][2] = ve_rn2 * (vd+vn*tl);
	Mvp[2][0] = 2*ve*omega_sl/RM;
	Mvp[2][1] = 0.0;
	double g = INS::get_gravity(phi, h);
	Mvp[2][2] = -ve_rn2*ve - vn_rm2*vn + 2*g/(sqrt(RM*RN)+h);
	// Mpv
	Mpv[0][0] = 1.0;
	Mpv[0][1] = 0.0;
	Mpv[0][2] = 0.0;
	Mpv[1][0] = 0.0;
	Mpv[1][1] = 1.0;
	Mpv[1][2] = 0.0;
	Mpv[2][0] = 0.0;
	Mpv[2][1] = 0.0;
	Mpv[2][2] = 1.0;
	// Mpp
	Mpp[0][0] = -vd / RM;
	Mpp[0][1] = 0.0;
	Mpp[0][2] = vn_rm;
	Mpp[1][0] = ve_rn * tl;
	Mpp[1][1] = -(vd + vn*tl)/RN;
	Mpp[1][2] = ve_rn;
	Mpp[2][0] = 0.0;
	Mpp[2][1] = 0.0;
	Mpp[2][2] = 0.0;

	// add K
	Vec omega_b = pre_info->dgm/pre_info->dt;
	double MKomega[3][3], MKf[3][3];
	for(int i = 0;i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			MKomega[i][j] = Cbn[i][j] * omega_b[j];
			MKf[i][j] = Cbn[i][j] * fb[j];
		}
	}

	// F
	for (int i = 0;i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			F(i  , j)   = Maa[i][j];
			F(i  , j+3) = Mav[i][j];
			F(i  , j+6) = Map[i][j];
			F(i+3, j)   = Mva[i][j];
			F(i+3, j+3) = Mvv[i][j];
			F(i+3, j+6) = Mvp[i][j];
			F(i+6, j)   = 0.0;
			F(i+6, j+3) = Mpv[i][j];
			F(i+6, j+6) = Mpp[i][j];
			F(i  , j+9) = -Cbn[i][j];
			F(i+3, j+12)=  Cbn[i][j];
			F(i   , j+15)= -MKomega[i][j];
			F(i+3,  j+18)= MKf[i][j];
		}
		F(i+9,  i+9) = -1/rel_t;
		F(i+12, i+12)= -1/rel_t;
		F(i+15, i+15)= -1/rel_t;
		F(i+18, i+18)= -1/rel_t;
	}

	// UKF to KF
	F = INS_INTERVAL * F;
	for(int i = 0; i<para_sum; i++)
		F(i, i) += 1.0;


	// G
	for(int i = 0;i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			G(i,    j)    = Cbn[i][j];
			G(i+3,  j+3)  = -Cbn[i][j];
		}
		G(i+9,  i+9)  = 1.0;
		G(i+12, i+12) = 1.0;
		G(i+15, i+15) = 1.0;
		G(i+18, i+18) = 1.0;
	}
	// check
	/*
	for( int i = 0; i< para_sum; i++)
	{
		for (int j = 0;j< para_sum; j++)
		{
			if(fabs(F(i, j)) > 10e-21)
				printf("%3d",10);
			else
				printf("%3d", 0);
		}
		cout<<endl;
	}*/
}

void update_HR_matrix(Vec att, Paras *init_para, GnssData *gnss_result, int iepoch)
{
	// arm correct
	double Cbn[3][3];
	double armk[3][3];
	INS::att2m(att, Cbn);
	INS::asym(init_para->antenna_arm, armk);
	dMatrix MCbn = dMatrix(3, 3);
	dMatrix Marm = dMatrix(3, 3);
	for(int i = 0;i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			MCbn(i, j) = Cbn[i][j];
			Marm(i, j) = armk[i][j];
		}
	}
	Marm = MCbn * Marm;
	for(int i = 0; i<3; i++)
	{
		for(int j = 0;j<3;j++)
		{
			H(i, j) = Marm(i, j);
		}
	}

	R(0, 0) = gnss_result->Nstd[iepoch] * gnss_result->Nstd[iepoch] * MRk;
	R(1, 1) = gnss_result->Estd[iepoch] * gnss_result->Estd[iepoch] * MRk;
	R(2, 2) = gnss_result->Ustd[iepoch] * gnss_result->Ustd[iepoch] * MRk;

}

void update_Z(Vec &gps_pos, Paras *init_paras, Info *cur_info)
{
	EarthPara eth;
	INS::earth_parameter(cur_info->pos[0], cur_info->pos[2], cur_info->vel, &eth);
	correct_arm(gps_pos, init_paras, cur_info);
	Z(0, 0) = (cur_info->pos[0] - gps_pos[0]) * (eth.RM + cur_info->pos[2]);
	Z(1, 0) = (cur_info->pos[1] - gps_pos[1]) * (eth.RN + cur_info->pos[2]) * cos(cur_info->pos[0]);
	Z(2, 0) = (cur_info->pos[2] - gps_pos[2]) * -1;
}

void kalman_filter(char mode)
{
	// one step prediction update
	if(mode == 'T')
	{
		// don't update X, as X is set to zero after each observation update
		dMatrix Qk = dMatrix(para_sum, para_sum);
		Qk = F * G * Q * ~G * ~F  + G * Q * ~G;
		for(int i = 0;i<para_sum;i++)
			for(int j = 0;j<para_sum;j++)
				Qk(i, j) = Qk(i, j) * INS_INTERVAL /2.0;
		DX0 = F * DX0 * ~F + Qk;
	}
	if(mode == 'M')
	{
		// X - att, vel, pos
		X = F * X;
		dMatrix K = dMatrix(3, 3);
		dMatrix V = dMatrix(3, 1);
		dMatrix tmp = dMatrix(3, 3);
		K = H * DX0 * ~H + R;
		matinv(K, 3);
		K = DX0 * ~H * K;
		V = Z - H * X;
		X = X + K * V;
		dMatrix I = dMatrix(para_sum, para_sum);
		for (int i = 0;i<para_sum; i++) I(i, i) = 1.0;
		DX0 = (I - K * H) * DX0;
	}
}

int get_gps_data(Vec &gps_pos, GnssData *gnss_result, double t)
{
	static double et;
	//match time
	for (int iepoch = 0;iepoch<gnss_result->nepoch;iepoch++)
	{
		if(et > t)
			return 0;
		et = gnss_result->gpstime[iepoch];
		if(fabs(et - t) < 1e-3)
		{
			gps_pos[0] = gnss_result->lat[iepoch];
			gps_pos[1] = gnss_result->lon[iepoch];
			gps_pos[2] = gnss_result->h[iepoch];
			return iepoch;
		}
	}
	return -1;
}
	
int get_imu_data(Info *pre_info, Info *cur_info, InsData *imu_data, InitState *state0, int iepoch, double dt)
{
	static bool first_flag = true;

	// imu state
	if(first_flag)
	{
		pre_info->att[0] = state0->roll;
		pre_info->att[1] = state0->pitch;
		pre_info->att[2] = state0->yaw;
		pre_info->vel[0] = state0->v_n;
		pre_info->vel[1] = state0->v_e;
		pre_info->vel[2] = state0->v_u;
		pre_info->pos[0] = state0->phi;
		pre_info->pos[1] = state0->lamda;
		pre_info->pos[2] = state0->h;
		cur_info->dt = dt;
		first_flag = false;
	}
	else
	{
		pre_info->att = cur_info->att;
		pre_info->vel = cur_info->vel;
		pre_info->pos = cur_info->pos;
	}

	// imu data
	double tr;
	if(fabs(dt - INS_INTERVAL) > 0.0001)
	{
		tr = dt / (imu_data->t[iepoch] - imu_data->t[iepoch-1]);
	}
	else
	{
		tr = 1.0;
	}
	pre_info->dam = cur_info->dam;
	pre_info->dgm = cur_info->dgm;
	for (int i = 0; i< 3; i++)
	{
		cur_info->dam[i]   = imu_data->da[iepoch*3    + i] * tr;
		cur_info->dgm[i]   = imu_data->dg[iepoch*3    + i] * tr;
		//pre_info->dam[i] = imu_data->da[iepoch*3 -3 + i] * tr;
		//pre_info->dam[i] = imu_data->dg[iepoch*3 -3 + i] * tr;
	}
	
	pre_info->dt = cur_info->dt;
	cur_info->dt = dt;
	return 1;
}

void correct_att(Vec datt, Vec &att)
{
	// rv2m
	static dMatrix MCbn = dMatrix(3, 3);
	static dMatrix MCnn = dMatrix(3, 3);
	double Cbn[3][3], Cnn[3][3];
	INS::att2m(att, Cbn);
	INS::rv2m(datt, Cnn);
	for(int i = 0; i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			MCbn(i, j) = Cbn[i][j];
			MCnn(i, j) = Cnn[i][j];
		}
	}
	matinv(MCnn, 3);
	MCbn = MCnn * MCbn;
	for(int i = 0; i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			Cbn[i][j] = MCbn(i, j);
		}
	}
	INS::m2att(Cbn, att);
}

void correct_arm(Vec &gps_pos, Paras *init_paras, Info *pre_info)
{	
	EarthPara eth;
	double Cbn[3][3];
	INS::earth_parameter(pre_info->pos[0], pre_info->pos[2], pre_info->vel, &eth);
	INS::att2m(pre_info->att, Cbn);
	Vec armb = init_paras->antenna_arm;
	Vec l_cor;

	double D_inv[3];
	D_inv[0] = 1.0/(eth.RM + pre_info->pos[2]);
	D_inv[1] = 1.0/(eth.RN + pre_info->pos[2])/cos(pre_info->pos[0]);
	D_inv[2] = -1;
	for(int i = 0;i<3;i++)
	{
		l_cor[i] = D_inv[i]* (Cbn[i][0]*armb[0] + Cbn[i][1]*armb[1] + Cbn[i][2]*armb[2]);
	}
	gps_pos[0] -= l_cor[0];
	gps_pos[1] -= l_cor[1];
	gps_pos[2] -= l_cor[2];
}

void correct_imu(Info * info, CorrImu * corr)
{
	for(int i = 0;i<3;i++)
	{
		info->dam[i] = (info->dam[i] - corr->dacc[i]*info->dt)/(1 + corr->Kacc[i]);
		info->dgm[i] = (info->dgm[i] - corr->dgym[i]*info->dt)/(1 + corr->Kgym[i]);
	}
}

void feedback(CorrImu *corr, Info *cur_info)
{
	Vec datt;
	for(int i = 0;i<3;i++)
	{
		datt[i]  = X(i,   0);
		cur_info->vel[i] -= X(i+3, 0);
		corr->dgym[i] += X(i+9, 0);
		corr->dacc[i] += X(i+12,0);
		corr->Kgym[i] += X(i+15,0);
		corr->Kacc[i] += X(i+18,0);
	}
	EarthPara eth;
	INS::earth_parameter(cur_info->pos[0], cur_info->pos[2], cur_info->vel, &eth);
	cur_info->pos[0] = cur_info->pos[0] - X(6, 0)/eth.RM;
	cur_info->pos[1] = cur_info->pos[1] - X(7, 0)/eth.RN/cos(cur_info->pos[0]);
	cur_info->pos[2] = cur_info->pos[2] - X(8, 0)/(-1);

	correct_att(datt, cur_info->att);
}

int integration(GnssData *gnss_result, InsData *imu_data, Paras *init_paras, InsState *integrate_result)
{
	InitState imu_state0;
	init_ins_state(init_paras, &imu_state0);
	if(INTERVALS == 0.0)
	{
		cout<<"gGPS interval is zero!"<<endl;
		return 0;
	}
	init_matrix(init_paras);

	int nepoch = imu_data->ncount;
	int s_epoch = 0;
	double last_time = init_paras->start_time, cur_time, dt;
	Vec gps_pos;
	Info pre_info, cur_info;
	CorrImu corr;
	for (int iepoch = 0; iepoch<nepoch; iepoch++)
	{
		// ins algorithm
		cur_time = imu_data->t[iepoch];
		if(cur_time <= init_paras->start_time)
			continue;
		// check if is integer time
		if (fmod(last_time + INS_INTERVAL, 1) > INS_INTERVAL || fabs(fmod(last_time, 1)) < 1e-5) // pure ins
		{
			dt = cur_time - last_time;
			if(!get_imu_data(&pre_info, &cur_info, imu_data, &imu_state0, iepoch, dt))
				continue;
			correct_imu(&cur_info, &corr);
			if(! INS::ins_update(cur_info.dam, pre_info.dam, cur_info.dgm, pre_info.dgm, 
				pre_info.att, pre_info.vel, pre_info.pos, cur_info.att, cur_info.vel, cur_info.pos, dt))
				return 0;
			// time update
			update_transform_matrix(init_paras->rel_t, &pre_info);
			kalman_filter('T');
		}
		else  // GPS KF
		{
			int gps_time = (int)(cur_time);
			cur_time = gps_time;
			dt = cur_time - last_time;
			// ins algorithm
			if(!get_imu_data(&pre_info, &cur_info, imu_data, &imu_state0, iepoch, dt))
				continue;
			correct_imu(&cur_info, &corr);
			if(! INS::ins_update(cur_info.dam, pre_info.dam, cur_info.dgm, pre_info.dgm, 
				pre_info.att, pre_info.vel, pre_info.pos, cur_info.att, cur_info.vel, cur_info.pos, dt))
				return 0;
			// time update
			update_transform_matrix(init_paras->rel_t, &pre_info);
			kalman_filter('T');
			// gps result
			int flag;
			if(!(flag = get_gps_data(gps_pos, gnss_result, gps_time*1.0)))
			{
				if(flag == -1)
					break;
				last_time = cur_time;
				iepoch--;
				continue;
			}

			// measurement update
			update_Z(gps_pos, init_paras, &cur_info);
			update_HR_matrix(cur_info.att, init_paras, gnss_result, flag);
			kalman_filter('M');
			// feedback
			feedback(&corr, &cur_info);
			// update time & X
			iepoch--;
			update_matrix(init_paras);
		}
		
		// result
		integrate_result->t[s_epoch] = cur_time;
		integrate_result->att[3*s_epoch + 0] = cur_info.att[0];
		integrate_result->att[3*s_epoch + 1] = cur_info.att[1];
		integrate_result->att[3*s_epoch + 2] = cur_info.att[2];
		integrate_result->vel[3*s_epoch + 0] = cur_info.vel[0];
		integrate_result->vel[3*s_epoch + 1] = cur_info.vel[1];
		integrate_result->vel[3*s_epoch + 2] = cur_info.vel[2];
		integrate_result->pos[3*s_epoch + 0] = cur_info.pos[0];
		integrate_result->pos[3*s_epoch + 1] = cur_info.pos[1];
		integrate_result->pos[3*s_epoch + 2] = cur_info.pos[2];
		s_epoch++;
		last_time = cur_time;
		// output process percentage
		if (s_epoch%10000 == 0)
			cout<<">";
	}
	integrate_result->nepoch = s_epoch;
	return 1;
}

void compare_result(InsState *real_state, InsState *com_result)
{
	FILE *fp, *fp0;
	if((fp = fopen("../../data/res.txt", "w")) == NULL)
		return;
	if((fp0 = fopen("../../data/result.txt", "w")) == NULL)
		return;

	// output residual
	double datt[3], dvel[3], dpos[3];
	double att_rms[3], vel_rms[3], pos_rms[3];
	for(int i = 0;i<3;i++)
	{
		att_rms[i] = 0.0;
		vel_rms[i] = 0.0;
		pos_rms[i] = 0.0;
	}
	int c_epoch0 = 0, r_epoch0 = 0, nepoch = 0;
	for(;c_epoch0 < com_result->nepoch;)
	{
		if(r_epoch0 >= real_state->nepoch)
			break;
		if(fabs(com_result->t[c_epoch0] - real_state->t[r_epoch0]) > INS_INTERVAL * 0.5)
		{
			if(com_result->t[c_epoch0] > real_state->t[r_epoch0])
			{
				r_epoch0++;
				continue;
			}
			else if(com_result->t[c_epoch0] < real_state->t[r_epoch0])
			{
				c_epoch0++;
				continue;
			}
		}
		if(fabs(com_result->t[c_epoch0] - real_state->t[r_epoch0]) > INS_INTERVAL * 0.5)
			cout<<"ERROR(compare_result)"<<endl;
		for(int i = 0;i<3;i++)
		{
			if(i < 2)
				dpos[i] = arc2degree(com_result->pos[c_epoch0*3 + i]) - arc2degree(real_state->pos[r_epoch0*3 + i]);
			else
				dpos[i] = com_result->pos[c_epoch0*3 + i] - real_state->pos[r_epoch0*3 + i];
			datt[i] = arc2degree(com_result->att[c_epoch0*3 + i]) - arc2degree(real_state->att[r_epoch0*3 + i]);
			dvel[i] = com_result->vel[c_epoch0*3 + i] - real_state->vel[r_epoch0*3 + i];
			datt[i] = norm_arc(datt[i]);
			att_rms[i] = att_rms[i] + datt[i] * datt[i];
			vel_rms[i] = vel_rms[i] + dvel[i] * dvel[i];
			pos_rms[i] = pos_rms[i] + dpos[i] * dpos[i];
		}
		fprintf(fp, "%14.5lf%14.8lf%14.8lf%10.3lf%10.3lf%10.3lf%10.3lf%14.8lf%14.8lf%14.8lf\n",
			com_result->t[c_epoch0], datt[0], datt[1], datt[2], dvel[0], dvel[1], dvel[2], dpos[0], dpos[1], dpos[2]);
		c_epoch0++;
		r_epoch0++;
		nepoch++;
	}
	for(int i = 0;i<3;i++)
	{
		att_rms[i] = sqrt(att_rms[i] / nepoch);
		vel_rms[i] = sqrt(vel_rms[i] / nepoch);
		pos_rms[i] = sqrt(pos_rms[i] / nepoch);
	}
	printf("att std(degree, degree, degree):\n %8.3lf%8.3lf%8.3lf\n", att_rms[0], att_rms[1], att_rms[2]);
	printf("vel std(m/s, m/s, m/s):\n %8.3lf%8.3lf%8.3lf\n", vel_rms[0], vel_rms[1], vel_rms[2]);
	printf("pos std(10-5 degree, 10-5 degree, m):\n %8.3lf%8.3lf%8.3lf\n", 1E5*pos_rms[0], 1E5*pos_rms[1], pos_rms[2]);

	// output result
	for(c_epoch0 = 0;c_epoch0 < com_result->nepoch;c_epoch0++)
	{
		if(com_result->t[c_epoch0] < 442336)
			continue;
		double att[3], vel[3], pos[3];
		for(int i = 0;i<3;i++)
		{
			att[i] = arc2degree(com_result->att[c_epoch0*3 + i]);
			vel[i] = com_result->vel[c_epoch0*3 + i];
			if(i<2)
				pos[i] = arc2degree(com_result->pos[c_epoch0*3 + i]);
			else
				pos[i] = com_result->pos[c_epoch0*3 + i];
		}
		fprintf(fp0, "%14.5lf%14.8lf%14.8lf%10.3lf%10.3lf%10.3lf%10.3lf%14.8lf%14.8lf%14.8lf\n",
			com_result->t[c_epoch0], pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], att[0], att[1], att[2]);
	}
	fclose(fp);
	fclose(fp0);
}

int algorithm()
{
	// initial
	Paras paras;
	init_parameters(&paras);
	string imufile = "../../data/d1_imu.bin";
	string gnssfile = "../../data/gps.txt";

	// read GPS result
	GnssData gnss_result;
	if(!read_gnssfile(gnssfile, &gnss_result))
		return 0;

	// read imu data and calculate
	InsData imu_data;
	if(!read_imufile(imufile, &imu_data))
		return 0;

	// GPS INU integration
	InsState integrate_result;
	cout<<endl<<"<<< begin to integrate GPS/IMU "<<endl;
	if(!integration(&gnss_result, &imu_data, &paras, &integrate_result))
		return 0;

	// result assessment
	cout<<"\n\nbegin to assess result"<<endl;
	InsState real_result;
	if(!read_result("../../data/d1_real.txt", &real_result))
		return 0;
	compare_result(&real_result, &integrate_result);

	return 1;
}

