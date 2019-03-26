#include "pro.h"

int read_resultfile(char* filename, InsState *ins_state)
{
	FILE *fp;
	int ncount = 0;
	double buf[10];
	int i;
	// 以"r"模式打开时，读到0x1A会意外终止
	if((fp = fopen(filename, "rb")) == NULL)
	{
		printf("can't open file %s\n", filename);
		return 0;
	}

	while(fread(&buf, sizeof(double), 10, fp))
	{
		ins_state->t = buf[0];
		for (i = 0;i<3;i++)
		{
			if(i==2)
				ins_state->pos[i] = buf[i + 1];
			else
				ins_state->pos[i] = degree2arc(buf[i + 1]);
			ins_state->vel[i] = buf[i + 4];
			ins_state->att[i] = degree2arc(buf[i + 7]);
		}
		ncount++;
		ins_state++;
		ins_state->nepoch = ncount;
	}
	fclose(fp);
	printf("read file successfully, the epoch number is %d\n", ncount);
	return 1;
}

void compare_result(InsState *test, InsState *real, InsState *res, double rms[9])
{
	int iepoch, jepoch;
	int i;
	double att_res, vel_res, pos_res;
	iepoch = jepoch = 0;

	for (i = 0; i< 9;i++) rms[i] = 0.0;
	while((test+iepoch)->nepoch > 0 || (test+iepoch)->t < (real+jepoch)->t)
	{
		if((test+iepoch)->t < (real+jepoch)->t)
		{
			iepoch++;
			continue;
		}
		if((test+iepoch)->t != (real+jepoch)->t)
		{
			printf("epoch match error!");
		}
		for (i = 0; i<3; i++)
		{
			att_res = (test+iepoch)->att[i] - (real+jepoch)->att[i];
			vel_res = (test+iepoch)->vel[i] - (real+jepoch)->vel[i];
			pos_res = (test+iepoch)->pos[i] - (real+jepoch)->pos[i];
			
			att_res = arc2degree(att_res);
			if(i<2) pos_res = arc2degree(pos_res);

			(res+jepoch)->att[i] = att_res;
			(res+jepoch)->vel[i] = vel_res;
			(res+jepoch)->pos[i] = pos_res;
			(res+jepoch)->t = (real+jepoch)->t;

			rms[i  ] += (att_res * att_res);
			rms[i+3] += (vel_res * vel_res);
			rms[i+6] += (pos_res * pos_res);
		}
		jepoch++;
		iepoch++;
	}
	for(i = 0;i<9;i++)
	{
		rms[i] = sqrt(rms[i] / jepoch);
	}
	printf("residual of attitude is (degree):\n");
	printf("%10.8lf  %10.8lf  %10.8lf\n", rms[0], rms[1], rms[2]);
	printf("residual of velocity is (m/s):\n");
	printf("%10.8lf  %10.8lf  %10.8lf\n", rms[3], rms[4], rms[5]);
	printf("residual of position is (degree, degree, m):\n");
	printf("%10.8lf  %10.8lf  %10.8lf\n", rms[6], rms[7], rms[8]);
}

void compute_residual(InsState *res, char* outfile, double acc)
{
	FILE *fp;
	int iepoch = 0;
	if((fp = fopen(outfile, "w")) == NULL)
	{
		printf("can't open output file!\n");
		return;
	}
	while((res+iepoch)->t > 0)
	{
		if((res+iepoch)->t >= 92281)
			break;
		fprintf(fp, "%12.3lf %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf\n", (res+iepoch)->t,
			(res+iepoch)->att[0]*acc, (res+iepoch)->att[1]*acc, (res+iepoch)->att[2]*acc,
			(res+iepoch)->vel[0]*acc, (res+iepoch)->vel[1]*acc, (res+iepoch)->vel[2]*acc,
			(res+iepoch)->pos[0]*acc, (res+iepoch)->pos[1]*acc, (res+iepoch)->pos[2]*acc);
		iepoch++;
	}
	fclose(fp);
}

void compute_result(InsState *test)
{
	char *outfile = "../../result.txt";
	FILE *fp;
	int iepoch = 0;
	if((fp = fopen(outfile, "w")) == NULL)
	{
		printf("can't open output file!\n");
		return;
	}
	while((test+iepoch)->t > 0)
	{
		fprintf(fp, "%12.3lf %12.7lf  %12.7lf  %12.7lf  %12.7lf  %12.7lf  %12.7lf  %12.7lf  %12.7lf  %12.7lf\n", (test+iepoch)->t,
			(test+iepoch)->att[0], (test+iepoch)->att[1], (test+iepoch)->att[2],
			(test+iepoch)->vel[0], (test+iepoch)->vel[1], (test+iepoch)->vel[2],
			(test+iepoch)->pos[0], (test+iepoch)->pos[1], (test+iepoch)->pos[2]);
		iepoch++;
	}
	fclose(fp);
}

int ins_process(char* filename, char* filename2, int test_flag)
{
	InitState state0;
	InsData *ins_data = (InsData*)malloc(sizeof(InsData) * DATASUM);
	InsState *ins_state = (InsState*)malloc(sizeof(InsState) * DATASUM);

	InsState *ins_state_r = (InsState*)malloc(sizeof(InsState) * DATASUM);
	InsState *res = (InsState*)malloc(sizeof(InsState) * DATASUM);

	double rms[9];

	int nepoch;
	char *outfile = "../../residual.txt";

	memset(ins_data, 0, sizeof(InsData) * DATASUM);
	memset(ins_state, 0, sizeof(InsState) * DATASUM);
	memset(ins_state_r, 0, sizeof(InsState) * DATASUM);
	memset(res, 0, sizeof(InsState) * DATASUM);

	test_flag = 0;
	// ****************原来的数据做测试*******************************
	if(test_flag == 1) 
	{
		init_state(&state0);
		//state0.roll -= degree2arc(0.2);
		//state0.yaw -=degree2arc(2);

		// read ins file
		if(!read_file_2(filename, ins_data))
		{
			return 0;
		}
		// ion algorithm
		if(!ins_algorithm(state0, ins_data, ins_state))
		{
			return 0;
		}
		compute_result(ins_state);
		free(ins_data);

		// read result file
		if(!read_resultfile(filename2, ins_state_r))
		{
			return 0;
		}

		//compare result
		compare_result(ins_state, ins_state_r, res, rms);

		// output residual
		// 为了提高精度，残差是以1e5输出的
		compute_residual(res, outfile, 1);
	}

	// *********************要计算的数据************************************
	if(test_flag == 0)
	{
		char filename_0[FILEMAX] = "../../data/imu_1208.txt";
		init_state(&state0);
		state0.yaw = degree2arc(77.9240439610316);
		state0.pitch = degree2arc(0.000686740684620727);
		state0.roll = degree2arc(0.00322623712156884);

		//state0.yaw = 0.0;
		//state0.pitch = 0.0;
		//state0.roll = 0.0;

		state0.phi = degree2arc(30.527817162);
		state0.lamda = degree2arc(114.356736560);
		state0.h = 76.3604;
		state0.t = 15060.0;

		// read ins file
		if(!read_file(filename_0, ins_data))
		{
			return 0;
		}
		// ion algorithm
		if(!ins_algorithm(state0, ins_data, ins_state))
		{
			return 0;
		}
		compute_result(ins_state);
		free(ins_data);

		// output residual
		// 为了提高精度，残差是以1e5输出的
		//compute_residual(res, outfile, 1);
	}

	return 1;
}