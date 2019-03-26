#include "pro.h"

int read_resultfile(string filename, InsState *ins_state)
{
	FILE *fp;
	// 以"r"模式打开时，读到0x1A会意外终止
	if((fp = fopen(filename.c_str(), "rb")) == NULL)
	{
		cout<<"can't open file "<<filename.c_str()<<endl;
		return 0;
	}

	int ncount = 0;
	double buf[10];
	while(!feof (fp))
	{
		try
		{
			fread(&buf, sizeof(double), 10, fp);
		}
		catch(...)
		{
			cout<<"file format error!"<<endl;
			return 0;
		}
		ins_state->t[ncount] = buf[0];
		for (int i = 0;i<3;i++)
		{
			if(i==2)
				ins_state->pos[ncount*3 + i] = buf[i + 1];
			else
				ins_state->pos[ncount*3 + i] = degree2arc(buf[i + 1]);
			ins_state->vel[ncount*3 + i] = buf[i + 4];
			ins_state->att[ncount*3 + i] = degree2arc(buf[i + 7]);
		}
		ncount++;
	}
	ins_state->nepoch = ncount;
	fclose(fp);
	cout<<"read file successfully, the epoch number is "<<ncount<<endl;
	return 1;
	return 2;
}

void compare_result(InsState *test, InsState *real, InsState *res, double rms[9])
{
	int jepoch = 0;
	for (int i = 0; i< 9;i++) rms[i] = 0.0;
	for( int iepoch = 0; iepoch<real->nepoch; )
	{
		if(test->t[jepoch] < real->t[iepoch])
		{
			jepoch++;
			continue;
		}
		if(test->t[jepoch] != real->t[iepoch])
		{
			printf("epoch match error!");
		}

		double att_res, vel_res, pos_res;
		for (int i=0; i<3; i++)
		{
			att_res = test->att[jepoch*3 + i] - real->att[iepoch*3 + i];
			vel_res = test->vel[jepoch*3 + i] - real->vel[iepoch*3 + i];
			pos_res = test->pos[jepoch*3 + i] - real->pos[iepoch*3 + i];
			
			att_res = arc2degree(att_res);
			if(i<2) pos_res = arc2degree(pos_res);

			res->att[iepoch*3 + i] = att_res;
			res->vel[iepoch*3 + i] = vel_res;
			res->pos[iepoch*3 + i] = pos_res;
			res->t[iepoch] = real->t[iepoch];

			rms[i  ] += (att_res * att_res);
			rms[i+3] += (vel_res * vel_res);
			rms[i+6] += (pos_res * pos_res);
		}
		jepoch++;
		iepoch++;
	}
	for(int i = 0;i<9;i++)
	{
		rms[i] = sqrt(rms[i] / real->nepoch);
	}
	printf("residual of attitude is (degree):\n");
	printf("%10.8lf  %10.8lf  %10.8lf\n", rms[0], rms[1], rms[2]);
	printf("residual of velocity is (m/s):\n");
	printf("%10.8lf  %10.8lf  %10.8lf\n", rms[3], rms[4], rms[5]);
	printf("residual of position is (degree, degree, m):\n");
	printf("%10.8lf  %10.8lf  %10.8lf\n", rms[6], rms[7], rms[8]);
}

void compute_residual(InsState *res, int endepoch, string outfile)
{
	FILE *fp;
	if((fp = fopen(outfile.c_str(), "w")) == NULL)
	{
		printf("can't open output file!\n");
		return;
	}
	for(int iepoch = 0;iepoch <endepoch; iepoch++)
	{
		fprintf(fp, "%12.3lf %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %7.4lf\n", res->t[iepoch],
			res->att[iepoch*3 + 0]*1e5, res->att[iepoch*3 + 1]*1e5, res->att[iepoch*3 + 2]*1e5,
			res->vel[iepoch*3 + 0]*1e5, res->vel[iepoch*3 + 1]*1e5, res->vel[iepoch*3 + 2]*1e5,
			res->pos[iepoch*3 + 0]*1e5, res->pos[iepoch*3 + 1]*1e5, res->pos[iepoch*3 + 2]*1e5);
	}
	fclose(fp);
}

int ins_process(string filename, string filename2)
{
	InsData *ins_data;
	ins_data = new InsData();
	// read ins file
	if(!INS::read_file(filename, ins_data))
	{
		return 0;
	}
	// ion algorithm
	InitState state0;
	InsState ins_state;
	if(!INS::ins_algorithm(state0, ins_data, &ins_state))
	{
		return 0;
	}
	delete ins_data;

	// read result file
	InsState ins_state_r;
	if(!read_resultfile(filename2, &ins_state_r))
	{
		return 0;
	}
	//compare result
	InsState res;
	double rms[9];
	compare_result(&ins_state, &ins_state_r, &res, rms);

	// output results
	string outfile = "../../residual.txt";
	compute_residual(&res, ins_state_r.nepoch, outfile);
	return 1;
}