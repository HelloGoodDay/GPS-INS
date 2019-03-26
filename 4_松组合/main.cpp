#include "process.h"
#include "time.h"

int main()
{
	double a = EARTH_e;
	clock_t start,finish;
	start = clock();

	algorithm();
	
	finish = clock();
	double totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"\nusing time: "<<totaltime<<" second"<<endl;
		return 1;}