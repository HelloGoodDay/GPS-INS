//#include "pro.h"
#include "vec.h"
#include <string.h>


int main()
{
	char filename1[FILEMAX] = "../../data/Data1.bin";
	char filename2[FILEMAX] = "../../data/Data1_PureINS.bin";

	printf("注意读文件时，读入的应该是角速度计和陀螺仪的增量\n");
	printf("检查时间间隔是否正确\n");

	if(!ins_process(filename1, filename2, 1))
	{
		printf("error!\n");
		return 0;
	}
	
	return 1;

}