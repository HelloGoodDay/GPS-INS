
#ifndef __GNSSCRDCLASS__H
#define __GNSSCRDCLASS__H

#include "coordinate_system_difine.h"
#include <math.h>


class CGNSSCRDClass  
{
public:

	void	RadianToDuFenMiao(double dRadian, int &dDeg, int &dMin, double &dSec);
	void	CartesianToGeodetic(CRDCARTESIAN cc, CRDGEODETIC &cg);//迪卡尔坐标转大地坐标
	void	GeodeticToCartesian(CRDGEODETIC cg, CRDCARTESIAN &cc);//大地坐标转迪卡尔坐标
	void	CartesianToTopocentric(CRDTOPOCENTRIC &ct, CRDCARTESIAN cc, CRDCARTESIAN ccCenter);//迪卡尔坐标转站心坐标
	void	TopocentricToTopocentricPolar(CRDTOPOCENTRIC ct, CRDTOPOCENTRICPOLAR &ctp);//站心坐标转站心极坐标
	double a;
	double e2;
	CGNSSCRDClass();
	virtual ~CGNSSCRDClass();
	
};

#endif