// GNSSCRDClass.cpp: implementation of the CGNSSCRDClass class.
//
//////////////////////////////////////////////////////////////////////
#include "coordinate_transformation.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CGNSSCRDClass::CGNSSCRDClass()
{
	a = WGS84_a;
	double f = WGS84_f;
	e2 = 2*f - f*f;
}

CGNSSCRDClass::~CGNSSCRDClass()
{

}

/************************************************************************/
/*	由大地坐标转换为笛卡尔坐标
说明：
cg：指向待转换的大地坐标的指针；
pcc：指向所转换出的笛卡尔坐标的指针；
dSemiMajorAxis：参考椭球的长半轴；
dFlattening：参考椭球的扁率。											*/
/************************************************************************/
void CGNSSCRDClass::GeodeticToCartesian (CRDGEODETIC cg, CRDCARTESIAN &cc)
{
	double B = cg.latitude;
	double L = cg.longitude;
	double H = cg.height;

	double N = a/sqrt(1-e2*sin(B)*sin(B));

	cc.x = (N + H)*cos(B)*cos(L);
	cc.y = (N + H)*cos(B)*sin(L);
	cc.z = (N*(1-e2) + H)*sin(B);
}

/************************************************************************/
/* 由笛卡尔坐标转换为大地坐标
说明：
cc：指向待转换的笛卡尔坐标的指针；
pcg：指向所转换出的大地坐标的指针；
dSemiMajorAxis：参考椭球的长半轴；
dFlattening：参考椭球的扁率。											*/
/************************************************************************/
void CGNSSCRDClass::CartesianToGeodetic (CRDCARTESIAN cc, CRDGEODETIC &cg)
{
//	double a = dSemiMajorAxis;
//	double f = dFlattening;
//	e2 = 2*f - f*f;
	
	cg.longitude = atan2(cc.y, cc.x);
	cg.latitude  = atan2(cc.z, (cc.x*cc.x + cc.y*cc.y));

	double tanB = 0.0;
	double tmp = 0.0;
	int i=0;
	do
	{
		i++;
		tmp = cg.latitude;
		tanB = (cc.z + a*e2*sin(tmp)/sqrt(1-e2*sin(tmp)*sin(tmp))) / sqrt(cc.x*cc.x + cc.y*cc.y);
		cg.latitude = atan(tanB);
	}
	while(fabs((cg.latitude-tmp)/PI*180*3600)>1.0e-10 && i<50);
    cg.height = sqrt(cc.x*cc.x + cc.y*cc.y)/cos(cg.latitude) -a/sqrt(1-e2*sin(cg.latitude)*sin(cg.latitude));
}


/************************************************************************/
/*由笛卡尔坐标转换为站心地平坐标
说明：
pct：指向所转换出的站心地平坐标的指针；
pcc：指向待转换的笛卡尔坐标的指针；
pccCenter：指向站心的笛卡尔坐标的指针；
dSemiMajorAxis：参考椭球的长半轴；
dFlattening：参考椭球的扁率。                                           */
/************************************************************************/
void CGNSSCRDClass::CartesianToTopocentric(CRDTOPOCENTRIC &ct, 
							 CRDCARTESIAN cc, 
							 CRDCARTESIAN ccCenter)
{
	double Xab = cc.x - ccCenter.x;
	double Yab = cc.y - ccCenter.y;
	double Zab = cc.z - ccCenter.z;

	CRDGEODETIC cgCenter;
	CartesianToGeodetic(ccCenter, cgCenter);

	ct.northing	=	-sin(cgCenter.latitude)*cos(cgCenter.longitude)*Xab
						-sin(cgCenter.latitude)*sin(cgCenter.longitude)*Yab
						+cos(cgCenter.latitude)*Zab;
	ct.easting	=	-sin(cgCenter.longitude)*Xab
						+cos(cgCenter.longitude)*Yab;
	ct.upping		=	 cos(cgCenter.latitude)*cos(cgCenter.longitude)*Xab
						+cos(cgCenter.latitude)*sin(cgCenter.longitude)*Yab
						+sin(cgCenter.latitude)*Zab;

}

/************************************************************************/
/* 由站心地平直角坐标转换为站心地平极坐标
说明：
pct：指向待转换的站心地平坐标的指针；     
pctp：指向所转换出的站心地平极坐标的指针；                              */
/************************************************************************/
void CGNSSCRDClass::TopocentricToTopocentricPolar(CRDTOPOCENTRIC ct, CRDTOPOCENTRICPOLAR &ctp)
{
	ctp.range	  = sqrt( ct.easting*ct.easting + ct.northing*ct.northing + ct.upping*ct.upping );
	ctp.azimuth	  = atan( ct.easting/ct.northing );
	if (ct.easting>=0)
	{
		if (ct.northing<0)
			ctp.azimuth+=PI;
	}		
	else
	{
		if (ct.northing<0)
			ctp.azimuth+=PI;
		else
			ctp.azimuth+=2*PI;
	}		
	ctp.elevation = asin( ct.upping/ctp.range );
}

void CGNSSCRDClass::RadianToDuFenMiao(double dRadian, int &dDeg, int &dMin, double &dSec)
{
	double dTmp = dRadian*180/PI;

	dDeg = int(dTmp);
	dTmp = (dTmp - dDeg)*60.0;
	dMin = int(dTmp);
	dSec = (dTmp - dMin)*60.0;

}

