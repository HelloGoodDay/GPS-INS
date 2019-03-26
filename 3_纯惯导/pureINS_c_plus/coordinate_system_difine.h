
#ifndef __GNSSCRDSTRUCTDEF__H
#define __GNSSCRDSTRUCTDEF__H


#ifndef PI
#define PI 3.1415926535898
#endif

#define WGS84_f  1/298.257223563 //CGCS2000_f  1/298.257222101
#define WGS84_a  6378137.0       //CGCS2000_a  6378137.0


typedef struct tagCRDCARTESIAN	//笛卡尔坐标 Coordinate Cartesian 
{
	double	x;
	double	y;
	double	z;
	tagCRDCARTESIAN()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
} CRDCARTESIAN, *PCRDCARTESIAN;



typedef struct tagCRDGEODETIC  //大地坐标 Coordinate Geodetic
{
	double	latitude;	//纬度
	double	longitude;	//经度
	double	height;		//高度
	tagCRDGEODETIC()
	{
		latitude = 0.0;
		longitude = 0.0;
		height = 0.0;
	}
} CRDGEODETIC, *PCRDGEODETIC;



typedef	struct tagCRDTOPOCENTRIC //站心地平坐标(线坐标形式) Coordinate TOPOCentric
{
	double	northing;
	double	easting;
	double	upping;
	tagCRDTOPOCENTRIC()
	{
		northing = 0.0;
		easting = 0.0;
		upping = 0.0;
	}
} CRDTOPOCENTRIC, *PCRDTOPOCENTRIC;


typedef	struct tagCRDTOPOCENTRICPOLAR	//站心地平坐标(极坐标形式) Coordinate TOPOCentric Polar
{
	double	range;			//距离
	double	azimuth;		//方位角
	double	elevation;		//高度角 = 90 - zenith distance
	tagCRDTOPOCENTRICPOLAR()
	{
		range = 0.0;
		azimuth = 0.0;
		elevation = 0.0;
	}
} CRDTOPOCENTRICPOLAR, *PCRDTOPOCENTRICPOLAR;



#endif