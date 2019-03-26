function [ g ] = get_gravity( B, h )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    a = 6378137.;
	b = 6356752.3141;
	f = 1 - b / a  ;
    EARTH_V = 7.292115E-5;
    GM = 3.986005e+14;
	m = EARTH_V*EARTH_V*a*a*b/GM;
	ga = 9.7803267715;
	gb = 9.8321863685;
	gl = (a*ga*cos(B)*cos(B) + b*gb*sin(B)*sin(B))/ ...
		sqrt(a*a*cos(B)*cos(B) + b*b*sin(B)*sin(B));
	g = gl*(1 - 2*h/a*(1+f+m-2*f*sin(B)*sin(B)) + 3*h*h/a/a);

end

