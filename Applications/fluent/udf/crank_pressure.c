
#include "udf.h"

DEFINE_PROFILE(crankcase_pressure,t,i)
{
	real crankAngle,Initial_CA,rpm,Pa,c1,c2,c3,c4,c5,c6;
	face_t f;

	rpm = RP_Get_Real("dynamesh/in-cyn/crank-rpm");		
	Initial_CA = RP_Get_Real("dynamesh/in-cyn/crank-start-angle");	

	crankAngle = 6.0*rpm*CURRENT_TIME + Initial_CA;	/* crank angle in degrees */
	c1 = crankAngle;
	c2 = c1*crankAngle;
	c3 = c2*crankAngle;
	c4 = c3*crankAngle;
	c5 = c4*crankAngle;
	c6 = c5*crankAngle;

	Pa = 1.01325e5;
	
	if(c1 < 140.5){
		begin_f_loop(f,t)
                { 
                        F_PROFILE(f,t,i) = (-2.3009735653e-7*c6 + 9.6941073775e-5*c5 - 1.5626002197e-2*c4 + 1.1980887092*c3 - 4.2285233312e1*c2 + 7.8048089568e2*c1 - 3.7857526556e3);
                }
                end_f_loop(f,t)
	}else{
		 begin_f_loop(f,t)
                {
                        F_PROFILE(f,t,i) = (6.5698976597e-8*c6 - 9.2510994385e-5 *c5 + 5.3633104688e-2*c4 - 1.6378246757e1*c3 + 2.7789178757e3*c2 - 2.4885966632e5*c1 + 9.2356941243e6);
                }
                end_f_loop(f,t)
	}
}
