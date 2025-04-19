#ifndef _MATCH_CONST_H_
#define _MATCH_CONST_H_
#include <string>
#include <cmath>
using namespace std;

typedef unsigned char       UInt8;
typedef unsigned short int  UInt16;
typedef short int  Int16;

const int PROC_SUCCEED = 0;
const int PROC_FAIL = 1;
const int MAX_DIM_SIZE = 6; /* assume maximum dimension size as 4 */

const unsigned short int   VIIRS_MISVAL_MIN_INT = 65528; 
const float VIIRS_MISVAL_MAX_FLT = -999.2; 


/****** Numerical  Constants ******/
const float PI = acos(-1.);
const float DEG2RAD = PI/180.;
const float RAD2DEG = 180./PI;
const double EARTH_RADIUS = 6371.2;     /* mean earth radius in km */

const double MOIST_AIR_LAPSE_RATE = 6.5/1000;
const double GRAVITY = 9.80665;
const float  DRYGAS = 287.05;

const float MISVAL_FLOAT = -999.9;
const unsigned char  MISVAL_UINT8 = 255;
const short MISVAL_INT16 = -9999;

#endif /*_MATCH_CONST_H_*/
