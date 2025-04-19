#ifndef _AER_VIIRS_EPS_MATCH_H_
#define _AER_VIIRS_EPS_MATCH_H_
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include "netcdf.h"
#include "hdf5.h"
#include "Match_Const.h"
using namespace std;

#ifndef EXTERN
#define EXTERN extern
#endif

/****** Data  Paths ******/
const string AERONET_PATH="../Data/AERONET/";
//const string AOD_PATH="/data/data314/hliu/DATA/VIIRS_NASA/AERDT_L2_VIIRS_SNPP/";  // YYYY/DDD/AERDT_L2_VIIRS_SNPP.A2012333.2306.011.2020200031009.nc
//const string MATCH_PATH="../Data/MATCHUP_DT/";
//const string AOD_PATH="/data/data314/hliu/DATA/VIIRS_NASA/AERDB_L2_VIIRS_SNPP/";  // YYYY/DDD/AERDB_L2_VIIRS_SNPP.A2012333.2306.011.2020200031009.nc
//const string MATCH_PATH="../Data/MATCHUP_DB_SNPP/";
const string AOD_PATH="/data/data314/hliu/DATA/VIIRS_NASA/AERDB_L2_VIIRS_NOAA20/";  // YYYY/DDD/AERDB_L2_VIIRS_NOAA20.A2023004.2124.002.2023080160611.nc
const string MATCH_PATH="../Data/MATCHUP_DB_NOAA20/";

/****** AERONET  Constants (V3) *******/
const int NUM_AERONET_WAVELENGTH = 22;
const float AERONET_WAVELENGTH[NUM_AERONET_WAVELENGTH] = 
              {1.640, 1.020, 0.870, 0.865, 0.779, 0.675, 0.667, 0.620, 
               0.560, 0.555, 0.551, 0.532, 0.531, 0.510, 0.500, 0.490, 
               0.443, 0.440, 0.412, 0.400, 0.380, 0.340};



/****** Matching  Constants ******/
const int NVLDAER = 2;  /* minimum number of valid AERONET measurements used for averaging */
const float MATCH_TIME_WINDOW = 1.0;  /* AERONET matching time window in hour (centered on satellite overpass time) */ 
const float MATCH_RADIUS = 27.5;      /* Satellite retrievals matching spatial domain (radius [km] of circle centered on station)  */ 
const int MIN_N_MATCH_PIX = 2;      /* Minimum number of pixels within matching domain */
const float MAXLAT = 80.;

/****** Data Structures for AERONET *******/

typedef struct {  
   float time;  /* fractional hour (0-1)*/
   float aods[NUM_AERONET_WAVELENGTH];
} AerAod;

typedef struct {  
   string  staName; 
   float   staLon;
   float   staLat;
   int     nmeas;   /* number measurements */
   AerAod* meas;    /* measurments of AODs at multiple wavelengths */
} OneDayAerAod;
  

/****** Data Structures for Match-ups *******/

typedef struct {  
   int   nMeas;     /* number of AERONET measurements */
   int   nPixs;     /* number of AOD retrieval pixels */
   float satTime;   /* fractional hour [0-24] */
   
   // AERONET data
   char   staName[32]; /* station info */
   float  staLon;
   float  staLat;  
   float  *meas;     /* (NUM_AERONET_WAVELENGTH+1)*nmeas
                        time (of day [0-1] + AOD at 22 wavelengths */
   
   // AOD outputs (nPixs)
   float *lon;   
   float *lat;
   float *solzen;
   float *satzen;
   float *relazi;  // Gordon convention: ABS(!PI-relazi)
   float *sctang;
   float *aod550;  // best estimate (QF>1 moderate or good)
   float *ae;      // 0.41/0.48 over arid land; 0.48/0.69 over vegetated/mixed land; 0.55/0.865 over ocean
   Int16 *lndSea;  // 0:Ocean 1:Land 

} MatchupRecord;


/****** Function Prototypes ******/
string int_to_str(int num);
int calendar_date (int yyyyddd);
int dayofyear (int yyyymmdd);
int getAerInfo(map<string,float*>& aerInfo);
int getDailyAerAod(map<string,float*>& aerInfo, string yyyyddd, vector<OneDayAerAod>& aerData);
int getVIIRSgranInfo(string yyyymmdd, map<float,string>&granInfo);
int matchup(string yyyymmdd, vector<OneDayAerAod>& aerData, map<float,string>& granInfo, 
            vector<MatchupRecord>& matchData);

void freeMemForMatch(MatchupRecord mr);

#endif /*_AER_VIIRS_EPS_MATCH_H_*/
