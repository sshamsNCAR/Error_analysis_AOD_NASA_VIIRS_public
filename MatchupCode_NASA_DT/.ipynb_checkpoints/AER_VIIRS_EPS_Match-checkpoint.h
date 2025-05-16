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
const string AERONET_PATH="/glade/campaign/acom/acom-da/SERVIR/ind-obs/aeronet/AOD_Level15_All_Points_V3/for/";
//const string AOD_PATH="/data/data314/hliu/DATA/VIIRS_NASA/AERDT_L2_VIIRS_SNPP/";  // YYYY/DDD/AERDT_L2_VIIRS_SNPP.A2012333.2306.011.2020200031009.nc
//const string MATCH_PATH="../Data/MATCHUP_DT/";
const string AOD_PATH="/glade/campaign/acom/acom-da/SERVIR/VIIRS/python_download_buffer/data-ingest/viirs_data/";  // YYYY/DDD/AERDT_L2_VIIRS_NOAA20.A2024003.0824.002.2024003210456.nc
// const string MATCH_PATH="/glade/u/home/sshams/code/AFRICA_SERVIR/Error_analysis_AOD_NASA/Data/MATCHUP_DT/";
const string MATCH_PATH="/glade/campaign/acom/acom-da/SERVIR/match_aeronet_viirs/matched_DT_1hours_TW/"; // we are using differnet time window for matching and save the data in different folder



/****** AERONET  Constants (V3) *******/
const int NUM_AERONET_WAVELENGTH = 22;
const float AERONET_WAVELENGTH[NUM_AERONET_WAVELENGTH] = 
              {1.640, 1.020, 0.870, 0.865, 0.779, 0.675, 0.667, 0.620, 
               0.560, 0.555, 0.551, 0.532, 0.531, 0.510, 0.500, 0.490, 
               0.443, 0.440, 0.412, 0.400, 0.380, 0.340};

 // Constants for the required wavelengths for 550 nm interpolation
const int INDEX_440 = 17;  // 440 nm is the 18th element (0-based index is 17)
const int INDEX_500 = 14;  // 500 nm is the 15th element (0-based index is 14)
const int INDEX_675 = 5;   // 675 nm is the 6th element (0-based index is 5)              

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
   float aod550;  // interpoalted AOD at 0.55um */
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
   float aerMean550;  // Mean AERONET AOD at 550 nm
   float aerStd550;   // Standard deviation of AERONET AOD at 550 nm
   
   // AOD outputs (nPixs)
   float *lon;   
   float *lat;
   float *solzen;
   float *solazi;
   float *satzen;
   float *satazi;
   float *aod550; // AOD at 550 nm all the pixels as array output
   float satMean550;  // Mean satellite AOD at 550 nm
   float satStd550;   // Standard deviation of satellite AOD at 550 nm
   float *ae1;    // 0.55 vs 0.86um
   float *ae2;    // 0.86 vs 2.13um
   Int16 *qf;     // 0: Bad 1:Marginal 2:Good 3:Very Good
   Int16 *lndSea; // 0:Ocean 1:Land and ephemeral water 2:Coastal

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
