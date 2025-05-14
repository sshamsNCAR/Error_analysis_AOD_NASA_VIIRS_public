#define EXTERN
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iomanip>
#include "AER_VIIRS_EPS_Match.h"
using namespace std;
  

/*************************** M A I N   P R O G R A M **************************/ 
/*
  Main code used to match VIIRS EPS AOD and ground AERONET measurements
  
  Command line argument: year month day
  
  Example:
     match_viirs_dt 2013 1 1
*/

int main(int argc, char** argv)
{    
   int status;
   
   /****** get command line arguments ******/
   if (argc != 4) {
      cerr << "match_main.cpp: Wrong number of arguments!" << endl;
      return 1;
   }
   stringstream ss;
   ss << std::setw(4) << std::setfill('0') << atoi(argv[1]);
   string year = ss.str();
   ss.str("");
   ss << std::setw(2) << std::setfill('0') << atoi(argv[2]);
   string month = ss.str();
   ss.str("");
   ss << std::setw(2) << std::setfill('0') << atoi(argv[3]);
   string day = ss.str();
   ss.str("");
   
   string yyyymmdd = year+month+day;
   string yyyyddd = int_to_str(dayofyear( atoi(yyyymmdd.c_str()) ));
 
   
   /****** get AERONET station info ******/
   map<string,float*> aerInfo;  /* staName; lon/lat */
   status = getAerInfo(aerInfo);
   if (status == PROC_FAIL) {
      cout << "main(): Cannot get AERONET station information" << endl;
      return 1;
   }
   
   /****** get daily AERONET data ******/
   vector<OneDayAerAod> aerData;
   status = getDailyAerAod(aerInfo, yyyyddd, aerData);
   if (status == PROC_FAIL) {
      cout << "main(): Cannot get AERONET data" << endl;
      return 1;
   }
   
   /****** get VIIRS daytime granule time/domain ******/
   map<float, string> granInfo;  /* key: granule midtime; value: GMTCO filename */
   status = getVIIRSgranInfo(yyyyddd, granInfo);
   if (status == PROC_FAIL) {
      cout << "main(): Cannot get granule info" << endl;
      return 1;
   }
   
   /****** match up ******/
   vector<MatchupRecord> matchData;
   status = matchup(yyyymmdd, aerData, granInfo, matchData);
   if (status == PROC_FAIL) {
      cout << "main(): Cannot find match-ups" << endl;
      return 1;
   }
   
   /****** output the match-ups ******/
   string outFile = MATCH_PATH+"Match_DB_"+yyyyddd+".csv";
   // ofstream output(outFile.c_str(), ios::out | ios::binary); 
   int nMatch = matchData.size();
   cout << "Writing match-up results to: " << outFile << endl;
   ofstream output(outFile.c_str());
   if (!output.is_open()) {
      cerr << "Error: Unable to open file " << outFile << endl;
      return 1;
   }
   cout << "*******************************************" << outFile << endl;
   
   // Write CSV header
   output << "date,nMeas,nPixs,satTime,staName,staLon,staLat,aerMean550,aerSTD550,satMean550,satSTD550,aod550,qf\n";
   // Write match-up data
   for (vector<MatchupRecord>::iterator im = matchData.begin(); im != matchData.end(); ++im) {
      // creating array to have all the single aod of sat pixels to be in the csv file
      std::ostringstream aod550_values, qf_values;
      // Create an empty string stream

      // Append all values in the aod550 array
      for (int i = 0; i < im->nPixs; ++i) {
         aod550_values << im->aod550[i];  // Add each AOD value
         if (i < im->nPixs - 1) {
            aod550_values << ";";  // Separate values with a semicolon
         }

            // Append qf values
         qf_values << im->qf[i];
         if (i < im->nPixs - 1) {
               qf_values << ";";  // Separate values with a semicolon
         }

         
      }
      output << yyyymmdd << ","  // Add the date
            << im->nMeas << ","  // Number of AERONET measurements
            << im->nPixs << ","  // Number of satellite pixels
            << im->satTime << ","  // Satellite overpass time
            << im->staName << ","  // AERONET station name
            << im->staLon << "," 
            << im->staLat << ","  // AERONET station coordinates
            << im->aerMean550 << ","  // Mean AERONET AOD at 550 nm
            << im->aerStd550 << ","   // Standard deviation of AERONET AOD at 550 nm
            << im->satMean550 << ","  // Mean satellite AOD at 550 nm
            << im->satStd550 << "," // Standard deviation of satellite AOD at 550 nm
            << "\"" << aod550_values.str() << "\","  // AOD at 550 nm for the pixel
            << "\"" << qf_values.str() << "\"," ; // Quality flag for all pixels (enclosed in quotes);
            
   }
   // output.write(reinterpret_cast <const char*>(&nMatch),sizeof(int));
   // for (vector<MatchupRecord>::iterator im=matchData.begin(); im!=matchData.end(); ++im) {
   //    int npix = im->nPixs;
   //    output.write(reinterpret_cast <const char*>(&im->nMeas),sizeof(int));
   //    output.write(reinterpret_cast <const char*>(&im->nPixs),sizeof(int));
   //    output.write(reinterpret_cast <const char*>(&im->satTime),sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->staName),32);
   //    output.write(reinterpret_cast <const char*>(&im->staLon),sizeof(float)*2);
   //    output.write(reinterpret_cast <const char*>(im->meas),(NUM_AERONET_WAVELENGTH+1)*sizeof(float)*im->nMeas);
   //    output.write(reinterpret_cast <const char*>(im->lon),npix*sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->lat),npix*sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->solzen),npix*sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->satzen),npix*sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->relazi),npix*sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->sctang),npix*sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->aod550),npix*sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->ae),npix*sizeof(float));
   //    output.write(reinterpret_cast <const char*>(im->lndSea),npix*sizeof(Int16));
   // }
   output.close();
   
   /****** free up memories ******/
  
   for (map<string,float*>::iterator it=aerInfo.begin(); it!=aerInfo.end(); ++it) {
      delete[] it->second;
   }
  
   for (vector<OneDayAerAod>::iterator ia=aerData.begin(); ia!=aerData.end(); ++ia) {
      delete[] ia->meas;
   }

/*    
   for (map<string,float*>::iterator it=granInfo.begin(); it!=granInfo.end(); ++it) {
      delete[] it->second;
   }
*/
   
   for (vector<MatchupRecord>::iterator im=matchData.begin(); im!=matchData.end(); ++im)
       freeMemForMatch(*im);
         
   return 0;
}
