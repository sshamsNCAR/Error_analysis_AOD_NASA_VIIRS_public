
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cstring>
#include <unistd.h>
#include <cmath>
#include <vector>
#include <stdexcept>
#include "AER_VIIRS_EPS_Match.h"
using namespace std;


/*******************************************************************************
            P R I V A T E   F U N C T I O N S
*******************************************************************************/


/* compute the two-point distance assuming a spherical Earth
   Inputs:
     lon1: longitude of the first point
     lat1: latitude of the first point
     lon2: longitude of the second point
     lat2: latitude of the second point
   Output:
     distance in km   
*/

float distance(float lon1, float lat1, float lon2, float lat2)
{
    double coslat1 = cos(lat1 * DEG2RAD);
    double coslat2 = cos(lat2 * DEG2RAD);
    double sinlat1 = sin(lat1 * DEG2RAD);
    double sinlat2 = sin(lat2 * DEG2RAD);
    double temp = coslat1*coslat2*cos((lon1-lon2)*DEG2RAD) + sinlat1*sinlat2;
                         
    if (abs(temp-1.) < 1.E-10) 
       return 0.0;
    else {  
       double angle = acos(temp);
       return float(EARTH_RADIUS * angle);            
    }            
}

/******************************************************************************/
/* return true if c is between a and b
   idl: true for handling cross International Date Line longitude pair (a,b)
*/    
bool isWithin(float a, float b, float c, bool checkIDL=false)
{
   bool idl = false;
   if (checkIDL)
       if (((a<0.) != (b<0.)) && (abs(a) > 100.))  idl = true;
   
   if (idl) return !((c >= a && c <= b) || (c >= b && c <= a));
   else     return ((c >= a && c <= b) || (c >= b && c <= a));
}

/******************************************************************************/

bool isLeapYear(int year) 
{
  if (((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0))
     return true;
  else
     return false;
}
/******************************************************************************/
/*  Day shift  */

int dayShift(int yyyyddd, int shift)
{
   int year = yyyyddd/1000;
   int day = yyyyddd%1000;
   
   day += shift;
   while (true) {
      int ydays;
      isLeapYear(year) ? ydays=366 : ydays=365;
      if (day > ydays) {
         day -= ydays;
         year++;
      } 
      else if (day < 1) {
         year--;
         isLeapYear(year) ? ydays=366 : ydays=365;
         day += ydays;
      }
      else break;
   }
   
   return year*1000+day;
}
/******************************************************************************/

bool file_exists(const string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

/******************************************************************************/

string GetStdoutFromCommand(string cmd) {
   string data;
   FILE * stream;
   const int max_buffer = 256;
   char buffer[max_buffer];
   //cmd.append(" 2>&1");

   stream = popen(cmd.c_str(), "r");
   if (stream) {
      while (!feof(stream))
         if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
      pclose(stream);
   }
   return data;
}
/******************************************************************************/

// sshams added for interpolation to get AOD at 550 nm
float interpolateAOD550(const float wavelengths[], const float aods[], int numWavelengths) {
   float targetWavelength = 0.550;

   // interpolation of 550 nm AOD using three wavenumbers 440, 500, and 675 nm

   // Extract the three wavelengths and corresponding AOD values
   float x[3] = {wavelengths[0], wavelengths[1], wavelengths[2]};
   float y[3] = {aods[0], aods[1], aods[2]};

   // Check for invalid AOD values (-999)
   for (int i = 0; i < 3; i++) {
      if (y[i] == -999) {
         return -999.0;  // Return -999 if any AOD value is invalid
      }
   }

   // Perform cubic spline interpolation
   // Calculate coefficients for the cubic polynomial
   float a = y[0];
   float b = (y[1] - y[0]) / (x[1] - x[0]);
   float c = ((y[2] - y[1]) / (x[2] - x[1]) - b) / (x[2] - x[0]);
   float d = ((y[2] - y[1]) / (x[2] - x[1]) - b) / ((x[2] - x[0]) * (x[2] - x[1]));

   // Interpolate for 550 nm
   float dx = targetWavelength - x[0];
   float interpolatedAOD = a + b * dx + c * dx * dx + d * dx * dx * dx;

   return interpolatedAOD;

   
}
/******************************************************************************/
/* 
   Get AERONET data within a give time window
     -- used to expand the daily data with previous/next day data
 */
int getAerAod(string fileName, int matchDay, int dayShift, float timeWindow[2], vector<AerAod>& obs)
{
   ifstream staFile(fileName.c_str(), ios::in);
   if (staFile.is_open()) {
      AerAod oneObs;
      string oneLine;
      int start, end;
      float time;
      // int numPassField = 4;
      // bool dataBegin = false;
      while(getline(staFile, oneLine)) {   /* loop over measurements for one station */
         end = oneLine.find(',', 0);
         if (end == string::npos) continue;  /* bypass the header */
         // if (!dataBegin) {  /* bypass the field naming line */
         //    if (oneLine.compare(0,12,"AERONET_Site") == 0) numPassField=5;  // for V1.5
         //    dataBegin = true;
         //    continue;
         // }
         
         /* Bypass the fields:
            V2.0 - Date(dd:mm:yyyy),Time(hh:mm:ss),Day_of_Year,Day_of_Year(Fraction) 
            V1.5 - AERONET_Site,Date(dd:mm:yyyy),Time(hh:mm:ss),Day_of_Year,Day_of_Year(Fraction) */
         // end = -1;
         // for (int i=0; i<numPassField; i++) {
         //    start = end+1;
         //    end = oneLine.find(',', start);
         // } 
         // Extract the fractional time
         start = end + 1;
         end = oneLine.find(',', start);
         
         start = end + 1;
         end = oneLine.find(',', start);
         // Extract the time field (hh:mm:ss)

         time = (float)atof((oneLine.substr(start, end-start)).c_str());
         if ((int)(time) != matchDay) {
            cout << "Skipping line (time does not match matchDay): time = " << time << ", matchDay = " << matchDay << endl;
            continue;
         }


         time = time - (int)(time);
         //if (abs(time) < 1E-5) time = 1.0; 
         if (time >= timeWindow[0] && time <= timeWindow[1]) {
            oneObs.time = time + dayShift;
            for (int i=0; i<NUM_AERONET_WAVELENGTH; i++) {
               start = end+1;
               end = oneLine.find(',', start);
               oneObs.aods[i] = (float)atof((oneLine.substr(start, end-start)).c_str());
            }

            // Interpolate AOD at 0.550 µm
            const float interp_550_wavelengths[3] = {0.440, 0.500, 0.675};
            // Function to extract AOD values for the required wavelengths
            const float interp_selectedAODs[3] = {
               oneObs.aods[INDEX_440],  // AOD at 440 nm (18th element, 0-based index is 17)
               oneObs.aods[INDEX_500],  // AOD at 500 nm (15th element, 0-based index is 14)
               oneObs.aods[INDEX_675]    // AOD at 675 nm (6th element, 0-based index is 5)
            };

            oneObs.aod550 = interpolateAOD550(interp_550_wavelengths, interp_selectedAODs, 3);

            obs.push_back(oneObs);
         }
      }  // end of measurements loop for one station

      staFile.close();
      return PROC_SUCCEED;
   }
   else {
     //cout << "getAerAod(): Cannot open AERONET file "+fileName << endl;
     cout << "getAerAod(): Cannot open AERONET file: " << fileName << endl;

     return PROC_FAIL;
   }
} 

/******************************************************************************/

/* Check whether the granule is at night
   Inputs:
      granFile: granule AERO EDR HDF5 filename 
      variable: variable path as "/group1/group2/varname"
   Return:
      true if night granule, false otherwise    
*/

bool isNightGranule(string granFile, string variable)
{   
   hid_t   file, dataset, grp1, grp2, attr, type;
   herr_t  status;
   char    dayNightFlag[32];
   size_t  pos1, pos2;

   file = H5Fopen(granFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   
   pos1 = variable.find_first_of("/");
   if (pos1 != 0) {
      pos1 = -1;
      pos2 = pos1;
   }
   else {
      pos2 = variable.find_first_of("/", pos1+1);
   }
   grp1 = H5Gopen1(file, (variable.substr(pos1+1, pos2-pos1-1)).c_str());   
   
   pos1 = pos2;
   pos2 = variable.find_first_of("/", pos1+1);
   grp2 = H5Gopen1(grp1, (variable.substr(pos1+1, pos2-pos1-1)).c_str());
   
   pos1 = pos2;
   pos2 = variable.find_first_of("/", pos1+1);
   dataset = H5Dopen (grp2, (variable.substr(pos1+1, pos2-pos1-1)).c_str(), H5P_DEFAULT);
   
   attr = H5Aopen_name(dataset, "N_Day_Night_Flag");
   type = H5Aget_type(attr);
   status = H5Aread(attr, type, dayNightFlag);   
   status = H5Tclose(type);
   status = H5Aclose(attr);
   status = H5Dclose(dataset);
   status = H5Gclose (grp2);
   status = H5Gclose (grp1);
   status = H5Fclose (file);
   //printf("The value of the attribute is: %s \n", dayNightFlag); 
   
   if ((strcmp(dayNightFlag, "Night") == 0) ||
       (strcmp(dayNightFlag, "Both") == 0)) 
      return true;
   else
      return false;   
}
/******************************************************************************/
/* read hdf5 data */

int getH5DataDim2(string fileName, string dataSetName, int& num_row, int& num_col)
{
   string   message;
   hid_t    file, dataset, dataspace; 
   hsize_t  dims[2];
   int      status;
   
   /* open file */
   file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) {
      cout << "getH5DataDim2(): Cannot open HDF5 file: "+fileName << endl;
      return PROC_FAIL;
   }
  
   /* open dataset */
   dataset = H5Dopen2(file, dataSetName.c_str(), H5P_DEFAULT);
   if (dataset < 0) {
      cout << "getH5DataDim2(): Cannot open HDF5 dateset: "+fileName+" "+dataSetName << endl;
      return PROC_FAIL;
   }
   
   /* open data space */
   dataspace = H5Dget_space(dataset);
   if (dataspace < 0) {
      cout << "getH5DataDim2(): Cannot open HDF5 dataspace: "+fileName+" "+dataSetName << endl;
      return PROC_FAIL;
   }
   
   /* get dimension */
   status  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
   if (status < 0) {
      cout << "getH5DataDim2(): Cannot get HDF5 data dimension: "+fileName+" "+dataSetName << endl;
      return PROC_FAIL;
   }
   
   num_row = dims[0];
   num_col = dims[1];
   H5Sclose(dataspace);
   H5Dclose(dataset);
   H5Fclose(file);
   return PROC_SUCCEED;
}

int readH5Data(string fileName, string dataSetName, 
               hid_t type_id, void* data, int dimension=1, bool hyperslab=false, 
               const int *start_in=NULL, const int *stride_in=NULL, 
               const int *count_in=NULL, const int *block_in=NULL)
{
   string   message;
   hid_t    file, dataset, dataspace, memspace; 
   hsize_t  dims[MAX_DIM_SIZE]; /* assume maximum dimension size as 4 */
   int      status;
   
   /* open file */
   file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) {
      cout << "readH5Data(): Cannot open HDF5 file: "+fileName << endl;
      return PROC_FAIL;
   }
  
   /* open dataset */
   dataset = H5Dopen2(file, dataSetName.c_str(), H5P_DEFAULT);
   if (dataset < 0) {
      cout << "readH5Data(): Cannot open HDF5 dateset: "+fileName+" "+dataSetName << endl;
      return PROC_FAIL;
   }
   
   /* open data space */
   dataspace = H5Dget_space(dataset);
   if (dataspace < 0) {
      cout << "readH5Data(): Cannot open HDF5 dataspace: "+fileName+" "+dataSetName << endl;
      return PROC_FAIL;
   }
   
   /* get dimension */
   if (dimension > 1) {
      status  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
      if (status < 0) {
         cout << "readH5Data(): Cannot get HDF5 data dimension: "+fileName+" "+dataSetName << endl;
         return PROC_FAIL;
      }
   }
   
   /* read data */
   if (hyperslab) {
      hsize_t start[2] = {start_in[0], start_in[1]};
      hsize_t stride[2] = {stride_in[0], stride_in[1]};
      hsize_t count[2] = {count_in[0], count_in[1]};
      hsize_t block[2] = {block_in[0], block_in[1]};
      /* read a hyperslab region within the data */
      memspace = H5Screate_simple (dimension, count, NULL); 
      status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, start,
                                    stride, count, block);
      if (status < 0) {
         cout << "readH5Data(): Cannot select hyperslab: "+fileName+" "+dataSetName << endl;
         return PROC_FAIL;
      }
      status = H5Dread (dataset, type_id, memspace, dataspace, H5P_DEFAULT, data);
   }
   else {
      /* read all data */
      status = H5Dread(dataset, type_id, H5S_ALL, dataspace, H5P_DEFAULT, data);
   }
   
   if (status < 0) {
      cout << "readH5Data(): Cannot read HDF5 dataset: "+fileName+" "+dataSetName << endl;
      return PROC_FAIL;
   }
   
   if (hyperslab) H5Sclose (memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);
   H5Fclose(file);
   return PROC_SUCCEED;
}

/******************************************************************************/
/* read dark-target VIIRS AOD data */

// sshams added
int readDBAod(string aodFile, MatchupRecord& mr, int offset, const int *start,
              const int *stride, const int *count, const int *block)
{
   int status;
   string dataSetName;
   int npix = count[0] * count[1];  // Total number of pixels

   // Allocate a single temporary buffer
   float *tmpData = new float[npix];

   // Read land AOD
   dataSetName = "Aerosol_Optical_Thickness_550_Land_Best_Estimate";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }

   // Process land AOD
   for (int i = 0; i < npix; i++) {
      if (tmpData[i] > -0.001) {
         *(mr.aod550 + offset + i) = tmpData[i];
         *(mr.lndSea + offset + i) = 1;  // Land
      }
   }

   // Read ocean AOD
   dataSetName = "Aerosol_Optical_Thickness_550_Ocean_Best_Estimate";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }

   // Process ocean AOD
   for (int i = 0; i < npix; i++) {
      if (*(mr.aod550 + offset + i) <= -0.001) {  // Only update if land AOD is invalid
         *(mr.aod550 + offset + i) = tmpData[i];
         *(mr.lndSea + offset + i) = 0;  // Ocean
      }
   }

   // Read land Angstrom Exponent
   dataSetName = "Angstrom_Exponent_Land_Best_Estimate";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }

   // Process land Angstrom Exponent
   for (int i = 0; i < npix; i++) {
      if (*(mr.lndSea + offset + i) == 1) {  // Only update for land pixels
         *(mr.ae + offset + i) = tmpData[i];
      }
   }

   // Read ocean Angstrom Exponent
   dataSetName = "Angstrom_Exponent_Ocean_Best_Estimate";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }

   // Process ocean Angstrom Exponent
   for (int i = 0; i < npix; i++) {
      if (*(mr.lndSea + offset + i) == 0) {  // Only update for ocean pixels
         *(mr.ae + offset + i) = tmpData[i];
      }
   }

   // Read Solar Zenith Angle
   dataSetName = "Solar_Zenith_Angle";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }
   memcpy(mr.solzen + offset, tmpData, npix * sizeof(float));

   // Read Viewing Zenith Angle
   dataSetName = "Viewing_Zenith_Angle";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }
   memcpy(mr.satzen + offset, tmpData, npix * sizeof(float));

   // Read Relative Azimuth Angle
   dataSetName = "Relative_Azimuth_Angle";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }
   memcpy(mr.relazi + offset, tmpData, npix * sizeof(float));

   // Read Scattering Angle
   dataSetName = "Scattering_Angle";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }
   memcpy(mr.sctang + offset, tmpData, npix * sizeof(float));

   // Read Quality Flag
   dataSetName = "Aerosol_Optical_Thickness_QA_Flag_Land";
   status = readH5Data(aodFile, dataSetName, H5T_NATIVE_FLOAT, tmpData, 2, true, start, stride, count, block);
   if (status < 0) {
      delete[] tmpData;
      return PROC_FAIL;
   }
   memcpy(mr.qf + offset, tmpData, npix * sizeof(Int16));

   // Clean up temporary buffer
   delete[] tmpData;

   return PROC_SUCCEED;
}

// int readDBAod(string aodFile, MatchupRecord& mr, int offset, const int *start,
//               const int *stride, const int *count, const int *block)
// {
//    int    status;
//    string dataSetName;
//    int npix = count[0]*count[1];   
   
//    float *tmpLnd = new float[npix];
//    float *tmpOcn = new float[npix];
//    float *tmpLndAe = new float[npix];
//    float *tmpOcnAe = new float[npix];
   
//    dataSetName = "Aerosol_Optical_Thickness_550_Land_Best_Estimate";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpLnd, 2,
//                        true, start, stride, count, block);
//    dataSetName = "Aerosol_Optical_Thickness_550_Ocean_Best_Estimate";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpOcn, 2,
//                        true, start, stride, count, block);
//    dataSetName = "Angstrom_Exponent_Land_Best_Estimate";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpLndAe, 2,
//                        true, start, stride, count, block);
//    dataSetName = "Angstrom_Exponent_Ocean_Best_Estimate";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpOcnAe, 2,
//                        true, start, stride, count, block);
//    for (int i=0; i<npix; i++) {
//       if (tmpLnd[i] > -0.001) {
//          *(mr.aod550+offset+i) = tmpLnd[i];
//          *(mr.ae+offset+i) = tmpLndAe[i];
//          *(mr.lndSea+offset+i) = 1;
//       }
//       else {
//          *(mr.aod550+offset+i) = tmpOcn[i];
//          *(mr.ae+offset+i) = tmpOcnAe[i];
//          *(mr.lndSea+offset+i) = 0;
//       }
//    }
      
      
//    dataSetName = "Solar_Zenith_Angle";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpLnd, 2,
//                        true, start, stride, count, block);
//    memcpy(mr.solzen+offset, tmpLnd, npix*sizeof(float));
      
//    dataSetName = "Viewing_Zenith_Angle";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpLnd, 2,
//                        true, start, stride, count, block);
//    memcpy(mr.satzen+offset, tmpLnd, npix*sizeof(float));
      
//    dataSetName = "Relative_Azimuth_Angle";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpLnd, 2,
//                        true, start, stride, count, block);
//    memcpy(mr.relazi+offset, tmpLnd, npix*sizeof(float));
      
//    dataSetName = "Scattering_Angle";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpLnd, 2,
//                        true, start, stride, count, block);
//    memcpy(mr.sctang+offset, tmpLnd, npix*sizeof(float));

//    dataSetName ="Aerosol_Optical_Thickness_QA_Flag_Land";
//    status = readH5Data(aodFile,  dataSetName, H5T_NATIVE_FLOAT, tmpLnd, 2,
//       true, start, stride, count, block);
//    memcpy(mr.qf+offset, tmpLnd, npix*sizeof(Int16));
   
//    delete[] tmpLnd;
//    delete[] tmpOcn;
//    delete[] tmpLndAe;
//    delete[] tmpOcnAe;
//    return PROC_SUCCEED;
// }

/******************************************************************************/
/* 
    Get the granule boundaries - west/east longitude; south/north latitude
    Input:
       geoInMem:  true if lon/lat data in memory
       gmtcoFile: name of the VIIRS granule GMTCO file
       granLon:   pointer to the granule longitude (if geoInMem is true)
       granLon:   pointer to the granule latitude (if geoInMem is true)
    Outputs:
       cornerLon:  longtitude of four corners (array index: [0,0],[0,n],[m,n],[m,0]) 
       cornerLat:  latitude of four corners
    
    return PROC_FAIL if any granule corners of is beyond MAXLAT   
*/
int getGranGeo(bool& geoInMem, string aodFile, float* granLon, float* granLat, 
               float cornerLon[4], float cornerLat[4], int& num_row, int& num_col)
{ 
   int status;
   
   status = getH5DataDim2(aodFile, "Longitude", num_row, num_col);
   if (status == PROC_FAIL)  {
      cout << "getGranGeo(): Cannot get 2D dimension from "<< aodFile << endl;
      return PROC_FAIL;
   }
   
   /* read lon/lat if needed */
   if (!geoInMem) {
      status = readH5Data(aodFile, "Longitude", H5T_NATIVE_FLOAT, granLon, 2);
      if (status == PROC_FAIL)  {
         cout << "getGranGeo(): Cannot read longitude from "<< aodFile << endl;
         return PROC_FAIL;
      }
      status = readH5Data(aodFile, "Latitude", H5T_NATIVE_FLOAT, granLat, 2);
      if (status == PROC_FAIL) {
         cout << "getGranGeo(): Cannot read latitude from "<< aodFile << endl;
         return PROC_FAIL;
      }
      geoInMem = true;
   }
   
   /* granule corner longitude/latitude ([0,0],[0,n],[m,n],[m,0]) 
      return fail if any one corner is beyond the MAXLAT */
   for (int i=0; i<num_row; i++) {
      int j = i*num_col;
      cornerLon[0] = *(granLon+j);
      cornerLat[0] = *(granLat+j);
      cornerLon[1] = *(granLon+j+num_col-1);
      cornerLat[1] = *(granLat+j+num_col-1);
      if (cornerLon[0] > -181. && cornerLon[1] > -181. &&
          cornerLat[0] >  -91. && cornerLat[1] > -91.) break;
   }
   if (abs(cornerLat[0]) > MAXLAT || abs(cornerLat[1]) > MAXLAT) return PROC_FAIL;
   for (int i=num_row-1; i>0; i--) {
      int j = i*num_col;
      cornerLon[2] = *(granLon+j+num_col-1);
      cornerLat[2] = *(granLat+j+num_col-1);
      cornerLon[3] = *(granLon+j);
      cornerLat[3] = *(granLat+j);
      if (cornerLon[2] > -181. && cornerLon[3] > -181. &&
          cornerLat[2] >  -91. && cornerLat[3] > -91.) break;
   }
   if (abs(cornerLat[2]) > MAXLAT || abs(cornerLat[3]) > MAXLAT) return PROC_FAIL;   
   
   return PROC_SUCCEED;
}
   

/******************************************************************************/
/* 
    Get the granule boundaries - west/east longitude; south/north latitude
    Input:
       cornerLon:  longtitude of four corners (array index: [0,0],[0,n],[m,n],[m,0]) 
       cornerLat:  latitude of four corners
    Outputs:
       westGranLon,  eastGranLon, northGranLat, southGranLat: granule boundaries
       asndNode: true if granule is in ascending node 
       crsIDL:   true if granule cross the International Date Line (+/- 180)
    return PROC_FAIL if    
*/
void getGranBoundary(float cornerLon[4], float cornerLat[4], float& westGranLon, float& eastGranLon, 
                     float& northGranLat, float& southGranLat, bool& asndNode, bool& crsIDL)
{ 
   /* get the side corner lon/lat */   
   float eastBndLon[2];
   float westBndLon[2];
   float northBndLat[2];
   float southBndLat[2];
   if (cornerLat[1] < cornerLat[2]) {
      /* ascending node: four corners are in SE->SW->NW->NE sequence */
      asndNode = true;
      eastBndLon[0] = cornerLon[0];
      eastBndLon[1] = cornerLon[3];
      westBndLon[0] = cornerLon[1];
      westBndLon[1] = cornerLon[2];
      northBndLat[0] = cornerLat[2];
      northBndLat[1] = cornerLat[3];
      southBndLat[0] = cornerLat[0];
      southBndLat[1] = cornerLat[1];
   }
   else {
      /* descending node: four corners are in NW->NE->SE->SW sequence */
      asndNode = false;
      eastBndLon[0] = cornerLon[1];
      eastBndLon[1] = cornerLon[2];
      westBndLon[0] = cornerLon[0];
      westBndLon[1] = cornerLon[3];
      northBndLat[0] = cornerLat[0];
      northBndLat[1] = cornerLat[1];
      southBndLat[0] = cornerLat[2];
      southBndLat[1] = cornerLat[3];
   }   
   
   /* get the latitude of north/south boundaries */
   northGranLat = (northBndLat[0] > northBndLat[1] ? northBndLat[0] : northBndLat[1]);
   southGranLat = (southBndLat[0] < southBndLat[1] ? southBndLat[0] : southBndLat[1]);
   
   
   /* get the longitude of east boundary  */
   if ((eastBndLon[0]<0.) == (eastBndLon[1]<0.)) { /* same sign */
      eastGranLon = (eastBndLon[0] > eastBndLon[1] ? eastBndLon[0] : eastBndLon[1]);
   }
   else { /* one positive and one negative */
      if (abs(eastBndLon[0]) > 100.) /* crsIDL=true */
         eastGranLon = (eastBndLon[0] < eastBndLon[1] ? eastBndLon[0] : eastBndLon[1]);
      else 
         eastGranLon = (eastBndLon[0] > eastBndLon[1] ? eastBndLon[0] : eastBndLon[1]);
   } 
   
   /* get the longitude of west boundary  */
   if ((westBndLon[0]<0.) == (westBndLon[1]<0.)) { /* same sign */
      westGranLon = (westBndLon[0] < westBndLon[1] ? westBndLon[0] : westBndLon[1]);
   }
   else { /* one positive and one negative */
      if (abs(westBndLon[0]) > 100.) /* crsIDL=true */
         westGranLon = (westBndLon[0] > westBndLon[1] ? westBndLon[0] : westBndLon[1]);
      else 
         westGranLon = (westBndLon[0] < westBndLon[1] ? westBndLon[0] : westBndLon[1]);
   }
   
   /* check crossing IDL */
   if ((eastGranLon<0.) != (westGranLon<0.)) {
      if (abs(eastGranLon) > 100.) crsIDL = true;
      else crsIDL = false;
   }
   else crsIDL = false;
} 

/******************************************************************************/

/* check whether station match with granule both in time and space 
   inputs: 
     sta: structure of {OneDayAerAod}
     granTime: granule middle time (start+end/2) in fractional hour (0-24)
     crsIDL: true if granule cross International Date Line
     eastGranLon, westGranLon, northGranLat, southGranLat: granule boundaries
   output:
     westStaLon, eastStaLon, northStaLat, southStaLat: station matching boundaries
     true if station matched the granule 
*/
bool staMatchGranule(OneDayAerAod& sta, float granTime, bool crsIDL, 
                     float westGranLon, float eastGranLon, float northGranLat, float southGranLat,
                     float& westStaLon, float& eastStaLon, float& northStaLat, float& southStaLat)
{
   /* get the matching square domain over station */
   float disLon = RAD2DEG*MATCH_RADIUS/(EARTH_RADIUS*cos(sta.staLat*DEG2RAD))+0.01;
   float disLat = RAD2DEG*MATCH_RADIUS/EARTH_RADIUS+0.01;
   westStaLon = sta.staLon - disLon;
   if (westStaLon < -180. )  westStaLon += 360.;
   eastStaLon = sta.staLon + disLon;
   if (eastStaLon > 180. )   eastStaLon -= 360.;
   northStaLat = sta.staLat + disLat;
   southStaLat = sta.staLat - disLat;
   
   /* spatical match */
   if (northStaLat < southGranLat || southStaLat > northGranLat)
      return false;
   if (crsIDL) {
      if ((westStaLon > eastGranLon && westStaLon < westGranLon) &&
          (eastStaLon > eastGranLon && eastStaLon < westGranLon))
         return false;
   }
   else{
      if (westStaLon > eastGranLon || eastStaLon < westGranLon)
         return false;
   }
   
   /* temporal match */
   float timeWindow[2] = {(granTime-MATCH_TIME_WINDOW*0.5)/24., (granTime+MATCH_TIME_WINDOW*0.5)/24.};
   int ntm = 0;
   for (int i=0; i<sta.nmeas-1; i++)
      if (sta.meas[i].time >= timeWindow[0] && sta.meas[i].time <= timeWindow[1]) ntm++;
   if (ntm < NVLDAER)
      return false;
      
   return true;
}                     

/******************************************************************************/
/* given granule lon/lat and matching domain, find the row/column indices 
   such the rectangular granule row/column domain enclose the matching domain 
    
   Inputs:
      lon:  granule longitude
      lat:  granule latitude
      asndNode: true if granule is in ascending node
      crsIDL: true if granule cross International Date Line
      westLon, eastLon, northLat, southLat: rectangular matching domain  
   Outputs:
      row[2]: row indices of the granule rectangular matching domain
      col[2]: column indices
   Return:
      -1: no pixels found within the matching lon/lat domain
      0:  pixels found, no need for next granule
      1:  pixels found, need next granule for the complete domain
*/ 

int findMatchPixels(float *lon, float *lat, int num_row, int num_col, bool asndNode, bool crsIDL, 
                    float westLon, float eastLon, float northLat, float southLat, 
                    int row[2], int col[2])
{                    
   int i, firstRow, lastRow;
   for (i=0; i<num_row; i++) {
      if (*(lon+i*num_col) > -200) {
         firstRow = i; 
         break;
      }
   }
   for (i=num_row-1; i>=0; i--) {
      if (*(lon+i*num_col) > -200) {
         lastRow = i; 
         break;
      }
   }
   
   int south, north, east, west, incS2N, incE2W;
   if (asndNode) {
      south = firstRow;
      north = lastRow;
      incS2N = 1;
      east = 0;
      west = num_col-1;
      incE2W = 1;
   }
   else {
      south = lastRow;
      north = firstRow;
      incS2N = -1;
      east = num_col-1;
      west = 0;
      incE2W = -1;
   }
   
   
   bool search = true;
   bool found;
   while (search) {
      search = false;
      
      // search for south bound
      i = south;
      while (i != north) {
         if (*(lat+i*num_col+east) >= southLat || *(lat+i*num_col+west) >= southLat) {
            if (i != south) search = true;
            south = i;
            break;
         }
         i += incS2N;
      }
      if (i == north) return -1;
   
      // search for north bound
      i = north;
      while (i != south) {
         if (*(lat+i*num_col+east) <= northLat || *(lat+i*num_col+west) <= northLat) {
            if (i != north) search = true;
            north = i;
            break;
         }
         i -= incS2N;
      }
      if (i == south) return -1;
      
      // search for east bound
      found = false;
      i = east;
      while (i != west) {
         if (crsIDL) {
            if (eastLon > 0.)
               found = ((*(lon+north*num_col+i) > 0. &&  *(lon+north*num_col+i) <= eastLon) ||
                        (*(lon+south*num_col+i) > 0. &&  *(lon+south*num_col+i) <= eastLon));
            else  
               found = ((*(lon+north*num_col+i) > 0. ||  *(lon+north*num_col+i) <= eastLon) ||
                        (*(lon+south*num_col+i) > 0. ||  *(lon+south*num_col+i) <= eastLon));           
         }
         else
            found = (*(lon+north*num_col+i) <= eastLon || *(lon+south*num_col+i) <= eastLon);
         
         if (found) {
            if (i != east) search = true;
            east = i;
            break;
         }
         i += incE2W;
      }
      if (i == west) return -1;
      
      // search for west bound
      found = false;
      i = west;
      while (i != east) {
         if (crsIDL) {
            if (westLon < 0.)
               found = ((*(lon+north*num_col+i) < 0. &&  *(lon+north*num_col+i) >= westLon) ||
                        (*(lon+south*num_col+i) < 0. &&  *(lon+south*num_col+i) >= westLon));
            else  
               found = ((*(lon+north*num_col+i) < 0. ||  *(lon+north*num_col+i) >= westLon) ||
                        (*(lon+south*num_col+i) < 0. ||  *(lon+south*num_col+i) >= westLon));           
         }
         else
            found = (*(lon+north*num_col+i) >= westLon || *(lon+south*num_col+i) >= westLon);
            
         if (found) {
            if (i != west) search = true;
            west = i;
            break;
         }
         i -= incE2W;
      }
      if (i == east) return -1;
   }
   
   if (asndNode) {
      row[0] = south;
      row[1] = north;
      col[0] = east;
      col[1] = west;
      if (north == lastRow) return 1;
      else return 0;
   }
   else {
      row[0] = north;
      row[1] = south;
      col[0] = west;
      col[1] = east;
      if (south == lastRow) return 1;
      else return 0;
   }
   
}

/******************************************************************************/

void allocMemForMatch (MatchupRecord &mr, int np) 
{
   mr.nPixs = np;
   
   mr.solzen = new float[np];
   mr.satzen = new float[np];
   mr.relazi = new float[np];
   mr.sctang = new float[np];
   mr.lon = new float[np];  
   mr.lat = new float[np];
   mr.aod550 = new float[np];
   mr.ae = new float[np];
   mr.lndSea = new short int[np];
   mr.qf = new short int[np];
   for (int i=0; i<np; i++) {
      *(mr.aod550+i) = MISVAL_FLOAT;
      *(mr.ae+i) = MISVAL_FLOAT;
      
   }
}
/******************************************************************************/
/* only collect the pixels with AOD retrievals */
    

int collectAod (MatchupRecord &mr)
{
   int np = 0;
   float dist;
   for (int i=0; i<mr.nPixs; i++) {
      if (*(mr.aod550+i) >-1.) {

         dist = distance(mr.staLon, mr.staLat, *(mr.lon+i), *(mr.lat+i));
         if (dist > MATCH_RADIUS) 
            *(mr.aod550+i) = -9.999;
         else
            np++;
      }   
   }
   if (np >= MIN_N_MATCH_PIX) { 
      float *lon = new float[np];  
      float *lat = new float[np];
      float *aod550 = new float[np];
      float *ae = new float[np];
      Int16*lndSea = new Int16[np];
      float *solzen = new float[np];
      float *relazi = new float[np];
      float *satzen = new float[np];
      float *sctang = new float[np];
      int n = 0;
      for (int i=0; i<mr.nPixs; i++) {
         if (*(mr.aod550+i) >-1.) {
            *(lon+n) = *(mr.lon+i);
            *(lat+n) = *(mr.lat+i);
            *(aod550+n) = *(mr.aod550+i);
            *(ae+n) = *(mr.ae+i);
            *(lndSea+n) = *(mr.lndSea+i);
            *(solzen+n) = *(mr.solzen+i);
            *(relazi+n) = *(mr.relazi+i);
            *(satzen+n) = *(mr.satzen+i);
            *(sctang+n) = *(mr.sctang+i);
      
            n++;
         }   
      } 
      delete[] mr.lon;    
      delete[] mr.lat;
      delete[] mr.aod550;
      delete[] mr.ae;
      delete[] mr.lndSea;
      delete[] mr.solzen;
      delete[] mr.relazi;
      delete[] mr.satzen;
      delete[] mr.sctang;
      
      mr.nPixs = np;
      mr.lon = lon;   
      mr.lat = lat;
      mr.aod550 = aod550;
      mr.ae = ae;
      mr.lndSea = lndSea;
      mr.solzen = solzen;
      mr.relazi = relazi;
      mr.satzen = satzen;
      mr.sctang = sctang;
      return PROC_SUCCEED;
   } 
   else { 
      delete[] mr.lon;    
      delete[] mr.lat;
      delete[] mr.aod550;
      delete[] mr.ae;
      delete[] mr.lndSea;
      delete[] mr.solzen;
      delete[] mr.relazi;
      delete[] mr.satzen;
      delete[] mr.sctang;
      mr.nPixs = 0;
      return PROC_FAIL;
   }
} 


/******************************************************************************/
/* check whether staion has been matched with previous and extended to current granules
    
   Inputs:
      matchData:  collected match records
      staName: AERONET station name
      currGranTime: current granule station time [0,24]
   Return:
      true if this station has been matched with pervious and current granules
*/ 
bool prevMatched(vector<MatchupRecord>& matchData, string staName, float currGranTime)
{
   for (vector<MatchupRecord>::iterator im=matchData.begin(); im!=matchData.end(); ++im) {
      if (staName.compare(im->staName) == 0) {
         float timeDiff = currGranTime - im->satTime;
         if (timeDiff > 0 && timeDiff < 8./60.) return true; 
      }
   }
   return false; 

}
/*******************************************************************************
            P U B L I C   F U N C T I O N S
*******************************************************************************/
/* convert integer to string */
string int_to_str(int num)
{
    stringstream ss;
    ss << num;
    return ss.str();
}
/******************************************************************************/
/* given YYYYDDD return YYYYMMDD */

int calendar_date (int yyyyddd)
{
   int MONDAY[] = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
   int year = yyyyddd/1000;
   int dayofyear = yyyyddd%1000;
   if (isLeapYear(year)) 
      for (int i=1; i<12; i++) MONDAY[i]++; 
      
   int month, day;
   for (int i=0; i<12; i++) {
      if (dayofyear <= MONDAY[i]) {
         month = i+1;
         if (i == 0) day=dayofyear;
         else day=dayofyear-MONDAY[i-1];
         break;
      }
   }
   return year*10000+month*100+day;
}
/******************************************************************************/
/* given YYYYMMDD return YYYYDDD */

int dayofyear (int yyyymmdd)
{
   int MONDAY[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
   int year = yyyymmdd/10000;
   int month =  (yyyymmdd%10000)/100;
   int day = yyyymmdd%100;
   if (isLeapYear(year)) 
      for (int i=2; i<month; i++) MONDAY[i]++; 
    
   return year*1000+ MONDAY[month-1]+day;
}

/******************************************************************************/
/* get AERONET station info (name, longitude, latitude) into a map<staName,[lon,lat]> */

int getAerInfo(map<string,float*>& aerInfo)
{
   string oneLine;
   string aerStaListFile = AERONET_PATH+"aeronet_locations.txt";

   ifstream aerInfoFile(aerStaListFile.c_str(), ios::in);
   if (aerInfoFile.is_open()) {
      // for (int i=0; i<2; i++) getline(aerInfoFile, oneLine);
      // Skip the first line (header)
      getline(aerInfoFile, oneLine);
      int start, end;
      string staName;
      float  *lonLat;
      while(getline(aerInfoFile, oneLine)) {
         start = 0;
         lonLat = new float[2];
         for (int i=0; i<3; i++) {
            end = oneLine.find(',', start);
            switch (i) {
               case 0:  /* station name */
                  staName = oneLine.substr(start, end-start);
                  break;
               case 1:  /* longitude */
                  lonLat[0] = (float)atof((oneLine.substr(start, end-start)).c_str());
                  break;
               case 2:  /* latitude */
                  lonLat[1] = (float)atof((oneLine.substr(start, end-start)).c_str());
                  break;
            }
            start = end+1;
            aerInfo.insert(pair<string,float*>(staName, lonLat));
         }
      }
   }
   else {
      cout << "getAerInfo(): Could not open AERONET station list file" << endl;
      return PROC_FAIL;
   }
   aerInfoFile.close();
   return PROC_SUCCEED;
}   
 
/******************************************************************************/
 /* get one day (plust half of the matching time window from previous/next day) 
    AERONET V3 L1.5 Direct Sun AOD measurements */ 

int getDailyAerAod(map<string,float*>& aerInfo, string yyyyddd, vector<OneDayAerAod>& aerData)
{
   int status;
   int day = atoi(yyyyddd.c_str());
   int preDay = dayShift(day, -1);
   int nextDay = dayShift(day, 1);

   // update by sshams to follow our filename and folder format
   // Extract year and day-of-year for current, previous, and next days
   string currYear = yyyyddd.substr(0, 4);
   string currDay = yyyyddd.substr(4, 3);

   string preYear = to_string(preDay / 1000);
   string preDayOfYear = to_string(preDay % 1000);
   if (preDayOfYear.length() < 3) preDayOfYear = string(3 - preDayOfYear.length(), '0') + preDayOfYear;

   string nextYear = to_string(nextDay / 1000);
   string nextDayOfYear = to_string(nextDay % 1000);
   if (nextDayOfYear.length() < 3) nextDayOfYear = string(3 - nextDayOfYear.length(), '0') + nextDayOfYear;

   // ostringstream str1;
   // str1 << preDay;
   // string preDayStr = str1.str();
   // str1.str("");
   // str1 << nextDay;
   // string nextDayStr = str1.str();
   
   float halfWindow = MATCH_TIME_WINDOW*0.5/24.;
   float timeWindow[3][2] = {{1.0-halfWindow, 1.0}, {0.0, 1.0}, {0.0, halfWindow}};
   int   matchDate[3] = {preDay%1000, day%1000, nextDay%1000};  
                           ;
   for (map<string,float*>::iterator it=aerInfo.begin(); it!=aerInfo.end(); ++it) {
      string staName = it->first;

      // sshams updated
      string aerFileName[3] = {
         AERONET_PATH + preYear + "/" + preDayOfYear + "/" + staName + ".dat",
         AERONET_PATH + currYear + "/" + currDay + "/" + staName + ".dat",
         AERONET_PATH + nextYear + "/" + nextDayOfYear + "/" + staName + ".dat"};

      // Debug print: Print the file paths being checked
      cout << "Looking for files for station: " << staName << endl;
      for (int i = 0; i < 3; i++) {
         cout << "  Checking file: " << aerFileName[i] << endl;
      }
      // string aerFileName[3] = { AERONET_PATH+preDayStr.substr(0,4)+"/"+preDayStr+"/"+staName+".dat", 
      //                           AERONET_PATH+yyyyddd.substr(0,4)+"/"+yyyyddd+"/"+staName+".dat", 
      //                           AERONET_PATH+nextDayStr.substr(0,4)+"/"+nextDayStr+"/"+staName+".dat"}; 
      //  sshams updated the code to work if at least one day it is available and can find the data in there
      /* get AERONET AODs */
      vector<AerAod> obs;
      bool fileAvailable = false;
      for (int i=0; i<3; i++) {
         if(file_exists(aerFileName[i])) {
            fileAvailable = true;
            status = getAerAod(aerFileName[i], matchDate[i], (i-1), timeWindow[i], obs);
            // break;
         }
         else {
            // Debug print: File not found
            cout << "  File not found: " << aerFileName[i] << endl;
         }
      }
      if (!fileAvailable) continue;
      
      
      // for (int i=0; i<3; i++) {
      //    status = getAerAod(aerFileName[i], matchDate[i], (i-1), timeWindow[i], obs);
      // }
      int nobs = obs.size();   /* number of measurements collected for current station */
      if (nobs == 0) {
         cout << "No valid observations for station: " << staName << endl;
          continue;}
      
      /* add current station measurements to the vector */
      OneDayAerAod oneSta;
      oneSta.staName = staName;
      oneSta.staLon = it->second[0];
      oneSta.staLat = it->second[1];
      oneSta.nmeas = nobs;
      oneSta.meas = new AerAod[nobs];
      memcpy(oneSta.meas, obs.data(), nobs*sizeof(AerAod));  
      aerData.push_back(oneSta);
      obs.clear();
       
   }  // end of station loop
   
   if (aerData.size() == 0) return PROC_FAIL;
   else return PROC_SUCCEED;
}

/******************************************************************************/
 /* get VIIRS daytime granule info (time, EPS AOD filename) */ 

int getVIIRSgranInfo(string yyyyddd, map<float,string>&granInfo)
{  
   /* get AOD file names:  YYYY/DDD/AERDT_L2_VIIRS_SNPP.A2012333.2306.011.2020200031009.nc*/
   string year = yyyyddd.substr(0,4);
   string dayOfYear = yyyyddd.substr(4,3);
   string command, pout;
   // command = "ls "+AOD_PATH+year+"/"+dayOfYear+"/AER*.A"+yyyyddd+".*.nc";
   command = "ls "+AOD_PATH+year+"/AERDB*.A"+yyyyddd+".*.nc";
   pout = GetStdoutFromCommand(command);
   
   /* check each granules */
   string aodFile, basename;
   int status, start, end, found, startTime, endTime;
   float granTime;
   start = 0;
   while ((end = pout.find('\n',start)) != string::npos) {
      /* get EPS aod file name */
      aodFile = pout.substr(start,end-start);
      start = end+1;

      /* get the file name without full path */
      found = aodFile.find_last_of('/');
      basename = aodFile.substr(found+1);

 
      /* get satellite overpass time = startingTime? */
      found = basename.find(yyyyddd);
      startTime = atoi((basename.substr(found+8,4)).c_str());
      granTime = startTime/100 + (startTime%100)/60.;
      
      /* insert the valid granule info into Map */
      granInfo.insert(pair<float,string>(granTime, aodFile));
   }
   
   if (granInfo.size() == 0) return PROC_FAIL;
   else  return PROC_SUCCEED;
}

  
/******************************************************************************/
  
int matchup(string yyyymmdd, vector<OneDayAerAod>& aerData, map<float,string>& granInfo, 
            vector<MatchupRecord>& matchData)
{ 
   map<float,string>::iterator nt; 
   bool  asndNode[2];      /* true if granule is in ascending node */
   bool  crsIDL[2];        /* true if granule cross the International Date Line */
   bool  granNextAdj;      /* true if next granule is adjacent (time difference < 2 min) */
   bool  granGeoInMem[2];  /* true if granule lon/lat is in memeory for current and next granule */
   int   curr, next, status;  
   float westGranLon, eastGranLon, northGranLat, southGranLat;   /* granule boundaries */
   float westStaLon, eastStaLon, northStaLat, southStaLat;       /* station matching boundaries */
   
   float *granLon[2];       /* granule lon for two consecutive granules */
   float *granLat[2];       /* granule lat for two consecutive granules */
   for (int i=0; i<2; i++) {
      granLon[i] = new float[250000];  //500x500
      granLat[i] = new float[250000]; 
   }
   int num_row[2];
   int num_col[2];
   
   /* loop over granules */
   curr = 0;
   granGeoInMem[curr] = false;
   for (map<float,string>::iterator it=granInfo.begin(); it!=granInfo.end(); ++it) {
 
      /* index of current and next granule in arrays of granLon/granLat/granGeoInMem */
      curr = curr % 2;          /* 0/1 */
      next = (curr+1) % 2;
      granGeoInMem[next] = false;
      
      /* get lon/lat boundaries [westGranLon, eastGranLon, northGranLat, southGranLat] for current granule */
      float cornerLon[4]; 
      float cornerLat[4];
      status = getGranGeo(granGeoInMem[curr], it->second, granLon[curr], granLat[curr],
                          cornerLon, cornerLat, num_row[curr], num_col[curr]);
      if (status == PROC_FAIL) {
         curr++; 
         continue;
      }
      getGranBoundary(cornerLon, cornerLat, westGranLon, eastGranLon, 
                      northGranLat, southGranLat, asndNode[curr], crsIDL[curr]);
   
      /* check the availability of next adjacent granule (less than 8 minutes) */
      nt = it;
      nt++;
      granNextAdj = false;  
      if (nt != granInfo.end()) {
         if ((nt->first - it->first) < 8./60.) granNextAdj = true;
      }
     
      /* loop over each station */
      for (vector<OneDayAerAod>::iterator ia=aerData.begin(); ia!=aerData.end(); ++ia) {
         /* no need to proceed if current station has been matched with previous granule */
         if (prevMatched(matchData, (*ia).staName, it->first)) continue;
        
         /* no need to proceed if current station is out of the granule domain (spatially and temporally) */ 
         if (!staMatchGranule(*ia, it->first, crsIDL, 
                              westGranLon, eastGranLon, northGranLat, southGranLat,
                              westStaLon, eastStaLon, northStaLat, southStaLat))
            continue;
         
         /* find the match-up VIIRS pixels from current (and possibly next) granule */
         int granRowDmn[2][2], granColDmn[2][2];
         bool pixAvlFromNext = false;
         status = findMatchPixels(granLon[curr], granLat[curr], num_row[curr], num_col[curr], asndNode[curr], crsIDL[curr], 
                                  westStaLon, eastStaLon, northStaLat, southStaLat, 
                                  granRowDmn[0], granColDmn[0]);
         if (status < 0) continue;  /* no matched pixel found */
         if (status == 1 && granNextAdj) { /* cross granule */
            /* read the lon/lat for next granule, get the boundaries */
            if (!granGeoInMem[next]) {
               status = getGranGeo(granGeoInMem[next], nt->second, granLon[next], granLat[next],
                                   cornerLon, cornerLat, num_row[next], num_col[next]);
               if (status == PROC_SUCCEED) {
                  float dummy1, dummy2, dummy3, dummy4;
                  getGranBoundary(cornerLon, cornerLat, dummy1, dummy2, dummy3, dummy4,
                                  asndNode[next], crsIDL[next]);
               }
               else {
                  /* invalid granule */
                  granNextAdj = false;
                  granGeoInMem[next] = false;
                  granInfo.erase(nt);
               }
               
            }
            if (granNextAdj) {  
               status = findMatchPixels(granLon[next], granLat[next], num_row[next], num_col[next], asndNode[next], crsIDL[next],
                                       westStaLon, eastStaLon, northStaLat, southStaLat, 
                                       granRowDmn[1], granColDmn[1]);
               if (status >= 0) pixAvlFromNext = true;
            }
         }

         /* collect the VIIRS data for the matched pixels */    
         MatchupRecord mr;
         mr.nMeas = 0;
         mr.nPixs = 0;  
         memset(mr.staName, 0, 32);
         (*ia).staName.copy(mr.staName, (*ia).staName.length(), 0);
         mr.staLon = (*ia).staLon;
         mr.staLat = (*ia).staLat;
         
         int nRow = 0;
         int nCol = 0;
         if (pixAvlFromNext) { 
            mr.satTime = 0.5*(it->first + nt->first);
            nRow = (granRowDmn[0][1]-granRowDmn[0][0]+1) + 
                   (granRowDmn[1][1]-granRowDmn[1][0]+1);
            if (granColDmn[0][0] < granColDmn[1][0]) granColDmn[0][0] = granColDmn[1][0];
            else granColDmn[1][0] = granColDmn[0][0];
            if (granColDmn[1][1] > granColDmn[0][1]) granColDmn[0][1] = granColDmn[1][1];
            else granColDmn[1][1] = granColDmn[0][1];
         } 
         else {
            mr.satTime = it->first;
            nRow = granRowDmn[0][1]-granRowDmn[0][0]+1;
         }  
         nCol = granColDmn[0][1]-granColDmn[0][0]+1;
         int np = nRow*nCol;
         if (np < MIN_N_MATCH_PIX) continue;
         
         /* collect matchup VIIRS pixel data */
         allocMemForMatch(mr, np);
         
         string aodFile;
         int offset, ncol;
         int start[2] = {0, granColDmn[0][0]};
         int count[2] = {0, nCol};
         int stride[2] = {1,1};
         int block[2] = {1,1};
         //hsize_t block[2] = {1,1};
         for (int ig=0; ig<2; ig++) {
            float *glon, *glat;
            if (ig == 0) {
               offset = 0;
               start[0] = granRowDmn[0][0];
               count[0] = granRowDmn[0][1]-granRowDmn[0][0]+1;
               glon = granLon[curr];
               glat = granLat[curr];
               ncol = num_col[curr];
            }
            else if (pixAvlFromNext) {
               offset = count[0]*count[1];
               start[0] = granRowDmn[1][0];
               count[0] = granRowDmn[1][1]-granRowDmn[1][0]+1;
               glon = granLon[next];
               glat = granLat[next];
               ncol = num_col[next];
            }
            else
               continue;
            
            /* lon/lat already in memory */
            int k = offset;
            for (int i=granRowDmn[ig][0]; i<=granRowDmn[ig][1]; i++) {
               for (int j=granColDmn[ig][0]; j<=granColDmn[ig][1]; j++) {
                  *(mr.lon+k) = *(glon+i*ncol+j);
                  *(mr.lat+k) = *(glat+i*ncol+j);
                  k++;
               }
            }
            
            /* get AOD data (NetCDF) */
            if (ig == 0) aodFile = it->second;
            else aodFile = nt->second;
            status = readDBAod(aodFile, mr, offset, start, stride, count, block);
            //if (status == PROC_FAIL) continue;
            if (status == PROC_FAIL) {
               freeMemForMatch(mr);
               break;
            }
            
         }
         if (status == PROC_FAIL) continue;
         
         /* only collect pixels with AOD retrievals */
         status = collectAod (mr);
         if (status == PROC_FAIL) continue;
         
         //sshams added
         /* Calculate satellite AOD statistics */
         float sumSatAOD = 0.0;       // Sum of satellite AOD values
         float sumSatAODSquared = 0.0; // Sum of squared satellite AOD values
         int nPixs = 0;

         for (int i = 0; i < mr.nPixs; i++) {
            float aod = mr.aod550[i];
            if (aod > -1.0) {  // Only include valid AOD values
               sumSatAOD += aod;
               sumSatAODSquared += aod * aod;
               nPixs++;
              //  cout << mr.aod550[i] << " ";  // Debug: Print AOD values
            }
         }
        //  // Debug: Print intermediate values
        //  cout << "Number of valid pixels: " << nPixs << endl;
        //  cout << "Sum of AOD values: " << sumSatAOD << endl;
        //  cout << "Sum of squared AOD values: " << sumSatAODSquared << endl;

         if (nPixs > 0) {
            float meanSatAOD = sumSatAOD / nPixs;
            float stdSatAOD = sqrt((sumSatAODSquared / nPixs) - (meanSatAOD * meanSatAOD));
           //  cout << "mean: " << meanSatAOD << endl;
            // Debug: Print the number of valid measurements and the required threshold
            cout << "***Overlap is found***" <<endl;
            cout << "Station: " << mr.staName <<endl;
            cout << "  Satellite Overpass Time: " << mr.satTime << " hours" << endl;
           // << ", Required: " << NVLDAER << endl;
            mr.satStd550 = stdSatAOD;  // Store standard deviation of satellite AOD at 550 nm
            mr.satMean550 = meanSatAOD;  // Store mean satellite AOD at 550 nm
         } else {
            mr.satStd550 = -999.0;  // Assign an invalid value if no valid pixels are found
         }


         /* collect the AERONET data for the matched pixels */  
         float timeWindow[2] = {(mr.satTime-MATCH_TIME_WINDOW*0.5)/24., 
                                (mr.satTime+MATCH_TIME_WINDOW*0.5)/24.};
         int nMeas = 0;
         int startIdx = -1;
         float sumAOD = 0.0;  // **Added: Sum of AOD values**
         float sumAODSquared = 0.0;  // **Added: Sum of squared AOD values**
         for (int i=0; i<(*ia).nmeas; i++) {
            if (((*ia).meas+i)->time >= timeWindow[0] && ((*ia).meas+i)->time <= timeWindow[1]) {
               float aod = ((*ia).meas + i)->aod550;  // Use the interpolated AOD at 0.550 µm
               if (aod != -999.0 && !std::isnan(aod)) {  // Skip invalid values
                  sumAOD += aod;  // **Added: Accumulate AOD values**
                  sumAODSquared += aod * aod;  // **Added: Accumulate squared AOD values**
                  nMeas++;
                  if(startIdx < 0) startIdx = i;
               }
            }
            else if (((*ia).meas+i)->time > timeWindow[1]) break;
         }
         if (nMeas == 0) {
            cout << "matchup(): AERONET nMeas == 0 for "+(*ia).staName << endl;
            freeMemForMatch(mr);
            continue;
         }
         
         /* Calculate average and standard deviation */
         float avgAOD = sumAOD / nMeas;  // **Added: Calculate average AOD**
         float stdAOD = sqrt((sumAODSquared / nMeas) - (avgAOD * avgAOD));  // **Added: Calculate standard deviation**

         mr.aerMean550 = avgAOD;  // **Store average AERONET AOD at 550 nm**
         mr.aerStd550 = stdAOD;   // **Store standard deviation of AERONET AOD at 550 nm**   
         mr.nMeas = nMeas;
         mr.meas = new float[nMeas*(NUM_AERONET_WAVELENGTH+1)]; 
         memcpy(mr.meas, (*ia).meas+startIdx, sizeof(float)*nMeas*(NUM_AERONET_WAVELENGTH+1));       
         
         /* insert to the matching vector */
         matchData.push_back(mr);
        
      } /* end of station loop */
      
      curr++; 
      
   }  /* end of granule loop */
 
   for (int i=0; i<2; i++) {
      delete[] granLon[i];
      delete[] granLat[i];
   }
   
   if (matchData.size() == 0) return PROC_FAIL;
   else return PROC_SUCCEED;
}                   


/******************************************************************************/
  
void freeMemForMatch (MatchupRecord mr)
{
   if (mr.nPixs > 0) {
      delete[] mr.solzen;
      delete[] mr.relazi;
      delete[] mr.satzen;
      delete[] mr.sctang;
      
      delete[] mr.lon;    
      delete[] mr.lat;
      delete[] mr.aod550;
      delete[] mr.ae;
      delete[] mr.lndSea;
   } 
   
   if (mr.nMeas > 0)  delete[] mr.meas; 
} 
/******************************************************************************/
