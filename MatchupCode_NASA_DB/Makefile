
HDF5_INC=/data/starfs1/libs/hdf5-1.8.13/include
HDF5_LIB=/data/starfs1/libs/hdf5-1.8.13/lib
NETCDF4_INC=/data/starfs1/libs/netcdf-4.3.2/include
NETCDF4_LIB=/data/starfs1/libs/netcdf-4.3.2/lib
LIB64=-L /usr/lib64 -lm
#LIB64=-L /data/data087/hliu/tools/lib64 -lnsl -lm

LOCAL_INC=./
INCS =  -I$(LOCAL_INC) -I${HDF5_INC} -I${NETCDF4_INC}
HDFLIBS = -L${HDF5_LIB} -Wl,-rpath ${HDF5_LIB} -lhdf5
NCLIBS = -L${NETCDF4_LIB} -Wl,-rpath,${NETCDF4_LIB} -lnetcdf
LIBS = ${HDFLIBS} ${LIB64} ${NCLIBS}
#CFLAGS = -Wall
CFLAGS =
DBGFLAG = -g -w
OPTFLAGS =
LDFLAGS =

#========================================
# Compile and Link
#
#CXX		= pgCC
CXX		= g++
CXX_COMPILE	= $(CXX) -c $(OPTFLAGS) $(DBGFLAG) $(CFLAGS) $(INCS) 
CXX_LINK 	= $(CXX) $(LDFLAGS) $(LIBS)

SRCS = match_functions.cpp Error_Handler.cpp match_main.cpp
OBJS = ${SRCS:.cpp=.o}
TARGET = ../Bin/match_viirs_db_n20

$(TARGET): $(OBJS)
	$(CXX) -o $@  $(OBJS) $(LIBS) 
.cpp.o:
	$(CXX) $(OPTFLAGS) $(DBGFLAG) $(CFLAGS) -o $@ $(INCS) -c $< 

clean:
	rm -f $(OBJS) $(MAINOBJ1) $(MAINOBJ2) $(MAINOBJ3) core 

veryclean: clean
	rm -f $(TARGET)
#************************************************

