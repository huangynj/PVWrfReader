# Cmake file to build the paraview wrf reader
# Dave Semeraro - Texas Advanced Computing Center
#
# 
cmake_minimum_required(VERSION 2.8)
# Append the package module path
# 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake")
# Find the paraview package
# 
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
# Find GDAL stuff - GDAL is a georef package
# Note - CMake comes with a FindGDAL module
# and that module looks for the GDAL_DIR 
# environmnt variable. Set it to point to 
# the directory you installed GDAL in
FIND_PACKAGE(GDAL REQUIRED)
INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
# Find the Netcdf packages
# I am going to grab the FindNetCDF.cmake module
# from VTK for this. The module is in CMake dir
FIND_PACKAGE(NetCDF REQUIRED)
INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIR})
# Find the PROJ stuff. Going to use
# Find_path and find_library here cause there
# is no module for PROJ that I know of.
#FIND_LIBRARY(PROJ_LIBRARY proj )
SET(PROJ_LIBRARY "/work/01197/semeraro/maverick/Packages/PROJ/lib/libproj.so")

ADD_PARAVIEW_PLUGIN(WRFReader "1.0" SERVER_MANAGER_XML WRFReader.xml SERVER_MANAGER_SOURCES WRFReader.cxx ${GDAL_LIBRARY} ${NETCDF_LIBRARY} ${PROJ_LIBRARY} REQUIRED_ON_SERVER)
TARGET_LINK_LIBRARIES(WRFReader PRIVATE  ${GDAL_LIBRARY} ${NETCDF_LIBRARY} ${PROJ_LIBRARY})


