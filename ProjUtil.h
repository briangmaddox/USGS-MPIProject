/***************************************************************************
                          ProjUtil.h  -  contains misc useful projection
                                         related functions
                             -------------------
    begin                : Fri Oct 13 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/

#ifndef PROJUTIL_H_
#define PROJUTIL_H_

#include <vector>
#include <algorithm>
#include "DRect.h"
#include "ProjectionMesh/ProjectionMesh.h"
#include "ProjectionIO/ProjectionReader.h"
#include "ProjectionIO/ProjectionWriter.h"
#include "ProjectorException.h"
#include "ProjectionParams.h"

//get projectionparams from projection
ProjectionParams getParams(ProjLib::Projection * proj) throw(ProjectorException);

//get a projection from ProjectionParams
ProjLib::Projection * SetProjection(const ProjectionParams & inParams) throw();

//get the datum from a string
ProjLib::DATUM GetDatum(std::string indatum) throw();

//get the unit from a
ProjLib::UNIT GetUnit(std::string inunit) throw();

//convert from one scale to another
MathLib::Point GetConvertedScale(double oldwidth, double oldheight,
                                   DRect outRec) throw ();

//convert to Packed Degree Minute seconds
double ConvertToDMS(double degrees) throw();

//convert from DMS
double ConvertFromDMS(double dms) throw();

//get minMax of two vecors
void getMinMax(std::vector<double>& array, double& min, double& max) throw();

//get SameScale converts linear units to perserve linear length (note:
//if the input scale is angluar then the lat and long at equator are used
MathLib::Point getSameScale(MathLib::Point inoldscale, ProjLib::Projection * in,
               ProjLib::Projection * out) throw(ProjectionException);

//Get the projection from a input file
Projection * SetProjection(std::string parameterfile) throw();

#endif

