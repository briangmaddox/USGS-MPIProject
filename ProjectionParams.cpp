/***************************************************************************
                          ProjectionParams.cpp  -  description
                             -------------------
    begin                : Thu Oct 12 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/


/***************************************************************************
                          ProjectionParams.cpp  -  description
                             -------------------
    begin                : Wed Oct 11 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/


#ifndef PROJECTIONPARAMS_CPP
#define PROJECTIONPARAMS_CPP


#include "ProjectionParams.h"


//**********************************************************************
ProjectionParams::ProjectionParams()
{
  projtype = ProjLib::GEO;
  datum = ProjLib::NAD83;
  unit = ProjLib::ARC_DEGREES;
  StdParallel1 = 0.0;
  StdParallel2 = 0.0;
  NatOriginLong = 0.0;
  NatOriginLat = 0.0;
  FalseOriginLong = 0.0;
  FalseOriginLat = 0.0;
  FalseOriginEasting = 0.0;
  FalseOriginNorthing = 0.0;
  CenterLong = 0.0;
  CenterLat = 0.0;
  CenterEasting = 0.0;
  CenterNorthing = 0.0;
  ScaleAtNatOrigin = 0.0;
  AzimuthAngle = 0.0;
  StraightVertPoleLong = 0.0;
  zone = 0;
  FalseEasting = 0.0;
  FalseNorthing = 0.0;
}

#endif
