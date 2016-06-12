/***************************************************************************
                          ProjUtil.cpp  -  Implementation for the misc
                                           projection functions
                             -------------------
    begin                : Fri Oct 13 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/

#ifndef PROJUTIL_CPP_
#define PROJUTIL_CPP_

#include "ProjUtil.h"

//****************************************************************************
ProjectionParams getParams(ProjLib::Projection * proj) 
  throw(ProjectorException)
{
  ProjectionParams ret; //the return structure

  try
  {
    //check to see that the projection is valid
    if (!proj)
      throw ProjectionException(PROJECTOR_INVALID_PROJECTION); //invalid projection


    //first get the projsys and datum
    ret.projtype = proj->getProjectionSystem();
    ret.datum = proj->getDatum();
    ret.unit =  proj->getUnit();
    switch(ret.projtype)
    {
    case UTM:
      ret.zone = dynamic_cast<UTMProjection *>(proj)->getZone();
      break;
    case SPCS:
      ret.zone = dynamic_cast<StatePlaneProjection *>(proj)->getZone();
      break;
    case ALBERS:
      {
        AlbersConicProjection* albers =
          dynamic_cast<AlbersConicProjection*>(proj);
        ret.StdParallel1 =
          ConvertFromDMS(albers->getFirstStandardParallel());
        ret.StdParallel2 =
          ConvertFromDMS(albers->getSecondStandardParallel());
        ret.CenterLong = ConvertFromDMS(albers->getCentralMeridian());
        ret.NatOriginLat = ConvertFromDMS(albers->getOriginLatitude());
        ret.FalseEasting = albers->getFalseEasting();
        ret.FalseNorthing = albers->getFalseNorthing();
        break;
      }
    case AZMEQD:
    case GNOMON:
    case LAMAZ:
    case ORTHO:
    case STEREO:
      {
        AzimuthalProjection* azim =
          dynamic_cast<AzimuthalProjection *>(proj);
        ret.CenterLong = ConvertFromDMS(azim->getCentralLongitude());
        ret.CenterLat = ConvertFromDMS(azim->getCentralLatitude());
        ret.FalseEasting = azim->getFalseEasting();
        ret.FalseNorthing = azim->getFalseNorthing();
        break;
      }
    case MILLER:
    case ROBIN:
    case SNSOID:
      {
        PseudocylindricalProjection* pc =
          dynamic_cast<PseudocylindricalProjection*>(proj);
        ret.CenterLong = ConvertFromDMS(pc->getCentralMeridian());
        ret.FalseEasting = pc->getFalseEasting();
        ret.FalseNorthing = pc->getFalseNorthing();
        break;
      }
    case EQUIDC:
      {
        EquidistantConicProjection* edc =
          dynamic_cast<EquidistantConicProjection*>(proj);
        ret.CenterLat = ConvertFromDMS(edc->getStandardParallel());
        ret.CenterLong = ConvertFromDMS(edc->getCentralMeridian());
        ret.NatOriginLat = ConvertFromDMS(edc->getOriginLatitude());
        ret.FalseEasting = edc->getFalseEasting();
        ret.FalseNorthing = edc->getFalseNorthing();
        ret.StdParallel1 = ConvertFromDMS(edc->getFirstStandardParallel());
        ret.StdParallel2 =
          ConvertFromDMS(edc->getSecondStandardParallel());
        break;
      }
    case EQRECT:
      {
        EquirectangularProjection* eqrect =
          dynamic_cast<EquirectangularProjection*>(proj);
        ret.CenterLat = ConvertFromDMS(eqrect->getTrueScaleLatitude());
        ret.CenterLong = ConvertFromDMS(eqrect->getCentralMeridian());
        ret.FalseEasting = eqrect->getFalseEasting();
        ret.FalseNorthing = eqrect->getFalseNorthing();
        break;
      }
    case HOM:
      {
        HotineObliqueMercatorProjection* homerc =
          dynamic_cast<HotineObliqueMercatorProjection*>(proj);

        // Make sure this is a type B HOM projection
        if ( 0.0 != homerc->getLatitude1()
             || 0.0 != homerc->getLatitude2()
             || 0.0 != homerc->getLongitude1()
             || 0.0 != homerc->getLongitude2()
             || 0.0 != homerc->getAzimuthPoint())
          return ret;

        ret.ScaleAtNatOrigin = homerc->getProjectionCenterScaleFactor();
        ret.AzimuthAngle = homerc->getAzimuthAngle();
        ret.CenterLong = ConvertFromDMS(homerc->getCentralMeridian());
        ret.CenterLat = ConvertFromDMS(homerc->getOriginLatitude());
        ret.FalseEasting = homerc->getFalseEasting();
        ret.FalseNorthing = homerc->getFalseNorthing();
        break;
      }
    case LAMCC:
      {
        LambertConformalConicProjection* lamcc =
          dynamic_cast<LambertConformalConicProjection*>(proj);
        ret.StdParallel1 =
          ConvertFromDMS(lamcc->getFirstStandardParallel());
        ret.StdParallel2 =
          ConvertFromDMS(lamcc->getSecondStandardParallel());
        ret.NatOriginLong = ConvertFromDMS(lamcc->getCentralMeridian());
        ret.FalseOriginLat = ConvertFromDMS(lamcc->getOriginLatitude());
        ret.FalseEasting = lamcc->getFalseEasting();
        ret.FalseNorthing = lamcc->getFalseNorthing();
        break;
      }
    case MERCAT:
      {
        MercatorProjection* merc =
          dynamic_cast<MercatorProjection*>(proj);
        ret.NatOriginLong = ConvertFromDMS(merc->getCentralMeridian());
        ret.NatOriginLat = ConvertFromDMS(merc->getOriginLatitude());
        ret.CenterEasting = merc->getFalseEasting();
        ret.CenterNorthing = merc->getFalseNorthing();
        break;
      }
    case POLYC:
      {
        PolyconicProjection* polycon =
          dynamic_cast<PolyconicProjection*>(proj);
        ret.CenterLong = ConvertFromDMS(polycon->getCentralMeridian());
        ret.CenterLat = ConvertFromDMS(polycon->getOriginLatitude());
        ret.FalseEasting = polycon->getFalseEasting();
        ret.FalseNorthing = polycon->getFalseNorthing();
        break;
      }
    case PS:
      {
        PolarStereographicProjection* polar =
          dynamic_cast<PolarStereographicProjection*>(proj);
        ret.StraightVertPoleLong =
          ConvertFromDMS(polar->getLongitudePole());
        ret.NatOriginLat = ConvertFromDMS(polar->getTrueScaleLatitude());
        ret.FalseEasting = polar->getFalseEasting();
        ret.FalseNorthing = polar->getFalseNorthing();
        break;
      }
    case ALASKA:
    case TM:
      {
        TransverseMercatorProjection* trans =
          dynamic_cast<TransverseMercatorProjection*>(proj);
        ret.ScaleAtNatOrigin = trans->getCentralMeridianScaleFactor();
        ret.CenterLong = ConvertFromDMS(trans->getCentralMeridian());
        ret.NatOriginLat = ConvertFromDMS(trans->getOriginLatitude());
        ret.FalseEasting = trans->getFalseEasting();
        ret.FalseNorthing = trans->getFalseNorthing();
        break;
      }
    case VGRINT:
      {
        VanDerGrintenProjection* vandg =
          dynamic_cast<VanDerGrintenProjection*>(proj);
        ret.CenterLat = ConvertFromDMS(vandg->getOriginLatitude());
        ret.CenterLong =ConvertFromDMS(vandg->getCentralMeridian());
        ret.FalseEasting = vandg->getFalseEasting();
        ret.FalseNorthing = vandg->getFalseNorthing();
        break;
      }
	  case GOOD: //good has no real params
      break;
    default:
      throw ProjectionException(PROJECTOR_INVALID_PROJECTION); //invalid projection
    }
    return ret;
  }
  catch (ProjectorException & temp)
  {
    throw temp;                                        //just rethrow
  }
  catch(...)
  {
    throw ProjectorException(PROJECTOR_ERROR_UNKOWN); //unkown error
  }

}

//****************************************************************************
ProjLib::Projection * SetProjection(const ProjectionParams & inParams) throw()
{
  Projection * proj = NULL;                   //the return projection


  try
    {
      if (inParams.projtype == ProjLib::GEO)
      {
        if (inParams.unit != ProjLib::ARC_DEGREES)
          return NULL;

        if (!(proj = new (std::nothrow) GeographicProjection
              (inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }

      if (inParams.projtype == ProjLib::UTM)
      {
        if (!(proj = new (std::nothrow) UTMProjection(inParams.zone,
                                                      inParams.datum,
                                                      inParams.unit)))
          throw std::bad_alloc();
      }

      if (inParams.projtype == ProjLib::SPCS)
      {
        if (inParams.datum == NAD83)
          StatePlaneProjection::setNAD83ParameterFilename
            (std::string("./nad83sp"));
        else
        {
          if (inParams.datum == NAD27)
            StatePlaneProjection::setNAD83ParameterFilename
            (std::string("./nad27sp"));
          else
          {
            return NULL;                      //not a valid state plane
          }
        }

        //create the state plane projection
        //here I am assuming that the zone will be
        //a correct state plane zone
        if (!(proj = new (std::nothrow) StatePlaneProjection
              (inParams.zone, inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }

      if (inParams.projtype == ProjLib::ALBERS)
      {
        if(!(proj =  new (std::nothrow) AlbersConicProjection
             ( ConvertToDMS(inParams.StdParallel1),
                 ConvertToDMS(inParams.StdParallel2),
               0.0, 0.0, ConvertToDMS(inParams.CenterLong),
               ConvertToDMS(inParams.NatOriginLat),
               inParams.FalseEasting, inParams.FalseNorthing,
               inParams.datum,
               inParams.unit)))
          throw std::bad_alloc();


      }
      if (inParams.projtype == ProjLib::AZMEQD)
      {
        if(!(proj = new (std::nothrow) AzimuthalEquidistantProjection
             ( ConvertToDMS(inParams.CenterLong),
               ConvertToDMS(inParams.CenterLat),
               inParams.FalseEasting, inParams.FalseNorthing,
               0.0, inParams.datum,
               inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::GNOMON)
      {

        if(!(proj = new(std::nothrow) GnomonicProjection
             ( ConvertToDMS(inParams.CenterLong),
               ConvertToDMS(inParams.CenterLat),
               inParams.FalseEasting,inParams.FalseNorthing,
               0.0, inParams.datum, inParams.unit)))
            throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::LAMAZ)
      {
        if(!(proj = new(std::nothrow) LambertAzimuthalProjection
             ( ConvertToDMS(inParams.CenterLong),
               ConvertToDMS(inParams.CenterLat),
               inParams.FalseEasting, inParams.FalseNorthing, 0.0,
               inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::ORTHO)
      {
        if(!(proj = new(std::nothrow) OrthographicProjection
             ( ConvertToDMS(inParams.CenterLong),
               ConvertToDMS(inParams.CenterLat),
               inParams.FalseEasting,inParams.FalseNorthing, 0.0,
               inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::STEREO)
      {
        if(!(proj = new(std::nothrow) StereographicProjection
             (ConvertToDMS(inParams.CenterLong),
              ConvertToDMS(inParams.CenterLat),
              inParams.FalseEasting, inParams.FalseNorthing,
              0.0, inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::MILLER)
      {
        if(!(proj = new (std::nothrow) MillerCylindricalProjection
             ( 0.0, ConvertToDMS(inParams.CenterLong),
               inParams.FalseEasting, inParams.FalseNorthing,
               inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::ROBIN)
      {
        if(!(proj = new(std::nothrow) RobinsonProjection
             ( 0.0, ConvertToDMS(inParams.CenterLong),
               inParams.FalseEasting, inParams.FalseNorthing,
               inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::SNSOID)
      {
        if(!(proj = new(std::nothrow) SinusoidalProjection
             ( 0.0, ConvertToDMS(inParams.CenterLong),
               inParams.FalseEasting,
               inParams.FalseNorthing, inParams.datum,
               inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::EQUIDC)
      {
        if ( inParams.StdParallel1 == inParams.StdParallel2 )
        {
          if(!(proj =  new(std::nothrow) EquidistantConicProjection
               ( ConvertToDMS(inParams.CenterLat), 0.0, 0.0,
                 ConvertToDMS(inParams.CenterLong),
                 ConvertToDMS(inParams.NatOriginLat),
                 inParams.FalseEasting,
                 inParams.FalseNorthing,
                 inParams.datum, inParams.unit)))
            throw std::bad_alloc();

        }
        else
        {
          if(!(proj = new(std::nothrow) EquidistantConicProjection
               ( ConvertToDMS(inParams.StdParallel1),
                 ConvertToDMS(inParams.StdParallel2),
                 0.0, 0.0,
                 ConvertToDMS(inParams.CenterLong),
                 ConvertToDMS(inParams.NatOriginLat),
                 inParams.FalseEasting,
                 inParams.FalseNorthing,
                 inParams.datum, inParams.unit)))
            throw std::bad_alloc();
        }
      }
      if (inParams.projtype == ProjLib::EQRECT)
      {
        if(!(proj =  new(std::nothrow) EquirectangularProjection
             ( ConvertToDMS(inParams.CenterLat), 0.0,
               ConvertToDMS(inParams.CenterLong),
               inParams.FalseEasting,inParams.FalseNorthing,
               inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::HOM)
      {
        if(!(proj = new (std::nothrow) HotineObliqueMercatorProjection
             ( inParams.ScaleAtNatOrigin, inParams.AzimuthAngle,
               0.0, 0.0, ConvertToDMS(inParams.CenterLong),
               ConvertToDMS(inParams.CenterLat), inParams.FalseEasting,
               inParams.FalseNorthing, inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::LAMCC)
      {
        if(!(proj = new(std::nothrow) LambertConformalConicProjection
             ( ConvertToDMS(inParams.StdParallel1),
               ConvertToDMS(inParams.StdParallel2),
               0.0, 0.0, ConvertToDMS(inParams.NatOriginLong),
               ConvertToDMS(inParams.FalseOriginLat),
               inParams.FalseEasting,
               inParams.FalseNorthing, inParams.datum,
               inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::MERCAT)
      {
        if(!(proj = new(std::nothrow) MercatorProjection
             ( 0.0, 0.0, ConvertToDMS(inParams.NatOriginLong),
               ConvertToDMS(inParams.NatOriginLat),
               inParams.CenterEasting, inParams.CenterNorthing,
               inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::POLYC)
      {
        if(!(proj = new(std::nothrow) PolyconicProjection
             ( 0.0, 0.0, ConvertToDMS(inParams.CenterLong),
               ConvertToDMS(inParams.CenterLat),
               inParams.FalseEasting,
               inParams.FalseNorthing, inParams.datum,
               inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::PS)
      {
        if(!(proj = new(std::nothrow) PolarStereographicProjection
             (ConvertToDMS(inParams.StraightVertPoleLong),
              ConvertToDMS(inParams.NatOriginLat), 0.0, 0.0,
              inParams.FalseEasting, inParams.FalseNorthing,
              inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::ALASKA)
      {
        if(!(proj =  new(std::nothrow) AlaskaConformalProjection
             (0.0, 0.0, inParams.FalseEasting,
              inParams.FalseNorthing,
              inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::TM)
      {
        if(!(proj = new(std::nothrow) TransverseMercatorProjection
             (inParams.ScaleAtNatOrigin, 0.0, 0.0,
              ConvertToDMS(inParams.CenterLong),
              ConvertToDMS(inParams.NatOriginLat),
              inParams.FalseEasting,
              inParams.FalseNorthing,
              inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::VGRINT)
      {
        if(!(proj = new(std::nothrow) VanDerGrintenProjection
             ( ConvertToDMS(inParams.CenterLat), 0.0,
               ConvertToDMS(inParams.CenterLong),
               inParams.FalseEasting,
               inParams.FalseNorthing,
               inParams.datum, inParams.unit)))
          throw std::bad_alloc();
      }
      if (inParams.projtype == ProjLib::GOOD)
      {
        if (!(proj = new (std::nothrow)
              GoodeHomolosineProjection(0.0, inParams.datum,
                                        inParams.unit)))
          throw std::bad_alloc();
      }
      return proj;
    }
  catch(...)
  {
    return NULL;
  }
}


//********************************************************************
ProjLib::DATUM GetDatum(std::string indatum) throw()
{

   if (indatum == std::string("ADINDAN"))       //see which datum it is
     return ADINDAN;

   if (indatum == std::string("ARC1950"))
     return ARC1950;

   if (indatum == std::string("ARC1960"))
     return ARC1960;

   if (indatum == std::string("AUSTRALIAN_GEODETIC_1966"))
     return AUSTRALIAN_GEODETIC_1966;

   if (indatum == std::string("AUSTRALIAN_GEODETIC_1984"))
     return AUSTRALIAN_GEODETIC_1984;

   if (indatum == std::string("CAPE"))
     return CAPE;

   if (indatum == std::string("EUROPEAN_DATUM_1950"))
     return EUROPEAN_DATUM_1950;

   if (indatum == std::string("HU_TZU_SHAN"))
     return HU_TZU_SHAN;

   if (indatum == std::string("INDIAN"))
     return INDIAN;

   if (indatum == std::string("NAD27"))
     return NAD27;

   if (indatum == std::string("NAD83"))
     return NAD83;

   if (indatum == std::string("ORDNANCE_SURVEY_1936"))
     return ORDNANCE_SURVEY_1936;

   if (indatum == std::string("PULKOVO_1942"))
     return PULKOVO_1942;

   if (indatum == std::string("PROVISIONAL_S_AMERICAN_1956"))
     return PROVISIONAL_S_AMERICAN_1956;

   if (indatum == std::string("TOKYO"))
     return TOKYO;

   if (indatum == std::string("WGS_72"))
     return WGS_72;

   if (indatum == std::string("WGS_84"))
     return WGS_84;

   return UNKNOWN_DATUM;                    //error finding datum
}

//******************************************************************
ProjLib::UNIT GetUnit(std::string inunit) throw()
{
  if (inunit == std::string("METERS"))          //see which unit it is in
    return METERS;

  if (inunit == std::string("ARC_DEGREES"))
    return ARC_DEGREES;

  if ((inunit == std::string("FEET")) || (inunit == std::string("US_FEET")))
    return US_FEET;

  return UNKNOWN_UNIT;                      //error finding projection

}



//********************************************************************
MathLib::Point GetConvertedScale(double oldwidth, double oldheight,
                                   DRect outRec) throw ()
{
  MathLib::Point ret;               //return value

  //try to maintain the same number of pixels in each direction
  ret.x = (outRec.right - outRec.left)/oldwidth;
  ret.y = (outRec.top - outRec.bottom)/oldheight;
  return ret;
}

// **************************************************************************
double ConvertToDMS(double degrees) throw()
{
  double temp;
  int deg, min;

  int sign = 1;

  temp = degrees;
  if (degrees < 0.0)             //make positive
  {
    sign = -1;
    temp = -temp;
  }

  deg = static_cast<int>(temp); //get the degrees
  temp -= deg;
  temp *= 60.0;
  min = static_cast<int>(temp); //get the minutes
  temp -= min;
  temp *= 60;                   //get the seconds

  temp = deg * 1000000 + min * 1000 + temp;
  return temp*sign;             //return approriate sign
}


//********************************************************
double ConvertFromDMS(double dms) throw()
{
  double sign = 1.0;
  double temp = dms;
  int deg;
  int min;
  if ( temp < 0.0 )
    {
      temp = -temp;
      sign = -1.0;
    }

  // Extract the degrees, minutes, and seconds
  deg = static_cast<int>(temp / 1000000);
  temp -= deg * 1000000;
  min = static_cast<int>(temp / 1000);
  temp -= min * 1000;
  temp /= 3600;
  temp += min / 60.0;
  temp += deg;
  return temp *sign;        //return apropriate sign
}

//***************************************************************************
void getMinMax(std::vector<double>& array, double& min, double& max) throw()
{
  std::sort(array.begin(), array.end());

  max = array[array.size() - 1];
  min = array[0];
}

//****************************************************************************
MathLib::Point getSameScale(MathLib::Point inoldscale,
               ProjLib::Projection * in,
               ProjLib::Projection * out) throw(ProjectionException)
{
  ProjLib::UNIT inunit;     //tempory input unit
  ProjLib::UNIT outunit;    //tempory output unit
  MathLib::Point ret;       //return point

  inunit = in->getUnit();
  outunit = out->getUnit();

  if (inunit != ProjLib::METERS &&
      inunit != ProjLib::US_FEET &&
      inunit != ProjLib::ARC_DEGREES)
    throw ProjectionException(PROJECTOR_INVALID_UNIT);

  if (outunit != ProjLib::METERS &&
      outunit != ProjLib::US_FEET   &&
      outunit != ProjLib::ARC_DEGREES)
    throw ProjectionException(PROJECTOR_INVALID_UNIT);



  //run through the cases
  if (inunit == ProjLib::METERS)
  {
    if (outunit == ProjLib::METERS)
      ret = inoldscale;

    if (outunit == ProjLib::US_FEET)
    {
      ret.x = inoldscale.x*0.3048;
      ret.y = inoldscale.y*0.3048;
    }

    if (outunit == ProjLib::ARC_DEGREES)
    {
      ret.x = inoldscale.x/110574.0;
      ret.y = inoldscale.y/111319.0;
    }
  }

  if (inunit == ProjLib::US_FEET)
  {
    if (outunit == ProjLib::METERS)
    {
      ret.x = inoldscale.x/0.3048;
      ret.y = inoldscale.y/0.3048;
    }
    if (outunit == ProjLib::US_FEET)
      ret = inoldscale;

    if (outunit == ProjLib::ARC_DEGREES)
    {
      ret.x = inoldscale.x/362776.0;
      ret.y = inoldscale.y/365221.0;
    }
  }

  if (inunit == ProjLib::ARC_DEGREES)
  {
    if (outunit == ProjLib::METERS)
    {
      ret.x = inoldscale.x * 110574.0;
      ret.y = inoldscale.y * 111319.0;
    }
    if (outunit == ProjLib::US_FEET)
    {
      ret.x = inoldscale.x*362776.0;
      ret.y = inoldscale.y*365221.0;
    }

    if (outunit == ProjLib::ARC_DEGREES)
    {
      ret = inoldscale;
    }
  }
  return ret;
}

//***********************************************************
Projection * SetProjection(std::string parameterfile) throw()
{

  std::ifstream in;                            //parameter file input
  std::string projtype, sdatum, sunit;         //projection type, datum and unit
  double StdParallel1 = 0.0;                   //projection options
  double StdParallel2 = 0.0;
  double NatOriginLong = 0.0;
  double NatOriginLat = 0.0;
  double FalseEasting = 0.0;
  double FalseNorthing = 0.0;
  //double FalseOriginLong = 0.0;
  double FalseOriginLat = 0.0;
  //double FalseOriginEasting = 0.0;
  //double FalseOriginNorthing = 0.0;
  double CenterLong = 0.0;
  double CenterLat = 0.0;
  double CenterEasting = 0.0;
  double CenterNorthing = 0.0;
  double ScaleAtNatOrigin = 1.0;
  double AzimuthAngle = 0.0;
  double StraightVertPoleLong = 0.0;
  int zone = 0;
  Projection * proj = NULL;                     //the projection to return

  try
    {
      in.open(parameterfile.c_str());           //open the parameter file


      in >> projtype;                           //read the projection type

      if (projtype == std::string("GEO"))
        {
          in >> sdatum;
          in >> sunit;
          if (!(proj = new (std::nothrow) GeographicProjection
                (GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }

      if (projtype == std::string("UTM"))
        {
          in >> zone;
          in >> sdatum;
          in >> sunit;


          if (!(proj = new (std::nothrow) UTMProjection(zone,
                                                        GetDatum(sdatum),
                                                        GetUnit(sunit))))
            throw std::bad_alloc();
        }

      if (projtype == std::string("SPCS"))
        {
          in >> zone;
          in >> sdatum;
          in >> sunit;

          if (sdatum == std::string("NAD83"))
            StatePlaneProjection::setNAD83ParameterFilename
              (std::string("./nad83sp"));
          else
          {
            if (sdatum == std::string("NAD27"))
              StatePlaneProjection::setNAD83ParameterFilename
                (std::string("./nad27sp"));
            else
            {
              return NULL;                      //not a valid state plane
            }
          }
          //create the state plane projection
          //here I am assuming that the zone will be
          //a correct state plane zone
          if (!(proj = new (std::nothrow) StatePlaneProjection
                (zone, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }

      if (projtype == std::string("ALBERS"))
        {
          in >> sdatum;
          in >> sunit;
          in >> StdParallel1
             >> StdParallel2
             >> CenterLong
             >> NatOriginLat
             >> FalseEasting
             >> FalseNorthing;

          if(!(proj =  new (std::nothrow) AlbersConicProjection
               ( ConvertToDMS(StdParallel1),
                 ConvertToDMS(StdParallel2),
                 0.0, 0.0, ConvertToDMS(CenterLong),
                 ConvertToDMS(NatOriginLat),
                 FalseEasting, FalseNorthing, GetDatum(sdatum),
                 GetUnit(sunit))))
            throw std::bad_alloc();


        }
      if (projtype == std::string("AZMEQD"))
        {
          in >> sdatum >> sunit >> CenterLong >> CenterLat
             >> FalseEasting >> FalseNorthing;


          if(!(proj = new (std::nothrow) AzimuthalEquidistantProjection
               ( ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat),
                 FalseEasting, FalseNorthing, 0.0, GetDatum(sdatum),
                 GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("GNOMON"))
        {
          in >> sdatum >> sunit >> CenterLong >> CenterLat
             >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) GnomonicProjection
               ( ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat),
                 FalseEasting,FalseNorthing,
                 0.0, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("LAMAZ"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> CenterLat >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) LambertAzimuthalProjection
               ( ConvertToDMS(CenterLong), ConvertToDMS(CenterLat),
                 FalseEasting, FalseNorthing, 0.0,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("ORTHO"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> CenterLat >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) OrthographicProjection
               ( ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat),
                 FalseEasting,FalseNorthing, 0.0,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("STEREO"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> CenterLat >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) StereographicProjection
               (ConvertToDMS(CenterLong),
                ConvertToDMS(CenterLat),
                FalseEasting, FalseNorthing,
                0.0, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("MILLER"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> FalseEasting >> FalseNorthing;

          if(!(proj = new (std::nothrow) MillerCylindricalProjection
               ( 0.0, ConvertToDMS(CenterLong),
                 FalseEasting, FalseNorthing,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("ROBIN"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) RobinsonProjection
               ( 0.0, ConvertToDMS(CenterLong),
                 FalseEasting, FalseNorthing,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("SNSOID"))
        {
          in >> sdatum >> sunit >> CenterLong
             >> FalseEasting >> FalseNorthing;
          if(!(proj = new(std::nothrow) SinusoidalProjection
               ( 0.0, ConvertToDMS(CenterLong),
                 FalseEasting,
                 FalseNorthing, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("EQUIDC"))
        {
          in >> sdatum >> sunit >> StdParallel1
             >> StdParallel2;

          if ( StdParallel1 == StdParallel2 )
            {
              in >> CenterLat >> CenterLong >> NatOriginLat
                 >> FalseEasting >> FalseNorthing;
              if(!(proj =  new(std::nothrow) EquidistantConicProjection
                   ( ConvertToDMS(CenterLat), 0.0, 0.0,
                     ConvertToDMS(CenterLong),
                     ConvertToDMS(NatOriginLat),
                     FalseEasting,
                     FalseNorthing,
                     GetDatum(sdatum), GetUnit(sunit))))
                throw std::bad_alloc();

            }
          else
            {
              in >> CenterLong >> NatOriginLat >> FalseEasting
                 >> FalseNorthing;
              if(!(proj = new(std::nothrow) EquidistantConicProjection
                   ( ConvertToDMS(StdParallel1),
                     ConvertToDMS(StdParallel2),
                     0.0, 0.0,
                     ConvertToDMS(CenterLong),
                     ConvertToDMS(NatOriginLat),
                     FalseEasting,
                     FalseNorthing,
                     GetDatum(sdatum), GetUnit(sunit))))
                throw std::bad_alloc();
            }
        }
      if (projtype == std::string("EQRECT"))
        {
          in >> sdatum >> sunit >> CenterLat >> CenterLong
             >> FalseEasting >> FalseNorthing;

          if(!(proj =  new(std::nothrow) EquirectangularProjection
               ( ConvertToDMS(CenterLat), 0.0,
                 ConvertToDMS(CenterLong),
                 FalseEasting,FalseNorthing,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("HOM"))
        {
          in >> sdatum >> sunit >> ScaleAtNatOrigin
             >> AzimuthAngle >> CenterLong
             >> CenterLat >> FalseEasting
             >> FalseNorthing;

          if(!(proj = new (std::nothrow) HotineObliqueMercatorProjection
               ( ScaleAtNatOrigin, AzimuthAngle,
                 0.0, 0.0, ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat), FalseEasting,
                 FalseNorthing, GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("LAMCC"))
        {
          in >> sdatum >> sunit >> StdParallel1 >> StdParallel2
             >> NatOriginLong >> FalseOriginLat
             >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) LambertConformalConicProjection
               ( ConvertToDMS(StdParallel1), ConvertToDMS(StdParallel2),
                 0.0, 0.0, ConvertToDMS(NatOriginLong),
                 ConvertToDMS(FalseOriginLat),
                 FalseEasting,FalseNorthing, GetDatum(sdatum),
                 GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("MERCAT"))
        {
          in >> sdatum >> sunit >> NatOriginLong
             >> NatOriginLat >> CenterEasting
             >> CenterNorthing;

          if(!(proj = new(std::nothrow) MercatorProjection
               ( 0.0, 0.0, ConvertToDMS(NatOriginLong),
                 ConvertToDMS(NatOriginLat),
                 CenterEasting, CenterNorthing,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("POLYC"))
        {

          in >> sdatum >> sunit >> CenterLong
             >> CenterLat >> FalseEasting
             >> FalseNorthing;

          if(!(proj = new(std::nothrow) PolyconicProjection
               ( 0.0, 0.0, ConvertToDMS(CenterLong),
                 ConvertToDMS(CenterLat), FalseEasting,
                 FalseNorthing, GetDatum(sdatum),
                 GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("PS"))
        {
          in >> sdatum >> sunit >> StraightVertPoleLong
             >> NatOriginLat >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) PolarStereographicProjection
               (ConvertToDMS(StraightVertPoleLong),
                ConvertToDMS(NatOriginLat), 0.0, 0.0,
                FalseEasting,FalseNorthing,
                GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("ALASKA"))
        {
          in >> sdatum >> sunit >> FalseEasting
             >> FalseNorthing;

          if(!(proj =  new(std::nothrow) AlaskaConformalProjection
                   (0.0, 0.0, FalseEasting, FalseNorthing,
                    GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }
      if (projtype == std::string("TM"))
        {
          in >> sdatum >> sunit >> ScaleAtNatOrigin
             >> CenterLong >> NatOriginLat
             >> FalseEasting >> FalseNorthing;

          if(!(proj = new(std::nothrow) TransverseMercatorProjection
               (ScaleAtNatOrigin, 0.0, 0.0,
                ConvertToDMS(CenterLong),
                ConvertToDMS(NatOriginLat), FalseEasting,
                FalseNorthing,
                GetDatum(sdatum), GetUnit(sunit))))
              throw std::bad_alloc();
        }
      if (projtype == std::string("VGRINT"))
        {
          in >> sdatum >> sunit >> CenterLat
             >> CenterLong >> FalseEasting
             >> FalseNorthing;
          if(!(proj = new(std::nothrow) VanDerGrintenProjection
               ( ConvertToDMS(CenterLat), 0.0,
                 ConvertToDMS(CenterLong),
                 FalseEasting,FalseNorthing,
                 GetDatum(sdatum), GetUnit(sunit))))
            throw std::bad_alloc();
        }

	    if (projtype == std::string("GOOD"))
	    {
		    in >> sdatum >> sunit;
		    if (!(proj = new (std::nothrow) GoodeHomolosineProjection(0.0,
			  GetDatum(sdatum), GetUnit(sunit))))
			    throw std::bad_alloc();
	    }


      in.close();
      if (!proj)  //unable to create the projection
      {
        return NULL;
      }

      return proj;
    }
  catch(...)
    {
      in.close();
      if (proj)
        delete proj;
      return NULL;
    }
}

#endif
