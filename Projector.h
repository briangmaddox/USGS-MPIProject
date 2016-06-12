/***************************************************************************
                          Projector.h  -  description
                             -------------------
    begin                : Thu Oct 12 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/


#ifndef PROJECTOR_H
#define PROJECTOR_H

#include <iostream>
#include "ProjUtil.h"
#include "ImageLib/DOQImageIFile.h"
#include "ImageLib/GeoTIFFImageOFile.h"
#include "ImageLib/GeoTIFFImageIFile.h"
#include "ImageLib/LRUCacheManager.h"
#include "ProjectionIO/ProjectionReader.h"
#include "ProjectionIO/ProjectionWriter.h"
#include "ProjectionMesh/ProjectionMesh.h"
#include "DRect.h"
#include "ProjectorException.h"
#include "ProjectionParams.h"
#include "BaseProgress.h"


#define CACHESIZE 100    //default is to try to cache 100 mbs of memory


//Reprojection object, converts one file to another
class Projector
{
public:
  //Constructors and Destructor
  Projector();
  Projector(std::string & ininfile, 
            ProjLib::Projection * inoutproject,
            const std::string & inoutfile = std::string("out.tif"))
            throw (ProjectorException, std::bad_alloc);
  Projector(std::string & ininfile, 
            const ProjectionParams& inParams,
            const std::string & inoutfile = std::string("out.tif"))
            throw (ProjectorException, std::bad_alloc);
  virtual ~Projector();

  //sets
  virtual void setInputFile(std::string & ininfile) throw(ProjectorException);
  void setOutputProjection(ProjLib::Projection * inoutproject) 
    throw(ProjectorException);
  void setOutputProjection(const ProjectionParams& inParams) 
    throw(ProjectorException);
  void setOutputFileName(const std::string& inoutfile) throw();
  void setPmeshName(const int & inpmeshname) throw();
  void setPmeshSize(const long int & inpmeshsize) throw();
  void setOutputScale(const MathLib::Point & innewscale) throw(); 
  
  //This function allows the user to set the cache size
  //Default is ~100mb of cached scanlines
  //add 3/6/01 by Chris Bilderback
  void setCacheSize(const unsigned int & incachesize) throw();

  //This function allows the user to specify that the output tiff
  //should be compressed with packbits compression.
  //Default is false.
  void setPackBits(const bool & inpackbits) throw();
  
  


  //This function allows user to have the same geographical 
  //output scale as input
  void setSameScale(bool insamescale) throw();
  
  //gets
  std::string getInputFile() const throw();
  ProjLib::Projection * getOutputProjection() const throw();
  ProjectionParams getOutputProjectionParams() const throw();
  int getPmeshName() const throw();
  int getPmeshSize() const throw();
  unsigned int getCacheSize() const throw();
  bool getPackBits() const throw();

  //main function which runs the projection
  virtual void
    project(BaseProgress * progress = NULL) 
    throw(ProjectorException);

protected:
  //setup input file and input projections
  void setupInput(std::string & ininfile) throw(ProjectorException);
  
  //setup the output file
  void setupOutput(std::string & inoutfile) throw(ProjectorException);

  //setup the pmesh 
  PmeshLib::ProjectionMesh * setupForwardPmesh() throw();
  PmeshLib::ProjectionMesh * setupReversePmesh() throw();

  //getExtents function gets the new bounding rectangle for the new image
  void getExtents(PmeshLib::ProjectionMesh * pmesh) throw(ProjectorException);


  ProjIOLib::ProjectionReader reader;
  ProjIOLib::ProjectionWriter writer;
  Projection * fromprojection, * toprojection;
  USGSImageLib::ImageIFile * infile;
  USGSImageLib::ImageOFile* out;                //the output file
  USGSImageLib::CacheManager* cache;            //cache

    
  //metrics
  long int oldheight, oldwidth, newheight, newwidth;
  DRect inRect, outRect;
  MathLib::Point oldscale, newscale;
  int pmeshsize;                                //pmesh metrics
  int pmeshname;
  int photo, spp, bps;
  std::string outfile;                          //outputfilename
  ProjectionParams Params;
  bool samescale;
  unsigned int cachesize;                       //the cache size in mb
  bool packbits;                                //whether to use packbits  
};

#endif





