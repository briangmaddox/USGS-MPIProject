/***************************************************************************
                          Projector.cpp  -  description
                             -------------------
    begin                : Thu Oct 12 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/


#ifndef PROJECTOR_CPP
#define PROJECTOR_CPP

#ifdef _WIN32
#pragma warning( disable : 4291 ) // Disable VC warning messages for
                                  // new(nothrow)
#endif


#include <strstream>
#include "Projector.h"
#include "ImageLib/RGBPalette.h"
#include <fstream>

//*********************************************************************
Projector::Projector() : fromprojection(NULL), toprojection(NULL),
infile(NULL), out(NULL), cache(NULL), 
oldheight(0), oldwidth(0), newheight(0), newwidth(0),
pmeshsize(4), pmeshname(0), outfile("out.tif"), 
samescale(false), cachesize(CACHESIZE), packbits(false) 
{
  //init the scales
  oldscale.x = newscale.x = 0;
  oldscale.y = newscale.y = 0;
  
}
  
//***********************************************************************
Projector::Projector(std::string & ininfile, 
            ProjLib::Projection * inoutproject,
            const std::string & inoutfile) 
  throw (ProjectorException, std::bad_alloc)
  : fromprojection(NULL),
    toprojection(NULL), infile(NULL), out(NULL), cache(NULL), 
    oldheight(0), 
    oldwidth(0), newheight(0), newwidth(0),
    pmeshsize(4), pmeshname(0), outfile("out.tif"), samescale(false),
    cachesize(CACHESIZE), packbits(false)
{
  oldscale.x = newscale.x = 0;                //initialize scale
  oldscale.y = newscale.y = 0;

  setInputFile(ininfile);                     //open the input file
  toprojection = inoutproject;                //set the to project
  Params = getParams(toprojection);           //get the Params

  outfile = inoutfile;                        //get the output filename
  
}
  
//********************************************************************
Projector::Projector(std::string & ininfile, 
            const ProjectionParams&  inParams,
            const std::string & inoutfile) 
  throw (ProjectorException, std::bad_alloc)
  : fromprojection(NULL),
    toprojection(NULL), infile(NULL), out(NULL), cache(NULL), 
    oldheight(0), 
    oldwidth(0), newheight(0), newwidth(0),
    pmeshsize(4), pmeshname(0), outfile("out.tif"), samescale(false),
    cachesize(CACHESIZE), packbits(false)
{
  oldscale.x = newscale.x = 0;                //initialize scale
  oldscale.y = newscale.y = 0;
  
  setInputFile(ininfile);                     //open the input file
  Params = inParams;                          //get the Params
  toprojection = SetProjection(Params);       //get the 
  outfile = inoutfile;                        //get the output filename
  
}

//***********************************************************************
Projector::~Projector()
{
  delete cache;
  delete toprojection;
  delete infile;
}
  

//**********************************************************************
void Projector::setInputFile(std::string & ininfile ) throw(ProjectorException)
{
  setupInput(ininfile); //setup the input file and params
}

//*************************************************************************
void Projector::
setOutputProjection(ProjLib::Projection * inoutproject) 
  throw(ProjectorException)
{
  if (!inoutproject)                 //check to see that its not null
    throw ProjectorException(PROJECTOR_INVALID_PROJECTION);

  //ADDED to fix memory passing problems -CBB 8/1/01
  //clone the projection
  toprojection = inoutproject->clone();
  Params = getParams(inoutproject);  //fill in the projection params struct
}

//*************************************************************************
void Projector::setOutputProjection(const ProjectionParams& inParams) 
  throw(ProjectorException)
{
  Params = inParams;
  toprojection = SetProjection(Params);  //create the projection from the 
                                         //params

  if (!toprojection)                     //check for valid projection
    throw ProjectorException(PROJECTOR_INVALID_PROJECTION);

}

//************************************************************************
void Projector::setOutputFileName(const std::string& inoutfile) throw()
{
  outfile = inoutfile; //get a different output filename
}

//************************************************************************
void Projector::setPmeshName(const int& inpmeshname) throw()
{
  pmeshname = inpmeshname; //set the pmesh name
}

//************************************************************************
void Projector::setPmeshSize(const long int & inpmeshsize) throw()
{
  pmeshsize = inpmeshsize;  //set the pmesh size
}

//***********************************************************************
void Projector::setOutputScale(const MathLib::Point & innewscale) throw()
{
  newscale = innewscale;   //set a new scale instead of apsect perserving
}


//***********************************************************************
void Projector::setSameScale(bool insamescale) throw()
{
  samescale = insamescale;
}

//*********************************************************
void Projector::setCacheSize(const unsigned int & incachesize) throw()
{
  cachesize = incachesize;
}


//**********************************************************************
ProjLib::Projection * Projector::getOutputProjection() const throw()
{
  return toprojection;    //give up the output projection
}

//**********************************************************************
ProjectionParams Projector::getOutputProjectionParams() const throw()
{
  return Params;         //give up the output Params
}

//**********************************************************************
int Projector::getPmeshName() const throw()
{
  return pmeshname;      //give the pmesh type
}

//**************************************************************
int Projector::getPmeshSize() const throw()
{
  return pmeshsize;     //return the square of the number of nodes
}

//**************************************************************
unsigned int Projector::getCacheSize() const throw()
{
  return cachesize;
}

//**************************************************************
void Projector::setPackBits(const bool & inpackbits) throw()
{
  packbits = inpackbits;
}


//**************************************************************
bool Projector::getPackBits() const throw()
{
  return packbits;
}
  
//**************************************************************
void Projector::project(BaseProgress * progress)
throw(ProjectorException)
{
  double x, y;                                 //temp variables
  long int _x, _y;                             //for world to pixel translation
  tdata_t scanline = NULL;                     //output scanline
  tdata_t inscanline = NULL;                   //input scanline
  const unsigned char * cachescanline = NULL;  //scanline from cache
  int sppcounter;                              //for colored images
  PmeshLib::ProjectionMesh * pmesh = NULL;     //projection mesh
  long int xcounter, ycounter;                 //counters for each direction
  try
  {
    if (!fromprojection || !toprojection)      //check for projections
      throw ProjectorException(PROJECTOR_PROJECTION);
    
    
    pmesh = setupForwardPmesh();               //setup the forward mesh

    getExtents(pmesh);                         //try to get the extents
    
    setupOutput(outfile);                      //setup the output file
      
    if (pmesh)                                 //check for existing pmesh
    {
      delete pmesh;                             
      pmesh = setupReversePmesh();             //setup the reverse mesh
    }
    
    if (!(scanline = new (std::nothrow) unsigned char [newwidth*spp*(bps/8)]))
      throw std::bad_alloc();
            
    if (!cache)
    {
      if (!(inscanline = new (std::nothrow) 
            unsigned char[oldwidth*spp*(bps/8)]))
        throw std::bad_alloc();
    }

    //init the status progress
    if (progress)
    {
      std::strstream tempstream;

      tempstream << "Reprojecting " << newheight << " lines." << std::ends;
      tempstream.freeze(0);
      progress->init(tempstream.str(),
                      NULL,
                      "Done.",
                      newheight, 
                      29);
      progress->start();  //start the progress
    }

                                               //now start adding pixels
    for (ycounter =0; ycounter < newheight; ycounter++)
    {
      if (progress && !(ycounter % 29))     //check for output status func
        progress->update(ycounter);
      
      
      for (xcounter = 0; xcounter < newwidth; xcounter++)
      {   
                                               //get the new x and y pixel
        x = outRect.left + newscale.x * xcounter;
        y = outRect.top  - newscale.y * ycounter;
        
        
        if (pmesh)                             //now get the old pixel
          pmesh->projectPoint(x, y);
        else
        {
          toprojection->projectToGeo(x, y, y, x);
          fromprojection->projectFromGeo(y, x, x, y);
        }
        
        _x = static_cast<long int>((x - inRect.left) / 
                                   (oldscale.x) + 0.5);
        _y = static_cast<long int>((inRect.top - y) / 
                                   (oldscale.y) + 0.5);
        /*if (_x >= oldwidth)
          _x = _x - 1;
          else
          if (_x < 0)
          _x = 0;
          if (_y >= oldheight)
          _y = _y -1;
          else
          if (_y < 0)
          _y = 0;*/
        if ((_x >= oldwidth) || (_x < 0) || (_y >= oldheight) 
            || (_y < 0))
        {
          if (bps != 16) //8bits per sample
          {
            for (sppcounter = 0; sppcounter < spp; sppcounter++)
              reinterpret_cast<unsigned char *>(scanline)
                [xcounter*spp + sppcounter] = 0; //out of bounds pixel
          }
          else//16 bits per sample
          {
             for (sppcounter = 0; sppcounter < spp; sppcounter++)
              reinterpret_cast<uint16 *>(scanline)
                [xcounter*spp + sppcounter] = 0; //out of bounds pixel
          }
        }
        else
        {
          if (cache) //8bits per sample
          {
            cachescanline = cache->getRawScanline(_y); //get a pointer to cache
            for (sppcounter = 0; sppcounter < spp; sppcounter++)
              reinterpret_cast<unsigned char *>(scanline)
                [xcounter*spp + sppcounter ] = 
                cachescanline[_x*spp + sppcounter];       //copy pixels
          }
          else//no cache
          {
            if (bps == 16)
            {
              dynamic_cast<USGSImageLib::TIFFImageIFile*>(infile)
                ->getRawScanline(_y, inscanline);
              for (sppcounter = 0; sppcounter < spp; sppcounter++)
                reinterpret_cast<uint16*>(scanline)
                  [xcounter*spp + sppcounter] = 
                  reinterpret_cast<uint16 *>(inscanline)
                  [_x*spp + sppcounter];
            }
            else
            {
              infile->
                getRawScanline(_y, reinterpret_cast<unsigned char *>
                               (inscanline));
              for (sppcounter = 0; sppcounter < spp; sppcounter++)
                reinterpret_cast<unsigned char *>(scanline)
                  [xcounter*spp + sppcounter ] = 
                  reinterpret_cast<unsigned char*>(inscanline)
                  [_x*spp + sppcounter];       //copy pixels
            }
          }
            
        }
        
      }
      out->putRawScanline(ycounter, scanline);       //write out scanlines
    }
    
    //finsh the progress
    if (progress)
      progress->done();


    writer.removeImage(0);                           //flush the output file
    out = NULL;
    delete [] reinterpret_cast<unsigned char *>(scanline);  
    delete [] reinterpret_cast<unsigned char *>(inscanline);
    //delete the scanline
    delete pmesh;                                    //delete the pmesh
  }
  catch(ProjectorException & temp)
  {
    throw temp;                                      //just rethrow
  }
  catch(...)
  {
    delete [] reinterpret_cast<unsigned char *>(scanline);
                             //delete the scanline
    delete [] reinterpret_cast<unsigned char *>(inscanline);
    delete pmesh;                                    //delete the pmesh
    writer.removeImage(0);                           //flush output file
    throw ProjectorException(PROJECTOR_ERROR_UNKOWN);
  }      
}

//**********************************************************
void Projector::setupInput(std::string & ininfile) throw(ProjectorException)
{
  double scale[3] = {0};                       //used to get the scale
  double tp[6] = {0};                          //used to get the tiepoints
  short unsigned int tpnum;                    //stores the number of tiepoints
  float xscale = 0.0;                          //scale in x dir
  USGSImageLib::GeoTIFFImageIFile* ingeo;      //used to recast as a geotiff
  USGSImageLib::DOQImageIFile * indoq;         //used to recast as a doq
  // unsigned int nlines;                         //the number of lines to cache

  try
  {
    //check for prexisting input file
    if (infile)
    {
      delete infile;
      infile = NULL;
    }
    
    //try to open the file as a doq
    if (!(infile = new(std::nothrow)USGSImageLib::DOQImageIFile(ininfile)))
      throw std::bad_alloc();
    
    if (!infile->good())
    {
      delete infile;
      
      if (!(infile = new(std::nothrow)
            USGSImageLib::GeoTIFFImageIFile(ininfile, false)))
        throw std::bad_alloc();
      
      if (!infile->good())
      {
        throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
      }
      
      //check for preexisting fromprojection
      if (fromprojection)
      {
        reader.removeProjection(0);
      }
      
      //create the from projection
      fromprojection = reader.createProjection
        (dynamic_cast<USGSImageLib::GeoTIFFImageIFile*>(infile));
      
      if (!fromprojection)
        throw ProjectorException(PROJECTOR_ERROR_BADINPUT);

      
      //get the scale and such
      ingeo = dynamic_cast<USGSImageLib::GeoTIFFImageIFile *>(infile);
      ingeo->getPixelScale(scale);
      oldscale.x = scale[0];
      oldscale.y = scale[1];
      ingeo->getTiePoints(tp, tpnum);
      inRect.left = tp[3];
      inRect.top  = tp[4];
      
    }
    else
    {
      fromprojection = reader.createProjection
        (dynamic_cast<USGSImageLib::DOQImageIFile*>(infile));
      
      if (!fromprojection)
      {
        throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
      }
      
      indoq = dynamic_cast<USGSImageLib::DOQImageIFile * > (infile);
      
      indoq->getHorizontalResolution(xscale);
      oldscale.x = xscale;
      oldscale.y = xscale;
      indoq->getXOrigin(inRect.left);
      indoq->getYOrigin(inRect.top);
      
    }
    
    //get some image metrics
    infile->getHeight(oldheight);
    infile->getWidth(oldwidth);
    infile->getSamplesPerPixel(spp);
    infile->getBitsPerSample(bps);
    infile->getPhotometric(photo);    
    
    //check for preexisting cache
    if (cache)
      delete cache;
    cache = NULL;
    
  
    if ((bps == 8) && cachesize) //check for 16 bit tiffs that can't use cache
    {
      /*
      //figure out the number lines to cache
      nlines = static_cast<unsigned int>((static_cast<double>(cachesize))/
        (static_cast<double>(bps*spp*oldwidth)/8388608.0));
      //check to see if the cache is larger than image size
      if (nlines > static_cast<unsigned int>(oldheight))
        nlines = oldheight;


      //check to see if the cache is now zero (really big image!) 
      if (nlines != 0)
      {
        //create the cache for the file
        if (!(cache = new(std::nothrow) USGSImageLib::CacheManager
              (infile, nlines/2, true)))
          throw std::bad_alloc();
        if (!cache->good())
        {
          //can't open the cache so delete it and set default cache
          delete cache;
          cache = NULL; //set it to NULL
        }
        }*/
      if (!(cache = new (std::nothrow) USGSImageLib::LRUCacheManager
            (infile, 1000)))
        throw std::bad_alloc();
    } 
      
    //set the dimensions
    inRect.right = inRect.left + oldscale.x * oldwidth;
    inRect.bottom = inRect.top - oldscale.y * oldheight;

  }
  catch(...)
  {
    delete infile;                  //delete stuff and rethrow        
    delete cache;
    throw ProjectorException(PROJECTOR_UNABLE_INPUT_SETUP);
  }
}

//***********************************************************************
void Projector::setupOutput(std::string & inoutfile) throw(ProjectorException)
{
  double tp[6] = {0};
  double res[3] = {0};
  
  try
  {
    //perform some checks
    if (!newwidth || !newheight || !infile)
      throw ProjectorException();
    
    if (out)                 //check for preexisting output file
    {
      writer.removeImage(0); //keep only one geotiff at a time
      out = NULL;
    }
    //Create the output file
    tp[3] = outRect.left;
    tp[4] = outRect.top;
    res[0] = newscale.x;
    res[1] = newscale.y;
    res[2] = 0;
    out = writer.create(toprojection, outfile, newwidth, 
                        newheight, photo,tp, res);
    

    if (!out)
      throw ProjectionException();

    
    //check for palette
    if (photo == PHOTO_PALETTE)
    {
      USGSImageLib::RGBPalette temppal;
      infile->getPalette(&temppal);
      out->setPalette(&temppal);
    }
    
    //see if we want to compress output
    if (packbits)
    {
      dynamic_cast<USGSImageLib::TIFFImageOFile*>(out)->setCompression
        (COMPRESSION_PACKBITS);
    }
    out->setSamplesPerPixel(spp);
    out->setBitsPerSample(bps);

  }
  catch(...)
  {
    throw ProjectorException(PROJECTOR_UNABLE_OUTPUT_SETUP);
  }
}



//**********************************************************************
PmeshLib::ProjectionMesh * Projector::setupForwardPmesh() throw()
{
  PmeshLib::ProjectionMesh * ret = NULL;          //return pmesh
  
  try
  {
    if (pmeshname != 0)
    {
      if(!(ret = new (std::nothrow) PmeshLib::ProjectionMesh))
        throw std::bad_alloc();
      
      //setup the foward mesh
      ret->setSourceMeshBounds(inRect.left, inRect.bottom, 
                               inRect.right, inRect.top);
      ret->setMeshSize(pmeshsize,pmeshsize); 
      ret->setInterpolator(pmeshname);
      //now project all the points onto the mesh
      ret->calculateMesh((*fromprojection), (*toprojection));
    }
    return ret;
  }
  catch(...)
  {
    delete ret;
    ret = NULL;
    return ret;
  }
}

//*********************************************************************
PmeshLib::ProjectionMesh * Projector::setupReversePmesh() throw()
{
  PmeshLib::ProjectionMesh * ret = NULL;          //return pmesh
  
  try
  {
    if (pmeshname != 0)
    {
      if(!(ret = new (std::nothrow) PmeshLib::ProjectionMesh))
        throw std::bad_alloc();
      
      //setup the reverse mesh
      ret->setSourceMeshBounds(outRect.left, outRect.bottom, 
                               outRect.right, outRect.top);
      ret->setMeshSize(pmeshsize,pmeshsize); 
      ret->setInterpolator(pmeshname);
      //now project all the points onto the mesh
      ret->calculateMesh((*toprojection), (*fromprojection));
    }
    return ret;
    
  }
  catch(...)
  {
    delete ret;
    ret = NULL;
    return ret;
  }
}



//**********************************************************************
void Projector::getExtents(PmeshLib::ProjectionMesh * pmesh) throw(ProjectorException)
{
  std::vector<double> xarr, yarr;     //vectors for sorting
  long int xcounter(0), ycounter(0);  //counters for loops
  double tempx(0), tempy(0);          //temp-o-vars

  try
  {
    //do some checks
    if (!infile)
      throw ProjectorException(PROJECTOR_INPUT_NOT_SET);
    
    if (!toprojection || !fromprojection)
      throw ProjectorException(PROJECTOR_PROJECTION);
    
    //make room in the vector for the edges
    xarr.resize(2*(oldwidth+oldheight));
    yarr.resize(2*(oldheight + oldwidth));
    
    //now we should be good to go
    for (xcounter = 0; xcounter < oldwidth; xcounter++)
    {
      //reproject the top line
      tempx = inRect.left + oldscale.x*xcounter;
      tempy = inRect.top;
      if (pmesh)
        pmesh->projectPoint(tempx, tempy);
      else
      {
        fromprojection->projectToGeo(tempx, tempy, tempy, tempx);
        toprojection->projectFromGeo(tempy, tempx, tempx, tempy);
      }
      xarr[xcounter] = tempx;
      yarr[xcounter] = tempy;
      
      //reproject the bottom line
      tempx = inRect.left + oldscale.x*xcounter;
      tempy = inRect.bottom;
      if (pmesh)
        pmesh->projectPoint(tempx, tempy);
      else
      {
        fromprojection->projectToGeo(tempx, tempy, tempy, tempx);
        toprojection->projectFromGeo(tempy, tempx, tempx, tempy);
      }
      xarr[xcounter + oldwidth] = tempx;
      yarr[xcounter + oldwidth] = tempy;
    }
    
    //now do the height
    for (ycounter = 0; ycounter < oldheight; ycounter++)
    {
      //reproject the left line
      tempy = inRect.top - oldscale.y * ycounter;
      tempx = inRect.left;
      if (pmesh)
        pmesh->projectPoint(tempx, tempy);
      else
      {
        fromprojection->projectToGeo(tempx, tempy, tempy, tempx);
        toprojection->projectFromGeo(tempy, tempx, tempx, tempy);
      }
      xarr[ycounter + 2*oldwidth] = tempx;
      yarr[ycounter + 2*oldwidth] = tempy;
      
      //reproject the right line
      tempy = inRect.top - oldscale.y * ycounter;
      tempx = inRect.right;
      if (pmesh)
        pmesh->projectPoint(tempx, tempy);
      else
      {
        fromprojection->projectToGeo(tempx, tempy, tempy, tempx);
        toprojection->projectFromGeo(tempy, tempx, tempx, tempy);
      }
      xarr[ycounter + 2*oldwidth + oldheight] = tempx;
      yarr[ycounter + 2*oldwidth + oldheight] = tempy;
    }
    
    
    //get the min and max of the coordinates
    getMinMax(xarr, outRect.left, outRect.right);
    getMinMax(yarr, outRect.bottom, outRect.top);
    
    
    if (samescale)
      newscale = getSameScale(oldscale, fromprojection, toprojection);
    
    //get the new scale
    if (!newscale.x || !newscale.y)
      newscale = GetConvertedScale(oldwidth,oldheight,outRect);
    
    //get the new pixel width and height
    newwidth  = static_cast<long int>
      ((outRect.right - outRect.left)/(newscale.x) + 0.5);
    newheight = static_cast<long int>
      ((outRect.top - outRect.bottom)/(newscale.y) + 0.5);

  }
  catch(ProjectorException & temp)
  {
    throw temp;            //just rethrow
  }
  catch(...)
  {
    throw ProjectorException(PROJECTOR_UNABLE_OUTPUT_SETUP);
  }
}

#endif




