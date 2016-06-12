#ifndef PROJECTORSTITCHER_CPP
#define PROJECTORSTITCHER_CPP



#include "ProjectorStitcher.h"


//**********************************************************************
ProjectorStitcher::ProjectorStitcher() : ACE_Task<ACE_MT_SYNCH>(), 
scanlinecond(listmutex), out(NULL), numberofscanlines(0), 
scanlinelist(NULL),current(0)  
{
}

//********************************************************************
ProjectorStitcher::~ProjectorStitcher()
{
  int counter;

  if (scanlinelist)
    {
      for (counter = 0; counter < numberofscanlines; counter++)
        delete [] scanlinelist[counter];
      delete [] scanlinelist;
    }
}

  
//***********************************************************************
int ProjectorStitcher::open(void*)
{
  //do a check for setup
  if (!scanlinelist || !out)
    return 0;

  activate(THR_NEW_LWP,1);
  return 1;
}

//**********************************************************************
int ProjectorStitcher::close(u_long flags)
{
  return 0;
}

//*********************************************************************
int ProjectorStitcher::svc()
{
  
  checkScanlines();                 //check the scanline cache
  done.set(true);                   //set the future item
  return 0;                         //done
}
      
//**********************************************************************
void ProjectorStitcher::
setOutputFile(USGSImageLib::ImageOFile * inout) throw()
{
  out = inout;
}
 
 
//*********************************************************************
void ProjectorStitcher::
setupScanlineList(const long int & inNumberOfScanlines) throw()
{
  int counter;                     //counter for loop

  try
    {
      if (scanlinelist)            //preexisting scanline cache
        {
          for (counter = 0; counter < numberofscanlines; counter++)
            delete [] scanlinelist[counter];
          delete [] scanlinelist;
          scanlinelist = NULL;
        }
      numberofscanlines = inNumberOfScanlines;

	                                //create the scanline list
      if (!(scanlinelist = new (std::nothrow) unsigned char 
            *[numberofscanlines]))
        throw std::bad_alloc();

      for (counter = 0; counter < numberofscanlines; counter++)
        scanlinelist[counter] = NULL;//init the scanline list

    }
  catch(...)
    {
      if (scanlinelist)     //clean up memory
        {
          for (counter = 0; counter < numberofscanlines; counter++)
            delete [] scanlinelist[counter];
          delete [] scanlinelist;
          scanlinelist = NULL;
        }
    }
}
      
//*******************************************************************
void ProjectorStitcher::insertScanline(const long int & inScanlineNumber, 
                      unsigned char * inBits) throw()
{
  //Since the stitcher deletes memory, just copy the pointer
  scanlinelist[inScanlineNumber] = inBits;
}

//*********************************************************************
float ProjectorStitcher::getPercentFull() const throw()
{
  float pfull = 0;              //return value
  int counter;                  //counter for loop
  try
    {
	                            //check for fullness of scanline list
      for (counter = 0; counter < numberofscanlines; counter++)
        if (scanlinelist[counter] != NULL)
          pfull++;

                                 //divide to get the averge
	  pfull = pfull/numberofscanlines; 
      return pfull;
    }
  catch(...)
    {
      return -1;                 //error
    }
}

//************************************************************************
ACE_Future<bool> ProjectorStitcher::getFuture() const throw()
{
 return done;                    //return future item
}

//*************************************************************************
void ProjectorStitcher::checkScanlines() throw()
{
  long int height;
  out->getHeight(height);        //get the height
 
  while (current < height)       //proccess scanlines
    {

      if (scanlinelist[current] != NULL)
      {                          //write out the current scanline
        
        out->putRawScanline(current, scanlinelist[current]);
        delete [] scanlinelist[current];
        scanlinelist[current] = NULL;
        current++;             //increment the current line
        
       }
      else
       {
         //std::cout << "Sticher Sleeping " << std::endl;
         ACE_OS::sleep(10);       //sleep then check again
       }
    }
}

#endif







