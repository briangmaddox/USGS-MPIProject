#ifndef PROJECTORSTITCHER_H
#define PROJECTORSTITCHER_H


#include "ace/OS.h"
#include "ace/Task.h"
#include "ace/Synch.h"
#include "ace/Future.h"
#include "ImageLib/ImageOFile.h"

class ProjectorStitcher : public ACE_Task<ACE_MT_SYNCH>
{
  
public:
  //constructor and destructor
  ProjectorStitcher();
  virtual ~ProjectorStitcher();

  
  //ace active object overloads
  virtual int open(void*);
  virtual int close(u_long flags = 0);
  virtual int svc();


  //setOutput file sets the output file being written to
  void setOutputFile(USGSImageLib::ImageOFile * inout) throw();
  
  //setupScanlineList sets up the scanline ptr array
  void setupScanlineList(const long int & inNumberOfScanlines) throw();

  //insertScanline puts a scanline in the list,
  //It locks the mutex and signals the condition variable when done
  void insertScanline(const long int & inScanlineNumber, 
                      unsigned char * inBits) throw();

  //getPercentFull returns the percent full of the list
  float getPercentFull() const throw();

  //getFuture lets the main thread know that the stitcher is done
  ACE_Future<bool> getFuture() const throw();

protected:
  
  //checkScanlines checks the list for the current scanline and writes
  //it out
  void checkScanlines() throw();

  ACE_Future<bool> done;              //is the stitcher done?
  ACE_Thread_Mutex listmutex;         //the scanline list mutex
  ACE_Condition<ACE_Thread_Mutex> scanlinecond;
  USGSImageLib::ImageOFile * out;
  long int numberofscanlines;
  unsigned char ** scanlinelist;
  long int current;
};


#endif
  
  




