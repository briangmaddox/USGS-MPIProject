/**
 * Impelementation file for the CLineProgress class
 **/


#ifndef CLINEPROGRESS_CPP
#define CLINEPROGRESS_CPP

#include "CLineProgress.h"


//************************************************
CLineProgress::CLineProgress() : BaseProgress()
{}

//************************************************
CLineProgress::~CLineProgress()
{
}

//************************************************
bool CLineProgress::init(const char * instartmesg,
                         const char * intagmesg,
                         const char * indonemesg,
                         unsigned int innumthings, 
                         unsigned int instepsize) throw()
{
  
  if (instartmesg)
  {
    startmesg = instartmesg;
  }

  if (intagmesg)
  {
    tagmesg = intagmesg;
  }

  if (indonemesg)
  {
    donemesg = indonemesg;
  }
  

  numthings = innumthings;
  stepsize = instepsize; 
  current = 0;
  
  if (!numthings || !stepsize)
  {
    //makes sense to give a output error here
    std::cout << "CLineProgress: Number of steps or stepsize invalid" 
              << std::endl;
    valid = false;                          //class is not valid
  }
  else
  {
    valid = true;                           //class is valid
  }

  return valid;                             //return that it is valid
}

//************************************************************
void CLineProgress::start() throw()
{
  //display the startup message
  if (startmesg.size() > 0 && valid)
  {
    std::cout << startmesg.c_str() << std::endl;
  }
}

//************************************************************
void CLineProgress::update() throw()
{
  if (valid)
  {
    current = current + stepsize;      //increment the progress
    
    spinner.update(static_cast<int>(current));   //update the spinner
  }
}

//*************************************************************  
void CLineProgress::update(int incurrent) throw()
{
  if (valid)
  {
    current = incurrent;
    spinner.update(static_cast<int>(current));
  }
}

//**************************************************************
void CLineProgress::done()  throw()
{
  if (valid)
  {
    valid = false;                     //no longer valid
    current = 0;                       //current is back to zero
    if (donemesg.size())
      spinner.done(donemesg.c_str());  //output the done mesg if it exists
  }
}


#endif


