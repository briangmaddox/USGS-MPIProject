/**
 * Implementation file for the BaseProgress class
 **/

#ifndef BASEPROGRESS_CPP
#define BASEPROGRESS_CPP

#include "BaseProgress.h"


//*******************************************************
BaseProgress::BaseProgress() :  numthings(0), stepsize(0),
                                current(0), valid(false)
                                
{
}

//********************************************************
BaseProgress::~BaseProgress()
{
}


//********************************************************
bool BaseProgress::getValid() const throw()
{
  return valid;
}

#endif

