/***************************************************************************
                          ProjectorException.cpp  -  description
                             -------------------
    begin                : Thu Oct 12 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/

#ifndef PROJECTOREXCEPTION_CPP
#define PROJECTOREXCEPTION_CPP


#include "ProjectorException.h"


//**********************************************************************
ProjectorException::ProjectorException() : myexception(PROJECTOR_ERROR_UNKOWN)
{}
  
//**********************************************************************
ProjectorException::
ProjectorException(const short unsigned int & inexception)
{
  myexception = inexception;
}

//**********************************************************************
ProjectorException::~ProjectorException()
{}

//*********************************************************************
void ProjectorException::
setException(const short unsigned int & inexception) throw()
{
  myexception = inexception;
}

//**********************************************************************
short unsigned int ProjectorException::getException() const throw()
{
  return myexception;
}

//**********************************************************************
std::string ProjectorException::getExceptionMessage() const throw()
{
  switch(myexception)
    {
    case 0: 
      return std::string("Projector Error: Bad Input");
    case 1:
      return std::string("Projector Error: Invalid Projection");
    case 2:
      return std::string("Projector Error: Invalid Unit");
    case 3:
      return std::string("Projector Error: Invalid Datum");
    case 4:
      return std::string("Projector Error: Unable To Setup Input File");
    case 5:
      return std::string("Projector Error: Unable To Setup Output File");
    case 6:
      return std::string("Projector Error: Input File Not Set");
    case 7:
      return std::string("Projector Error: Input or Output Projection is Invalid");
    default:
      return std::string("Projector Error: Unkown Error");
    }
}

#endif
