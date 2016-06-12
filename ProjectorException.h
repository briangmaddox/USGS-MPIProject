/***************************************************************************
                          ProjectorException.h  -  description
                             -------------------
    begin                : Thu Oct 12 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/


#ifndef PROJECTOREXCEPTION_H
#define PROJECTOREXCEPTION_H

#include <new>
#include <string>


//defines for errors

#define PROJECTOR_ERROR_BADINPUT       0
#define PROJECTOR_INVALID_PROJECTION   1
#define PROJECTOR_INVALID_UNIT         2
#define PROJECTOR_INVALID_DATUM        3
#define PROJECTOR_UNABLE_INPUT_SETUP   4
#define PROJECTOR_UNABLE_OUTPUT_SETUP  5
#define PROJECTOR_INPUT_NOT_SET        6
#define PROJECTOR_PROJECTION           7
#define PROJECTOR_ERROR_UNKOWN         255


//ProjectorException Class
class ProjectorException
{
public:
  //Constructors and Destructor
  ProjectorException();
  ProjectorException(const short unsigned int & inexception);
  ~ProjectorException();

  void setException(const short unsigned int & inProjection) throw();
  short unsigned int getException() const throw();
  std::string getExceptionMessage() const throw();

protected:
  short unsigned int myexception;
};

#endif
