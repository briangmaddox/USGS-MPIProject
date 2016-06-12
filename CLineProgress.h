/** ClineProgress is a command line progress that 
 *  uses SpinnerThingy to output progress to a command line
 **/


#ifndef CLINEPROGRESS_H
#define CLINEPROGRESS_H

#include "MiscUtils/SpinnerThingy.h"
#include "BaseProgress.h"


/**
 * CLineProgress class is a command line progress class
 * that is really just a wrapper around spinnerthing.
 **/
class CLineProgress : public BaseProgress
{
 public:
  
  /**
   * Main constructor for the class.
   **/
  CLineProgress();
  
  /**
   * Destructor for the class
   **/
  virtual ~CLineProgress();


  /**
   * init sets up the progress display
   * innumthings = the number of things to do,
   * instepsize is the step size that the progress will be updated in.
   * instartmesg is the message to be displayed at start.
   * intagmesg  is a message to be displayed with each step
   * indonemesg is a message to be displayed when done
   **/
  bool init(const char * instartmesg,
            const char * intagmesg,
            const char * indonemesg,
            unsigned int innumthings, 
            unsigned int instepsize) throw();
  /**
   * start puts the progress in the start state nomatter where it was before
   **/
  void start() throw();


  /**
   * update updates the progress by one step.
   **/
  void update() throw();
  
  /**
   * overloaded update allows you to set the current step
   **/
  void update(int incurrent) throw();

  /**
   * done puts the progress in the done state no matter where it was before.
   **/
  void done()  throw();

private:
  
  MiscUtils::SpinnerThingy  spinner;   //the spinner that does the update
};

#endif
