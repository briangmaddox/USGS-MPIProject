/**
 * This forms the base class of that displays the progress
 * of a work process.
 * By Chris Bilderback
 **/


#ifndef BASEPROGRESS_H
#define BASEPROGRESS_H

#include <string>

/**
 * BaseProgress is a abstract class from
 * displaying progress.
 **/
class BaseProgress
{
public:
  
  /** 
   * Main constructor for the class
   **/
  BaseProgress();

  /**
   * Main destructor for the class
   **/
  virtual ~BaseProgress();

  
  /**
   * init sets up the progress display
   * innumthings = the number of things to do,
   * instepsize is the step size that the progress will be updated in.
   * instartmesg is the message to be displayed at start.
   * intagmesg  is a message to be displayed with each step
   * indonemesg is a message to be displayed when done
   **/
  virtual bool init(const char * instartmesg,
                    const char * intagmesg,
                    const char * indonemesg,
                    unsigned int innumthings, 
                    unsigned int instepsize) throw()=0;
  /**
   * start puts the progress in the start state nomatter where it was before
   **/
  virtual void start() throw() = 0;


  /**
   * update updates the progress by one step.
   **/
  virtual void update() throw()=0;
  
  /**
   * overloaded update allows you to set the current step
   **/
  virtual void update(int incurrent) throw()=0;

  /**
   * done puts the progress in the done state no matter where it was before.
   **/
  virtual void done()  throw()=0;

  /**
   * get valid returns whether the progress is valid
   **/
  bool getValid() const throw();

protected:
  
  std::string startmesg, tagmesg, donemesg;  //the messages
  unsigned int numthings, stepsize;          //the number of things and the
                                             //step size
  unsigned current;                          //the current state
  bool valid;                                //is the progress valid
};


#endif
