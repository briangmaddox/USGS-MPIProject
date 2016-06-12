//This is the sticher class that 
//stitches the image together that does not use
//Ace but instead uses pthread.

#ifndef STITCHER_H_
#define STITCHER_H_

#include "ImageLib/ImageOFile.h"
#include <ace/OS.h>
#include <ace/Synch.h>
#include <queue>
#include "StitcherNode.h"




//This is the function that starts a thread that the lives in
//the class
void * start_func(void * class_instance);


//The stitcher class
class Stitcher
{
 public:
  
  /**
   * This is the constructor for the stitcher class 
   * that takes a pointer to the output file.
   * Immediately starts the interal thread
   * DON'T MESS WITH THE FILE AFTER CREATING
   * THE THREAD UNTIL THE THREAD IS TERMINATED
   * (IMAGELIB is NOT thead safe)
   **/
  Stitcher(USGSImageLib::ImageOFile * inout);

  
  /**
   * Destructor: It is up to the user to close the output file
   **/
  virtual ~Stitcher();

  
  /**
   * Add function adds a stitcher node to the stitchers interal
   * queue to be processed.
   * (The stitcher will delete the node when it processes it)
   **/
  void add(StitcherNode * temp) throw();
  
  /**
   * Wait function will wait (block) until the stitcher completely 
   * finishes writing the image.
   **/
  void wait() throw();

  /**
   * run function is where the the thread actually runs and should
   * not ever be called by any outside thread.
   **/
  void run() throw();
  
 private:

  USGSImageLib::ImageOFile * out;
  ACE_Thread_Mutex donemutex;              //the mutex to see if we are done
  ACE_Thread_Mutex workmutex;              //the mutex for the working queue
  ACE_Condition<ACE_Thread_Mutex> workcond;//the condition if there is 
                                           //something in the queue
  bool done;                               //the done tag tells the stitcher
                                           //that it is done
  ACE_Condition<ACE_Thread_Mutex> waitcond;//the condition that the waiter
                                           //waits on
  std::queue<StitcherNode *> workqueue;    //the working que

};


#endif
