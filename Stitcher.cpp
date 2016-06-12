/**
 *This is the implementation file for the stitcher node
 **/


#ifndef STITCHER_CPP
#define STITCHER_CPP


#include "Stitcher.h"


//*************************************************************
void * start_func(void * class_instance)
{
  //run the class
  reinterpret_cast<Stitcher*>(class_instance)->run();
  return 0;

}


//*************************************************************
Stitcher::Stitcher(USGSImageLib::ImageOFile * inout)
  : donemutex(), workmutex(), workcond(workmutex), done(false),
    waitcond(donemutex)
{
  
  //copy the output file
  out = inout;
  
  //start the thread
  ACE_Thread::spawn_n (1, (ACE_THR_FUNC)start_func, 
                       reinterpret_cast<void *>(this));

}
  
  
//*************************************************************
Stitcher::~Stitcher()
{
  //This is really just for premature termination by the calling thread
  StitcherNode * tempnode(0);
  
  //lock the que
  workmutex.acquire();
  
  //check the queue
  while(workqueue.size())
  {
    tempnode = workqueue.front();
    workqueue.pop();
    delete tempnode;
  }
  
  //create the termination node
  if (!(tempnode = new (std::nothrow) StitcherNode(-1, -1, 0)))
    throw std::bad_alloc();
  
  //put it in the queue
  workqueue.push(tempnode);

  //signal the thread
  workcond.signal();
  workmutex.release();
  
  //lock the done mutex and wait
  donemutex.acquire();
  while(!done)
    waitcond.wait();
  donemutex.release();
  

  
}

  
//***********************************************************
void Stitcher::add(StitcherNode * temp) throw()
{
  //aquire the workmutex
  workmutex.acquire();
  //put the work in the que
  workqueue.push(temp);
  //signal the thread
  workcond.signal();
  workmutex.release();
}

//**********************************************************
void Stitcher::wait() throw()
{
  //try to aquire the done mutex
  donemutex.acquire();
  while(!done)
    waitcond.wait();
  donemutex.release();
}

//**************************************************************
void Stitcher::run() throw()
{
  
  long int counter(0), newheight(0), newwidth(0);
  int spp(0), bps(0);
  StitcherNode * temp(0);
  int start(0), stop(0);
  unsigned char * data(0);
  bool ldone(false);
  //check the output image
  if (!out)
  {
    donemutex.acquire();
    done = true;
    waitcond.signal();
    donemutex.release();
    return; //exit
  }

  //get the image info
  out->getHeight(newheight);
  out->getWidth(newwidth);
  out->getSamplesPerPixel(spp);
  out->getBitsPerSample(bps);
  //loop
  while(!ldone)
  {
    
    //check the queue
    workmutex.acquire();
    if (!workqueue.size())
    {
      workcond.wait();
      workmutex.release();
    }
    else
    {
      temp = workqueue.front();
      workqueue.pop();
      workmutex.release();
      start = temp->getStart();
      stop = temp->getEnd();
      data = temp->getData();
      
      //check early termination
      if (start < 0)
      {
        ldone = true;
        delete temp;
      }
      else
      {
      
        //write the scanlines to disk
        for (counter = start; counter <= stop;  ++counter)
        {
          //write it
          out->putRawScanline(counter, 
                              &(data[newwidth*spp*(counter-start)]));
        } 
        
        //delete the scanline chunk
        delete temp;
      

        //check termination
        if (stop == newheight-1)
        {
          //we are done so set the done flag
          ldone = true;
        }
      
      }      
    }
  }


  //termination
  donemutex.acquire();
  done = true;
  waitcond.signal();
  donemutex.release();

}


#endif







