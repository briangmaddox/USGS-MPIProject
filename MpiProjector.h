/**
 * This forms the core of the master mpi projection program.
 * By Chris Bilderback.
 **/


#ifndef MPIPROJECTOR_H_
#define MPIPROJECTOR_H_

#include "Projector.h"
#include "MessageTags.h"
#include <mpi.h>
#include <queue>
#include "Stitcher.h"

//The master pvm projector
class MpiProjector : public Projector
{
 public:
  //Constructor and Destructor
  MpiProjector();
  virtual ~MpiProjector();
  
  //allows user to set the number of slaves spawned
  void setNumberOfSlaves(const int & innumslaves) throw();
  int getNumberOfSlaves() const throw();
  
  //allows user to set the chunksize which is by default 1
  //added 3/5/01 by Chris Bilderback
  void setChunkSize(const int & inchunksize) throw();
  int getChunkSize() const throw();

  //allows the user to set the sequence of chunksizes that the projector
  //sends chunks to the slaves 
  void setSequence(const int * insequence,
                   const int & insequencesize) throw(std::bad_alloc);
  
  //allows the user to set whether to use a random sequence of 
  //chunksizes between minchunk and maxchunk
  void setRandomSequence(const int & inminchunk,
                         const int & inmaxchunk) throw();
  
  //This resets the sequencing method to the default
  void resetSequencing() throw();


  //setEvenChunks divides up the scanlines in the file
  //addd 3/5/01 by Chris Bilderback
  //This automatically sets slavelocal
  void setEvenChunks(bool inevenchunks) throw();
  bool getEvenChunks() const throw();

  //allows the user to say whether the slave should store
  //all its data locally and then send back its work at end
  //FOR RIGHT NOW: the slave assumes that if you want to store
  // the work locally it will write it's chunk to disk no
  // matter the work size. DEFAULT is true
  void setSlaveStoreLocal(const bool & inslavelocal) throw();
  bool getSlaveStoreLocal() const throw();


  //this allows the user to set were the path to store the 
  //slave locally is stored.
  //This is initally set to "./" indicating that the working
  //director should be used (this may not be what we want)
  void setSlaveLocalDir(const std::string & inslavelocalpath) throw();
  std::string getSlaveStoreLocalDir() const throw();

  //This tells the master to use the threaded stitcher when writting 
  //output
  void setStitcher(bool institcher) throw();
  bool getStitcher() const throw();

  //overloaded to save the filename
  virtual void setInputFile(std::string & ininfile) throw(ProjectorException);
  
  //main function which runs the projection
  virtual void 
    project(BaseProgress * progress = NULL) 
    throw(ProjectorException);
  
protected:
  //This function sends the stitcher a chunk
  long int sendStitcher(Stitcher * mystitch, unsigned char * buffer, 
                        long int buffersize) throw();
  
  
  //If the slaves were set to store the data locally this function
  //will handle that projection case
  bool projectslavelocal(BaseProgress * progress = NULL)
    throw(ProjectorException);
  
  //If the slaves don't store the data locally
  bool projectnoslavelocal(BaseProgress * progress = NULL)
    throw(ProjectorException);


  //sendSlaveSetup function does all the pvm packing to send to the slave
  //for setup
  bool sendSlaveSetup(int rank) throw();

  //unpackScanline unpacks a scanline from the slave and writes it
  long int unpackScanline(unsigned char * buffer, 
                          long int buffersize) throw();

  int numofslaves;                   //the number of slaves
  bool evenchunks;                   //are we using even chunks
  bool slavelocal;                   //should the slaves store and then send?
  int sequencemethod;                //the method of chunksize sequence
                                     //0 is uniform, 1 is sequence, 2 random
  int minchunk, maxchunk;            //used for the random sequencing and to
                                     //tell the slave what the maxchunk size is
  int * sequence;
  int sequencesize;

  std::string inputfilename;
  std::string slavelocalpath;        //the base directory were the slaves
                                     //should store there files
  bool stitcher;                     //wheather the master should use the
                                     //sticher

};

#endif
  



