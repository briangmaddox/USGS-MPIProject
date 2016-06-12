#ifndef MPIPROJECTORSLAVE_H
#define MPIPROJECTORSLAVE_H



#include "Projector.h"
#include "MessageTags.h"
#include <mpi.h>
#include <queue>
#include <fstream>
#include <strstream>

//This is the node that descibes the file when the slave stores
//a chunk locally
class ChunkFile
{
 public:
  //Default constructor for the class.
  ChunkFile() {}
  //Secondary constructor for the class
  ChunkFile(const std::string & inpath,
            const int & instarty,
            const int & inendy) : path(inpath), starty(instarty),
                                  endy(inendy) {}


  std::string path;   //the path the local file
  long int starty;    //starting scanline
  long int endy;      //the endy
};

 


//Slave pvm projector
class MpiProjectorSlave : public Projector
{
 public:
  //Constructor and Destructor
  MpiProjectorSlave();
  virtual ~MpiProjectorSlave();

  //connect function attempts to connect to the master 
  //and runs the projection
  bool connect() throw();

 protected:
  //unpackSetup function unpacks setup info from the master
  void unpackSetup() throw();

  //storelocal function handles when the master tells the slave
  //to store its information locally.
  bool storelocal() throw();
  
  
  bool slavelocal;                 //default is false
  int mastertid, mytid;            //pvm name
  unsigned int maxchunk;           //maximum chunksize
  std::string basepath;            //the path to the local file directory
 

};

#endif
