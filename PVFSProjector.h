/**
 *This is a extension onto the pvm project class
 *that is designed to take advantage of pvfs partitioning
 *By Chris Bilderback
 **/

#ifndef PVFSPROJECTOR_H
#define PVFSPROJECTOR_H

#include "MpiProjector.h"
#include <vector>
#include <hash_map>


//PVFSProject class partitions slaves to
//work on only a specific portion of the file
//which can reduce network traffic as less slaves
//are talking to the master to get their data.
class PVFSProjector : public MpiProjector
{
 public:

  //Constructor and destructor
  PVFSProjector();
  ~PVFSProjector();

  //This function indicates the number of pvfs partition that
  //is being worked with. (Default is 2)
  void setPartitionNumber(const unsigned int & inPart) throw();
  unsigned int getPartitionNumber() const throw();

  //main function which runs the projection
  virtual void 
    project(BaseProgress * progress = NULL) 
    throw(ProjectorException);

 protected:
  
  //This function actually projects the file
  void projectPVFS(BaseProgress * progress = NULL)
    throw(ProjectorException);

  //This function tells the slave to terminate
  long int PVFSProjector::terminateSlave(MPI_Status status,
                          Stitcher * mystitch,  
                          unsigned char * buffer,
                          long int buffersize) throw(ProjectorException);



  std::vector<long int> mcounters;    //membership counters
  std::vector<long int> mstop;        //where to stop each membership
  std::hash_map<int, unsigned int> 
    membership;                       //nodal membership map.

};

#endif
