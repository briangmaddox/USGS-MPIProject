#ifndef MPIPROJECTOR_CPP_
#define MPIPROJECTOR_CPP_


#include "MpiProjector.h"
#include <stdlib.h>
#include <time.h>
#include <cmath>

//*******************************************************************
MpiProjector::MpiProjector() : Projector(),
                               numofslaves(0),
                               evenchunks(false),slavelocal(false), 
                               sequencemethod(0),
                               minchunk(0), maxchunk(0),
                               sequence(0), sequencesize(0),
                               slavelocalpath("./"), stitcher(false)
{}

//*******************************************************************
MpiProjector::~MpiProjector()
{
  delete [] sequence;
}

//*******************************************************************
void MpiProjector::setNumberOfSlaves(const int & innumslaves) throw()
{
  
  numofslaves = innumslaves;                //set the number of slaves

}

//**********************************************************************
int MpiProjector::getNumberOfSlaves() const throw()
{
  return numofslaves;
}


//**********************************************************************
void MpiProjector::setChunkSize(const int & inchunksize) throw()
{
  resetSequencing();

  if (inchunksize == 0)                //check to see if it is zero
    maxchunk = 1;
  else
    maxchunk = inchunksize;           //set the chunksize

  
}

//**********************************************************************
int MpiProjector::getChunkSize() const throw()
{
  return maxchunk;
}

//**********************************************************************
void MpiProjector::setEvenChunks(bool inevenchunks) throw()
{
  evenchunks = inevenchunks;
}

//**********************************************************************
bool MpiProjector::getEvenChunks() const throw()
{
  return evenchunks;
}

//***********************************************************************
void MpiProjector::setSlaveStoreLocal(const bool & inslavelocal) throw()
{
  slavelocal = inslavelocal;   //set the slave local
}

//***********************************************************************
bool MpiProjector::getSlaveStoreLocal() const throw()
{
  return slavelocal;
}

  
//*********************************************************************
void MpiProjector::project(BaseProgress * progress) 
    throw(ProjectorException)
{
  PmeshLib::ProjectionMesh * pmesh(0);           //projection mesh
  int rank(0);                                   //my rank in mpi
 
  try
  {
    if (!numofslaves)                            //must have a slave pvmproject
    {
      Projector::project(progress);
      return;
    }
    
    if (!fromprojection || !toprojection)        //check for projection
    {
      throw ProjectorException(PROJECTOR_ERROR_UNKOWN);
    }
    
    pmesh = setupForwardPmesh();                 //try setup the forward
                                                 //pmesh
    
    getExtents(pmesh);                           //get the extents
    
    
    if(cache)                                    //delete the cache
    {
      delete cache;
      cache = NULL;
    }
      
    setupOutput(outfile);                        //create the output file
        
    
    if (pmesh)                                   //delete uneeded mesh
    {
      delete pmesh;
      pmesh = NULL;
    }

    //check the rank in MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank != 0)
    {
      //master must be rank zero
      throw ProjectorException(PROJECTOR_ERROR_UNKOWN);
    }


        
    //branch on whether to have a slave store locally or not
    if (!slavelocal)
    {
      if (!projectnoslavelocal(progress))
          throw ProjectorException(PROJECTOR_ERROR_UNKOWN);
    }
    else
    {
      if (!projectslavelocal(progress))
        throw ProjectorException(PROJECTOR_ERROR_UNKOWN);
    }
  }
  catch(...)
  {
    if (pmesh)
      delete pmesh;
    throw ProjectorException(PROJECTOR_ERROR_UNKOWN);
  }

}
    

//*******************************************************************
void MpiProjector::setInputFile(std::string & ininfile) 
throw(ProjectorException)
{
  Projector::setInputFile(ininfile);
  inputfilename = ininfile;
}


//********************************************************************
bool MpiProjector::sendSlaveSetup(int rank) throw()
{
  char tempbuffer[100];
  int temp;
  unsigned char * buf(0);
  int bufsize(0), tempsize(0);
  int position(0);
  try
  {
    //calculate the buffersize
    MPI_Pack_size(200, MPI_CHAR, MPI_COMM_WORLD, &tempsize);
    bufsize += tempsize;
    MPI_Pack_size(2, MPI_LONG, MPI_COMM_WORLD, &tempsize);
    bufsize += tempsize;
    MPI_Pack_size(23, MPI_DOUBLE, MPI_COMM_WORLD, &tempsize);
    bufsize += tempsize;
    MPI_Pack_size(8, MPI_INT, MPI_COMM_WORLD, &tempsize);
    bufsize += tempsize;
    
    if (!(buf = new (std::nothrow) unsigned char[bufsize]))
      throw std::bad_alloc();

    position = 0;
    //pack the stuff
    strcpy(tempbuffer, inputfilename.c_str());
    //pack the input file name
    MPI_Pack(tempbuffer, 100, MPI_CHAR,
            buf, bufsize, &position, MPI_COMM_WORLD);
    //image metrics
    MPI_Pack(&newheight, 1, MPI_LONG,
            buf, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&newwidth, 1, MPI_LONG,
            buf, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&newscale.x, 1, MPI_DOUBLE,
            buf, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&newscale.y, 1, MPI_DOUBLE,
            buf, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&outRect.left, 1, MPI_DOUBLE,
            buf, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&outRect.top, 1, MPI_DOUBLE,
            buf, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&outRect.bottom, 1, MPI_DOUBLE,
            buf, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&outRect.right, 1, MPI_DOUBLE,
            buf, bufsize, &position, MPI_COMM_WORLD);
   
    if (slavelocal)
      temp = 1;
    else
      temp = 0;

    //pack slave local
    MPI_Pack(&temp, 1, MPI_INT,
             buf, bufsize, &position, MPI_COMM_WORLD);

    //pack the max chunksize
    MPI_Pack(&maxchunk, 1, MPI_INT,
            buf, bufsize, &position, MPI_COMM_WORLD);

    strcpy(tempbuffer, slavelocalpath.c_str());
    MPI_Pack(tempbuffer, 100, MPI_CHAR,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    //pack pmesh info
    MPI_Pack(&pmeshsize, 1, MPI_INT,
            buf, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&pmeshname, 1, MPI_INT,
            buf, bufsize, &position, MPI_COMM_WORLD);
    
    //pack the projection parameters
    MPI_Pack(reinterpret_cast<int *>(&Params.projtype), 1, MPI_INT,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(reinterpret_cast<int *>(&Params.datum), 1, MPI_INT,
             buf, bufsize, &position, MPI_COMM_WORLD);

    MPI_Pack(reinterpret_cast<int *>(&Params.unit), 1, MPI_INT,
            buf, bufsize, &position, MPI_COMM_WORLD);
   
    MPI_Pack(&Params.StdParallel1, 1, MPI_DOUBLE,
              buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.StdParallel2, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);

    MPI_Pack(&Params.NatOriginLong, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);

    MPI_Pack(&Params.NatOriginLat, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.FalseOriginLong, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.FalseOriginLat, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.FalseOriginEasting, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);   
    
    MPI_Pack(&Params.FalseOriginNorthing, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.CenterLong, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.CenterLat, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.CenterEasting, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.CenterNorthing, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);

    MPI_Pack(&Params.ScaleAtNatOrigin, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.AzimuthAngle, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.StraightVertPoleLong, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.zone, 1, MPI_INT,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.FalseEasting, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    MPI_Pack(&Params.FalseNorthing, 1, MPI_DOUBLE,
             buf, bufsize, &position, MPI_COMM_WORLD);
    
    
    MPI_Send(buf, position, MPI_PACKED, rank,
                 WORK_MSG, MPI_COMM_WORLD);

    delete [] buf;

    return true;
  }
  catch(...)
  {
    delete [] buf;
    //don't do anything
    return false;
  }
}

//*******************************************************
bool MpiProjector::projectslavelocal(BaseProgress * progress)
  throw(ProjectorException)
{
  //not supported yet..
  return false;
  /*
  int  bufferid(0), len(0), tag(0),    //pvm tags
    temptid(0);
  long int ycounter(1), endofchunk(0), beginofchunk(0);
  long int countmax(0);      //this is the total number of scanlines sent
  long int chunkssent(0);    //this is the number of chunks sent
  int chunkdif(maxchunk-minchunk);
  long int * chunknode(0);   //represents a node in the chunk queue
  std::queue<long int *> chunkqueue;
  unsigned char * scanline(0);
  std::strstream tempstream;

  try
  {
     //init the status progress
    if (progress)
    {
      tempstream << "Reprojecting " << newheight << " lines." << std::ends;
      tempstream.freeze(0);
      progress->init(tempstream.str(),
                      NULL,
                      "Done.",
                      newheight, 
                      29);
      tempstream.seekp(0);
      progress->start();  //start the progress
    }

    //init the random number generator
    if (sequencemethod == 2)
      srand48(time(NULL));
    

    while (countmax < newheight)
    {
      beginofchunk = countmax;    //set the beginning line
      //build the message
      switch(sequencemethod)
      {
      case 0:
        endofchunk = countmax + maxchunk-1;
        break;
      case 1:
        ++ycounter;
        if (ycounter >= sequencesize)
          ycounter = 0;
        endofchunk = countmax + sequence[ycounter]-1;
        break;
      case 2:
        endofchunk = countmax + static_cast<int>(drand48()*chunkdif
                                                 + minchunk) -1;
        break;
      }
      
                                                            
      //check the newheight
      if (endofchunk >= newheight)
      {
        endofchunk = newheight-1;
      }

      //update the countmax
      countmax += (endofchunk - beginofchunk) + 1;
      
      //create the chunk node
      if (!(chunknode = new (std::nothrow) long int [3]))
        throw std::bad_alloc();
      
      //store chunk dimensions
      chunknode[0] = beginofchunk; 
      chunknode[1] = endofchunk;

      //update the number of chunks sent
      ++chunkssent;
              
      //update the output
      if (progress && !((beginofchunk) % 11))
        progress->update(beginofchunk);
      
      bufferid = pvm_recv(-1, -1);      //blocking wait for all messages
      pvm_bufinfo(bufferid, &len, &tag, &temptid);
      

      switch(tag)
      {
      case SETUP_MSG:
        //pack the info in
        sendSlaveSetup();
        //add the first bit of work
        pvm_pklong(&(beginofchunk), 1, 1);
        pvm_pklong(&(endofchunk), 1, 1);
        //store that this node has this chunk
        chunknode[2] = temptid;
        chunkqueue.push(chunknode);

        //send the message
        pvm_send(temptid, SETUP_MSG);
        break;
      case WORK_MSG:
        //pack the next work
        pvm_initsend(PvmDataDefault);
        pvm_pklong(&(beginofchunk), 1, 1);
        pvm_pklong(&endofchunk, 1, 1);     //last chunk
        
        chunknode[2] = temptid;
        chunkqueue.push(chunknode);
        pvm_send(temptid, WORK_MSG);
        break;
      case ERROR_MSG:
      default:
        //slave sent error
        throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
      }
    }

    //first tell the slaves to projecting loop
    for (ycounter = 0; ycounter < numofslaves; ++ycounter)
    {
      bufferid = pvm_recv(-1, -1);      //blocking wait for all msgs
      pvm_bufinfo(bufferid, &len, &tag, &temptid);
      if (tag == ERROR_MSG)
      {
        throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
      }

      pvm_initsend(PvmDataDefault);
      pvm_send(childtid[ycounter], EXIT_MSG);   //tell slave to exit

    }

    
    if (progress)
      progress->done();

    //create the scaline
    if (!(scanline = new unsigned char[newwidth*spp]))
      throw std::bad_alloc();
    
    if (progress)
    {
      tempstream << "Writing " << newheight << " lines." << std::ends;
      tempstream.freeze(0);
      progress->init(tempstream.str(),
                     NULL,
                     "Done.",
                     newheight, 
                     29);
      tempstream.seekp(0);
      progress->start();  //start the progress
    }

    //finish writting scanlines
    while (chunkqueue.size())
    {
      chunknode = chunkqueue.front();
      
      //get the number of scanlines in this chunk
      beginofchunk = chunknode[0];
      endofchunk = chunknode[1];

      //update the output
      if (progress && !((beginofchunk) % 11))
        progress->update(beginofchunk);

      //get the slave id
      temptid = chunknode[2];
      
      for (ycounter = beginofchunk; ycounter <= endofchunk; ++ycounter)
      {
        //ask for a scanline
        pvm_initsend(PvmDataDefault);
        pvm_send(temptid, WORK_MSG);
        
        bufferid = pvm_recv(-1, -1);      //blocking wait for all msgs
        pvm_bufinfo(bufferid, &len, &tag, &temptid);
        
        if (tag == ERROR_MSG)
        {
          throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
        }
        
        if(pvm_upkbyte(reinterpret_cast<char*>(scanline),
                       newwidth*spp, 1) < 0)
          throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
        //write it
        out->putRawScanline(ycounter, scanline); 
      }
      
      //break it out of the loop
      pvm_initsend(PvmDataDefault);
      pvm_send(temptid, WORK_MSG);
      

      chunkqueue.pop();         //remove the chunk
      delete []  chunknode;     //delete the chunk node
    }

    if (progress)
      progress->done();

    writer.removeImage(0);                      //flush the output image
    out = NULL;
    return true;
  }
  catch(...)
  {
    //see if there is anthing left in the queue
    while (chunkqueue.size())
    {
      chunknode = chunkqueue.front();
      chunkqueue.pop();
      delete [] chunknode;
    }
      

    return false;
  }
  */
}

//*******************************************************
bool MpiProjector::projectnoslavelocal(BaseProgress * progress)
  throw(ProjectorException)
{
  int msize(0), membersize(0), position(0); 
  long int buffersize(0);
  MPI_Status status;         
  unsigned char * buffer(0); //the buffer for sending data
  long int ycounter(1), endofchunk(0), beginofchunk(0);
  long int countmax(0);      //this is the total number of scanlines sent
  long int chunkssent(0);    //this is the number of chunks sent
  long int chunksgot(0);     //this is the number of chunks that we have got
  unsigned int chunkdif(maxchunk-minchunk);
  Stitcher * mystitch(0);    //this is the sticher pointer (if we use it)
    

  try
  {
     //init the status progress
    if (progress)
    {
      std::strstream tempstream;

      tempstream << "Reprojecting " << newheight << " lines." << std::ends;
      tempstream.freeze(0);
      progress->init(tempstream.str(),
                      NULL,
                      "Done.",
                      newheight, 
                      29);
      progress->start();  //start the progress
    }

    if (stitcher)  //see if we want a stitcher
    {
      //creates the stitcher thread
      if (!(mystitch = new (std::nothrow) Stitcher(out)))
        
        throw std::bad_alloc();
    }

    //init the random number generator
    if (sequencemethod == 2)
      srand48(time(NULL));
    
    
    //figure out the maximum buffer size based on the system;
    MPI_Pack_size(2, MPI_LONG, MPI_COMM_WORLD, &membersize);
    buffersize+=membersize;
    MPI_Pack_size(maxchunk*newwidth*spp, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD,
                  &membersize);
    buffersize+=membersize;

    
    //ask for the buffer
    if (!(buffer = new unsigned char[buffersize]))
      throw std::bad_alloc();
    


    while (countmax < newheight)
    {
      beginofchunk = countmax;    //set the beginning line
      //build the message
      switch(sequencemethod)
      {
      case 0:
        endofchunk = countmax + maxchunk -1;
        break;
      case 1:
        ++ycounter;
        if (ycounter >= sequencesize)
          ycounter = 0;
        endofchunk = countmax + sequence[ycounter]-1;
        break;
      case 2:
        endofchunk = countmax + static_cast<int>(drand48()*chunkdif
                                                 + minchunk) -1;
        break;
      }
      
                                                            
      //check the newheight
      if (endofchunk >= newheight)
      {
        endofchunk = newheight-1;
      }

      //update the countmax
      countmax += (endofchunk - beginofchunk) + 1;
      
      //update the number of chunks sent
      ++chunkssent;
              
      //update the output
      if (progress && !((beginofchunk) % 11))
        progress->update(beginofchunk);
      
      //do a blocking wait for any message.
      MPI_Recv(buffer, buffersize, MPI_PACKED, MPI_ANY_SOURCE,
               MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      
      

      switch(status.MPI_TAG)
      {
      case SETUP_MSG:
        //pack the info in and send to slave
        sendSlaveSetup(status.MPI_SOURCE);
        //add the first bit of work
        position = 0;
        MPI_Pack(&(beginofchunk), 1, MPI_LONG, buffer, buffersize,
                 &position, MPI_COMM_WORLD);
        MPI_Pack(&(endofchunk), 1, MPI_LONG, buffer, buffersize,
                 &position, MPI_COMM_WORLD);
        //send to the slave
        MPI_Send(buffer, position, MPI_PACKED, status.MPI_SOURCE,
                 WORK_MSG, MPI_COMM_WORLD);
        break;
      case WORK_MSG:
        //unpack the scanline
        ++chunksgot; //got a chunk
        
        MPI_Get_count(&status, MPI_PACKED, &msize);
        if (stitcher)
        {
          sendStitcher(mystitch, buffer, msize);
        }
        else
          unpackScanline(buffer, msize);
        //pack the next work
        position = 0;
        MPI_Pack(&(beginofchunk), 1, MPI_LONG, buffer, buffersize,
                 &position, MPI_COMM_WORLD);
        MPI_Pack(&(endofchunk), 1, MPI_LONG, buffer, buffersize,
                 &position, MPI_COMM_WORLD);
        //send to the slave
        MPI_Send(buffer, position, MPI_PACKED, status.MPI_SOURCE,
                 WORK_MSG, MPI_COMM_WORLD);
        break;
      case ERROR_MSG:
      default:
         throw ProjectorException(PROJECTOR_ERROR_BADINPUT);
      }
    }

    //finish writting scanlines
    for (ycounter = chunksgot; ycounter < chunkssent; ++ycounter)
    {
      MPI_Recv(buffer, buffersize, MPI_PACKED, MPI_ANY_SOURCE,
               MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (status.MPI_TAG == WORK_MSG)
      {
        if (stitcher)
        {
          sendStitcher(mystitch, buffer, buffersize);
        }
        else
          unpackScanline(buffer, buffersize);
      }
      else
      {
         throw ProjectorException(PROJECTOR_ERROR_BADINPUT); 
         //should never happen
      }
    }
    
    //slave termination
    for (ycounter = 0; ycounter < numofslaves; ++ycounter)
    {
      MPI_Send(buffer, buffersize, MPI_PACKED, ycounter+1,
               EXIT_MSG, MPI_COMM_WORLD);
    }

    if (progress)
      progress->done();

    if (stitcher)
    {
      mystitch->wait();
      //remove the stitcher
      delete mystitch;
    }

    writer.removeImage(0);                      //flush the output image

    out = NULL;
    return true;
  }
  catch(...)
  {
    if (stitcher)
    {
      delete mystitch;                          //should stop the stitcher
      mystitch = NULL;
    }

    return false;
  }
}

//******************************************************
void MpiProjector::setStitcher(bool institcher) throw()
{
  stitcher = institcher;
}

//******************************************************
bool MpiProjector::getStitcher() const throw()
{
  return stitcher;
}


//***************************************************************
void MpiProjector::setSequence(const int * insequence,
                 const int & insequencesize) throw(std::bad_alloc)
{
  int counter(0);
  int max(0);                        //finds the max of the squence
  resetSequencing();                 //reset the sequencing
  
  sequencesize = insequencesize;     //set the size
  
  //create the sequence
  if(!(sequence = new (std::nothrow) int [sequencesize]))
    throw std::bad_alloc();
  
  //copy the sequence
  for (; counter < insequencesize; ++counter)
  {
    sequence[counter] = insequence[counter];
    if (sequence[counter] > max)
      max = sequence[counter];
  }

  maxchunk = max;   //set the maxchunksize
  //set the sequence method
  sequencemethod = 1;
}

  
//****************************************************************
void MpiProjector::setRandomSequence(const int & inminchunk,
                       const int & inmaxchunk) throw()
{
  resetSequencing(); 

  //set the chunksizes and sequence method
  minchunk = inminchunk;
  maxchunk = inmaxchunk;
  sequencemethod = 2;
}


  
//****************************************************************
void MpiProjector::resetSequencing() throw()
{
  delete [] sequence;
  sequence = 0;
  sequencesize = 0;
  minchunk = 0;
  maxchunk = 1;
  sequencemethod = 0;
}

//****************************************************************
void MpiProjector::
setSlaveLocalDir(const std::string & inslavelocalpath) throw()
{
  slavelocalpath = inslavelocalpath;
}

//**************************************************************
std::string MpiProjector::getSlaveStoreLocalDir() const throw()
{
  return slavelocalpath;
}

//*******************************************************
long int MpiProjector::unpackScanline(unsigned char * buffer,
                                      long int buffersize) throw()
{
  unsigned char * tempscanline=NULL;
  long int scanlinenumber(0);
  long int endscanline(0);
  long int counter(0); 
  int position(0);
  try
  {
    MPI_Unpack(buffer, buffersize, &position,
               &scanlinenumber, 1, MPI_LONG, MPI_COMM_WORLD);
    MPI_Unpack(buffer, buffersize, &position,
               &endscanline, 1, MPI_LONG, MPI_COMM_WORLD);
    
    //create enough space to hold the entire chunk
    if (!(tempscanline = new (std::nothrow) unsigned char 
          [(endscanline-scanlinenumber + 1)*newwidth*spp]))
      throw std::bad_alloc();

    //unpack all of the scanlines
    MPI_Unpack(buffer, buffersize, &position,
               tempscanline, (endscanline-scanlinenumber + 1)*newwidth*spp, 
               MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);

    for (counter = scanlinenumber; counter <= endscanline; ++counter)
    {
      //write it
      out->putRawScanline(counter, 
              &(tempscanline[newwidth*spp*(counter-scanlinenumber)]) ); 
    }

    delete [] tempscanline;
  
    return endscanline-scanlinenumber + 1;

  }
  catch(...)
  {
    return -1;
    //something bad happened 
    if (tempscanline)
      delete [] tempscanline;
  }
}

//*********************************************************************
long int MpiProjector::sendStitcher(Stitcher * mystitch, 
                                    unsigned char * buffer,
                                    long int buffersize) throw()
{
  unsigned char * tempscanline(0);
  long int scanlinenumber(0);
  long int endscanline(0);
  StitcherNode * temp(0);
  int position(0);
  try
  {
    MPI_Unpack(buffer, buffersize, &position,
               &scanlinenumber, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, buffersize, &position,
               &endscanline, 1, MPI_INT, MPI_COMM_WORLD);
    
    //create enough space to hold the entire chunk
    if (!(tempscanline = new (std::nothrow) unsigned char 
          [(endscanline-scanlinenumber + 1)*newwidth*spp]))
      throw std::bad_alloc();

    //unpack all of the scanlines
    MPI_Unpack(buffer, buffersize, &position,
               tempscanline, (endscanline-scanlinenumber + 1)*newwidth*spp, 
               MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);

    //create the stitcher node
    if (!(temp = new (std::nothrow) StitcherNode(scanlinenumber,
                                                 endscanline,
                                                 tempscanline)))
      throw std::bad_alloc();
    
    
    //add the chunk to the stitchers queue
    mystitch->add(temp);
    
    return endscanline-scanlinenumber + 1;
    
  }
  catch(...)
  {
    return -1;
    delete temp;
    delete [] tempscanline;
  }
}


#endif





