
#ifdef _WIN32
#define __GNU_LIBRARY__
#endif

#ifdef _WIN32
#pragma warning( disable : 4291 ) // Disable VC warning messages for
                                  // new(nothrow)
#endif


#include <iostream>
#include <fstream>
#include <time.h>
#include "PVFSProjector.h"
#include "CLineProgress.h"
#include "ProjUtil.h"
#include "inparms.h"
using namespace ProjLib;


int main(int argc, char *argv[])
{
  MpiProjector * projector(0);   //if this the master then use this
  CLineProgress progress;        //progress for command line output
  clock_t start = 0, finish = 0; //for time files
  std::ofstream out;             //for output time files
  Projection * outproj = NULL;   //output projection
  inputparm inparms;             //for master setup and parameters



  try
  {
    if (argc == 1)
    {
      //we are the master with out a flat file so get input
      inparms.getinput();
    }
    else
    {
      if (!inparms.read_parm_file(argv[1]))
      {
        std::cout << "Unable to read the configuration file" << std::endl;
        return 0;
      }
      
    }

    //generate the output projection based on the flat file
    outproj = SetProjection(inparms.parameterfile);
    if (!outproj)
    {
      std::cout << "Could not create the output projection!" << std::endl;
      return 0;
    }
    
    //check for the time file
    if (inparms.timefile)
    {
      out.open(inparms.logname.c_str()); //open the ouput file
      start = time(NULL);                //get the start
      out << start << std::endl;         //output the start
    }

    
    if (inparms.numPartitions > 0)
    {
      if (!(projector = new (std::nothrow) PVFSProjector))
        throw std::bad_alloc();
      
      dynamic_cast<PVFSProjector *>(projector)->setPartitionNumber
        (inparms.numPartitions);
    }
    else
      //create the projector
      if(!(projector = new (std::nothrow) MpiProjector))
        throw std::bad_alloc();

    projector->setOutputFileName(inparms.outfile_name);
    projector->setOutputProjection(outproj);
    projector->setPmeshName(inparms.pmeshname);
    projector->setPmeshSize(inparms.pmeshsize);
    if (inparms.chunksize > 0)
      projector->setChunkSize(inparms.chunksize);
    else
      projector->setChunkSize(1);

    projector->setSlaveStoreLocal(inparms.storelocal);

    projector->setStitcher(inparms.stitcher);

    if (!inparms.samescale)
      projector->setOutputScale(inparms.newscale);
    else
      projector->setSameScale(true);
    
    projector->setInputFile(inparms.filename);
    projector->setNumberOfSlaves(inparms.numofslaves);

    if(inparms.numofslaves > 0)
      MPI_Init(&argc, &argv);           //start MPI

    projector->project(&progress);

    //see if we need to write the time file
    if (inparms.timefile)
    {
      finish = time(NULL);
      out << finish << std::endl;
      out << (finish - start);
      out.close();
    }
    
    //done with the output projection
    delete outproj;
    delete projector;
    outproj = NULL;
    
    if(inparms.numofslaves > 0)
      MPI_Finalize();           //start MPI


  }
  catch(ProjectorException & pe)
  {
    std::cout << pe.getExceptionMessage() << std::endl;
    delete projector;
    delete outproj;
    outproj = NULL;
  }
  catch(USGSImageLib::ImageException & ie)
  {
    std::string error;

    ie.getString(error);
    std::cout << error << std::endl;
    delete projector;
    delete outproj;
    outproj = NULL;
  }
  catch(...)
  {
    std::cout << "You are the weakest link. Good-bye." << std::endl;
    delete projector;
    delete outproj;
    outproj = NULL;
  }

  return 0;  //done

}










