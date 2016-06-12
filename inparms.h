//*************************************************************
//Header File
//Inputparm class reads input parameters for the projector
//Author: Keith Athmer
//*************************************************************

#include <string>
#include "MathLib/Point.h"

class inputparm
{
public:

  inputparm();         //Default Constructor:
   
  ~inputparm();        //Destructor:

  void getinput();  
  
  bool write_parm_file(const std::string& outfilename);

  bool read_parm_file(const std::string& infilename);


	
  int numofslaves;
  bool samescale;
  bool timefile;
  MathLib::Point newscale;
  int pmeshsize, pmeshname;
  std::string logname, filename, parameterfile, outfile_name;
  /** These were added for the updated options in the new projector program
      CBB 5/8/2001 **/ 
  int chunksize;                  //the chunksize to process with (default 1)
  bool storelocal;                //whether or not to store chunks locally
                                  //(default no)
  bool stitcher;                  //whether or not to use the stitcher on
                                  //as the master (default no)
  int numPartitions;              //the pvfs partitions

protected:

  //this funciton reads characters (until the enter is pressed) and
  //fills a buffer (with the enter stripped).
  //If the buffer is empty the false is returned (user hit enter)
  bool inputparm::read(std::string & inRead) throw();

    
};  //end of class
