//**************************************************************************
//Implementation File
//**************************************************************************
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
#include <cstdlib>
#include <queue>
#include "inparms.h"
#include "MiscUtils/cmp_nocase.h"


inputparm::inputparm()
{
  numofslaves = 0;
  samescale = false;
  timefile = false;
  pmeshname = 0;
  pmeshsize = 4;
  newscale.x = 0;
  newscale.y = 0;  
  chunksize = 0;  
  storelocal = false;      
  stitcher = false; 
  numPartitions = 0;
}//constructor

inputparm::~inputparm()
{
}//default destrutor


//***************************************************************************
//Get input parameters
//
//This function ask the user information to set up the pameter list
//
//***************************************************************************
void inputparm::getinput()
{
  std::string parm_file;
  std::string inbuf;        //input buffer
  int pmsh_ans = 0;

  //first clear the stdin buffer if there is anthing in it
  std::cin.ignore(std::cin.rdbuf()->in_avail());

  //Get input 
  std::cout << "Enter the number of slaves to use (default 0)" << std::endl;
  std::getline(std::cin, inbuf);

  if (!inbuf.size())
  {
    numofslaves = 0;
  }
  else
  {
    numofslaves = std::atoi(inbuf.c_str());
  }

  std::cout << "Please enter the chunksize (1 for default)." << std::endl;
  std::getline(std::cin, inbuf);

  if(!inbuf.size())
  {
    chunksize = 1;
  }
  else
  {
    chunksize = std::atoi(inbuf.c_str());;
  }

  std::cout << "Do you want to have the slaves store data locally? (Y/N)"
            << " (default N)" << std::endl;
  std::getline(std::cin, inbuf);

  if (!inbuf.size())
  {
    storelocal = false;
  }
  else
  {
    if (!MiscUtils::cmp_nocase(inbuf, "Y"))
    {
      storelocal = true;
    }
    else
      storelocal = false;
  }

  std::cout << "Do you wish to use the multithreaded master? (Y/N)" 
            << " (Default N)" << std::endl;
  std::getline(std::cin, inbuf);

  if (!inbuf.size())
  {
    stitcher = false;
  }
  else
  {
    if (!MiscUtils::cmp_nocase(inbuf, "Y"))
    {
      stitcher = true;
    }
    else
      stitcher = false;
  }

  std::cout << "Do you want output in the same scale? (y or n) (default y)"
            << std::endl;
  std::getline(std::cin, inbuf);

  if (!inbuf.size())
  {
    samescale = true;
    newscale.x = 0;
    newscale.y = 0;
  }
  else
  {
    if (!MiscUtils::cmp_nocase(inbuf,"N"))
    {
      samescale = false;
      std::cout << "Enter the new x scale factor." << std::endl;
      std::cin >> newscale.x;
    
      std::cout << "Enter the new y scale factor." << std::endl;
      std::cin >> newscale.y;

      //get rid of the excess std::cin buffer stuff
      std::cin.ignore(std::cin.rdbuf()->in_avail());
    }
    else
    {
      samescale = true;
      newscale.x = 0;
      newscale.y = 0;
    }
  }

  std::cout << "Do you want to log the processing time? (y or n) (default n)" 
            << std::endl;
  std::getline(std::cin, inbuf);
  if (inbuf.size())
  {
    if(!MiscUtils::cmp_nocase(inbuf,"Y"))
    {
     
      timefile = true;
      std::cout << "Enter the name of the logfile for timings." << std::endl;
      std::getline(std::cin, inbuf);
      logname = inbuf;
    }
    else
    {
      timefile = false;
      logname = "logname";
    }
  }
  else
  {
    timefile = false;
    logname = "logname";
  }

  std::cout << "Do you want to use PVFS partitioning? {y/n} (default y)"
            << std::endl;
  std::getline(std::cin, inbuf);
  if (inbuf.size())
  {
    if (!MiscUtils::cmp_nocase(inbuf, "Y"))
    {
      std::cout << "Enter the number of partitions:" << std::endl;
      std::cin >> numPartitions;
      std::cin.ignore(std::cin.rdbuf()->in_avail());
    }
  }
  
  std::cout << "Choose one of the following for the Pmesh name:" << std::endl;
  std::cout << "0=None(Default),1=LeastSqrs, 2=Bilinear, 3=bicubic." 
            << std::endl;
  std::getline(std::cin, inbuf);
  if(inbuf.size())
  {
    if (!MiscUtils::cmp_nocase(inbuf, "0"))
      pmeshname = 0;
    if (!MiscUtils::cmp_nocase(inbuf, "1"))
      pmeshname = 6;
    if (!MiscUtils::cmp_nocase(inbuf, "2"))
      pmeshname = 8;
    if (!MiscUtils::cmp_nocase(inbuf, "3"))
      pmeshname = 10;
  }
  else
    pmeshname = 0;

  if (pmeshname > 0)
  {
    std::cout << "Enter the Pmesh Size." << std::endl;
    std::cin >> pmsh_ans;
  
    if ( pmsh_ans < 4)
    {
      pmeshsize = 4;
    }
    else
    {
      pmeshsize = pmsh_ans;
    }
  }
  else
    pmeshsize = 4;


  std::cout << "Enter the name of the input file." << std::endl;
  std::cin >> filename;

  std::cout << "Enter the name of parameter file." << std::endl;
  std::cin >> parameterfile;

  std::cout << "Enter the name of the output file." << std::endl;
  std::cin >> outfile_name;

  std::cout << "Do you want to save these parameters? (y or n)" << std::endl;
  std::cin >> inbuf;

  if (!MiscUtils::cmp_nocase(inbuf, "Y"))
  {
    std::cout << "Enter the name of the parameter file to save the settings" 
              << std::endl;
    std::cin >> parm_file;
    write_parm_file(parm_file);
  }


}


//****************************************************************************
//Write Function
//
//This function take a file name in to write the input parameters to.
//It is called by the getinputs function
//
//****************************************************************************
bool inputparm::write_parm_file(const std::string& outfilename)
{
  
  std::ofstream outfile;

  outfile.open(outfilename.c_str());

   if (! outfile)
  {
    std::cerr << "can't open the output file" << std::endl;
    return false;
  }

  outfile << numofslaves << std::endl;
  outfile << samescale << std::endl;
  outfile << newscale.x << std::endl;
  outfile << newscale.y << std::endl;
  outfile << timefile << std::endl;
  outfile << logname << std::endl;
  outfile << pmeshsize << std::endl;
  outfile << pmeshname << std::endl;
  outfile << filename << std::endl;
  outfile << parameterfile << std::endl;
  outfile << outfile_name << std::endl;
  outfile << chunksize << std::endl;
  outfile << storelocal << std::endl;
  outfile << stitcher << std::endl;
  outfile << numPartitions << std::endl;
  outfile.close();

  return true;

}
  

//****************************************************************************
//Read Function
//
//This function reads the input parmeters from an input file 
//****************************************************************************

bool inputparm::read_parm_file(const std::string &infilename)
{
 
  std::ifstream infile;
 
  infile.open (infilename.c_str());
  
  if (! infile)
  {
    std::cerr << "can't open input file" << std::endl;
    return false;
  }

  infile >> numofslaves;
  infile >> samescale;
  infile >> newscale.x;
  infile >> newscale.y;
  infile >> timefile;
  infile >> logname;
  infile >> pmeshsize;
  infile >> pmeshname;
  infile >> filename;
  infile >> parameterfile;
  infile >> outfile_name;
  infile >> chunksize;
  infile >> storelocal;
  infile >> stitcher;
  infile >> numPartitions;
  infile.close();
  
  return true;
    
}







