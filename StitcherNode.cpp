/**
 *This is the implementation file for the stitchernode
 **/

#ifndef STITCHERNODE_CPP
#define STITCHERNODE_CPP

#include "StitcherNode.h"


//********************************************************************
StitcherNode::
StitcherNode(const int & instart, const int & inend, 
             unsigned char * indata) : start(instart),
                                             end(inend),
                                             data(indata)
{}

//*******************************************************************
StitcherNode::~StitcherNode()
{
  //as promised delete the data
  delete [] data;
}

//********************************************************************
void StitcherNode::
setData(const int & instart, const int & inend,
        unsigned char * indata) throw()
{
  start = instart;
  end = inend;
  data = indata;
}

#endif




