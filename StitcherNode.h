/**
 *This is the node that goes in the stitcher queue that 
 *represents a chunk of output image to be written to disc.
 **/

#ifndef STITCHERNODE_H
#define STITCHERNODE_H




//This represents a node in the stitcher's que
class StitcherNode
{
public:
  
  /**
   * This is the main constructor for the sticher node
   * Start chunk is the starting scanline number of the queue and
   * end chunk is the ending scanline of the que
   * and data is the actual data to be written.
   **/
  StitcherNode(const int & instart, const int & inend, 
               unsigned char * indata);

  /**
   *Destructor deletes the chunk of memory that is passed into it
   **/
  virtual ~StitcherNode();

  /**
   *This sets the chunk data
   **/
  void setData(const int & instart, const int & inend,
                unsigned char * indata) throw();
  
  /**
   *This returns a pointer to the scanline data
   **/
  unsigned char * getData() throw();
  
  /**
   *This returns the start scanline number
   **/
  int getStart() const throw();


  /**
   *This returns the end scanline number
   **/
  int getEnd() const throw();
  
  
 private:
  
  int start, end;
  unsigned char * data;
};


//inline functions

//***************************************************
inline unsigned char * StitcherNode::getData() throw()
{
  return data;
}

//**************************************************
inline int StitcherNode::getStart() const throw()
{
  return start;
}

//**************************************************
inline int StitcherNode::getEnd() const throw()
{
  return end;
}


#endif

