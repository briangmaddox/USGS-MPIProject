


#include <iostream>
#include "MpiProjectorSlave.h"

int main(int argc, char **argv)
{
  MpiProjectorSlave slave;
  int temp(0);

  //start Mpi
  MPI_Init(&argc, &argv);

  try
  {
    temp = slave.connect();
    MPI_Finalize();
    return temp;
  }
  catch(...)
  {
    MPI_Finalize();
    return 0;
  }
}
