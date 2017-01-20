#include <phasta.h>
#include <phastaIO.h>
#include <phstream.h>
void testReadSolution(grstream grs) {
      FILE* f = openGRStreamRead(grs,"restart");
      int params[128];
      readHeader(f,"solution",params,3,"binary");
      int size = params[0]*params[1];
      double* data = (double*) calloc(size,sizeof(double));
      readDataBlock(f,data,size,"double","binary");
      free(data);
      fclose(f);
}
