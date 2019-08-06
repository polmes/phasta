#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
//#define OMPI_SKIP_MPICXX 1 //Added in the CMakeList.txt file
#include <mpi.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "phastaIO.h"
enum {
  DIR_FANOUT = 2048
};

inline int
cscompare( const char teststring[],
	   const char targetstring[] )
{
  char* s1 = const_cast<char*>(teststring);
  char* s2 = const_cast<char*>(targetstring);

  while( *s1 == ' ') s1++;
  while( *s2 == ' ') s2++;
  while( ( *s1 )
	 && ( *s2 )
	 && ( *s2 != '?')
	 && ( tolower( *s1 )==tolower( *s2 ) ) ) {
    s1++;
    s2++;
    while( *s1 == ' ') s1++;
    while( *s2 == ' ') s2++;
  }
  if ( !( *s1 ) || ( *s1 == '?') ) return 1;
  else return 0;
}

inline int 
computenitems(const int localpartid, const int fieldid, const int myrank, const char *fieldName, int ***para, const int intHeader, const int numVariables) {
// This routine computes the number of items in the data block based on 
// - the name of the fields 
// - the integers read from the header

  int nItems = -1;

  if (cscompare("nbc values",fieldName))
    nItems = para[localpartid][fieldid][0] * (numVariables+1);
  else if  (cscompare("nbc codes",fieldName))
    nItems = para[localpartid][fieldid][0] * 2;
  else if ( intHeader==1)
    nItems = para[localpartid][fieldid][0];
  else
    nItems = para[localpartid][fieldid][0] * para[localpartid][fieldid][1];
  return nItems;
}


int main(int argc, char *argv[]) {

  MPI_Init(&argc,&argv);

  int myrank, N_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &N_procs);

  FILE * pFile,* AccessoryFileHandle;
  char target[1024], pool[256], tempString[128],numpe[8],numstart[8];
  char * temp, * token;
  int fieldCompareMark;
  bool readField;
  int i, j, k, N_geombc_double, N_geombc_integer, N_restart_double;
  int N_restart_integer, N_steps, N_parts, N_files;

  pFile = fopen("./IO.O2N.input","r");
  if (pFile == NULL)
    printf("Error openning\n");

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_geombc_double = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_geombc_integer = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_restart_double = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_restart_integer = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  strncpy(numstart,temp,1);
  N_steps = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  strncpy(numpe,temp,1);
  N_parts = atoi(temp);

  if(myrank==0){
    printf("numpe is %d and start is %d\n",N_parts,N_steps);
  }

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_files = atoi(temp);

  double ***Dfield; 
  double ***FixedDfield; 
  double **CoordsOld; 
  double **CoordsNew;
  int ***Ifield;
  int ***paraD, ***paraI, *expectD, *expectI;
  char **fieldNameD, **fileTypeD, **dataTypeD, **headerTypeD;
  char **fieldNameI, **fileTypeI, **dataTypeI, **headerTypeI;

  int* WriteLockD = new int[N_geombc_double];
  int* WriteLockI = new int[N_geombc_integer];

  int nppp = N_parts/N_procs;
  int startpart = myrank * nppp +1;
  int endpart = startpart + nppp - 1;
  char gfname[64], numTemp[128];
  int iarray[10], igeom, isize,TempFileHandle;

  ////////////////////////////////////////
  // Test if the user has given the right  parameters related to the number of parts, procs and SyncIO files
  ////////////////////////////////////////

  if (N_parts%N_files!=0)
    {
      printf("Input error: number of parts should be a multiple of number of files!\n");
      printf("Please modify the IO.O2N.input file!\n");
      return 0;
    }
  if (N_procs%N_files!=0)
    {
      printf("Input error: number of procs should be a multiple of number of files!\n");
      printf("Please modify the IO.O2N.input file!\n");
      return 0;
    }
  if (N_parts%N_procs!=0)
    {
      printf("Input error: number of parts should be a multiple of number of procs!\n");
      printf("Please modify the IO.O2N.input file!\n");
      return 0;
    }

  /////////////////////////////////////
  // Create numpe.in and numstart.dat
  ////////////////////////////////////////

  if (myrank==0)
    {

//MR CHANGE
      bzero((void*)gfname,64);
      sprintf(gfname,"./%d-procs_case-SyncIO-%d-FM",N_parts,N_files);
      if(0<mkdir(gfname,0777)) { printf("ERROR - Could not create procs_case-SyncIO directory\n"); return 1; }
//MR CHANGE END

      bzero((void*)gfname,64);
      sprintf(gfname,"./%d-procs_case-SyncIO-%d-FM/numstart.dat",N_parts,N_files);
      AccessoryFileHandle=fopen(gfname,"w");
      fprintf(AccessoryFileHandle,"%d",N_steps);
      fclose(AccessoryFileHandle);

      bzero((void*)gfname,64);
      sprintf(gfname,"./%d-procs_case-SyncIO-%d-FM/numpe.in",N_parts,N_files);
      AccessoryFileHandle=fopen(gfname,"w");
      fprintf(AccessoryFileHandle,"%d",N_parts);
      fclose(AccessoryFileHandle);
    }

  ///////////////////////////////////////
  Dfield = new double**[nppp];  // first rank is parts per processor
  CoordsOld = new double*[nppp];
  CoordsNew = new double*[nppp];
  Ifield = new int**[nppp];

  paraD = new int**[nppp];
  paraI = new int**[nppp];
  /////////////////////////////////////

  int* interiorMark = new int[nppp];
  int* boundaryMark = new int[nppp];
  int* codesMark = new int[nppp];
  int* valuesMark = new int[nppp];
  int* numVariables = new int[nppp];
//  int numBoundaryFields[nppp], numInteriorFields[nppp];

  if(myrank==0){
    printf("Starting to read some headers in the restart.##.## and geombc.dat.## files\n");
  }

  // The number of variables does not vary from one part to another. Do not read this info from every part.
  // The ideal would be to ask only rank 0 to read only part 0 and broadcast the information to the other ranks. It saves a lot of time for very big meshes.
  // Right now, every rank will read the number of variables from the solution field of one part only, which is already better than reading this info from every part.
    int ithree = 3;
    j = 0;

    // Test if file exist in the procs_case directory
    bzero((void*)gfname,64);
    sprintf(gfname,"./%d-procs_case-1PPP-BM/restart.%d.%d", N_parts, N_steps, startpart+j);
    int subdir = -1;
    bool existsubdir = false;
    FILE *pfiletest;
    pfiletest = fopen(gfname, "r");
    if (pfiletest == NULL ) {
      // Test if file exist in the procs_case/subdir directory
      subdir = (startpart+j-1) / DIR_FANOUT;
      bzero((void*)gfname,64);
      sprintf(gfname,"./%d-procs_case-1PPP-BM/%d/restart.%d.%d",N_parts, subdir, N_steps, startpart+j);
      pfiletest = fopen(gfname, "r");
      if (pfiletest == NULL ) {
        printf("[%d] File %s does not exit - abort\n",myrank, gfname);
        abort();
      }
      else{
        existsubdir = true;
        fclose(pfiletest);
      }
    }
    else {
      fclose(pfiletest);
    }

    // Debug
    //printf("Rank: %d - subdir: %d - path: %s\n",myrank, subdir, gfname);

    openfile(gfname,"read",&TempFileHandle);
    readheader( &TempFileHandle,
		   "solution",
		   (void*)iarray,
		   &ithree,
		   "double",
		   "binary" );
    closefile(&TempFileHandle, "read");
    for ( j = 0; j < nppp; j++ )
      numVariables[j] = iarray[1]; //iarray[i] contains the number of variables from the header of the solution field

    MPI_Barrier(MPI_COMM_WORLD); //added by MR


  /////////////////////////
  for ( i = 0; i < nppp; i++ )
    {
      Dfield[i] = new double*[N_geombc_double]; // Space for the datablock of double format  Second rank is geomdouble
      paraD[i] = new int*[N_geombc_double];     // Integers in the header of each field of double format
      Ifield[i] = new int*[N_geombc_integer];   // Space for the datablock of integer format
      paraI[i] = new int*[N_geombc_integer];    // Integers in the header of each field of integer format

    }

  expectD = new int[N_geombc_double];          // Expected number of integers in the header for each field of double format
  expectI = new int[N_geombc_integer];         // Expected number of integers in the header for each field of integer format

  fieldNameD = new char*[N_geombc_double];     // Name of the field in double format
  fileTypeD = new char*[N_geombc_double];      // geombc or restart (useless if associated with geombc file but read)
  dataTypeD = new char*[N_geombc_double];      // Integer or double (useless if associated with double data but read)
  headerTypeD = new char*[N_geombc_double];    // block (means with data block) or header (just a header with no block)

  fieldNameI = new char*[N_geombc_integer];
  fileTypeI = new char*[N_geombc_integer];
  dataTypeI = new char*[N_geombc_integer];
  headerTypeI = new char*[N_geombc_integer];

  /////////////////////////

  for ( i = 0; i < N_geombc_double; i++ )
    {
      WriteLockD[i]=0;

      fieldNameD[i] = new char[128];
      fileTypeD[i] = new char[128];
      dataTypeD[i] = new char[128];
      headerTypeD[i] = new char[128];
    }

  for ( i = 0; i < N_geombc_integer; i++ )
    {
      WriteLockI[i]=0; //This take value 1 if the field requested in IO.input is not found (typically for tpblocks)

      fieldNameI[i] = new char[128];
      fileTypeI[i] = new char[128];
      dataTypeI[i] = new char[128];
      headerTypeI[i] = new char[128];
    }

  ////////////////////////////////////////////////////////////////
  // Reading IO.input        
  // temporary fix: in the new version the double and integer
  //                can mix and match to avoid the order confusion
  ////////////////////////////////////////////////////////////////

  char S1[128],S2[128],S3[128],S4[128],S5[128];
  int double_counter=0,int_counter=0;

  for ( i = 0; i < N_geombc_double+N_geombc_integer; i++ )
    {
      fgets( target, 1024, pFile );
      temp = strtok( target, ";" );
      token = strtok( temp, "," );
      strcpy( S1, token );
      token = strtok ( NULL, "," );
      strcpy( S2, token );
      token = strtok ( NULL, "," );
      strcpy( S3, token );
      token = strtok ( NULL, "," );
      strcpy( S4, token );
      token = strtok ( NULL, "," );
      strcpy( S5, token );

      //if (cscompare(S3,"double"))
      if (cscompare("double",S3))
	{
	  strcpy( fileTypeD[double_counter], S1 );
	  strcpy( fieldNameD[double_counter], S2 );
	  strcpy( dataTypeD[double_counter], S3 );
	  strcpy( headerTypeD[double_counter], S4 );
	  strcpy( numTemp, S5 );
	  expectD[double_counter] = atoi (numTemp);
	  double_counter++;
	}

      //if (cscompare(S3,"integer"))
      if (cscompare("integer",S3))
	{
	  strcpy( fileTypeI[int_counter], S1 );
	  strcpy( fieldNameI[int_counter], S2 );
	  strcpy( dataTypeI[int_counter], S3 );
	  strcpy( headerTypeI[int_counter], S4 );
	  strcpy( numTemp, S5 );
	  expectI[int_counter] = atoi (numTemp);
	  int_counter++;
	}
    }

  //////////////////////////////////////////////////////////////

  //for ( i = 0; i < N_geombc_double; i++) {
    //printf("%d %s %s %s %s %d\n", myrank, fileTypeD[i], fieldNameD[i], dataTypeD[i], headerTypeD[i], expectD[i]);
  //}


//  printf("rank %d is waiting\n",myrank);
  MPI_Barrier(MPI_COMM_WORLD); //already there

  if(myrank==0){
    printf("Starting to read some blocks (doubles) in the geombc.dat.## files\n");
  }

  for ( i = 0; i < nppp; i++ )
    {
      if(existsubdir) { 
        subdir = (startpart+i-1) / DIR_FANOUT;
        sprintf(gfname,"./%d-procs_case-1PPP-BM/%d/geombc.dat.%d",N_parts, subdir, startpart+i);
      }
      else
        sprintf(gfname,"./%d-procs_case-1PPP-BM/geombc.dat.%d",N_parts,startpart+i);

      openfile(gfname,"read",&igeom);

//      MPI_Barrier(MPI_COMM_WORLD);

      for ( j = 0; j < N_geombc_double; j++ )
	{

//	  printf("rank: %d - read double: %s\n",myrank,fieldNameD[j]);
	  paraD[i][j] = new int[expectD[j]];

	  for ( k = 0; k < 10; k++ )
	    iarray[k]=0;

	  iarray[0]=-1;
	  WriteLockD[j]=0;
	  readheader( &igeom,
		       fieldNameD[j],
		       (void*)iarray,
		       &expectD[j],
		       "double",
		       "binary" );
	  if ( iarray[0]==-1 )  // The field requested in IO.O2N.input has not been found (should be a tpblocks)
	      WriteLockD[j]=1;

//          printf("rank: %d - part: %d - field: %s - iarray: %d\n",myrank,startpart+i,fieldNameD[j],iarray[0]);

//	  MPI_Barrier(MPI_COMM_WORLD);

	  // get the parameter list in data header ...
	  // different fields have different ways to get this list ...
	  for ( k = 0; k < expectD[j]; k++ )
	    paraD[i][j][k] = iarray[k];

	  if ( WriteLockD[j]==1) // Put the value of the expected integers in the header to 0 when field not present
	    for ( k = 0; k < expectD[j]; k++ )
	      paraD[i][j][k] = 0;

/*          int iproc;
          for(iproc=0; iproc<N_procs; iproc++){
            MPI_Barrier(MPI_COMM_WORLD);
            if(iproc == myrank){
              printf(" iproc: %d ", iproc);
              printf("part: %d myrank: %d field: %s header: ",startpart+i,myrank,fieldNameD[j]);
              for ( k = 0; k < expectD[j]; k++ )
                printf(" %d ",iarray[k]);
              printf("\n");
            }
            usleep(100);
            MPI_Barrier(MPI_COMM_WORLD);
          }*/


	  if ( cscompare("block",headerTypeD[j]) )
	    {
/*	      if ( expectD[j]==1)
		isize = paraD[i][j][0];
	      else
		isize = paraD[i][j][0] * paraD[i][j][1];

	      if (cscompare("nbc values",fieldNameD[j]))
		isize = paraD[i][j][0] * (numVariables[i]+1);
*/

//              int test;
//              test = computenitems(i,j,myrank,fieldNameD[j],paraD,expectD[j],numVariables[i]);
//              printf("irank: %d fieldname: %s ParaD: %d\n",myrank,fieldNameD[j], paraD[0][0][0], numVariables[i]);
//              if(test != isize)
//                printf("PROBLEM fieldname: %s part: %d isize: %d test: %d\n",fieldNameD[j],startpart+i,isize,test); 
//              else 
//                printf("fieldname: %s part: %d isize: %d test: %d\n",fieldNameD[j],startpart+i,isize,test); 

              isize = computenitems(i,j,myrank,fieldNameD[j],paraD,expectD[j],numVariables[i]);
              //printf("fieldname: %s part: %d isize: %d\n",fieldNameD[j],startpart+i,isize); 

	      Dfield[i][j] = new double[isize]; // third rank is nVert*nvars
	      readdatablock( &igeom,
			      fieldNameD[j],
			      (void*)Dfield[i][j],
			      &isize,
			      "double",
			      "binary" );
               if(j==0) {
                 CoordsOld[i]=new double[isize];
                 for(k=0; k<isize; k++) CoordsOld[i][k]=Dfield[i][j][k];
               } 
	    }
	}
//      MPI_Barrier(MPI_COMM_WORLD);
      closefile(&igeom, "read");
    }

  if(myrank==0){
    printf("Starting to read coordinates (doubles) in the NEW geombc.dat.## files\n");
  }

  for ( i = 0; i < nppp; i++ )
    {
      if(existsubdir) { 
        subdir = (startpart+i-1) / DIR_FANOUT;
        sprintf(gfname,"./%d-procs_case-1PPP-GM/%d/geombc.dat.%d",N_parts, subdir, startpart+i);
      }
      else
        sprintf(gfname,"./%d-procs_case-1PPP-GM/geombc.dat.%d",N_parts,startpart+i);

      openfile(gfname,"read",&igeom);

         // hack that exploits the fact that j=1 from old data matches j=1 from new since it is coord
          j=0;


	  for ( k = 0; k < 10; k++ )
	    iarray[k]=0;

	  iarray[0]=-1;
	  WriteLockD[j]=0;
	  readheader( &igeom,
		       fieldNameD[j],
		       (void*)iarray,
		       &expectD[j],
		       "double",
		       "binary" );

	  if ( cscompare("block",headerTypeD[j]) )
	    {

              isize = computenitems(i,j,myrank,fieldNameD[j],paraD,expectD[j],numVariables[i]);
              //printf("fieldname: %s part: %d isize: %d\n",fieldNameD[j],startpart+i,isize); 

	      CoordsNew[i] = new double[isize];
	      readdatablock( &igeom,   // this was changed to point to new direccory and should get new coords
			      fieldNameD[j],
			      (void*)CoordsNew[i],
			      &isize,
			      "double",
			      "binary" );
	    }
	}
//      MPI_Barrier(MPI_COMM_WORLD);
      closefile(&igeom, "read");





  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("Starting to read some blocks (integers) in the geombc.dat.## files\n");
  }

  // Count the number of interior and boundary tpblocks for 2 new headers named
  // 'total number of different interior tpblocks' and
  // 'total number of different boundary tpblocks'
  int interiorCounter, boundaryCounter;
  interiorCounter=0;  
  boundaryCounter=0;  
  for ( j = 0; j < N_geombc_integer; j++ )
    {
       if (cscompare("connectivity interior",fieldNameI[j]))
	{
                  //printf("part: %d, fieldNameI[j]: %s\n",GPID,fieldNameI[j]);
		  interiorCounter++;
        }
        else if (cscompare("connectivity boundary",fieldNameI[j]))
	{
                  //printf("part: %d, fieldNameI[j]: %s\n",GPID,fieldNameI[j]);
		  boundaryCounter++;
        }
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("There are %d total connectivity interior and %d total connectivity boundary\n", interiorCounter, boundaryCounter);
  }
 


  // Now, start to read the integer fields
  for ( i = 0; i < nppp; i++ )
    {
      if(existsubdir) { 
        subdir = (startpart+i-1) / DIR_FANOUT;
        sprintf(gfname,"./%d-procs_case-1PPP-BM/%d/geombc.dat.%d",N_parts, subdir, startpart+i);
      }
      else
        sprintf(gfname,"./%d-procs_case-1PPP-BM/geombc.dat.%d", N_parts, startpart+i);

      openfile(gfname,"read",&igeom);

//      MPI_Barrier(MPI_COMM_WORLD);

//      printf("gfname is %s and nppp is %d myrank %d\n",gfname,i,myrank);

      for ( j = 0; j < N_geombc_integer; j++ )
	{

//	  printf("Writing integer ... %s\n",fieldNameI[j]);
	  paraI[i][j] = new int[expectI[j]];

	  for ( k = 0; k < 10; k++ )
	    iarray[k]=0;

	  //	  printf("myrank %d and i %d j %d numBou is %d\n",myrank,i,j,numBoundaryFields[i]);

	  WriteLockI[j]=0;
	  iarray[0]=-1;

          if ( cscompare("total number of interior tpblocks",fieldNameI[j] ) )
          { 
             iarray[0] = interiorCounter; //New header that does not exist in the posix file 
          }
          else if ( cscompare("total number of boundary tpblocks",fieldNameI[j] ) )
          { 
             iarray[0] = boundaryCounter; //New header that does not exist in the posix file 
          }
          else
          {
	     readheader( &igeom,
		          fieldNameI[j],
		          (void*)iarray,
		          &expectI[j],
		          "integer",
		          "binary" );
	     if ( iarray[0]==-1)
	       WriteLockI[j]=1; // The field was not found in the posix geombc file
          }

	  //MPI_Barrier(MPI_COMM_WORLD);

/*          int iproc;
          for(iproc=0; iproc<N_procs; iproc++){
            MPI_Barrier(MPI_COMM_WORLD);
            if(iproc == myrank){
              printf(" iproc: %d ", iproc);
              printf("part: %d myrank: %d field: %s header: ",startpart+i,myrank,fieldNameI[j]);
              for ( k = 0; k < expectI[j]; k++ )
                printf(" %d ",iarray[k]);
              printf("\n");
            }
            usleep(100);
            MPI_Barrier(MPI_COMM_WORLD);
          }*/

	  for ( k = 0; k < expectI[j]; k++ )
	    paraI[i][j][k] = iarray[k];

	  if ( WriteLockI[j]==1) //The field is not present but SyncIO needs it to read collectively. Put 0.
	    for ( k = 0; k < expectI[j]; k++ )
	      paraI[i][j][k] = 0;

	  if ( cscompare("block",headerTypeI[j]) )
	    {
/*	      if ( expectI[j]==1)
		isize = paraI[i][j][0];
	      else
		isize = paraI[i][j][0] * paraI[i][j][1];

	      if (cscompare("nbc codes",fieldNameI[j]))
		isize = paraI[i][j][0] * 2;
*/
//              int test;
//              test = computenitems(i,j,myrank,fieldNameI[j],paraI,expectI[j],numVariables[i]);
//              printf("irank: %d fieldname: %s ParaI: %d\n",myrank,fieldNameI[j], parapI[0][0][0], numVariables[i]);
//              if(test != isize)
//                printf("PROBLEM fieldname: %s part: %d isize: %d test: %d\n",fieldNameI[j],startpart+i,isize,test); 
//              else 
//                printf("fieldname: %s part: %d isize: %d test: %d\n",fieldNameI[j],startpart+i,isize,test); 

              isize = computenitems(i,j,myrank,fieldNameI[j],paraI,expectI[j],numVariables[i]);
              //printf("fieldname: %s part: %d isize: %d\n",fieldNameI[j],startpart+i,isize); 

	      Ifield[i][j] = new int[isize];
	      readdatablock( &igeom,
			      fieldNameI[j],
			      (void*)Ifield[i][j],
			      &isize,
			      "integer",
			      "binary" );
	    }
	}
//      MPI_Barrier(MPI_COMM_WORLD);
      closefile(&igeom, "read");
    }

  MPI_Barrier(MPI_COMM_WORLD); //added by MR



  ///////////////////// Writing geometry SYNCIO Removed but need a few declarations ///////

  int nppf = N_parts/N_files;
  int N_geombc = N_geombc_double + N_geombc_integer;
  int writeHandle, GPID;
  char fname[255],fieldtag[255];
  char dirname[1024];


  bzero((void*)fname,255);

  /////////////////////// restart data ////////////////////////////

  int irestart;

  Dfield = new double**[N_restart_double];  // inexplicably, restarts transpose first two indexes...first is restartDoubles??
  FixedDfield = new double**[N_restart_double];  
  Ifield = new int**[N_restart_integer];

  paraD = new int**[N_restart_double];
  paraI = new int**[N_restart_integer];

  expectD = new int[N_restart_double];
  expectI = new int[N_restart_integer];

  fieldNameD = new char*[N_restart_double];
  fileTypeD = new char*[N_restart_double];
  dataTypeD = new char*[N_restart_double];
  headerTypeD = new char*[N_restart_double];

  fieldNameI = new char*[N_restart_integer];
  fileTypeI = new char*[N_restart_integer];
  dataTypeI = new char*[N_restart_integer];
  headerTypeI = new char*[N_restart_integer];

  if (N_restart_double>0)
    for ( i = 0; i < N_restart_double; i++ )
      {
	WriteLockD[i]=0;
	Dfield[i] = new double*[nppp]; // second rank is parts per process
	FixedDfield[i] = new double*[nppp]; 

	paraD[i] = new int*[nppp];

	fieldNameD[i] = new char[128];
	fileTypeD[i] = new char[128];
	dataTypeD[i] = new char[128];
	headerTypeD[i] = new char[128];
      }

  if (N_restart_integer>0)
    for ( i = 0; i < N_restart_integer; i++ )
      {
	WriteLockI[i]=0;
	Ifield[i] = new int*[nppp];

	paraI[i] = new int*[nppp];

	fieldNameI[i] = new char[128];
	fileTypeI[i] = new char[128];
	dataTypeI[i] = new char[128];
	headerTypeI[i] = new char[128];
      }

  ////////////////////////////////////////////////////////////////
  // temporary fix: in the new version the double and integer
  //                can mix and match to avoid the order confusion
  ////////////////////////////////////////////////////////////////

  double_counter=0,int_counter=0;

  for ( i = 0; i < N_restart_double+N_restart_integer; i++ )
    {
      fgets( target, 1024, pFile );
      temp = strtok( target, ";" );
      token = strtok( temp, "," );
      strcpy( S1, token );
      token = strtok ( NULL, "," );
      strcpy( S2, token );
      token = strtok ( NULL, "," );
      strcpy( S3, token );
      token = strtok ( NULL, "," );
      strcpy( S4, token );
      token = strtok ( NULL, "," );
      strcpy( S5, token );

      if (cscompare(S3,"double"))
	{
	  strcpy( fileTypeD[double_counter], S1 );
	  strcpy( fieldNameD[double_counter], S2 );
	  strcpy( dataTypeD[double_counter], S3 );
	  strcpy( headerTypeD[double_counter], S4 );
	  strcpy( numTemp, S5 );
	  expectD[double_counter] = atoi (numTemp);
	  double_counter++;
	}

      if (cscompare(S3,"integer"))
	{
	  strcpy( fileTypeI[int_counter], S1 );
	  strcpy( fieldNameI[int_counter], S2 );
	  strcpy( dataTypeI[int_counter], S3 );
	  strcpy( headerTypeI[int_counter], S4 );
	  strcpy( numTemp, S5 );
	  expectI[int_counter] = atoi (numTemp);
	  int_counter++;
	}
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("Starting to read some blocks (doubles) in the restart.dat.##.## files\n");
  }

  for ( i = 0; i < N_restart_double; i++ )
    {
      for ( j = 0; j < nppp; j++ )
	{

          if(existsubdir) { 
            subdir = (startpart+j-1) / DIR_FANOUT;
	    sprintf(gfname,"./%d-procs_case-1PPP-BM/%d/restart.%d.%d",N_parts, subdir, N_steps,startpart+j);
          }
          else
	    sprintf(gfname,"./%d-procs_case-1PPP-BM/restart.%d.%d",N_parts, N_steps, startpart+j);

	  openfile(gfname,"read",&irestart);

	  for ( k = 0; k < 10; k++ )
	    iarray[k]=0;

	  paraD[i][j] = new int[expectD[i]];

	  iarray[0]=-1;
	  readheader( &irestart,
		       fieldNameD[i],
		       (void*)iarray,
		       &expectD[i],
		       "double",
		       "binary" );

	  for ( k = 0; k < expectD[i]; k++ )
	    paraD[i][j][k] = iarray[k];

	  if ( iarray[0]==-1 )
	      WriteLockD[i]=1;
	  if ( WriteLockD[i]==0 )
	    {
	      if ( cscompare("block",headerTypeD[i]) )
		{
		  if ( expectD[i]==1)
		    isize = paraD[i][j][0];
		  else
		    isize = paraD[i][j][0] * paraD[i][j][1];

		  Dfield[i][j] = new double[isize]; // third rank is nVerts*nvars_per_vert
		  FixedDfield[i][j] = new double[isize]; 
		  readdatablock( &irestart,
				  fieldNameD[i],
				  (void*)Dfield[i][j],
				  &isize,
				  "double",
				  "binary" );

		}
	    }
	  closefile(&irestart, "read");
	}
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("Starting to read some blocks (integers) in the restart.dat.##.## files\n");
  }


  for ( i = 0; i < N_restart_integer; i++ )
    {
      for ( j = 0; j < nppp; j++ )
	{

          if(existsubdir) { 
            subdir = (startpart+j-1) / DIR_FANOUT;
	    sprintf(gfname,"./%d-procs_case-1PPP-BM/%d/restart.%d.%d",N_parts, subdir, N_steps, startpart+j);
          }
          else
	    sprintf(gfname,"./%d-procs_case-1PPP-BM/restart.%d.%d",N_parts, N_steps, startpart+j);

	  openfile(gfname,"read",&irestart);

	  for ( k = 0; k < 10; k++ )
	    iarray[k]=0;

	  paraI[i][j] = new int[expectI[i]];

	  iarray[0]=-1;
	  readheader( &irestart,
		       fieldNameI[i],
		       (void*)iarray,
		       &expectI[i],
		       "integer",
		       "binary" );

	  for ( k = 0; k < expectI[i]; k++ )
	    paraI[i][j][k] = iarray[k];

	  if ( iarray[0]==-1 )
	      WriteLockI[i]=1;
	  if ( WriteLockI[i]==0 )
	    {

	      if ( cscompare("block",headerTypeI[i]) )
		{
		  if ( expectI[i]==1)
		    isize = paraI[i][j][0];
		  else
		    isize = paraI[i][j][0] * paraI[i][j][1];

		  Ifield[i][j] = new int[isize];
		  readdatablock( &irestart,
				  fieldNameI[i],
				  (void*)Ifield[i][j],
				  &isize,
				  "integer",
				  "binary" );
		}
	    }
	  closefile(&irestart, "read");
	}
    }

  fclose(pFile);
  ///////////////////// Fixing ///////////////////////////////

  double etol=1e-7;
  for ( i = 0; i < nppp; i++ ) {
//DEBUG     bzero((void*)fname,255);
//DEBUG     sprintf(fname,"./%d-procs_case-SyncIO-%d-FM/fixlog.%d.%d",N_parts,N_files,myrank,i);
//DEBUG     pFile = fopen(fname,"w");
//DEBUG     fprintf(pFile,"start\n");
//DEBUG     MPI_Barrier(MPI_COMM_WORLD);

    j=0;
//    if ( expectD[i]==1)
//      isize = paraD[i][j][0];
//    else
//      isize = paraD[i][j][0] * paraD[i][j][1];
//FAIL    isize = computenitems(i,j,myrank,fieldNameD[j],paraD,expectD[j],numVariables[i]);
    int nverts=paraD[0][i][0]; //FAIL isize/numVariables[i];  // probably paraD[i][0]
    int k,k2,kf,kd,ipm,ipm2;
    int iskip=nverts;
    int iskip2=2*nverts;
    int ifound;
    double xO,yO,zO,xN,yN,zN,pO,pN;
//DEBUG    fprintf(pFile,"start loop k over %d verts \n ",nverts);
//DEBUG    MPI_Barrier(MPI_COMM_WORLD);
    for ( k = 0; k < nverts; k++ ) { // k counts old coordinate
      ipm=k; // pointer to the x for kth coordinate
      xO=CoordsOld[i][ipm];    
      yO=CoordsOld[i][(ipm+iskip)]; 
      zO=CoordsOld[i][(ipm+iskip2)]; 
//DEBUG      fprintf(pFile,"start loop looking for %f %f %f \n ",xO,yO,zO);
      for ( k2 = 0; k2 < nverts; k2++ ) { // k2 counts new coordinates
        ipm2=k2; // pointer
        xN=CoordsNew[i][ipm2];
        yN=CoordsNew[i][(ipm2+iskip)];
        zN=CoordsNew[i][(ipm2+iskip2)];
//        if(abs(CoordsOld[i][ipm]           - CoordsNew[i][ipm2]          ) <etol &&
//           abs(CoordsOld[i][(ipm+iskip)]   - CoordsNew[i][(ipm2+iskip)] ) <etol &&
//           abs(CoordsOld[i][(ipm+iskip2)] - CoordsNew[i][(ipm2+iskip2)]) < etol ) {
        ifound=0;
        if(abs(xO -xN) <etol &&
           abs(yO -yN) < etol &&
           abs(zO -zN) < etol ) {
//DEBUG           fprintf(pFile,"%d %d \n",k,k2);
           ifound=1;
           if(k==nverts/10 && myrank==0 ) //  && i==0)
              printf("10 percent complete of first ppp: k=%d k2=%d\n", k, k2);
           for ( kf=0; kf < N_restart_double; kf++) {  // kf counts the double fields being transferred
             int nvarsIR=paraD[kf][i][1];
//DEBUG             fprintf(pFile,"field number %d attempting to map %d variables \n",kf+1,nvarsIR);
             for ( kd=0; kd < nvarsIR; kd++)  {  // note paraD AND Dfield areflipped relative to Geom??
               FixedDfield[kf][i][k2+kd*nverts]= Dfield[kf][i][k+kd*nverts];
//fordeug               pN=FixedDfield[kf][i][k2+kd*nverts]= Dfield[kf][i][k+kd*nverts];
//fordeug               pO=Dfield[kf][i][k+kd*nverts];
             }
           }
           break;
        }
      }
//DEBUG      if(ifound==0) fprintf(pFile," k=%d no match \n",k);
      if(ifound==0) printf(" k=%d no match \n",k);
    }
//DEBUG    fprintf(pFile,"all %d verts found for this part \n",nverts);
//DEBUG    fclose(pFile);
  }       

  if(myrank==0 ) //  && i==0)
     printf("finished fixing \n");

  ///////////////////// Writing SyncIO///////////////////////////////

  int N_restart = N_restart_double + N_restart_integer;

  bzero((void*)fname,255);
  sprintf(fname,"./%d-procs_case-SyncIO-%d-FM/restart-dat.%d.%d",N_parts,N_files,N_steps,((int)(myrank/(N_procs/N_files))+1));
  initphmpiio(&N_restart, &nppf, &N_files,&writeHandle, "write");
  openfile(fname, "write", &writeHandle);

  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("Starting to write some blocks (doubles) in the restart-dat.##.## files\n");
  }

  for ( i = 0; i < N_restart_double; i++ )
    {
      for (  j = 0; j < nppp; j++  )
	{

	  if (WriteLockD[i]==0)
	    {
	      GPID = startpart + j;
	      bzero((void*)fieldtag,255);
	      sprintf(fieldtag,"%s@%d",fieldNameD[i],GPID);

	      if ( expectD[i]==1)
		isize = paraD[i][j][0];
	      else
		isize = paraD[i][j][0] * paraD[i][j][1];

	      for ( k = 0; k < expectD[i]; k++ )
		iarray[k] = paraD[i][j][k];

	      if ( cscompare("header",headerTypeD[i]) )
		isize = 0;

	      writeheader( &writeHandle,
			    fieldtag,
			    (void*)iarray,
			    &expectD[i],
			    &isize,
			    "double",
			    "binary");

	      writedatablock( &writeHandle,
			       fieldtag,
			       (void*)FixedDfield[i][j],
			       &isize,
			       "double",
			       "binary" );
	      if ( cscompare("block",headerTypeD[i]) )
// wait... will write posix too		delete [] FixedDfield[i][j];
		delete [] Dfield[i][j];
	    }
// wait... will write posix too	  delete [] paraD[i][j];
	}
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("Starting to write some blocks (integers) in the restart.dat.##.## files\n");
  }

  for ( i = 0; i < N_restart_integer; i++ )
    {
      for (  j = 0; j < nppp; j++  )
	{

	  if (WriteLockI[i]==0)
	    {
	      GPID = startpart + j;
	      bzero((void*)fieldtag,255);
	      sprintf(fieldtag,"%s@%d",fieldNameI[i],GPID);

	      if ( expectI[i]==1)
		isize = paraI[i][j][0];
	      else
		isize = paraI[i][j][0] * paraI[i][j][1];

	      for ( k = 0; k < expectI[i]; k++ )
		iarray[k] = paraI[i][j][k];

	      if ( cscompare("header",headerTypeI[i]) )
		isize = 0;

	      writeheader( &writeHandle,
			    fieldtag,
			    (void*)iarray,
			    &expectI[i],
			    &isize,
			    "integer",
			    "binary");

	      writedatablock( &writeHandle,
			       fieldtag,
			       (void*)Ifield[i][j],
			       &isize,
			       "integer",
			       "binary" );

// wait... will write posix too	      if ( cscompare("block",headerTypeI[i]) )
// wait... will write posix too		delete [] Ifield[i][j];
	    }
// wait... will write posix too	  delete [] paraI[i][j];
	}

    }

  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("Closing restart-dat.##.## files\n");
  }

  closefile(&writeHandle, "write");
  finalizephmpiio(&writeHandle);



  //////////////////////////writing Fixed Posix////////////////////////////

  int irstou;
  int magic_number = 362436;
  int* mptr = &magic_number;
  int nitems = 1;
  int ret;

  bzero((void*)fname,255);
  sprintf(fname,"./%d-procs_case-1PPP-FM",N_parts);
  if(myrank == 0) {
    ret=mkdir(fname,0777);
    if(ret<0) {
      if(errno == EEXIST) {
       // acceptable
      } else {
        printf("ERROR - Could not create procs_case-1PPP-FM directory\n");
        return 1;
      }
    }
    printf("mkdir 1PPP succeeded on %d \n", myrank);

  }
  MPI_Barrier(MPI_COMM_WORLD);
  ret=chdir(fname);
  if(ret<0) {
      printf("ERROR - Could not chdir procs_case-1PPP-FM directory\n");
      return 1;
  }
  if(myrank == 0) printf("chdir to 1PPP succeeded on %d \n", myrank);

  // newer versions of phastaIO do the fanout for posix but assume that calling program is already in the 
  // posix part_case directory. 

  if(myrank == 0) {
    bzero((void*)dirname,sizeof(dirname));
    getcwd(dirname, sizeof(dirname));
    printf("in a fanout created in phastaIO underneath %s \n", dirname);
  }

/* Doubtful that we need to come back to top dir but this would do it.
    ret=chdir("..");
      if(ret<0) {
          printf("ERROR - Could not chdir back to cwd  \n"); 
          return 1; 
      }
      bzero((void*)dirname,sizeof(dirname));
      getcwd(dirname, sizeof(dirname));
      printf("current working director is %s \n", fname);
*/
  for ( j = 0; j < nppp; j++ ) {
  bzero((void*)fname,255);

  int mypart=1+myrank*nppp+j;  // or  1+myrank+j*N_parts
  sprintf(fname,"restart.%d.%d",N_steps,mypart);
  openfile(fname,"write", &irstou);

  /* writing the top ascii header for the restart file */

  writestring( &irstou,"# PHASTA Input File Version 2.0\n");
  writestring( &irstou,
                "# format \"keyphrase : sizeofnextblock usual headers\"\n");

  bzero( (void*)fname, 255 );
  writestring( &irstou, fname );

  writestring( &irstou, fname );
  writestring( &irstou,"\n");


  isize = 1;
  nitems = 1;
  iarray[ 0 ] = 1;
  writeheader( &irstou, "byteorder magic number ",
                (void*)iarray, &nitems, &isize, "integer", "binary" );

  nitems = 1;
  writedatablock( &irstou, "byteorder magic number ",
                   (void*)mptr, &nitems, "integer", "binary" );

  for ( i = 0; i < N_restart_double; i++ )
    {
//wrong      for ( j = 0; j < nppp; j++ )
//        {
          if ( WriteLockD[i] == 0 )
            {
              if ( cscompare("header",headerTypeD[i]) )
                {
                  bzero( (void*)fname, 255 );
                  sprintf(fname,"%s : < 0 > %d\n", fieldNameD[i],paraD[i][j][0]);                 
                  writestring( &irstou, fname );
                }

              if ( cscompare("block",headerTypeD[i]) )
                {
                  if ( expectD[i]==1 )
                    isize = paraD[i][j][0];
                  else
                    isize = paraD[i][j][0] * paraD[i][j][1];

                  for ( k = 0; k < expectD[i]; k++ )
                    iarray[k] = paraD[i][j][k];

                  if ( cscompare("header",headerTypeD[i]) )
                    isize = 0;

                  writeheader( &irstou,
                                fieldNameD[i],
                                (void*)iarray,
                                &expectD[i],
                                &isize,
                                "double",
                                "binary");
                  writedatablock( &irstou,
                                   fieldNameD[i],
                                   (void*)FixedDfield[i][j],
                                   &isize,
                                   "double",
                                   "binary");
                }

              
              if ( cscompare("block",headerTypeD[i]) )
                delete [] FixedDfield[i][j];
 //           } //closes wrong loop over parts per process....posix does NOT write multiple parts in a file
          delete [] paraD[i][j];
        }
    }


  for ( i = 0; i < N_restart_integer; i++ )
    {
//wrong      for ( j = 0; j < nppp; j++ )
//        {

          if ( WriteLockI[i] == 0 )
            {

              if ( cscompare("header",headerTypeI[i]) )
                {
                  bzero( (void*)fname, 255 );
                  sprintf(fname,"%s : < 0 > %d\n", fieldNameI[i],paraI[i][j][0]);                 
                  writestring( &irstou, fname );
                }
    
              if ( cscompare("block",headerTypeI[i]) )
                {
                  if ( expectI[i]==1 )
                    isize = paraI[i][j][0];
                  else
                    isize = paraI[i][j][0] * paraI[i][j][1];
    
                  for ( k = 0; k < expectI[i]; k++ )
                    iarray[k] = paraI[i][j][k];
    
                  writeheader( &irstou,
                                fieldNameI[i],
                                (void*)iarray,
                                &expectI[i],
                                &isize,
                                "integer",
                                "binary");
                  writedatablock( &irstou,
                                   fieldNameI[i],
                                   (void*)Ifield[i][j],
                                   &isize,
                                   "integer",
                                   "binary");
                }
              
              if ( cscompare("block",headerTypeI[i]) )
                delete [] Ifield[i][j];
//            } //closes wrong loop over parts per process....posix does NOT write multiple parts in a file
          delete [] paraI[i][j];
        }
    }


  closefile( &irstou, "write" );
  } // proper closing of the parts per process loop...a separate file for each part


  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    printf("Free memory related to restart-dat.##.## files\n");
  }


  for ( i = 0; i < N_restart_double; i++ )
    {
      delete [] FixedDfield[i];
      delete [] Dfield[i];
      delete [] paraD[i];

      delete [] fieldNameD[i];
      delete [] fileTypeD[i];
      delete [] dataTypeD[i];
      delete [] headerTypeD[i];
    }

  for ( i = 0; i < N_restart_integer; i++ )
    {
      delete [] Ifield[i];
      delete [] paraI[i];

      delete [] fieldNameI[i];
      delete [] fileTypeI[i];
      delete [] dataTypeI[i];
      delete [] headerTypeI[i];
    }

  delete [] Dfield;
  delete [] FixedDfield;
  delete [] Ifield;

  delete [] paraD;
  delete [] paraI;

  delete [] expectD;
  delete [] expectI;

  delete [] fieldNameD;
  delete [] fileTypeD;
  delete [] dataTypeD;
  delete [] headerTypeD;

  delete [] fieldNameI;
  delete [] fileTypeI;
  delete [] dataTypeI;
  delete [] headerTypeI;

  delete [] WriteLockD;
  delete [] WriteLockI;

  delete [] interiorMark;
  delete [] boundaryMark;
  delete [] codesMark;
  delete [] valuesMark;
  delete [] numVariables;

  if (myrank==0)
    {
      printf("\nFinished transfer, please check data using:\n");
      printf(" grep -a ': <' filename \n\n");
      printf("Note that the size of the fields is computed based on previous geombc files\n");
      printf("Check the routine 'computenitems' if you have any reason to think it has changes for the fields you are interested in\n\n");
    }

  MPI_Finalize();

}


