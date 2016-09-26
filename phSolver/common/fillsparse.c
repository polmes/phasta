#ifdef HAVE_PETSC
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>
#include <petsc.h>

#include "common_c.h"
#include "FCMangle.h"

#define fillsparsecpetscs FortranCInterface_GLOBAL_(fillsparsecpetscs, FILLSPARSECPETSCS)
#define fillsparsecpetscc FortranCInterface_GLOBAL_(fillsparsecpetscc, FILLSPARSECPETSCC)
#define fillsparsecpetsci FortranCInterface_GLOBAL_(fillsparsecpetsci, FILLSPARSECPETSCI)

#define COLMAJ2D(row,col,numrow) (row-1)+(col-1)*numrow
#define COLMAJ3D(a,b,c,amax,bmax,cmax) (a-1)+amax*((b-1)+bmax*(c-1))
#define ROWMAJ2D_ONE(row,col,numcol) (row-1)*numcol+(col-1)
typedef long long int gcorp_t;

void fillsparsecpetscs(gcorp_t* ieng, double* EGmass, Mat* lhsP)
{
        int npro = propar.npro;
        int nshl = shpdat.nshl;
        double* mb = (double*) malloc(sizeof(double)*nshl*nshl); //block to insert
        int e,i,j,aa; //following along with fillsparse.f
        PetscInt* locat = (PetscInt*) malloc(sizeof(PetscInt)*nshl);
        for(e=0;e<npro;e++)
        {
         for(aa=0;aa<nshl;aa++) locat[aa]=ieng[e+npro*aa]-1;
//         for(aa=0;aa<nshl;aa++) assert(locat[aa]>=0);
         for (i=0; i<nshl; i++)  {   // fill up Ke with respective egmass 
           for (j=0; j<nshl; j++)  {
            mb[nshl*i + j] = EGmass[e + npro*(i + nshl*j)];
           }
         }
         //MatSetValuesBlocked(*lhsP, nshl , locat, nshl, locat, mb, ADD_VALUES);
         PetscInt petsc_nshl;
         petsc_nshl = (PetscInt) nshl;
         MatSetValues(*lhsP, petsc_nshl , locat, petsc_nshl, locat, mb, ADD_VALUES);
        }
        free(mb);
	free(locat);
}
void fillsparsecpetscc(gcorp_t* ieng, double* EGmass, Mat* lhsP)
{
        int npro = propar.npro;
        int nshl = shpdat.nshl;
        int nedof = conpar.nedof;
        double* mb = (double*) malloc(sizeof(double)*nedof*nedof); //block to insert
        int e,i,j,aa; //following along with fillsparse.f
        //int* locat = (int*) malloc(sizeof(int)*nshl);
        PetscInt* locat = (PetscInt*) malloc(sizeof(PetscInt)*nshl);
        for(e=0;e<npro;e++)
        {
         for(aa=0;aa<nshl;aa++) locat[aa]=ieng[e+npro*aa]-1;
//         for(aa=0;aa<nshl;aa++) assert(locat[aa]>=0);
         for (i=0; i<nedof; i++)  {   /* fill up Ke with respective egmass */
           for (j=0; j<nedof; j++)  {
            mb[nedof*i + j] = EGmass[e + npro*(i + nedof*j)];
           }
         }
         //MatSetValuesBlocked(*lhsP, nshl , locat, nshl, locat, mb, ADD_VALUES);
         PetscInt petsc_nshl;
         petsc_nshl = (PetscInt) nshl;
         MatSetValuesBlocked(*lhsP, petsc_nshl , locat, petsc_nshl, locat, mb, ADD_VALUES);
        }
        free(mb);
	free(locat);
}
void fillsparsecpetsci(gcorp_t* ieng, double* xlhs, Mat* lhsP)
{
        int bsz = genpar.bsz;
        int npro = propar.npro;
        int nshl = shpdat.nshl;
        int nflow = conpar.nflow;
//        double* mbt = (double*) malloc(sizeof(double)*nflow*nflow); //sub-block to insert
        double* mb = (double*) malloc(sizeof(double)*nflow*nflow*nshl*nshl); //block to insert
        int nshlnfl,jhjmp,ihjmp,bsznf,id,e,iv,ih,jv,jh,nfsq,nfsqnsh,i00,i04; 
        //int* locat = (int*) malloc(sizeof(int)*nshl);
        PetscInt* locat = (PetscInt*) malloc(sizeof(PetscInt)*nshl);
        nfsq=nflow*nflow;
        nfsqnsh=nfsq*nshl;
        bsznf=bsz*nflow;
        nshlnfl=nflow*nshl;
        ihjmp=nfsq*bsz;
        jhjmp=ihjmp*nshl;
        for(e=0;e<npro;e++)
        {
         for(ih=0;ih<nshl;ih++) locat[ih]=ieng[e+npro*ih]-1;
//         for(aa=0;aa<nshl;aa++) assert(locat[aa]>=0);

         for(ih=0; ih<nshl;  ih++) {
           for(jh=0; jh<nshl;  jh++) {
             i00=e+ih*ihjmp+jh*jhjmp;
             for(iv=0; iv<4; iv++) {
               for(jv=0; jv<4; jv++) {
                 id=jv+iv*nshlnfl+jh*nflow+ih*nfsqnsh; //16*nshl;
                 mb[id]=xlhs[i00+jv*bsznf+iv*bsz]; // mbt[iv*4+jv];
               }
             }
           } 
         }
         //MatSetValuesBlocked(*lhsP, nshl , locat, nshl, locat, mb, ADD_VALUES);
         PetscInt petsc_nshl;
         petsc_nshl = (PetscInt) nshl;
         MatSetValuesBlocked(*lhsP, petsc_nshl , locat, petsc_nshl, locat, mb, ADD_VALUES);
        }
        free(mb);
	free(locat);
}
#endif

