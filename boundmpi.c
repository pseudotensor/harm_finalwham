#include "decs.h"

#define DEBUG 0

// boundtype specifies whether to bound scalar or to bound vector that is only needed to be bound along that direction
int bound_mpi(int boundstage, int boundtype, FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR])
{
  FTYPE (*prim2bound)[N2M][N3M][NPR];
  int dir;
#if(DEBUG)
  int i,j,k,pl;
#endif

#if(USEMPI)
  /* These arrays contain designations that identify 
   * each recv and send */
  static MPI_Request requests[COMPDIM * 2 * 2];
  // format of map for requests[dir*2+recv/send(0/1)]
#endif

#if(USEMPI)

  /*
   *
   * 1. do outflow/inflow/etc. boundary conditions (in bounds)
   * 2. pack data into workbc arrays and do send/recv's
   * 3. for each transfer do a wait/unpack
   * 4. do all send waits
   *
   * NOTE: order is important so corner zones are assigned correctly
   *
   * workbc[PACK][][] is data that is being sent
   * workbc[UNPACK][][] is data that is being recvd
   *
   */

  /*
   *  -> higher values of x1
   * 2 -> lower values of x1
   *
   */

  
  /* bounds has already done non-MPI boundary conditions; now do MPI
   * boundary conditions */

  // must go in this order (0,1) then (2,3) then (4,5) or visa versa,
  // a single directions order doesn't matter.  This method of
  // ordering is as opposed to directly transfering the corner zones
  // to the corner CPU.

  // This may be faster since all transfers can proceed at once,
  // although this may be slower since no transfers can occur until
  // packing is completed.  This way packing and transfering occur
  // simultaneously.

  // Although l/r are packed together, since in the end we have to
  // wait for both l/r to complete, so equal time completion is
  // favored //over asynch completion.

  // Also, transfering corner zones with small message sizes increases
  // the importance of latency.

  // for 2D:
  // I choose left-right N1M/N2M first, then up/down N1M/N2M.  Could
  // just do N1/N2 for interior for L/R, but true boundary needs full
  // N1M/N2M exchanged since cpu sets boundary using normal bc code
  // which needs to get transfered to neight(i.e. currently if corner
  // has bctype 99/?  then doesn't do corner)

  // GODMARK: Make sure 3D makes sense (no extra things to do)

#if(DEBUG)
  PBOUNDLOOP(pl){
    FULLLOOP{
      prim[i][j][k][pl]=-1-pl*100; // should be turned into myid-k*100
    }
    ZLOOP {
      prim[i][j][k][pl]=myid-pl*100; // differentiates but clear per pr
      //      fprintf(log_file,"%d %d %d %d %21.15g\n",i,j,k,pl,prim[i][j][k][pl]);
    }
  }
#endif


  // this is designed to copy corners by indirectly copying them.  The
  // corners eventually get to corner-related CPUs by passing through
  // another cpu.  This required ordering the left/right and up/down
  // procedures.

  // one could copy the corners directly and get more bandwidth since
  // would transfers 2X as much data, but corners would transfer very
  // slowly alone, and basically has the same number of operations
  // required as does edge transfers.
  //  fprintf(fail_file,"innerboundhere1\n"); fflush(fail_file);

  ///////////////////
  //
  // x -dir
  //
  /////////////////
  if(boundtype==BOUNDPRIMTYPE) prim2bound=prim;
  else if(boundtype==BOUNDFLUXTYPE) prim2bound=F1;

  if((boundstage==STAGEM1)||(boundstage==STAGE0)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundtype][dir][DIRIF]) pack(dir,boundtype,prim2bound,workbc);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundtype][dir][DIRIF]) sendrecv(dir,boundtype,workbc,requests);
  }
  //fprintf(fail_file,"innerboundhere2\n"); fflush(fail_file);
  if((boundstage==STAGEM1)||(boundstage==STAGE1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundtype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundtype][dir][DIRIF]) unpack(dir,boundtype,workbc,prim2bound);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundtype][dir][DIRIF]) sendwait(dir,requests);
  }
  ///////////////////
  //
  // y -dir
  //
  /////////////////
  //fprintf(fail_file,"innerboundhere3\n"); fflush(fail_file);
  if(boundtype==BOUNDPRIMTYPE) prim2bound=prim;
  else if(boundtype==BOUNDFLUXTYPE) prim2bound=F2;

  if((boundstage==STAGEM1)||(boundstage==STAGE2)){
    // now dir=0,1(X1UP,X1DN) is done, so can start 2,3(X2UP,X2DN)
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundtype][dir][DIRIF]) pack(dir,boundtype,prim2bound,workbc);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundtype][dir][DIRIF]) sendrecv(dir,boundtype,workbc,requests);
  }
  //fprintf(fail_file,"innerboundhere4\n"); fflush(fail_file);
  if((boundstage==STAGEM1)||(boundstage==STAGE3)){
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundtype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundtype][dir][DIRIF]) unpack(dir,boundtype,workbc,prim2bound);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundtype][dir][DIRIF]) sendwait(dir,requests);
  }
  //fprintf(fail_file,"innerboundhere5\n"); fflush(fail_file);

  ///////////////////
  //
  // z -dir
  //
  /////////////////
  //fprintf(fail_file,"innerboundhere3\n"); fflush(fail_file);
  if(boundtype==BOUNDPRIMTYPE) prim2bound=prim;
  else if(boundtype==BOUNDFLUXTYPE) prim2bound=F3;

  if((boundstage==STAGEM1)||(boundstage==STAGE4)){
    // now dir=0,1,2,3(X1UP,X1DN,X2UP,X2DN) is done, so can start 4,5(X3UP,X3DN)
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundtype][dir][DIRIF]) pack(dir,boundtype,prim2bound,workbc);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundtype][dir][DIRIF]) sendrecv(dir,boundtype,workbc,requests);
  }
  //fprintf(fail_file,"innerboundhere4\n"); fflush(fail_file);
  if((boundstage==STAGEM1)||(boundstage==STAGE5)){
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundtype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundtype][dir][DIRIF]) unpack(dir,boundtype,workbc,prim2bound);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundtype][dir][DIRIF]) sendwait(dir,requests);
  }
  //fprintf(fail_file,"innerboundhere5\n"); fflush(fail_file);

#if(DEBUG)
    fprintf(log_file,"\n\nafter\n\n");
    FULLLOOP{
      PBOUNDLOOP(pl){
        fprintf(log_file,"%d %d %d %d %21.15g\n",i,j,k,pl,prim[i][j][k][pl]);
      }
    }
    myexit(0);
#endif



  // end if mpi
#endif


    return(0);

}	
// end function


// PACKLOOP allows one to alter which i or j or k is faster iterated
#define PACKLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,pr,num) GENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop) PLOOPMPI(pr,num)

// packs data for shipment
void pack(int dir, int boundtype, FTYPE prim[][N2M][N3M][NPR],FTYPE workbc[][COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM])
{
  // dir=direction sending
  int i,j,k;
  int bci,pr;

  bci=0;
  PACKLOOP(i,j,k
	   ,dirset[boundtype][dir][DIRPSTART1]
	   ,dirset[boundtype][dir][DIRPSTOP1]
	   ,dirset[boundtype][dir][DIRPSTART2]
	   ,dirset[boundtype][dir][DIRPSTOP2]
	   ,dirset[boundtype][dir][DIRPSTART3]
	   ,dirset[boundtype][dir][DIRPSTOP3]
	   ,pr
	   ,dirset[boundtype][dir][DIRNUMPR]
	   ){
    /*
    if(bci>=dirset[boundtype][dir][DIRSIZE]){
      dualfprintf(fail_file,"pack memory leak: bci: %d dirset[%d][DIRSIZE]: %d\n",bci,dirset[boundtype][dir][DIRSIZE]);
      myexit(10);
    }
    */
    workbc[PACK][dir][bci++] = prim[i][j][k][pr];
  }
}

#if(USEMPI)
void sendrecv(int dir, int boundtype, FTYPE workbc[][COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM],MPI_Request *requests)
{
  MPI_Irecv(workbc[UNPACK][dir],
	    dirset[boundtype][dir][DIRSIZE],
	    MPI_FTYPE,
	    dirset[boundtype][dir][DIROTHER],
	    dirset[boundtype][dir][DIRTAGR],
	    MPI_COMM_WORLD,
	    &requests[dir*2+REQRECV]);

  MPI_Isend(workbc[PACK][dir],
	    dirset[boundtype][dir][DIRSIZE],
	    MPI_FTYPE,
	    dirset[boundtype][dir][DIROTHER],
	    dirset[boundtype][dir][DIRTAGS],
	    MPI_COMM_WORLD,
	    &requests[dir*2+REQSEND]);

}
#endif


#if(USEMPI)
void recvwait(int dir,MPI_Request *requests)
{
  MPI_Wait(&requests[dir*2+REQRECV], &mpichstatus);

}
#endif


void unpack(int dir, int boundtype, FTYPE workbc[][COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM],FTYPE prim[][N2M][N3M][NPR])
{
  // dir is direction receiving from
  int i,j,k;
  int bci,pr;

  bci=0;
  PACKLOOP(i,j,k
	   ,dirset[boundtype][dir][DIRUSTART1]
	   ,dirset[boundtype][dir][DIRUSTOP1]
	   ,dirset[boundtype][dir][DIRUSTART2]
	   ,dirset[boundtype][dir][DIRUSTOP2]
	   ,dirset[boundtype][dir][DIRUSTART3]
	   ,dirset[boundtype][dir][DIRUSTOP3]
	   ,pr
	   ,dirset[boundtype][dir][DIRNUMPR]
	   ){
    prim[i][j][k][pr]=workbc[UNPACK][dir][bci++];
  }
}

#if(USEMPI)

void sendwait(int dir,MPI_Request *requests)
{
  MPI_Wait(&requests[dir*2+REQSEND], &mpichstatus);
}
#endif
