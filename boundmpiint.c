#include "decs.h"


int bound_mpi_int(int boundstage, int prim[][N2M][N3M][NUMPFLAGS])
{
  int dir;

#if(USEMPI)
  /* These arrays contain designations that identify 
   * each recv and send */
  static MPI_Request requests[COMPDIM * 2 * 2];
  // format of map for requests[dir*2+recv/send(0/1)]
#endif

#if(USEMPI)

  
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) pack_int(dir,prim,workbc_int);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) sendrecv_int(dir,workbc_int,requests);
  }
  if((boundstage==STAGE1)||(boundstage==STAGEM1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) unpack_int(dir,workbc_int,prim);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) sendwait(dir,requests);
  }

  if((boundstage==STAGE2)||(boundstage==STAGEM1)){
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) pack_int(dir,prim,workbc_int);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) sendrecv_int(dir,workbc_int,requests);
  }
  if((boundstage==STAGE3)||(boundstage==STAGEM1)){
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) unpack_int(dir,workbc_int,prim);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) sendwait(dir,requests);
  }

  if((boundstage==STAGE4)||(boundstage==STAGEM1)){
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) pack_int(dir,prim,workbc_int);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) sendrecv_int(dir,workbc_int,requests);
  }
  if((boundstage==STAGE5)||(boundstage==STAGEM1)){
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) unpack_int(dir,workbc_int,prim);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[BOUNDPRIMTYPE][dir][DIRIF]) sendwait(dir,requests);
  }



  // end if mpi
#endif

  return(0);

}	
// end function


// PACKLOOP allows one to alter which i or j is faster iterated
#define PACKLOOP_INT(i,j,k,istart,istop,jstart,jstop,kstart,kstop) GENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop) FLOOP(pl)

// packs data for shipment
void pack_int(int dir,int prim[][N2M][N3M][NUMPFLAGS],int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM])
{
  // dir=direction sending
  int i,j,k;
  int pl;
  int bci;


  bci=0;
  PACKLOOP_INT(i,j,k
	   ,dirset[BOUNDPRIMTYPE][dir][DIRPSTART1]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRPSTOP1]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRPSTART2]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRPSTOP2]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRPSTART3]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRPSTOP3]
	       ){
    /*
    if(bci>=dirset[BOUNDPRIMTYPE][dir][DIRSIZE]){
      dualfprintf(fail_file,"pack memory leak: bci: %d dirset[%d][DIRSIZE]: %d\n",bci,dirset[BOUNDPRIMTYPE][dir][DIRSIZE]);
      myexit(10);
    }
    */
    workbc_int[PACK][dir][bci++] = prim[i][j][k][pl];
  }
}

#if(USEMPI)
void sendrecv_int(int dir,int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM],MPI_Request *requests)
{
  MPI_Irecv(workbc_int[UNPACK][dir],
	    dirset[BOUNDPRIMTYPE][dir][DIRSIZE]/NPR, // GODMARK: what's /NPR for?
	    MPI_INT,
	    dirset[BOUNDPRIMTYPE][dir][DIROTHER],
	    dirset[BOUNDPRIMTYPE][dir][DIRTAGR],
	    MPI_COMM_WORLD,
	    &requests[dir*2+REQRECV]);

  MPI_Isend(workbc_int[PACK][dir],
	    dirset[BOUNDPRIMTYPE][dir][DIRSIZE]/NPR,
	    MPI_INT,
	    dirset[BOUNDPRIMTYPE][dir][DIROTHER],
	    dirset[BOUNDPRIMTYPE][dir][DIRTAGS],
	    MPI_COMM_WORLD,
	    &requests[dir*2+REQSEND]);

}
#endif


void unpack_int(int dir,int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM],int prim[][N2M][N3M][NUMPFLAGS])
{
  // dir is direction receiving from
  int i,j,k;
  int pl;
  int bci;

  bci=0;
  PACKLOOP_INT(i,j,k
	   ,dirset[BOUNDPRIMTYPE][dir][DIRUSTART1]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRUSTOP1]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRUSTART2]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRUSTOP2]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRUSTART3]
	   ,dirset[BOUNDPRIMTYPE][dir][DIRUSTOP3]
	       ){
    prim[i][j][k][pl]=workbc_int[UNPACK][dir][bci++];
  }
}

