/* 
   To start mpi run, first do:

   mpirun -np 4 ./grmhd

   e.g. 4 cpus using mpi:

   rm nohup.out ; nohup sh -c 'mpirun -np 4 ./grmhd' &

   Note: cannot have any cpu write to same file or pipe to same file
   without parallel I/O or blocking through each CPU

 */

#define USEMPI 0

// whether to use MPI over GM
// This forces no forks to be called.  Can be used for non-gm too
#if(USEMPI==1)
#include "mpi.h"
// can't use ROMIO unless file system is shared across all CPUs (i.e. not BH)
// ROMIO for tungsten since problems with minmem method.
#define USEROMIO 0 // choice, whether to use ROMIO parallel I/O package
#define USEGM 0	   // choice (avoids system/fork/etc calls)
#else
#define USEROMIO 0 // no choice
#define USEGM 0			// always 0, can't have GM without MPI
#endif



// whether to simultaneously compute and transfer bc (i.e. actually use non-blocking with a purpose).
// -1: use old loop even, super NONONONO!
// 0: no
// 1: yes with if type loops over fewer loops
// 2: yes with more loops without if, more blocks
#define SIMULBCCALC -1


#define MPIEQUALNONMPI 1
// 0= relax condition for MPI to equal non-MPI
// 1= guarentee not only that MPI boundaries are transfered but that an MPI run is identical to non-MPI run.




// this is the total stencil half width, where sts is full stencil size as: (sts-1)/2 which represents the one-sided safetey size in number of zones
// 2 even for parabolic, due to magnetic field stencil!
#define SAFESIZE (2)








///////////////////
//
// MNEMOICS
//
//////////////////////

#define INITROMIO 0
#define WRITECLOSEROMIO 1
#define READROMIO 2
#define READFREEROMIO 3


// still use MPI data type for communicating types to functions
#if(USEMPI==0)
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)2)
#define MPI_BYTE           ((MPI_Datatype)3)
#define MPI_SHORT          ((MPI_Datatype)4)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)5)
#define MPI_INT            ((MPI_Datatype)6)
#define MPI_UNSIGNED       ((MPI_Datatype)7)
#define MPI_LONG           ((MPI_Datatype)8)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)9)
#define MPI_FLOAT          ((MPI_Datatype)10)
#define MPI_DOUBLE         ((MPI_Datatype)11)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)12)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)13)
#endif




// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define MPI_FTYPE MPI_FLOAT
#elif(REALTYPE==DOUBLETYPE)
#define MPI_FTYPE MPI_DOUBLE
#elif(REALTYPE==LONGDOUBLETYPE)
#define MPI_FTYPE MPI_LONG_DOUBLE
#endif


#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define MPI_SFTYPE MPI_FLOAT
#elif(SENSITIVE==DOUBLETYPE)
#define MPI_SFTYPE MPI_DOUBLE
#elif(SENSITIVE==LONGDOUBLETYPE)
#define MPI_SFTYPE MPI_LONG_DOUBLE
#endif

#if(COUNTTYPE==DOUBLETYPE)
#define MPI_CTYPE MPI_DOUBLE
#elif(COUNTTYPE==LONGLONGINTTYPE)
#define MPI_CTYPE MPI_LONG_LONG_INT
#endif


#define BUFFERMAP (bufferoffset+(k*N1*N2+j*N1+i)*numcolumns+nextbuf++)
#define BUFFERMAP2 (k*N1*N2+j*N1+i)
#define BUFFERINIT0 bufferoffset=0
// mpi uses BUFFERINIT in global.h as well


#define PACK 1
#define UNPACK 2

#define REQRECV 0
#define REQSEND 1

#define DIRNUMVARS 19
#define DIRIF      0
#define DIRSIZE    1
#define DIROTHER   2
#define DIRTAGS    3
#define DIRTAGR    4
#define DIROPP     5
#define DIRPSTART1 6
#define DIRPSTOP1  7
#define DIRUSTART1 8
#define DIRUSTOP1  9
#define DIRPSTART2 10
#define DIRPSTOP2  11
#define DIRUSTART2 12
#define DIRUSTOP2  13
#define DIRPSTART3 14
#define DIRPSTOP3  15
#define DIRUSTART3 16
#define DIRUSTOP3  17
#define DIRNUMPR   18


#define TEXTOUTPUT 0
#define BINARYOUTPUT 1
#define MIXEDOUTPUT 2 // means header is text and dump is binary (handled by dump_gen()

#define UNSORTED 0
#define SORTED 1


// simple algorithm, but eats alot of memory on cpu=0 (unbounded) if doing sorted output
#define MPICOMBINESIMPLE 0
// more parallel:
#define MPICOMBINEMINMEM 1 // homebrew, but buggy on tungsten/mako, no problem on BH cluster -- ever.
#define MPICOMBINEROMIO 2 // requires romio package

// for various uses (dumping and special boundary/comp interchange routine
#define STAGEM1 (-1)
#define STAGE0 0
#define STAGE1 1
#define STAGE2 2
#define STAGE3 3
#define STAGE4 4
#define STAGE5 5
#define STAGE6 6
#define STAGE7 7

// #define DATADIR "./"
#define DATADIR ""

// extention for data files
#define DATEXT ".dat"
#define PAREXT ".par"
#define INEXT ".in"
#define OUTEXT ".out"
#define PPEXT ".pp"

#define CPUTXT ".%04d"

// also includes 3 flags: 2 b^2 and 1 utoprim fail flag
#define PLOOPMPI(pr,num) for(pr=0;pr<num;pr++)


#if(SIMULBCCALC==1)
// GODMARK : not setup for 3D and didn't yet work for 2D

// we introduce a stage-dependent loop (non-interfering stages and parts of stages)
// these are completely general and change shape depending upon absolute specified ranges.  All ranges are as specified as if normal non-MPI
#define MIDDLEI(istop) (istop-(N1-1-SAFESIZE)/2)
// safe left i across j
#define STAGECONDITION0(istart,istop,jstart,jstop) ((j>=jstart+SAFESIZE)&&(j<=jstop-SAFESIZE)&&(i>=istart+SAFESIZE)&&(i<=MIDDLEI(istop)))
// safe right i across j
#define STAGECONDITION1(istart,istop,jstart,jstop) ((j>=jstart+SAFESIZE)&&(j<=jstop-SAFESIZE)&&(i>=MIDDLEI(istop)+1)&&(i<=istop-SAFESIZE))
// left unsafe across j
#define STAGECONDITION20(istart,istop,jstart,jstop) ((j>=jstart)&&(j<=jstop)&&(i>=istart)&&(i<=istart+SAFESIZE-1))
// right unsafe across j
#define STAGECONDITION21(istart,istop,jstart,jstop) ((j>=jstart)&&(j<=jstop)&&(i>=istop+1-SAFESIZE)&&(i<=istop))
// upper j unsafe across i
#define STAGECONDITION22(istart,istop,jstart,jstop) ((j>=jstart)&&(j<=SAFESIZE-1)&&(i>=istart+SAFESIZE)&&(i<=istop-SAFESIZE))
// lower j unsafe across i
#define STAGECONDITION23(istart,istop,jstart,jstop) ((j>=N2-SAFESIZE)&&(j<=jstop)&&(i>=istart+SAFESIZE)&&(i<=istop-SAFESIZE))


#define STAGECONDITION2(istart,istop,jstart,jstop) (STAGECONDITION20(istart,istop,jstart,jstop)||STAGECONDITION21(istart,istop,jstart,jstop)||STAGECONDITION22(istart,istop,jstart,jstop)||STAGECONDITION23(istart,istop,jstart,jstop))

#define STAGECONDITION(istart,istop,jstart,jstop) ((stage==STAGEM1)||(stage==STAGE0)&&(STAGECONDITION0(istart,istop,jstart,jstop))||(stage==STAGE1)&&(STAGECONDITION1(istart,istop,jstart,jstop))||(stage==STAGE2)&&(STAGECONDITION2(istart,istop,jstart,jstop)))


#define CZLOOP ZLOOP if(STAGECONDITION(0,N1-1,0,N2-1))
#define CZSLOOP(istart,istop,jstart,jstop) ZSLOOP(istart,istop,jstart,jstop) if(STAGECONDITION(istart,istop,jstart,jstop))
/*
#define CZLOOP ZLOOP
#define CZSLOOP(istart,istop,jstart,jstop) ZSLOOP(istart,istop,jstart,jstop)
*/
// need another set of condition loops for flux calc since need fluxes a bit farther out than values.
// note that flux doesn't need more safe zones! (not even flux_ct)
// flux_ct requires F1 flux to be computed an extra i+1, and j-1,j,j+1
// flux_ct requires F2 flux to be computed an extra j+1, and i-1,i,i+1
// again, these outer fluxes don't use more primitive variables outside the safe zone

// fluxes stored ok
#define FZLOOP(istart,jstart) CZSLOOP(istart,N1,jstart,N2)

// emf requires its own storage since don't want to recalculate emf for each stage, only do necessary calculations
// this is already provided by the static array in flux_ct
// i+1 and j+1
#define EMFZLOOP CZSLOOP(0,N1,0,N2)
// used inside flux_ct()
// just stores fluxes, which we have F1 and F2 for storage
// i+1
#define F1CTZLOOP CZSLOOP(0,N1,0,N2-1,N3-1)
// j+1
#define F2CTZLOOP CZSLOOP(0,N1-1,0,N2,N3-1)
// k+1
#define F2CTZLOOP CZSLOOP(0,N1-1,0,N2-1,N3)

// also need a loop for dq-like calculations (i.e. for get_bsqflags()) i-1,i+1  and j-1,j+1  this is safe too since stencil is no bigger than loop size
// must keep new dq1 and dq2 in storage since each stage gets new set but calculations require previous stage values
#define DQZLOOP CZSLOOP(-1,N1,-1,N2)

// goes over only primitive quantities for rescale()
// affects pr, but undone, so no special storage needed
//#define PREDQZLOOP CZSLOOP(-N1BND,N1-1+N1BND,-N2BND,N2-1+N2BND)
#define PREDQZLOOP FULLLOOP


#elif(SIMULBCCALC==2)
#define STAGESETUPM1(istart,istop,jstart,jstop,is,ie,js,je) js=jstart;je=jstop;is=istart;ie=istop;
#define MIDDLEI(istop) (istop-(N1-1-SAFESIZE)/2)
// safe left i across j
#define STAGESETUP0(istart,istop,jstart,jstop,is,ie,js,je) js=jstart+SAFESIZE;je=jstop-SAFESIZE;is=istart+SAFESIZE;ie=MIDDLEI(istop);
// safe right i across j
#define STAGESETUP1(istart,istop,jstart,jstop,is,ie,js,je) js=jstart+SAFESIZE;je=jstop-SAFESIZE;is=MIDDLEI(istop)+1;ie=istop-SAFESIZE;
// left unsafe across j
#define STAGESETUP20(istart,istop,jstart,jstop,is,ie,js,je) js=jstart;je=jstop;is=istart;ie=istart+SAFESIZE-1;
// right unsafe across j
#define STAGESETUP21(istart,istop,jstart,jstop,is,ie,js,je) js=jstart;je=jstop;is=istop+1-SAFESIZE;ie=istop;
// upper j unsafe across i
#define STAGESETUP22(istart,istop,jstart,jstop,is,ie,js,je) js=jstart;je=SAFESIZE-1;is=istart+SAFESIZE;ie=istop-SAFESIZE;
// lower j unsafe across i
#define STAGESETUP23(istart,istop,jstart,jstop,is,ie,js,je) js=N2-SAFESIZE;je=jstop;is=istart+SAFESIZE;ie=istop-SAFESIZE;

#define STAGECONDITION(istart,istop,jstart,jstop,is,ie,js,je) if(stage==STAGEM1){ STAGESETUPM1(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE0){ STAGESETUP0(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE1){ STAGESETUP1(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE2){ STAGESETUP20(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE3){ STAGESETUP21(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE4){ STAGESETUP22(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE5){ STAGESETUP23(istart,istop,jstart,jstop,is,ie,js,je) }

#define STAGELOOP(is,ie,js,je) for(i=is;i<=ie;i++) for(j=js;j<=je;j++)

#define TYPE2 (0)
#define CZSLOOP(istart,istop,jstart,jstop,is,ie,js,je) STAGECONDITION(istart,istop,jstart,jstop,is,ie,js,je) STAGELOOP(is,ie,js,je)
#define CZLOOP CZSLOOP(0,N1-1,0,N2-1,isc,iec,jsc,jec)
#define FZLOOP(istart,jstart) CZSLOOP(istart,N1,jstart,N2,isc,iec,jsc,jec)
#define EMFZLOOP CZSLOOP(0,N1,0,N2,isc,iec,jsc,jec)
#define F1CTZLOOP CZSLOOP(0,N1,0,N2-1,N3-1,isc,iec,jsc,jec)
#define F2CTZLOOP CZSLOOP(0,N1-1,0,N2,N3-1,isc,iec,jsc,jec)
#define F3CTZLOOP CZSLOOP(0,N1-1,0,N2-1,N3,isc,iec,jsc,jec)
#define DQZLOOP CZSLOOP(-1,N1,-1,N2,isc,iec,jsc,jec)
#define PREDQZLOOP CZSLOOP(-N1BND,N1-1+N1BND,-N2BND,N2-1+N2BND,isc,iec,jsc,jec)


/*
#define TYPE2 (1)
#define CZLOOP STAGELOOP(isc,iec,jsc,jec)
// ief1 and jef1 are equal to ief2 and jef2, so ok
#define FZLOOP(istart,jstart) STAGELOOP(istart,ief1,jstart,jef1)
#define EMFZLOOP STAGELOOP(ise,iee,jse,jee)
#define F1CTZLOOP STAGELOOP(isf1ct,ief1ct,jsf1ct,jef1ct)
#define F2CTZLOOP STAGELOOP(isf2ct,ief2ct,jsf2ct,jef2ct)
#define DQZLOOP STAGELOOP(isdq,iedq,jsdq,jedq)
#define PREDQZLOOP STAGELOOP(ispdq,iepdq,jspdq,jepdq)
*/










#else



////////////////////////////////////////////
//
// normal non-MPI or standard MPI
//
////////////////////////////////////////////



#define CZLOOP ZLOOP
// can use CZSLOOP if don't care about extra calculations
#define CZSLOOP(istart,istop,jstart,jstop,kstart,kstop) ZSLOOP(istart,istop,jstart,jstop,kstart,kstop)
// these are special loops with very careful calculation ranges to avoid extra calculations
// but assumes are global variables which are assigned, so remains intact under whatever circumstances needed next
#define FZLOOP(istart,jstart,kstart) ZSLOOP(istart,OUTM1,jstart,OUTM2,kstart,OUTM3)
#define EMFZLOOP ZSLOOP(0,OUTM1,0,OUTM2,0,OUTM3)
#define PREEMFZLOOP ZSLOOP(INM1,OUTM1,INM2,OUTM2,INM3,OUTM3)
#define F1CTZLOOP ZSLOOP(0,OUTM1,0,N2-1,0,N3-1)
#define F2CTZLOOP ZSLOOP(0,N1-1,0,OUTM2,0,N3-1)
#define F3CTZLOOP ZSLOOP(0,N1-1,0,N2-1,0,OUTM3)
#define PREDQZLOOP FULLLOOP
#define DQZLOOP ZSLOOP(INM1,OUTM1,INM2,OUTM2,INM3,OUTM3)



#endif





////////////////////////////////////////////
//
// function declarations
//
////////////////////////////////////////////



extern void init_mpi(int argc, char *argv[]);
extern void init_genfiles(int gopp);
extern int init_MPI(int argc, char *argv[]);
extern void init_placeongrid(void);
extern int myexit(int call_code);

extern int bound_mpi(int boundstage, int boundtype, FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]);

extern int bound_mpi_int(int boundstage, int prim[][N2M][N3M][NUMPFLAGS]);
extern void pack(int dir, int boundtype, FTYPE prim[][N2M][N3M][NPR],FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM]);
extern void unpack(int dir, int boundtype, FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM],FTYPE prim[][N2M][N3M][NPR]);
extern void pack_int(int dir,int prim[][N2M][N3M][NUMPFLAGS],int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM]);
extern void unpack_int(int dir,int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM],int prim[][N2M][N3M][NUMPFLAGS]);
#if(USEMPI)
extern void sendrecv(int dir,int boundtype,FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM],MPI_Request *requests);
extern void sendrecv_int(int dir,int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM],MPI_Request *requests);
extern void recvwait(int dir,MPI_Request *request);
extern void sendwait(int dir,MPI_Request *request);
#endif

extern void mpimax(SFTYPE*maxptr);
extern void mpimin(SFTYPE*minptr);
extern void mpiisum(int*maxptr);
extern void mpiisum0(int*sumptr, int recvid);
extern void mpildsum0(long int*sumptr, int recvid);
extern void mpifmin(FTYPE*minptr);


extern int getsizeofdatatype(MPI_Datatype datatype);

extern void mpiio_init(int bintxt, int sorted, FILE ** fp, long headerbytesize, int which, char *filename, int numcolumns,
		       MPI_Datatype datatype, void **jonio, void **writebuf);
extern void mpiio_combine(int bintxt, int sorted,
			  int numcolumns, MPI_Datatype datatype,
			  FILE ** fp, void *jonio, void *writebuf);
extern void mpiio_seperate(int bintxt, int sorted, int stage,
			   int numcolumns,
			   MPI_Datatype datatype, FILE ** fp, void *jonio,
			   void *writebuf);
extern void mpiios_init(int bintxt, int sorted, FILE ** fp, int which, int headerbytesize, char *filename, int numcolumns,
			MPI_Datatype datatype, void **jonio, void **writebuf);
extern void mpiiomin_final(int numcolumns,FILE **fp, void *jonio, void *writebuf);
extern void mpiio_minmem(int readwrite, int whichdump, int i, int j, int k, int bintxt, int sorted,
			      int numcolumns, MPI_Datatype datatype,
			      FILE ** fpptr, void *jonio, void *writebuf);


extern void mpiioromio_init_combine(int operationtype, int which, long headerbytesize, char *filename, int numcolumns,MPI_Datatype datatype, void **writebufptr, void *writebuf);

#if(USEMPI)
extern void mpiios_combine(int bintxt, MPI_Datatype datatype, int numcolumns,
			   FILE ** fp, void *jonio, void *writebuf);
extern void mpiios_seperate(int bintxt, int stage, MPI_Datatype datatype, int numcolumns,
			    FILE ** fp, void *jonio,
			    void *writebuf);
extern void mpiiotu_combine(MPI_Datatype datatype, int numcolumns,
			    FILE ** fp, void *writebuf);
extern void mpiiotu_seperate(int stage, MPI_Datatype datatype, int numcolumns,
			     FILE ** fp,void *writebuf);
#endif

#define MYOUT stderr		// normally stderr
