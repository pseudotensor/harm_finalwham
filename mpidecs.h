extern int romiocoliter;
extern int periodicx1, periodicx2, periodicx3;
extern int mpiperiodicx1, mpiperiodicx2, mpiperiodicx3;
extern int skipix1, reflectix1, reflectox1;
extern int skipix2, reflectix2, reflectox2;
extern int skipix3, reflectix3, reflectox3;
extern int intix1, intox1, intix2, intox2, intix3, intox3;
extern int skipintix1, skipintix2, skipintix3;
extern int ncpux1, ncpux2, ncpux3;
extern int myid, numprocs;
extern FILE* log_file;
extern FILE* fail_file;
extern FILE* logfull_file;
extern char myidtxt[MAXFILENAME];
extern int totalzones, realtotalzones;
extern int rtotalzones;
extern int itotalzones;
extern int sizes[COMPDIM + 1][MAXCPUS];
extern int isizes[COMPDIM + 1][MAXCPUS];
extern int totalsize[COMPDIM + 1];
extern int itotalsize[COMPDIM + 1];
extern int mycpupos[COMPDIM + 1];		// my position amongst the cpus
extern int dirset[NUMBOUNDTYPES][COMPDIM*2][DIRNUMVARS];
extern int srdir[6];			// which direction this cpu
				// sends/receives normal interior data
extern int startpos[COMPDIM + 1];
extern int endpos[COMPDIM + 1];		// startj and endj are where this CPU
				// located on full grid 
extern int *startpos0[COMPDIM+1];
extern int *endpos0[COMPDIM+1];
extern int *mycpupos0[COMPDIM+1];

extern int procnamelen;
#if(USEMPI)
extern char processor_name[MPI_MAX_PROCESSOR_NAME];
extern MPI_Status mpichstatus;

extern FTYPE workbca[2][COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM];	// [1=out/2=in][0=right,2=up,1=left,3=down,4=out,5=in][datawidth]
extern int workbc_inta[2][COMPDIM * 2][NBIGBND * NBIGSM];	// [1=out/2=in][0=right,2=up,1=left,3=down,4=out,5=in][datawidth]
#endif
extern FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM];
extern int (*workbc_int)[COMPDIM * 2][NBIGBND * NBIGSM];

// MPI transmit vars, so minimum local code changes
extern FTYPE ndtsend, bsq_maxsend;

// for data output
extern int nextbuf,numcolumns;
extern int bufferoffset;
extern int joniosize,writebufsize;
