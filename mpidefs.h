int romiocoliter;
int periodicx1, periodicx2, periodicx3;
int mpiperiodicx1, mpiperiodicx2, mpiperiodicx3;
int skipix1, reflectix1, reflectox1;
int skipix2, reflectix2, reflectox2;
int skipix3, reflectix3, reflectox3;
int intix1, intox1, intix2, intox2, intix3, intox3;
int skipintix1, skipintix2, skipintix3;
int ncpux1, ncpux2, ncpux3;
int myid, numprocs;
FILE *log_file;
FILE *fail_file;
FILE *logfull_file;
char myidtxt[MAXFILENAME];
int totalzones, realtotalzones;
int rtotalzones;
int itotalzones;
int sizes[COMPDIM + 1][MAXCPUS];
int isizes[COMPDIM + 1][MAXCPUS];
int totalsize[COMPDIM + 1];
int itotalsize[COMPDIM + 1];
int mycpupos[COMPDIM + 1];		// my position amongst the cpus
int dirset[NUMBOUNDTYPES][COMPDIM*2][DIRNUMVARS];
int srdir[6];			// which direction this cpu
				// sends/receives normal interior data
int startpos[COMPDIM + 1];
int endpos[COMPDIM + 1];		// startj and endj are where this CPU
				// located on full grid 
int *startpos0[COMPDIM+1];
int *endpos0[COMPDIM+1];
int *mycpupos0[COMPDIM+1];

int procnamelen;
#if(USEMPI)
char processor_name[MPI_MAX_PROCESSOR_NAME];
MPI_Status mpichstatus;

FTYPE workbca[2][COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM];	// [1=out/2=in][0=right,2=up,1=left,3=down,4=out,5=in][datawidth]
int workbc_inta[2][COMPDIM * 2][NBIGBND * NBIGSM];	// [1=out/2=in][0=right,2=up,1=left,3=down,4=out,5=in][datawidth]
#endif
FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM];
int (*workbc_int)[COMPDIM * 2][NBIGBND * NBIGSM];

// MPI transmit vars, so minimum local code changes
FTYPE ndtsend, bsq_maxsend;

// for data output
int nextbuf,numcolumns;
int bufferoffset;
int joniosize,writebufsize;
