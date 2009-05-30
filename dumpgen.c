#include "decs.h"



int dump_gen(int readwrite, long dump_cnt, int bintxt, int whichdump, MPI_Datatype datatype, char *fileprefix, char *fileformat, char *filesuffix, int (*headerfun) (int bintxt, FILE*headerptr),int (*setgetcontent) (int i, int j, int k, MPI_Datatype datatype, void*setbuf))
{
  int i = 0, j = 0, k = 0, l = 0, col = 0;
  FILE **fpp;
  char dfnam[MAXFILENAME];
  char dfnamreal[MAXFILENAME];
  char localfileformat[MAXFILENAME];
  void *jonio;
  void *writebuf;
  char truemyidtxt[MAXFILENAME];
  char filerw[MAXFILENAME];
  FILE *headerptr;
  int numfiles,coliter;
  void *setbuf;
  int sizeofdatatype;
  int romiocloopend;
  int headerbintxt;
  int dumpbintxt;
  fpos_t headerendpos;
  long headerbytesize;
  int binextension;


  ////////////
  //
  // setup file format for header and dump
  //
  ////////////

  if(readwrite==READFILE){
	strcpy(filerw,"rb");  //atch
	//if(bintxt==BINARYOUTPUT)
    //else strcpy(filerw,"rt");
  }
  else if(readwrite==WRITEFILE){
	strcpy(filerw,"wb");  //atch
    //if(bintxt==BINARYOUTPUT) strcpy(filerw,"w");
	//else strcpy(filerw,"wt");
  }

  if(bintxt==BINARYOUTPUT) headerbintxt=dumpbintxt=BINARYOUTPUT;
  else if(bintxt==TEXTOUTPUT) headerbintxt=dumpbintxt=TEXTOUTPUT;
  else if(bintxt==MIXEDOUTPUT){
    headerbintxt=TEXTOUTPUT;
    dumpbintxt=BINARYOUTPUT;
  }



  numcolumns=dnumcolumns[whichdump];
  docolsplit=DOCOLSPLIT[whichdump]; // docolsplit global var for now


  ////////////////////
  //
  // See if enough HD space
  //
  ////////////////////

  if(mpicombine){
    if(dumpbintxt==BINARYOUTPUT){
      if(myid==0) isenoughfreespace((unsigned long long)(sizeof(FTYPE))*(unsigned long long)numcolumns*(unsigned long long)(totalsize[1])*(unsigned long long)(totalsize[2])*(unsigned long long)(totalsize[3]) );
      else isenoughfreespace(0);
    }
    else{// text
      if(myid==0) isenoughfreespace((unsigned long long)(22)*(unsigned long long)numcolumns*(unsigned long long)(totalsize[1])*(unsigned long long)(totalsize[2])*(unsigned long long)(totalsize[3]) );
      else isenoughfreespace(0);
    }
  }
  else{
    if(dumpbintxt==BINARYOUTPUT){
      isenoughfreespace((unsigned long long)(sizeof(FTYPE))*(unsigned long long)numcolumns*(unsigned long long)(N1)*(unsigned long long)(N2)*(unsigned long long)(N3) );
    }
    else{// text
      isenoughfreespace((unsigned long long)(21)*(unsigned long long)numcolumns*(unsigned long long)(N1)*(unsigned long long)(N2)*(unsigned long long)(N3) );
    }
  }


  /////////////////////
  //
  // Allocate memory for setting up setbuf
  //
  //////////////////////

  sizeofdatatype=getsizeofdatatype(datatype);
  setbuf=malloc(numcolumns*sizeofdatatype);
  if(setbuf==NULL){
    dualfprintf(fail_file,"cannot allocate memory to setbuf in %s %s with numcolumns=%d and sizeofdatatype=%d\n",fileprefix,filesuffix,numcolumns,sizeofdatatype);
    myexit(1);
  }


  //  trifprintf("numcolumns=%d sizeofdatatype=%d setbuf=%d total=%d\n",numcolumns,sizeofdatatype,setbuf,numcolumns*sizeofdatatype);

  //////////////////////////////////
  //
  //  Set up DOCOLSPLIT for normal and ROMIO loop
  //
  ///////////////////////////////////



  if(docolsplit){
    numfiles=numcolumns;
    if(mpicombine&&USEMPI&&USEROMIO) romiocloopend=numfiles;
    else romiocloopend=1;
  }
  else{
    numfiles=1;
    romiocloopend=1;
  }


  //////////////////////////////////
  //
  //  Define file output and open it
  //
  ///////////////////////////////////

  // say whether .bin is allowed or not if binary
  if(fileprefix[0]=='i') binextension=0; // images don't need binary extension
  else binextension=1;


  // sometimes all CPUs need to know filename (e.g. ROMIO)
  // setup file suffix
  if((dumpbintxt==BINARYOUTPUT)&&(binextension)){
    if(USEMPI&&(mpicombine==0)&&(numprocs>1)) sprintf(truemyidtxt,".bin.%04d",myid);
    else strcpy(truemyidtxt,".bin");
  }
  else{
    if(USEMPI&&(mpicombine==0)&&(numprocs>1)) sprintf(truemyidtxt,".%04d",myid);
    else strcpy(truemyidtxt,"");
  }

  // setup filename
  if(dump_cnt>=0){
    strcpy(localfileformat,"%s");
    strcat(localfileformat,fileformat);
    strcat(localfileformat,"%s");
    strcat(localfileformat,"%s");
    sprintf(dfnam, localfileformat, fileprefix, dump_cnt, filesuffix, truemyidtxt);
  }
  else{ // then no file number wanted (i.e. for gdump())
    sprintf(dfnam, "%s%s%s", fileprefix, filesuffix, truemyidtxt);
  }



  ////////////////
  //
  // open files, or open files for header if mpicombine==1, which for mpicombine==1 gets reopened later by MPI routines
  //
  ///////////////

  if((USEMPI&&(myid==0)&&(mpicombine==1))||(mpicombine==0)){// for mpicombine==1 even with ROMIO, real filename and header not needed

    // only one CPU does header if mpicombine==1, header+dump done in all CPUs if mpicombine==0
    // create files for each column, or each column's header if mpicombine==1
    if((fpp=(FILE**)malloc(sizeof(FILE*)*numfiles))==NULL){
      dualfprintf(fail_file,"couldn't open fpp in dump()\n");
      myexit(2);
    }// now fpp[i] indexes a list of file pointers


    // setup each file corresponding to each column
    COLLOOP(coliter){
      if(docolsplit&&(numfiles>1)){
	sprintf(dfnamreal,"%s-col%04d",dfnam,coliter);
      }
      else strcpy(dfnamreal,dfnam);
      
      if ((fpp[coliter] = fopen(dfnamreal, filerw)) == NULL) {
	dualfprintf(fail_file, "error opening %s %s file\n",fileprefix,filesuffix);
	myexit(2);
      }
      //////////////////////////////////
      //
      //  read or write header: header is read/written in whatever style chosen to the top of each dump file created
      //
      ///////////////////////////////////
      headerfun(headerbintxt,fpp[coliter]); // outputs header to each column file (or just one file, or all CPU files, etc.)
      headerbytesize=ftell(fpp[coliter]);
    }
    
    // don't close if mpicombine==0, since used in a moment, else mpicombine==1 it's reopened by MPI routines
    if(USEMPI&&(myid==0)&&(mpicombine==1)) COLLOOP(coliter) fclose(fpp[coliter]); // will get reopened later by MPI routines


  }
  // need to broadcast the header size to other CPUs for ROMIO
#if(USEMPI&&USEROMIO)
  MPI_Bcast(&headerbytesize,1,MPI_LONG,0,MPI_COMM_WORLD);
#endif

    
  ///////////////////////////////////////////////////////////
  //
  // loop over columns for per-column buffer ROMIO dump
  //
  //
  ////////////////////////////////////////////////////////////


  ROMIOCOLLOOP(romiocoliter) { // only loop if mpicombine&&USEMPI&&USEROMIO&&docolsplit==1
    if(romiocloopend>1) trifprintf("romiocoliter=%d of romiocloopend=%d\n",romiocoliter,romiocloopend);


    
    // setup MPI buffer if mpicombine==1
    if( mpicombine == 0 ) { // then one file per CPU if USEMPI or just normal file writing on 1CPU
      writebuf=NULL;
    }
    else mpiio_init(dumpbintxt,sortedoutput, fpp, headerbytesize, readwrite, dfnam, numcolumns, datatype, &jonio, &writebuf);
    // if USEROMIO==1 then numcolumns interpreted properly for docolsplit


    if(readwrite==READFILE){
      //////////////////////////////////
      //
      // read DUMP 
      //
      //////////////////////////////////
      
      if (mpicombine == 1) {
#if(USEMPI)
	mpiio_seperate(binaryoutput,sortedoutput, STAGE1, numcolumns, datatype, fpp, jonio, writebuf);
#endif
      }
    }



    //////////////////
    //
    // DUMP LOOP
    //
    //////////////////



    if(readwrite==READFILE){
      BUFFERINIT0;
      DUMPGENLOOP {
	// buffer init starts the parallel index
	BUFFERINIT;
	// initialize to 0th column
	COLINIT;

	///////////////////////////////////////
	//
	// READFILE
	//
	//////////////////////
	if((mpicombine)&&(truempicombinetype==MPICOMBINEMINMEM)) mpiio_minmem(READFILE,whichdump,i,j,k,dumpbintxt,sortedoutput,numcolumns,datatype, fpp,jonio,writebuf);

	
	// read all at once
	myfread(dumpbintxt,datatype,setbuf,0,numcolumns,i,j,k,fpp,writebuf);
	
	// check
	if(nextbuf!=numcolumns){
	  dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns/buffers attempted (nextbuf=%d)\n",numcolumns,nextbuf);
	  myexit(1);
	}
	
	// get the content of 1 row
	setgetcontent(i,j,k,datatype,setbuf);
	
	// check
	if(nextcol!=numcolumns){
	  dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns attempted (nextcol=%d)\n",numcolumns,nextcol);
	  myexit(1);
	}
      }// end DUMPGENLOOP
    }// end readwrite==READFILE
    else if(readwrite==WRITEFILE){
      BUFFERINIT0;
      DUMPGENLOOP {



	// buffer init starts the parallel index
	BUFFERINIT;
	// initialize to 0th column
	COLINIT;
	///////////////////////////////////////
	//
	// WRITEFILE
	//
	//////////////////////

	// set the content of 1 row
	setgetcontent(i,j,k,datatype,setbuf);

	// check
	if(nextcol!=numcolumns){
	  dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns attempted (nextcol=%d)\n",numcolumns,nextcol);
	  myexit(1);
	}

	// write all at once
	myfwrite(dumpbintxt,datatype,setbuf,0,numcolumns,i,j,k,fpp,writebuf);

	
	// check
	if(nextbuf!=numcolumns){
	  dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns/buffers attempted (nextbuf=%d)\n",numcolumns,nextbuf);
	  myexit(1);
	}


	// finish up this row
	if((mpicombine==0)&&(dumpbintxt==TEXTOUTPUT)) COLLOOP(coliter) fprintf(fpp[coliter],"\n");
	if((mpicombine)&&(truempicombinetype==MPICOMBINEMINMEM)) mpiio_minmem(WRITEFILE,whichdump,i,j,k,dumpbintxt,sortedoutput,numcolumns,datatype, fpp, jonio,writebuf);



      }// end DUMPGENLOOP
    }//end readwrite==WRITEFILE
  



    //////////////////
    //
    // Close dump file and write/close file if mpicombine==1
    //
    //////////////////


    if (mpicombine == 0){
      COLLOOP(coliter) if (fpp[coliter] != NULL) fclose(fpp[coliter]);
    }
    else{
#if(USEMPI)
      if(readwrite==WRITEFILE) mpiio_combine(dumpbintxt, sortedoutput, numcolumns, datatype, fpp, jonio, writebuf);
      else if(readwrite==READFILE) mpiio_seperate(binaryoutput,sortedoutput, STAGE2, numcolumns, datatype, fpp, jonio, writebuf);
#endif
    }

  }// end column loop for ROMIO&&docolsplit



  // free the set/get buffer
  if(setbuf!=NULL) free(setbuf);



  return (0);
}


