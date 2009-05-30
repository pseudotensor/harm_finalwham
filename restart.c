
/* restart functions; restart_init and restart_dump */

#include "decs.h"




int restart_init(int which)
{
  char ans[100];
  struct of_geom geom;
  int failreturn;
  FTYPE ucon[NDIM];
  FTYPE utmax=0.0;
  int i,j,k,pl;
  int failflag=0;
  extern int set_dt(FTYPE (*prim)[N2M][N3M][NPR], SFTYPE *dt);
  extern int init_grid_post_set_grid(void);
  extern int post_par_set(void);



  trifprintf("begin restart init\n");

  ranc(7);
  restart_read(which);
  // set all CPUs to have same restart header stuff (needed in mpicombine==1 mode)
  restart_read_defs();
  read_coord_parms(); // assumes doesn't overwrite what just set (i.e. what's written to restart file)

  ////////////////
  //
  // write header to screen to all CPUs
  //
  ////////////////
  fprintf(log_file,"header contents below\n"); fflush(log_file);
  write_restart_header(TEXTOUTPUT,log_file);

  ////////////////
  //
  // check if failed inside rdump file
  //
  ////////////////
  if(failed!=0){
    dualfprintf(log_file,"WARNING: failed=%d in rdump file.  Must clense the way between us.\n",failed);
    dualfprintf(log_file,"Setting failed=0\n");
    failed=0;
  }

  // certain global.h parameters may override restart parameters.  This happens, say, when the user decides not to DOAVG after restart, but DOAVGDIAG is still set.  Without this the code will segfault in diag.c
  // corresponds to "if" statements in initbase.c as well, where user has either choice or no choice.
  if(DOAVG==0) DOAVGDIAG=0;



  ////////////////////////
  //
  // test read by looking at written file
  //
  ////////////////////////
  restart_write(3);


  trifprintf("proc: %d t=%21.15g failed=%d\n", myid, t, failed);

  /* set metric functions */
  set_grid();

  trifprintf("proc: %d grid restart completed: failed=%d\n", myid,failed);

  // user post_set_grid function
  init_grid_post_set_grid();
  
  trifprintf("proc: %d init_grid_post_set_grid completed: failed=%d\n", myid,failed);

  // set parameters dependent on user parameters
  post_par_set();

  trifprintf("proc: %d post_par_set completed: failed=%d\n", myid,failed);



  ////////////////////////////////////////////////////////////////
  //
  // NOW have all parameters and data.
  //
  // NOW CHECK THE DATA and apply "fixups" if necessary
  //
  ////////////////////////////////////////////////////////////////




#if(CHECKRHONEGZERORESTART)

  // report if data has rho<0 or u<0
  ZLOOP{
    if(STEPOVERNEGRHO!=-1){
      if(p[i][j][k][RHO]<=0.0){
	dualfprintf(fail_file,"restart data has negative mass density at i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
	failflag++;
      }
    }
    if(STEPOVERNEGU!=-1){
      if(p[i][j][k][UU]<=0.0){
	dualfprintf(fail_file,"restart data has negative ie density at i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
	failflag++;
      }
    }

    // check 3-velocity
    if(WHICHVEL==VEL3){
      if(jonchecks){
	get_geometry(i,j,k,CENT,&geom);
	failreturn=check_pr(p[i][j][k],p[i][j][k],&geom,-2,-1);
	if(failreturn){
	  dualfprintf(fail_file,"restart data has large or imaginary u^t=%21.15g at i=%d j=%d k=%d.  I will attempt to correct.\n",1.0/sqrt(uttdiscr),startpos[1]+i,startpos[2]+j,startpos[3]+k);
	}
	if(1.0/sqrt(uttdiscr)>utmax) utmax=1.0/sqrt(uttdiscr);
	// need to settle over limit u^t's
	failreturn=check_pr(p[i][j][k],p[i][j][k],&geom,-1,-1);
	if(failreturn){
	  dualfprintf(fail_file,"restart data has imaginary u^t at i=%d j=%d k=%d.  Unable to correct.\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
	  return(1);
	}
      }
    }
  }
  if(WHICHVEL==VEL3){
    if(jonchecks){
      dualfprintf(fail_file,"max u^t of restart data=%21.15g\n",utmax);
    }
  }
#endif

  if(failflag>0){
    dualfprintf(fail_file,"Restart data has at least %d failures, please correct\n",failflag);
    return(1);
  }

#if(FIXUPAFTERRESTART)
  if(fixup(STAGEM1,p,-1)>=1)
    FAILSTATEMENT("restart.c:restart_init()", "fixup()", 1);

  trifprintf( "proc: %d fixup restart completed: failed=%d\n", myid,failed);
#endif

  /* bound */
  if (bound_prim(STAGEM1,p) >= 1) {
    fprintf(fail_file, "restart_init:bound_prim: failure\n");
    fflush(fail_file);
    return (1);
  }

  trifprintf( "proc: %d bound restart completed: failed=%d\n", myid,failed);
  
  if(pre_fixup(STAGEM1,p)>=1)
    FAILSTATEMENT("init.c:init()", "postbc_fixup()", 1);

#if(CHECKRHONEGZERORESTART)
  // make sure all zones are good now
  if(STEPOVERNEGRHO!=-1){
    FULLLOOP{
      if(p[i][j][k][RHO]<=0.0){
	dualfprintf(fail_file,"restart data has negative mass density at i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
	return(1);
	
      }
    }
    if(STEPOVERNEGU!=-1){
      if(p[i][j][k][UU]<=0.0){
	dualfprintf(fail_file,"restart data has negative ie density at i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
	return(1);
	
      }
    }
  }
#endif



  ////////////////////////
  //
  // Now adjust timestep since may be dt=0 or whatever if last run ended on tf
  //
  ///////////////////////
  set_dt(p,&dt);




  trifprintf( "proc: %d restart completed: failed=%d\n", myid,failed);
  
  trifprintf("end restart init\n");

  /* done! */
  return (0);

}





int restart_write(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
  int truedump_cnt;

  trifprintf("begin dumping rdump# %ld ... ",dump_cnt);

  whichdump=RDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/rdump");
  if(dump_cnt>=0) {
    strcpy(fileformat,"-%01ld"); // very quick restart
    truedump_cnt=dump_cnt;
  }
  else {
    strcpy(fileformat,"--%04ld"); // assumed at some longer cycle and never overwritten .. must rename this to normal format to use as real rdump.
    truedump_cnt=-dump_cnt-1;
  }
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,truedump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,write_restart_header,rdump_content)>=1) return(1);

  trifprintf("end dumping rdump# %ld ... ",dump_cnt);


  return(0);

}


int rdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  myset(datatype,p[i][j][k],0,NPRDUMP,writebuf);
  myset(datatype,unew[i][j][k],0,NPRDUMP,writebuf);

  return(0);
}








int restart_read(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin reading rdump# %ld ... ",dump_cnt);

  whichdump=RDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/rdump");
  strcpy(fileformat,"-%01ld");
  strcpy(filesuffix,"");
  
  if(dump_gen(READFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,read_restart_header,rdump_read_content)>=1) return(1);

  trifprintf("end reading rdump# %ld ... ",dump_cnt);


  return(0);

}


int rdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  myget(datatype,p[i][j][k],0,NPRDUMP,writebuf);
  myget(datatype,unew[i][j][k],0,NPRDUMP,writebuf);


  return(0);
}





// headerptr created and only used here OR passed a given pointer
int read_restart_header(int bintxt, FILE*headerptr)
{
  int ii;
  int k,dir;
  int enerregion, floor, tscale;
  int idum1,idum2,idum3;

  trifprintf("begin reading header of restart file\n");

  if(bintxt==BINARYOUTPUT){

    /* read in global variables, in binary */
    fread(&idum1, sizeof(int), 1, headerptr);
    fread(&idum2, sizeof(int), 1, headerptr);
    fread(&idum3, sizeof(int), 1, headerptr); // extra 3D thing

    // all cpus read the rest of header the same
    fread(&t, sizeof(SFTYPE), 1, headerptr);
    fread(&tf, sizeof(SFTYPE), 1, headerptr);
    fread(&nstep, sizeof(long), 1, headerptr);
    fread(&a, sizeof(FTYPE), 1, headerptr);
    fread(&gam, sizeof(FTYPE), 1, headerptr);
    fread(&cour, sizeof(FTYPE), 1, headerptr);
      
    fread(&DTd, sizeof(SFTYPE), 1, headerptr);
    fread(&DTavg, sizeof(SFTYPE), 1, headerptr);
    fread(&DTener, sizeof(SFTYPE), 1, headerptr);
    fread(&DTi, sizeof(SFTYPE), 1, headerptr);
    // fread(&DTr, sizeof(SFTYPE), 1, headerptr) ;
    fread(&DTr, sizeof(long), 1, headerptr);
    fread(&DTdebug, sizeof(SFTYPE), 1, headerptr);
    fread(&dump_cnt, sizeof(long), 1, headerptr);
    fread(&image_cnt, sizeof(long), 1, headerptr);
    fread(&rdump_cnt, sizeof(long), 1, headerptr);
    fread(&avg_cnt, sizeof(long), 1, headerptr);
    fread(&debug_cnt, sizeof(long), 1, headerptr);
      
      
    fread(&dt, sizeof(SFTYPE), 1, headerptr);
    fread(&lim, sizeof(int), 1, headerptr);
    fread(&TIMEORDER, sizeof(int), 1, headerptr);
    fread(&fluxmethod, sizeof(int), 1, headerptr);
    fread(&FLUXB, sizeof(int), 1, headerptr);
    fread(&UTOPRIMVERSION, sizeof(int), 1, headerptr);
    fread(&failed, sizeof(int), 1, headerptr);
      
    fread(&R0, sizeof(FTYPE), 1, headerptr);
    fread(&Rin, sizeof(FTYPE), 1, headerptr);
    fread(&Rout, sizeof(FTYPE), 1, headerptr);
    fread(&hslope, sizeof(FTYPE), 1, headerptr);
    fread(&defcoord, sizeof(FTYPE), 1, headerptr);

    fread(&BCtype[X1UP],sizeof(int), 1, headerptr);
    fread(&BCtype[X1DN],sizeof(int), 1, headerptr);
    fread(&BCtype[X2UP],sizeof(int), 1, headerptr);
    fread(&BCtype[X2DN],sizeof(int), 1, headerptr);
    fread(&BCtype[X3UP],sizeof(int), 1, headerptr); // 3D thing
    fread(&BCtype[X3DN],sizeof(int), 1, headerptr);

    // new May 6, 2003
    fread(&realnstep, sizeof(long), 1, headerptr);
    fread(&debugfail,sizeof(int), 1, headerptr);
    fread(&whichrestart,sizeof(int), 1, headerptr);
    fread(&cooling,sizeof(int), 1, headerptr);
    fread(&restartsteps[0],sizeof(long), 1, headerptr);
    fread(&restartsteps[1],sizeof(long), 1, headerptr);
    fread(&GAMMIEDUMP,sizeof(int), 1, headerptr);
    fread(&GAMMIEIMAGE,sizeof(int), 1, headerptr);
    fread(&GAMMIEENER,sizeof(int), 1, headerptr);
    fread(&DODIAGS,sizeof(int), 1, headerptr);
    fread(&DOENERDIAG,sizeof(int), 1, headerptr);
    fread(&DOGDUMPDIAG,sizeof(int), 1, headerptr);
    fread(&DORDUMPDIAG,sizeof(int), 1, headerptr);
    fread(&DODUMPDIAG,sizeof(int), 1, headerptr);
    fread(&DOAVGDIAG,sizeof(int), 1, headerptr);
    fread(&DOIMAGEDIAG,sizeof(int), 1, headerptr);
    fread(&DOAREAMAPDIAG,sizeof(int), 1, headerptr);
    fread(&POSDEFMETRIC,sizeof(int), 1, headerptr);
    fread(&periodicx1,sizeof(int), 1, headerptr);
    fread(&periodicx2,sizeof(int), 1, headerptr);
    fread(&periodicx3,sizeof(int), 1, headerptr); // 3D thing
    fread(&binaryoutput,sizeof(int), 1, headerptr);
    fread(&sortedoutput,sizeof(int), 1, headerptr);
    fread(&defcon,sizeof(FTYPE), 1, headerptr);
    fread(&SAFE,sizeof(FTYPE), 1, headerptr);
    fread(&RHOMIN,sizeof(FTYPE), 1, headerptr);
    fread(&UUMIN,sizeof(FTYPE), 1, headerptr);
    fread(&RHOMINLIMIT,sizeof(FTYPE), 1, headerptr);
    fread(&UUMINLIMIT,sizeof(FTYPE), 1, headerptr);
    fread(&BSQORHOLIMIT,sizeof(FTYPE), 1, headerptr);
    fread(&BSQOULIMIT,sizeof(FTYPE), 1, headerptr);
    fread(&GAMMAMAX,sizeof(FTYPE), 1, headerptr);
    fread(&GAMMADAMP,sizeof(FTYPE), 1, headerptr);
    fread(&GAMMAFAIL,sizeof(FTYPE), 1, headerptr);      
    // end new May 6, 2003

    PDUMPLOOP(k) fread(&prMAX[k],sizeof(FTYPE),1,headerptr);

    fread(&rescaletype,sizeof(int),1,headerptr);

    // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
    fread(&avgscheme,sizeof(int),1,headerptr);
    fread(&do_transverse_flux_integration,sizeof(int),1,headerptr);
    fread(&do_source_integration,sizeof(int),1,headerptr);
    fread(&do_conserved_integration,sizeof(int),1,headerptr);
    fread(&INVERTFROMAVERAGEIFFAILED,sizeof(int),1,headerptr);
    fread(&LIMIT_AC_PRIM_FRAC_CHANGE,sizeof(int),1,headerptr);
    fread(&MAX_AC_PRIM_FRAC_CHANGE,sizeof(FTYPE),1,headerptr);
    fread(&DOENOFLUX,sizeof(int),1,headerptr);


    // new June 6, 2003 (cumulatives)
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(k) fread(&pcumreg_tot[enerregion][dir][k],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fread(&fladdreg_tot[enerregion][k],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fread(&sourceaddreg_tot[enerregion][k],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fread(&Ureg_init_tot[enerregion][k],sizeof(FTYPE),1,headerptr);
    TSCALELOOP(tscale) FLOORLOOP(floor) fread(&failfloorcountlocal_tot[tscale][floor],sizeof(CTYPE),1,headerptr);
    // end new June 6,2003

    if(DODISS) ENERREGIONLOOP(enerregion) fread(&dissreg_tot[enerregion],sizeof(FTYPE),1,headerptr);
    
    if(DOLUMVSR) fread(&lumvsr_tot,sizeof(SFTYPE),ncpux1*N1,headerptr);
    if(DODISSVSR) fread(&dissvsr_tot,sizeof(SFTYPE),ncpux1*N1,headerptr);


 
  }
  // assumes headerptr is normal file pointer already defined
  else{
    fscanf(headerptr,RESTARTHEADER,
	   &idum1,&idum2,&idum3,
	   &t,&tf,&nstep,&a,&gam,&cour,
	   &DTd,&DTavg,&DTener,&DTi,&DTr,&DTdebug,&dump_cnt,&image_cnt,&rdump_cnt,&avg_cnt,&debug_cnt,
	   &dt,&lim,&TIMEORDER,&fluxmethod,&FLUXB,&UTOPRIMVERSION,&failed,
	   &R0,&Rin,&Rout,&hslope,&defcoord,
	   &BCtype[X1UP],&BCtype[X1DN],&BCtype[X2UP],&BCtype[X2DN],&BCtype[X3UP],&BCtype[X3DN],
	   &realnstep,&debugfail,&whichrestart,&cooling,&restartsteps[0],&restartsteps[1],&GAMMIEDUMP,&GAMMIEIMAGE,&GAMMIEENER,&DODIAGS,&DOENERDIAG,&DOGDUMPDIAG,&DORDUMPDIAG,&DODUMPDIAG,&DOAVGDIAG,&DOIMAGEDIAG,&DOAREAMAPDIAG,&POSDEFMETRIC,&periodicx1,&periodicx2,&periodicx3,&binaryoutput,&sortedoutput,&defcon,&SAFE,&RHOMIN,&UUMIN,&RHOMINLIMIT,&UUMINLIMIT,&BSQORHOLIMIT,&BSQOULIMIT,&GAMMAMAX,&GAMMADAMP,&GAMMAFAIL
	   );
    PDUMPLOOP(k) fscanf(headerptr,HEADERONEIN,&prMAX[k]);

    fscanf(headerptr,"%d",&rescaletype);

    // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
    fscanf(headerptr,"%d",&avgscheme);
    fscanf(headerptr,"%d",&do_transverse_flux_integration);
    fscanf(headerptr,"%d",&do_source_integration);
    fscanf(headerptr,"%d",&do_conserved_integration);
    fscanf(headerptr,"%d",&INVERTFROMAVERAGEIFFAILED);
    fscanf(headerptr,"%d",&LIMIT_AC_PRIM_FRAC_CHANGE);
    fscanf(headerptr,HEADERONEIN,&MAX_AC_PRIM_FRAC_CHANGE);
    fscanf(headerptr,"%d",&DOENOFLUX);


    // new June 6, 2003 (cumulatives)
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(k) fscanf(headerptr,HEADERONEIN,&pcumreg_tot[enerregion][dir][k]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fscanf(headerptr,HEADERONEIN,&fladdreg_tot[enerregion][k]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fscanf(headerptr,HEADERONEIN,&sourceaddreg_tot[enerregion][k]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fscanf(headerptr,HEADERONEIN,&Ureg_init_tot[enerregion][k]);
#if(COUNTTYPE==LONGLONGINTTYPE)
    TSCALELOOP(tscale) FLOORLOOP(floor) fscanf(headerptr,"%lld",&failfloorcountlocal_tot[tscale][floor]);
#elif(COUNTTYPE==DOUBLETYPE)
    TSCALELOOP(tscale) FLOORLOOP(floor) fscanf(headerptr,"%lf",&failfloorcountlocal_tot[tscale][floor]);
#endif
    // end new June 6,2003
    
    if(DODISS) ENERREGIONLOOP(enerregion) fscanf(headerptr,HEADERONEIN,&dissreg_tot[enerregion][0]);
    
    // assumes same CPU geometry during restart
    if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) fscanf(headerptr,HEADERONEIN,&lumvsr_tot[ii]);

    // assumes same CPU geometry during restart
    if(DODISSVSR) for(ii=0;ii<ncpux1*N1;ii++) fscanf(headerptr,HEADERONEIN,&dissvsr_tot[ii]);



    // flush to just after the header line in case binary read of data
    while(fgetc(headerptr)!='\n');
  }

  /////////////////
  //
  // some checks
  if (idum1 != totalsize[1]) {
    dualfprintf(fail_file, "error reading restart file; N1 differs\n");
    dualfprintf(fail_file, "got totalsize[1]=%d needed totalsize[1]=%d\n",idum1,totalsize[1]);
    myexit(3);
  }
  if (idum2 != totalsize[2]) {
    dualfprintf(fail_file, "error reading restart file; N2 differs\n");
    dualfprintf(fail_file, "got totalsize[2]=%d needed totalsize[2]=%d\n",idum2,totalsize[2]);
    myexit(4);
  }
  if (idum3 != totalsize[3]) {
    dualfprintf(fail_file, "error reading restart file; N3 differs\n");
    dualfprintf(fail_file, "got totalsize[3]=%d needed totalsize[3]=%d\n",idum3,totalsize[3]);
    myexit(4);
  }


  trifprintf("end reading header of restart file\n");

  return(0);
}


// all cpus should do this
int restart_read_defs(void)
{
  int enerregion;
  int floor,tscale;
  int dir,k;
  int ii;

  if(myid==0){
    ////////////
    //
    // Define for cpu=0 only, which will continue to keep track of the total after restart
    //
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(k) pcumreg[enerregion][dir][k]=pcumreg_tot[enerregion][dir][k];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fladdreg[enerregion][k]=fladdreg_tot[enerregion][k];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) sourceaddreg[enerregion][k]=sourceaddreg_tot[enerregion][k];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) Ureg_init[enerregion][k]=Ureg_init_tot[enerregion][k];

    TSCALELOOP(tscale) FLOORLOOP(floor) failfloorcountlocal[tscale][floor]=failfloorcountlocal_tot[tscale][floor];


    if(DODISS) ENERREGIONLOOP(enerregion) dissreg[enerregion][0]=dissreg_tot[enerregion][0];
    // assume dissfunpos[][][] not restored since zeroed out each dump

    if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) lumvsr[ii]=lumvsr_tot[ii];
    if(DODISSVSR) for(ii=0;ii<ncpux1*N1;ii++) dissvsr[ii]=dissvsr_tot[ii];

  }



  // now that CPU=0 has the restart header, pass to other CPUs
#if(USEMPI)
  MPI_Bcast(&t, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tf, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nstep, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gam, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cour, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&DTd, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DTavg, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DTener, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DTi, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DTr, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DTdebug, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dump_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&image_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rdump_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&avg_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&debug_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&dt, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&lim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TIMEORDER, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&fluxmethod, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&FLUXB, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&UTOPRIMVERSION, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&failed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&R0, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Rin, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Rout, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hslope, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&defcoord, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);

  MPI_Bcast(&BCtype[X1UP], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BCtype[X1DN], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BCtype[X2UP], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BCtype[X2DN], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BCtype[X3UP], 1, MPI_INT, 0, MPI_COMM_WORLD); // 3D thing
  MPI_Bcast(&BCtype[X3DN], 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&realnstep, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&debugfail, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&whichrestart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cooling, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&restartsteps[0], 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&restartsteps[1], 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GAMMIEDUMP, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GAMMIEIMAGE, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GAMMIEENER, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DODIAGS, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DOENERDIAG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DOGDUMPDIAG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DORDUMPDIAG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DODUMPDIAG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DOAVGDIAG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DOIMAGEDIAG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DOAREAMAPDIAG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&POSDEFMETRIC, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&periodicx1, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&periodicx2, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&periodicx3, 1, MPI_INT, 0, MPI_COMM_WORLD); // 3D thing
  MPI_Bcast(&binaryoutput, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sortedoutput, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&defcon, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&SAFE, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&RHOMIN, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&UUMIN, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&RHOMINLIMIT, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&UUMINLIMIT, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BSQORHOLIMIT, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BSQOULIMIT, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GAMMAMAX, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GAMMADAMP, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GAMMAFAIL, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);

  // new June 6, 2003 (cumulatives)
  // don't broadcast cumulatives or initial conditions, only needed by cpu=0  
  // end new June 6,2003

  PDUMPLOOP(k) MPI_Bcast(&prMAX[k], 1, MPI_FTYPE, 0, MPI_COMM_WORLD);

  MPI_Bcast(&rescaletype, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
  MPI_Bcast(&avgscheme,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&do_transverse_flux_integration,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&do_source_integration,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&do_conserved_integration,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&INVERTFROMAVERAGEIFFAILED,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&LIMIT_AC_PRIM_FRAC_CHANGE,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&MAX_AC_PRIM_FRAC_CHANGE,1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DOENOFLUX,1, MPI_INT, 0, MPI_COMM_WORLD);


#endif



  return(0);
}



int write_restart_header(int bintxt,FILE*headerptr)
{
  int k,dir,floor,tscale;
  int enerregion;
  int ii;

  trifprintf("begin writing header of restart file\n");

  if(bintxt==BINARYOUTPUT){

    /* write out key global variables, in binary */
    fwrite(&totalsize[1], sizeof(int), 1, headerptr);
    fwrite(&totalsize[2], sizeof(int), 1, headerptr);
    fwrite(&totalsize[3], sizeof(int), 1, headerptr); // extra 3D thing
      
    fwrite(&t, sizeof(SFTYPE), 1, headerptr);
    fwrite(&tf, sizeof(SFTYPE), 1, headerptr);
    fwrite(&nstep, sizeof(long), 1, headerptr);
    fwrite(&a, sizeof(FTYPE), 1, headerptr);
    fwrite(&gam, sizeof(FTYPE), 1, headerptr);
    fwrite(&cour, sizeof(FTYPE), 1, headerptr);
      
    fwrite(&DTd, sizeof(SFTYPE), 1, headerptr);
    fwrite(&DTavg, sizeof(SFTYPE), 1, headerptr);
    fwrite(&DTener, sizeof(SFTYPE), 1, headerptr);
    fwrite(&DTi, sizeof(SFTYPE), 1, headerptr);
    // fwrite(&DTr, sizeof(SFTYPE),1, headerptr) ;
    fwrite(&DTr, sizeof(long), 1, headerptr);
    fwrite(&DTdebug, sizeof(SFTYPE), 1, headerptr);
    fwrite(&dump_cnt, sizeof(long), 1, headerptr);
    fwrite(&image_cnt, sizeof(long), 1, headerptr);
    fwrite(&rdump_cnt, sizeof(long), 1, headerptr);
    fwrite(&avg_cnt, sizeof(long), 1, headerptr);
    fwrite(&debug_cnt, sizeof(long), 1, headerptr);
      
    fwrite(&dt, sizeof(SFTYPE), 1, headerptr);
    fwrite(&lim, sizeof(int), 1, headerptr);
    fwrite(&TIMEORDER, sizeof(int), 1, headerptr);
    fwrite(&fluxmethod, sizeof(int), 1, headerptr);
    fwrite(&FLUXB, sizeof(int), 1, headerptr);
    fwrite(&UTOPRIMVERSION, sizeof(int), 1, headerptr);
    fwrite(&failed, sizeof(int), 1, headerptr);
      
    fwrite(&R0, sizeof(FTYPE), 1, headerptr);
    fwrite(&Rin, sizeof(FTYPE), 1, headerptr);
    fwrite(&Rout, sizeof(FTYPE), 1, headerptr);
    fwrite(&hslope, sizeof(FTYPE), 1, headerptr);
    fwrite(&defcoord, sizeof(FTYPE), 1, headerptr);

    fwrite(&BCtype[X1UP],sizeof(int), 1, headerptr);
    fwrite(&BCtype[X1DN],sizeof(int), 1, headerptr);
    fwrite(&BCtype[X2UP],sizeof(int), 1, headerptr);
    fwrite(&BCtype[X2DN],sizeof(int), 1, headerptr);
    fwrite(&BCtype[X3UP],sizeof(int), 1, headerptr); // 3d thing
    fwrite(&BCtype[X3DN],sizeof(int), 1, headerptr);

    fwrite(&realnstep, sizeof(long), 1, headerptr);
    fwrite(&debugfail,sizeof(int), 1, headerptr);
    fwrite(&whichrestart,sizeof(int), 1, headerptr);
    fwrite(&cooling,sizeof(int), 1, headerptr);
    fwrite(&restartsteps[0],sizeof(long), 1, headerptr);
    fwrite(&restartsteps[1],sizeof(long), 1, headerptr);
    fwrite(&GAMMIEDUMP,sizeof(int), 1, headerptr);
    fwrite(&GAMMIEIMAGE,sizeof(int), 1, headerptr);
    fwrite(&GAMMIEENER,sizeof(int), 1, headerptr);
    fwrite(&DODIAGS,sizeof(int), 1, headerptr);
    fwrite(&DOENERDIAG,sizeof(int), 1, headerptr);
    fwrite(&DOGDUMPDIAG,sizeof(int), 1, headerptr);
    fwrite(&DORDUMPDIAG,sizeof(int), 1, headerptr);
    fwrite(&DODUMPDIAG,sizeof(int), 1, headerptr);
    fwrite(&DOAVGDIAG,sizeof(int), 1, headerptr);
    fwrite(&DOIMAGEDIAG,sizeof(int), 1, headerptr);
    fwrite(&DOAREAMAPDIAG,sizeof(int), 1, headerptr);
    fwrite(&POSDEFMETRIC,sizeof(int), 1, headerptr);
    fwrite(&periodicx1,sizeof(int), 1, headerptr);
    fwrite(&periodicx2,sizeof(int), 1, headerptr);
    fwrite(&periodicx3,sizeof(int), 1, headerptr); // 3d thing
    fwrite(&binaryoutput,sizeof(int), 1, headerptr);
    fwrite(&sortedoutput,sizeof(int), 1, headerptr);
    fwrite(&defcon,sizeof(FTYPE), 1, headerptr);
    fwrite(&SAFE,sizeof(FTYPE), 1, headerptr);
    fwrite(&RHOMIN,sizeof(FTYPE), 1, headerptr);
    fwrite(&UUMIN,sizeof(FTYPE), 1, headerptr);
    fwrite(&RHOMINLIMIT,sizeof(FTYPE), 1, headerptr);
    fwrite(&UUMINLIMIT,sizeof(FTYPE), 1, headerptr);
    fwrite(&BSQORHOLIMIT,sizeof(FTYPE), 1, headerptr);
    fwrite(&BSQOULIMIT,sizeof(FTYPE), 1, headerptr);
    fwrite(&GAMMAMAX,sizeof(FTYPE), 1, headerptr);
    fwrite(&GAMMADAMP,sizeof(FTYPE), 1, headerptr);
    fwrite(&GAMMAFAIL,sizeof(FTYPE), 1, headerptr);

    PDUMPLOOP(k) fwrite(&prMAX[k],sizeof(FTYPE),1,headerptr);

    fwrite(&rescaletype,sizeof(int),1,headerptr);

    // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
    fwrite(&avgscheme,sizeof(int),1,headerptr);
    fwrite(&do_transverse_flux_integration,sizeof(int),1,headerptr);
    fwrite(&do_source_integration,sizeof(int),1,headerptr);
    fwrite(&do_conserved_integration,sizeof(int),1,headerptr);
    fwrite(&INVERTFROMAVERAGEIFFAILED,sizeof(int),1,headerptr);
    fwrite(&LIMIT_AC_PRIM_FRAC_CHANGE,sizeof(int),1,headerptr);
    fwrite(&MAX_AC_PRIM_FRAC_CHANGE,sizeof(FTYPE),1,headerptr);
    fwrite(&DOENOFLUX,sizeof(int),1,headerptr);


    // new June 6, 2003 (cumulatives)
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(k) fwrite(&pcumreg_tot[enerregion][dir][k],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fwrite(&fladdreg_tot[enerregion][k],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fwrite(&sourceaddreg_tot[enerregion][k],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fwrite(&Ureg_init_tot[enerregion][k],sizeof(FTYPE),1,headerptr);

    TSCALELOOP(tscale) FLOORLOOP(floor) fwrite(&failfloorcountlocal_tot[tscale][floor],sizeof(CTYPE),1,headerptr);
    // end new June 6,2003

    if(DODISS) ENERREGIONLOOP(enerregion) fwrite(&dissreg_tot[enerregion][0],sizeof(FTYPE),1,headerptr);

    if(DOLUMVSR) fwrite(lumvsr_tot,sizeof(SFTYPE),ncpux1*N1,headerptr);
    if(DODISSVSR) fwrite(dissvsr_tot,sizeof(SFTYPE),ncpux1*N1,headerptr);


  }
  else if(bintxt==TEXTOUTPUT){
    fprintf(headerptr,WRITERESTARTHEADER,
	    totalsize[1],totalsize[2], totalsize[3], //3
	    t,tf,nstep,a,gam,cour, // 6
	    DTd,DTavg,DTener,DTi,DTr,DTdebug,dump_cnt,image_cnt,rdump_cnt,avg_cnt,debug_cnt, // 11
	    dt,lim,TIMEORDER,fluxmethod,FLUXB,UTOPRIMVERSION,failed,  // 3
	    R0,Rin,Rout,hslope,defcoord, // 5
	    BCtype[X1UP],BCtype[X1DN],BCtype[X2UP],BCtype[X2DN],BCtype[X3UP],BCtype[X3DN], // 6
	    realnstep,debugfail,whichrestart,cooling,restartsteps[0],restartsteps[1],GAMMIEDUMP,GAMMIEIMAGE,GAMMIEENER,DODIAGS,DOENERDIAG,DOGDUMPDIAG,DORDUMPDIAG,DODUMPDIAG,DOAVGDIAG,DOIMAGEDIAG,DOAREAMAPDIAG,POSDEFMETRIC,periodicx1,periodicx2,periodicx3,binaryoutput,sortedoutput,defcon,SAFE,RHOMIN,UUMIN,RHOMINLIMIT,UUMINLIMIT,BSQORHOLIMIT,BSQOULIMIT,GAMMAMAX,GAMMADAMP,GAMMAFAIL // 23+10+1
	    // total = 65 (not updated)
	    );

    // end new June 6,2003
    PDUMPLOOP(k) fprintf(headerptr,HEADERONEOUT,prMAX[k]);
    fprintf(headerptr,"%d ",rescaletype);


    // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
    fprintf(headerptr,"%d ",avgscheme);
    fprintf(headerptr,"%d ",do_transverse_flux_integration);
    fprintf(headerptr,"%d ",do_source_integration);
    fprintf(headerptr,"%d ",do_conserved_integration);
    fprintf(headerptr,"%d ",INVERTFROMAVERAGEIFFAILED);
    fprintf(headerptr,"%d ",LIMIT_AC_PRIM_FRAC_CHANGE);
    fprintf(headerptr,HEADERONEOUT,MAX_AC_PRIM_FRAC_CHANGE);
    fprintf(headerptr,"%d ",DOENOFLUX);


      // new June 6, 2003 (cumulatives)
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(k) fprintf(headerptr,HEADERONEOUT,pcumreg_tot[enerregion][dir][k]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fprintf(headerptr,HEADERONEOUT,fladdreg_tot[enerregion][k]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fprintf(headerptr,HEADERONEOUT,sourceaddreg_tot[enerregion][k]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fprintf(headerptr,HEADERONEOUT,Ureg_init_tot[enerregion][k]);

#if(COUNTTYPE==LONGLONGINTTYPE)
    TSCALELOOP(tscale) FLOORLOOP(floor) fprintf(headerptr,"%lld ",failfloorcountlocal_tot[tscale][floor]);
#elif(COUNTTYPE==DOUBLETYPE)
    TSCALELOOP(tscale) FLOORLOOP(floor) fprintf(headerptr,"%21.15g ",failfloorcountlocal_tot[tscale][floor]);
#endif
    
    if(DODISS) ENERREGIONLOOP(enerregion) fprintf(headerptr,HEADERONEOUT,dissreg_tot[enerregion][0]);

    // assumes same CPU geometry during restart
    if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) fprintf(headerptr,HEADERONEOUT,lumvsr_tot[ii]);

    // assumes same CPU geometry during restart
    if(DODISSVSR) for(ii=0;ii<ncpux1*N1;ii++) fprintf(headerptr,HEADERONEOUT,dissvsr_tot[ii]);


    fprintf(headerptr,"\n");
  }
  fflush(headerptr);


  trifprintf("end writing header of restart file\n");
  return(0);
}
