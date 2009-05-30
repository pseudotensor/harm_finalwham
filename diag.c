
#include "decs.h"

/* diagnostics subroutine */

#define DIAGREPORT {trifprintf("t=%21.15g to do: tener=%21.15g (dt=%21.15g): dump_cnt=%ld @ t=%21.15g (dt=%21.15g) : avg_cnt=%ld @ t=%21.15g (dt=%21.15g) : debug_cnt=%ld @ t=%21.15g (dt=%21.15g) : image_cnt=%ld @ t=%21.15g (dt=%21.15g): restart=%d @ nstep=%ld (dt=%ld)\n",t,tener,DTener,dump_cnt,tdump,DTd,avg_cnt,tavg,DTavg,debug_cnt,tdebug,DTdebug,image_cnt,timage,DTi,whichrestart,nrestart,DTr);}


int diag(int call_code)
{
  //////////////////
  //
  // BEGIN DIAG VARS
  //
  static int firsttime = 1;
	int pl;
	FTYPE asym[NPR], norm[NPR], maxasym[NPR];

  FILE *imagecnt_file, *dumpcnt_file,*avgcnt_file,*debugcnt_file,*fieldlinecnt_file;

  static SFTYPE tlastener,tlastdump,tlastimage,tlastareamap,tlastavg,tlastdebug,tlastfieldline;
  static int dumpc, avgc, imagec, restartc, enerc,debugc,fieldlinec;
  static SFTYPE tdump, tavg, timage, tener,tdebug,tfieldline;
  static long nlastrestart,nrestart;
  int dogdump,doavg, dordump,dodump,doener,doimagedump,doareamap,dodebug,dofieldlinedump;

  int floor,i,j,k;
  int dir,interpi,enodebugi;
  int enodebugdump(long dump_cnt);
  int asym_compute_1(FTYPE (*prim)[N2M][N3M][NPR]);


  ///////////////////////////
  //
  // setup timing of writes to files
  //

  if ((call_code == INIT_OUT) || (firsttime == 1)) {
    
    // no need to repeat if restarting
    tlastener  = DTener*(SFTYPE)((int)(t/DTener))-SMALL;
    nlastrestart = (long) (DTr*(SFTYPE)((realnstep/DTr))-1);
    tlastdump  = DTd*(SFTYPE)((int)(t/DTd))-SMALL;
    tlastavg  = DTavg*(SFTYPE)((int)(t/DTavg))-SMALL;
    tlastimage = DTi*(SFTYPE)((int)(t/DTi))-SMALL;
    tlastdebug  = DTdebug*(SFTYPE)((int)(t/DTdebug))-SMALL;


    tlastareamap = t-SMALL;
    
    dumpc = avgc = imagec = restartc = enerc = debugc = 0;

    if (RESTARTMODE == 0) {
      tdump = timage = tener = t;
      tavg=t+DTavg; // do next time
      nrestart = nstep;

      dump_cnt = 0;
      image_cnt = 0;
      rdump_cnt = 0;
      avg_cnt = 0;
      debug_cnt = 0;

      appendold = 0;
    } else {
      setrestart(&appendold);

      // assuming started at t=0 and nstep=0 for original run
      // time to dump NEXT
      tener  = DTener*(SFTYPE)((int)(t/DTener)+1);
      nrestart = (long)(DTr*(SFTYPE)((realnstep/DTr)+1));
      tdump  = DTd*(SFTYPE)((int)(t/DTd)+1);
      timage = DTi*(SFTYPE)((int)(t/DTi)+1);
      tavg  = DTavg*(SFTYPE)((int)(t/DTavg)+1);
      tdebug  = DTdebug*(SFTYPE)((int)(t/DTdebug)+1);
      DIAGREPORT;
    }
  }



  ///////////////////////
  //
  // setup what we will dump this call to diag()
  //

  // output grid (probaly want both fullgrid (to make sure ok) and compute grid to compare with data dumps
  if((DOGDUMPDIAG)&&(!GAMMIEDUMP)&&(firsttime&&(RESTARTMODE==0))){
    dogdump=1;
  }
  else dogdump=0;


  if( ((DODUMPDIAG)&&(DODIAGEVERYSUBSTEP||((t!=tlastdump)&&(t >= tdump || (RESTARTMODE&&dofaildump&&(nstep>=steptofaildump)) || call_code==FINAL_OUT ))) )  ){
    dodump=1;
  }
  else dodump=0;

  if( ((DODEBUG)&&(DOENODEBUGEVERYSUBSTEP||DODIAGEVERYSUBSTEP||((t!=tlastdebug)&&(t >= tdebug || call_code==FINAL_OUT))) )  ){
    dodebug=1;
  }
  else dodebug=0;

  
  if((DOAVGDIAG)&&((t!=tlastavg)&&(t >= tavg))){
    doavg=1;
  }
  else doavg=0;
    
  if((DORDUMPDIAG)&&( ((nlastrestart!=nrestart)&&(failed == 0) && (realnstep >= nrestart ))||(call_code==FINAL_OUT) ) ){
    dordump=1;
  }
  else dordump=0;

  // t!=tlast avoids duplicate entries
  if(DOENERDIAG&&(DODIAGEVERYSUBSTEP||((t!=tlastener)&&( (t >= tener)||(call_code==INIT_OUT)||(call_code==FINAL_OUT)||firsttime)))){
    doener=1;
  }
  else doener=0;

  if((DOAREAMAPDIAG)&&(t != tlastareamap)&&(dofailmap)&&(nstep>=steptofailmap)){
    doareamap=1;
  }
  else doareamap=0;

  /* image dump at regular intervals */
  if((DOIMAGEDIAG)&&(DODIAGEVERYSUBSTEP||((t!=tlastimage)&&(t >= timage))) ){
    doimagedump=1;
  }
  else doimagedump=0;

  /* fieldline dump at regular intervals */
  // use image time period
  if((DOFIELDLINE)&&(DODIAGEVERYSUBSTEP||((t!=tlastfieldline)&&(t >= tfieldline))) ){
    dofieldlinedump=1;
    fieldline_cnt=image_cnt; // force to go with image dump (needed also so restart() knows the count)
  }
  else dofieldlinedump=0;

  //////////////////////////////////////////////////////////////////////
  //
  ////////////////// NOW WRITE TO FILES
  //
  //////////////////////////////////////////////////////////////////////


#if(ASYMDIAGCHECK)
  asym_compute_1(pdump);
#endif



  //////////////////////
  //
  // Grid dump
  //
  /////////////////////

  if(dogdump) gdump();



  // if doing simulbccalc type loop in step_ch.c then need to bound when doing diagnostics since not done yet
  if(SIMULBCCALC>=1){
    if(dodump||dordump||doener){
      // for dump, rdump, and divb in ener
      bound_prim(-1,p);
    }
  }


  //////////////////////
  //
  // ener dump (integrated quantities: integrate and possibly  dump them too (if doener==1))
  //
  /////////////////////

  if(doener||dordump||(call_code==FINAL_OUT)||(call_code==INIT_OUT)){ // need integratd quantities for restart dump, but don't write them to ener file.

    dump_ener(doener,dordump,call_code);
    if(doener||(call_code==FINAL_OUT)||(call_code==INIT_OUT)){
      frdotout(); // need to include all terms and theta fluxes on horizon/outer edge at some point GODMARK
    }
    while (t >= tener) {
      enerc++;
      tener = enerc * DTener;
    }
    tlastener=t;
  }




  ///////////////////////
  //
  // RESTART DUMP
  //
  ///////////////////////
  

  if(dordump){
    DIAGREPORT;
    trifprintf("dumping: restart: %d nstep: %ld nlastrestart: %ld nrestart: %ld restartc: %d\n", whichrestart,realnstep,nlastrestart,nrestart,restartc);
    restart_write((long)whichrestart);	// 0 1 0 1 0 1 ...
    restartsteps[whichrestart] = realnstep;
    whichrestart = !whichrestart;
    while (realnstep >= nrestart) {
      restartc++;
      nrestart = restartc * DTr;
    }
    nlastrestart=realnstep;
  }
      
  
  ///////////////////////
  //
  // DUMP
  //
  ///////////////////////
  
  if(dodump){
    DIAGREPORT;
    trifprintf("dumping: dump_cnt=%ld : t=%21.15g tlastdump=%21.15g tdump=%21.15g dumpc=%d\n", dump_cnt,t,tlastdump,tdump,dumpc);
    
    /* make regular dump file */
    if (dump(dump_cnt) >= 1){
      dualfprintf(fail_file,"unable to print dump file\n");
      return (1);
    }
    restart_write(-(long)dump_cnt-1);// so can restart at a dump without reconstructing the rdump from a dump.


    if(DODISS){
      /* make dissdump file */
      if (dissdump(dump_cnt) >= 1){
	dualfprintf(fail_file,"unable to print dissdump file\n");
	return (1);
      }
    }


    // iterate counter
    dump_cnt++;
    while (t >= tdump) {
      dumpc++;
      tdump = dumpc * DTd;
    }
    // output number of dumps
    myfopen("dumps/0_numdumps.dat","w","error opening dump count file\n",&dumpcnt_file);      
    myfprintf(dumpcnt_file, "# Number of dumps\n%ld\n", dump_cnt);
    myfclose(&dumpcnt_file,"Couldn't close dumpcnt_file");

    tlastdump=t;
  }

  ///////////////////////
  //
  // DEBUG DUMP
  //
  ///////////////////////
  
  if(dodebug){
    DIAGREPORT;
    trifprintf("debug dumping: debug_cnt=%ld : t=%21.15g tlastdebug=%21.15g tdebug=%21.15g debugc=%d\n", debug_cnt,t,tlastdebug,tdebug,debugc);
    
    /* make regular dump file */
    if (debugdump(debug_cnt) >= 1){
      dualfprintf(fail_file,"unable to print debug dump file\n");
      return (1);
    }

    if(DOENODEBUG){
      /* make regular dump file */
      if (enodebugdump(debug_cnt) >= 1){
	dualfprintf(fail_file,"unable to print enodebug dump file\n");
	return (1);
      }
    }

    // iterate counter
    debug_cnt++;
    while (t >= tdebug) {
      debugc++;
      tdebug = debugc * DTdebug;
    }
    // output number of dumps
    myfopen("dumps/0_numdebug.dat","w","error opening debug dump count file\n",&debugcnt_file);      
    myfprintf(debugcnt_file, "# Number of debug dumps\n%ld\n", debug_cnt);
    myfclose(&debugcnt_file,"Couldn't close debugcnt_file");

    tlastdebug=t;
  }


  ///////////////////////
  //
  // AVG
  //
  ///////////////////////

  if(DOAVGDIAG){
    // do every time step
    // assume can't fail, but can apparently
    if(average_calc(doavg)>=1) return(1);
  }
  
  if(doavg){
    DIAGREPORT
    trifprintf("avging dump: avg_cnt=%ld : t=%21.15g tlastavg=%21.15g tavg=%21.15g avgc=%d\n", avg_cnt,t,tlastavg,tavg,avgc);
    
    /* make avg dump file */
    if (avgdump(avg_cnt) >= 1){
      dualfprintf(fail_file,"unable to print avg file\n");
      return (1);
    }
#if(DOAVG2)
    /* make avg dump file */
    if (avgdump2(avg_cnt) >= 1){
      dualfprintf(fail_file,"unable to print avg2 file\n");
      return (1);
    }
#endif

    // iterate counter
    avg_cnt++;
    while (t >= tavg) {
      avgc++;
      tavg = avgc * DTavg;
    }
    // output number of avgs
    myfopen("dumps/0_numavgs.dat","w","error opening avg count file\n",&avgcnt_file);      
    myfprintf(avgcnt_file, "# Number of avgs\n%ld\n", avg_cnt);
    myfclose(&avgcnt_file,"Couldn't close avgcnt_file");

    tlastavg=t;
  }

  ///////////////////////
  //
  // AREA MAP
  //
  ///////////////////////
  if(doareamap){
    if(area_map(call_code, TIMESERIESAREAMAP, 20, ifail, jfail, kfail, pdump)>=1) return(1);
    tlastareamap=t;
  }

  ///////////////////////
  //
  // IMAGE
  //
  ///////////////////////

  if(doimagedump){
    DIAGREPORT;
    trifprintf("image dump %ld : t=%21.15g tlastimage=%21.15g timage=%21.15g imagec=%d\n", image_cnt, t,tlastimage,timage,imagec);
    
    
    /* make regular image file */
    if(image_dump(image_cnt)>=1) return(1);

    // iterate counter
    image_cnt++;
    while (t >= timage) {
      imagec++;
      timage = imagec * DTi;
    }

    // output number of images
    myfopen("images/0_numimages.dat","w","error opening image count file\n",&imagecnt_file);      
    myfprintf(imagecnt_file, "# Number of images\n%ld\n", image_cnt);
    myfclose(&imagecnt_file,"Couldn't close imagecnt_file");
    
    tlastimage=t;
  }


  ///////////////////////
  //
  // FIELDLINE
  //
  ///////////////////////


  if(dofieldlinedump){
    DIAGREPORT;
    trifprintf("fieldline dump %ld : t=%21.15g tlastfieldline=%21.15g tfieldline=%21.15g fieldlinec=%d\n", fieldline_cnt, t,tlastfieldline,tfieldline,fieldlinec);
    
    // (after processing) equivalent to image in interest
    if(fieldlinedump(fieldline_cnt)>=1) return(1);

    // iterate counter
    fieldline_cnt++;
    while (t >= tfieldline) {
      fieldlinec++;
      tfieldline = fieldlinec * DTi;
    }

    // output number of fieldlines
    myfopen("dumps/0_numfieldlines.dat","w","error opening fieldline count file\n",&fieldlinecnt_file);      
    myfprintf(fieldlinecnt_file, "# Number of fieldlines\n%ld\n", fieldline_cnt);
    myfclose(&fieldlinecnt_file,"Couldn't close fieldlinecnt_file");
    
    tlastfieldline=t;
  }


  ///////////////////////
  //
  // DIAG clearings
  //
  // finally some clearing based upon what called
  //
  ////////////////////////

  if(DODEBUG){ // shouldn't clear these till after ener and debug dump done so both have all timescale data.
    if(doener){
      // cleanse the ener time scale for the failure diag
      ZLOOP FLOORLOOP(floor) failfloorcount[i][j][k][ENERTS][floor]=0;
    }
    if(dodebug){
      // clense failure diag
      ZLOOP FLOORLOOP(floor) failfloorcount[i][j][k][DEBUGTS][floor]=0;
    }
    if(doimagedump){
      // clense the failure counts for this time scale
      ZLOOP FLOORLOOP(floor) failfloorcount[i][j][k][IMAGETS][floor]=0;
    }
  }

  if(DOENODEBUG){
    if(dodebug){
      FULLLOOP DIMENLOOP(dir) INTERPLOOP(interpi) PLOOP(pl) ENODEBUGLOOP(enodebugi){
	if(dir<=2 && pl<=U2){
	  enodebugarray[i][j][k][dir-1][interpi][pl][enodebugi]=0;
	}
      }
    }
  }

  if(DODISS){
    // clear diss
    if(dodump){
      FULLLOOP{
	dissfunpos[i][j][k][0]=0.0;
	dissfunpos[i][j][k][1]=0.0; // clear failures too
      }
    }
  }


/*	//symmetry check
	j = 0;
	k = 0;
	PLOOP(pl) maxasym[pl] = 0.0;

	PLOOP(pl) {
		if( pl >= U1 && pl <= U3 ) {
			norm[pl] = coordparams.timescalefactor;
		}
		else if( pl == UU ) {
			norm[pl] = coordparams.timescalefactor * coordparams.timescalefactor;
		}
		else {
			norm[pl] = 1.;
		}

		norm[pl] = 1 / norm[pl];
	}
	
	for( i = 0; i < N1; i++ ) {
		PLOOP( pl ) {
			if( pl == U1 ) {
				asym[pl] = fabs( p[i][j][k][pl] + p[N1 - 1 - i][j][k][pl] );
			}
			else {
				asym[pl] = fabs( p[i][j][k][pl] - p[N1 - 1 - i][j][k][pl] );
			}


			asym[pl] /= norm[pl];

			maxasym[pl] = MAX( asym[pl], maxasym[pl] );
		}
	}

	dualfprintf( fail_file, "asym %ld ", realnstep );
		
	PLOOP(pl) if( pl <= U1 ) dualfprintf( fail_file, " %22.16g", maxasym[pl] );

	dualfprintf( fail_file, "\n" );
*/
  firsttime = 0;
  return (0);
}

/** some diagnostic routines **/





// for implosion problem
int asym_compute_1(FTYPE (*prim)[N2M][N3M][NPR])
{
  int i,j,k;

  FULLLOOP{

    if(prim[i][j][k][RHO]!=prim[j][i][k][RHO]){
      dualfprintf(fail_file,"ASYM in RHO %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][RHO],prim[j][i][k][RHO]);
    }

    if(prim[i][j][k][UU]!=prim[j][i][k][UU]){
      dualfprintf(fail_file,"ASYM in UU %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][UU],prim[j][i][k][UU]);
    }


    if(prim[i][j][k][U1]!=prim[j][i][k][U2]){
      dualfprintf(fail_file,"ASYM in U1 %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][U1],prim[j][i][k][U2]);
    }


      
  }


  return(0);
}

// for implosion problem
int asym_compute_2(FTYPE (*prim)[N2M][N3M][NPR])
{
  int i,j,k;

  ZLOOP{

    if(prim[i][j][k][RHO]!=prim[j][i][k][RHO]){
      dualfprintf(fail_file,"ASYM in RHO %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][RHO],prim[j][i][k][RHO]);
    }

    if(prim[i][j][k][UU]!=prim[j][i][k][UU]){
      dualfprintf(fail_file,"ASYM in UU %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][UU],prim[j][i][k][UU]);
    }


    if(prim[i][j][k][U1]!=prim[j][i][k][U2]){
      dualfprintf(fail_file,"ASYM in U1 %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][U1],prim[j][i][k][U2]);
    }


      
  }


  return(0);
}



// 2D only for now since really only useful for 2D imaging

/* map out region around failure point */
int area_map(int call_code, int type, int size, int i, int j, int k, FTYPE prim[][N2M][N3M][NPR])
{
  int pl;
  int l,m,ll,mm;
  FTYPE vmin1, vmax1, vmin2, vmax2;
  int ignorecourant;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM],V[NDIM];
  FTYPE divb;
  FTYPE tens_em[NDIM][NDIM], tens_matt[NDIM][NDIM], b[NDIM],
    ucon[NDIM];
  FTYPE U[NPR];
  int lowersizex1,uppersizex1;
  int lowersizex2,uppersizex2;

  static FILE* fileptr;
  static int firsttime=1;
  static int domap=0;
  static int doclose=0;



  trifprintf("\nStart area_map function ... ");


  k=N3/2+SHIFT3; // middle cell // GODMARK : 2D only below


  if(i-(-N1BND)<size/2) lowersizex1=i-(-N1BND);
  else lowersizex1=size/2;
  if((N1-1+N1BND)-i<size/2) uppersizex1=(N1-1+N1BND)-i;
  else uppersizex1=size/2;
  
  if(j-(-N2BND)<size/2) lowersizex2=j-(-N2BND);
  else lowersizex2=size/2;
  if((N2-1+N2BND)-j<size/2) uppersizex2=(N2-1+N2BND)-j;
  else uppersizex2=size/2;


  if(firsttime){
    if((type==TIMESERIESAREAMAP)&&(dofailmap)){
      if((fileptr=fopen("areamap","wt"))==NULL){
	dualfprintf(fail_file,"Cannot open ./areamap on proc=%d\n",myid);
	domap=0;
      }
      else domap=1;
    }
  }

  if((type==TIMESERIESAREAMAP)&&domap&&(call_code==2)){
    doclose=1;
  }
  else doclose=0;

  if(type==FINALTDUMPAREAMAP){
    dualfprintf(fail_file, "area map\n");
    dualfprintf(fail_file, "failure at: i=%d j=%d k=%d\n",i+startpos[1],j+startpos[2],k+startpos[3]);
    coord(i,j,k,CENT,X);
    dualfprintf(fail_file, "failure at: i=%d j=%d k=%d\n",i+startpos[1],j+startpos[2],k+startpos[3]);



    PLOOP(pl) {// all vars
      
      dualfprintf(fail_file, "variable %d \n", pl);
      
      dualfprintf(fail_file, "i = \t ");
      for(l=i-lowersizex1;l<=i+uppersizex1;l++){
	ll=l+startpos[1];
	dualfprintf(fail_file, "%21d", ll);
      }
      dualfprintf(fail_file, "\n");
      for(m=j-lowersizex2;m<=j+uppersizex2;m++){
	mm=m+startpos[2];
	dualfprintf(fail_file, "j = %d \t ",mm);
	for(l=i-lowersizex1;l<=i+lowersizex1;l++){
	  ll=l+startpos[1];
	  dualfprintf(fail_file, "%21.15g ",prim[l][m][k][pl]);
	}
	dualfprintf(fail_file, "\n");
      }
    }
  }
  else if((type==TIMESERIESAREAMAP)&&(domap)){
    if(firsttime){
      // GODMARK: 2D only, not outputting z-related stuff so function remains consistent with SM macro
      fprintf(fileptr,"%21.15g %d %d %21.15g %21.15g %21.15g %21.15g %d %d %d %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
	      t,totalsize[1],totalsize[2],startx[1],startx[2],dx[1],dx[2],lowersizex1,uppersizex1,lowersizex2,uppersizex2,startpos[1]+i,startpos[2]+j,gam,a,R0,Rin,Rout,hslope);
      fflush(fileptr);
    }
    for(m=j-size/2;m<=j+size/2;m++){
      if((m<-N2BND)||(m>N2-1+N2BND)) continue;
      mm=m+startpos[2];
      for(l=i-size/2;l<=i+size/2;l++){
	if((l<-N1BND)||(l>N1-1+N1BND)) continue;
	ll=l+startpos[1];
	
	
	coord(l, m, k, CENT, X);
	bl_coord(X, V); 
	get_geometry(l, m, k, CENT, &geom);
	if (!failed) {
	  if (get_state(p[l][m][k], &geom, &q) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "get_state() dir=0", 1);
	  if (vchar(p[l][m][k], &q, 1, &geom, &vmax1, &vmin1,&ignorecourant) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "vchar() dir=1or2", 1);
	  if (vchar(p[l][m][k], &q, 2, &geom, &vmax2, &vmin2,&ignorecourant) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "vchar() dir=1or2", 2);
	  // GODMARK: no 3-direction char.
	}
	if((l>=-1)&&(l<=N1+1)&&(m>=-1)&&(m<=N2+1)&&(k>=-1)&&(k<=N3+1) ){ SETFDIVB(divb, p, l, m, k);}
	else divb=0.0;

	// same order as dump.c for first columns (easy sm read)
	fprintf(fileptr,
		"%d %d "
		"%21.15g %21.15g "
		"%21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "
		"%21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g "
		"%21.15g %ld\n",
		ll,mm,
		X[1],X[2],
		V[1],V[2],
		p[l][m][k][0],
		p[l][m][k][1],
		p[l][m][k][2],
		p[l][m][k][3],
		p[l][m][k][4],
		p[l][m][k][5],
		p[l][m][k][6],
		p[l][m][k][7],
		divb,
		q.ucon[0],q.ucon[1],q.ucon[2],q.ucon[3],
		q.ucov[0],q.ucov[1],q.ucov[2],q.ucov[3],
		q.bcon[0],q.bcon[1],q.bcon[2],q.bcon[3],
		q.bcov[0],q.bcov[1],q.bcov[2],q.bcov[3],
		vmin1,vmax1,vmin2,vmax2,
		geom.g,
		t,realnstep);
      }
    }
    fflush(fileptr);
  }

  if(doclose) if(fileptr!=NULL) fclose(fileptr);


  /* print out other diagnostics here */

  firsttime=0;
  trifprintf("end area_map function.\n");  
  return(0);
}



/* evaluate fluxed based diagnostics; put results in global variables */

#define JETBSQORHO (3.162)


// notice that F1 and F2 have arbitrary eomfunc that must be divided out to get real flux
int diag_flux(FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR], SFTYPE Dt)
{
  int fluxdir;
  int i, j, k, pl, l, dir,fl,enerregion;
  SFTYPE surface,surf2;
  SFTYPE surgdet;
  FTYPE (*flux)[N2M][N3M][NPR];
  int start1,start2,start3,stop1,stop2,stop3;
  int gpos;
  int ii;
  FTYPE ftemp;
  FTYPE ftemp0,ftemp1,ftemp2,ftemp3,ftemp4,ftemp5,ftemp6;
  FTYPE pgas,bsq,bsqorho;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM],V[NDIM];
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];
  int condjet2;
  FTYPE Ftemp[NPR],Ftempdiag[NPR];
  int firstinloop;
  FTYPE pr[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);



  // initialize
  ENERREGIONLOOP(enerregion){
    doflux=dofluxreg[enerregion];
    enerpos=enerposreg[enerregion];
    pcum=pcumreg[enerregion];
    pdot=pdotreg[enerregion];
    pdotterms=pdottermsreg[enerregion];
    

    DIRLOOP(dir) PLOOP(pl){
      pdot[dir][pl]=0;
      FLLOOP(fl) pdotterms[dir][fl][pl]=0;
      if(enerregion==0) FLLOOP(fl) pdottermsjet2[dir][fl][pl]=0;
    }
    // true outer boundary surface fluxes (per direction, per conserved variable)
    DIRLOOP(dir) {
      if (doflux[dir] >= 0) {// then this cpu is involved in flux BOUNDARY calculation
	// otherwise don't add to it


	///////
	// assumes rectangular region
	//
	if((dir==X1UP)||(dir==X1DN)){

	  surface = dx[2]*dx[3];

	  start1=stop1=doflux[dir];

	  start2=enerpos[X2DN];
	  stop2=enerpos[X2UP];

	  start3=enerpos[X3DN];
	  stop3=enerpos[X3UP];

	  flux=F1;
	  gpos=FACE1;
	  fluxdir=1;
	}
	else if((dir==X2UP)||(dir==X2DN)){

	  surface = dx[1]*dx[3];

	  start2=stop2=doflux[dir];

	  start1=enerpos[X1DN];
	  stop1=enerpos[X1UP];

	  start3=enerpos[X3DN];
	  stop3=enerpos[X3UP];

	  flux=F2;
	  gpos=FACE2;
	  fluxdir=2;
	}
	else if((dir==X3UP)||(dir==X3DN)){

	  surface = dx[1]*dx[2];

	  start3=stop3=doflux[dir];

	  start1=enerpos[X1DN];
	  stop1=enerpos[X1UP];

	  start2=enerpos[X2DN];
	  stop2=enerpos[X2UP];

	  flux=F3;
	  gpos=FACE3;
	  fluxdir=3;
	}
	else{
	  dualfprintf(fail_file,"no such direction: %d\n",dir);
	  myexit(1);
	}

	// zero out summation quantities
	PLOOP(pl){
	  pdot[dir][pl]=0;
	  FLLOOP(fl) pdotterms[dir][fl][pl]=0;
	  if(enerregion==0) FLLOOP(fl) pdottermsjet2[dir][fl][pl]=0;
	}
	GENLOOP(i,j,k,start1,stop1,start2,stop2,start3,stop3){
	  // now add up all zones for each conserved quantity (PLOOP)

	  ////////////
	  //
	  // do standard flux for accounting
	  //
	  ////////////

	  get_geometry(i, j, k, gpos, &geom);
	  PLOOP(pl) Ftemp[pl]=flux[i][j][k][pl]*surface; // in UEVOLVE form
	  // GODMARK: for finite volume method, below doesn't change the result if eomfunc=gdet.
	  // Otherwise flux would have to be completely recomputed for gdet case JUST for diagnostic to be consistent at higher order
	  UtoU(UEVOLVE,UDIAG,&geom,Ftemp,Ftempdiag); // convert to diag form
	  PLOOP(pl) pdot[dir][pl]  += Ftempdiag[pl];

	  // DEBUG
	  //	  PLOOP(pl) if(!finite(pdot[dir][pl])){
	  //	    dualfprintf(fail_file,"not finite: i=%d j=%d k=%d :: dir=%d pl=%d %g : %g\n",i,j,k,dir,pl,pdot[dir][pl],flux[i][j][k][pl]);
	  //	  }


	
	  ////////////
	  //
	  // do term-by-term flux for physics accounting
	  // somewhat like dumps and somewhat like avg of terms
	  // GODMARK: For finite volume method this does not differentiate between point and average values!
	  //
	  ////////////
	  get_geometry(i, j, k, CENT, &geom);
	  PLOOP(pl) pr[pl]=prim[i][j][k][pl];

	  //coord(i, j, k, CENT, X);
	  //bl_coord(X, V);
	    // if failed, then data output for below invalid, but columns still must exist    
	    // CENT since p at center
	  if (!failed) {
	    if (get_state(pr, &geom, &q) >= 1)
	      FAILSTATEMENT("diag.c:diag_flux()", "get_state() dir=0", 1);
	  }
	  else {// do a per zone check, otherwise set to 0
	    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
	    if (get_state(pr, &geom, &q) >= 1){
	      for (l = 0; l < NDIM; l++)
		q.ucon[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.ucov[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.bcon[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.bcov[l]=0;
	    }
	    whocalleducon=0; // return to normal state
	  }
	  
	  // somewhat like function which averages stress terms
	  pgas = pressure_rho0_u(pr[RHO],pr[UU]);
	  bsq=0; for(l=0;l<NDIM;l++) bsq+=(q.bcon[l])*(q.bcov[l]);
	  if(enerregion==0){
	    bsqorho=bsq/pr[RHO]; // b^2/\rho
	    // we assume user will check if this condition makes sense for a particular simulation.
	    // condition answer for recording for jet2 region
	    //
	    // a real analysis would trace a field line from the horizon to the outer edge in the funnel-type region and only include the region inside, where eventually we have at the outer edge (or will have) unbound/outbound flow.
	    if(dir==X1DN){
	      condjet2=(bsqorho>JETBSQORHO); // assumes that plunging region never develops such large values.  Can occur, but generally not so.  Can raise if problems.
	    }
	    else if(dir==X1UP){ // assumes in jet2 region, but non-jet2 region could have this property.
	      condjet2=((q.ucon[1]>0.0)&&(-q.ucov[0]-1.0>0.0)); // outgoing and unbound at outer edge
	    }
	    else{
	      condjet2=0;
	    }
	  }
	  else condjet2=0; // never touches pdottermsjet2 then
	  
	  
	  surgdet=(geom.g)*surface;

	  // loop and if's since some calculations are redundantly simple for similar pl
	  PLOOP(pl){
	    if(pl==RHO){
	      ftemp0=pr[pl]*(q.ucon[fluxdir])*surgdet;
	      pdotterms[dir][0][pl]+=ftemp0; // only one part
	      if(condjet2) pdottermsjet2[dir][0][pl]+=ftemp0; // only one part
	      //	  pdot[dir][pl]  += ftemp0 * surgdet;
	    }
	    // part0-6
	    else if((pl>=UU)&&(pl<=U3)){
	      l=pl-UU;
	      // we currently DO NOT add flux[RHO] to flux[UU], just assume reader knows this is from the native stress
	      ftemp0=pgas*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      pdotterms[dir][0][pl]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][pl]+=ftemp0;

	      ftemp1=p[i][j][k][RHO]*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      pdotterms[dir][1][pl]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][pl]+=ftemp1;


	      ftemp2=p[i][j][k][UU]*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      pdotterms[dir][2][pl]+=ftemp2;
	      if(condjet2)	    pdottermsjet2[dir][2][pl]+=ftemp2;


	      ftemp3=bsq*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      pdotterms[dir][3][pl]+=ftemp3;
	      if(condjet2)	    pdottermsjet2[dir][3][pl]+=ftemp3;


	      ftemp4=pgas*delta(fluxdir,pl-UU)*surgdet;
	      pdotterms[dir][4][pl]+=ftemp4;
	      if(condjet2)	    pdottermsjet2[dir][4][pl]+=ftemp4;


	      ftemp5=0.5*bsq*delta(fluxdir,pl-UU)*surgdet;
	      pdotterms[dir][5][pl]+=ftemp5;
	      if(condjet2)	    pdottermsjet2[dir][5][pl]+=ftemp5;


	      ftemp6=-(q.bcon[fluxdir])*(q.bcov[l])*surgdet;
	      pdotterms[dir][6][pl]+=ftemp6;
	      if(condjet2)	    pdottermsjet2[dir][6][pl]+=ftemp6;

	    }
	    else if(pl==B1){
	      ftemp0=(q.bcon[1])*(q.ucon[fluxdir])*surgdet; // flux_b1 term1
	      pdotterms[dir][0][pl]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][pl]+=ftemp0;


	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[1])*surgdet; // flux_b1 term2
	      pdotterms[dir][1][pl]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][pl]+=ftemp1;

	    }
	    else if(pl==B2){
	      ftemp0=(q.bcon[2])*(q.ucon[fluxdir])*surgdet; // flux_b2 term1
	      pdotterms[dir][0][pl]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][pl]+=ftemp0;

	    
	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[2])*surgdet; // flux_b2 term2
	      pdotterms[dir][1][pl]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][pl]+=ftemp1;

	    }
	    else if(pl==B3){
	      ftemp0=(q.bcon[3])*(q.ucon[fluxdir])*surgdet; // flux_b3 term1
	      pdotterms[dir][0][pl]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][pl]+=ftemp0;

	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[3])*surgdet; // flux_b3 term2
	      pdotterms[dir][1][pl]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][pl]+=ftemp1;

	    }

	  }// end PLOOP over term-by-term fluxes

	}// end GENLOOP

	// cumulative only based upon pdot
	// pdot contains entire sum over grid of relevant surface integral value for this direction and ener-region.
	PLOOP(pl) pcum[dir][pl]+=pdot[dir][pl]*Dt;

      }// end if doflux
    }// end DIRLOOP
  }// end ENERloop
#if(COMPUTEFRDOT)
  // radial flux vs. radius
  flux=F1;
  surface = dx[2]*dx[3];
  for(i=0;i<N1;i++){

    PLOOP(pl) frdot[i][pl]=0;

    for(j=0;j<N2;j++) for(k=0;k<N3;k++){
      get_geometry(i, j, k, FACE1, &geom);
      PLOOP(pl) Ftemp[pl]=flux[i][j][k][pl]*surface; // UEVOLVE form
      // GODMARK: for finite volume method, below doesn't change the result if eomfunc=gdet.
      // Otherwise flux would have to be completely recomputed for gdet case JUST for diagnostic to be consistent at higher order
      UtoU(UEVOLVE,UDIAG,&geom,Ftemp,Ftempdiag); // convert to diag form
      PLOOP(pl) frdot[i][pl]+=Ftempdiag[pl];
    }

  }
#endif
  // GODMARK
  // want all fluxes vs theta on horizon


  return(0);

}


#define DEBUGFRLOOP 1

// write the flux vs. radius
void frdotout(void)
{
  int i,j,k,pl,l;
  SFTYPE ftemp;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  SFTYPE frdottemp[N1][NPR];
  SFTYPE *frtot;
  int ospos1;
  FILE*frout;
  static int firsttime=1;

  if(numprocs==1){
    frtot=(SFTYPE (*))(&frdot[0][0]);
  }
  else{
#if(USEMPI)
    if(myid==0){
      frtot=(SFTYPE*) malloc(sizeof(SFTYPE)*totalsize[1]*NPR);
      if(frtot==NULL){
	dualfprintf(fail_file,"Cannot get frtot memory\n");
	myexit(1);
      }
      else{
	for(i=0;i<totalsize[1];i++) PLOOP(pl){
	  frtot[i*NPR+pl]=0;
	}
      }
      for(l=0;l<numprocs;l++){ // just go over all cpus and assume only added to frdot per cpu for correct cpus.
	ospos1=(l%ncpux1)*N1;
	if(l==0){ // assumes cpu=0 is main cpu and is on horizon
	  for(i=0;i<N1;i++) PLOOP(pl){
	    frdottemp[i][pl]=frdot[i][pl];
	  }
	}
	else{
	  MPI_Irecv(frdottemp,N1*NPR,MPI_SFTYPE,l,l,MPI_COMM_WORLD,&rrequest);
	  MPI_Wait(&rrequest,&mpichstatus);
	}
	for(i=0;i<N1;i++) PLOOP(pl){
	  frtot[(ospos1+i)*NPR+pl]+=frdottemp[i][pl];
#if(DEBUGFRLOOP)
	  if((ospos1+i)*NPR+pl>=totalsize[1]*NPR){
	    dualfprintf(fail_file,"outside bounds: %d\n",(ospos1+i)*NPR+pl);
	    myexit(1);
	  }
#endif
	}
      }
    }
    else{
      MPI_Isend(frdot,N1*NPR,MPI_SFTYPE,0,myid,MPI_COMM_WORLD,&srequest);
      MPI_Wait(&srequest,&mpichstatus);
    }
#endif
  }
  // now we have frtot with full fluxes vs. radius (totalsize[1]), so output

  if(myid==0){
    frout=fopen("frdot.out","at");
    if(frout==NULL){
      dualfprintf(fail_file,"Cannot open frdot.out\n");
      myexit(1);
    }
    if(firsttime){
      fprintf(frout,"%21.15g %ld %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
	      t,realnstep,totalsize[1],totalsize[2],totalsize[3],startx[1],startx[2],startx[3],dx[1],dx[2],dx[3]);
      fflush(frout);
    }

    for(i=0;i<totalsize[1];i++){
      fprintf(frout,"%21.15g %d ",t,i);
      PDUMPLOOP(pl){// dump only dump prims
	fprintf(frout,"%21.15g ",frtot[i*NPR+pl]);
      }
      fprintf(frout,"\n");
    }
    
    if(frout!=NULL) fclose(frout);
#if(USEMPI)
		if( numprocs != 1 ) free(frtot);  //atch corrected: without the "if" on 1 proc. in MPI this leads to freeing the stack memory
#endif
  }
  firsttime=0;
}





void init_varstavg(void)
{
  int i,j,k,ii;

  ZLOOP{
    for(ii=0;ii<NUMNORMDUMP;ii++){
      normalvarstavg[i][j][k][ii]=0.0;
      anormalvarstavg[i][j][k][ii]=0.0;
    }      
    for(ii=0;ii<NDIM;ii++){
#if(CALCFARADAYANDCURRENTS)
      jcontavg[i][j][k][ii]=0.0;
      jcovtavg[i][j][k][ii]=0.0;
      ajcontavg[i][j][k][ii]=0.0;
      ajcovtavg[i][j][k][ii]=0.0;
#endif
      massfluxtavg[i][j][k][ii]=0.0;
      amassfluxtavg[i][j][k][ii]=0.0;
    }
    for(ii=0;ii<NUMOTHER;ii++){
      othertavg[i][j][k][ii]=0.0;      
      aothertavg[i][j][k][ii]=0.0;      
    }
#if(CALCFARADAYANDCURRENTS)
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][k][ii]=0.0;
      fcovtavg[i][j][k][ii]=0.0;
      afcontavg[i][j][k][ii]=0.0;
      afcovtavg[i][j][k][ii]=0.0;
    }
#endif
    for(ii=0;ii<NUMSTRESSTERMS;ii++){
      tudtavg[i][j][k][ii]=0.0;
      atudtavg[i][j][k][ii]=0.0;
    }      
    
  }
}

void final_varstavg(FTYPE IDT)
{
  int i,j,k,ii;

  ZLOOP{
    for(ii=0;ii<NUMNORMDUMP;ii++){
      normalvarstavg[i][j][k][ii]=normalvarstavg[i][j][k][ii]*IDT;
      anormalvarstavg[i][j][k][ii]=anormalvarstavg[i][j][k][ii]*IDT;
    }      
    for(ii=0;ii<NDIM;ii++){
#if(CALCFARADAYANDCURRENTS)
      jcontavg[i][j][k][ii]=jcontavg[i][j][k][ii]*IDT;
      jcovtavg[i][j][k][ii]=jcovtavg[i][j][k][ii]*IDT;
      ajcontavg[i][j][k][ii]=ajcontavg[i][j][k][ii]*IDT;
      ajcovtavg[i][j][k][ii]=ajcovtavg[i][j][k][ii]*IDT;
#endif
      massfluxtavg[i][j][k][ii]*=IDT;
      amassfluxtavg[i][j][k][ii]*=IDT;
    }
    for(ii=0;ii<NUMOTHER;ii++){
      othertavg[i][j][k][ii]*=IDT;
      aothertavg[i][j][k][ii]*=IDT;
    }
#if(CALCFARADAYANDCURRENTS)
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][k][ii]=fcontavg[i][j][k][ii]*IDT;
      fcovtavg[i][j][k][ii]=fcovtavg[i][j][k][ii]*IDT;
      afcontavg[i][j][k][ii]=afcontavg[i][j][k][ii]*IDT;
      afcovtavg[i][j][k][ii]=afcovtavg[i][j][k][ii]*IDT;
    }
#endif
    for(ii=0;ii<NUMSTRESSTERMS;ii++){
      tudtavg[i][j][k][ii]=tudtavg[i][j][k][ii]*IDT;
      atudtavg[i][j][k][ii]=atudtavg[i][j][k][ii]*IDT;
    }      
  }
}


int set_varstavg(FTYPE tfrac)
{
  int i,j,k;
  int iii;
  int ll;
  int l,ii,aii;
  FTYPE ftemp;
  FTYPE ftemp0,ftemp1,ftemp2,ftemp3,ftemp4,ftemp5,ftemp6;
  FTYPE pgas,bsq;
  FTYPE jcov[NDIM];
  FTYPE fcov[NUMFARADAY];
  FTYPE V[NDIM], vmin[NDIM], vmax[NDIM];
  int ignorecourant;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];


  ZLOOP{

    // just like dumps
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    // if failed, then data output for below invalid, but columns still must exist    
    get_geometry(i, j, k, CENT, &geom);
    if (!failed) {
      if (get_state(p[i][j][k], &geom, &q) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "get_state() dir=0", 1);

      if (vchar(p[i][j][k], &q, 1, &geom, &vmax[1], &vmin[1],&ignorecourant) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "vchar() dir=1or2or3", 1);
      if (vchar(p[i][j][k], &q, 2, &geom, &vmax[2], &vmin[2],&ignorecourant) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "vchar() dir=1or2or3", 2);
      if (vchar(p[i][j][k], &q, 3, &geom, &vmax[3], &vmin[3],&ignorecourant) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "vchar() dir=1or2or3", 3);
    }
    else {// do a per zone check, otherwise set to 0
      whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
      if (get_state(p[i][j][k], &geom, &q) >= 1){
	for (iii = 0; iii < NDIM; iii++)
	  q.ucon[iii]=0;
	for (iii = 0; iii < NDIM; iii++)
	  q.ucov[iii]=0;
	for (iii = 0; iii < NDIM; iii++)
	  q.bcon[iii]=0;
	for (iii = 0; iii < NDIM; iii++)
	  q.bcov[iii]=0;
      }
      if (vchar(p[i][j][k], &q, 1, &geom, &vmax[1], &vmin[1],&ignorecourant) >= 1){
	vmax[1]=vmin[1]=0;
      }
	
      if (vchar(p[i][j][k], &q, 2, &geom, &vmax[2], &vmin[2],&ignorecourant) >= 1){
	vmax[2]=vmin[2]=0;
      }

      if (vchar(p[i][j][k], &q, 2, &geom, &vmax[3], &vmin[3],&ignorecourant) >= 1){
	vmax[3]=vmin[3]=0;
      }


      whocalleducon=0; // return to normal state

    }

    SETFDIVB(divb, p, i, j, k);


    ii=0;
    aii=0;

    for(iii=0;iii<NPR;iii++){
      normalvarstavg[i][j][k][ii++]+=p[i][j][k][iii]*tfrac;
      anormalvarstavg[i][j][k][aii++]+=fabs(p[i][j][k][iii])*tfrac;
    }
    normalvarstavg[i][j][k][ii++]+=divb*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(divb)*tfrac;

    for (iii = 0; iii < NDIM; iii++) normalvarstavg[i][j][k][ii++]+=q.ucon[iii]*tfrac;
    for (iii = 0; iii < NDIM; iii++) anormalvarstavg[i][j][k][aii++]+=fabs(q.ucon[iii])*tfrac;

    for (iii = 0; iii < NDIM; iii++) normalvarstavg[i][j][k][ii++]+=q.ucov[iii]*tfrac;
    for (iii = 0; iii < NDIM; iii++) anormalvarstavg[i][j][k][aii++]+=fabs(q.ucov[iii])*tfrac;

    for (iii = 0; iii < NDIM; iii++) normalvarstavg[i][j][k][ii++]+=q.bcon[iii]*tfrac;
    for (iii = 0; iii < NDIM; iii++) anormalvarstavg[i][j][k][aii++]+=fabs(q.bcon[iii])*tfrac;

    for (iii = 0; iii < NDIM; iii++) normalvarstavg[i][j][k][ii++]+=q.bcov[iii]*tfrac;
    for (iii = 0; iii < NDIM; iii++) anormalvarstavg[i][j][k][aii++]+=fabs(q.bcov[iii])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmin[1]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmin[1])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmax[1]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmax[1])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmin[2]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmin[2])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmax[2]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmax[2])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmin[3]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmin[3])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmax[3]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmax[3])*tfrac;


#if(CALCFARADAYANDCURRENTS)
    lower_vec(jcon[i][j][k],&geom,jcov);
    for(ii=0;ii<NDIM;ii++){
      jcontavg[i][j][k][ii]+=jcon[i][j][k][ii]*tfrac;
      jcovtavg[i][j][k][ii]+=jcov[ii]*tfrac;
      ajcontavg[i][j][k][ii]+=fabs(jcon[i][j][k][ii])*tfrac;
      ajcovtavg[i][j][k][ii]+=fabs(jcov[ii])*tfrac;
    }
#endif
    
    for(ii=0;ii<NDIM;ii++){
      ftemp=(geom.g)*p[i][j][k][RHO]*(q.ucon[ii]);
      massfluxtavg[i][j][k][ii]+=ftemp*tfrac;
      amassfluxtavg[i][j][k][ii]+=fabs(ftemp)*tfrac;
    }

    ii=0;
    aii=0;
    ftemp=(q.ucon[3])/(q.ucon[0]);
    othertavg[i][j][k][ii++]=ftemp*tfrac;
    aothertavg[i][j][k][aii++]=fabs(ftemp)*tfrac;

#if(CALCFARADAYANDCURRENTS)
    lowerf(fcon[i][j][k],&geom,fcov);
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][k][ii]+=fcon[i][j][k][ii]*tfrac;
      fcovtavg[i][j][k][ii]+=fcov[ii]*tfrac;
      afcontavg[i][j][k][ii]+=fabs(fcon[i][j][k][ii])*tfrac;
      afcovtavg[i][j][k][ii]+=fabs(fcov[ii])*tfrac;
    }
#endif

    pgas = pressure_rho0_u(p[i][j][k][RHO],p[i][j][k][UU]);
    bsq=0; for(iii=0;iii<NDIM;iii++) bsq+=(q.bcon[iii])*(q.bcov[iii]);

    // part0
    ii=0;
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp0=pgas*(q.ucon[iii])*(q.ucov[l]);
      tudtavg[i][j][k][ii]+=ftemp0*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp0)*tfrac;
      ii++;
    }
    // part1
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp1=p[i][j][k][RHO]*(q.ucon[iii])*(q.ucov[l]);
      tudtavg[i][j][k][ii]+=ftemp1*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp1)*tfrac;
      ii++;
    }
    // part2
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp2=p[i][j][k][UU]*(q.ucon[iii])*(q.ucov[l]);
      tudtavg[i][j][k][ii]+=ftemp2*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp2)*tfrac;
      ii++;
    }
    // part3
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp3=bsq*(q.ucon[iii])*(q.ucov[l]);
      tudtavg[i][j][k][ii]+=ftemp3*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp3)*tfrac;
      ii++;
    }
    // part4
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp4=pgas*delta(iii,l);
      tudtavg[i][j][k][ii]+=ftemp4*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp4)*tfrac;
      ii++;
    }
    // part5
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp5=0.5*bsq*delta(iii,l);
      tudtavg[i][j][k][ii]+=ftemp5*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp5)*tfrac;
      ii++;
    }
    // part6
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp6=-(q.bcon[iii])*(q.bcov[l]);
      tudtavg[i][j][k][ii]+=ftemp6*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp6)*tfrac;
      ii++;
    }

  }

  return(0);

}


// if doavg==1, then assume this is call before dumping
int average_calc(int doavg)
{
  static FTYPE lastdt;
  static int calls=0;
  static FTYPE tavgi,tavgf;
  static int tavgflag=1;
 
  if(calls>0){ // since need 2 times

    if(tavgflag){
      // gets reached on next call after dump call or first time
      init_varstavg();
      tavgflag=0;
      tavgi=t;
    }

    // always do
    if(set_varstavg(0.5*(lastdt+dt))>=1) return(1);

    if(doavg==1){
      tavgflag=1;
      tavgf=t;
      final_varstavg(1.0/(tavgf-tavgi));
      // expect to dump after this function ends and before next call to this function
    }
  }
  calls++;
  lastdt=dt;

  return(0);
}


void diag_source_all(struct of_geom *ptrgeom, FTYPE *dU,SFTYPE Dt)
{
  int pl,enerregion;
  FTYPE ftemp[NPR];
  FTYPE ftempdiag[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);


  // does not matter what stage
  if(Dt>0.0){

    //    dualfprintf(fail_file,"got here: i=%d j=%d t=%21.15g\n",ptrgeom->i,ptrgeom->j,t);


    ENERREGIONLOOP(enerregion){
      enerpos=enerposreg[enerregion];
      sourceaddterms=sourceaddtermsreg[enerregion];
      sourceadd=sourceaddreg[enerregion];

      if(WITHINENERREGION(ptrgeom->i,ptrgeom->j,ptrgeom->k)){
	PLOOP(pl) ftemp[pl]=Dt*dVF*dU[pl]; // in UEVOLVE form
	// GODMARK: for finite volume method, below doesn't change the result if eomfunc=gdet.
	// Otherwise source would have to be completely recomputed for gdet case JUST for diagnostic to be consistent at higher order
	UtoU(UEVOLVE,UDIAG,ptrgeom,ftemp,ftempdiag); // convert to diag form

	// now assign diagnostic form of source
	PLOOP(pl){
	  sourceadd[pl]+=ftempdiag[pl];
#if(DOLUMVSR)
	  // GODMARK: only correct for diagonal coordinate Jacobian in which each i is same radius for all j
	  if(pl==UU) if(enerregion==0) lumvsr[startpos[1]+ptrgeom->i]+=ftempdiag[pl];
#endif
	} // end PLOOP on diag
      }
    }
  }

}



void diag_source_comp(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt)
{
  int sc,pl,enerregion;
  FTYPE ftemp[NPR];
  FTYPE ftempdiag[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  // does not matter what stage
  if(Dt>0.0){

    ENERREGIONLOOP(enerregion){
      enerpos=enerposreg[enerregion];
      sourceaddterms=sourceaddtermsreg[enerregion];
      sourceadd=sourceaddreg[enerregion];

      if(WITHINENERREGION(ptrgeom->i,ptrgeom->j,ptrgeom->k)){
	SCLOOP(sc){
	  PLOOP(pl) ftemp[pl]=Dt*dVF*dUcomp[sc][pl]; // in UEVOLVE form
	  // GODMARK: for finite volume method, below doesn't change the result if eomfunc=gdet.
	  // Otherwise source would have to be completely recomputed for gdet case JUST for diagnostic to be consistent at higher order
	  UtoU(UEVOLVE,UDIAG,ptrgeom,ftemp,ftempdiag); // convert to diag form

	  // now assign diagnostic form of source
	  PLOOP(pl){
	    sourceaddterms[sc][pl]+=ftempdiag[pl];
	  } // end PLOOP on diag
	} // end SCLOOP
      }
    }
  }

}



// compute dissipated energy due to (e.g.) shocks and reconnection
int diss_compute(int evolvetype, int inputtype, FTYPE *U, struct of_geom *ptrgeom, FTYPE *pr)
{
  extern int Utoprimdiss(int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, int *otherfail);
  FTYPE prother[NPR];
  int otherfail;
  int enerregion;
  FTYPE dissenergy;
  int pl;
  struct of_state q;
  FTYPE Unew[NPR];
  FTYPE primtoUcons;

  if(DOENTROPY==DOEVOLVECOMPAREENTROPY){

    PLOOP(pl) prother[pl]=pr[pl]; // guess

    // invert with entropy evolution
    Utoprimdiss(evolvetype, inputtype, U,  ptrgeom, prother,&otherfail);

    // result now contains internal energy (prother[UU,ENTROPY]) as derived from entropy evolution
    // notice that all other quantities are could also be different (if doentropy==evolvefullentropy), hence the prother variable for temporary storage.
    
    // at this point, ie version of entropy
    // just overwrite entropy primitive, leave rest same as from full energy equation
    pr[ENTROPY]=prother[UU];
    // now compare pr[UU] and pr[ENTROPY] with some kind of diagnostic?
    if(evolvetype==EVOLVEUTOPRIM){
      // then during evolution and pr[UU]-pr[ENTROPY] is relevant to physical calculation
      // store difference
      if(otherfail==UTOPRIMNOFAIL){ // only use if inversion succeeded (otherwise assume entropy evolution wanted negative internal energy and so not a good solution)
	// only correct for diagonal coordinate Jacobian in which each i is same radius for all j                         

	dissenergy=pr[UU]-pr[ENTROPY];

	// only for enerregion==0
	if(DODISSVSR) dissvsr[startpos[1]+ptrgeom->i]+=dissenergy*ptrgeom->g * dVF;

	if(DODISS){
	  // local integral
	  ENERREGIONLOOP(enerregion){
	    diss=dissreg[enerregion];
	    if( WITHINENERREGION(ptrgeom->i,ptrgeom->j,ptrgeom->k) ){
	      diss[0]+=dissenergy*ptrgeom->g * dVF; // actual energy
	    }
	  }

	  // function over all space to be written as dump file
	  // energy density, which can be integrated in SM since grid is given
	  dissfunpos[ptrgeom->i][ptrgeom->j][ptrgeom->k][0]+=dissenergy;
	}
	//	    dualfprintf(fail_file,"diss=%g\n",diss[ptrgeom->i][ptrgeom->j][ptrgeom->k]);

      }// end if didn't fail

      // report failure to invert
      if(DODISS){
	if(otherfail!=UTOPRIMNOFAIL) dissfunpos[ptrgeom->i][ptrgeom->j][ptrgeom->k][1]+=1.0;
      }
    }// end if evolving
  }// end if doing comparison


  // now must redefine U[ENTROPY] so consistent with p(U[normalgrmhd])
  if(DOENOFLUX==ENOFINITEVOLUME){
    //    primtoUcons=UENTROPY; // not UENTROPY since this doesn't have gdet factor!
    if (get_state(pr, ptrgeom, &q) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
    if (primtoU(UEVOLVE, pr, &q, ptrgeom, Unew) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 1);
    U[ENTROPY] = Unew[ENTROPY]; // now conserved entropy is consistent with real primitive state
  }
  // this is done automatically if doing NOENOFLUX since U is obtained again from new primitives.


  return(0);

}
