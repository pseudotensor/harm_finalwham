
/**
 *
 * this contains the generic piece of code for advancing
 * the primitive variables 
 *
**/


#include "decs.h"


/** end algorithmic choices **/


int step_ch(void)
{
  int step_ch_simplempi(void);
  int step_ch_supermpi(void);


#if(SIMULBCCALC==-1)
  MYFUN(step_ch_simplempi(),"step_ch.c:step_ch()", "step_ch_simplempi()", 1);
#else
  MYFUN(step_ch_supermpi(),"step_ch.c:step_ch()", "step_ch_supermpi()", 1);
#endif

  /* done! */
  return (0);
}




int step_ch_simplempi()
{
  int boundstage;
  SFTYPE mydt;
  int stage, stagei,stagef;
  FTYPE ndt,lastndt;
  FTYPE (*pi)[N2M][N3M][NPR];
  FTYPE (*pb)[N2M][N3M][NPR];
  FTYPE (*pf)[N2M][N3M][NPR];
  FTYPE (*prevpf)[N2M][N3M][NPR];
  FTYPE (*pii[4])[N2M][N3M][NPR];
  FTYPE (*pbb[4])[N2M][N3M][NPR];
  FTYPE (*pff[4])[N2M][N3M][NPR];
  FTYPE (*uii[4])[N2M][N3M][NPR];
  FTYPE (*uff[4])[N2M][N3M][NPR];
  FTYPE (*ucum[4])[N2M][N3M][NPR];
  //  FTYPE alphaik[MAXSTAGES][MAXSTAGES],betaik[MAXSTAGES];
  FTYPE Cunew[MAXTIMEORDER][4],CUf[MAXTIMEORDER][4];
  int i, j, k;
  int numtimeorders;
  int timeorder;
  SFTYPE dt4diag;
  int finalstep;
  extern int diag_flux(FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],SFTYPE Dt);
  void setup_rktimestep(int *numtimeorders,
			FTYPE (*p)[N2M][N3M][NPR],FTYPE (*pk)[N1M][N2M][N3M][NPR],
			FTYPE (*pii[4])[N2M][N3M][NPR],FTYPE (*pbb[4])[N2M][N3M][NPR],FTYPE (*pff[4])[N2M][N3M][NPR],
			FTYPE (*uii[4])[N2M][N3M][NPR],FTYPE (*uff[4])[N2M][N3M][NPR],FTYPE (*ucum[4])[N2M][N3M][NPR],
			FTYPE (*Cunew)[4],FTYPE (*CUf)[4]);
  extern int advance(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
		     FTYPE ui[][N2M][N3M][NPR],FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
		     FTYPE *CUf,FTYPE *Cunew,int timeorder, int numtimeorders, FTYPE *ndt);

  int asym_compute_2(FTYPE (*prim)[N2M][N3M][NPR]);



#if(CALCFARADAYANDCURRENTS)

  if((WHICHCURRENTCALC==CURRENTCALC0)||(WHICHCURRENTCALC==CURRENTCALC2)){
    // for this method, all J components are located back in time
    current_doprecalc(CURTYPEX,p);
    current_doprecalc(CURTYPEY,p);
    current_doprecalc(CURTYPEZ,p);
    // get faraday in case used for time averaging or dumping it
    current_doprecalc(CURTYPEFARADAY,p);
  }

#endif

  // setup time-stepping
  setup_rktimestep(&numtimeorders,p,pk,pii,pbb,pff,uii,uff,ucum,Cunew,CUf);

  // global debug tracking var
  numstepparts=numtimeorders;



  // general temporal loop
  for(timeorder=0;timeorder<numtimeorders;timeorder++){
  // global debug tracking var
    steppart=timeorder;
    //    if(timeorder==numtimeorders-1) dt4diag=dt; else dt4diag=-1.0; // set dt for accounting diagnostics
    if(timeorder==numtimeorders-1) finalstep=1; else finalstep=0;

#if(PRODUCTION==0)
    trifprintf("|%ds",timeorder);
#endif


    //the starting and the ending times of the current substep
    if(timeorder!=numtimeorders-1){
      tstepparti = t + CUf[timeorder][3] * dt;
      tsteppartf = t + CUf[timeorder][2] * dt +  CUf[timeorder][3] * dt;
    }
    else{
      tstepparti=t;
      tsteppartf=t+dt;
    }


    // prefixup
#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD )
    // force-free and cold GRMHD don't use pre_fixup, currently, even if could
    MYFUN(pre_fixup(-1, pii[timeorder]),"step_ch.c:step_ch_simple()", "pre_fixup()", 1);
#endif

    // advance
    MYFUN(advance(-1, pii[timeorder], pbb[timeorder], pff[timeorder], uii[timeorder], uff[timeorder], ucum[timeorder], CUf[timeorder], Cunew[timeorder],timeorder,numtimeorders,&ndt),"step_ch.c:step_ch_simplempi()", "advance()", 1);


    // force-free and cold GRMHD don't use post_fixup
#if( (EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD) && (1 == CHECKSOLUTION || UTOPRIMADJUST == UTOPRIMAVG) )
    // if CHECKSOLUTION==1, then need values to be bounded right now, since use them to check whether even good solutions are really good.
    // post_fixup() will use previous time step pff boundary values to fixup_utoprim() if this is not called.
    // bound advanced values before post_fixup() so fixup_utoprim() has updated boundary values to base fixup on.

    MYFUN(bound_prim(-1,pff[timeorder]),"step_ch.c:step_ch_simplempi()", "bound_prim()", 1);


#if(PRODUCTION==0)
    trifprintf("b");
#endif

    // done when all timeorders are completed, so stencil used doesn't matter

    // postfixup
    // post_fixup: If this modifies values and then uses the modified values for other modified values, then must bound_prim() after this
    // if one doesn't care about MPI being same as non-MPI, then can move bound_prim() above to below error_check() and then remove prc<-pv in fixup_utoprim()
    MYFUN(post_fixup(-1,pff[timeorder],pbb[timeorder],finalstep),"step_ch.c:advance()", "post_fixup()", 1);

#if(PRODUCTION==0)
    trifprintf( "x");
#endif


#endif
    // must check before MPI operation (since asymmetries would
    // desynchronize cpus
    MYFUN(error_check(1),"step_ch.c", "error_check", 1);

    // bound final values (comes after post_fixup() since changes made by post_fixup)
    //#if(MPIEQUALNONMPI)

#if(ASYMDIAGCHECK)
    dualfprintf(fail_file,"1before bound\n");
    asym_compute_2(pff[timeorder]);
    dualfprintf(fail_file,"2before bound\n");
#endif


    MYFUN(bound_prim(-1,pff[timeorder]),"step_ch.c:step_ch_simplempi()", "bound_prim()", 2);

    //#endif


#if(CALCFARADAYANDCURRENTS)
    if((WHICHCURRENTCALC==CURRENTCALC0)||(WHICHCURRENTCALC==CURRENTCALC2)){
      // puts J at the time center, but hard to know if RK is at mid point in time except for midpoint method
      // compute current_doprecalc if near half-step in time
      if(
	 ((numtimeorders>=3)&&(timeorder==1))
	 ||((numtimeorders<=2)&&(timeorder==0))
	 )
	current_doprecalc(CURTYPET,pff[timeorder]); // should be called using half-time step data
    }
#endif

    /* perform diagnostics */
    // no error check since assume if step_ch passed, diag(1) will pass
    if (DODIAGS && DODIAGEVERYSUBSTEP ){ //SASMARK -- moved the diags calls here
      pdump = pff[timeorder];
      diag(1);
#if(PRODUCTION==0)
      trifprintf( "D");
#endif
    }

  }// end timeorder
  
  /////////////////
  //
  // FOR ANY TIME ORDER
  //
  //////////////////

#if(CALCFARADAYANDCURRENTS)
  // compute final faradays
  if(WHICHCURRENTCALC==CURRENTCALC1){
    // compute faraday components needed for time centering of J
    current_doprecalc(CURTYPET,p);
    // J is located at this time
    current_doprecalc(CURTYPEX,p);
    current_doprecalc(CURTYPEY,p);
    current_doprecalc(CURTYPEZ,p);
    current_doprecalc(CURTYPEFARADAY,p); // for diagnostics
  }
  // compute current
  current_calc(cfaraday);
#endif

#if(ACCURATEDIAG==0)
  // compute flux diagnostics (uses global F1/F2)
  // this doesn't exactly make conservation work -- should have in middle step point using full step.  But if PARA, no middle point that's exact.
  // think about where to put this
  // GODMARK
  diag_flux(p,F1, F2, F3, dt); // should use REAL dt, not within a timeorderd RK integration step
#endif

#if(PRODUCTION==0)
  trifprintf( "d");
#endif

  /* check timestep
     if (dt < 1.e-9) {
     trifprintf( "timestep too small\n");
     myexit(11);
     }
  */

  // increment time
  t += dt;
  tstepparti = tsteppartf = t;

  realnstep++;

  // set next timestep
  // find global minimum value of ndt over all cpus
  mpifmin(&ndt);
  if (ndt > SAFE * dt)    ndt = SAFE * dt;
  dt = ndt;
  // don't step beyond end of run
  if (t + dt > tf) dt = tf - t;


  /* done! */
  return (0);
}




int step_ch_supermpi()
{
  int boundstage;
  SFTYPE mydt;
  int stage, stagei,stagef;
  int timeorder;
  FTYPE ndt,lastndt;
  FTYPE (*pi)[N2M][N3M][NPR];
  FTYPE (*pb)[N2M][N3M][NPR];
  FTYPE (*pf)[N2M][N3M][NPR];
  FTYPE (*prevpf)[N2M][N3M][NPR];
  FTYPE (*pii[4])[N2M][N3M][NPR];
  FTYPE (*pbb[4])[N2M][N3M][NPR];
  FTYPE (*pff[4])[N2M][N3M][NPR];
  FTYPE (*uii[4])[N2M][N3M][NPR];
  FTYPE (*uff[4])[N2M][N3M][NPR];
  FTYPE (*ucum[4])[N2M][N3M][NPR];
  //  FTYPE alphaik[MAXSTAGES][MAXSTAGES],betaik[MAXSTAGES];
  FTYPE Cunew[MAXTIMEORDER][4],CUf[MAXTIMEORDER][4];
  int i, j, k, pl;
  int numtimeorders;
  SFTYPE dt4diag;
  int finalstep;
  extern int diag_flux(FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR],SFTYPE Dt);
  void setup_rktimestep(int *numtimeorders,
			FTYPE (*p)[N2M][N3M][NPR],FTYPE (*pk)[N1M][N2M][N3M][NPR],
			FTYPE (*pii[4])[N2M][N3M][NPR],FTYPE (*pbb[4])[N2M][N3M][NPR],FTYPE (*pff[4])[N2M][N3M][NPR],
			FTYPE (*uii[4])[N2M][N3M][NPR],FTYPE (*uff[4])[N2M][N3M][NPR],FTYPE (*ucum[4])[N2M][N3M][NPR],
			FTYPE (*Cunew)[4],FTYPE (*CUf)[4]);
  extern int advance(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
		     FTYPE ui[][N2M][N3M][NPR],FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
		     FTYPE *CUf,FTYPE *Cunew,int timeorder, int numtimeorders, FTYPE *ndt);



  // setup time-stepping
  setup_rktimestep(&numtimeorders,p,pk,pii,pbb,pff,uii,uff,ucum,Cunew,CUf);


  // SPECIAL BOUNDARY/COMPUTATION MPI METHOD (out of date, and doesn't yet work right even if essentially complete code)
  /* check timestep */
  if (dt < MINDT) {
    trifprintf( "timestep too small\n");
    myexit(11);
  }
  
  lastndt=1E30; // initialize lastndt
  for(timeorder=1;timeorder<=numtimeorders;timeorder++){
#if(PRODUCTION==0)
    trifprintf("-to%d/%d-",timeorder,numtimeorders);
#endif
    if(numtimeorders==2){
      // note that pb (used for flux calc which has a stencil
      // calculation) must be different from pf so new stencils in
      // different stages won't affect stencil calculations -- must
      // use old values, not new from most previous temporary stage
      //
      // pi however can be the same as pf since each pi is replaced 1
      // zone at a time with a 0 stencil.
      if(timeorder==1){
	pi=p;
	pb=p;
	pf=pk[0]; // different already, so good for simulbccalc
	prevpf=p; // previous final true array
	mydt=0.5*dt;
      }
      else if(timeorder==2){
	pi=p;
	pb=pk[0];
	pf=p;
	prevpf=pk[0];
	mydt=dt;
      }
    }
    else if(numtimeorders==1){
      pi=p;
      pb=p;
      if(SIMULBCCALC<=0) pf=p; else pf=pk[0]; // need to be different if doing simulbccalc
      prevpf=p;
      mydt=dt;
    }
    if(SIMULBCCALC<=0){ stagei=STAGEM1; stagef=STAGEM1; }
    else if(SIMULBCCALC==1) { stagei=STAGE0; stagef=STAGE2;}
    else if(SIMULBCCALC==2) { stagei=STAGE0; stagef=STAGE5;}

    
    // initialize bound stage
    if(SIMULBCCALC) boundstage=STAGE0;
    else boundstage=STAGEM1;
    for(stage=stagei;stage<=stagef;stage++){
#if(PRODUCTION==0)
      if(SIMULBCCALC) trifprintf("!s%d!",stage);
#endif
      // setup stage loop
#if(SIMULBCCALC==2)
#if(TYPE2==1)
      // GODMARK: isf1, etc. are NOT defined?!
	STAGECONDITION(0,N1-1,0,N2-1,isc,iec,jsc,jec);
	STAGECONDITION(0,N1,-1,N2,isf1,ief1,jsf1,jef1);
	STAGECONDITION(-1,N1,0,N2,isf2,ief2,jsf2,jef2);
	STAGECONDITION(0,N1,0,N2,ise,iee,jse,jee);
	STAGECONDITION(0,N1,0,N2-1,isf1ct,ief1ct,jsf1ct,jef1ct);
	STAGECONDITION(0,N1-1,0,N2,isf2ct,ief2ct,jsf2ct,jef2ct);
	STAGECONDITION(-1,N1,-1,N2,isdq,iedq,jsdq,jedq);
	STAGECONDITION(-2,N1+1,-2,N2+1,ispdq,iepdq,jspdq,jepdq);
	// GODMARK : probably not right for general boundary condition size
#endif
#endif

      // only bounding if safe zones, unsafe needs bz complete
      if(stage<STAGE2){
	bound_prim(boundstage, prevpf);
	if(stage!=STAGEM1) boundstage++;
      }

      // done here instead of local since pseudo-complicated calculation that might slow the dq calculation if done locally per zone
      MYFUN(pre_fixup(stage, prevpf),"step_ch.c:advance()", "pre_fixup()", 1);

      // go from previous solution to new solution
      partialstep=timeorder;      
      // not right for numtimeorders==4 // GODMARK
      // advance
      MYFUN(advance(-1, pii[timeorder], pbb[timeorder], pff[timeorder], uii[timeorder], uff[timeorder], ucum[timeorder], CUf[timeorder], Cunew[timeorder],timeorder,numtimeorders,&ndt),"step_ch.c:step_ch_supermpi()", "advance()", 1);
      //      MYFUN(advance(-1,pii[timeorder], ulast, pbb[timeorder], CUf[timeorder], pff[timeorder], Cunew[timeorder],unew,timeorder,numtimeorders,&ndt),"step_ch.c:step_ch()", "advance()", 1);
      //      MYFUN(advance(stage, pi, pb, mydt, pf, 0.0, unew, stage, numtimeorders, &ndt),"step_ch.c:step_ch()", "advance()", 1);

      // must check before MPI operation (since asymmetries would desynchronize cpus)
      if(stage<STAGE2){
	MYFUN(error_check(1),"step_ch.c", "error_check", 1);
      }
      if(stage!=STAGEM1){
	if(stage<STAGE2){
	  bound_prim(boundstage, prevpf);
	  boundstage++;
	}
      }
      if(timeorder==numtimeorders){
	if(ndt>lastndt) ndt=lastndt; // don't change if last was lower
	else lastndt=ndt; // new is lower, keep it
      }
    }
    if(timeorder==numtimeorders){// only do on full step
      // find global minimum value of ndt over all cpus
      mpifmin(&ndt);
    }
    // done when all stages are completed, so stencil used doesn't matter
    MYFUN(post_fixup(-1,pf,pb,1),"step_ch.c:advance()", "post_fixup()", 1);
  }


#if(ACCURATEDIAG==0)
  /* evaluate diagnostics based on fluxes on second pass (Dt=dt)*/
  // need to do this every timestep, but only after all stages are complete and on full timestep
  diag_flux(p, F1, F2, F3, dt); // should use REAL dt, not within a timeorderd RK integration step
#endif


  // copy the contents to the final working array
  if((numtimeorders==1)&&(SIMULBCCALC)) FULLLOOP PLOOP(pl) p[i][j][k][pl]=pf[i][j][k][pl];
  
  
  /* increment time */
  t += dt;
  realnstep++;
  
  // new timestep
  if (ndt > SAFE * dt)    ndt = SAFE * dt;
  dt = ndt;    
  /* but don't step beyond end of run */
  if (t + dt > tf)    dt = tf - t;

  

  /* done! */
  return (0);
}




// for the ith stage:

// Uf^i = ulast^i = CUf^{i0} Ui^i + CUf^{i1} ulast^i + CUf^{i2} dU^i

// unew^i = Cunew^{i0} Ui^i + Cunew^{i1} dU^i + Cunew^{i2} Uf^i
void setup_rktimestep(int *numtimeorders,
		      FTYPE (*p)[N2M][N3M][NPR],FTYPE (*pk)[N1M][N2M][N3M][NPR],
		      FTYPE (*pii[4])[N2M][N3M][NPR],FTYPE (*pbb[4])[N2M][N3M][NPR],FTYPE (*pff[4])[N2M][N3M][NPR],
		      FTYPE (*uii[4])[N2M][N3M][NPR],FTYPE (*uff[4])[N2M][N3M][NPR],FTYPE (*ucum[4])[N2M][N3M][NPR],
		      FTYPE (*Cunew)[4],FTYPE (*CUf)[4])
{


  // to avoid special copying of final pff->p, always use p as final pff
  if(TIMEORDER==4){
    // RK4 stepping
    *numtimeorders=TIMEORDER;

    // Ui ulast dU(pb)
    CUf[0][0]=1.0;  CUf[0][1]=0.0;      CUf[0][2]=0.5;  CUf[0][3] = 0.0;
    CUf[1][0]=1.0;  CUf[1][1]=0.0;      CUf[1][2]=0.5;  CUf[1][3] = 0.0;
    CUf[2][0]=1.0;  CUf[2][1]=0.0;      CUf[2][2]=1.0;  CUf[2][3] = 0.0;
    CUf[3][0]=1.0;  CUf[3][1]=0.0;      CUf[3][2]=1.0;  CUf[3][3] = 0.0;

    // Ui dU(Ub) Uf
    Cunew[0][0]=1.0;  Cunew[0][1]=1.0/6.0;      Cunew[0][2]=0.0;
    Cunew[1][0]=0.0;  Cunew[1][1]=1.0/3.0;      Cunew[1][2]=0.0;
    Cunew[2][0]=0.0;  Cunew[2][1]=1.0/3.0;      Cunew[2][2]=0.0;
    Cunew[3][0]=0.0;  Cunew[3][1]=1.0/6.0;      Cunew[3][2]=0.0;

    //primitive values used for initial state, fluxes, final state (where you output)
    pii[0]=p;    pbb[0]=p;       pff[0]=pk[0]; // produces U1
    pii[1]=p;    pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
    pii[2]=p;    pbb[2]=pk[1];   pff[2]=pk[0]; // produces U3
    pii[3]=p;    pbb[3]=pk[0];   pff[3]=p; // produces U4 (only dU part used)
    
    uii[0]=uinitial;  uff[0]=ulast; ucum[0]=unew;
    uii[1]=uinitial;  uff[1]=ulast; ucum[1]=unew;
    uii[2]=uinitial;  uff[2]=ulast; ucum[2]=unew;
    uii[3]=uinitial;  uff[3]=ulast; ucum[3]=unew;

  }
  else if(TIMEORDER==3){
    // TVD optimal RK3 method as in Shu's report
    *numtimeorders=3;
    
    // Ui ulast dU(pb)
    CUf[0][0]=1.0;      CUf[0][1]=0.0;      CUf[0][2]=1.0;      CUf[0][3] = 0.0;
    CUf[1][0]=3.0/4.0;  CUf[1][1]=1.0/4.0;  CUf[1][2]=1.0/4.0;  CUf[1][3] = 0.0;
    CUf[2][0]=1.0/3.0;  CUf[2][1]=2.0/3.0;  CUf[2][2]=2.0/3.0;  CUf[2][3] = 0.0;
    
    // Ui dU(Ub) Uf
    // unew=U3
    Cunew[0][0]=0.0;   Cunew[0][1]=0.0;      Cunew[0][2]=0.0;
    Cunew[1][0]=0.0;   Cunew[1][1]=0.0;      Cunew[1][2]=0.0;
    Cunew[2][0]=0.0;   Cunew[2][1]=0.0;      Cunew[2][2]=1.0;
    
    //always starting the substeps from the initial time
    pii[0]=p;      pbb[0]=p;       pff[0]=pk[0]; // produces U1
    pii[1]=p;      pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
    pii[2]=p;      pbb[2]=pk[1];   pff[2]=p; // produces U3

    uii[0]=uinitial;  uff[0]=ulast; ucum[0]=unew;
    uii[1]=uinitial;  uff[1]=ulast; ucum[1]=unew;
    uii[2]=uinitial;  uff[2]=ulast; ucum[2]=unew;

  }
  else if(TIMEORDER==2){
#if(1)
    // midpoint method

    *numtimeorders=TIMEORDER;

    // old unew not used for this method (i.e. [?][1]=0)
    CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=0.5; CUf[0][3] = 0.0;
    CUf[1][0]=1.0; CUf[1][1]=0.0; CUf[1][2]=1.0; CUf[1][3] = 0.0;

    // Ui dU(Ub) Uf
    // unew=U2
    Cunew[0][0]=0.0;   Cunew[0][1]=0.0;      Cunew[0][2]=0.0;
    Cunew[1][0]=0.0;   Cunew[1][1]=0.0;      Cunew[1][2]=1.0;

    pii[0]=p;    pbb[0]=p;       pff[0]=pk[0];
    pii[1]=p;    pbb[1]=pk[0];   pff[1]=p;

    uii[0]=uinitial;  uff[0]=ulast; ucum[0]=unew;
    uii[1]=uinitial;  uff[1]=ulast; ucum[1]=unew;

#else
    *numtimeorders=TIMEORDER;
    // TVD RK2 (Chi-Wang Shu 1997 - eq 4.10)
    // actually less robust than generic midpoint method above

    CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=1.0; CUf[0][3] = 0.0;
    CUf[1][0]=0.5; CUf[1][1]=0.5; CUf[1][2]=0.5; CUf[1][3] = 0.0;

    // Ui dU(Ub) Uf
    // unew=U2
    Cunew[0][0]=0.0;   Cunew[0][1]=0.0;      Cunew[0][2]=0.0;
    Cunew[1][0]=0.0;   Cunew[1][1]=0.0;      Cunew[1][2]=1.0;

    pii[0]=p;    pbb[0]=p;       pff[0]=pk[0];
    pii[1]=p;    pbb[1]=pk[0];   pff[1]=p;

    uii[0]=uinitial;  uff[0]=ulast; ucum[0]=unew;
    uii[1]=uinitial;  uff[1]=ulast; ucum[1]=unew;

#endif
  }
  else if(TIMEORDER==1){
    // Euler method
    *numtimeorders=TIMEORDER;

    CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=1.0; CUf[0][3] = 0.0;

    // Ui dU(Ub) Uf
    // unew=U1
    Cunew[0][0]=0.0;   Cunew[0][1]=0.0;      Cunew[0][2]=1.0;

    pii[0]=p;    pbb[0]=p;    pff[0]=p;

    uii[0]=uinitial;  uff[0]=ulast; ucum[0]=unew;

  }

}





int bound_prim(int boundstage, FTYPE prim[][N2M][N3M][NPR])
{
  extern int bound_prim_user(int boundstage, FTYPE prim[][N2M][N3M][NPR]);
  extern int bound_prim_user_after_mpi(int boundstage, FTYPE prim[][N2M][N3M][NPR]);

  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){
    
    MYFUN(bound_prim_user(boundstage,prim),"step_ch.c:bound_prim()", "bound_prim_user()", 1);
    
  }// end if stage0 or stagem1
  
  if(USEMPI){

    MYFUN(bound_mpi(boundstage, BOUNDPRIMTYPE, prim, NULL, NULL, NULL),"step_ch.c:bound_prim()", "bound_mpi()", 1);

  }

  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    MYFUN(bound_prim_user_after_mpi(boundstage,prim),"step_ch.c:bound_prim()", "bound_prim_user_after_mpi()", 1);

  }// end if stage0 or stagem1

  
  return(0);
}


int bound_flux(int boundstage, FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR])
{
  extern int bound_flux_user(int boundstage, FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]);


 // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    MYFUN(bound_flux_user(boundstage,F1,F2,F3),"step_ch.c:bound_flux()", "bound_flux_user()", 1);

  }// end if stage0 or stagem1


  if(USEMPI){

    MYFUN(bound_mpi(boundstage, BOUNDFLUXTYPE, NULL, F1, F2, F3),"step_ch.c:bound_flux()", "bound_mpi()", 1);

  }

  return(0);
}


int bound_pflag(int boundstage, int prim[][N2M][N3M][NUMPFLAGS])
{
  extern int bound_pflag_user(int boundstage, int prim[][N2M][N3M][NUMPFLAGS]);

  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    MYFUN(bound_pflag_user(boundstage,prim),"step_ch.c:bound_pflag()", "bound_pflag_user()", 1);

  }// end bound stage


  if(USEMPI){

    MYFUN(bound_mpi_int(boundstage, prim),"step_ch.c:bound_pflag()", "bound_mpi_int()", 1);

  }

  return(0);

}
