
#include "decs.h"


int init(int argc, char *argv[])
{
  extern int pre_init(int argc, char *argv[]);
  extern int init_defgrid(void);
  extern int init_defglobal(void);
  extern int init_defconsts(void);
  extern int init_grid_post_set_grid(void);
  extern int post_par_set(void);
  extern int init_primitives(FTYPE p[][N2M][N3M][NPR]);
  extern void filterffde(int i, int j, int k, FTYPE *pr);
  extern void filter_coldgrmhd(int i, int j, int k, FTYPE *pr);
  extern int set_dt(FTYPE (*prim)[N2M][N3M][NPR], SFTYPE *dt);  
  //
  extern int init_conservatives(FTYPE p[][N2M][N3M][NPR], FTYPE Utemp[][N2M][N3M][NPR], FTYPE U[][N2M][N3M][NPR]);
  extern int post_init(void);
  int i,j,k;
  int pl;



  fprintf(stderr,"Start init\n"); fflush(stderr);

  pre_init(argc,argv);
  // user specific pre_init
  pre_init_specific_init();


  if(RESTARTMODE==1){
    // loads primitives and if doing ENO method then also fills conserved quantities.  Always assume conserved quantities are in the file, regardless of whether used
    if (restart_init(WHICHFILE) >= 1) {
      dualfprintf(fail_file, "main:restart_init: failure\n");
      return(1);
    }
  }
  else{
    // rest is for fresh start, and thus the above restart must include all these parameters since not otherwise set!


    // define parameters
    init_defglobal(); // init default global parameters
    init_defconsts(); // init default physical constants
    init_global(); // request choices for global parameters

    init_defgrid(); // init default grid
    init_grid(); // request choices for grid/coordinate/metric parameters

    set_coord_parms(); // requires correct defcoord at least
    write_coord_parms(); // output coordinate parameters to file

    // get grid
    set_grid();

    // user post_set_grid function
    init_grid_post_set_grid();

    trifprintf("MCOORD=%d\n",MCOORD);
    trifprintf("COMPDIM=%d\n",COMPDIM);
    trifprintf("MAXBND=%d\n",MAXBND);

    // once all parameters are set, now can set dependent items that may be used to set primitives or conservatives
    post_par_set();



    // user function that should fill p with primitives
    MYFUN(init_primitives(p),"initbase.c:init()", "init_primitives()", 0);

    // dump analytic solution
    pdump=panalytic;
    if (dump(9999) >= 1){
      dualfprintf(fail_file,"unable to print dump file\n");
      return (1);
    }
    pdump=p;



    //    LOOPF{
    // p[i][j][k][UU]=0.0;
    // }

    // mandatory filters on user supplied primitives

    /////////////////////////////
    //
    // Filter to get correct degenerate FFDE solution
    //
    /////////////////////////////// 
    
#if(EOMTYPE==EOMFFDE)
    trifprintf("System filtered to FFDE\n");
    // filter to get force-free
    FULLLOOP{
      filterffde(i,j,k,p[i][j][k]);
    }
#endif

#if(EOMTYPE==EOMCOLDGRMHD)
    trifprintf("System filtered to cold GRMHD\n");
    // filter to get cold GRMHD
    FULLLOOP{
      filter_coldgrmhd(i,j,k,p[i][j][k]);
    }
#endif


    /////////////////////////////
    //
    // Fixup and Bound variables since field may have changed
    // Also setup pre_fixup() type quantities
    //
    /////////////////////////////// 
    
    trifprintf("System Fixup and Bound\n");

#if(FIXUPAFTERINIT)
    if(fixup(STAGEM1,p,0)>=1)
      FAILSTATEMENT("initbase.c:init()", "fixup()", 1);
#endif
    
    if (bound_prim(STAGEM1,p) >= 1)
      FAILSTATEMENT("initbase.c:init()", "bound_prim()", 1);
    
    if(pre_fixup(STAGEM1,p)>=1)
      FAILSTATEMENT("initbase.c:init()", "pre_fixup()", 1);



    // set initial dt correctly rather than random choice
    set_dt(p,&dt);
    trifprintf("dt=%21.15g at t=%21.15g at nstep=%ld at realnstep=%ld\n",dt,t,nstep,realnstep);
    

  }

  // after all parameters and primitives are set, then can set these items
  post_init();
  // user post_init function
  post_init_specific_init();


  if(RESTARTMODE==0){
    // then set initial conserved quantities
    if(DOENOFLUXMEMORY){
      init_conservatives(p, ulast, uinitial);
    }
    else{ // then no difference between point and average and no memory for ulast/unitial
      init_conservatives(p, ulast, unew);
    }
  }
  else{
    // unew should be read in, now assign to unitial for finite volume method
    if(DOENOFLUXMEMORY){
      FULLLOOP PLOOP(pl){
	uinitial[i][j][k][pl] =unew[i][j][k][pl];
      }
    }
  }


  trifprintf("FLUXB=%d\n",FLUXB);

  trifprintf("end init.c\n");
  return (0);

}








int pre_init(int argc, char *argv[])
{
  int ii;
  int dir,pl,sc,fl,floor,enerregion;
  int tscale;
  int i,j,k;
  extern void set_arrays(void);
  int checki;


  // things initialized whether restarting or init fresh

  ranc(0);  // power up random number generator in case used without init

#if(CHECKONINVERSION)
  checki=0;
  strcpy(globalinvtext[checki++],"Qdotnp");
  strcpy(globalinvtext[checki++],"Qtsq");
  strcpy(globalinvtext[checki++],"Qtsqorig");
  strcpy(globalinvtext[checki++],"Bsq");
  strcpy(globalinvtext[checki++],"DD");
  strcpy(globalinvtext[checki++],"QdotB");
  strcpy(globalinvtext[checki++],"WWp");
  strcpy(globalinvtext[checki++],"Qtcon1");
  strcpy(globalinvtext[checki++],"Qtcon2");
  strcpy(globalinvtext[checki++],"Qtcon3");
  strcpy(globalinvtext[checki++],"Bcon1");
  strcpy(globalinvtext[checki++],"Bcon2");
  strcpy(globalinvtext[checki++],"Bcon3");
  if(checki>NUMGLOBALINV){
    dualfprintf(fail_file,"Not enough memory for globalinvtext: checki=%d NUMGLOBALINV=%d\n",checki,NUMGLOBALINV);
  }
#endif


  // below 2 now determined at command line.  See init_mpi() in init_mpi.c (myargs and init_MPI).
  //  RESTARTMODE=0;// whether restarting from rdump or not (0=no, 1=yes)
  //WHICHFILE=0; // see diag.c for dump_cnt and image_cnt starts
  // user defined parameter
  restartonfail=0; // whether we are restarting on failure or not and want special diagnostics

  if(WHICHVEL==VEL3){
    jonchecks=1; // whether to include jon's checks to make sure u^t real and some rho/u limits throughout code
    jonchecks=0;
  }
  else jonchecks=0; // not relevant

  // choice// GODMARK: not convenient location, but needed for init_mpi()
  periodicx1=0;
  periodicx2=0;
  periodicx3=0;// GODMARK: periodic in \phi for 3D spherical polar

  if(USEMPI&&USEROMIO){
    binaryoutput=MIXEDOUTPUT; // choice: mixed or binary
    sortedoutput=SORTED; // no choice
  }
  else{
    // choice
    binaryoutput=TEXTOUTPUT;
    sortedoutput=SORTED;
  }


  // init MPI (assumes nothing in set_arrays.c used here)
  init_mpi(argc, argv);

  /////////////////
  //
  // setup files for writing and reading (must come after init_mpi())
  //
  makedirs();


  // init arrays
  set_arrays();


  // set default variable to dump (must come before init() where if failed or other reasons can dump output)
  udump=unew;
  ubound=unew;
  pdump = p;



  // get EOS for a given EOMTYPE
  pickeos_eomtype(EOMTYPE);


  init_dumps();

  // must go here b4 restart if restarting
  ENERREGIONLOOP(enerregion){
    // used for each region, related to global quantities
    // no need to initialize _tot quantities, they are overwritten during MPI sum in diag.c
    fladd=fladdreg[enerregion];
    fladdterms=fladdtermsreg[enerregion];
    U_init=Ureg_init[enerregion];
    pcum=pcumreg[enerregion];
    pdot=pdotreg[enerregion];
    pdotterms=pdottermsreg[enerregion];
    sourceaddterms=sourceaddtermsreg[enerregion];
    sourceadd=sourceaddreg[enerregion];
    diss=dissreg[enerregion];

    PDUMPLOOP(pl){
      fladd[pl] = 0;
      FLOORLOOP(floor) fladdterms[floor][pl]=0;
      U_init[pl] = 0;
      DIRLOOP(dir){
	pcum[dir][pl]=0;
	pdot[dir][pl]=0;
	FLLOOP(fl) pdotterms[dir][fl][pl]=0;
	if(enerregion==0) FLLOOP(fl) pdottermsjet2[dir][fl][pl]=0; // needed for other not-flux cpus!
      }
      sourceadd[pl] = 0;
      SCLOOP(sc) sourceaddterms[sc][pl] = 0;
    }
    diss[0] = 0;

    if(DOLUMVSR) if(enerregion==0) for(ii=0;ii<ncpux1*N1;ii++) lumvsr[ii]=0;
    if(DODISSVSR) if(enerregion==0) for(ii=0;ii<ncpux1*N1;ii++) dissvsr[ii]=0;
  }
  
  // start counter
  // full loop since user may choose to count something in boundary zones
  if(DODEBUG) FULLLOOP TSCALELOOP(tscale) FLOORLOOP(floor) failfloorcount[i][j][k][tscale][floor]=0;
#if(CALCFARADAYANDCURRENTS)
  // zero out jcon since outer boundaries not set ever since j^\mu involves spatial derivatives that don't exist outside a certain point
  for(pl=0;pl<NDIM;pl++) FULLLOOP jcon[i][j][k][pl]=0.0;
#endif

  return(0);
}

int init_defgrid(void)
{
  // sets metric
  a=0.0;
  // set coordinates
  defcoord=LOGRSINTH;
  // sets parameters of coordinates, default changes
  R0 = 0.0;
  Rin = 0.98 * Rhor;
  Rout = 40.;
  hslope = 0.3;

  return(0);
}


int init_defglobal(void)
{
  int i;

#if(!PRODUCTION)
  debugfail=1;
#else
  debugfail=0; // no messages in production -- assumes all utoprim-like failures need not be debugged
#endif
  // whether to show debug into on failures.  Desirable to turn off if don't care and just want code to burn away using given hacks/kludges
  // 0: no messages
  // 1: critical messages
  // 2: all failure messages

  // included in rdump
  defcon = 1.0;
  /* maximum increase in timestep */
  SAFE=1.3;
  whichrestart = 0;
  restartsteps[0] = 0;
  restartsteps[1] = 0;
  nstep = realnstep = 0;
  failed = 0;
  cour = 0.5;  //atch: modified the courant factor from 0.9

  avgscheme=WENO5BND;
  //  avgscheme=DONOR;
  
  do_transverse_flux_integration = 1;
  do_source_integration = 1;
  do_conserved_integration = 1;

  INVERTFROMAVERAGEIFFAILED = 1;
  LIMIT_AC_PRIM_FRAC_CHANGE = 1;
  MAX_AC_PRIM_FRAC_CHANGE = 0.1;

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  //lim = WENO5FLAT;
  lim = WENO5BND;
  //  lim = WENO3;
  //lim = DONOR;
  //lim = MINM;
  //  lim = PARA;
  //lim = MC;
  //lim = PARA;
  //  lim = PARAFLAT;
  //lim = MC;
  TIMEORDER=4;
  // whether/which ENO used to interpolate fluxes
  DOENOFLUX = ENOFINITEVOLUME;
  //DOENOFLUX= NOENOFLUX;
  //DOENOFLUX=ENOFLUXRECON;
  //  fluxmethod=MUSTAFLUX;
  //fluxmethod=FORCEFLUX;
  fluxmethod=HLLFLUX;
  //fluxmethod=HLLLAXF1FLUX;
  //fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM5D1;  //UTOPRIM2DFINAL;
  //UTOPRIMVERSION=UTOPRIM5D2;
  //  UTOPRIMVERSION=UTOPRIM2DFINAL;
#elif(EOMTYPE==EOMFFDE)
  // PARA and TO=4 and HLL not trustable in FFDE so far
  lim = MC;
  TIMEORDER=2;
  fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM2DFINAL;
  // whether/which ENO used to interpolate fluxes
  DOENOFLUX = ENOFINITEVOLUME;
  //  DOENOFLUX= NOENOFLUX;
  //DOENOFLUX=ENOFLUXRECON;
#endif


  t = 0.;
  tstepparti = t;
  tsteppartf = t;

  tf = 1.0;
  DTd = DTavg = DTdebug = 1.0;
  DTener=1.0;
  DTi=1.0;
  DTr=1;



  GAMMIEDUMP=0;// whether in Gammie output types (sets filename with 3 numbers and doesn't output i,j)
  GAMMIEIMAGE=0; // Gammie filename and single density output
  GAMMIEENER=0; // currently doing gener as well as ener, but this would also output some messages in gammie form

  // DOCOLSPLIT
  //
  // 0: don't ..
  // 1: split dump files into 1 column per file with ending number being column number
  // works in MPI-mode as well.  ROMIO+DOCOLSPLIT is useful for tungsten with low memory and small files to avoid diskspace and memory limits.

  // default
  for(i=0;i<NUMDUMPTYPES;i++){
    DOCOLSPLIT[i]=0;
  }
  // otherwise specify for each dump type


  DODIAGEVERYSUBSTEP = 0;
  DOENODEBUGEVERYSUBSTEP = 0;
 

  DODIAGS=1; // whether to do diagnostics
  // specify individual diagnostics to be done
  DOENERDIAG=1;
  DOGDUMPDIAG=1;
  DORDUMPDIAG=1;
  DODUMPDIAG=1;
  if(DOAVG){
    DOAVGDIAG=1; // choice
  }
  else DOAVGDIAG=0; // no choice
  DOIMAGEDIAG=1;
  DOAREAMAPDIAG=1;

  POSDEFMETRIC=0; // see metric.c, bounds.c, and coord.c

  rescaletype=1;
  // 0: normal
  // 1: extended b^2/rho in funnel
  // 2: conserve E and L along field lines
  //   1 or 2 required to conserve E and L along field line

  /** FIXUP PARAMETERS **/
  RHOMIN=1.e-4;
  UUMIN =1.e-6;
  RHOMINLIMIT=1.e-20;
  UUMINLIMIT =1.e-20;

  // limit of B^2/rho if using that flag
  BSQORHOLIMIT=20.0;
  BSQOULIMIT=100.0;
  GAMMADAMP=5.0;

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  // GODMARK -- unstable beyond about 25, but can sometimes get away with 1000
  GAMMAMAX=25.0; // when we think gamma is just too high and may cause unstable flow, but solution is probably accurate.
#else
  GAMMAMAX=2000.0;
#endif

  GAMMAFAIL=100.0*GAMMAMAX; // when we think gamma is rediculous as to mean failure and solution is not accurate.
  prMAX[RHO]=20.0;
  prMAX[UU]=20.0;
  prMAX[U1]=100.0;
  prMAX[U2]=100.0;
  prMAX[U3]=100.0;
  prMAX[B1]=100.0;
  prMAX[B2]=100.0;
  prMAX[B3]=100.0;

  // some physics
  gam=4./3.; // ultrarelativistic gas, assumes pgas<prad and radiation
	     // doesn't escape
  // gam=5/3 for non-relativistic gas, such as neucleons in collapsar model
  cooling=0;
  // cooling: 0: no cooling 1: adhoc thin disk cooling 2: neutrino cooling for collapsar model

  // boundary conditions (default is for 3D spherical polar grid -- full r,pi,2pi)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;
  BCtype[X3UP]=PERIODIC;
  BCtype[X3DN]=PERIODIC;

  return(0);
}


int init_defconsts(void)
{


  // some constants that define neutrino cooling and spin evolution
  // cgs
  msun=1.989E33;
  lsun=3.89E33;
  G=6.672E-8;
  H=6.6262E-27;
  C=2.99792458E10;
  mn=1.675E-24;
  me=9.10938188E-28;
  kb=1.3807E-16;
  arad=8.*pow(M_PI,5.0)*pow(kb,4.0)/(15*C*C*C*H*H*H);
  sigmasb=arad*C/4.0;
  sigmamat=6.652E-29*100*100;
  mevocsq=1.783E-27;
  ergPmev=1.602E-6;

  // GRB
  mb=mn;
  //  M=3.*msun;
  Mdot=0.1*msun;
  Mdotc=0.35; // approx code units mass accretion rate
  
  // units
  Lunit=G*M/(C*C);
  Tunit=G*M/(C*C*C);
  rhounit=Mdot/(Mdotc*Lunit*Lunit*C);
  Munit=rhounit*pow(Lunit,3.0);
  mdotunit=Munit/Tunit;
  energyunit=Munit*C*C;
  edotunit=energyunit/Tunit;
  Pressureunit=rhounit*C*C;
  Tempunit=Pressureunit*mb/(rhounit*kb);
  Bunit=C*sqrt(rhounit);
  massunitPmsun=Munit/msun;
  
  // physics stuff
  ledd=4.*M_PI*C*G*M*mb/sigmamat;
  leddcode=ledd/edotunit;

 
  return(0);
}






// assumes normal field p
int vpot2field(FTYPE (*A)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE p[][N2M][N3M][NPR])
{
  int i,j,k;
  struct of_geom geom;



  /* flux-ct */

  // A[1] located at CORN1
  // A[2] located at CORN2
  // A[3] located at CORN3

  // F_{\mu\nu} \equiv A_{\nu,\mu} - A_{\mu,\nu}

  // B^i \equiv \dF^{it}

  // F_{ij} = \detg B^k [ijk]

  // F_{\theta\phi} = \detg B^r

  // F_{\phi r} = \detg B^\theta

  // F_{r\theta} = \detg B^\phi


  // \detg B^x = A_{z,y} - A_{y,z}
  // \detg B^y = A_{x,z} - A_{z,x}
  // \detg B^z = A_{y,x} - A_{x,y}

  // loop optimized for filling cells rather than for speed


  FULLLOOP{

    // can do full loop since A[1,2,3] are defined even on their plus parts (redundant but simpler than messing with B's post-hoc by linear extrapolation or something)
    // puts burden on computing A[1,2,3] (such as having radius or anything at plus part)
    // this burden is easier since coord() and bl_coord() have no limits on the i,j,k they are given.
    // so in principle one can give an analytic expression
    // however, depends on metric and such things if converting expression from one basis to another.
    // however, this way is more correct than post-hoc extrapolation.
    // only an issue for fixed boundary conditions, where one could just specify the flux directly.

    get_geometry(i, j, k, CENT, &geom);

    p[i][j][k][B1]  = +(AVGCORN_1(A[3],i,jp1mac(j),k)-AVGCORN_1(A[3],i,j,k))/(geom.g*dx[2]);
    p[i][j][k][B1] += -(AVGCORN_1(A[2],i,j,kp1mac(k))-AVGCORN_1(A[2],i,j,k))/(geom.g*dx[3]);

    p[i][j][k][B2]  = +(AVGCORN_2(A[1],i,j,kp1mac(k))-AVGCORN_2(A[1],i,j,k))/(geom.g*dx[3]);
    p[i][j][k][B2] += -(AVGCORN_2(A[3],ip1mac(i),j,k)-AVGCORN_2(A[3],i,j,k))/(geom.g*dx[1]);

    p[i][j][k][B3]  = +(AVGCORN_3(A[2],ip1mac(i),j,k)-AVGCORN_3(A[2],i,j,k))/(geom.g*dx[1]);
    p[i][j][k][B3] += -(AVGCORN_3(A[1],i,jp1mac(i),k)-AVGCORN_3(A[1],i,j,k))/(geom.g*dx[2]);

  }

  



  

  return(0);
}


int post_init(void)
{


  trifprintf("begin: post_init\n");

  // in synch always here
  if (error_check(3)) {
    dualfprintf(fail_file, "error_check detected failure at main:1\n");
    dualfprintf(fail_file, "Bad initial conditions\n");
    myexit(1);
  }

  // all calculations that do not change for any initial conditions or problem setup or restart conditions

  // these 2 below are used prior, but not initialized otherwised on restart
  // some calculations
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);


  // Setup diagnostics  
  find_horizon();
  setflux();
  if(DOJETDIAG)  setjetflux();
 

  if(RESTARTMODE!=0){
    setfailresponse(restartonfail);
  }


#if(CALCFARADAYANDCURRENTS)

#if(CURRENTST0)
  // setup faraday so t=0 dump has diagnostics
  // may want currents for dump=0 (time-derivative terms will be wrong)
  current_doprecalc(CURTYPET,p);
  current_doprecalc(CURTYPEX,p);
  current_doprecalc(CURTYPEY,p);
  current_doprecalc(CURTYPEZ,p);
  // compute current
  current_calc(cfaraday);

#else

  // need to at least compute t=0 faraday so currents can start to be computed before reach step_ch.c (just an ordering issue with how step_ch.c does this case)
  if(WHICHCURRENTCALC==1){
    // compute faraday components needed for time centering of J
    current_doprecalc(CURTYPET,p);
  }
#endif

#if(FARADAYT0)
  current_doprecalc(CURTYPEFARADAY,p);
#endif

#endif



  trifprintf("end: post_init\n");
  return(0);
}



// after all parameters are set, can call this
int post_par_set(void)
{
  int boundaries_interp_set(void);
  int interp_loop_set(void);
  int orders_set(void);
  void check_bnd_num(void);

  trifprintf("Setting orders\n");
  orders_set();

  trifprintf("Boundary number checks\n");
  check_bnd_num();

  // setup boundary conditions related to interpolation
  trifprintf("Setting boundaries\n");
  boundaries_interp_set();

  trifprintf("Setting interp loop\n");
  interp_loop_set();

  return(0);
}




int orders_set(void)
{
  int l;

  // maximum number of points *potentially* used for a given interpolation type
  // for a point (i), we assume that the first point potentially used is:
  // si = i - (int)(interporder/2)
  // for WENO4, then if i=2, then si=0 and list of potentially used points is 0,1,2,3 .  This is offset such that the i=2 point is correctly centered at the edge with 2 points on edge side.
  // notice that WENO4 only gives correctly centered left state.  Using WENO4, right state is found from right neighbor.
  // WENO4 useful for avg->point at an interface, but not for point->point at different location, which requires odd WENO.
  // for WENO5, then if i=2, then si=0 and list is 0,1,2,3,4 .  Here, point is assumed to be centered of middle (2) cell.
  interporder[DONOR]=1;
  interporder[VANL]=3;
  interporder[MINM]=3;
  interporder[MC]=3;
  interporder[PARA]=5;
  interporder[PARAFLAT]=7;
  interporder[CSSLOPE]=5;
  interporder[WENO3]=3;
  interporder[WENO4]=4;
  interporder[WENO5]=5;
  interporder[WENO5BND]=11;   //replace the 5 for WENO5 with 11 = 5 + (2+2 for reduction) + (1+1 for dp/p reduction) to have enough boundary zones for l&r stencil reduction; also see interpline.c:slope_lim_linetype()
  interporder[WENO5FLAT]=7;
  interporder[WENO6]=6;
  interporder[WENO7]=7;
  interporder[WENO8]=8;
  interporder[WENO9]=9;
  interporder[ENO3]=5;
  interporder[ENO5]=9;

  for(l=0;l<NUMINTERPS;l++){
    if(interporder[l]>MAXSPACEORDER){
      dualfprintf(fail_file,"MAXSPACEORDER=%d while interporder[%d]=%d\n",MAXSPACEORDER,l,interporder[l]);
      myexit(13);
    }
  }

  // check to make sure interpolation type is correctly set
  // checked in flux.c from now on for each type of interpolation
  //  if(DOEXTRAINTERP){
  //    if( (VARTOINTERP!=PRIMTOINTERP_VSQ)||(RESCALEINTERP!=1) ){
  //      dualfprintf(fail_file,"Must set VARTOINTERP==PRIMTOINTERP_VSQ and RESCALEINTERP==1 in definit.h/init.h to use DOEXTRAINTERP=%d\n",DOEXTRAINTERP);
  //      myexit(77);
  //    }
  //  }


  return(0);
}


// check that there are enough boundary zones for interpolation order used
void check_bnd_num(void)
{
  int totalo;
  int get_num_bnd_zones_used(void);
  
  totalo=get_num_bnd_zones_used();

  
  // MAXBND is the number of boundary zones to be set by boundary routines and to be passed by MPI routines
  if(totalo>MAXBND){
    dualfprintf(fail_file,"Not enough MAXBND=%d for avgscheme interporder[%d]=%d or lim interporder[%d]=%d\n",MAXBND,avgscheme,interporder[avgscheme],lim,interporder[lim]);
    myexit(1);
  }

  if(FULLOUTPUT&&USEMPI){
    dualfprintf(fail_file,"Cannot use FULLOUTPUT!=0 when USEMPI=1\n");
    myexit(200);
  }

}


int get_num_bnd_zones_used(void)
{
  int avgo;
  int interpo;
  int totalo;

  if(DOENOFLUX==ENOFINITEVOLUME){
    // number of zones one way for finite volume scheme to convert Uavg -> Upoint
    avgo=(interporder[avgscheme]-1)/2;
  }
  else avgo=0; // no need for extra zones

  // number of zones one way to have for interpolation to get fluxes
  // need to get flux at i,j,k=-1 in any case, and need boundary zones from there, so 1 extra effective boundary zone for interpolation
  interpo=(interporder[lim]-1)/2+1;

  totalo=avgo+interpo;

  return(totalo);
}



// define range of allowed data to access as boundary information
// these should define the maximum amount of information that can be used
int boundaries_interp_set(void)
{
  int itemp;
  int l;



  /////////////////////////
  // x1

  // x1 lower
  if(OUTFLOWAVOIDBC&&(mycpupos[1]==0)&&(BCtype[X1DN]==OUTFLOW)){
    itemp=0;
    ijkminmax[NONENOINTERPTYPE][1][0]=ijkminmaxud[NONENOINTERPTYPE][X1DN]=INFULL1; // no method to avoid boundary data
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][1][0]   =ijkminmaxud[l][X1DN]   =itemp; // avoid boundary data
    }
  }
  else{ // then normal or MPI boundary
    itemp=INFULL1;
    ijkminmax[NONENOINTERPTYPE][1][0]=ijkminmaxud[NONENOINTERPTYPE][X1DN]=itemp;
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][1][0]   =ijkminmaxud[l][X1DN]   =itemp;
    }
  }

  // x1 upper
  if(OUTFLOWAVOIDBC&&(mycpupos[1]==ncpux1-1)&&(BCtype[X1UP]==OUTFLOW)){
    itemp=N1-1;
    ijkminmax[NONENOINTERPTYPE][1][1]=ijkminmaxud[NONENOINTERPTYPE][X1UP]=OUTFULL1; // no method to avoid boundary data
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][1][1]   =ijkminmaxud[l][X1UP]   =itemp; // avoid boundary data
    }
  }
  else{ // then normal or MPI boundary
    itemp=OUTFULL1;
    ijkminmax[NONENOINTERPTYPE][1][1]=ijkminmaxud[NONENOINTERPTYPE][X1UP]=itemp;
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][1][1]   =ijkminmaxud[l][X1UP]   =itemp;
    }
  }

  /////////////////////////////
  // x2

  // x2 lower
  if(OUTFLOWAVOIDBC&&(mycpupos[2]==0)&&(BCtype[X2DN]==OUTFLOW)){
    itemp=0;
    ijkminmax[NONENOINTERPTYPE][2][0]=ijkminmaxud[NONENOINTERPTYPE][X2DN]=INFULL2; // no method to avoid boundary data
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][2][0]   =ijkminmaxud[l][X2DN]   =itemp; // avoid boundary data
    }
  }
  else{ // then normal or MPI boundary
    itemp=INFULL2;
    ijkminmax[NONENOINTERPTYPE][2][0]=ijkminmaxud[NONENOINTERPTYPE][X2DN]=itemp;
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][2][0]   =ijkminmaxud[l][X2DN]   =itemp;
    }
  }

  // x2 upper
  if(OUTFLOWAVOIDBC&&(mycpupos[2]==ncpux2-1)&&(BCtype[X2UP]==OUTFLOW)){
    itemp=N2-1;
    ijkminmax[NONENOINTERPTYPE][2][1]=ijkminmaxud[NONENOINTERPTYPE][X2UP]=OUTFULL2; // no method to avoid boundary data
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][2][1]   =ijkminmaxud[l][X2UP]   =itemp; // avoid boundary data
    }
  }
  else{ // then normal or MPI boundary
    itemp=OUTFULL2;
    ijkminmax[NONENOINTERPTYPE][2][1]=ijkminmaxud[NONENOINTERPTYPE][X2UP]=itemp;
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][2][1]   =ijkminmaxud[l][X2UP]   =itemp;
    }
  }

  /////////////////////////////
  // x3

  // x3 lower
  if(OUTFLOWAVOIDBC&&(mycpupos[3]==0)&&(BCtype[X3DN]==OUTFLOW)){
    itemp=0;
    ijkminmax[NONENOINTERPTYPE][3][0]=ijkminmaxud[NONENOINTERPTYPE][X3DN]=INFULL3; // no method to avoid boundary data
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][3][0]   =ijkminmaxud[l][X3DN]   =itemp; // avoid boundary data
    }
  }
  else{ // then normal or MPI boundary
    itemp=INFULL3;
    ijkminmax[NONENOINTERPTYPE][3][0]=ijkminmaxud[NONENOINTERPTYPE][X3DN]=itemp;
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][3][0]   =ijkminmaxud[l][X3DN]   =itemp;
    }
  }

  // x3 upper
  if(OUTFLOWAVOIDBC&&(mycpupos[3]==ncpux3-1)&&(BCtype[X3UP]==OUTFLOW)){
    itemp=N3-1;
    ijkminmax[NONENOINTERPTYPE][3][1]=ijkminmaxud[NONENOINTERPTYPE][X3UP]=OUTFULL3; // no method to avoid boundary data
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][3][1]   =ijkminmaxud[l][X3UP]   =itemp; // avoid boundary data
    }
  }
  else{ // then normal or MPI boundary
    itemp=OUTFULL3;
    ijkminmax[NONENOINTERPTYPE][3][1]=ijkminmaxud[NONENOINTERPTYPE][X3UP]=itemp;
    for(l=1;l<NUMINTERPTYPES;l++){
      ijkminmax[l][3][1]   =ijkminmaxud[l][X3UP]   =itemp;
    }
  }


  // can add conditions based upon user choice of boundary conditions.
  // the above will determine how far interpolation can grab points
  // needs to consider MPI boundary and normal boundary

	
  return(0);

}


// define range over which various loops go
int interp_loop_set(void)
{
  int dir;
  int avgo;


  // the fluxloop[dir][FIJKDEL]'s can be used for general purpose to get idel, jdel, kdel

  if(DOENOFLUX==ENOFINITEVOLUME){

    // (interporder[avgscheme]-1)/2 is the number of points to the left and the number of ponits to the right that are needed for the finite volume scheme

    // scheme used to convert Uavg -> Upoint requires extra zones
    avgo=(interporder[avgscheme]-1)/2;
    // i.e. if avgscheme=WENO5, then interporder[WENO5]=5 and the Uconsloop goes from -2 .. N+1 inclusive  and fluxes go from -2 to N+2 inclusive
    // same for WENO4
    // this is -avgo .. N-1+avgo for Uconsloop and -avgo .. N+avgo for fluxes

    dir=1;
    fluxloop[dir][FIDEL]=SHIFT1;
    fluxloop[dir][FJDEL]=0;
    fluxloop[dir][FKDEL]=0;
    fluxloop[dir][FFACE]=FACE1;
    fluxloop[dir][FIS]=-avgo*N1NOT1;
    fluxloop[dir][FIE]=N1-1+(avgo+1)*N1NOT1;
    fluxloop[dir][FJS]=INFULL2;
    fluxloop[dir][FJE]=OUTFULL2;
    fluxloop[dir][FKS]=INFULL3;
    fluxloop[dir][FKE]=OUTFULL3;

    dir=2;
    fluxloop[dir][FIDEL]=0;
    fluxloop[dir][FJDEL]=SHIFT2;
    fluxloop[dir][FKDEL]=0;
    fluxloop[dir][FFACE]=FACE2;
    fluxloop[dir][FIS]=INFULL1;
    fluxloop[dir][FIE]=OUTFULL1;
    fluxloop[dir][FJS]=-avgo*N2NOT1;
    fluxloop[dir][FJE]=N2-1+(avgo+1)*N2NOT1;
    fluxloop[dir][FKS]=INFULL3;
    fluxloop[dir][FKE]=OUTFULL3;

    dir=3;
    fluxloop[dir][FIDEL]=0;
    fluxloop[dir][FJDEL]=0;
    fluxloop[dir][FKDEL]=SHIFT3;
    fluxloop[dir][FFACE]=FACE3;
    fluxloop[dir][FIS]=INFULL1;
    fluxloop[dir][FIE]=OUTFULL1;
    fluxloop[dir][FJS]=INFULL2;
    fluxloop[dir][FJE]=OUTFULL2;
    fluxloop[dir][FKS]=-avgo*N3NOT1;
    fluxloop[dir][FKE]=N3-1+(avgo+1)*N3NOT1;

    // loop over averaged U to get Uf
    //    Uconsloop[FIDEL]=0;
    //    Uconsloop[FJDEL]=0;
    //    Uconsloop[FKDEL]=0;
    Uconsloop[FFACE]=CENT;
    Uconsloop[FIS]=-avgo*N1NOT1;
    Uconsloop[FIE]=N1-1+avgo*N1NOT1;
    Uconsloop[FJS]=-avgo*N2NOT1;
    Uconsloop[FJE]=N2-1+avgo*N2NOT1;
    Uconsloop[FKS]=-avgo*N3NOT1;
    Uconsloop[FKE]=N3-1+avgo*N3NOT1;

    // inversion for this method just involves CZLOOP

  }
  else if((DOENOFLUX==NOENOFLUX)||(DOENOFLUX==ENOFLUXRECON)||(DOENOFLUX==ENOFLUXSPLIT)){

    // Uconsloop for these methods just involve normal CZLOOP
    // inversion for this method just involves CZLOOP
    dir=1;
    fluxloop[dir][FIDEL]=SHIFT1;
    fluxloop[dir][FJDEL]=0;
    fluxloop[dir][FKDEL]=0;
    fluxloop[dir][FFACE]=FACE1;
    fluxloop[dir][FIS]=0;
    fluxloop[dir][FIE]=OUTM1;
    fluxloop[dir][FJS]=-SHIFT2;  //atch: loop over additional row to provide enough fluxes for FLUXCT, etc. to operate near the boundary
    fluxloop[dir][FJE]=N2-1+SHIFT2; // " " 
    fluxloop[dir][FKS]=-SHIFT3;     // " "
    fluxloop[dir][FKE]=N3-1+SHIFT3; // " "

    dir=2;
    fluxloop[dir][FIDEL]=0;
    fluxloop[dir][FJDEL]=SHIFT2;
    fluxloop[dir][FKDEL]=0;
    fluxloop[dir][FFACE]=FACE2;
    fluxloop[dir][FIS]=-SHIFT1;   //atch: loop over additional row to provide enough fluxes for FLUXCT, etc. to operate near the boundary
    fluxloop[dir][FIE]=N1-1+SHIFT1; // " "
    fluxloop[dir][FJS]=0; 
    fluxloop[dir][FJE]=OUTM2;
    fluxloop[dir][FKS]=-SHIFT3;    // " "
    fluxloop[dir][FKE]=N3-1+SHIFT3;// " "

    dir=3;
    fluxloop[dir][FIDEL]=0;
    fluxloop[dir][FJDEL]=0;
    fluxloop[dir][FKDEL]=SHIFT3;
    fluxloop[dir][FFACE]=FACE3;
    fluxloop[dir][FIS]=-SHIFT1;   //atch: loop over additional row to provide enough fluxes for FLUXCT, etc. to operate near the boundary
    fluxloop[dir][FIE]=N1-1+SHIFT1;  // " "
    fluxloop[dir][FJS]=-SHIFT2;      // " "
    fluxloop[dir][FJE]=N2-1+SHIFT2;  // " "
    fluxloop[dir][FKS]=0;
    fluxloop[dir][FKE]=OUTM3;

  }

  return(0);

}










int find_horizon(void)
{
  int i, j, k, ii;
  FTYPE r1, r2;
  FTYPE X[NDIM],V[NDIM];
  int horizoncpupos1, gotit;
  FTYPE horizonvalue;
  // called after grid is setup for all cpus


  trifprintf("begin: find_horizon ... ");

  // find cpu column that brackets the horizon and determine the
  // i-offset of horizon

  horizonvalue = Rhor;
  horizoni = -100;
  gotit = 0;
  for (ii = numprocs - 1; ii >= 0; ii--) { // should get done by first row
    if (ii == myid) {
      for (i = N1 - 1; i >= 0; i--) {
        j = N2 / 2;             // doesn't matter
	k = N3 / 2; // doesn't matter
        coord(i, j, k, CENT, X);
        bl_coord(X, V);
	r1=V[1];
        coord(ip1, j, k, CENT, X);
        bl_coord(X, V);
	r2=V[1];
        if (fabs(r1 - horizonvalue) <= (r2 - r1)) {     // find horizon
          horizoni = i;
          horizoncpupos1 = mycpupos[1];
          break;
        }
      }
    }
    if (numprocs > 0) {
#if(USEMPI)
      MPI_Bcast(&horizoni, 1, MPI_INT, ii, MPI_COMM_WORLD);
      MPI_Bcast(&horizoncpupos1, 1, MPI_INT, ii, MPI_COMM_WORLD);
#endif
    }
    if (horizoni >= 0)
      gotit = 1;                // can stop entire process
    if (mycpupos[1] != horizoncpupos1) {
      horizoni = -100;
    }                           // reset if not right cpu group
    if (gotit)
      break;
  }
  trifprintf("horizoni: %d horizoncpupos1: %d\n", horizoni,
             horizoncpupos1);
  // just a check
  dualfprintf(log_file,"horizoni: %d mycpupos[1]: %d horizoncpupos1: %d\n", horizoni, mycpupos[1], horizoncpupos1);

  trifprintf("end: find_horizon\n");
  return(0);
}


// determine if this cpu is doing what flux
int setflux(void)
{
  int dir;

  // setup pointers
  enerpos=enerposreg[GLOBALENERREGION];
  doflux=dofluxreg[GLOBALENERREGION];

  // all CPUs , total space for global region
  // inclusive loop range (0..N-1)
  enerpos[X1DN]=0;
  enerpos[X1UP]=N1-1;
  enerpos[X2DN]=0;
  enerpos[X2UP]=N2-1;
  enerpos[X3DN]=0;
  enerpos[X3UP]=N3-1;


  // only 0 through N-1 mean do flux

  // doflux<0 means ignore this cpu in flux calculation
  // doflux>=0 is used to set where flux is computed.
  // If dimension exists (i.e. N>1), then OUTM = N is outer flux on boundary
  // If dimension doesn't exist, then no such outer flux or outer boundary, so stay on same boundary (e.g., i=0) rather than i=N.

  // x1
  if((N1>1)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    trifprintf("proc: %d doing flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if((N1>1)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    trifprintf("proc: %d doing flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  // x2
  if((N2>1)&&(mycpupos[2]==0)){
    doflux[X2DN]=0;
    trifprintf("proc: %d doing flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if((N2>1)&&(mycpupos[2]==ncpux2-1)){
    doflux[X2UP]=OUTM2;
    trifprintf("proc: %d doing flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;

  // x3
  if((N3>1)&&(mycpupos[3]==0)){
    doflux[X3DN]=0;
    trifprintf("proc: %d doing flux X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  if((N3>1)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    trifprintf("proc: %d doing flux X3UP\n",myid);
  }
  else doflux[X3UP]=-100;

  // fluxes are on edges of zone, so 0 and N are on edge fluxes

  DIRLOOP(dir) trifprintf("proc: %d %d doflux[%d]=%d enerpos[%d]=%d\n",myid,GLOBALENERREGION,dir,doflux[dir],dir,enerpos[dir]);


  return(0);
}

// GODMARK
// the below only works for a grid with full Pi.  Crashes code at runtime otherwise!  should fix.

// assume this is an intrinsically axisymmetric function, so k (\phi) is just carried along -- no truncation in \phi

// determine the flux positions for each CPU for the jet region (jetedge)
// AND the range of volume integration for cons. variables in jet (jetpos).
int setjetflux(void)
{
  FTYPE X[NDIM],V[NDIM],r,th;
  int i,j,k,pl;
  int dir;
  FTYPE startth,endth,thetajet;
  int jetedge[NUMJETS];


  if(defcoord==EQMIRROR){
    dualfprintf(fail_file,"setjetflux() not setup to work for non-full-Pi grids\n");
    myexit(1);
  }


  // jet region is assumed to be within a constant theta slice
  // this is theta w.r.t. polar axis
  thetajet=M_PI*0.5-h_over_r_jet;
  // find j for which theta is at our prespecified point

  i=0;j=0;k=0;
  coord(i, j, k, FACE2, X);
  bl_coord(X, V);
  startth=V[2];
  i=0;j=N2;k=0;
  coord(i, j, k, FACE2, X);
  bl_coord(X, V);
  endth=V[2];

  // assumes 0<thetajet<Pi/2
  if((fabs(startth-thetajet)/thetajet<1E-8)||
     (fabs(endth-thetajet)/thetajet<1E-8)||
     (fabs(startth-(M_PI-thetajet))/thetajet<1E-8)||
     (fabs(endth-(M_PI-thetajet))/thetajet<1E-8)
     ){
    dualfprintf(fail_file,"thetajet is on top of grid, move h_over_r_jet a bit\n");
    myexit(1);
  }


  ////////////////////
  //
  // INNERJET
  //

  // setup pointers
  enerpos=enerposreg[INNERJETREGION];
  doflux=dofluxreg[INNERJETREGION];


  // see if jet edge is related to this CPU
  // assumes increasing j is increasing th
  if((startth<=thetajet)&&(endth<=thetajet)){
    // if cpu entirely within inner theta jet
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
    jetedge[INNERJET]=-100;
  }
  else if((startth<thetajet)&&(endth>thetajet)){
    // if inner jet edge is on this CPU but not on boundary
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    i=0;
    for(j=0;j<=OUTM2;j++){
      coord(i, j, k, FACE2, X);
      bl_coord(X,V);
      r=V[1];
      th=V[2];
      // look for switch from below to above thetajet at inner theta jet edge
      if(th>thetajet){
	enerpos[X2UP]=jm1;
	jetedge[INNERJET]=j;
	break;
      }
    }
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
  }
  else if((startth>=thetajet)&&(endth>=thetajet)){
    // if cpu is entirely not contained in inner jet
    enerpos[X1DN]=-100;
    enerpos[X1UP]=-100;
    enerpos[X2DN]=-100;
    enerpos[X2UP]=-100;
    enerpos[X3DN]=-100;
    enerpos[X3UP]=-100;
    jetedge[INNERJET]=-100;
  }
  else{
    trifprintf("problem with INNERJET setjetflux()\n");
    myexit(1);
  }



  // left edge (any directional condition would do)
  if((N1>1)&&(enerpos[X1DN]!=-100)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    trifprintf("proc: %d doing inner jet flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  // right edge (any directional condition would do)
  if((N1>1)&&(enerpos[X1DN]!=-100)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    trifprintf("proc: %d doing inner jet flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  // lower theta boundary
  if((N2>1)&&(mycpupos[2]==0)){
    doflux[X2DN]=0;
    trifprintf("proc: %d doing inner jet flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;
  
  // upper theta boundary
  if((N2>1)&&(jetedge[INNERJET]!=-100)){ // only get flux if CPU has edge
    doflux[X2UP]=jetedge[INNERJET];
    trifprintf("proc: %d doing inner jet flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;

  if((N3>1)&&(enerpos[X3DN]!=-100)&&(mycpupos[3]==0)){
    doflux[X3DN]=0;
    trifprintf("proc: %d doing inner jet flux X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  // right edge (any directional condition would do)
  if((N3>1)&&(enerpos[X3DN]!=-100)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    trifprintf("proc: %d doing inner jet flux X3UP\n",myid);
  }
  else doflux[X3UP]=-100;



  DIRLOOP(dir) trifprintf("proc: %d %d doflux[%d]=%d enerpos[%d]=%d\n",myid,INNERJETREGION,dir,doflux[dir],dir,enerpos[dir]);


  /////////////////////
  //
  // OUTERJET
  //

  // setup pointers
  enerpos=enerposreg[OUTERJETREGION];
  doflux=dofluxreg[OUTERJETREGION];


  // see if outer jet edge is related to this CPU
  if((startth<=M_PI-thetajet)&&(endth<=M_PI-thetajet)){
    // if cpu entirely not within outer jet region
    enerpos[X1DN]=-100;
    enerpos[X1UP]=-100;
    enerpos[X2DN]=-100;
    enerpos[X2UP]=-100;
    enerpos[X3DN]=-100;
    enerpos[X3UP]=-100;
    jetedge[OUTERJET]=-100;
  }
  else if((startth<M_PI-thetajet)&&(endth>M_PI-thetajet)){
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    // if outer jet edge is on this CPU but not on boundary
    i=0;k=0;
    for(j=0;j<=OUTM2;j++){
      coord(i, j, k, FACE2, X);
      bl_coord(X, V);
      th=V[2];
      // look for switch from below to above thetajet at inner theta jet edge
      if(th>M_PI-thetajet){
	enerpos[X2DN]=jm1;
	jetedge[OUTERJET]=jm1;
	break;
      }
    }
    enerpos[X2UP]=N2-1;

    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;

  }
  else if((startth>=M_PI-thetajet)&&(endth>=M_PI-thetajet)){
    // if cpu is entirely containe within outer jet
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
    jetedge[OUTERJET]=-100;
  }
  else{
    trifprintf("problem with OUTERJET setjetflux()\n");
    myexit(1);
  }

  if((N1>1)&&(enerpos[X1DN]!=-100)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    trifprintf("proc: %d doing outer jet flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if((N1>1)&&(enerpos[X1DN]!=-100)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    trifprintf("proc: %d doing outer jet flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  if((N2>1)&&(jetedge[OUTERJET]!=-100)){
    doflux[X2DN]=jetedge[OUTERJET];
    trifprintf("proc: %d doing outer jet flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if((N2>1)&&(mycpupos[2]==ncpux2-1)){
    doflux[X2UP]=OUTM2;
    trifprintf("proc: %d doing outer jet flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;
  // fluxes are on edges of zone, so 0 and N are on edge fluxes

  if((N3>1)&&(enerpos[X3DN]!=-100)&&(mycpupos[3]==0)){
    doflux[X3DN]=0; 
    trifprintf("proc: %d doing outer jet flux X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  if((N3>1)&&(enerpos[X3DN]!=-100)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    trifprintf("proc: %d doing outer jet flux X3UP\n",myid);
  }
  else doflux[X3UP]=-100;


  DIRLOOP(dir) trifprintf("proc: %d %d doflux[%d]=%d enerpos[%d]=%d\n",myid,OUTERJETREGION,dir,doflux[dir],dir,enerpos[dir]);


  return(0);
}


int init_dumps(void)
{
  struct blink * blinkptr;
  struct blink * cpulinkptr;
  int i,numlists,numcells;
  int maxnumcolumns;


  trifprintf("begin: init_dumps\n");


  ///////////////////////////
  //
  // setup number of columns per dump file (see dumpgen.c or dump.c for how used)
  //
  /////////////////////////

  // now setup the data output/input organization for chunking method for each number of columns
  dnumcolumns[IMAGECOL]=1;

  dnumcolumns[RDUMPCOL]=NPRDUMP*2; // primitives and conservatives
  //dnumcolumns[RDUMPCOL]=NPRDUMP; // primitives only


  if(GAMMIEDUMP)  dnumcolumns[DUMPCOL]=2*3 + NPRDUMP*2 + 1 + NDIM * NDIM + 6 + 1
#if(CALCFARADAYANDCURRENTS)
		    + NDIM*2
		    + 2*6
#endif
		    ;
  else  dnumcolumns[DUMPCOL]=3*3 + NPRDUMP*2 + 1 + NDIM * NDIM + 6 + 1 
#if(CALCFARADAYANDCURRENTS)
	  + NDIM*2
	  + 2*6
#endif
	  ;    // 61 total if also doing currents and faraday, 41 otherwise


  // 205+4+4*4 currently
  //dnumcolumns[GDUMPCOL]=3*3+NDIM*NDIM*NDIM+NPG*NDIM*NDIM*2+NPG+4+4*4;
  //NPG was replaced with unity in order to avoid excessive dumping of info (only center info now available)
  dnumcolumns[GDUMPCOL]=3*3  +   NDIM*NDIM*NDIM  +   1*NDIM*NDIM*2   +   1  +  NDIM   +   NDIM*NDIM;
  //t^i x^i V^i,     \Gamma^\mu_{\nu\tau},     g^{\mu\nu} g_{\mu\nu}, \sqrt{-g}, \gamma_\mu, dx^\mu/dx^\nu
  


  // 36+29+8*2+4*2+2+12*2+96*2=339
  dnumcolumns[AVGCOL]=3*3 + 1 + NUMNORMDUMP  // (6+1+29=36)
    + NUMNORMDUMP // |normal terms| (29)
#if(CALCFARADAYANDCURRENTS)
    + NDIM*2 // jcon/jcov (8)
    + NDIM*2 // |jcon|/|jcov| (8)
#endif
    + NDIM*2 // massflux/|massflux|
    + NUMOTHER*2 // other stuff and fabs of each
#if(CALCFARADAYANDCURRENTS)
    +6*2 // fcon/fcov (12)
    +6*2 // |fcon|,|fcov| (12)
#endif
    +7*16 // Tud all 7 parts, all 4x4 terms (112)
    +7*16 // |Tud| all 7 parts, all 4x4 terms (112)
    ;


  if(DOAVG2){
    dnumcolumns[AVGCOL]-=224;
    dnumcolumns[AVG2COL]=10 + 224; // otherwise doesn't exist so don't need to set
  }
  else dnumcolumns[AVG2COL]=0;


  
  if(DODEBUG){
    dnumcolumns[DEBUGCOL]=NUMFAILFLOORFLAGS*NUMTSCALES;
  }
  else dnumcolumns[DEBUGCOL]=0;

  if(DOENODEBUG){
    //dnumcolumns[ENODEBUGCOL]=NUMENODEBUGS;
    dnumcolumns[ENODEBUGCOL]=(3-1)* NUMINTERPTYPES * (NPR-4) * NUMENODEBUGS;  //SASMARK2
  }
  else dnumcolumns[ENODEBUGCOL]=0;

  if(DOFIELDLINE){
    dnumcolumns[FIELDLINECOL]=NUMFIELDLINEQUANTITIES;
  }
  else dnumcolumns[FIELDLINECOL]=0;


  if(DODISS){
    dnumcolumns[DISSDUMPCOL]=NUMDISSFUNPOS;
  }
  else dnumcolumns[DISSDUMPCOL]=0;


  trifprintf("dump number of columns\n\n0=IMAGE, 1=RDUMP, 2=DUMP, 3=GDUMP, 4=AVG, 5=AVG2, 6=DEBUG 7=FIELDLINE 8=ENODEBUG 9=DISSDUMPCOL (see global.h)\n");
  for(i=0;i<NUMDUMPTYPES;i++) trifprintf("dnumcolumns[%d]=%d\n",i,dnumcolumns[i]);








  /////////////////////
  //
  // Below shouldn't be modified for user purposes
  //
  // setup number of buffers
  //
  ///////////////////////

  maxnumcolumns=0;
  for(i=0;i<NUMDUMPTYPES;i++){
    if(maxnumcolumns<dnumcolumns[i]) maxnumcolumns=dnumcolumns[i];
  }
  // buffer must at least hold maxcolumns of data, and since buffer is only N1*N2*N3 big, make sure that at least NUMBUFFERS*N1*N2*N3>maxnumcolumns
  if(N1*N2*N3<maxnumcolumns) NUMBUFFERS=(int)ceil((FTYPE)maxnumcolumns/((FTYPE)(N1*N2*N3)));
  else NUMBUFFERS=1;


  for(i=0;i<NUMDUMPTYPES;i++) if(dnumcolumns[i]>0) setuplinklist(dnumcolumns[i],i);


  trifprintf("end setuplinklists: %d\n",NUMDUMPTYPES);



  ///////////////////////////
  //
  // setup link lists for setuplinklist()
  //
  ///////////////////////////

  trifprintf("start per cpu lists\n");
  // check link lists
  for(i=0;i<NUMDUMPTYPES;i++){
    if(dnumcolumns[i]>0){
      fprintf(log_file,"i=%d\n",i); fflush(log_file);
      blinkptr=blinkptr0[i];
      numlists=0;
      numcells=0;
      while(blinkptr!=NULL){
	numcells+=blinkptr->num;
	//      fprintf(log_file,"i=%d num=%d, numtotal=%d\n",i,blinkptr->num,numcells); fflush(log_file);
	numlists++;
	blinkptr=blinkptr->np; // next one
      }
      fprintf(log_file,"i=%d numlists=%d numcells=%d\n",i,numlists,numcells);
      numlists=0;
    }
  }

  // check cpu=0 link list
  if(myid==0){
    trifprintf("start cpu==0 lists\n");
    for(i=0;i<NUMDUMPTYPES;i++){
      if(dnumcolumns[i]>0){
	fprintf(log_file,"i=%d\n",i); fflush(log_file);
	cpulinkptr=cpulinkptr0[i];
	numlists=0;
	numcells=0;
	while(cpulinkptr!=NULL){
	  numcells+=cpulinkptr->num;
	  //	fprintf(log_file,"i=%d num=%d, cpu=%d, li=%d, lj=%d, lk=%d, col=%d, numtotal=%d\n",i,cpulinkptr->num,cpulinkptr->cpu,cpulinkptr->i,cpulinkptr->j,cpulinkptr->k,cpulinkptr->col,numcells); fflush(log_file);
	  numlists++;
	  cpulinkptr=cpulinkptr->np; // next one
	}
	fprintf(log_file,"i=%d numlists=%d numcells=%d\n",i,numlists,numcells);
	numlists=0;
      }
    }
  }


  trifprintf("end: init_dumps\n");


  return(0);
}


// setuplinklist() not a user function

// use for mpiminio() functions (see init_mpi.c)
int setuplinklist(int numcolumns,int which)
{
  int gcount,lcount,numlinks;
  int i,j,k,col,li,lj,lk,pi,pj,pk,pid,firstlink;
  struct blink * clinkptr0, *clinkptr;
  struct blink * linkptr0for0, *linkptrfor0;
  int *lcountfor0;
  int firstlinkfor0;
  int *firstlijk,*li0,*lj0,*lk0,*lcol0;
  int ri,rj,rk,rcol;
  int *cpulist0;
  int numcpusinlist0,lcpu,itercpu,buffersize;
  int maxnumsize;


  maxnumsize=(int)(ceil(ceil((FTYPE)(N1*N2*N3*NUMBUFFERS)/(FTYPE)numcolumns)*(FTYPE)(numcolumns)));

  if(myid==0){
    // cpulist0's size is maximum possible number of cpus in a list due to buffer size
    //    buffersize=(int)(ceil(ceil((FTYPE)(N1*N2*N3*NUMBUFFERS)/(FTYPE)numcolumns)*(FTYPE)(numcolumns)/(FTYPE)N1));
    buffersize=numprocs;
    fprintf(stderr,"max cpus in a list=%d\n",buffersize); fflush(stderr);
    if((cpulist0=(int*)malloc(sizeof(int)*buffersize))==NULL){
      dualfprintf(fail_file,"can't allocate cpulist0\n");
      myexit(10000);
    }
    if((lcountfor0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lcountfor0\n");
      myexit(10000);
    }
    if((firstlijk=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate firstlijk\n");
      myexit(10000);
    }
    if((li0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate li0\n");
      myexit(10000);
    }
    if((lj0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lj0\n");
      myexit(10000);
    }
    if((lk0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lk0\n");
      myexit(10000);
    }
    if((lcol0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lcol0\n");
      myexit(10000);
    }
    for(i=0;i<buffersize;i++){
      cpulist0[i]=0;
    }
    for(i=0;i<numprocs;i++){
      lcountfor0[i]=firstlijk[i]=li0[i]=lj0[i]=lk0[i]=lcol0[i]=0;
    }
  }



  numcpusinlist0=0;

  clinkptr0=NULL;
  gcount=0;
  lcount=0;
  numlinks=0;
  firstlink=1;
  if(myid==0){
    for(itercpu=0;itercpu<numprocs;itercpu++){  firstlijk[itercpu]=1; }
    linkptr0for0=NULL;
    firstlinkfor0=1;
  }

  /////////////////////////
  // general loop
  for(k=0;k<ncpux3*N3;k++)  for(j=0;j<ncpux2*N2;j++)  for(i=0;i<ncpux1*N1;i++) for(col=0;col<numcolumns;col++){
    // relative local index
    li=i%N1;
    lj=j%N2;
    lk=k%N3;
    // cpu position number
    pi=(int)(i/N1);
    pj=(int)(j/N2);
    pk=(int)(k/N3);
    // cpu id for this data
    pid=pk*ncpux2*ncpux1+pj*ncpux1+pi;
    if(myid==pid) lcount++;
    if(myid==0){
      lcountfor0[pid]++;
      // below is if we have need this cpu's data (pid) and need to mark starting point on full grid
      if(firstlijk[pid]){
	cpulist0[numcpusinlist0++]=pid;
	li0[pid]=i;
	lj0[pid]=j;
	lk0[pid]=k;
	lcol0[pid]=col;
	if(col!=0){
	  dualfprintf(fail_file,"col!=0 col=%d, so chunking bad\n",col);
	  myexit(10000);
	}
	firstlijk[pid]=0;
      }
    }
    gcount++;
    //    if(myid==0){
    //  fprintf(fail_file,"%d %d %d %d\n",numcpusinlist0,gcount,pid,cpulist0[numcpusinlist0]); fflush(fail_file);
    // }
    //    fprintf(log_file,"%d %d %d %d %d %d %d %d\n",li,lj,lk,pi,pj,pk,pid,lcount,gcount); fflush(log_file);
    // 1st below if is to catch every buffer amount, while 2nd if part is needed to account for when the number of buffers is such that the last buffer isn't completely needed
    // this should work for any numcolumns or NUMBUFFERS, even at very last zone no matter what
    // chunk in minimum size of numcolumns
    if((gcount%((int)ceil((double)(N1*N2*N3*NUMBUFFERS/numcolumns))*numcolumns)==0)||(gcount==totalzones*numcolumns)){
      // ok, so numcolumns can't exceed the buffer size, highly unlikely to happen, and checked for!
      if(myid==0){
	// must do in order determined to have data, not numerical order
	for(itercpu=0;itercpu<numcpusinlist0;itercpu++){
	  lcpu=cpulist0[itercpu];
	  if(lcountfor0[lcpu]>0){
	    if(itercpu==0){ // first cpu in list
	      ri=li0[lcpu];
	      rj=lj0[lcpu];
	      rk=lk0[lcpu];
	      rcol=lcol0[lcpu];
	    }
	    if(firstlinkfor0){
	      linkptrfor0=linkptr0for0=addlink(NULL);
	      firstlinkfor0=0;
	    }
	    else{
	      linkptrfor0=addlink(linkptrfor0);
	    }
	    linkptrfor0->cpu=lcpu;
	    linkptrfor0->num=lcountfor0[lcpu];
	    linkptrfor0->i=li0[lcpu];
	    linkptrfor0->j=lj0[lcpu];
	    linkptrfor0->k=lk0[lcpu];
	    linkptrfor0->col=lcol0[lcpu];
	    linkptrfor0->ri=ri;
	    linkptrfor0->rj=rj;
	    linkptrfor0->rk=rk;
	    linkptrfor0->rcol=rcol;
	    linkptrfor0->end=0;
	    
	    lcountfor0[lcpu]=0; // reset counter for this id
	    firstlijk[lcpu]=1; // reset starting value
	  }
	  else{
	    fprintf(fail_file,"wtf: shoudn't be here\n");
	    myexit(10000);
	  }
	}
	// the last link is here identified as the last in the series of cpus to communicate with.  There's at least one new link here!
	linkptrfor0->end=1;
	numcpusinlist0=0; // reset list of cpus for this list
      }
      if(lcount>0){
	fprintf(log_file,"numcolumns=%d lcount=%d\n",numcolumns,lcount); fflush(log_file);
        // initialize another structure
        // set previous structure value to this structure, set this next one to NULL
        if(firstlink){
          clinkptr=clinkptr0=addlink(NULL);
	  clinkptr->num=lcount;
          firstlink=0;
        }
        else{
          clinkptr=addlink(clinkptr);
	  clinkptr->num=lcount;
        }
        lcount=0;
      }
    }// otherwise continue
  }      // now we have a link list for each cpu that determines how long each next buffer is that needs to be sent to cpu=0
  blinkptr0[which]=clinkptr0;
  cpulinkptr0[which]=linkptr0for0;

  return(0);
}

// add link for forward-only link list
struct blink * addlink(struct blink * clinkptr)
{
  struct blink *pb;

  pb=(struct blink *)malloc(sizeof(struct blink));
  pb->np=NULL; // terminate list
  // set last link's pointer to this new structure
  if(clinkptr!=NULL) clinkptr->np=pb;

  return(pb);
}



int pi2Uavg(FTYPE (*prim)[N2M][N3M][NPR], FTYPE (*Upoint)[N2M][N3M][NPR], FTYPE (*Uavg)[N2M][N3M][NPR])
{
  struct of_geom geom;
  struct of_state q;
  int i,j,k;
  extern void avg2cen_interp(int whichquantity, int interporflux,FTYPE (*prims_from_avg_cons)[N2M][N3M][NPR], FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR]);
  int pl;
  FTYPE (*Ucomputed)[N2M][N3M][NPR];
  FTYPE (*Ustored)[N2M][N3M][NPR];



  FULLLOOP{
    // find dU(pb)
    // only find source term if non-Minkowski and non-Cartesian
    // set geometry for centered zone to be updated
    get_geometry(i, j, k, CENT, &geom);
    

    // find Ui(pi)
    MYFUN(get_state(prim[i][j][k], &geom, &q),"initbasec:pi2Uavg()", "get_state()", 1);
    MYFUN(primtoU(UEVOLVE,prim[i][j][k], &q, &geom, Upoint[i][j][k]),"initbase.c:pi2Uavg()", "primtoU()", 1);

  }

  // volume integrate dUgeom
  if(DOENOFLUX == ENOFINITEVOLUME && avgscheme > 3 && do_conserved_integration ) {  //SASMARK: added do_conserved_integration so do not do c2a on conserved when conserved integratin is off
    // c2a_1 c2a_2 c2a_3
    avg2cen_interp(ENOCONSERVED, ENOCENT2AVGTYPE, prim, Upoint, Uavg);
  }
  else {
    FULLLOOP PLOOP(pl) {
      Uavg[i][j][k][pl] = Upoint[i][j][k][pl];
    }
  }
				    

  
  // volume integrate dUgeom  
  // c2a_1 c2a_2 c2a_3
  //avg2cen_interp(ENOCONSERVED, ENOCENT2AVGTYPE, prim, Upoint, Uavg);  //SASMARK: repetition?? commented out

  // need to copy Uavg -> unew for diagonistics
  FULLLOOP PLOOP(pl){
    unew[i][j][k][pl] = Uavg[i][j][k][pl];
  }

  // now both unew and uinitial set correctly

  return(0);
}


void makedirs(void)
{

  if ((USEMPI == 0) || (USEMPI && (!USEGM))) {
    if ((mpicombine && (myid == 0)) || (mpicombine == 0)) {
      system("mkdir dumps");
      system("mkdir images");
    }
#if(USEMPI)
    MPI_Barrier(MPI_COMM_WORLD);	// all cpus wait for directory
    // to be created
#endif
  }
}

