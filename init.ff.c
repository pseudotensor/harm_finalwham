//B3 apparently gets overwritten by the expression from differentiating the vector potential

/*

to use this, need to probably change following elsewhere:

1) AVOIDCS in phys.ffde.c
2) EOMTYPE in global.h
3) MCOORD KSCOORDS in metric.h
4) probably defcoord==0's entries in coord.c
5) maybe show boundary zones by modifying dump.c
6) ZEROOUTFLOWFLUX 1 in step_ch.h ?


*/

#include "decs.h"


#define SLOWFAC 1.0		/* reduce u_phi by this amount */

static SFTYPE rhomax=0,umax=0,bsq_max=0,beta,rin;

static FTYPE mm;

#define NTHETAMAX 10005

static FTYPE Thetavstheta[NTHETAMAX],dThetadthetavstheta[NTHETAMAX],theta[NTHETAMAX],dtheta,myfloati,myTheta,myThetac,mydThetac;
static int NTHETA;

static FTYPE  INDEXN;  //SASMARKX:  changed to FTYPE from int, otherwise INDEXN = (1.0/4.0) would be truncated to zero


// 0 : split monopole
// 1 : monopole
// 2 : Wald
// 3 : BZ paraboloidal
// 4 : GRMHD nearly-paraboloidal
// 5 : BZ para but monopole near horizon
// 6 : GRMHD nearly-paraboloidal but monopole near horizon
// 7 : constant velocity Ramesh disk // can do cylindrical coordinatse also
#define PROBLEMTYPE 7

#define RAMESHTYPE 1  // 1 -- nu = 1, M = 0.25, s = 0.5 solution



int pre_init_specific_init(void)
{
  void get_nearlypara_data(void);
  void get_ramesh_data(void);
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.2;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  Omegastar=0;




  if((PROBLEMTYPE==4)||(PROBLEMTYPE==6)){
    get_nearlypara_data();
  }
  else if(PROBLEMTYPE==7){
    get_ramesh_data();
  }


  return(0);
}

int init_conservatives(FTYPE p[][N2M][N3M][NPR], FTYPE Utemp[][N2M][N3M][NPR], FTYPE U[][N2M][N3M][NPR])
{
  extern int pi2Uavg(FTYPE (*prim)[N2M][N3M][NPR], FTYPE (*Upoint)[N2M][N3M][NPR], FTYPE (*Uavg)[N2M][N3M][NPR]);

  trifprintf("begin init_conservatives\n");
  pi2Uavg(p, Utemp, U);
  trifprintf("end init_conservatives\n");

  return(0);

}


int post_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *after* all other calculations
	extern int set_dt(FTYPE (*prim)[N2M][N3M][NPR], SFTYPE *dt);  //atch update
	
	set_dt( p, &dt);  //atch update: always use dt determined by the courant factor

  pdump = p;
  // pre-bound dump file
  if (dump(9999) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }

  return(0);
}

int init_grid(void)
{
  SFTYPE rh;


#if(MCOORD==CYLMINKMETRIC)

  //  Zin=-1E3;
  //  Zout=1E3;
  //  Rin=1.0;
  //  Rout=1E3;

  if(0){
    //Zeqin=1E-5;
    Zeqin=0.0;
    Zin=-1E2;
    Zout=1E2;
    Rin=20.0;
    //  Rin=1E-5;
    //  Rin=1E-2;
    //  Rin=1.0;
    Rout=1E3;
  }

  if(1){
    Zeqin=0;
    Zin=-1E3;
    Zout=1E3;
    Rin=1E-4;
    Rout=1E3;
  }


  defcoord=666; // see coord.c for npow parameter


#elif(MCOORD==SPCMINKMETRIC)


  if(0){
    // metric stuff first
    a = 0.0 ;

    // make changes to primary coordinate parameters R0, Rin, Rout, hslope
    R0 = -1.0;
    Rhor=rhor_calc(0);
    Rout = 1E3;

    //  Rin=setRin(setihor());
    Rin = 0.98 * Rhor;

    hslope = 0.1;

    // define coordinate type
    defcoord = 9;
  }
  else if(1){
    // make changes to primary coordinate parameters R0, Rin, Rout, hslope
    Rhor=rhor_calc(0);

    if(0){
      R0 = -1.0;
      Rout = 1E3;
    }
    else if(1){
      R0 = -3.0;
      Rout = 1E4;
    }


    //  Rin=setRin(setihor());
    Rin = 1.0;

    hslope=1.0;

    // define coordinate type
    //    defcoord = 0;
    //    defcoord = 9;
    defcoord = 11;
  }


#endif

  /* output choices */
  tf = 10.0*Rout;

  DTd = 100;			/* dumping frequency, in units of M */
  DTavg = DTd;
  DTener = 2;			/* logfile frequency, in units of M */
  DTi = tf/100.0;			/* image file frequ., in units of M */
  DTdebug = DTd; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;			/* restart file period in steps */

  return(0);
}

int init_global(void)
{

  SAFE=1.3;
  //  cour = 0.9;
  cour=0.8;
  //  cour = 0.5;

  ///////////////////////
  //
  // ENO-RELATED STUFF
  //
  ///////////////////////
  //  avgscheme=WENO5BND;
  avgscheme=DONOR;
  
  do_transverse_flux_integration = 1;
  do_source_integration = 1;
  do_conserved_integration = 1;

  INVERTFROMAVERAGEIFFAILED = 1;
  LIMIT_AC_PRIM_FRAC_CHANGE = 1;
  MAX_AC_PRIM_FRAC_CHANGE = 0.1;


#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  //lim = WENO5BND;
  lim = PARA;
  TIMEORDER=4;
  DOENOFLUX = ENOFINITEVOLUME;
  //  DOENOFLUX = NOENOFLUX;
  fluxmethod=HLLFLUX;
  //  fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  //  UTOPRIMVERSION=UTOPRIM5D1;

#elif(EOMTYPE==EOMFFDE)
  // PARA and TO=4 and HLL not trustable in FFDE so far
  lim = MC;
  TIMEORDER=2;
  fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM2DFINAL;
  // whether/which ENO used to interpolate fluxes
  //DOENOFLUX = ENOFINITEVOLUME;
  DOENOFLUX= NOENOFLUX;
  //DOENOFLUX=ENOFLUXRECON;
#endif



  ranc(7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=0;

  //  BCtype[X1UP]=FIXEDOUTFLOW;
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;
  BCtype[X3UP]=PERIODIC;
  BCtype[X3DN]=PERIODIC;


  rescaletype=1;  //atch: determines the way floors are set in fixup.c:set_density_floors()


  if(BCtype[X1UP]==FIXEDOUTFLOW){ // then doing bondi inflow
    // avoids constant floor activation -- trying to be more physical
    prfloorcoef[RHO]=RHOMIN/100.0;
    prfloorcoef[UU]=UUMIN/100.0;
  }
  else{
    prfloorcoef[RHO]=RHOMIN;
    prfloorcoef[UU]=UUMIN;
  }


#if(MCOORD==CYLMINKMETRIC)

  //  BCtype[X1UP]=FIXED;
  //  BCtype[X1DN]=FIXED;
  //  BCtype[X2UP]=FIXED;
  //  BCtype[X2DN]=FIXED;
  //  BCtype[X1UP]=FIXED;
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=POLARAXIS;
  //BCtype[X1DN]=FIXED;
  BCtype[X2UP]=OUTFLOW;
  BCtype[X2DN]=OUTFLOW;


#else


  if(0){
    BCtype[X1UP]=FIXED;
    BCtype[X1DN]=FIXED;
  }
  else{
    BCtype[X1UP]=RAMESHOUTFLOW;
    BCtype[X1DN]=FIXED;
  }

  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;

#endif

  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];


  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}


int init_grid_post_set_grid(void)
{
  return(0);
}

int check_grid_position(void)
{
  extern int get_num_bnd_zones_used(void);
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  if( MCOORD == KSCOORDS || MCOORD == BLCOORDS ) {
    //computed on all cpus
    Rhor=rhor_calc(0);
    Risco=rmso_calc(PROGRADERISCO);
  }
  
  // diagnostic
  // determine nature of inner radial edge (assumes myid==0 is always there)
  if(myid==0){
    i= - get_num_bnd_zones_used();  //actual number of boundary zones required by the scheme
    j=k=0;
    coord(i,j,k, FACE1, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    trifprintf("rmin: %21.15g\n", r);
    if( MCOORD == KSCOORDS || MCOORD == BLCOORDS ) {
      trifprintf("rmin/rh: %21.15g\n", r / Rhor );
      //    trifprintf("rmin/rsing: %21.15g\n", r / (a+SMALL));
      if(r/Rhor<=1.0){
        trifprintf("inner grid is inside horizon\n");
      }
      else{
        trifprintf("inner grid is outside horizon\n");
        myexit(1);  //do not allow grid outside of horizon
      }
    }
    else if( MCOORD == SPCMINKMETRIC ) {
      if( r <= 0.0 ) {
        trifprintf("rmin < 0 inner grid goes through the point singularity\n");
        myexit(1);
      }
      //else {
      //  trifprintf("rmin > 0 inner grid avoids the point singularity\n");
      //}
    }
  }
  
  return(0);

}


int init_primitives(FTYPE p[][N2M][N3M][NPR])
{
  int whichvel, whichcoord;
  int initreturn;
  int i = 0, j = 0, k = 0, l;
  struct of_geom geom;
  FTYPE r,th,X[NDIM],V[NDIM];
  FTYPE (*A)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3];
  int normalize_densities(FTYPE p[][N2M][N3M][NPR]);
  int init_vpot(int l, int i, int j, int k, FTYPE *A);
  int normalize_field(FTYPE p[][N2M][N3M][NPR]);
  int init_dsandvels(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *p);
  int init_atmosphere(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *pr);
  int init_vpot2field(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE pr[][N2M][N3M][NPR]);
  int init_postfield(FTYPE pr[][N2M][N3M][NPR]);
  extern int check_grid_position(void);
  int res;





  res = check_grid_position();
  if(res) return(res);



  ///////////////////////////////////
  //
  // Assign primitive variables
  //
  ///////////////////////////////////
  trifprintf("Assign primitives\n");

  // assume we start in bl coords and convert to KSprim
  FULLLOOP{
    initreturn=init_dsandvels(&whichvel, &whichcoord,i,j,k,p[i][j][k]); // request densities for all computational centers
    if(initreturn>0) return(1);
    else{
      // transform from whichcoord to MCOORD
      if (bl2met2metp2v(whichvel,whichcoord,p[i][j][k], i,j,k) >= 1)
	FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);
    }
  }

  /////////////////////////////
  //
  // normalize density if wanted
  //
  /////////////////////////////// 
  // at this point densities are still standard, so just send "p"
  trifprintf("Normalize densities\n");
  normalize_densities(p);


  /////////////////////////////
  //
  // Define an atmosphere if wanted
  //
  /////////////////////////////// 

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  // normalized atmosphere
  trifprintf("Add atmosphere\n");
  ZLOOP{
    initreturn=init_atmosphere(&whichvel, &whichcoord,i,j,k,p[i][j][k]);
    if(initreturn>0) return(1);
    else{
      // transform from whichcoord to MCOORD
      if (bl2met2metp2v(whichvel, whichcoord,p[i][j][k], i,j,k) >= 1)
	FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);
    }
  }
#endif


  trifprintf("User filtered to FFDE\n");
  // filter to get force-free
  FULLLOOP{
    filterffde(i,j,k,p[i][j][k]);
  }

  // copy over initial solution as analytic solution
  // SET ANALYTIC SOLUTION FROM vector potential-based solution
  // NEEDED FOR BOUND in case uses panalytic
  ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND, -N3BND, N3-1+N3BND) PLOOP(pl){
    panalytic[i][j][k][pl]=p[i][j][k][pl];
  }


  /////////////////////////////
  //
  // Fixup and Bound variables since some primitive quantities may have changed
  // These may be used to define vector potential below
  // Also setup pre_fixup() type quantities
  //
  /////////////////////////////// 
  trifprintf("Fixup and Bound #1\n");

#if(EOMTYPE!=EOMFFDE)
  // assume EOMFFDE doesn't use "density/ie" to set field, so no need to bound, and no field definition is bad for EOMFFDE
#if(FIXUPAFTERINIT)
  if(fixup(STAGEM1,p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 1);
#endif

  if (bound_prim(STAGEM1,p) >= 1)
    FAILSTATEMENT("init.c:init()", "bound_prim()", 1);

  if(pre_fixup(STAGEM1,p)>=1)
    FAILSTATEMENT("init.c:init()", "pre_fixup()", 1);
#endif



  /////////////////////////////
  //
  // Initialize field from vector potential
  //
  /////////////////////////////// 
  A=emf; // dummy memory space not used till computation so safe.


  trifprintf("Initialize field from vector potential\n");
  FULLLOOPP1{
    for(l=1;l<=3;l++) A[l][i][j][k] = 0.;
  }

  FULLLOOPP1{
    // GODMARK: Caution: Possible to use quantity off grid
    // (e.g. density) to define lower corner value of A, which then
    // defines B at center for lower cells.
    // Do have *grid* quantities for everywhre A is.
    for(l=1;l<=3;l++) 
      init_vpot(l,i,j,k,&A[l][i][j][k]); // request vector potential for all computational corners
  }
  trifprintf("Initialize field from vector potential assign\n");

  // pre-bound dump file (no velocity information)
  pdump = p;
  if (dump(9995) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }

  init_vpot2field(A,p);

  normalize_field(p);
  
  pdump = p;

  // pre-bound dump file (with field from vector potential)
  if (dump(9997) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }

  // set velocity information from field from vector potential
  // may want to use the field to adjust the solution
  
  trifprintf("init_postfield\n");
  //    MYFUN(init_postfield(p),"initbase.c:init()","init_postfield()",1);
  init_postfield(p);

  
  pdump = p;

  // pre-bound dump file
  if (dump(9998) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }

  


  return(0);


}




int init_dsandvels(int *whichvel, int*whichcoord, int ii, int jj, int kk, FTYPE *pr)
{
  FTYPE X[NDIM], V[NDIM], r, th;
  struct of_geom geom;

  FTYPE Fcov[NDIM][NDIM];
  FTYPE Mcon[NDIM][NDIM];
  FTYPE Mcov[NDIM][NDIM];
  FTYPE etacov[NDIM],etacon[NDIM];
  FTYPE Ecov[NDIM],Econ[NDIM],Bcov[NDIM],Bcon[NDIM];
  FTYPE  alpha;
  int j,k;

  void Fcov_numerical(FTYPE *X, FTYPE (*Fcov)[NDIM]);
  extern void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  extern void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);
  extern int EBtopr(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern int EBtopr_2(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern int OBtopr(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);

  struct of_state q;
  FTYPE faradaytest[NDIM][NDIM];


  if(EOMTYPE!=EOMFFDE){
    dualfprintf(fail_file,"Are you sure?\n");
    myexit(1);
  }


  coord(ii, jj, kk, CENT, X);
  bl_coord(X, V);
  r = V[1];
  th = V[2];
  get_geometry(ii, jj, kk, CENT, &geom); 

#define B0 1.0

  pr[RHO]=pr[UU]=0.0;
  pr[U1]=pr[U2]=pr[U3]=0.0;
  pr[B2]=pr[B3]=0;
  if(PROBLEMTYPE==0){
    if(th<M_PI*0.5)  pr[B1]=B0*gdet[horizoni][jj][kk][CENT]/(gdet[ii][jj][kk][CENT]);
    else pr[B1]=-B0*gdet[horizoni][jj][kk][CENT]/(gdet[ii][jj][kk][CENT]);
  }
  else if(PROBLEMTYPE==1){
    // Ruben's talk says they set $\dF^{tr} = C\sin{\theta}/\detg$.
    pr[B1]=B0*gdet[horizoni][N2-1-jj][kk][CENT]/(gdet[ii][N2-1-jj][kk][CENT]);
  }
  else if(PROBLEMTYPE==2){
    //    first get F_{\mu\nu}
    Fcov_numerical(X, Fcov);
    //    dualfprintf(fail_file,"Fcov\n");
    //DLOOP(j,k) dualfprintf(fail_file,"Fcov[%d][%d]=%21.15g\n",j,k,Fcov[j][k]);
    //    dualfprintf(fail_file,"%21.15g %21.15g\n",j,k,Fcov[0][3],Fcov[3][0]);

    // lapse
    // define \eta_\alpha
    // assume always has 0 value for space components
    alpha = 1./sqrt(- geom.gcon[0][0]);

    etacov[TT]=-alpha; // any constant will work.
    SLOOPA(j) etacov[j]=0.0; // must be 0

    // shift
    // define \eta^\beta
    raise_vec(etacov,&geom,etacon);
    //    dualfprintf(fail_file,"raise\n");

    //    DLOOPA(j) dualfprintf(fail_file,"etacon[%d]=%21.15g etacov[%d]=%21.15g\n",j,etacon[j],j,etacov[j]);


    // then get E^\alpha and B^\alpha
    DLOOPA(j) Ecov[j]=0.0;
    DLOOP(j,k) Ecov[j]+=etacon[k]*Fcov[j][k];
    raise_vec(Ecov,&geom,Econ);
    //    dualfprintf(fail_file,"Ecov[3]=%2.15g\n",Ecov[3]);
    //DLOOPA(j) dualfprintf(fail_file,"Econ[%d]=%2.15g\n",j,Econ[j]);


    MtoF(3,Fcov,&geom,Mcon);
    //    dualfprintf(fail_file,"MtoF\n");
    //    DLOOP(j,k) dualfprintf(fail_file,"Mcon[%d][%d]=%21.15g\n",j,k,Mcon[j][k]);

    DLOOPA(j) Bcon[j]=0.0;
    DLOOP(j,k) Bcon[j]+=etacov[k]*Mcon[k][j];

    //    DLOOPA(j) dualfprintf(fail_file,"Econ[%d]=%21.15g Bcon[%d]=%21.15g\n",j,Econ[j],j,Bcon[j]);

    EBtopr(Econ,Bcon,&geom,pr);
    //    dualfprintf(fail_file,"EBtopr\n");

#if(0)
    // check where faraday changed
    get_state(pr,&geom,&q);

    faraday_calc(0,q.bcon,q.ucon,&geom,faradaytest);
    //    DLOOP(j,k) dualfprintf(fail_file,"%21.15g  %21.15g\n",faradaytest[j][k],Fcov[j][k]);
    DLOOP(j,k){
      if(fabs(faradaytest[j][k]-Fcov[j][k])>1E-10){
        dualfprintf(fail_file,"1 %d %d : %21.15g  %21.15g\n",ii,jj,faradaytest[j][k],Fcov[j][k]);
      }
    }
    if(fabs(faradaytest[0][3])>1E-10) dualfprintf(fail_file,"1 Fcov=%21.15g faraday=%21.15g\n",Fcov[0][3],faradaytest[0][3]);
#endif

  }
  else if(PROBLEMTYPE==3){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }
  else if(PROBLEMTYPE==4){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }
  else if(PROBLEMTYPE==5){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }
  else if(PROBLEMTYPE==6){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }
  else if(PROBLEMTYPE==7){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }




  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);
}



//#define NTHETAMAX 10000


// assumes normal field in pr
int init_vpot(int l, int ii, int jj, int kk, FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE X[NDIM],r,th;
  FTYPE Xc[NDIM],rc,thc;
  FTYPE V[NDIM],Vc[NDIM];
  struct of_geom geom;
  struct of_geom geomc;
  FTYPE mcov[NDIM],mcon[NDIM],kcov[NDIM],kcon[NDIM];
  FTYPE th2;
  //  static FTYPE Thetavstheta[NTHETAMAX],dThetadthetavstheta[NTHETAMAX],theta[NTHETAMAX],dtheta,myfloati,myTheta,myThetac;
  int i,k;
  static int firsttime=1;
  FTYPE rprime,zprime;
  FTYPE B0para;
  FTYPE rtrans;
  //  int NTHETA;
  FTYPE Ac;
  FTYPE theu;
  void interpfun(int i, FTYPE th2, FTYPE *theta, FTYPE *fun, FTYPE *answer);
  FTYPE R_smooth;
  void set_ramesh_solution(int ii, int jj, int kk, FTYPE *A);




  if(l!=3){// NOT A_\phi -- atch
    *A = 0;
    return(0);
  }



  if(PROBLEMTYPE<=1) return(0); // otherwise setup poloidal components using vector potential


  if(PROBLEMTYPE==2){

    coord(ii, jj, kk, CORN3, X);
    bl_coord(X, V);
    r = V[1];
    th = V[2];
    get_geometry(ii,jj,kk,CORN3,&geom);


    mcon[TT]=0;
    mcon[RR]=0;
    mcon[TH]=0;
    mcon[PH]=1.0;

    kcon[TT]=1.0;
    kcon[RR]=0;
    kcon[TH]=0;
    kcon[PH]=0;

    lower_vec(mcon,&geom,mcov);
    lower_vec(kcon,&geom,kcov);


    // A_\phi
    *A = -B0*(mcov[PH]+2.0*a*kcov[PH]);
  }
  else if((PROBLEMTYPE==3)||(PROBLEMTYPE==5)){


    coord(ii, jj, kk, CORN3, X);
    bl_coord(X, V);
    r = V[1];
    th = V[2];
    get_geometry(ii,jj,kk,CORN3,&geom);


    if(th<=M_PI*0.5){
      th2=th;
    }
    else{
      th2=M_PI-th;
    }

    // BZ paraboloidal
    if(PROBLEMTYPE==3) rprime=r;
    // BZ paraboloidal with monopole near horizon
    else if(PROBLEMTYPE==5) rprime=sqrt(pow(r,2)+pow(2.0*Rhor,2));

    // setup so A=0 on poles
    // B0para setup so B0 is value of aphi on horizon
    B0para=2.0*B0/(Rhor+4.0*log(2.0)-2.0);
    *A = 0.5*B0para*(rprime*(1.0-cos(th2))+2.0*(1.0+cos(th2))*(1.0-log(1.0+cos(th2))) - 4.0*(1.0-log(2.0)) );

  }
  else if((PROBLEMTYPE==4)||(PROBLEMTYPE==6)){
    // GRMHD nearly-paraboloidal

    // already have data

    coord(ii, jj, kk, CORN3, X);
    bl_coord(X, V);
    r = V[1];
    th = V[2];
    get_geometry(ii,jj,kk,CORN3,&geom);


    dtheta=theta[1]-theta[0]; // uniform grid


    if(th<0.0){
      // function is even
      th2=-th;
    }
    else if(th<=M_PI*0.5){
      // theta is itself
      th2=th;
    }
    else if(th<=M_PI){
      th2=M_PI-th;
    }
    else{
      // function is even
      th2=-(M_PI-th);
    }

    // linearly interpolate \Theta using \theta
    myfloati=(th2-0.0)/dtheta;

    myTheta=Thetavstheta[(int)myfloati]
    +(Thetavstheta[(int)myfloati+1]-Thetavstheta[(int)myfloati])/(theta[(int)myfloati+1]-theta[(int)myfloati])*(th2-theta[(int)myfloati]);

    rtrans=Rhor;

    // n=1/4 strict
    if(PROBLEMTYPE==4) rprime=r;
    // n=1/4 but relaxes to n=1 near horizon
    else if(PROBLEMTYPE==6) rprime=sqrt(pow(r,2)+pow(rtrans,2));

    *A = 0.5*B0*pow(rprime,-INDEXN+1.0)*myTheta;

  }



  else if(PROBLEMTYPE==7){

    set_ramesh_solution(ii, jj, kk, A);


  }



  //  dualfprintf(fail_file,"i=%d j=%d\n",ii,jj);
  //  PLOOP(pl) dualfprintf(fail_file,"%21.15g ",p[ii][jj][kk][pl]);
  //  dualfprintf(fail_file,"\n");


  return(0);

}





void get_nearlypara_data(void)
{
  int i;
  FILE * inTheta;
  int dumi;



  // for PROBLEMTYPE==4
  NTHETA=10000;
  INDEXN=(1.0/4.0); // requires Theta to be l=-n




  if(myid==0){
    if( (inTheta=fopen("myThetavstheta.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open myThetavstheta.txt\n");
      myexit(100);
    }

    // read in data file
    while(fgetc(inTheta)!='\n'); // skip first line, a comment
    for(i=0;i<NTHETA;i++){
      fscanf(inTheta,"%d %lf %lf %lf\n",&dumi,&theta[i],&Thetavstheta[i],&dThetadthetavstheta[i]);  //SASMARK:  actually, theta[i] is not \theta but u = Z/R = tan(\theta)
    }

  }

  // send data to all CPUs
#if(USEMPI)
  MPI_Bcast(&NTHETA,1,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&theta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Thetavstheta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&dThetadthetavstheta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
#endif

}





void get_ramesh_data(void)
{
  int i;
  FILE * inTheta;
  FTYPE ftemp,ftemp1,ftemp2;


  ////////////////////////
  // remove empty line at top of file and convert D to E
  // and set NTHETA
  //
  // CHOOSE which Ramesh problem
  //
  ///////////////////////////////

  if(myid==0){
#if(RAMESHTYPE==0)
    // for nu.75_m.5.txt
    Ttpow = 0; // not known right now
    NTHETA=1314;
    if( (inTheta=fopen("nu.75_m.5.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open nu.75_m.5.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==1)
    // for nu1.0_m.25.txt
    // T(\theta)\propto \theta^{Ttpow} -> \theta_j \propto r^{-1/(1+Ttpow/nu)}
    Ttpow = 1.0;
    NTHETA=2030;
    if( (inTheta=fopen("nu1.0_m.25.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open nu1.0_m.25.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==125)
    // or nu1.0_m.25_hres.txt
    Ttpow = 1.0;
    NTHETA=2030;
    if( (inTheta=fopen("nu1.0_m.25_hres.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open nu1.0_m.25_hres.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==2)
    // for nu.75_m.0005.txt
    Ttpow = 0.0; // not known right now
    NTHETA=1206;
    if( (inTheta=fopen("nu.75_m.0005.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open nu.75_m.0005.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==3)
    // for ramesh_nu1_m.25_s1.txt
    Ttpow = 1.0;
    NTHETA=2030;
    if( (inTheta=fopen("ramesh_nu1_m.25_s1.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu1_m.25_s1.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==4)
    // for ramesh_nu.75_m.1_s.5.txt
    Ttpow = 3.0;
    NTHETA=1505;
    if( (inTheta=fopen("ramesh_nu.75_m.1_s.5.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu.75_m.1_s.5.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==5)
    // for ramesh_nu.75_m.1_s.28146.txt
    Ttpow = 1.361;
    NTHETA=2032;
    if( (inTheta=fopen("ramesh_nu.75_m.1_s.28146.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu.75_m.1_s.28146.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==6)
    // for ramesh_nu1.25_m.4_soM2.5.txt
    Ttpow = 0.0; // special solution
    NTHETA=1502;
    if( (inTheta=fopen("ramesh_nu1.25_m.4_soM2.5.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu1.25_m.4_soM2.5.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==7)
    // for ramesh_nu1.25_m.4_soM1.6.txt
    Ttpow = 2.0-nu; // for \nu>1
    NTHETA=1865;
    if( (inTheta=fopen("ramesh_nu1.25_m.4_soM1.6.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu1.25_m.4_soM1.6.txt\n");
      myexit(100);
    }
#endif	

    if(NTHETA>=NTHETAMAX){
      dualfprintf(fail_file,"NTHETA=%d and NTHETAMAX=%d\n",NTHETA,NTHETAMAX);
      myexit(7);
    }


    ///////////////////////////////////////////
    //
    // READ IN DATA FILE
    //
    ////////////////////////////////////////////

    // skip first blank line
    //while(fgetc(inTheta)!='\n'); // skip first line, a comment

    // read in data file

    // read in first row (parameters)
    fscanf(inTheta,"%lf %lf %lf %lf",&nu,&mm,&ss,&ucrit);
    //dualfprintf(fail_file,"got0: %g %g %g %g\n",nu,mm,ss,ucrit);

    // NTHETA does not include first line



    for(i=0;i<NTHETA;i++){
      //	while(!feof(inTheta)){
      fscanf(inTheta,"%lf %lf %lf",&ftemp,&ftemp1,&ftemp2);
      // Z/R itself
      theta[i]=ftemp; // ftemp=u= Z/R=tan(\theta) //SASMARK:  actually, theta[i] is not \theta but u = Z/R = tan(\theta)
      Thetavstheta[i]=ftemp1;
      dThetadthetavstheta[i]=ftemp2;
#if(0)
      // theta
      theta[NTHETA-1-i]=atan(1.0/ftemp); // true theta since file has ftemp=u= Z/R=tan(\theta)
      Thetavstheta[NTHETA-1-i]=ftemp1;
      dThetadthetavstheta[NTHETA-1-i]=ftemp2;
#endif
      //	  dualfprintf(fail_file,"got1: %d %g %g %g\n",i,theta[NTHETA-1-i],Thetavstheta[NTHETA-1-i],dThetadthetavstheta[NTHETA-1-i]);
      //	  i++;
    }
    //	NTHETA=i;

    // \theta_j \propto r^{-jetalpha}
    jetalpha = 1.0/(1.0+(Ttpow/nu));

    trifprintf("Ramesh parameters: %g %g %g %g %g %g\n",nu,mm,ss,ucrit,Ttpow,jetalpha);


  }





  // send data to all CPUs
#if(USEMPI)
  MPI_Bcast(&nu,1,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&mm,1,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&ss,1,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&ucrit,1,MPI_FTYPE,0,MPI_COMM_WORLD);

  MPI_Bcast(&Ttpow,1,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&jetalpha,1,MPI_FTYPE,0,MPI_COMM_WORLD);


  MPI_Bcast(&NTHETA,1,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&theta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Thetavstheta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&dThetadthetavstheta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
#endif


}






////////////////////////////////////////////////////
//
// Ramesh constant velocity disk
//
///////////////////////////////////////////////////
void set_ramesh_solution(int ii, int jj, int kk, FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE X[NDIM],V[NDIM],r,th;
  FTYPE Xc[NDIM],Vc[NDIM],rc,thc;
  struct of_geom geom;
  struct of_geom geomc;
  FTYPE th2;
  FILE * inTheta;
  //  static FTYPE Thetavstheta[NTHETAMAX],dThetadthetavstheta[NTHETAMAX],theta[NTHETAMAX],dtheta,myfloati,myTheta,myThetac;
  int dumi;
  int i,k;
  static int firsttime=1;
  FTYPE rprime,zprime;
  FTYPE B0para;
  FTYPE rtrans;
  //  int NTHETA;
  FTYPE Ac;
  FTYPE theu;
  void interpfun(int i, FTYPE th2, FTYPE *theta, FTYPE *fun, FTYPE *answer);
  FTYPE R_smooth;



  ///////////////////////////////////////////
  //
  // Already have data, now interpolate data to HARM grid
  //
  ////////////////////////////////////////////


  coord(ii, jj, kk, CORN3, X);
  bl_coord(X, V);
  r = V[1];
  th = V[2];
  get_geometry(ii,jj,kk,CORN3,&geom);


  coord(ii, jj, kk, CENT, Xc);
  bl_coord(Xc, Vc);
  rc = Vc[1];
  thc = Vc[2];
  get_geometry(ii,jj,kk,CENT,&geomc);




#if(MCOORD==CYLMINKMETRIC)

  // choose smoothing cylindrical radius
  // only used for u=z/R
  // u = (|z|+R_smooth)/R
  R_smooth = 1.0;

  // u = (|z|+R_smooth)/R
  theu=th2=( fabs(th) + R_smooth )/fabs(r); // Z/R really, but solution symmetric around Z=0
  // R
  rprime=fabs(r);
#else
  // assume hidding singularity at r=0 with star, so no smoothing needed
  R_smooth=0.0;

  if(th<0.0){
    // function is even
    th2=-th;
  }
  else if(th<=M_PI*0.5){
    // theta is itself
    th2=th;
  }
  else if(th<=M_PI){
    th2=M_PI-th;
  }
  else{
    // function is even
    th2=-(M_PI-th);
  }
  // R
  rprime=fabs(r)*fabs(sin(th2));
  // z
  zprime=fabs(r)*fabs(cos(th2));
  //    theu = ( fabs(zprime)+R_smooth )/fabs(rprime);
  theu = fabs(cos(th2)/sin(th2));
#endif


  // \theta is not uniform, so loop to find first theta
  for(i=0;i<NTHETA;i++){
    if(theta[i]>theu) break;   //SASMARK:  actually, theta[i] is not \theta but u = Z/R = tan(\theta)
  }

  if(i==NTHETA) i=NTHETA-1;
  

  interpfun(i, theu, theta, Thetavstheta, &myTheta);

  // A_\phi = P(R,Z) = R^\nu T(u)
  *A = B0*pow(rprime,nu)*myTheta;

  if(MCOORD!=CYLMINKMETRIC){
    if(th<0.0 || th>M_PI ){
      *A *=-1.0;   //SASMARK:  B^r = dA_\phi/d\theta has to be symmetric, so flipping the sign is correct here
    }
  }

  // for testing purposes (dump999x), place omegafanalytic into internal energy.  This is set back to 0 before evolution.
  if( ii < N1 + N1BND && jj < N2 + N2BND && kk < N3 + N3BND )  //SASMARKx: to avoid getting out of bounds
    p[ii][jj][kk][RHO]=*A;


  //////////////////////////////////////////////////////////
  //
  /////////////////////////////// CENTERED QUANTITIES (only applies for diagonal metric)
  //
  //////////////////////////////////////////////////////////

  //////////////////////////////
  //
  // set analytic B^\phi.  With this choice and v^i chosen by stationarity, the analytic solution is set so FIXED boundaries are correct
  //
  ///////////////////////////////
#if(MCOORD==CYLMINKMETRIC)
  // u = (|Z|+R_smooth)/R
  // u=Z/R really, but solution symmetric around Z=0 and symmetric around R
  theu=th2=(fabs(thc)+R_smooth)/fabs(rc);
  // R
  rprime=fabs(rc);
#else
  if(thc<0.0){
    // function is even
    th2=-thc;
  }
  else if(thc<=M_PI*0.5){
    // theta is itself
    th2=thc;
  }
  else if(thc<=M_PI){
    th2=M_PI-thc;
  }
  else{
    // function is even
    th2=-(M_PI-thc);
  }
  // R
  // bit different than how setting vector potential, but maybe A_\phi is even and B^\phi not, across polar axes?
  rprime=fabs(rc)*fabs(sin(thc)); // R
  zprime=fabs(rc)*fabs(cos(thc)); // z
  //    theu=fabs( fabs(r*cos(thc)) + R_smooth )/fabs(r*sin(thc)) ; // u=(|Z|+R_smooth)/R
  theu = (zprime+R_smooth)/rprime;
#endif

  // linearly interpolate \Theta using \theta
  // \theta is not uniform, so loop to find first theta

  for(i=0;i<NTHETA;i++){
    if(theta[i]>theu) break;  //SASMARK:  actually, theta[i] is not \theta but u = Z/R = tan(\theta)
  }

  if(i==NTHETA) i=NTHETA-1;

  interpfun(i, theu, theta, Thetavstheta, &myThetac);
  interpfun(i, theu, theta, dThetadthetavstheta, &mydThetac);



  // GODMARK: in this code, A is over same range as p, so ok to do this, else need check on i,j,k


  // in orthonormal basis
  if(R_smooth==0.0 && (MCOORD==CYLMINKMETRIC)){
    p[ii][jj][kk][B3]=B0*(nu*ss*pow(rprime,nu-2.0)*pow(myThetac,(nu-1.0)/nu));
    // transform to coordinate basis for diagonal metric (at least \phi not connected to any other coordinate -- always true so far in HARM)
    p[ii][jj][kk][B3]/=sqrt(geomc.gcov[PH][PH]);

    // just do directly (not general)
    // original Ramesh formula is super divergent!?
    //    p[ii][jj][B3]=B0*(-nu*ss*pow(rprime,nu-3.0)*pow(myTheta,(nu-1.0)/nu));



    // B3
    // below not general
    //    p[ii][jj][B3]=B0*(-nu*ss*pow(rprime,nu-3.0)*pow(myThetac,(nu-1.0)/nu));

    if(p[ii][jj][kk][B3]>1E30) p[ii][jj][kk][B3]=1E30;
    else if(p[ii][jj][kk][B3]<-1E30) p[ii][jj][kk][B3]=-1E30;

  }
  else{
    //      p[ii][jj][B3]=B0*(nu*ss*pow(rprime+R_smooth,nu-2.0)*pow(myThetac,(nu-1.0)/nu));
    // transform to coordinate basis for diagonal metric (at least \phi not connected to any other coordinate -- always true so far in HARM)
    //      p[ii][jj][B3]/=sqrt(geomc.gcov[PH][PH]);
    p[ii][jj][kk][B3]=0.0;
  }


  if(MCOORD==CYLMINKMETRIC){
    // B1
    // orthonormal : B^{\hat{R}} = dA_\phi/d\theta 1/\detg \sqrt{g_{RR}} = R^{\nu-1} dT/d\theta = R^{\nu-2} dT/du
    // where for CYL mydThetac = dT/du
    p[ii][jj][kk][B1]=B0*(pow(rprime,nu-2.0)*mydThetac);
    // coordinate
    p[ii][jj][kk][B1]/=sqrt(geomc.gcov[1][1]);  //SASMARKx
    if(p[ii][jj][kk][B1]>1E30) p[ii][jj][kk][B1]=1E30;
    else if(p[ii][jj][kk][B1]<-1E30) p[ii][jj][kk][B1]=-1E30;

    // B2
    // 
    //    p[ii][jj][B2]=B0*pow(rprime,nu-2.0)*(nu*myThetac  - theu*mydThetac);
    p[ii][jj][kk][B2]=-B0*pow(rprime,nu-2.0)*(nu*myThetac  - theu*mydThetac);
    p[ii][jj][kk][B2]/=sqrt(geomc.gcov[2][2]);  //SASMARKx
    if(p[ii][jj][kk][B2]>1E30) p[ii][jj][kk][B2]=1E30;
    else if(p[ii][jj][kk][B2]<-1E30) p[ii][jj][kk][B2]=-1E30;

#if(MCOORD==CYLMINKMETRIC)
    // B^\phi is antisymmetric across equator
    // thc=Z.  Normal range is Z>0
    if(thc<0.0){
      p[ii][jj][kk][B3]*=-1.0;
      p[ii][jj][kk][B1]*=-1.0;
    }
    if(rc<0.0){
      p[ii][jj][kk][B1]*=-1.0;
    }
#else
    p[ii][jj][kk][B1]*=-1.0;

    // B^\phi is antisymmetric across equator
    if(thc>M_PI*0.5){
      p[ii][jj][kk][B3]*=-1.0;
      p[ii][jj][kk][B1]*=-1.0;
    }
    if(rc<0.0){
      p[ii][jj][kk][B1]*=-1.0;
    }
#endif

  }
  else{
    p[ii][jj][kk][B1]=0.0;
    p[ii][jj][kk][B2]=0.0;
  }





  ////////////////
  //
  // need omegaf analytic! (only used at center)
  //
  ////////////////


  Ac=B0*pow(rprime,nu)*myThetac;

  ////////////////////
  //
  // Here, R_smooth is used for \Omega as well!
  //
#if(MCOORD==CYLMINKMETRIC)
  if(R_smooth==0.0){
    omegafanalytic[ii][jj][kk][0]=mm*pow(Ac/Thetavstheta[0]+R_smooth,-1.0/nu);
    //omegafanalytic[ii][jj][kk][0]=mm*pow(Ac/Thetavstheta[0],-1.0/nu);
    // test below
    //omegafanalytic[ii][jj][kk][0]=mm*pow(B0*pow(rprime,nu)*1E-8,-1.0/nu);
  }
  else{
    //      dualfprintf(fail_file,"got here assigning smoothed omega_f: %g %g\n",thc,5.0*R_smooth);
    // GODMARK: assumes vertical direction is as resolved on the R_smooth scale as the cylindrical radial direction
    // otherwise the disk fixed \Omega_F will be contaminated too much by the smoothing
    if(fabs(thc)<2.0*R_smooth){
      omegafanalytic[ii][jj][kk][0]=mm/(fabs(rc)+R_smooth);
    }
    else{
      omegafanalytic[ii][jj][kk][0]=mm*pow(Ac/Thetavstheta[0]+R_smooth,-1.0/nu)/(fabs(thc)+R_smooth);
    }

  }
#else
  if(rc<=Rin){
    // fixed \omega inside star (surface)
    omegafanalytic[ii][jj][kk][0]=mm*pow(B0*pow(Rin,nu),-1.0/nu); // equatorial value for all of star
  }
//  else if(fabs(thc-M_PI*0.5)<0.2){ // GODMARK arbitrary 
  else if( (totalsize[2]/2-(startpos[2]+jj+1)<2)&&((startpos[2]+jj)-totalsize[2]/2<2)){ // GODMARK arbitrary 
    // within disk but not in star
    omegafanalytic[ii][jj][kk][0]=mm*pow(Ac/Thetavstheta[0],-1.0/nu);
  }
  else omegafanalytic[ii][jj][kk][0]=0.0;
#endif


  // TRYING boundary condition with \omega antisymemtric across cyl R=0
#if(MCOORD==CYLMINKMETRIC)
  if(rc<0.0){
    omegafanalytic[ii][jj][kk][0]*=-1.0;
  }
#else
  if(rc<0.0){
    omegafanalytic[ii][jj][kk][0]*=-1.0;
  }

  if(thc<0.0 || thc>M_PI ){ // assuming consistently using boundary condition where u^\phi is antisymmetric
    omegafanalytic[ii][jj][kk][0]*=-1.0;
  }

#endif


  // for testing purposes (dump999x), place omegafanalytic into internal energy.  This is set back to 0 before evolution.
  p[ii][jj][kk][UU]=omegafanalytic[ii][jj][kk][0];

  //     dualfprintf(fail_file,"%d %d : th2=%21.15g : theta[%d]=%21.15g myTheta=%21.15g A=%21.15g\n",ii,jj,th2,i,theta[i],myTheta,*A);




}






#define LINEARTYPE 0
#define LOGTYPE 1
#define QUADRATICTYPE 2

#define INTERPTYPE LOGTYPE
//#define INTERPTYPE LINEARTYPE
//#define INTERPTYPE QUADRATICTYPE

void interpfun(int i, FTYPE th2, FTYPE *theta, FTYPE *fun, FTYPE *answer)
{
  FTYPE slope,intercept;
  FTYPE slope1,slope2,xminusx0;
  FTYPE f0,f1,f2,x0,x1,x2;



  if(INTERPTYPE==LINEARTYPE){
    // linearly interpolate \Theta using \theta
    //      *answer = fun[i-1] + (fun[i]-fun[i-1])/(theta[i]-theta[i-1])*(th2-theta[i-1]);    
    slope = (fun[i]-fun[i-1])/(theta[i]-theta[i-1]);
    intercept = fun[i-1];
    *answer = slope*(th2-theta[i-1]) + intercept;

  }
  else if(INTERPTYPE==QUADRATICTYPE){
    // quadratically interpolate \Theta using \theta
    if(i-1<0){
      f0=fun[i];
      f1=fun[i+1];
      f2=fun[i+2];
      x0=theta[i];
      x1=theta[i+1];
      x2=theta[i+2];	
    }
    else if(i+1>=NTHETA){
      f0=fun[i-2];
      f1=fun[i-1];
      f2=fun[i];
      x0=theta[i-2];
      x1=theta[i-1];
      x2=theta[i];	
    }
    else{
      f0=fun[i-1];
      f1=fun[i];
      f2=fun[i+1];
      x0=theta[i-1];
      x1=theta[i];
      x2=theta[i+1];
    }

    slope2 = ((f0-f2)/(x0-x2) - (f2-f1)/(x2-x1))/(x0-x1);
    slope1 = (f0-f1)/(x0-x1) + (f0-f2)/(x0-x2) - (f2-f1)/(x2-x1);
    xminusx0 = (th2-x0);

    *answer = slope2*pow(xminusx0,2.0) + slope1*xminusx0 + f0;
  }
  else if(INTERPTYPE==LOGTYPE){
    // log interpolate \Theta using \theta
    slope = log(fun[i]/fun[i-1])/log(theta[i]/theta[i-1]);
    if(fabs(slope)<1E-10 || !isfinite(slope)) *answer=fun[0];  //SASMARKx
    else if(fun[i-1]<0.0){
      // assume bi-log
      *answer=-exp( slope*log(th2/theta[i-1])+log(-fun[i-1]) );
    }
    else *answer=exp( slope*log(th2/theta[i-1])+log(fun[i-1]) );

    //dualfprintf(fail_file,"ii=%d jj=%d slope=%g myTheta=%g\n",ii,jj,slope,myTheta);
  }



}




int init_vpot2field(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE pr[][N2M][N3M][NPR])
{
  extern int vpot2field(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE p[][N2M][N3M][NPR]);

#if( PROBLEMTYPE<=1 )
  //#if( (PROBLEMTYPE<=1)||(PROBLEMTYPE==7) ) // GODMARK : test of errors
  return(0);
  // didn't need vector potential formulation
#else
  return(vpot2field(A,pr));
#endif
}



int init_postfield(FTYPE pr[][N2M][N3M][NPR])
{
  FTYPE Bccon[NDIM];
  FTYPE prnew[NPR];
  int i,j,k;
  extern int OBtopr(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  struct of_geom geom;
  FTYPE r,th,X[NDIM], V[NDIM];
  FTYPE get_omegastar(struct of_geom *geom, FTYPE r, FTYPE th);


  ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND, -N3BND, N3-1+N3BND){

    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r = V[1];
    th = V[2];
    get_geometry(i, j, k, CENT, &geom);


    Bccon[0]=0;
    Bccon[1]=pr[i][j][k][B1];
    Bccon[2]=pr[i][j][k][B2];
    Bccon[3]=pr[i][j][k][B3];

    Omegastar=get_omegastar(&geom,r,th);

    if(OBtopr(Omegastar,Bccon,&geom,prnew)>=1){
      dualfprintf(fail_file, "OBtopr(vcon2pr): space-like error in init_postfield()\n");
      dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
      failed=1;
      return(1);

    }

    pr[i][j][k][U1]=prnew[U1];
    pr[i][j][k][U2]=prnew[U2];
    pr[i][j][k][U3]=prnew[U3];

  }

  return(0);

}

FTYPE get_omegastar(struct of_geom *geom, FTYPE r, FTYPE th)
{
  FTYPE omegak,omegaf,omegah,omegadiskbh;
  FTYPE ftemp1,ftemp2;
  FTYPE rtrans,ntrans;
  FTYPE CON,CON2;
  FTYPE omegatrans;
  FTYPE vtrans;



#if(PROBLEMTYPE<7)

  omegak=1.0/(pow(r,3.0/2.0)+a);
  omegah=(a/(2.0*Rhor));

  //  omegadiskbh=0.6*omegah;
  // does quite well till about dump 11, then has problems and jump to above equator is omegaf2/omegah~0.27.  Leads to crispy omegaf by dump=20

  //  omegadiskbh=0.2652*omegah; // goes uu0~30 before dump 1
  // above chooses sharp changed to omegaf~1 just outside equator by 2nd dump
  // by dump=8 or 11, polar region has settled to standard form with 0.65*omegah just above equator

  //  omegadiskbh=1.0*omegah; // lasts till dump 8 before uu0~30 , but close to disk omegaf~0.27*omegah after settles to non-force-free

  omegadiskbh=0.4*omegah;

#if(0)
  // Just Keplerian for all radii
  omegaf=omegak;
#elif(0)
  if(r>Risco){
    omegaf=omegak;
  }
  else{
    omegaf=1.0/(pow(Risco,3.0/2.0)+a);
  }
#elif(0)
  omegaf = omegak/(1.0 + pow(Rhor/r,2.0));
#elif(0)
  // goes from \Omega_K at large radii to (omegadiskbh) on horizon
  if(r>2.0*Rhor){
    ftemp1=omegak*(1.0-pow(2.0*Rhor/r,3.0));
    ftemp2=(omegadiskbh)*pow(2.0*Rhor/r,3.0);
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh;
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see omegaf_vs_omegah_bz77_grffe_grmhd.nb
  // math1
  rtrans=Rhor;
  ntrans=6.21;
  omegadiskbh=(0.5/0.85)*omegah;

  if(r>rtrans){
    ftemp1=omegak*(1.0-pow(rtrans/r,3.0));
    ftemp2=(omegadiskbh)*pow(rtrans/r,ntrans);
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh;
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see omegaf_vs_omegah_bz77_grffe_grmhd.nb
  // math2
  rtrans=Rhor;
  ntrans=10.744;
  omegadiskbh=(0.34)*omegah;

  if(r>rtrans){
    ftemp1=omegak*(1.0-pow(rtrans/r,3.0));
    ftemp2=(omegadiskbh)*pow(rtrans/r,ntrans);
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh;
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see omegaf_vs_omegah_bz77_grffe_grmhd.nb
  // math3
  rtrans=2.0*Rhor;
  //  omegadiskbh=(0.34)*omegah; // works ok
  //  omegadiskbh=(0.315)*omegah;
  omegadiskbh=(0.27)*omegah;
  //  omegadiskbh=(0.4)*omegah;
  //  omegadiskbh=(0.32)*omegah;
  ntrans=3.0/(omegadiskbh*(a+pow(rtrans,3.0/2.0)));

  if(r>rtrans){
    ftemp1=omegak*(1.0-pow(rtrans/r,3.0));
    ftemp2=(omegadiskbh)*pow(rtrans/r,ntrans);
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh;
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see omegaf_vs_omegah_bz77_grffe_grmhd.nb
  // math4
  rtrans=2.0*Rhor;
  omegadiskbh=(0.5/0.85)*omegah;
  CON=-0.0479309;
  CON2=0.0612062;

  if(r>rtrans){
    ftemp1=omegak;
    ftemp2=0;
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh+CON*pow(r-Rhor,2.0)+CON2*(r-Rhor);
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see para_bz_omegaf_compute_plot.nb
  omegadiskbh=0.27*omegah;
  rtrans=2.0*Rhor;
  CON=-3.0*omegadiskbh*pow(Rhor - rtrans,-2); // B
  CON2=2.0*omegadiskbh*pow(-Rhor + rtrans,-3); // C

  if(r>rtrans) omegaf=0;
  else if(r<Rhor) omegaf=omegadiskbh;
  else omegaf=omegadiskbh + CON*(r-Rhor)*(r-Rhor) + CON2*(r-Rhor)*(r-Rhor)*(r-Rhor);

#elif(1)

  // see para_bz_omegaf_compute_plot.nb
  omegadiskbh=0.27*omegah;
  rtrans=2.0*Rhor;
  vtrans=0.5;
  //  omegatrans=(vtrans/rtrans)/(sqrt(geom->gcov[3][3])/rtrans);
  omegatrans=(vtrans/rtrans)/(r/rtrans);

  // fsB
  CON = (-3.0*omegadiskbh*rtrans*rtrans - (Rhor - 4.0*rtrans)*vtrans)/( pow(rtrans-Rhor,2)*rtrans*rtrans );

  // fsC
  CON2 =(+2.0*omegadiskbh*rtrans*rtrans + (Rhor - 3.0*rtrans)*vtrans)/( pow(rtrans-Rhor,3)*rtrans*rtrans );


  if(r>=rtrans) omegaf=omegatrans;
  else if(r<=Rhor) omegaf=omegadiskbh;
  else omegaf=omegadiskbh + CON*(r-Rhor)*(r-Rhor) + CON2*(r-Rhor)*(r-Rhor)*(r-Rhor);

#endif


#elif(PROBLEMTYPE==7)

  // mm is global variable with local scope to this file

  omegaf=omegafanalytic[geom->i][geom->j][geom->k][0];
  // only true at equator
  //#if(MCOORD==CYLMINKMETRIC)
  //  omegaf=mm/r;
  //#else
  //  omegaf=mm/(r*sin(th));
  //#endif

  //  dualfprintf(fail_file,"omegaf=%21.15g mm=%21.15g\n",omegaf,mm);

#endif


  return(omegaf);



  // NON-ROTATING DISK
  //  return(0.0);

}

// assumes we are fed the true densities
int normalize_densities(FTYPE p[][N2M][N3M][NPR])
{

  return(0);

}


// assumes normal field definition
int normalize_field(FTYPE p[][N2M][N3M][NPR])
{
  int i,j,k;
  FTYPE bsq_ij;
  SFTYPE bsq_max, norm, beta_act;
  struct of_geom geom;
  FTYPE X[NDIM], V[NDIM];
  FTYPE r,th;

  if(PROBLEMTYPE==7) return(0); // don't worry about normalization

  bsq_max = 0.;
  ZLOOP {
    get_geometry(i, j, k, CENT, &geom);    
    coord(i, j, k, CENT, X);
    bl_coord(X, V);

    r = V[1];
    th = V[2];

    if((r>Rhor)&&(fabs(th-M_PI*0.5)<0.1)){
      if (bsq_calc(p[i][j][k], &geom, &bsq_ij) >= 1)
        FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);

      if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    }
  }

  mpimax(&bsq_max);
  trifprintf("initial bsq_max: %21.15g\n", bsq_max);

  /* finally, normalize to set field strength */
  norm = sqrt((B0*B0)/bsq_max);

  bsq_max = 0.;
  ZLOOP {
    p[i][j][k][B1] *= norm;
    p[i][j][k][B2] *= norm;
    p[i][j][k][B3] *= norm;

    get_geometry(i, j, k, CENT, &geom);
    coord(i, j, k, CENT, X);
    bl_coord(X, V);

    r = V[1];
    th = V[2];

    if((r>Rhor)&&(fabs(th-M_PI*0.5)<0.1)){
      if (bsq_calc(p[i][j][k], &geom, &bsq_ij) >= 1)
        FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);

      if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    }

  }
  mpimax(&bsq_max);
  trifprintf("new initial bsq_max: %21.15g\n", bsq_max);

  trifprintf("new bsq_max: %21.15g\n", bsq_max);




  return(0);
}

/* this is a little test calculation with a radial field, designed to 
make the scheme fail */

/*     Br0 = 1.0 ; ZLOOP { GSET(i,j,k,CENT) p[i][j][k][B1] = Br0/(rcurr*rcurr) 
; p[i][j][k][B2] = 0. ; p[i][j][k][B3] = 0. ; } */


#undef SLOWFAC



#define GAMMIEDERIVATIVE 0
#define NUMREC 1

//#define DXDERTYPE NUMREC 
#define FCOVDERTYPE GAMMIEDERIVATIVE

// see conn_func() for notes
#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
#define DXDELTA 1E-5
#elif(REALTYPE==LONGDOUBLETYPE)
#define DXDELTA 1E-6
#endif

void Fcov_numerical(FTYPE *X, FTYPE (*Fcov)[NDIM])
{
  int j,k,l;
  FTYPE Xhk[NDIM], Xlk[NDIM];
  FTYPE Xhj[NDIM], Xlj[NDIM];
  FTYPE mcovhj,mcovlj,kcovhj,kcovlj;
  FTYPE mcovhk,mcovlk,kcovhk,kcovlk;
  FTYPE mcov_func_mcoord(FTYPE* X, int i, int j); // i not used
  FTYPE kcov_func_mcoord(FTYPE* X, int i, int j); // i not used
  extern FTYPE dfridr(FTYPE (*func)(FTYPE*,int,int), FTYPE *X,int ii, int jj, int kk);

  if(FCOVDERTYPE==GAMMIEDERIVATIVE){

    for(k=0;k<NDIM;k++){
      for(j=0;j<NDIM;j++){

        for(l=0;l<NDIM;l++) Xlk[l]=Xhk[l]=Xlj[l]=Xhj[l]=X[l]; // location of derivative
        Xhk[k]+=DXDELTA; // shift up
        Xlk[k]-=DXDELTA; // shift down

        Xhj[j]+=DXDELTA; // shift up
        Xlj[j]-=DXDELTA; // shift down

        //	  dualfprintf(fail_file,"got here1: k=%d j=%d\n",k,j);


        mcovhj=mcov_func_mcoord(Xhk,0,j); // 0 not used
        //	  dualfprintf(fail_file,"got here1.1: k=%d j=%d\n",k,j);
        mcovlj=mcov_func_mcoord(Xlk,0,j); // 0 not used
        //	  dualfprintf(fail_file,"got here1.2: k=%d j=%d\n",k,j);
        mcovhk=mcov_func_mcoord(Xhj,0,k); // 0 not used
        //	  dualfprintf(fail_file,"got here1.3: k=%d j=%d\n",k,j);
        mcovlk=mcov_func_mcoord(Xlj,0,k); // 0 not used
        //	  dualfprintf(fail_file,"got here1.4: k=%d j=%d\n",k,j);

        kcovhj=kcov_func_mcoord(Xhk,0,j); // 0 not used
        //	  dualfprintf(fail_file,"got here1.5: k=%d j=%d\n",k,j);
        kcovlj=kcov_func_mcoord(Xlk,0,j); // 0 not used
        //	  dualfprintf(fail_file,"got here1.6: k=%d j=%d\n",k,j);
        kcovhk=kcov_func_mcoord(Xhj,0,k); // 0 not used
        //	  dualfprintf(fail_file,"got here1.7: k=%d j=%d\n",k,j);
        kcovlk=kcov_func_mcoord(Xlj,0,k); // 0 not used
        //	  dualfprintf(fail_file,"got here1.8: k=%d j=%d\n",k,j);

        //	  dualfprintf(fail_file,"got here2\n");

        Fcov[j][k] = B0*(
          +(mcovhj - mcovlj) / (Xhk[k] - Xlk[k])
          -(mcovhk - mcovlk) / (Xhj[j] - Xlj[j])
          +2.0*a*(
          +(kcovhj - kcovlj) / (Xhk[k] - Xlk[k])
          -(kcovhk - kcovlk) / (Xhj[j] - Xlj[j])
          )
          );
      }// j
    }// k
  }
  else if(FCOVDERTYPE==NUMREC){

    for(k=0;k<NDIM;k++) for(j=0;j<NDIM;j++){
      // 0 in dfridr not used
      Fcov[j][k]=B0*(
        +dfridr(mcov_func_mcoord,X,0,j,k)
        -dfridr(mcov_func_mcoord,X,0,k,j)
        +2.0*a*(+dfridr(kcov_func_mcoord,X,0,j,k)
        -dfridr(kcov_func_mcoord,X,0,k,j)
        )
        );
    }

  }
}

#undef GAMMIEDERIVATIVE
#undef NUMREC
#undef FCOVDERTYPE
#undef DXDELTA




// returns MCOORD m_\mu form of m^\mu={0,0,0,1} value for jth element
FTYPE mcov_func_mcoord(FTYPE* X, int ii, int jj) // i not used
{
  extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);  //atch

  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE mcon[NDIM];
  FTYPE mcov[NDIM];
  struct of_geom geom;
  int i,j,k;
  FTYPE gcovpert[NDIM];


  //  dualfprintf(fail_file,"got here3.1: %d %d\n",ii,jj);
  gcov_func(1,MCOORD,X,gcovmcoord,gcovpert);

  //  dualfprintf(fail_file,"got here3.2\n");
  DLOOP(j,k) gengcov[j][k]=gcovmcoord[j][k];
  geom.gcov=gengcov;

  //  dualfprintf(fail_file,"got here3.3\n");
  mcon[TT]=0.0;
  mcon[RR]=0.0;
  mcon[TH]=0.0;
  mcon[PH]=1.0;
  //  dualfprintf(fail_file,"got here3.4\n");

  // lower only needs geom->gcov
  lower_vec(mcon,&geom,mcov);
  //  dualfprintf(fail_file,"got here3.5\n");

  return(mcov[jj]);
}

// returns MCOORD k_\mu form of k^\mu={1,0,0,0} value for jth element
FTYPE kcov_func_mcoord(FTYPE* X, int ii, int jj) // i not used
{
  extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);  //atch

  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE kcon[NDIM];
  FTYPE kcov[NDIM];
  struct of_geom geom;
  int i,j,k;
  FTYPE gcovpert[NDIM];

  gcov_func(1,MCOORD,X,gcovmcoord,gcovpert);

  DLOOP(j,k) gengcov[j][k]=gcovmcoord[j][k];
  geom.gcov=gengcov;

  kcon[TT]=1.0;
  kcon[RR]=0.0;
  kcon[TH]=0.0;
  kcon[PH]=0.0;

  // lower only needs geom->gcov
  lower_vec(kcon,&geom,kcov);

  return(kcov[jj]);
}
