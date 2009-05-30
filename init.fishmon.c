
/* 
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


#define SLOWFAC 1.0		/* reduce u_phi by this amount */

SFTYPE rhomax=0,umax=0,bsq_max=0,beta,rin;

int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.2;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

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

  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;

  //  UTOPRIMVERSION =   UTOPRIM5D1;
  //UTOPRIMVERSION =   UTOPRIM2DFINAL;

  return(0);
}

int init_grid(void)
{
  SFTYPE rh;
  
  // metric stuff first
  a = 0.9375 ;
  

  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.0;
  Rhor=rhor_calc(0);
  Rout = 40.;
 
  Rin=setRin(setihor());
  //  Rin = 0.98 * Rhor;
  
  hslope = 0.3;

  // define coordinate type
  defcoord = 9;



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
  DOENOFLUX = ENOFINITEVOLUME;
  //  DOENOFLUX= NOENOFLUX;
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


  rescaletype=1;


  if(BCtype[X1UP]==FIXEDOUTFLOW){ // then doing bondi inflow
    // avoids constant floor activation -- trying to be more physical
    prfloorcoef[RHO]=RHOMIN/100.0;
    prfloorcoef[UU]=UUMIN/100.0;
  }
  else{
    prfloorcoef[RHO]=RHOMIN;
    prfloorcoef[UU]=UUMIN;
  }



  /* output choices */
  tf = 2000.0;

  DTd = 50.;			/* dumping frequency, in units of M */
  DTavg = 50.0;
  DTener = 2.0;			/* logfile frequency, in units of M */
  DTi = 2.0;			/* image file frequ., in units of M */
  DTdebug = 50.0; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;			/* restart file period in steps */

  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int pl;
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];


  get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
  set_atmosphere(0,WHICHVEL,&realgeom,pratm); // set velocity in chosen WHICHVEL frame in any coordinate system

  if(pr[RHO]<pratm[RHO]){
    PLOOP(pl) pr[pl]=pratm[pl];
  }
  

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}

int init_grid_post_set_grid(void)
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);

  
  // diagnostic
  // determine nature of inner radial edge (assumes myid==0 is always there)
  if(myid==0){
    i=-INFULL1;
    j=k=0;
    coord(i,j,k, FACE1, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    trifprintf("rmin: %21.15g\n", r);
    trifprintf("rmin/rh: %21.15g\n", r / Rhor );
    //    trifprintf("rmin/rsing: %21.15g\n", r / (a+SMALL));
    if(r/Rhor<=1.0){
      trifprintf("inner grid is inside horizon\n");
    }
    else{
      trifprintf("inner grid is outside horizon\n");
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
    for(l=1;l<=3;l++) init_vpot(l,i,j,k,&A[l][i][j][k]); // request vector potential for all computational corners
  }
  trifprintf("Initialize field from vector potential assign\n");

  init_vpot2field(A,p);

  normalize_field(p);




  return(0);


}



// unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr)
{
  SFTYPE randfact;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom realgeom,geom;
  

  /* for disk interior */
  SFTYPE l, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;
  SFTYPE rmax, lfish_calc(SFTYPE rmax);
  SFTYPE rh;
  //  FTYPE pratm[NPR];


  rin = 6. ;
  rmax = 12. ;
  l = lfish_calc(rmax) ;
  kappa = 1.e-3 ;
  beta = 1.e2 ;
  randfact = 4.e-2;
  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];



  sth = sin(th);
  cth = cos(th);

  /* calculate lnh */
  DD = r * r - 2. * r + a * a;
  AA = (r * r + a * a) * (r * r + a * a) - DD * a * a * sth * sth;
  SS = r * r + a * a * cth * cth;
  
  thin = M_PI / 2.;
  sthin = sin(thin);
  cthin = cos(thin);
  DDin = rin * rin - 2. * rin + a * a;
  AAin = (rin * rin + a * a) * (rin * rin + a * a)
    - DDin * a * a * sthin * sthin;
  SSin = rin * rin + a * a * cthin * cthin;
  
  if (r >= rin) {
    lnh = 0.5 * log((1. + sqrt(1. + 4. * (l * l * SS * SS) * DD /
			       (AA * sth * AA * sth))) / (SS * DD /
							  AA))
      - 0.5 * sqrt(1. +
		   4. * (l * l * SS * SS) * DD / (AA * AA * sth *
						  sth))
      - 2. * a * r * l / AA -
      (0.5 *
       log((1. +
	    sqrt(1. +
		 4. * (l * l * SSin * SSin) * DDin / (AAin * AAin *
						      sthin *
						      sthin))) /
	   (SSin * DDin / AAin))
       - 0.5 * sqrt(1. +
		    4. * (l * l * SSin * SSin) * DDin / (AAin *
							 AAin *
							 sthin *
							 sthin))
       - 2. * a * rin * l / AAin);
  } else
    lnh = 1.;
  

  
  /* regions outside torus */
  // this region is already in Kerr Schild prime in proper primitive quantity for velocity
  if (lnh < 0. || r < rin) {


    get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,&realgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  /* region inside magnetized torus; u^i is calculated in
     Boyer-Lindquist coordinates, as per Fishbone & Moncrief, so it
     needs to be transformed at the end */
  else {
    hm1 = exp(lnh) - 1.;
    rho = pow(hm1 * (gam - 1.) / (kappa * gam), 1. / (gam - 1.));
    u = kappa * pow(rho, gam) / (gam - 1.);
    ur = 0.;
    uh = 0.;
    
    /* calculate u^phi */
    expm2chi = SS * SS * DD / (AA * AA * sth * sth);
    up1 = sqrt((-1. + sqrt(1. + 4. * l * l * expm2chi)) / 2.);
    up = 2. * a * r * sqrt(1. + up1 * up1) / sqrt(AA * SS * DD) +
      sqrt(SS / AA) * up1 / sth;
    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0) - 0.5));
    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;

    *whichvel=VEL4;
    *whichcoord=BLCOORDS;
    return(0);
  }
}


#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2

#define FIELDTYPE DISKFIELD

// assumes normal field in pr
int init_vpot(int l, int i, int j, int k, FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom geom;





  *A=0;

  if(l==3){// A_\phi

    coord(i, j, k, CORN3, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];


    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISKVERT)){
      *A += 0.5*r*sin(th) ;
    }
    /* field-in-disk version */
    
    if((FIELDTYPE==DISKFIELD)||(FIELDTYPE==DISKVERT)){
      // average of density that lives on CORN3


      // since init_vpot() is called for all i,j,k, can't use
      // non-existence values, so limit averaging:
      if((i==-N1BND)&&(j==-N2BND)){
	rho_av = p[i][j][k][RHO];
      }
      else if(i==-N1BND){
	rho_av = AVGN_2(p,i,j,k,RHO);
      }
      else if(j==-N2BND){
	rho_av = AVGN_1(p,i,j,k,RHO);
      }
      else{ // normal cells
	rho_av = AVGN_for3(p,i,j,k,RHO);
      }

      q = rho_av / rhomax - 0.2;

      if (q > 0.)      *A += q;
    }
  }

  return(0);

}


int init_vpot2field(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE pr[][N2M][N3M][NPR])
{
  extern int vpot2field(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE p[][N2M][N3M][NPR]);


  return(vpot2field(A,pr));

}

// assumes we are fed the true densities
int normalize_densities(FTYPE p[][N2M][N3M][NPR])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;


  rhomax=0;
  umax=0;
  ZLOOP{
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];

    if (p[i][j][k][RHO] > rhomax)   rhomax = p[i][j][k][RHO];
    if (p[i][j][k][UU] > umax && r > rin)    umax = p[i][j][k][UU];
  }

  mpimax(&rhomax);
  mpimax(&umax);
  trifprintf("rhomax: %21.15g umax: %21.15g\n", rhomax, umax);


  ZLOOP{
    p[i][j][k][RHO] /= rhomax;
    p[i][j][k][UU] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;

  return(0);
}


// assumes normal field definition
int normalize_field(FTYPE p[][N2M][N3M][NPR])
{
  int i,j,k;
  FTYPE bsq_ij;
  SFTYPE bsq_max, norm, beta_act;
  struct of_geom geom;
  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;

  bsq_max = 0.;
  ZLOOP {
    get_geometry(i, j, k, CENT, &geom);    

    if(FIELDTYPE==VERTFIELD){
      coord(i, j, k, CENT, X);
      bl_coord(X, V);
      r=V[1];
      th=V[2];
      
      if((r>rin)&&(fabs(th-M_PI*0.5)<4.0*M_PI*dx[2]*hslope)){
	if (bsq_calc(p[i][j][k], &geom, &bsq_ij) >= 1)
	  FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
	
	if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
      }
    }
    else{
      if (bsq_calc(p[i][j][k], &geom, &bsq_ij) >= 1)
	FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      
      if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    }
  }

  mpimax(&bsq_max);
  trifprintf("initial bsq_max: %21.15g\n", bsq_max);

  /* finally, normalize to set field strength */
  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
  trifprintf("initial beta: %21.15g (should be %21.15g)\n", beta_act,beta);
  norm = sqrt(beta_act / beta);
  
  bsq_max = 0.;
  ZLOOP {
    p[i][j][k][B1] *= norm;
    p[i][j][k][B2] *= norm;
    p[i][j][k][B3] *= norm;

    get_geometry(i, j, k, CENT, &geom);
    if (bsq_calc(p[i][j][k], &geom, &bsq_ij) >= 1)
      FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
    if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    
  }
  mpimax(&bsq_max);
  trifprintf("new initial bsq_max: %21.15g\n", bsq_max);

  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);

  trifprintf("new bsq_max: %21.15g\n", bsq_max);
  trifprintf("final beta: %21.15g (should be %21.15g)\n", beta_act,beta);

  return(0);
}



#undef SLOWFAC

SFTYPE lfish_calc(SFTYPE r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
	   ((-2. * a * r * (pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2))) /
	    sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
	    ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) +
					  pow(a,
					      2) * (2. + r))) / sqrt(1 +
								     (2.
								      *
								      a)
								     /
								     pow
								     (r,
								      1.5)
								     -
								     3.
								     /
								     r)))
	  / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
	     (pow(a, 2) + (-2. + r) * r))
	  );
}
