//1. Created new outflow type: RAMESHOUTFLOW
//2. Corrected RAMESHOUTFLOW (mistyped k for rk)


#include "decs.h"

/* bound array containing entire set of primitive variables */

// GODMARK: something seriously wrong with EXTRAP=1 (EOMFFDE)

#define EXTRAP 0
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())

// to help protect the pole from death blows to the computational grid
// a sort of crushing regularization
#define POLEDEATH 0
// causes problems with stability at just beyond pole
// for field line plots, can just set B^\theta=0 along pole


// in order to avoid accessing undefined data, but still fill corner
// zones, the ORDER of boundary LOOPS is as follows:

// X1 in&out: LOOPN2 LOOPN3 LOOPBOUNDIN1 & LOOPBOUNDOUT1
// X2 in&out: LOOPF1 LOOPN3 LOOPBOUNDIN2 & LOOPBOUNDOUT2  // LOOPF1 ok if X1 dimension not there, then LOOPF1->LOOPN1
// X3 in&out: LOOPF1 LOOPF2 LOOPBOUNDIN3 & LOOPBOUNDOUT3  // as above


///Called before the MPI boundary routines
int bound_prim_user(int boundstage, FTYPE prim[][N2M][N3M][NPR])
{
  extern void bl_coord_2d(FTYPE *X, FTYPE *V1, FTYPE *V2);
  FTYPE rr, rth, r, th;

  int i, j, k, pl;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  FTYPE prescale[NPR];
  FTYPE X[NPR];







  ///////////////////////////
  //
  // X1 inner OUTFLOW/FIXEDOUTFLOW
  //
  ///////////////////////////

  if (mycpupos[1] == 0) {
    if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)){
      /* inner r boundary condition: u, just copy */
      LOOPN2 LOOPN3{
#if(EXTRAP==0)
        ri=0;
        rj=j;
        rk=k;
        LOOPBOUND1IN PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
#elif(EXTRAP==1)

        ri=0;
        rj=j;
        rk=k;
        LOOPBOUND1IN{
          for(pl=RHO;pl<=UU;pl++){
            prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][CENT]/gdet[i][j][k][CENT] ;
          }
          pl=U1;	    // treat U1 as special
          prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. - (i-ri)*dx[1]) ;
          for(pl=U2;pl<=U3;pl++){
            prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. + (i-ri)*dx[1]) ;
          }
          pl=B1; // treat B1 special
          prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][CENT]/gdet[i][j][k][CENT] ;
          for(pl=B2;pl<=B3;pl++){
            prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. + (i-ri)*dx[1]) ;
          }
        }
#elif(EXTRAP==2)
        ri=0;
        rj=j;
        rk=k;
        get_geometry(ri, rj, rk, CENT, &rgeom);
        rescale(1,1,prim[ri][rj][rk],&rgeom,prescale);
        LOOPBOUND1IN{
          // set guess
          PBOUNDLOOP(pl) prim[i][j][k][pl]=prim[ri][rj][k][pl];
          get_geometry(i, j, k, CENT, &geom);	    
          rescale(-1,1,prim[i][j][k],&geom,prescale);
        }
#endif
        LOOPBOUND1IN{
#if(WHICHVEL==VEL4)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_4vel(1,prim[i][j][k],&geom, 0) ;
#elif(WHICHVEL==VEL3)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_3vel(1,prim[i][j][k],&geom, 0) ;
          // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
          if(jonchecks){
            //fixup1zone(prim[i][j][k],&geom,0);
            failreturn=check_pr(prim[i][j][k],prim[i][j][k],&geom,-3);
            if(failreturn){
              dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
              if (fail(FAIL_BCFIX) >= 1) return (1);
            }
          }
#endif
#elif(WHICHVEL==VELREL4)
          get_geometry(i,j,k,CENT,&geom) ;
          inflow_check_rel4vel(1,prim[i][j][k],&geom,0) ;
          if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom,0)>=1)
            FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
        }
      }// end 2 3
    }
    else if(BCtype[X1DN]==FIXED){
      // inner radial fixed boundary conditions
      LOOPN2 LOOPN3{
        LOOPBOUND1IN PBOUNDLOOP(pl) prim[i][j][k][pl] = panalytic[i][j][k][pl];
      }
      //PBOUNDLOOP prim[i][j][k] = p[i][j][k];
    }// end if correct bound type
    else {
      dualfprintf( fail_file, "Unknown BCtype[X1DN] = %d\n", BCtype[X1DN] );
    }

  }// end if mycpupos[1]==0


  ///////////////////////////
  //
  // X1 outer OUTFLOW/FIXEDOUTFLOW
  //
  ///////////////////////////


  // outer r BC:
  if (mycpupos[1] == ncpux1 - 1) {
    if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)){
      /* outer r BC: outflow */

      LOOPN2 LOOPN3{
#if(EXTRAP==0)
        ri=N1-1;
        rj=j;
        rk=k;
        LOOPBOUND1OUT PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
#elif(EXTRAP==1)
        ri=N1-1;
        rj=j;
        rk=k;
        LOOPBOUND1OUT{
          for(pl=RHO;pl<=UU;pl++){
            prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][CENT]/gdet[i][j][k][CENT] ;
          }
          pl=U1; // treat U1 as special
          prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. - 2*(i-ri)*dx[1]) ;
          for(pl=U2;pl<=U3;pl++){
            prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. - (i-ri)*dx[1]) ;
          }
          pl=B1; // treat B1 special
          prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][CENT]/gdet[i][j][k][CENT] ;
          for(pl=B2;pl<=B3;pl++){
            prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. - (i-ri)*dx[1]) ;
          }
        }
#elif(EXTRAP==2)
        ri=N1-1;
        rj=j;
        rk=k;
        get_geometry(ri, rj, rk, CENT, &rgeom);
        rescale(1,1,prim[ri][rj][rk],&rgeom,prescale);
        LOOPBOUND1OUT{
          // set guess
          PBOUNDLOOP(pl) prim[i][j][k][pl]=prim[ri][rj][rk][pl];
          get_geometry(i, j, k, CENT, &geom);
          rescale(-1,1,prim[i][j][k],&geom,prescale);
        }
#endif

        LOOPBOUND1OUT{
#if(WHICHVEL==VEL4)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_4vel(1,prim[i][j][k],&geom,0) ;
#elif(WHICHVEL==VEL3)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_3vel(1,prim[i][j][k],&geom,0) ;
          // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
          if(jonchecks){
            //fixup1zone(prim[i][j][k],&geom,0);
            failreturn=check_pr(prim[i][j][k],prim[i][j][k],&geom,-3);
            if(failreturn){
              dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
              if (fail(FAIL_BCFIX) >= 1) return (1);
            }
          }
#endif
#elif(WHICHVEL==VELREL4)
          get_geometry(i,j,k,CENT,&geom) ;
          inflow_check_rel4vel(1,prim[i][j][k],&geom,0) ;
          if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom, 0)>=1)
            FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif	
        }
      }// end 2 3
    }// end if correct bound type
    else if( BCtype[X1UP]==RAMESHOUTFLOW ){
      LOOPN2 LOOPN3{
        ri=N1-1;
        rj=j;
        rk=k;
	      get_geometry(ri, rj, rk, CENT, &rgeom);
	      coord(ri, rj, rk, CENT, X);
	      bl_coord_2d(X, &rr, &rth);

        LOOPBOUND1OUT PBOUNDLOOP(pl) {
	        //	    prim[i][j][k] = (prim[ri][rj][k]-prim[ri-1][rj][k])*(i-ri)+prim[ri][rj][k];
	        // interpolator state
	        get_geometry(i, j, k, CENT, &geom);
	        coord(i, j, k, CENT, X);
	        bl_coord_2d(X, &r, &th);

	        // r B^{\hat{r}} is nearly constant (for nu=1.0)
	        // GODMARK
	        // Ramesh
	        if(pl==B1) prim[i][j][k][pl] = prim[ri][rj][rk][pl]* (sqrt(rgeom.gcov[1][1])*pow(rr,2.0-nu))/(sqrt(geom.gcov[1][1])*pow(r,2.0-nu));
	        else if(pl==B2) prim[i][j][k][pl] = prim[ri][rj][rk][pl]* (sqrt(rgeom.gcov[2][2])*pow(rr,2.0-nu))/(sqrt(geom.gcov[2][2])*pow(r,2.0-nu));
	        else if(pl==B3) prim[i][j][k][pl] = prim[ri][rj][rk][pl]* (sqrt(rgeom.gcov[3][3])*pow(rr,2.0-nu))/(sqrt(geom.gcov[3][3])*pow(r,2.0-nu));
	        else prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	      }
      }
    } //end RAMESHOUTFLOW
    else {
      dualfprintf( fail_file, "Unknown BCtype[X1UP] = %d\n", BCtype[X1UP] );
    }
  }// end if mycpu is correct


  ///////////////////////////
  //
  // X2 inner POLARAXIS
  //
  ///////////////////////////


  /* inner polar BC (preserves u^t rho and u) */
  if (mycpupos[2] == 0) {
    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM) ){
      LOOPF1 LOOPN3{
        ri=i;
        rj=0;
        rk=k;
        LOOPBOUND2IN PBOUNDLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj+(rj-j-1)][rk][pl];
      }
    }
    else {
      dualfprintf( fail_file, "Unknown BCtype[X2DN] = %d\n", BCtype[X2DN] );
    }

    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPN3{
        LOOPBOUND2IN {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            prim[i][j][k][U2] *= -1.;  //SASMARKx: corrected to be the same as in ramt38
            prim[i][j][k][U3] *= -1.;  //generally, U2 should be anti-symmetric; B3 can be either symmetric or antisymmetric -- does not matter if in 2d
            prim[i][j][k][B2] *= -1.;
            prim[i][j][k][B3] *= -1.;
          }
          else{
            prim[i][j][k][U2] *= -1.;
            prim[i][j][k][U3] *= 1.;
            prim[i][j][k][B2] *= -1.;
            prim[i][j][k][B3] *= 1.;
          }
        }
      }// end loop 13

#if(POLEDEATH)
      // fixup
      LOOPF1 LOOPN3 {
        for (j = -N2BND; j < 0+POLEDEATH; j++) {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            prim[i][j][k][U2] *= 0;
            prim[i][j][k][B2] *= 0.;
          }
          else{
            prim[i][j][k][U2] *= 0.;
            prim[i][j][k][B2] *= 0.;
          }
        }
      }// end loop 13
#endif

    } // end if POLARXIS or ASYMM

  }// end if mycpupos[2]==0


  ///////////////////////////
  //
  // X2 outer POLARAXIS
  //
  ///////////////////////////


  if (mycpupos[2] == ncpux2-1) {
    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM) ){
      LOOPF1 LOOPN3{
        ri=i;
        rj=N2-1;
        rk=k;
        LOOPBOUND2OUT PBOUNDLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj+(rj-j+1)][rk][pl];
      }
    }
    else {
      dualfprintf( fail_file, "Unknown BCtype[X2UP] = %d\n", BCtype[X2UP] );
    }

    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPF3{
        LOOPBOUND2OUT {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            prim[i][j][k][U2] *= -1.;   //SASMARKx:  corrected to be the same as in ramt38
            prim[i][j][k][U3] *= -1.;
            prim[i][j][k][B2] *= -1.;
            prim[i][j][k][B3] *= -1.;
          }
          else{
            prim[i][j][k][U2] *= -1.;
            prim[i][j][k][U3] *= 1.;
            prim[i][j][k][B2] *= -1.;
            prim[i][j][k][B3] *= 1.;
          }
        }
      }// end loop 13

#if(POLEDEATH)
      // fixup
      LOOPF1 LOOPN3 {
        for (j = N2-POLEDEATH; j <= N2-1+N2BND; j++) {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            prim[i][j][k][U2] *= 0;
            prim[i][j][k][B2] *= 0.;
          }
          else{
            prim[i][j][k][U2] *= 0.;
            prim[i][j][k][B2] *= 0.;
          }
        }
      }// end loop 13
#endif

    } // end if POLARXIS or ASYMM
  }// end if mycpupos[2]==ncpux2-1



  // periodic x3
  if ( (mycpupos[3] == 0)&&(mycpupos[3] == ncpux3 - 1) ) {
    if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){
      // just copy from one side to another

      LOOPF1 LOOPF2{

        // copy from upper side to lower boundary zones
        ri=i;
        rj=j;
        rk=N3;
        LOOPBOUND3IN PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk+k][pl];

        // copy from lower side to upper boundary zones
        ri=i;
        rj=j;
        rk=0;
        LOOPBOUND3OUT PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk+(k-N3)][pl];
      }
    }
    else {
      dualfprintf( fail_file, "Unknown BCtype[X3UP] = %d, BCtype[X3DN] = %d\n", BCtype[X3UP], BCtype[X3DN] );
    }

  }

  return (0);
}




// see interpline.c
int apply_bc_line(int doinverse, int iterglobal, int recontype, int bs, int be, FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM])
{
  int flip_y(int iterglobal, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM]);

  if(doinverse==0){
    flip_y(iterglobal, recontype, bs, be, yin);
  }
  else{
    flip_y(iterglobal, recontype, bs, be, yin);
    flip_y(iterglobal, recontype, bs, be, yout);
  }

  return(0);

}


#include "reconstructeno.h"

int flip_y(int iterglobal, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM])
{
  int pl,myi;


#if( WENO_DIR_FLIP_CONS_SIGN_DN )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterglobal == WENO_DIR_FLIP_CONS_SIGN_DN && (recontype == CVT_C2A || recontype == CVT_A2C) && mycpupos[iterglobal] == 0 ) { 
    PLOOP(pl) 
      for( myi = bs; myi < 0; myi++ ) {
        y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif

#if( WENO_DIR_FLIP_CONS_SIGN_UP )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterglobal == WENO_DIR_FLIP_CONS_SIGN_UP && (recontype == CVT_C2A || recontype == CVT_A2C)  && mycpupos[iterglobal] == numbercpu[iterglobal] - 1 ) { 
    PLOOP(pl) 
      for( myi = N1*(iterglobal==1) + N2*(iterglobal==2) + N3*(iterglobal==3); myi <= be; myi++ ) {
        y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif


  return(0);

}


///Called after the MPI boundary routines
int bound_prim_user_after_mpi(int boundstage, FTYPE prim[][N2M][N3M][NPR])
{
  int bound_disk(FTYPE prim[][N2M][N3M][NPR]);

  bound_disk(prim);  //disk boundary condition

  return(0);
}


int bound_disk(FTYPE prim[][N2M][N3M][NPR])
{
  int i, j, k, pl;
  struct of_geom geom,rgeom,r2geom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j
  FTYPE prescale[NPR];
  FTYPE primtest[NPR];

  FTYPE Bcon[NDIM],Bcov[NDIM],ucon[NDIM],Bsq;
  int sanitycheck;
  FTYPE X[NDIM];
  FTYPE rcent,thcent;
  FTYPE r,th;
  FTYPE rr,rth;
  FTYPE rr2,rth2;
  extern FTYPE get_omegastar(struct of_geom *geom, FTYPE r, FTYPE th);
  FTYPE Mcon[NDIM][NDIM],Fcov[NDIM][NDIM],romegaf2;
  extern void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  extern void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE Mcon[NDIM][NDIM]);
  void bl_coord_2d(FTYPE *X, FTYPE *V1, FTYPE *V2);  //wrapper for the 3d version
  FTYPE realgammamax;
  int tscale;
  FTYPE fun;
  static int firsttime=1;
  int limit_3vel(FTYPE *Bcon, struct of_geom *geom, FTYPE *vcon, FTYPE *pr);
  extern int OBtopr(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  FTYPE prnew[NPR];
  FTYPE outflowvar1,outflowvar2,outflowvar3;
  FTYPE ftempa,ftempb;

#define DELTAJ 2

  //#define RDISKINNER (2.0*Rhor)
#define RDISKINNER (-1E30)
#define RDISKOUTER (1E30)

#define RNEARBLACKHOLE (-1E30) // outflow \Omega_F for r<RNEARBLACKHOLE (doesn't work apparently)

  // time while static boundary conditions
  //#define TSTATIC (15.0)
#define TSTATIC (0.0)

  // whether to use new method to limit v^i to be time-like
#define LIMIT3VELBOUND 1


  // 0 = 0th order prims
  // 1 = linear prims
  // 2 = 0th order ramesh self-sim
  // 3 = linear ramesh self-sim
#define BOUNDARYINTERPTYPE 3





  // check boundary situation (can't copy from nothing)
  if( ((startpos[2]==totalsize[2]/2)||(startpos[2]==totalsize[2]/2-1)) && ( (DELTAJ>=N1BND)||(DELTAJ>=N2BND)) ){
    dualfprintf(fail_file,"DELTAJ out of bounds\n");
    myexit(1);
  }

#if(1)  //SYmmetrization stuff //atch
  if(!firsttime){
    k = 0; //no looping along phi direction
    // make the underlying evolved state exactly symmetric
    LOOPF2 LOOPF1{ // loop over entire domain looking for equator
      if((i>=-N1BND)&&(i<N1+N1BND)&&( (startpos[2]+j-totalsize[2]/2<0)&&(startpos[2]+j-totalsize[2]/2>=-DELTAJ) ) ){ 
        // then at equator for an even-sized grid with 2 zones around equator for 4 zones total
        // one side of equator

        coord(i, j, k, CENT, X);
        bl_coord_2d(X, &r, &th);

        if((r>RDISKINNER)&&(r<RDISKOUTER)){

          ri=i;
          // symmetric partner across equator (for even or odd sized grids)
          // global j
          rj=totalsize[2]/2-(startpos[2]+j) + (totalsize[2]-1)/2;
          // local j
          rj=rj - startpos[2];
          rk=0;

          // average rather than choose one
          // GODMARK: really doesn't apply at t=0

          PLOOP(pl){
            if((pl==B1)||(pl==B3)||(pl==U2)){
              // antisymmetric
              prim[i][j][k][pl]=0.5*(prim[i][j][k][pl]-prim[ri][rj][rk][pl]);
              prim[ri][rj][rk][pl]=0.5*(prim[ri][rj][rk][pl]-prim[i][j][k][pl]);
            }
            else{
              // symmetric
              prim[i][j][k][pl]=0.5*(prim[i][j][k][pl]+prim[ri][rj][rk][pl]);
              prim[ri][rj][rk][pl]=0.5*(prim[i][j][k][pl]+prim[ri][rj][rk][pl]);
            }
          }
          //	dualfprintf(fail_file,"i=%d j=%d rj=%d\n",i,j,rj);
        }
      }
    }
  }
#endif



  /////////////////////////////
  //
  // Set B^i
  //
  //////////////////////////////

  k = 0;

  LOOPF2 LOOPF1{ // loop over entire domain looking for equator


    coord(i, j, k, CENT, X);
    bl_coord_2d(X, &rcent, &thcent);

    // "lower"-j side (smaller j's)
    // one side of equator and outside horizon
    if((i>=-N1BND)&&(i<N1+N1BND)&&( (startpos[2]+j-totalsize[2]/2<0)&&(startpos[2]+j-totalsize[2]/2>=-DELTAJ) ) ){ 
      // then at equator for an even-sized grid with 2 zones around equator for 4 zones total
      // one side of equator

      //    coord(i, j, CENT, X);
      //    bl_coord(X, &r, &th);
      //    get_geometry(i, j, CENT, &geom);


      if((rcent>RDISKINNER)&&(rcent<RDISKOUTER)){
        // OUTFLOW some quantities into conductor that is the current sheet
        ri=i;
        rj=-(DELTAJ+1)-startpos[2]+totalsize[2]/2;
        rk=0;

        // outflow B^r
        //	  prim[i][j][B1]=prim[ri][rj][B1]+(prim[ri][rj][B1]-prim[ri][rj-1][B1])*(j-rj);
        //	  prim[i][j][B1]=prim[ri][rj][B1];
        // outflow B^\phi
        //prim[i][j][B3]=prim[ri][rj][B3]+(prim[ri][rj][B3]-prim[ri][rj-1][B3])*(j-rj);
        //	  prim[i][j][B3]=prim[ri][rj][B3];

        // define function that becomes 1.0 at TSTATIC
        //	if(t<=TSTATIC) fun = 1.0*sin(M_PI*0.5*t/TSTATIC);
        if(t<=TSTATIC) fun = 0.0;
        else fun=1.0;
        if(fun<0.0) fun=0.0;
        if(fun>1.0) fun=1.0;

        // smoothly turn on "outflow" boundary condition to avoid transients
        // fix B^r


        if(BOUNDARYINTERPTYPE==2){
          // 0th order with Ramesh dependence

          // reference state
          get_geometry(ri, rj, rk, CENT, &rgeom);
          coord(ri, rj, rk, CENT, X);
          bl_coord_2d(X, &rr, &rth);

          // interpolate state
          get_geometry(i, j, k, CENT, &geom);
          coord(i, j, k, CENT, X);
          bl_coord_2d(X, &r, &th);

          outflowvar1 = prim[ri][rj][rk][B1]* (sqrt(rgeom.gcov[1][1])*pow(rr,2.0-nu))/(sqrt(geom.gcov[1][1])*pow(r,2.0-nu));
          outflowvar3 = prim[ri][rj][rk][B3]* (sqrt(rgeom.gcov[3][3])*pow(rr,2.0-nu))/(sqrt(geom.gcov[3][3])*pow(r,2.0-nu));
        }
        else if(BOUNDARYINTERPTYPE==3){
          // 1st order with Ramesh dependence

          // reference state
          get_geometry(ri, rj, rk, CENT, &rgeom);
          coord(ri, rj, rk, CENT, X);
          bl_coord_2d(X, &rr, &rth);

          // other reference state for interpolation
          get_geometry(ri, rj-1, rk, CENT, &r2geom);
          coord(ri, rj-1, rk, CENT, X);
          bl_coord_2d(X, &rr2, &rth2);

          // interpolate state
          get_geometry(i, j, k, CENT, &geom);
          coord(i, j, k, CENT, X);
          bl_coord_2d(X, &r, &th);


          ftempa=prim[ri][rj][rk][B1]* (sqrt(rgeom.gcov[1][1])*pow(rr,2.0-nu));
          ftempb=prim[ri][rj-1][rk][B1]* (sqrt(r2geom.gcov[1][1])*pow(rr2,2.0-nu));

          outflowvar1 = ftempa + (ftempa - ftempb)*(j-rj); // (th-rth);
          outflowvar1 /=(sqrt(geom.gcov[1][1])*pow(r,2.0-nu));

          ftempa=prim[ri][rj][rk][B3]* (sqrt(rgeom.gcov[3][3])*pow(rr,2.0-nu));
          ftempb=prim[ri][rj-1][rk][B3]* (sqrt(r2geom.gcov[3][3])*pow(rr2,2.0-nu));

          outflowvar3 = ftempa + (ftempa - ftempb)*(j-rj); // (th-rth);
          outflowvar3 /=(sqrt(geom.gcov[3][3])*pow(r,2.0-nu));
        }
        else if(BOUNDARYINTERPTYPE==1){
          // linear interp
          outflowvar1 = prim[ri][rj][rk][B1]+(prim[ri][rj][rk][B1]-prim[ri][rj-1][rk][B1])*(j-rj);
          outflowvar3 = prim[ri][rj][rk][B3]+(prim[ri][rj][rk][B3]-prim[ri][rj-1][rk][B3])*(j-rj);
        }
        else if(BOUNDARYINTERPTYPE==0){
          // 0th order interp
          outflowvar1 = prim[ri][rj][rk][B1];
          outflowvar3 = prim[ri][rj][rk][B3];
        }

        // set B1
        prim[i][j][k][B1]=(1.0-fun)*panalytic[i][j][k][B1] + fun*outflowvar1;
        // set B3
        prim[i][j][k][B3]=(1.0-fun)*panalytic[i][j][k][B3] + fun*outflowvar3;

        //	if(!((rcent>RDISKINNER)&&(rcent<RDISKOUTER))){
        //	  prim[i][j][k][B2]=prim[ri][rj][rk][B2];
        //	}

#if(0)
        if(rcent<2.0*Rhor){
          // outflow B^\theta
          //	  prim[i][j][k][B2]=prim[ri][rj][rk][B2];
          prim[i][j][k][B2]=0; // monopole
        }
        else{
          // fix B^\theta
          prim[i][j][k][B2]=panalytic[i][j][k][B2];
        }
#endif


        // fix B^\theta
        prim[i][j][k][B2]=panalytic[i][j][k][B2];



#if(0)// outflow instead of evolved state
        // instead of letting evolution be default, let outflow'ed 4-velocity be default
        prim[i][j][k][U1]=prim[ri][rj][rk][U1];
        //	prim[i][j][k][U2]=prim[ri][rj][rk][U2];
        prim[i][j][k][U2]=0.0;
        prim[i][j][k][U3]=prim[ri][rj][rk][U3];
#endif

        if(rcent<=RNEARBLACKHOLE){ // outflow omegaf2 (use it below to set v^i)

          get_geometry(ri, rj, rk, CENT, &rgeom);

          prim[i][j][k][B2]=panalytic[i][j][k][B2];


          Max_con(prim[ri][rj][rk],&rgeom,Mcon);
          MtoF(3,Fcov,&rgeom,Mcon);
          romegaf2=Fcov[TT][TH]/Fcov[TH][PH];
        }

      }// end if within radius of disk
    }// end if one side of equator



    k = 0;

    // "upper"-j side (smaller j's)
    // other side of equator
    if((i>=-N1BND)&&(i<N1+N1BND)&&( (startpos[2]+j-totalsize[2]/2<=DELTAJ-1)&&(startpos[2]+j-totalsize[2]/2>=0) ) ){ 
      // then at equator for an even-sized grid with 2 zones around equator for 4 zones total



      if((rcent>RDISKINNER)&&(rcent<RDISKOUTER)){
        // OUTFLOW some quantities into conductor that is the current sheet
        ri=i;
        rj=(DELTAJ)-startpos[2]+totalsize[2]/2;
        rk=0;

        //    coord(i, j, k, CENT, X);
        //    bl_coord_2d(X, &r, &th);
        //    get_geometry(i, j, k, CENT, &geom);


        // outflow B^r
        //	  prim[i][j][k][B1]=prim[ri][rj][rk][B1]+(prim[ri][rj][rk][B1]-prim[ri][rj-1][rk][B1])*(j-rj);
        //	  prim[i][j][k][B1]=prim[ri][rj][rk][B1];
        // outflow B^\phi
        //prim[i][j][k][B3]=prim[ri][rj][rk][B3]+(prim[ri][rj][rk][B3]-prim[ri][rj-1][rk][B3])*(j-rj);
        //	  prim[i][j][k][B3]=prim[ri][rj][rk][B3];

        // define function that becomes 1.0 at TSTATIC
        //	if(t<=TSTATIC) fun = 1.0*sin(M_PI*0.5*t/TSTATIC);
        if(t<=TSTATIC) fun = 0.0;
        else fun=1.0;
        if(fun<0.0) fun=0.0;
        if(fun>1.0) fun=1.0;
        // smoothly turn on "outflow" boundary condition to avoid transients


        if(BOUNDARYINTERPTYPE==2){
          // 0th order with Ramesh dependence

          // reference state
          get_geometry(ri, rj, rk, CENT, &rgeom);
          coord(ri, rj, rk, CENT, X);
          bl_coord_2d(X, &rr, &rth);

          // interpolate state
          get_geometry(i, j, k, CENT, &geom);
          coord(i, j, k, CENT, X);
          bl_coord_2d(X, &r, &th);

          outflowvar1 = prim[ri][rj][rk][B1]* (sqrt(rgeom.gcov[1][1])*pow(rr,2.0-nu))/(sqrt(geom.gcov[1][1])*pow(r,2.0-nu));
          outflowvar3 = prim[ri][rj][rk][B3]* (sqrt(rgeom.gcov[3][3])*pow(rr,2.0-nu))/(sqrt(geom.gcov[3][3])*pow(r,2.0-nu));
        }
        else if(BOUNDARYINTERPTYPE==3){
          // 1st order with Ramesh dependence

          // reference state
          get_geometry(ri, rj, rk, CENT, &rgeom);
          coord(ri, rj, rk, CENT, X);
          bl_coord_2d(X, &rr, &rth);

          // other reference state for interpolation
          get_geometry(ri, rj+1, rk, CENT, &r2geom);
          coord(ri, rj+1, rk, CENT, X);
          bl_coord_2d(X, &rr2, &rth2);

          // interpolate state
          get_geometry(i, j, k, CENT, &geom);
          coord(i, j, k, CENT, X);
          bl_coord_2d(X, &r, &th);


          ftempa=prim[ri][rj][rk][B1]* (sqrt(rgeom.gcov[1][1])*pow(rr,2.0-nu));
          ftempb=prim[ri][rj+1][rk][B1]* (sqrt(r2geom.gcov[1][1])*pow(rr2,2.0-nu));

          outflowvar1 = ftempa + (ftempb - ftempa)*(j-rj); //(th-rth);
          outflowvar1 /=(sqrt(geom.gcov[1][1])*pow(r,2.0-nu));

          ftempa=prim[ri][rj][rk][B3]* (sqrt(rgeom.gcov[3][3])*pow(rr,2.0-nu));
          ftempb=prim[ri][rj+1][rk][B3]* (sqrt(r2geom.gcov[3][3])*pow(rr2,2.0-nu));

          outflowvar3 = ftempa + (ftempb - ftempa)*(j-rj); //(th-rth);
          outflowvar3 /=(sqrt(geom.gcov[3][3])*pow(r,2.0-nu));
        }
        else if(BOUNDARYINTERPTYPE==1){
          // linear interp
          outflowvar1 = prim[ri][rj][rk][B1]+(prim[ri][rj+1][rk][B1]-prim[ri][rj][rk][B1])*(j-rj);
          outflowvar3 = prim[ri][rj][rk][B3]+(prim[ri][rj+1][rk][B3]-prim[ri][rj][rk][B3])*(j-rj);
        }
        else if(BOUNDARYINTERPTYPE==0){
          // 0th order interp
          outflowvar1 = prim[ri][rj][rk][B1];
          outflowvar3 = prim[ri][rj][rk][B3];
        }

        // set B1
        prim[i][j][k][B1]=(1.0-fun)*panalytic[i][j][k][B1] + fun*outflowvar1;
        // set B3
        prim[i][j][k][B3]=(1.0-fun)*panalytic[i][j][k][B3] + fun*outflowvar3;



        //	if(!((rcent>RDISKINNER)&&(rcent<RDISKOUTER))){
        //	  prim[i][j][k][B2]=prim[ri][rj][rk][B2];
        //	}

#if(0)
        if(rcent<2.0*Rhor){
          // outflow B^\theta
          //	  prim[i][j][k][B2]=prim[ri][rj][rk][B2];
          prim[i][j][k][B2]=0; // monopole
        }
        else{
          // fix B^\theta
          prim[i][j][k][B2]=panalytic[i][j][k][B2];
        }
#endif

        // fix B^\theta
        prim[i][j][k][B2]=panalytic[i][j][k][B2];


#if(0)// outflow instead of evolved state
        // instead of letting evolution be default, let outflow'ed 4-velocity be default
        prim[i][j][k][U1]=prim[ri][rj][rk][U1];
        //	prim[i][j][U2]=prim[ri][rj][rk][U2];
        prim[i][j][k][U2]=0.0;
        prim[i][j][k][U3]=prim[ri][rj][rk][U3];
#endif


        if(rcent<=RNEARBLACKHOLE){ // outflow omegaf2 (use it below to set v^i)

          get_geometry(ri, rj, rk, CENT, &rgeom);

          prim[i][j][k][B2]=panalytic[i][j][k][B2];

          Max_con(prim[ri][rj][rk],&rgeom,Mcon);
          MtoF(3,Fcov,&rgeom,Mcon);
          romegaf2=Fcov[TT][TH]/Fcov[TH][PH];
        }


      }

#if(0)
      if((rcent>RDISKINNER)&&(rcent<RDISKOUTER)&&(i>=-N1BND)&&(i<N1+N1BND)&&( (abs(startpos[2]+j-(totalsize[2]-1)/2)<=DELTAJ) ) ){

        // OUTFLOW some quantities into conductor that is the current sheet

        // just average B^r
        prim[i][j][k][B1]=0.5*(prim[i][j+1][k][B1]+prim[i][j-1][k][B1]);
        prim[i][j][k][B3]=0.5*(prim[i][j+1][k][B3]+prim[i][j-1][k][B3]);

        if(!((rcent>RDISKINNER)&&(rcent<RDISKOUTER))){
          prim[i][j][k][B2]=0.5*(prim[i][j+1][k][B2]+prim[i][j-1][k][B2]);
        }

      }
#endif
    }

  }// end LOOPF2/1



  //////////
  ///  STAR
  /////////


  // loop must be independent because uses values defined from above (at least as setup right now -- could have loops above not go over star)
  //
  k = 0;
  LOOPF2 LOOPF1{ // loop over entire domain looking for star/axis

    coord(i, j, k, CENT, X);
    bl_coord_2d(X, &rcent, &thcent);


    //////////////////////////
    //
    // "fixed" axis or "fixed" stellar surface
    //
    ///////////////////////////
    if(1&&(BCtype[X1DN]==FIXED)){
      if(startpos[1]+i<0){

        ri=0;
        rj=j;
        rk=0;

        // define function that becomes 1.0 at TSTATIC
        //	if(t<=TSTATIC) fun = 1.0*sin(M_PI*0.5*t/TSTATIC);
        if(t<=TSTATIC) fun = 0.0;
        else fun=1.0;
        if(fun<0.0) fun=0.0;
        if(fun>1.0) fun=1.0;

        // smoothly turn on "outflow" boundary condition to avoid transients


        if(BOUNDARYINTERPTYPE==2){
          // 0th order with Ramesh dependence

          // reference state
          get_geometry(ri, rj, rk, CENT, &rgeom);
          coord(ri, rj, rk, CENT, X);
          bl_coord_2d(X, &rr, &rth);

          // interpolate state
          get_geometry(i, j, k, CENT, &geom);
          coord(i, j, k, CENT, X);
          bl_coord_2d(X, &r, &th);

          outflowvar2 = prim[ri][rj][rk][B2]* (sqrt(rgeom.gcov[2][2])*pow(rr,2.0-nu))/(sqrt(geom.gcov[2][2])*pow(r,2.0-nu));
          outflowvar3 = prim[ri][rj][rk][B3]* (sqrt(rgeom.gcov[3][3])*pow(rr,2.0-nu))/(sqrt(geom.gcov[3][3])*pow(r,2.0-nu));
        }
        else if(BOUNDARYINTERPTYPE==3){
          // 1st order with Ramesh dependence

          // reference state
          get_geometry(ri, rj, rk, CENT, &rgeom);
          coord(ri, rj, rk, CENT, X);
          bl_coord_2d(X, &rr, &rth);

          // other reference state for interpolation
          get_geometry(ri+1, rj, rk, CENT, &r2geom);
          coord(ri+1, rj, rk, CENT, X);
          bl_coord_2d(X, &rr2, &rth2);

          // interpolate state
          get_geometry(i, j, k, CENT, &geom);
          coord(i, j, k, CENT, X);
          bl_coord_2d(X, &r, &th);


          ftempa=prim[ri][rj][rk][B2]* (sqrt(rgeom.gcov[2][2])*pow(rr,2.0-nu));
          ftempb=prim[ri+1][rj][rk][B2]* (sqrt(r2geom.gcov[2][2])*pow(rr2,2.0-nu));

          outflowvar2 = ftempa + (ftempa - ftempb)*(i-ri); //(r-rr);
          outflowvar2 /=(sqrt(geom.gcov[2][2])*pow(r,2.0-nu));

          ftempa=prim[ri][rj][rk][B3]* (sqrt(rgeom.gcov[3][3])*pow(rr,2.0-nu));
          ftempb=prim[ri+1][rj][rk][B3]* (sqrt(r2geom.gcov[3][3])*pow(rr2,2.0-nu));

          outflowvar3 = ftempa + (ftempb - ftempa)*(i-ri); //(r-rr);
          outflowvar3 /=(sqrt(geom.gcov[3][3])*pow(r,2.0-nu));
        }
        else if(BOUNDARYINTERPTYPE==1){
          // linear interp
          outflowvar2 = prim[ri][rj][rk][B2]+(prim[ri+1][rj][rk][B2]-prim[ri][rj][rk][B2])*(i-ri);
          outflowvar3 = prim[ri][rj][rk][B3]+(prim[ri+1][rj][rk][B3]-prim[ri][rj][rk][B3])*(i-ri);
        }
        else if(BOUNDARYINTERPTYPE==0){
          // 0th order interp
          outflowvar2 = prim[ri][rj][rk][B2];
          outflowvar3 = prim[ri][rj][rk][B3];
        }

        // set B2
        prim[i][j][k][B2]=(1.0-fun)*panalytic[i][j][k][B2] + fun*outflowvar2;
        // set B3
        prim[i][j][k][B3]=(1.0-fun)*panalytic[i][j][k][B3] + fun*outflowvar3;



        //	if(!((rcent>RDISKINNER)&&(rcent<RDISKOUTER))){
        //	  prim[i][j][k][B2]=prim[ri][rj][rk][B2];
        //	}

        // fix B^r
        prim[i][j][k][B1]=panalytic[i][j][k][B1];


      }// end if near axis
    }// end if fixed near axis




  }// end LOOPF2/1





  ////////////////////////////////////
  //
  //      GET VELOCITY FROM FIELD AND OMEGAF
  //
  /////////////////////////////////////


  k = 0;

  LOOPF2 LOOPF1{ // loop over entire domain looking for equator

    coord(i, j, k, CENT, X);
    bl_coord_2d(X, &rcent, &thcent);
    get_geometry(i, j, k, CENT, &geom);

    //	if(abs(startpos[2]+j-totalsize/2)<=2){ // then at equator
    if(
      // then at equator for an even-sized grid with 2 zones around equator for 4 zones total
      ((i>=-N1BND)&&(i<N1+N1BND)&&( (startpos[2]+j-totalsize[2]/2<=DELTAJ-1)&&(startpos[2]+j-totalsize[2]/2>=-DELTAJ) ))||
      // then near axis
      (1&& (BCtype[X1DN]==FIXED)&&(startpos[1]+i<0) )
      ){


        if((rcent>RDISKINNER)&&(rcent<RDISKOUTER)){
          // outflow inner radial region till t=30, then try to fix
          /*
          if(
          ((rcent>2.0*Rhor) || (t>30.0))&&
          //	 ((rcent>3.0) || (t>30.0))&&
          ((rcent>RDISKINNER)&&(rcent<RDISKOUTER))
          ){
          */



          if((rcent>RDISKINNER)&&(rcent<RDISKOUTER)){ // assume done above in split loops
            // fix B^r, assume other field can do whatever they want.
            //	  prim[i][j][k][B1]=panalytic[i][j][k][B1];
            // fix B^\theta, assume other field can do whatever they want.
            //	  prim[i][j][k][B2]=panalytic[i][j][k][B2];
          }

          // fix E_\phi=0, which completely determines all components of v^i
          // try it, but don't have to force it except at late times.
          // since near horizon there is the light surface and the grid goes through this surface, we can't guarantee that our v^i is time-like.
          Bcon[0]=0;
          Bcon[1]=prim[i][j][k][B1];
          Bcon[2]=prim[i][j][k][B2];
          Bcon[3]=prim[i][j][k][B3];


          if(rcent>RNEARBLACKHOLE){ // outflow omegaf2 (use it below to set v^i)
            // must be same as in init.c
            Omegastar=get_omegastar(&geom,rcent,thcent);
          }
          else{
            Omegastar=romegaf2;
          }


#if(1)
          ///////////////////////////////
          //
          // new way to get velocity
          if(OBtopr(Omegastar,Bcon,&geom,prnew)>=1){
            dualfprintf(fail_file, "OBtopr(bounds_disk): space-like error in init_postfield()\n");
            dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
            failed=1;
            return(1);
          }
          // assign answer
          prim[i][j][k][U1]=prnew[U1];
          prim[i][j][k][U2]=prnew[U2];
          prim[i][j][k][U3]=prnew[U3];



          ///////////////////////
          //
          // old way to get velocity
#else 


          lower(Bcon,&geom,Bcov);
          Bsq=0.0;
          for(k=0;k<=3;k++) Bsq+=Bcon[k]*Bcov[k];

          ////////////////
          // damp the calculation of v^i by damping Bsq
          //	if(rcent<2.0*Rhor) Bsq*=1.0/(1.0-exp(-t/10)+SMALL);


          vcon[1]=-Bcon[1]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;
          if((rcent>RDISKINNER)&&(rcent<RDISKOUTER)) vcon[2]=-Bcon[2]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;
          else vcon[2]=0.0; // like AVOIDCS in phys.ffde.c (kinda stationary since v^\theta=0 on equator)

          vcon[3]=Omegastar-Bcon[3]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;

#if(LIMIT3VELBOUND)

          // guarenteed to return prim[i][j] that's time-like with Lorentz factor limited by GAMMAMAX
          limit_3vel(Bcon, &geom, vcon, prim[i][j][k]);

#else
          // check to make sure that chosen v^i is time-like.
          // allow modification of v^i since evolution toward stationary solution may not be adiabatic

          // get electric field from v and B ? Not sure how to get E^2(v^i,B^i)


          //	    for(k=1;k<=3;k++) dualfprintf(fail_file,"Bcon[%d]=%21.15g Bcov[%d]=%21.15g\n",k,Bcon[k],k,Bcov[k]);

          //	    for(k=1;k<=3;k++) dualfprintf(fail_file,"vcon[%d]=%21.15g Bsq=%21.15g\n",k,vcon[k],Bsq);

          //	    for(k=1;k<=3;k++) if(isnan(vcon[k])){
          //	      dualfprintf(fail_file,"nan encountered: k=%d vcon=%21.15g\n",k,vcon[k]);
          //	    }
          // sticks primitive velocity into prim[i][j][U1->U3]
          whocalleducon=1;
          sanitycheck=1;
          if(vcon2pr(WHICHVEL,vcon,&geom,primtest)>=1){
            failed=0; // don't fail
            //	    dualfprintf(fail_file,"v^\phi[fixed]=%g\n",vcon[PH]);
            //return(1);
            // keep old velocity
            sanitycheck=0;
            // mark when this type of problem
            if(DODEBUG){
              TSCALELOOP failfloorcount[i][j][k][tscale][COUNTUTOPRIMFAILRHONEG]++;
            }

            ///////////////////////
            //
            // keeping old velocity leads to non-stationary solution and noise from current sheet dissipating leading to large u^t in evolved solution
            //

          }
          else{ // then primtest is ok, but do sanity check first
            for(k=U1;k<=U3;k++) if(isnan(primtest[k])){
              dualfprintf(fail_file,"nan encountered: i=%d j=%d k=%d prim=%21.15g\n",i,j,k,primtest[k]);
              sanitycheck=0;
              // mark when this type of problem
              if(DODEBUG){
                TSCALELOOP failfloorcount[i][j][k][tscale][COUNTUTOPRIMFAILUNEG]++;
              }
            }
            if(sanitycheck){ // then primtest ok, so assign
              prim[i][j][k][U1]=primtest[U1];
              prim[i][j][k][U2]=primtest[U2];
              prim[i][j][k][U3]=primtest[U3];


            }
            // otherwise still avoid primtest

          }
          whocalleducon=0;

          // if no good stationary v^i, then reduce to AVOIDCS like behavior
          if(!sanitycheck){
            prim[i][j][U2]=0.0;
          }

          // if set boundary condition, then ignore failure from inversion
          // code in FFDE mode doesn't actually use pflag.
          if(sanitycheck){
            pflag[i][j][k][FLAGUTOPRIMFAIL]=UTOPRIMNOFAIL;
          }

#endif // end if old way of dealing with when v^i is not good (LIMIT3VELBOUND)


#if(0)
          /////////////////////////////
          //
          // additional limiting?
          //
          ////////////////////////////

          if(sanitycheck){ // make sure still ok gamma
            // force flow to not move too fast inside ergosphere
            if(rcent<2) realgammamax=3;
            else realgammamax=GAMMAMAX;

            // limit gamma
            limit_gamma(realgammamax,prim[i][j][k],&geom,-1.0);
          }
#endif

#endif // end if old way to get velocity

        }

    }// end if at equator
  } // end loop over domain


  firsttime=0;
  return(0);
}



//wrapper for the 3d version that provides the interface of a 2d version
void bl_coord_2d(FTYPE *X, FTYPE *V1, FTYPE *V2)
{
  FTYPE V[NDIM];

  bl_coord( X, V );
  *V1 = V[1];
  *V2 = V[2];

}
