
#include "decs.h"

/* bound array containing entire set of primitive variables */

// GODMARK: something seriously wrong with EXTRAP=1 (EOMFFDE)

#define EXTRAP 0
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())

// to help protect the pole from death blows to the computational grid
// a sort of crushing regularization
#define POLEDEATH 1
// causes problems with stability at just beyond pole
// for field line plots, can just set B^\theta=0 along pole


// in order to avoid accessing undefined data, but still fill corner
// zones, the ORDER of boundary LOOPS is as follows:

// X1 in&out: LOOPN2 LOOPN3 LOOPBOUNDIN1 & LOOPBOUNDOUT1
// X2 in&out: LOOPF1 LOOPN3 LOOPBOUNDIN2 & LOOPBOUNDOUT2  // LOOPF1 ok if X1 dimension not there, then LOOPF1->LOOPN1
// X3 in&out: LOOPF1 LOOPF2 LOOPBOUNDIN3 & LOOPBOUNDOUT3  // as above



int bound_prim_user(int boundstage, FTYPE prim[][N2M][N3M][NPR])
{
  int i, j, k, pl;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  FTYPE prescale[NPR];







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
    }// end if correct bound type

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

    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPN3{
	LOOPBOUND2IN {
	  if(POSDEFMETRIC==0){
	    // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	    // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
	    prim[i][j][k][U2] *= 1.;
	    prim[i][j][k][U3] *= 1.;
	    prim[i][j][k][B2] *= 1.;
	    prim[i][j][k][B3] *= 1.;
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

    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPF3{
	LOOPBOUND2OUT {
	  if(POSDEFMETRIC==0){
	    // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	    // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
	    prim[i][j][k][U2] *= 1.;
	    prim[i][j][k][U3] *= 1.;
	    prim[i][j][k][B2] *= 1.;
	    prim[i][j][k][B3] *= 1.;
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
