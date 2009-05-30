
#include "decs.h"

/* bound array containing entire set of primitive variables */

// For fluxes, e.g. F1, assume fluxes exist everywhere -- including j/k boundary zones.  Only i-boundary zones need to be bounded.
// This assumes ZSLOOP(is,ie,js,je,ks,ke) is over boundary zones in flux.c, which in general to be compatible with any flux method (including finite volume) this is how it should be.

// in order to avoid accessing undefined data, but still fill corner
// zones, the ORDER of boundary LOOPS is as follows:

// X1 in&out: LOOPF2 LOOPF3 LOOPFACEBOUNDIN1 & LOOPFACEBOUNDOUT1
// X2 in&out: LOOPF1 LOOPF3 LOOPFACEBOUNDIN2 & LOOPFACEBOUNDOUT2
// X3 in&out: LOOPF1 LOOPF2 LOOPFACEBOUNDIN3 & LOOPFACEBOUNDOUT3

// With fluxes, only need to bound each dir-flux along that direction (as presently used by ENO-type schemes)

// Assume flux at 0 through N are computed correctly. So only need fluxes in other boundary zones.
// Self-assigns for 0 or N for simplicity of coding

// OUTFLOW leaves true edge of boundary unchanged
// Therefore, if FIXEDOUTFLOW, then extrapolation is always ok.
// if OUTFLOW, then extrapolation is ok as long as flux is from active zones out of boundary

// order of outflow extrap
// 0: none/ copy
// 1: first order
#define EXTRAP 0 //atch




int bound_flux_user(int boundstage, FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR])
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

  // first allow extrapolation
  //
  if (mycpupos[1] == 0) {
    if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)||(BCtype[X1DN]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* inner r boundary condition: u, just copy */
      LOOPF2 LOOPF3{
	ri=0; // since ri=0 is boundary flux
	rj=j;
	rk=k;
	LOOPFACEBOUND1IN PBOUNDLOOP(pl){
#if(EXTRAP==0)
	  // zeroth orde extrap
	  F1[i][j][k][pl] = F1[ri][rj][rk][pl];
#else
	  // linear extrap
	  F1[i][j][k][pl] = F1[ri][rj][rk][pl]+(F1[ri+1][rj][rk]-F1[ri][rj][rk])*(FTYPE)(i-ri);
#endif
	  // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
	  if((pl==RHO)&&(BCtype[X1DN]==OUTFLOW)){ // only for RHO is it obvious what to do without primitives
	    if(F1[i][j][k][pl]>0.0) F1[i][j][k][pl]=0.0; // GODMARK: hope this is enough to shut-down inflow (what about energy flux?)
	  }
	}// end pl and face loop
      }// end 2 3
    }// end if correct bound type
    else if(BCtype[X1DN]== ASYMM) {  //atch: added this; however unsure if any symmetry exists for fluxes
      LOOPN2 LOOPN3 {
	LOOPFACEBOUND1IN PBOUNDLOOP(pl)  {
	  F1[i][j][k][pl] = - F1[-i][j][k][pl];
	  if( U1 == pl || B1 == pl ) F1[i][j][k][pl] *= -1.0;  //resymm them
	}
      }
    }
    else if(BCtype[X1DN]== BCEXTRAP_VEL3 || BCtype[X1DN] == BCEXTRAP) { 
      //extrapolates fluxes with 6th order
      LOOPN2 LOOPN3{
	LOOPFACEBOUND1IN PBOUNDLOOP(pl) {
	  //	    F1[i][j][k][pl] = interpn( 6, i,  
	  //				       0, F1[0][j][k][pl], 
	  //				       1, F1[1][j][k][pl], 
	  //				       2, F1[2][j][k][pl], 
	  //				       3, F1[3][j][k][pl],
	  //				       4, F1[4][j][k][pl], 
	  //				       5, F1[5][j][k][pl]	);
	}
      }
    }
  }// end if mycpupos[1]==0


  ///////////////////////////
  //
  // X1 outer OUTFLOW/FIXEDOUTFLOW
  //
  ///////////////////////////


  // outer r BC:
  if (mycpupos[1] == ncpux1 - 1) {
    if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)||(BCtype[X1UP]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* outer r BC: outflow */

      LOOPF2 LOOPF3{
	ri=N1; // correct for flux
	rj=j;
	rk=k;
	LOOPFACEBOUND1OUT PBOUNDLOOP(pl){
#if(EXTRAP==0)
	  // zeroth orde extrap
	  F1[i][j][k][pl] = F1[ri][rj][rk][pl];
#else
	  // linear extrap
	  F1[i][j][k][pl] = F1[ri][rj][rk][pl]+(F1[ri][rj][rk]-F1[ri-1][rj][rk])*(FTYPE)(i-ri);
#endif
	  // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
	  if((pl==RHO)&&(BCtype[X1UP]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
	    if(F1[i][j][k][pl]<0.0) F1[i][j][k][pl]=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
	  }
	}// end pl and face loop
      }// end 2 3
    }// end if correct bound type
    else if(BCtype[X1UP]== ASYMM) {   //atch: added this; however unsure if any symmetry exists for fluxes
      LOOPN2 LOOPN3 {
	LOOPFACEBOUND1OUT PBOUNDLOOP(pl)  {
	  F1[i][j][k][pl] = - F1[N1 + N1 - i][j][k][pl];
	  if( U1 == pl || B1 == pl ) F1[i][j][k][pl] *= -1.0;  //resymm them
	}
      }
    }
    else if(BCtype[X1UP]== BCEXTRAP_VEL3 || BCtype[X1UP] == BCEXTRAP) { 
      //extrapolates fluxes with 6th order
      LOOPN2 LOOPN3{
	LOOPFACEBOUND1OUT PBOUNDLOOP(pl) {
	  //	    F1[i][j][k][pl] = interpn( 6, i, 
	  //				       N1,   F1[N1][j][k][pl], 
	  //				       N1-1, F1[N1-1][j][k][pl], 
	  //				       N1-2, F1[N1-2][j][k][pl], 
	  //				       N1-3, F1[N1-3][j][k][pl],
	  //				       N1-4, F1[N1-4][j][k][pl], 
	  //				       N1-5, F1[N1-5][j][k][pl]	);
	}
      }
    }
  }// end if mycpu is correct


  ///////////////////////////
  //
  // X2 inner POLARAXIS
  //
  ///////////////////////////


  if (mycpupos[2] == 0) {
    if((BCtype[X2DN]==OUTFLOW)||(BCtype[X2DN]==FIXEDOUTFLOW)||(BCtype[X2DN]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* inner 2 BC: outflow */

      LOOPF1 LOOPF3{
	ri=i; // correct for flux
	rj=0;
	rk=k;
	LOOPFACEBOUND2IN PBOUNDLOOP(pl){
#if(EXTRAP==0)
	  // zeroth orde extrap
	  F2[i][j][k][pl] = F2[ri][rj][rk][pl];
#else
	  // linear extrap
	  F2[i][j][k][pl] = F2[ri][rj][rk][pl]+(F2[ri][rj][rk]-F2[ri][rj-l][rk])*(FTYPE)(j-rj);
#endif
	  // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
	  if((pl==RHO)&&(BCtype[X2DN]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
	    if(F2[i][j][k][pl]<0.0) F2[i][j][k][pl]=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
	  }
	}// end pl and face loop
      }// end 2 3
    }// end if correct bound type
    else if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM) ){
      LOOPF1 LOOPF3{
	ri=i;
	rj=0;
	rk=k;
	LOOPFACEBOUND2IN PBOUNDLOOP(pl)  F2[i][j][k][pl] = F2[ri][rj+(rj-j)][rk][pl]; // self-assigns for j=0
      }
    }
    // F2 of U2/B2 is symmetric, so no need to antisymmetrize
    // F2 antisymmetric with respect to all other quantities

    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPF3{
	LOOPFACEBOUND2IN {
	  //if(POSDEFMETRIC==0){  //SASMARK need to asymmetrize the fluxes; don't understand why should not asymm. the fluxes if POSDEFMETRIC is 0
	  //	// u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	  //	// for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
	  //}
	  //else{
	  PLOOP(pl) F2[i][j][k][pl]*=-1; // anti-sym all
	  F2[i][j][k][U2]*=-1; // re-sym U2
	  F2[i][j][k][B2]*=-1; // re-sym B2
	  //}
	}
      }// end loop 13
    } // end if POLARXIS or ASYMM
  }// end if mycpupos[2]==0


  ///////////////////////////
  //
  // X2 outer POLARAXIS
  //
  ///////////////////////////


  if (mycpupos[2] == ncpux2-1) {
    if((BCtype[X2UP]==OUTFLOW)||(BCtype[X2UP]==FIXEDOUTFLOW)||(BCtype[X2UP]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* outer 2 BC: outflow */

      LOOPF1 LOOPF3{
	ri=i; // correct for flux
	rj=N2;
	rk=k;
	LOOPFACEBOUND2OUT PBOUNDLOOP(pl){
#if(EXTRAP==0)
	  // zeroth orde extrap
	  F2[i][j][k][pl] = F2[ri][rj][rk][pl];
#else
	  // linear extrap
	  F2[i][j][k][pl] = F2[ri][rj][rk][pl]+(F2[ri][rj][rk]-F2[ri][rj-l][rk])*(FTYPE)(j-rj);
#endif
	  // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
	  if((pl==RHO)&&(BCtype[X2UP]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
	    if(F2[i][j][k][pl]>0.0) F2[i][j][k][pl]=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
	  }
	}// end pl and face loop
      }// end 2 3
    }// end if correct bound type
    else if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM) ){
      LOOPF1 LOOPF3{
	ri=i;
	rj=N2-1;
	rk=k;
	LOOPFACEBOUND2OUT PBOUNDLOOP(pl)  F2[i][j][k][pl] = F2[ri][rj+(rj-j+2)][rk][pl]; // self-assigns for j=N
      }
    }

    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPF3{
	LOOPFACEBOUND2OUT {
	  //if(POSDEFMETRIC==0){  //SASMARK need to asymmetrize the fluxes; don't understand why should not asymm. the fluxes if POSDEFMETRIC is 0
	  //	// u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	  //	// for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
	  //}
	  //else{
	  PLOOP(pl) F2[i][j][k][pl]*=-1; // anti-sym all
	  F2[i][j][k][U2]*=-1; // re-sym U2
	  F2[i][j][k][B2]*=-1; // re-sym B2
	  //}
	}
      }// end loop 13
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
	LOOPFACEBOUND3IN PBOUNDLOOP(pl) F3[i][j][k][pl] = F3[ri][rj][rk+k][pl]; // for k=0 -> k=N3

	// copy from lower side to upper boundary zones
	ri=i;
	rj=j;
	rk=0;
	LOOPFACEBOUND3OUT PBOUNDLOOP(pl) F3[i][j][k][pl] = F3[ri][rj][rk+(k-N3)][pl]; // for k=N3 -> k=0
      }
    }
  }

  //x3 inner
  if ( mycpupos[3] == 0 ) {
    if((BCtype[X3DN]==OUTFLOW)||(BCtype[X3DN]==FIXEDOUTFLOW)||(BCtype[X3DN]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* inner 3 BC: outflow */

      LOOPF1 LOOPF2{
	ri=i; // correct for flux
	rj=j;
	rk=0;
	LOOPFACEBOUND3IN PBOUNDLOOP(pl){
#if(EXTRAP==0)
	  // zeroth orde extrap
	  F3[i][j][k][pl] = F3[ri][rj][rk][pl];
#else
	  // linear extrap
	  F3[i][j][k][pl] = F3[ri][rj][rk][pl]+(F3[ri][rj][rk]-F3[ri][rj][rk-l])*(FTYPE)(k-rk);
#endif
	  // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
	  if((pl==RHO)&&(BCtype[X3DN]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
	    if(F3[i][j][k][pl]<0.0) F3[i][j][k][pl]=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
	  }
	}// end pl and face loop
      }// end 2 3
    }// end if correct bound type
  }

  //x3 outer
  if ( mycpupos[3] == ncpux3 - 1 ) {
    if((BCtype[X3UP]==OUTFLOW)||(BCtype[X3UP]==FIXEDOUTFLOW)||(BCtype[X3UP]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* outer 3 BC: outflow */

      LOOPF1 LOOPF2{
	ri=i; // correct for flux
	rj=j;
	rk=N3;
	LOOPFACEBOUND3OUT PBOUNDLOOP(pl){
#if(EXTRAP==0)
	  // zeroth orde extrap
	  F3[i][j][k][pl] = F3[ri][rj][rk][pl];
#else
	  // linear extrap
	  F3[i][j][k][pl] = F3[ri][rj][rk][pl]+(F3[ri][rj][rk]-F3[ri][rj][rk-l])*(FTYPE)(k-rk);
#endif
	  // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
	  if((pl==RHO)&&(BCtype[X3UP]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
	    if(F3[i][j][k][pl]>0.0) F3[i][j][k][pl]=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
	  }
	}// end pl and face loop
      }// end 2 3
    }// end if correct bound type
  }

  return (0);
}
