
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



int bound_prim_user(int boundstage, FTYPE prim[][N2M][N3M][NPR])
{
	int i, j, k, pl;
	struct of_geom geom,rgeom;
	FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
	int failreturn;
#endif
	int ri, rj, rk; // reference i,j,k
	//FTYPE prescale[NPR];
	FTYPE pr0[NPRBOUND]; 
	FTYPE pr1[NPRBOUND];
	FTYPE pr2[NPRBOUND];
	FTYPE pr3[NPRBOUND];
	FTYPE pr4[NPRBOUND];
	FTYPE pr5[NPRBOUND];
	#define BONDIEXTRAPORDER 1
	FTYPE BondiC1[BONDIEXTRAPORDER+1], BondiC2[BONDIEXTRAPORDER+1], BondiK[BONDIEXTRAPORDER+1]; //for storing Bondi intergrals
	FTYPE BondiC1x = -0.00675;
	FTYPE BondiC2x = 1.373125;
	FTYPE BondiKx = 1.;
	FTYPE BondiC1interp;
	FTYPE BondiC2interp;
	FTYPE BondiKinterp;
	int ind;

	int asym_compute_1(FTYPE (*prim)[N2M][N3M][NPR]);
	int	bondi_generate_solution( FTYPE C1, FTYPE C2, FTYPE K, int i1, int j1, int k1, int face, FTYPE *primout );
	int	bondi_compute_integrals( int i0, int j0, int k0, FTYPE *primin, FTYPE *C1out, FTYPE *C2out, FTYPE *Kout );


	// periodic x1
	if ( (mycpupos[1] == 0)&&(mycpupos[1] == ncpux1 - 1) ) {
		if( (BCtype[X1DN]==PERIODIC)&&(BCtype[X1UP]==PERIODIC) ){
			// just copy from one side to another

			LOOPF2 LOOPF3{

				// copy from upper side to lower boundary zones
				ri=N1;
				rj=j;
				rk=k;
				LOOPBOUND1IN PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri+i][rj][rk][pl];

				// copy from lower side to upper boundary zones
				ri=0;
				rj=j;
				k=k;
				LOOPBOUND1OUT PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri+(i-N1)][rj][rk][pl];
			}
		}
	}

	///////////////////////////
	//
	// X1 inner OUTFLOW/FIXEDOUTFLOW
	//
	///////////////////////////

	if (mycpupos[1] == 0) {
		if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)||(BCtype[X1DN]==OUTFLOWNOINFLOW)){
			/* inner r boundary condition: u, just copy */
			LOOPN2 LOOPN3{
				ri=0;
				rj=j;
				rk=k;

				//avoid reading from the MPI boundary which is not set and may have NAN's
				if( !finite( prim[ri][rj][rk][0] ) ) continue;

				LOOPBOUND1IN PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
				if( BCtype[X1DN]==OUTFLOWNOINFLOW ) {
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
		}// end if correct bound type
		else if(BCtype[X1DN] == SYMM || BCtype[X1DN] == ASYMM || BCtype[X1DN] == POLARAXIS) { 
			LOOPN2 LOOPN3 {
				LOOPBOUND1IN PBOUNDLOOP(pl)  {
					prim[i][j][k][pl] = prim[-1-i][j][k][pl];
					if( BCtype[X1DN]==ASYMM || BCtype[X1DN]== POLARAXIS )  {
						if( pl == U1 || pl == B1 )	prim[i][j][k][pl] *= -1.;
					}
				}
			}
		}
		else if(BCtype[X1DN]== BONDIINTOUTFLOW) { 
			LOOPN2 LOOPN3{
				//compute Bondi integrals in active zones
				for( ind = 0; ind < BONDIEXTRAPORDER; ind++ ) {
					bondi_compute_integrals( ind, j, k, prim[ind][j][k], &BondiC1[ind], &BondiC2[ind], &BondiK[ind] );
				}
				//extrapolate C1 at BONDIEXRAPORDERth order
				LOOPBOUND1IN {
					//the last two do not matter since order = 5
					BondiC1interp = interpn( BONDIEXTRAPORDER, i, 0, BondiC1[0], 1,  BondiC1[1], 2,  BondiC1[2], 3,  BondiC1[3], 4, BondiC1[4], 5, BondiC1[5]); 
					BondiC2interp = interpn( BONDIEXTRAPORDER, i, 0, BondiC2[0], 1,  BondiC2[1], 2,  BondiC2[2], 3,  BondiC2[3], 4, BondiC2[4], 5, BondiC2[5]); 
					BondiKinterp =  interpn( BONDIEXTRAPORDER, i, 0, BondiK[0],  1,  BondiK[1],  2,  BondiK[2],  3,  BondiK[3],  4, BondiK[4],  5, BondiK[5]); 
					//outflow C1,C2,K from the grid
					bondi_generate_solution( BondiC1interp, BondiC2interp, BondiKinterp, i, j, k, CENT, prim[i][j][k] );  
				}
			}
		}
		else if(BCtype[X1DN]== BCEXTRAP) { 
			LOOPN2 LOOPN3 {
				//PBOUNDLOOP(pl) pr0[pl] = prim[0][j][k][pl];
				//PBOUNDLOOP(pl) pr1[pl] = prim[1][j][k][pl];
				//LOOPBOUND1IN PBOUNDLOOP(pl) prim[i][j][k][pl]=pr1[pl]*(i-0)+pr0[pl]*(1-i); // extrapolate all primitive quantities, linearly
				PBOUNDLOOP(pl) pr0[pl] = prim[0][j][k][pl];
				PBOUNDLOOP(pl) pr1[pl] = prim[1][j][k][pl];
				PBOUNDLOOP(pl) pr2[pl] = prim[2][j][k][pl];
				PBOUNDLOOP(pl) pr3[pl] = prim[3][j][k][pl];
				PBOUNDLOOP(pl) pr4[pl] = prim[4][j][k][pl];
				PBOUNDLOOP(pl) pr5[pl] = prim[5][j][k][pl];
				LOOPBOUND1IN {
					PBOUNDLOOP(pl) prim[i][j][k][pl] = interpn( 6, i, 0, pr0[pl], 1, pr1[pl], 2, pr2[pl], 3, pr3[pl], 4, pr4[pl], 5, pr5[pl] );
				}
			}
		}
		else if(BCtype[X1DN]== RESCALEOUTFLOW) { 
			LOOPN2 LOOPN3 {
				LOOPBOUND1IN {
					bounds_generate( i, j, k, prim[i][j][k] );
				}
				//relax them to connect to the solution in the active grid cells smoothly
				bounds_generate( 0, j, k, pr0 );
				LOOPBOUND1IN PBOUNDLOOP(pl) {
					prim[i][j][k][pl] += (prim[0][j][k][pl] - pr0[pl]);
				}
			}
		}
		else if(BCtype[X1DN]== BCEXTRAP_VEL3) { 
			LOOPN2 LOOPN3 {
				PBOUNDLOOP(pl) pr0[pl] = prim[0][j][k][pl];
				PBOUNDLOOP(pl) pr1[pl] = prim[1][j][k][pl];
				PBOUNDLOOP(pl) pr2[pl] = prim[2][j][k][pl];
				PBOUNDLOOP(pl) pr3[pl] = prim[3][j][k][pl];
				PBOUNDLOOP(pl) pr4[pl] = prim[4][j][k][pl];
				PBOUNDLOOP(pl) pr5[pl] = prim[5][j][k][pl];
				metp2met2bl( VEL3, MCOORD, pr0, 0, j, k );
				metp2met2bl( VEL3, MCOORD, pr1, 1, j, k );
				metp2met2bl( VEL3, MCOORD, pr2, 2, j, k );
				metp2met2bl( VEL3, MCOORD, pr3, 3, j, k );
				metp2met2bl( VEL3, MCOORD, pr4, 4, j, k );
				metp2met2bl( VEL3, MCOORD, pr5, 5, j, k );
				LOOPBOUND1IN {
					PBOUNDLOOP(pl) prim[i][j][k][pl] = interpn( 6, i, 0, pr0[pl], 1, pr1[pl], 2, pr2[pl], 3, pr3[pl], 4, pr4[pl], 5, pr5[pl] );
					//PBOUNDLOOP(pl) prim[i][j][k][pl]=pr1[pl]*(i-0)+pr0[pl]*(1-i); // extrapolate all primitive quantities, linearly
					bl2met2metp2v( VEL3, MCOORD, prim[i][j][k], i, j, k );
				}
			}
		}
		else if(BCtype[X1DN]== FIXED) { 
			//calls init_dsandvels() to regenerate the solution
			LOOPN2 LOOPN3{
				LOOPBOUND1IN bounds_generate( i, j, k, prim[i][j][k] );  //SASMARK
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
		if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)||(BCtype[X1UP]==OUTFLOWNOINFLOW)){
			/* outer r BC: outflow */

			LOOPN2 LOOPN3{
				ri=N1-1;
				rj=j;
				rk=k;

				//avoid reading from the MPI boundary which is not set and may have NAN's
				if( !finite( prim[ri][rj][rk][0] ) ) continue;

				LOOPBOUND1OUT PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
				if( BCtype[X1UP]==OUTFLOWNOINFLOW ){
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
				}
			}// end 2 3
		}// end if correct bound type
		else if(BCtype[X1UP]== SYMM || BCtype[X1UP]== ASYMM || BCtype[X1UP]== POLARAXIS ) { 
			LOOPN2 LOOPN3 {
				LOOPBOUND1OUT PBOUNDLOOP(pl)  {
					prim[i][j][k][pl] = prim[N1 - 1 + N1 - i][j][k][pl];
					if( BCtype[X1UP]==ASYMM || BCtype[X1UP]== POLARAXIS )  {
						if( pl == U1 || pl == B1 )	prim[i][j][k][pl] *= -1.;
					}
				}
			}
		}
		else if(BCtype[X1UP]== BCEXTRAP) { 
			//for (j = 0; j < N2; j++) {
			//	for(i = N1; i < N1 + N1BND; i++) {
			//		PBOUNDLOOP(pl) pr0[pl] = prim[N1-2][j][k][pl];
			//		PBOUNDLOOP(pl) pr1[pl] = prim[N1-1][j][k][pl];
			//		LOOPBOUND1OUT PBOUNDLOOP(pl) prim[i][j][k][pl]=pr1[pl]*(i-N1+2)+pr0[pl]*(N1-1-i);
			//	}
			//}
			LOOPN2 LOOPN3{
				PBOUNDLOOP(pl) pr0[pl] = prim[N1-1][j][k][pl];
				PBOUNDLOOP(pl) pr1[pl] = prim[N1-2][j][k][pl];
				PBOUNDLOOP(pl) pr2[pl] = prim[N1-3][j][k][pl];
				PBOUNDLOOP(pl) pr3[pl] = prim[N1-4][j][k][pl];
				PBOUNDLOOP(pl) pr4[pl] = prim[N1-5][j][k][pl];
				PBOUNDLOOP(pl) pr5[pl] = prim[N1-6][j][k][pl];
				LOOPBOUND1OUT {
					PBOUNDLOOP(pl) prim[i][j][k][pl] = interpn( 6, i, N1-1, pr0[pl], N1-2, pr1[pl], N1-3, pr2[pl], N1-4, pr3[pl], N1-5, pr4[pl], N1-6, pr5[pl] );
				}
			}
		}
		else if(BCtype[X1UP]== BONDIMDOTOUTFLOW) { 
			LOOPN2 LOOPN3{
				//compute Bondi integrals in active zones
				for( ind = 0; ind < BONDIEXTRAPORDER; ind++ ) {
					bondi_compute_integrals( N1 - 1 - ind, j, k, prim[N1 - 1 - ind][j][k], &BondiC1[ind], &BondiC2[ind], &BondiK[ind] );
				}
				//extrapolate C1 at BONDIEXRAPORDERth order
				LOOPBOUND1OUT {
					//the last two do not matter since order = 5
					BondiC1interp = interpn( BONDIEXTRAPORDER, i, N1-1, BondiC1[0], N1-2,  BondiC1[1], N1-3,  BondiC1[2], N1-4,  BondiC1[3], N1-5, BondiC1[4], N1-6, BondiC1[5]); 
					//outflow Mdot ~ C1 from the grid, use default values for other integrals
					bondi_generate_solution( BondiC1interp, BondiC2x, BondiKx, i, j, k, CENT, prim[i][j][k] );  
				}
			}
		}
		else if(BCtype[X1UP]== RESCALEOUTFLOW) { 
			//FIX all quantities to the initial values; then rescale them to match the closest active grid zone assuming the trend is the same as in initial conditions
			LOOPN2 LOOPN3{
				//generate fixed boundary conditions
				LOOPBOUND1OUT {
					bounds_generate( i, j, k, prim[i][j][k] );
				}
				//relax them to connect to the solution in the active grid cells smoothly
				i = N1-1;
				bounds_generate( i, j, k, pr0 );
				LOOPBOUND1OUT PBOUNDLOOP(pl) {
					prim[i][j][k][pl] += (prim[N1-1][j][k][pl] - pr0[pl]);
				}
			}
		}
		else if(BCtype[X1UP]== RESCALEFIXEDALLOUTFLOWU1) { 
			//FIX all quantities to the initial values; then rescale them to match the closest active grid zone assuming the trend is the same as in initial conditions
			LOOPN2 LOOPN3{
				//generate fixed boundary conditions
				LOOPBOUND1OUT {
					bounds_generate( i, j, k, prim[i][j][k] );
				}
				//relax them to connect to the solution in the active grid cells smoothly
				i = N1-1;
				bounds_generate( i, j, k, pr0 );
				LOOPBOUND1OUT PBOUNDLOOP(pl) {
					if( pl == U1 ) prim[i][j][k][pl] += (prim[N1-1][j][k][pl] - pr0[pl]);
					else prim[i][j][k][pl] -= (prim[N1-1][j][k][pl] + pr0[pl])*(i-N1+1);
				}
			}
		}
		else if(BCtype[X1UP]== FIXED_RESCALEOUTFLOWU1) { 
			//FIX all quantities to the initial values; then rescale them to match the closest active grid zone assuming the trend is the same as in initial conditions
			LOOPN2 LOOPN3{
				//generate fixed boundary conditions
				LOOPBOUND1OUT {
					bounds_generate( i, j, k, prim[i][j][k] );
				}
				//relax them to connect to the solution in the active grid cells smoothly
				i = N1-1;
				bounds_generate( i, j, k, pr0 );
				pl = U1;
				LOOPBOUND1OUT {
					prim[i][j][k][pl] += (prim[N1-1][j][k][pl] - pr0[pl]);
				}
			}
		}
		else if(BCtype[X1UP]== BCU1EXTRAPOTHERFIXED) { 
			LOOPN2 LOOPN3{
				//fix everything except U1
				LOOPBOUND1OUT {
					bounds_generate( i, j, k, prim[i][j][k] );
				}
				//extrapolate U1
				pl = U1;
				pr0[pl] = prim[N1-1][j][k][pl];
				pr1[pl] = prim[N1-2][j][k][pl];
				pr2[pl] = prim[N1-3][j][k][pl];
				pr3[pl] = prim[N1-4][j][k][pl];
				pr4[pl] = prim[N1-5][j][k][pl];
				pr5[pl] = prim[N1-6][j][k][pl];
				LOOPBOUND1OUT {
					prim[i][j][k][pl] = interpn( 6, i, N1-1, pr0[pl], N1-2, pr1[pl], N1-3, pr2[pl], N1-4, pr3[pl], N1-5, pr4[pl], N1-6, pr5[pl] );
				}
			}
		}
		else if(BCtype[X1UP]== BCEXTRAPCONSTRAINED) { 
			LOOPN2 LOOPN3{
				//fix everything at 1st boundary zone
				i = N1; 
				bounds_generate( i, j, k, prim[i][j][k] );
				
				//extrapolate everything through the active zone values (i = N1-5..N1-1) and the 1st boundary zone value (i = N1)
				PBOUNDLOOP(pl)  {
					pr0[pl] = prim[N1-0][j][k][pl];
					pr1[pl] = prim[N1-1][j][k][pl];
					pr2[pl] = prim[N1-2][j][k][pl];
					pr3[pl] = prim[N1-3][j][k][pl];
					pr4[pl] = prim[N1-4][j][k][pl];
					pr5[pl] = prim[N1-5][j][k][pl];
			
					LOOPBOUND1OUT {
						prim[i][j][k][pl] = interpn( 6, i, N1, pr0[pl], N1-1, pr1[pl], N1-2, pr2[pl], N1-3, pr3[pl], N1-4, pr4[pl], N1-5, pr5[pl] );
					}
				}

				//extrapolate velocity (not caring about boundary zones)
				pl = U1;
				pr0[pl] = prim[N1-1][j][k][pl];
				pr1[pl] = prim[N1-2][j][k][pl];
				pr2[pl] = prim[N1-3][j][k][pl];
				pr3[pl] = prim[N1-4][j][k][pl];
				pr4[pl] = prim[N1-5][j][k][pl];
				pr5[pl] = prim[N1-6][j][k][pl];
				LOOPBOUND1OUT {
					prim[i][j][k][pl] = interpn( 6, i, N1-1, pr0[pl], N1-2, pr1[pl], N1-3, pr2[pl], N1-4, pr3[pl], N1-5, pr4[pl], N1-6, pr5[pl] );
				}
			}
		}
		else if(BCtype[X1UP]== BCEXTRAP_VEL3) { 
			//extrapolates velocities after converting them to VEL3
			LOOPN2 LOOPN3{
				PBOUNDLOOP(pl) pr0[pl] = prim[N1-1][j][k][pl];
				PBOUNDLOOP(pl) pr1[pl] = prim[N1-2][j][k][pl];
				PBOUNDLOOP(pl) pr2[pl] = prim[N1-3][j][k][pl];
				PBOUNDLOOP(pl) pr3[pl] = prim[N1-4][j][k][pl];
				PBOUNDLOOP(pl) pr4[pl] = prim[N1-5][j][k][pl];
				PBOUNDLOOP(pl) pr5[pl] = prim[N1-6][j][k][pl];
				metp2met2bl( VEL3, MCOORD, pr0, N1-1, j, k );
				metp2met2bl( VEL3, MCOORD, pr1, N1-2, j, k );
				metp2met2bl( VEL3, MCOORD, pr2, N1-3, j, k );
				metp2met2bl( VEL3, MCOORD, pr3, N1-4, j, k );
				metp2met2bl( VEL3, MCOORD, pr4, N1-5, j, k );
				metp2met2bl( VEL3, MCOORD, pr5, N1-6, j, k );
				//interp3grid_x1( i, j, prim[i][j], N1-1, prim[N1-1][j], N1-2, prim[N1-2][j], N1-3, prim[N1-3][j] );
				//interp3grid_x1( i, j, prim[i][j][k], N1-1, pr0, N1-2, pr1, N1-3, pr2 );
				LOOPBOUND1OUT {
					//bad?PBOUNDLOOP(pl) prim[i][j][k][pl]=pr1[pl]*(i-N1+2)+pr0[pl]*(N1-1-i);
					//PBOUNDLOOP(pl) prim[i][j][k][pl]=- pr1[pl]*(i-N1+1) - pr0[pl]*(N1-2-i);

					PBOUNDLOOP(pl) prim[i][j][k][pl] = interpn( 6, i, N1-1, pr0[pl], N1-2, pr1[pl], N1-3, pr2[pl], N1-4, pr3[pl], N1-5, pr4[pl], N1-6, pr5[pl] );
					bl2met2metp2v( VEL3, MCOORD, prim[i][j][k], i, j, k );
				}
			}
		}
		else if(BCtype[X1UP]== FIXED) { 
			//calls init_dsandvels() to regenerate the solution
			LOOPN2 LOOPN3{
				LOOPBOUND1OUT bounds_generate( i, j, k, prim[i][j][k] );  //SASMARK
			}
		}				


	}// end if mycpu is correct


	// periodic x2
	if ( (mycpupos[2] == 0)&&(mycpupos[2] == ncpux2 - 1) ) {
		if( (BCtype[X2DN]==PERIODIC)&&(BCtype[X2UP]==PERIODIC) ){
			// just copy from one side to another

			LOOPF1 LOOPF3{

				// copy from upper side to lower boundary zones
				ri=i;
				rj=N2;
				rk=k;
				LOOPBOUND2IN PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj+j][rk][pl];

				// copy from lower side to upper boundary zones
				ri=i;
				rj=0;
				rk=k;
				LOOPBOUND2OUT PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj+(j-N2)][rk][pl];
			}
		}
	}

	///////////////////////////
	//
	// X2 inner
	//
	///////////////////////////


	/* inner polar BC (preserves u^t rho and u) */
	if (mycpupos[2] == 0) {
		if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM) ){
			LOOPF1 LOOPN3{
				ri=i;
				rj=0;
				rk=k;

				LOOPBOUND2IN PBOUNDLOOP(pl) {
					prim[i][j][k][pl] = prim[ri][rj+(rj-j-1)][rk][pl];
					if( BCtype[X2DN]==ASYMM || BCtype[X2DN]== POLARAXIS )  {
						if( pl == U2 || pl == B2 )	prim[i][j][k][pl] *= -1.;
					}
				}
			}
		}
		else if( (BCtype[X2DN]==OUTFLOW) || (BCtype[X2DN]==OUTFLOWNOINFLOW) ){
			LOOPF1 LOOPN3{
				ri=i;
				rj=0;
				rk=k;

				//avoid reading from the MPI boundary which is not set and may have NAN's
				if( !finite( prim[ri][rj][rk][0] ) ) continue;

				LOOPBOUND2IN PBOUNDLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj][rk][pl];
				if( BCtype[X2DN]==OUTFLOWNOINFLOW ){
					LOOPBOUND2IN{
#if(WHICHVEL==VEL4)
						get_geometry(i, j, k, CENT, &geom);
						inflow_check_4vel(2,prim[i][j][k],&geom,0) ;
#elif(WHICHVEL==VEL3)
						get_geometry(i, j, k, CENT, &geom);
						inflow_check_3vel(2,prim[i][j][k],&geom,0) ;
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
						inflow_check_rel4vel(2,prim[i][j][k],&geom,0) ;
						if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom, 0)>=1)
							FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif	
					}
				}
			}
		}
		else if(BCtype[X2DN]== JETINJECTION) { 
			//calls init_dsandvels() to regenerate the solution after outflowing the solution first
			LOOPF1 LOOPN3{
				ri=i;
				rj=0;
				rk=k;
				LOOPBOUND2IN PBOUNDLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj][rk][pl];

				//populates the bounds with the jet solution only in the limited area where the jet is injected
				//other grid cells are left to their old values that just got assigned above
				LOOPBOUND2IN bounds_generate( i, j, k, prim[i][j][k] );  //SASMARK
			}
		}				
		else if(BCtype[X2DN]== FIXED) { 
			//calls init_dsandvels() to regenerate the solution
			LOOPF1 LOOPN3{
				LOOPBOUND2IN bounds_generate( i, j, k, prim[i][j][k] );  //SASMARK
			}
		}				

	}// end if mycpupos[2]==0


	///////////////////////////
	//
	// X2 outer
	//
	///////////////////////////


	if (mycpupos[2] == ncpux2-1) {
		if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM) ){
			LOOPF1 LOOPN3{
				ri=i;
				rj=N2-1;
				rk=k;
				LOOPBOUND2OUT PBOUNDLOOP(pl)  {
					prim[i][j][k][pl] = prim[ri][rj+(rj-j+1)][rk][pl];
					if( BCtype[X2UP]==ASYMM || BCtype[X2UP]== POLARAXIS ) {
						if( pl == U2 || pl == B2 )	prim[i][j][k][pl] *= -1.;
					}
				}
			}
		}
		else if( BCtype[X2UP] == OUTFLOW || BCtype[X2UP] == OUTFLOWNOINFLOW ){
			LOOPF1 LOOPN3{
				ri=i;
				rj=N2-1;
				rk=k;

				//avoid reading from the MPI boundary which is not set and may have NAN's
				if( !finite( prim[ri][rj][rk][0] ) ) continue;

				LOOPBOUND2OUT PBOUNDLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj][rk][pl];
				if( BCtype[X2UP]==OUTFLOWNOINFLOW ){
					LOOPBOUND2OUT{
#if(WHICHVEL==VEL4)
						get_geometry(i, j, k, CENT, &geom);
						inflow_check_4vel(2,prim[i][j][k],&geom,0) ;
#elif(WHICHVEL==VEL3)
						get_geometry(i, j, k, CENT, &geom);
						inflow_check_3vel(2,prim[i][j][k],&geom,0) ;
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
						inflow_check_rel4vel(2,prim[i][j][k],&geom,0) ;
						if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom, 0)>=1)
							FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif	
					}
				}
			}
		}
		else if(BCtype[X2UP]== FIXED) { 
			//calls init_dsandvels() to regenerate the solution
			LOOPF1 LOOPN3{
				LOOPBOUND2OUT bounds_generate( i, j, k, prim[i][j][k] );  //SASMARK
			}
		}				

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


	///////////////////////////
	//
	// X3 inner 
	//
	///////////////////////////


	/* inner polar BC (preserves u^t rho and u) */
	if (mycpupos[3] == 0) {
		if((BCtype[X3DN]==POLARAXIS)||(BCtype[X3DN]==SYMM)||(BCtype[X3DN]==ASYMM) ){
			LOOPF1 LOOPN2{
				ri=i;
				rj=j;
				rk=0;
				LOOPBOUND3IN PBOUNDLOOP(pl) {
					prim[i][j][k][pl] = prim[ri][rj][rk+(rk-k-1)][pl];
					if( BCtype[X3DN]==ASYMM || BCtype[X3DN]== POLARAXIS )  {
						if( pl == U3 || pl == B3 )	prim[i][j][k][pl] *= -1.;
					}
				}
			}
		}
		else if( BCtype[X3DN]==OUTFLOW || BCtype[X3DN] == OUTFLOWNOINFLOW ){
			LOOPF1 LOOPN2{
				ri=i;
				rj=j;
				rk=0;

				//avoid reading from the MPI boundary which is not set and may have NAN's
				if( !finite( prim[ri][rj][rk][0] ) ) continue;

				LOOPBOUND3IN PBOUNDLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj][rk][pl];
				if( BCtype[X3DN]==OUTFLOWNOINFLOW ){
					LOOPBOUND3IN{
#if(WHICHVEL==VEL4)
						get_geometry(i, j, k, CENT, &geom);
						inflow_check_4vel(3,prim[i][j][k],&geom,0) ;
#elif(WHICHVEL==VEL3)
						get_geometry(i, j, k, CENT, &geom);
						inflow_check_3vel(3,prim[i][j][k],&geom,0) ;
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
						inflow_check_rel4vel(3,prim[i][j][k],&geom,0) ;
						if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom, 0)>=1)
							FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif	
					}
				}
			}
		}

	}// end if mycpupos[3]==0


	///////////////////////////
	//
	// X3 outer
	//
	///////////////////////////


	if (mycpupos[3] == ncpux3-1) {
		if((BCtype[X3UP]==POLARAXIS)||(BCtype[X3UP]==SYMM)||(BCtype[X3UP]==ASYMM) ){
			LOOPF1 LOOPN2{
				ri=i;
				rj=j;
				rk=N3-1;
				LOOPBOUND3OUT PBOUNDLOOP(pl)  {
					prim[i][j][k][pl] = prim[ri][rj][rk+(rk-k+1)][pl];
					if( BCtype[X3UP]==ASYMM || BCtype[X3UP]== POLARAXIS ) {
						if( pl == U3 || pl == B3 )	prim[i][j][k][pl] *= -1.;
					}
				}
			}
		}
		else if( BCtype[X3UP] == OUTFLOW || BCtype[X3UP] == OUTFLOWNOINFLOW ){
			LOOPF1 LOOPN2{
				ri=i;
				rj=j;
				rk=N3-1;

				//avoid reading from the MPI boundary which is not set and may have NAN's
				if( !finite( prim[ri][rj][rk][0] ) ) continue;

				LOOPBOUND3OUT PBOUNDLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj][rk][pl];
				if( BCtype[X3UP]==OUTFLOWNOINFLOW ){
					LOOPBOUND3OUT{
#if(WHICHVEL==VEL4)
						get_geometry(i, j, k, CENT, &geom);
						inflow_check_4vel(3,prim[i][j][k],&geom,0) ;
#elif(WHICHVEL==VEL3)
						get_geometry(i, j, k, CENT, &geom);
						inflow_check_3vel(3,prim[i][j][k],&geom,0) ;
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
						inflow_check_rel4vel(3,prim[i][j][k],&geom,0) ;
						if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom, 0)>=1)
							FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif	
					}
				}
			}

		}

	}// end if mycpupos[3]==ncpux2-1


#if(ASYMDIAGCHECK)
	dualfprintf(fail_file,"1after bound\n");
	asym_compute_1(prim);
	dualfprintf(fail_file,"2after bound\n");
#endif



	return (0);
}

int	bondi_compute_integrals( int i0, int j0, int k0, FTYPE *primin, FTYPE *C1out, FTYPE *C2out, FTYPE *Kout )
{
	FTYPE C1, C2, K;
  FTYPE X[NDIM],V[NDIM];
	FTYPE r;
	FTYPE ucon[NDIM];
	FTYPE n = 1 / (gam - 1);
	FTYPE p;
	FTYPE rho;
	FTYPE uu1;
	FTYPE T;
	int pl;
  struct of_geom geom;
	FTYPE prbl[NPR];


	/////////
	//
	// compute the value fo C1 based on primin[] and compute the values of primout[] that correspond to (C1, C2x):
	//

	//convert primin[] form prime coords to Boyer-Lindquist prbl[]
	PLOOP(pl) prbl[pl] = primin[pl];
	//convert primitives from WHICHVEL MCOORD to VEL4 KSCOORD
	metp2met2bl( VEL4, KSCOORDS, prbl, i0, j0, k0 );

	//compute BL coordinates of primin[]
	coord(i0, j0, k0, CENT, X);  //get coordinates of (i0, j0, k0) point in prime (internal) coordinate system, X
  bl_coord(X, V);              //get these coords in the orthonormal basis, V

	r = V[1];
	p = (gam-1) * prbl[UU];
	rho = prbl[RHO];
	uu1 = prbl[U1];
	T = p / rho;      //"temperature" of the flow
	
	//compute the values of integrals
	C1 = pow( T, n ) * uu1 * r * r;  //C1 = rho u^r r^2
	C2 = pow( 1 + (1+n) * T, 2. ) * ( 1 - 2. / r + uu1 * uu1 );
	K  = p / pow( rho, gam );

	if( C1out ) *C1out = C1;
	if( C2out ) *C2out = C2;
	if( Kout ) *Kout = K;

	return( 0 );	
}


//generates the bondi solution at a given location (i0, j0, k0) at a given (face) in MCOORD
int	bondi_generate_solution( FTYPE C1, FTYPE C2, FTYPE K, int i1, int j1, int k1, int face, FTYPE *primout )
{
	//exact values of integrals for Bondi problem, taken from bondigrid4acc.nb
	FTYPE C1x = -0.00675;
	FTYPE C2x = 1.373125;
  FTYPE X[NDIM],V[NDIM];
	FTYPE ucon[NDIM];
	FTYPE n = 1 / (gam - 1);
	int pl;
	FTYPE relative_accuracy = 1.0e-15;
	int inversion_flag;
  struct of_geom geom;
	FTYPE prbl[NPR];

	int bondi_invert_integrals( FTYPE n, FTYPE C1, FTYPE C2, FTYPE K, FTYPE relative_accuracy, FTYPE r, FTYPE *primout );
	void mettometpface(int ii, int jj, int kk, int face, FTYPE*ucon);
	void gset_face(int getprim, int whichcoord, int i, int j, int k, int face, struct of_geom *ptrgeom);

	//compute BL coordinates of primout[]
	coord(i1, j1, k1, face, X);  //get coordinates of (i0, j0, k0) point in prime (internal) coordinate system, X
  bl_coord(X, V);              //get these coords in the orthonormal basis, V

	inversion_flag = bondi_invert_integrals( n, C1, C2, K, relative_accuracy, V[1], primout );
	if( inversion_flag ) { //inversion error occured
		myexit( inversion_flag );
	}

	ucon[0]=0;
  ucon[1]=primout[U1];
  ucon[2]=primout[U2];
  ucon[3]=primout[U3];

	//convert primitives from VEL4 KSCOORD to WHICHVEL MCOORD (internal coord system)
	//bl2met2metp2v( VEL4, KSCOORDS, primout, i0, j0, k0 );

	
  // pr is in whichcoord coordinates
  // get geometry (non-prime coords)
  gset_face(0,KSCOORDS,i1,j1,k1,face,&geom);
  // convert whichvel-pr in whichcoord coords to ucon in whichcoord coordinates
  if (pr2ucon(VEL4,primout, &geom, ucon) >= 1) FAILSTATEMENT("init.c:init()", "pr2ucon()", 1);

	//converts velocity from MCOORD non-prime to MCOORD prime (since time is not involved, no need for ucon[0] to be defined; not general SASMARKx
	mettometpface( i1, j1, k1, face, ucon );

  // get prime geometry
  get_geometry( i1, j1, k1, face, &geom ) ;
  // convert from MCOORD prime 4-vel to MCOORD prime WHICHVEL-vel(i.e. primitive velocity of evolution)
  ucon2pr(WHICHVEL,ucon,&geom,primout);

	return( 0 );
}

// find the con/cov forms of the chosen metric
void gset_face(int getprim, int whichcoord, int i, int j, int k, int face, struct of_geom *ptrgeom)
{
  FTYPE X[NDIM];
  struct of_geom tempgeom;
  extern void assign_eomfunc(struct of_geom *geom, FTYPE eomfuncgen);
  extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  extern void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
	extern void eomfunc_func(int getprim, int whichcoord, FTYPE *X, FTYPE *eomfunc);

  FTYPE eomfuncgen;

  ptrgeom->i=i;
  ptrgeom->j=j;
  ptrgeom->k=k;
  ptrgeom->p=face;
  icurr=i;
  jcurr=j;
  jcurr=k;
  pcurr=CENT;
  
  if(whichcoord>=0){
    coord(i, j, k, face, X);
    gcov_func(getprim,whichcoord,X,gengcov,gengcovpert);
    ptrgeom->g=gdet_func(gengcov); // must come after gcov_func() above
    gcon_func(getprim,whichcoord,X,gengcov,gengcon); // must come after gcov_func() above
    
    eomfunc_func(getprim, whichcoord,X,&eomfuncgen);
    assign_eomfunc(ptrgeom,eomfuncgen); // must come after assigning ptrgeom->g above
    
    ptrgeom->gcov=gengcov; // pointer
    ptrgeom->gcon=gengcon; // pointer
    ptrgeom->gcovpert=gengcovpert; // pointer
  }
  else if(whichcoord==PRIMECOORDS){ // special case
    get_geometry(i,j,k,face,ptrgeom);
  }
  else{
    dualfprintf(fail_file,"gset(): no such whichcoord=%d\n",whichcoord);
    myexit(1);
  }

}

// MCOORD -> prime MCOORD
void mettometpface(int ii, int jj, int kk, int face, FTYPE*ucon)
{
  int j,k;
  FTYPE X[NDIM], V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE tmp[NDIM];

  coord(ii, jj, kk, face, X);
  bl_coord(X, V);
  //  r=V[1]; th=V[2];

  dxdxprim(X, V, dxdxp);
  // actually gcon_func() takes inverse of first arg and puts result into second arg.
  matrix_inverse(dxdxp,idxdxp);

  /* transform ucon */
  // this is u^j = (iT)^j_k u^k, as in mettobl() above
  DLOOPA(j) tmp[j] = 0.;
  DLOOP(j,k) tmp[j] += idxdxp[j][k] * ucon[k];
  DLOOPA(j) ucon[j] = tmp[j];
  
  // note that u_{k,BL} = u_{j,KSP} (iT)^j_k  

  // note that u_{k,KSP} = u_{j,BL} T^j_k  

  // note that u^{j,BL} =  T^j_k u^{k,KSP}   // (T) called ks2bl in grmhd-transforms.nb

  // note that u^{j,KSP} = (iT)^j_k u^{k,BL} // (iT) called bl2ks in grmhd-transforms.nb

  // So T=BL2KSP for covariant components and KSP2BL for contravariant components
  // and (iT)=BL2KSP for contra and KSP2BL for cov

  // where here T=dxdxp and (iT)=idxdxp (not transposed, just inverse)

  /* done! */
}
//outflow Mdot from the grid, use default values for other integrals
//int	bondi_compute_mdot_outflow( int i0, int j0, int k0, FTYPE *primin, int i1, int j1, int k1, FTYPE *primout )
//{
//	//exact values of integrals for Bondi problem, taken from bondigrid4acc.nb
//	FTYPE C1x = -0.00675;
//	FTYPE C2x = 1.373125;
//	FTYPE C1;         //computed values of C1 (based on primin[])
//	FTYPE prbl[NPR];  //primitives in Boyer-Lindquist
//  FTYPE X[NDIM],V[NDIM];
//	FTYPE n = 1 / (gam - 1);
//	int pl;
//	FTYPE relative_accuracy = 1.0e-14;
//	int inversion_flag;
//	int bondi_invert_integrals( FTYPE n, FTYPE C1, FTYPE C2, FTYPE K, FTYPE relative_accuracy, FTYPE r, FTYPE *primout );
//
//	/////////
//	//
//	// compute the value fo C1 based on primin[] and compute the values of primout[] that correspond to (C1, C2x):
//	//
//
//	//convert primin[] form prime coords to Boyer-Lindquist prbl[]
//	PLOOP(pl) prbl[pl] = primin[pl];
//	//convert primitives from WHICHVEL MCOORD to VEL4 KSCOORD
//	metp2met2bl( VEL4, KSCOORDS, prbl, i0, j0, k0 );
//
//	//compute BL coordinates of primin[]
//	coord(i0, j0, k0, CENT, X);  //get coordinates of (i0, j0, k0) point in prime (internal) coordinate system, X
//  bl_coord(X, V);              //get these coords in the orthonormal basis, V
//
//
//	//compue the value of C1
//	C1 = prbl[RHO] * prbl[U1] * V[1] * V[1];  //C1 = rho u^r r^2
//
//	//if( C1 < 1.7 * C1x ) {  //no solutions if the value of C1 is too large
//	//	trifprintf( "C1 = %21.15lg > 1.7 C1x, therefore no solution will be found\n", C1 );
//	//	myexit( 1 );
//	//}
//
//	//compute BL coordinates of primout[]
//	coord(i1, j1, k1, CENT, X);  //get coordinates of (i1, j1, k1) point in prime (internal) coordinate system, X
//  bl_coord(X, V);              //get these coords in the orthonormal basis, V
//
//	inversion_flag = bondi_invert_integrals( n, C1, C2x, 1, relative_accuracy, V[1], primout );
//
//	if( inversion_flag ) { //inversion error occured
//		myexit( inversion_flag );
//	}
//
//	//convert primitives back to WHICHVEL MCOORD from VEL4 KSCOORD
//	bl2met2metp2v( VEL4, KSCOORDS, primout, i1, j1, k1 );
//
//	return( 0 );
//}

//finds the values of primtive quantities that correspond to the integrals; returns the primitives in KSCOORDS
int bondi_invert_integrals( FTYPE n, FTYPE C1, FTYPE C2, FTYPE K, FTYPE relative_accuracy, FTYPE r, FTYPE *primout )
{
	int pl;
	FTYPE T0;
	FTYPE T1;
	FTYPE u1;
	FTYPE rho;
	int iter;
	int max_iter = 250;


	//initialize Told = T0 and Tnew = T1 with values that are too different so that the algorithm picks up and starts iterating with T = T1
	T0 = 0.2;
	T1 = 0.1;

	if( r > 0.35 && r < 2.5 ){
		//this iteration form works well for 0.6 < r < 3
		//iterate to solve T = f(T) until the difference between Tnew and Told is smaller than relative_accuracy
		for( iter = 0; iter < max_iter && fabs((T0 - T1)/T1) > relative_accuracy; iter++ ) {
			T0 = T1;
			u1 = C1 / ( r * r * pow(T0, n) );  //the approximation to radial velocity
			T1 = T0 + 0.02 * ( T0 - ( sqrt(C2 / (1 - 2 / r + u1 * u1)) - 1 ) / ( 1 + n ) );  //new approximation to temperature, which is = p/rho
		}
	}
	else if( r > 9 ) {
		//this iteration form works well for 9 < r < 30
		//iterate to solve T = f(T) until the difference between Tnew and Told is smaller than relative_accuracy
		for( iter = 0; iter < max_iter && (T0 - T1)/T1 > relative_accuracy; iter++ ) {
			T0 = T1;
			u1 = C1 / ( r * r * pow(T0, n) );  //the approximation to radial velocity
			T1 = ( sqrt(C2 / (1 - 2 / r + u1 * u1)) - 1 ) / ( 1 + n );  //new approximation to temperature, which is = p/rho
		}

		if( T1 > T0 ) {
			trifprintf( "Inversion of Bondi integrals failed: root finding procedure expects temperature to be monotonically decreasing as it approaches the root;\nIt is not: Tnew = %21.15g > Told = %21.15g", T1, T0 );
			return( 1 );
		}
	}
	else {
		trifprintf( "Bad r = %lg, supported range: 0.6 < r < 3 and 9 < r < 30\n", r );
		return( 1 );
	}


	if( iter > max_iter ) {
		trifprintf( "Inversion of Bondi integrals failed: maximum number of iterations reached iter = %d = maxiter; current accuracy = (T0-T1)/T1 = %21.15lg; desired accuracy = %21.15lg\n", 
			iter, (T0-T1)/T1, relative_accuracy );
		return( 2 );
	}


	// compute other primitive quantities based on termperature T1
  PLOOP(pl) primout[pl] = 0.0;  //zero out the output array

	rho = pow( T1/K, n );									       //rho = T^n / K^n
	primout[RHO] = rho;
	primout[UU]  = rho * T1 * n;					       //ie = rho T / (gam - 1) = rho T n; n = 1/(gam - 1); T = p / rho = (gam-1) u / rho.
	primout[U1] =  C1 / ( pow(T1, n) * r * r );  //u^r = C1 / ( T^n * r^2 );

	return( 0 );  //ok
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

int bound_prim_user_after_mpi(int boundstage, FTYPE prim[][N2M][N3M][NPR])
{
   return( 0 );
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
