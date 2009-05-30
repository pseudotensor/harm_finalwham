
// GODMARK: Should redo flux's so that fluxes are accessed by F1[j][k][i] F2[k][i][j] F3[i][j][k] for faster differencing in advance.c

#include "decs.h"

// see fluxcompute.c for non-computer science, real physics calculations of flux
int fluxcalc(int stage, FTYPE pr[][N2M][N3M][NPR],
	     FTYPE F1[][N2M][N3M][NPR], 
	     FTYPE F2[][N2M][N3M][NPR], 
	     FTYPE F3[][N2M][N3M][NPR], 
	     FTYPE CUf,
	     FTYPE *ndt1,
	     FTYPE *ndt2,
	     FTYPE *ndt3
	     )
{
  int fluxcalc_1d(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt);
  extern int flux_ct(int stage, FTYPE pr[][N2M][N3M][NPR],FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]);
  void fix_flux(FTYPE (*pr)[N2M][N3M][NPR],FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]) ;
  extern void flux_interp(int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prims_guess)[N2M][N3M][NPR], FTYPE (*p2interpm)[N2M][N3M][NPR], FTYPE (*p2interpp)[N2M][N3M][NPR], FTYPE (*pleft)[N2M][N3M][NPR], FTYPE (*pright)[N2M][N3M][NPR]);
  extern int bound_flux(int boundstage, FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]);
  int dir;
  int idel, jdel, kdel, face, idel1, idel2, jdel1, jdel2, kdel1, kdel2;
  int is, ie, js, je, ks, ke;
  FTYPE (*dqvec[NDIM])[N2M][N3M][NPR2INTERP];
  FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR];
  FTYPE (*fluxveca[NDIM])[N2M][N3M][NPR];
  FTYPE (*fluxvecb[NDIM])[N2M][N3M][NPR];
  FTYPE *ndtvec[NDIM];
  int pl,i,j,k;
  int Nvec[NDIM],avgdirtype[NDIM],odir1,odir2;
	FTYPE limit_fluxc2a_prim_change( 
			 int dir, FTYPE pr[][N2M][N3M][NPR],
	     FTYPE fluxvec_point[][N2M][N3M][NPR],
			 FTYPE fluxvec_avg[][N2M][N3M][NPR]);   //atch


  /////////////////////////
  //
  // SETUP dimensionality
  //
  /////////////////////////

  fluxvec[1]=F1;
  fluxvec[2]=F2;
  fluxvec[3]=F3;

  // storage for UNSPLIT flux_interp method for FV method
  fluxveca[1]=Fa;
  fluxveca[2]=Fa;
  fluxveca[3]=Fa;
  fluxvecb[1]=Fb;
  fluxvecb[2]=Fb;
  fluxvecb[3]=Fb;

  dqvec[1]=dq1;
  dqvec[2]=dq2;
  dqvec[3]=dq3;

  ndtvec[1]=ndt1;
  ndtvec[2]=ndt2;
  ndtvec[3]=ndt3;

  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  avgdirtype[1]=ENOFLUXAVG1TYPE;
  avgdirtype[2]=ENOFLUXAVG2TYPE;
  avgdirtype[3]=ENOFLUXAVG3TYPE;

	#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )  //SUPERSASMARK atch
		//zero out the array with lower order fractions at every substep before c2e reconstructions are performed
		FULLLOOP {
		 DIMENLOOP(dir)	weno_prim_lower_order_fraction[dir][i][j][k] = 0.0;
		}
	#endif

  DIMENLOOP(dir){

    // set dimension having no influence on dt by default
    *(ndtvec[dir])=BIG;

    // skip to next dir if no such dimension
    if(Nvec[dir]==1) continue;

    // get loop details
    idel = fluxloop[dir][FIDEL];
    jdel = fluxloop[dir][FJDEL];
    kdel = fluxloop[dir][FKDEL];
    face = fluxloop[dir][FFACE];
    is=fluxloop[dir][FIS];  //loop over the interfaces where fluxes are computed -- atch, use ZSLOOP( is, ie, js, je, ks, ke ) { ... }
    ie=fluxloop[dir][FIE];
    js=fluxloop[dir][FJS];
    je=fluxloop[dir][FJE];
    ks=fluxloop[dir][FKS];
    ke=fluxloop[dir][FKE];


    MYFUN(fluxcalc_1d(stage, pr, dir, is, ie, js, je, ks, ke, idel, jdel, kdel, face, dqvec[dir], fluxvec[dir], CUf, ndtvec[dir]),"flux.c:fluxcalc()", "fluxcalc_1d", dir);

#if(PRODUCTION==0)
    trifprintf("%d",dir);
#endif
  }// end DIMENLOOP(dir)


  //////////////////////////////
  //
  // FIXFLUX
  //
  /////////////////////////////

  if(FIXUPFLUX){
    fix_flux(pr, F1, F2, F3);
#if(PRODUCTION==0)
    trifprintf( "x");
#endif
  }





#if(PRODUCTION==0)
  trifprintf( "c");
#endif

  //////////////////////////////
  //
  // ENOFY FLUXES
  //
  /////////////////////////////

  if((DOENOFLUX==ENOFLUXRECON)&&(avgscheme>3)){
    // convert quasi-averaged fluxes to quasi-point fluxes
    
    // bound the fluxes (only needed by this method)
    bound_flux(-1,F1,F2,F3);

#if(0)
    // GODMARK : test of bound_flux()
    if(nstep==30){
      // test of bound_flux().  Appears to work correctly.
      //
      // read in SM:
      // first add header from dump file and make size of grid as FULLLOOP zones
      //
      // jrdpheader3d 0_fail.out
      // da dumps/0_fail.out lines 2 1000000
      // read '%d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g' {i j k f10 f20 f11 f21 f12 f22 f13 f23 f14 f24 f15 f25 f16 f26 f17 f27}
      // set r=4+i set h=j set ph=k gsetup gammienew
      //
      FULLLOOP{
	fprintf(fail_file,"%d %d %d ",i,j,k);
	PLOOP(pl) fprintf(fail_file,"%g %g ",F1[i][j][k][pl],F2[i][j][k][pl]);
	fprintf(fail_file,"\n");
      }
      //		myexit(0);
    }
#endif

    // F1 input and F1 as output is ok since each full 1D line is operated on, and each 1D line never used information from other lines, and only 1 direction per F1/F2/F3, so completely done with that dataset
#if(N1>1)
    dir=1;
    idel = fluxloop[dir][FIDEL];
    jdel = fluxloop[dir][FJDEL];
    kdel = fluxloop[dir][FKDEL];
    flux_interp(ENOFLUXRECONTYPE, dir, idel, jdel, kdel, pr, F1, NULL, F1, NULL); 
#endif
#if(N2>1)
    dir=2;
    idel = fluxloop[dir][FIDEL];
    jdel = fluxloop[dir][FJDEL];
    kdel = fluxloop[dir][FKDEL];
    flux_interp(ENOFLUXRECONTYPE, dir, idel, jdel, kdel, pr, F2, NULL, F2, NULL);
#endif
#if(N3>1)
    dir=3;
    idel = fluxloop[dir][FIDEL];
    jdel = fluxloop[dir][FJDEL];
    kdel = fluxloop[dir][FKDEL];
    flux_interp(ENOFLUXRECONTYPE, dir, idel, jdel, kdel, pr, F3, NULL, F3, NULL);
#endif


#if(0)
    // GODMARK : test of bound_flux()
    if(nstep==30){
      // test of bound_flux().  Appears to work correctly.
      //
      // read in SM:
      // first add header from dump file and make size of grid as FULLLOOP zones
      //
      // jrdpheader3d 0_fail.out
      // da dumps/0_fail.out lines 2 1000000
      // read '%d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g' {i j k f10 f20 f11 f21 f12 f22 f13 f23 f14 f24 f15 f25 f16 f26 f17 f27}
      // set r=4+i set h=j set ph=k gsetup gammienew
      //
      FULLLOOP{
	fprintf(log_file,"%d %d %d ",i,j,k);
	PLOOP(pl) fprintf(log_file,"%g %g ",F1[i][j][k][pl],F2[i][j][k][pl]);
	fprintf(log_file,"\n");
      }
      myexit(0);
    }
#endif
#if(PRODUCTION==0)
    trifprintf( "e");
#endif
    
    
  }
  else if((DOENOFLUX==ENOFLUXSPLIT)&&(avgscheme>3)){
    // GODMARK: NOTICE: this will probably go inside fluxcalc_fluxspliteno() because the procedure is
    
    // -1) no need to bound fluxes since SHOULD exist everywhere (should be taken care of by new fluxcalc_fluxspliteno() function)
    // 0) assume primitives everywhere defined (all grid including ghost zones)
    // 1) Compute split F's (F^+ F^-) for all grid including ghost zones, store in F1m/p F2m/p F3m/p
    // 2) use flux_interp(ENOFLUXSPLITTYPE,...) for each dir and get F1l/r F2l/r F3l/r
    // 3) add together fluxes and store in F1,F2,F3 (can be same memory space as used for F1m/p F2m/p F3m/p since done with it, but can't be same memory space as F123l/r
    
    //    dir=1; flux_interp(ENOFLUXSPLITTYPE, dir, idel, jdel, kdel, F1m, F1p, F1l, F1r);
    //    dir=2; flux_interp(ENOFLUXSPLITTYPE, dir, idel, jdel, kdel, F2m, F2p, F2l, F2r);
    //    dir=3; flux_interp(ENOFLUXSPLITTYPE, dir, idel, jdel, kdel, F3m, F3p, F3l, F3r);
    // these should be read in correct form as Fr Fl in fluxcalc_fluxspliteno() such that for i-th index:
    // dir=1
    // Fr(i) = F1l(i)
    // Fl(i) = F1r(i-1)
  }
  else if((DOENOFLUX==ENOFINITEVOLUME)&&(avgscheme>3)&&do_transverse_flux_integration){
    // integrate point fluxes to surface averaged fluxes

    DIMENLOOP(dir){
      
      // skip to next dir if no such dimension
      if(Nvec[dir]==1) continue;

			#if( LIMIT_FLUXC2A_PRIM_CHANGE ) //atch
			  //otherwise fluxvectemp points to space in memory set in set_arrays.c
				//intitialize it with the flux given by the riemann solver in the dir direction
				FULLLOOP PLOOP(pl) fluxvectemp[i][j][k][pl] = fluxvec[dir][i][j][k][pl];
			#endif

			// other dimensions
      odir1=dir%3+1;
      odir2=(dir+1)%3+1;

      // get loop details
      idel1 = fluxloop[odir1][FIDEL];
      jdel1 = fluxloop[odir1][FJDEL];
      kdel1 = fluxloop[odir1][FKDEL];

      idel2 = fluxloop[odir2][FIDEL];
      jdel2 = fluxloop[odir2][FJDEL];
      kdel2 = fluxloop[odir2][FKDEL];

      // do quasi-Strang splitting if chosen or if <3D since in 1D or 2D there is no ordering issue
      if( (FLUXDIMENSPLIT==QUASISTRANG)||(!(N1NOT1*N2NOT1*N3NOT1)) ){

				// Quasi-Strang Splitting (alternate order of dimension) of ordering of c2a integration of flux over 2 orthogonal directions
				if((nstep*(long long)numstepparts+(long long)steppart)%2){
					// m%3+1 gives next 1->2,2->3,3->1
					// 3-(4-m)%3 = (dir+1)%3+1 gives previous 1->3,2->1,3->2
					if(Nvec[odir1]>1) flux_interp(avgdirtype[odir1], dir, idel1, jdel1, kdel1, pr, fluxvec[dir], NULL, fluxvec[dir], NULL);
					if(Nvec[odir2]>1) flux_interp(avgdirtype[odir2], dir, idel2, jdel2, kdel2, pr, fluxvec[dir], NULL, fluxvec[dir], NULL);
				}
				else{
					if(Nvec[odir2]>1) flux_interp(avgdirtype[odir2], dir, idel2, jdel2, kdel2, pr, fluxvec[dir], NULL, fluxvec[dir], NULL);
					if(Nvec[odir1]>1) flux_interp(avgdirtype[odir1], dir, idel1, jdel1, kdel1, pr, fluxvec[dir], NULL, fluxvec[dir], NULL);
				}
			}
			else if(FLUXDIMENSPLIT==UNSPLIT){
				// only relevant for 3D, and Fa and Fb are only used in this case
				if(Nvec[odir1]>1) flux_interp(avgdirtype[odir1], dir, idel1, jdel1, kdel1, pr, fluxvec[dir], NULL, fluxveca[dir], NULL);
				if(Nvec[odir2]>1) flux_interp(avgdirtype[odir2], dir, idel2, jdel2, kdel2, pr, fluxveca[dir], NULL, fluxveca[dir], NULL);

				if(Nvec[odir2]>1) flux_interp(avgdirtype[odir2], dir, idel2, jdel2, kdel2, pr, fluxvec[dir], NULL, fluxvecb[dir], NULL);
				if(Nvec[odir1]>1) flux_interp(avgdirtype[odir1], dir, idel1, jdel1, kdel1, pr, fluxvecb[dir], NULL, fluxvecb[dir], NULL);

				// symmetrized flux
				FULLLOOP PLOOP(pl){
					fluxvec[dir][i][j][k][pl]=0.5*(fluxveca[dir][i][j][k][pl]+fluxvecb[dir][i][j][k][pl]);
				}
			}

			#if( LIMIT_FLUXC2A_PRIM_CHANGE )   //atch
				//limits how different fluxvec[dir] (averaged fluxes) is from fluxvectemp (point fluxes) are based on how much rho and gamma would
			  //change due to such a difference in flux during the timestep, dt
				limit_fluxc2a_prim_change( dir, pr, fluxvectemp, fluxvec[dir] );
			#endif
		} // end over DIMENLOOP



#if(PRODUCTION==0)
    trifprintf( "e");
#endif
    
    
  }



  //////////////////////////////
  //
  // FLUXCT : must come after modifying fluxes so that divb=0 is preserved to machine error
  // i.e. flux_ct modifies fluxes in special way
  //
  /////////////////////////////

  MYFUN(flux_ct(stage, pr, F1, F2, F3),"step_ch.c:advance()", "flux_ct",1);


  return(0);
  
  
}


FTYPE limit_fluxc2a_prim_change( 
			 int dir, FTYPE pr[][N2M][N3M][NPR],
	     FTYPE fluxvec_point[][N2M][N3M][NPR],
			 FTYPE fluxvec_avg[][N2M][N3M][NPR])
{
	int is, ie, js, je, ks, ke;
	int i, j, k;
	int i1, j1, k1;
	int pl;
	int index;
	int idel, jdel, kdel;
	FTYPE delta_u_left[NPR];
	FTYPE delta_u_right[NPR];
	FTYPE delta_u[2][NPR];
	FTYPE Upoint[NPR], Upoint_updated[NPR];
	FTYPE pr_updated[NPR];
	int pflag_current;
	int pflag_backup;
	FTYPE frac_point_flux_used_array[2];
	FTYPE frac_point_flux_used;
	struct of_geom geom; //atch
	struct of_state q; //atch
	extern FTYPE limit_prim_correction( FTYPE fractional_difference_threshold, struct of_geom *geom, FTYPE *pin, FTYPE *pout );


	//unit vector in the flux direction
	idel = (dir == 1);
	jdel = (dir == 2);
	kdel = (dir == 3);

	//start and end indices for the loop over the fluxes in the dir direction that will be used for evolving the conserved quantities
	//since the conserved quantites are confined to the range of indices given by Uconsloop and there is one more flux
	//that conserved quantities, we expand the loop by one grid cell in the direction of the flux
	is = Uconsloop[FIS]; 
	ie = Uconsloop[FIE] + idel;
	js = Uconsloop[FJS];
	je = Uconsloop[FJE] + jdel;
	ks = Uconsloop[FKS];
	ke = Uconsloop[FKE] + kdel;

	ZSLOOP( is, ie, js, je, ks, ke ) {
		PLOOP( pl ) { //loop over the interfaces where the fluxes that determine the evolution of conserved quantities are
			//
			// !! Insert c2a limiting code here !! SASMARK atch
			//

			//
			// compute the changes that would be made to conserved quantities due to the c2a corrections for the fluxes at all interfaces during dt
			// fluxvec_avg contains the averaged fluxes
			//


			//pseudo-update to conserved quantity located to the left of the current interface, i.e. located at (i - idel, j - jdel, k - kdel),
			//due to c2a correction to the flux at (i,j,k) interface in the dir direction 
			delta_u_left[pl] = - dt * (fluxvec_avg[i][j][k][pl] - fluxvec_point[i][j][k][pl]) / dx[dir];

			//to the right:
			delta_u_right[pl] = - delta_u_left[pl];

			//unite them in one array
			delta_u[0][pl] = delta_u_left[pl];
			delta_u[1][pl] = delta_u_right[pl];
		}

		//loop through the two neighbouring cells that are affected by the fluxvec_avg[i][j][k];
		//index = 0 is the left and index = 1 is the right cells
		for( index = 0; index <= 1; index++ ) {

			i1 = i + (index-1) * idel;
			j1 = j + (index-1) * jdel;
			k1 = k + (index-1) * kdel;

			//1. get u_point_adjacent[0..1] from pr's -- inconsistent but this way avoid bounding

			// set geometry for zone to be updated
			get_geometry(i1, j1, k1, CENT, &geom);

			// find U(pr)
			MYFUN(get_state(pr[i1][j1][k1], &geom, &q),"flux.c:fluxcalc()", "get_state()", 1);
			MYFUN(primtoU(UEVOLVE,pr[i1][j1][k1], &q, &geom, Upoint),"step_ch.c:advance()", "primtoU()", 1);


			//2. add delta_u[0..1] to it
			//that's what Upoint would become were it evolve only due to c2a flux correction on one interface
			PLOOP(pl) Upoint_updated[pl] = Upoint[pl] + delta_u[index][pl]; 
			PLOOP(pl) pr_updated[pl] = pr[i1][j1][k1][pl];  //fill in the initial guess for inversion

			//3. invert & check if the difference between pr_adjacent_updated[0..1] and pr_adjacent_original[0..1]= pr is too large
			// invert point Upoint_updated-> point pr_updated
			pflag_backup = pflag[i1][j1][k1][FLAGUTOPRIMFAIL]; //back up the old inversion flag, just in case, probably not needed anyway
			
			MYFUN(Utoprimgen(EVOLVEUTOPRIM, UEVOLVE, Upoint_updated, &geom, pr_updated),"flux.c:fluxcalc()", "Utoprimgen", 1);

			pflag_current = pflag[i1][j1][k1][FLAGUTOPRIMFAIL];  //backup the new inversion flag
			pflag[i1][j1][k1][FLAGUTOPRIMFAIL] = pflag_backup;   //restore the inversion flag

			if( 
				//(pflag_backup==UTOPRIMFAILUNEG || pflag_backup==UTOPRIMNOFAIL) && //not sure if can use this since it will always fix up the failure in the end
				(pflag_current==UTOPRIMFAILUNEG || pflag_current==UTOPRIMNOFAIL)   //if only u < 0 or no failure it is OK (since only check rho & gamma)
				) {

					//4. limit the flux c2a correction (difference btw. pr & pr_updated based on this difference
					frac_point_flux_used_array[index] = limit_prim_correction(MAX_AC_PRIM_FRAC_CHANGE, &geom, pr[i1][j1][k1], pr_updated);

			}
			else {
				frac_point_flux_used_array[index] = 1.; //if inversion from new value failed, let there be no c2a correction done to flux
			}
		}

		frac_point_flux_used = MAX( frac_point_flux_used_array[0], 	frac_point_flux_used_array[1] );

		if( 0 < frac_point_flux_used ) { 
			dualfprintf( fail_file, "Limited flux integration, dir = %d, i = %d, j = %d\n", dir, i, j );
		}

		PLOOP(pl) fluxvec_avg[i1][j1][k1][pl] = frac_point_flux_used * fluxvec_point[i1][j1][k1][pl] + (1.0 - frac_point_flux_used) * fluxvec_avg[i1][j1][k1][pl];
		
	} //end ZLOOP

	return( 0 );

}



int fluxcalc_1d(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt)
{
  int fluxcalc_standard(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt);
  //  int fluxcalc_fluxspliteno(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE F[][N2M][N3M][NPR], FTYPE *ndt);
  
  if((DOENOFLUX==NOENOFLUX)||(DOENOFLUX==ENOFLUXRECON)||(DOENOFLUX==ENOFINITEVOLUME)){  //atch correct: added "||(DOENOFLUX==ENOFINITEVOLUME)"
    MYFUN(fluxcalc_standard(stage,pr,dir,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,CUf,ndt),"flux.c:fluxcalc_1d()", "fluxcalc_standard()", 1);
  }
  else if(DOENOFLUX==ENOFLUXSPLIT){
    //MYFUN(fluxcalc_fluxspliteno(stage,pr,dir,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,CUf,ndt),"flux.c:fluxcalc_1d()", "fluxcalc_fluxspliteno()", 1);
  }

  return(0);
  
}







int fluxcalc_standard(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt)
{
	int check_plpr(int dir, int i, int j, int k, int idel, int jdel, int kdel, struct of_geom *geom, FTYPE pr[][N2M][N3M][NPR], FTYPE *p_l, FTYPE *p_r);
	void slope_lim(int numprims, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP]);
	int rescale(int which, int dir, FTYPE *pr, struct of_geom *geom,FTYPE*newvar);
	int global_vchar(FTYPE pointspeed[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE wspeed[][2][N1M][N2M][N3M]);
	void getp2interplr(int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r);
	FTYPE (*p2interp)[N2M][N3M][NPR2INTERP];
	int i, j, k, pl;
	FTYPE p_l[NPR], p_r[NPR];
	FTYPE p2interp_l[NPR2INTERP],p2interp_r[NPR2INTERP];
	FTYPE dtij;
	FTYPE ctop;
	struct of_geom geom;
	int p2SFUevolve(int dir, FTYPE *p, struct of_geom *geom, struct of_state *state, FTYPE *F, FTYPE *U);
	int reallim;
	int locallim;
	extern int flux_compute_general(int i, int j, int k, int dir, struct of_geom *geom, FTYPE CUf, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *ctop);
	extern int get_global_wavespeeds(int dir,FTYPE *pr,struct of_geom *ptrgeom, FTYPE *pleft);



	//////////////////////////
	//
	// rescale before interpolation
	//
	////////////////////////////
#if(RESCALEINTERP)
	// assume if DOEXTRAINTERP==1, then must get here
	p2interp=prc; // it's different
#else
	p2interp=pr; // it's itself
#endif


	//////////////////////////
	//
	// rescale before interpolation AND/OR get wavespeeds
	//
	////////////////////////////

#if((RESCALEINTERP)|| ( STOREWAVESPEEDS) )

	FULLLOOP{
		// get geometry for center pre-interpolated values
		get_geometry(i, j, k, CENT, &geom); 

#if(RESCALEINTERP)
		// assume no need for a guess to p2interp to get pr
		rescale(1,dir,pr[i][j][k],&geom,p2interp[i][j][k]);
#endif

#if(STOREWAVESPEEDS)
		MYFUN(get_global_wavespeeds(dir,pr[i][j][k],&geom,&state,pleft[i][j][k]),"flux.c:fluxcalc_standard()", "get_global_wavespeeds()", 0);
#endif // end if STOREWAVESPEEDS
	}// end FULLLOOP

#endif // end if rescaling or storing wavespeeds




#if(STOREWAVESPEEDS)
	// get char. velocity estimate as some average over local or global zones
	global_vchar(pleft, dir, is, ie, js, je, ks, ke, idel, jdel, kdel, wspeed);
	// now wspeed contains left/right fastest wave speed (with sign preserved)
#endif  // otherwise use very local estimate



	/////////////////////////////////////
	//
	// evaluate slopes (dq) or get pleft/pright of (possibly rescaled) primitive variables
	// c2e reconstruction: p2interp -> pleft & pright (indexed by grid cell #) -- atch comment
	//
	/////////////////////////////////////
	slope_lim(NPR2INTERP,dir,idel,jdel,kdel,pr,p2interp,dq,pleft,pright);




	//////////////////////////////////////
	//
	// flux loop : Extra "expand" zone for the purpose of averaging flux to get emf at corner.  Only used by field components, see flux_ct().
	// This loop is over interfaces where fluxes are evaluated -- atch
	//
	////////////////////////////////////////

#if((SIMULBCCALC==2)&&(TYPE2==1))
	FZLOOP(is,js,ks)
#else
	ZSLOOP( is, ie, js, je, ks, ke ) 
#endif
	{


		//////////////////////////////////
		//
		// this avoids problems (bad fluxes) on the pole
		//
		///////////////////////////////////
#if(ZEROPOLEFLUX==1)
		if((ZEROPOLEFLUX==1)&&(dir == 2) && ( ((startpos[2]+j == 0)&&(BCtype[X2DN]==POLARAXIS)) || ((startpos[2]+j == totalsize[2])&&(BCtype[X2UP]==POLARAXIS))  )) {
			// designed to just outright avoid the calculation.
			PLOOP(pl) F[i][j][k][pl] = 0. ;
			// GODMARK: not correct for pl==U2, where still have pressure
		}
		else {
#endif

			//////////////////////////////////////
			//
			// interpolate primitive using slope (dq) or directly from pleft and pright
			// For FV: p_left, p_right (indexed by grid cell #) -> p2interp_l, p2interp_r (indexed by interface #) -- atch comment
			//
			/////////////////////////////////////
			getp2interplr(dir,idel,jdel,kdel,i,j,k,p2interp,dq,pleft,pright,p2interp_l,p2interp_r);


			////////////////////////
			//
			// get the geometry for the flux face
			//
			///////////////////////
			get_geometry(i, j, k, face, &geom);



			/////////////////////////////////////
			//      
			// after interpolation, unrescale from p2interp to normal primitive 
			//
			///////////////////////////////////
#if(RESCALEINTERP)
			// setup plausible p_l/p_r in case used for some reason (e.g. inversion starting guess)
			// this is sufficient for utoprim_1D...a better guess does no better (see interpU code
			PLOOP(pl){
				p_l[pl]=pr[i][j][k][pl];
				p_r[pl]=pr[i][j][k][pl];
			}
			rescale(-1,dir,p_l,&geom,p2interp_l);
			rescale(-1,dir,p_r,&geom,p2interp_r);
#else
			PLOOP(pl){
				p_l[pl]=p2interp_l[pl];
				p_r[pl]=p2interp_r[pl];
			}
#endif


			//////////////////////////////
			//
			// Must preserve divb in 1D Riemann problem, so B^{dir} must be continuous
			//
			// GODMARK: Should really interpolate SUCH THAT this is automatically satisfied?
			//
			// GODMARK: This does not enforce E_\perp to be continuous for stationary flow!
			//
			///////////////////////////////

#if(BDIRCONT)
			if(dir==1) p_l[B1]=p_r[B1]=0.5*(p_l[B1]+p_r[B1]);
			else if(dir==2) p_l[B2]=p_r[B2]=0.5*(p_l[B2]+p_r[B2]);
			else if(dir==3) p_l[B3]=p_r[B3]=0.5*(p_l[B3]+p_r[B3]);
#endif

			///////////////////////
			//
			// correct interpolated quantities
			// no fixup accounting for these intermediate quantities
			//
			//////////////////////
			MYFUN(check_plpr(dir, i, j, k, idel, jdel, kdel, &geom, pr, p_l, p_r),"step_ch.c:fluxcalc()", "check_plpr()", 1);


			//////////////////////////////////
			//
			// actually compute the flux
			//
			/////////////////////////////////

			MYFUN(flux_compute_general(i, j, k, dir, &geom, CUf,  p_l, p_r, F[i][j][k], &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);


			/////////////////////////////
			//
			// evaluate restriction on timestep
			//
			///////////////////////////////

			dtij = cour * dx[dir] / ctop;
			if (dtij < *ndt) *ndt = dtij;


#if(ZEROPOLEFLUX==1) // close crazy ifdef
		}
#endif
	}

	// GODMARK: Should unrescale p2interp or something so that can use it in Athena2 part of fluxct rather than regenerating interpolated quantities
	// no longer needed since acts on different stored quantity, not pr
	//  if(RESCALEINTERP){
	// rescale before interpolation
	//    PREDQZLOOP{
	//  get_geometry(i, j, k, CENT, &geom);
	//  rescale(-1,dir,pr[i][j][k],&geom);
	// }
	// }


	return (0);
}






int check_plpr(int dir, int i, int j, int k, int idel, int jdel, int kdel, struct of_geom *geom, FTYPE pr[][N2M][N3M][NPR], FTYPE *p_l, FTYPE *p_r)
{
  int pl;

#if(EVOLVECHECKS)
#if(WHICHVEL==VEL4)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_4vel(dir,p_l,geom,-1.0);
  inflow_check_4vel(dir,p_r,geom,-1.0);
#endif
#elif(WHICHVEL==VEL3)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_3vel(dir,p_l,geom,-1.0);
  inflow_check_3vel(dir,p_r,geom,-1.0);
#endif
#if(JONCHECKS2)
  // must verify if this p makes sense (u^t sense)
  MYFUN(check_pr(p_l,pr[i-idel][j-jdel][k-kdel],geom,1,-1.0),"flux.c:check_plpr()", "check_pr()", 1);	
  
  MYFUN(check_pr(p_r,pr[i][j][k],geom,2,-1.0),"flux.c:check_plpr()", "check_pr()", 2);
#endif
#elif(WHICHVEL==VELREL4)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_rel4vel(dir,p_l,geom,-1.0);
  inflow_check_rel4vel(dir,p_r,geom,-1.0);
#endif
  // need to limit gamma since gamma may be large for interpolated value and would lead to bad fluxes
  MYFUN(limit_gamma(GAMMAMAX,p_l,geom,-1),"flux.c:check_plpr()", "limit_gamma()", 1);  //jon corr, see email re: SPINT warnings from 4/24/2006 10:54 p.m.
  MYFUN(limit_gamma(GAMMAMAX,p_r,geom,-1),"flux.c:check_plpr()", "limit_gamma()", 2);  //jon corr
#endif// end if WHICHVEL==VEL4REL      
#endif


#if(PRODUCTION==0)
  PLOOP(pl) if(!isfinite(p_l[pl])){
    dualfprintf(fail_file,"p_l is a nan pl=%d\n",pl);
  }
  PLOOP(pl) if(!isfinite(p_r[pl])){
    dualfprintf(fail_file,"p_r is a nan pl=%d\n",pl);
  }
#endif

  return(0);
}







//////////////////////////////////////
//
// interpolate primitive using slope or just copy pleft/pright into the correct place
//
// |=interface
// i=zone center of ith zone
//
// |              |       dq(i)        |
// |         pl(i)|pr(i)    i          |
// |              |pleft(i)   pright(i)|
//
//
//
//////////////////////////////////////

void getp2interplr(int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r)
{
  FTYPE Xcent[NDIM],Xleft[NDIM];
  FTYPE Vcent[NDIM],Vleft[NDIM];
  FTYPE rleft,rcent,thleft,thcent;
  int pl;
  int locallim;
  int choose_limiter(int i, int j, int k, int pl);



	if((HORIZONSUPERFAST)&&((lim<PARA)&&(LIMADJUST==0))&&(dir==1)){
		// since this uses dq's to shift, must always have dq's everywhere.  Hence LIMADJUST==0 and lim<PARA must be true
		coord(im1, j, k, CENT, Xleft);
		coord(i, j, k, CENT, Xcent);
		bl_coord(Xleft, Vleft);
		rleft=Vleft[1];
		thleft=Vleft[2];
		bl_coord(Xcent, Vcent);
		rcent=Vcent[1];
		thcent=Vcent[2];

		PLOOPINTERP(pl){
			locallim=choose_limiter(i,j,k,pl);
			// get interpolated quantity
			if((locallim<PARA)&&(LIMADJUST==0)){
				if(rleft>Rhor) p2interp_l[pl] = p2interp[i - idel][j - jdel][k - kdel][pl] + 0.5 * dq[i - idel][j - jdel][k - kdel][pl];
				else p2interp_l[pl] = p2interp[i][j][k][pl] - 0.5 * dq[i][j][k][pl];
				if(rcent>Rhor) p2interp_r[pl] = p2interp[i][j][k][pl] - 0.5 * dq[i][j][k][pl];
				else p2interp_r[pl] = p2interp[i+idel][j+jdel][k+kdel][pl] - 1.5 * dq[i+idel][j+jdel][k+kdel][pl];
			}
			else{
				p2interp_l[pl] = pright[i-idel][j-jdel][k-kdel][pl];
				p2interp_r[pl] = pleft[i][j][k][pl];
			}
		} // end PLOOP
	}// if horizonsuperfast
	else{
		/////////////////////////////
		//
		// standard routine
		//
		/////////////////////////////
		PLOOPINTERP(pl){
			locallim=choose_limiter(i,j,k,pl);
			// get interpolated quantity
			if((locallim<PARA)&&(LIMADJUST==0)){
				p2interp_l[pl] = p2interp[i - idel][j - jdel][k - kdel][pl] + 0.5 * dq[i - idel][j - jdel][k - kdel][pl];
				p2interp_r[pl] = p2interp[i][j][k][pl] - 0.5 * dq[i][j][k][pl];
			}
			else{
				//p_l & p_r for the current interface -- atch comment
				p2interp_l[pl] = pright[i-idel][j-jdel][k-kdel][pl];
				p2interp_r[pl] = pleft[i][j][k][pl];
			}
		}
	}
}





// choose limiter
int choose_limiter(int i, int j, int k, int pl)
{
#if(LIMADJUST==LIMITERFIXED)
	return(lim);
#else

#if(HYDROLIMADJUSTONLY)
	if(pl<B1) return(pflag[i][j][k][FLAGREALLIM]);
	else return(lim);
#else
	return(pflag[i][j][k][FLAGREALLIM]);
#endif

#endif

}






// slope_lim() is provided p2interp and returns pleft/pright
//
// |=interface
// i=zone center of ith zone
//
// |              |      p2interp(i)   |
// |         pl(i)|pr(i)    i          |
// |         Fl(i)|Fr(i)    i          |
// |         Ul(i)|Ur(i)    i          |
// |              |pleft(i)   pright(i)|
// |              |F(i)                |
//
void slope_lim(int numprims, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP])
{
	extern void slope_lim_linetype_c2e(int numprims, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP]);
	extern void slope_lim_pointtype(int numprims, int pl, int dir, int idel, int jdel, int kdel, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP]);
	int pl;


	if( lim>=WENO3 ){ // this overrides lim, but lim must still be set properly
		slope_lim_linetype_c2e(numprims, ENOINTERPTYPE, dir, idel, jdel, kdel, primreal, p2interp, pleft, pright);
	}
	else{
		PLOOPINTERP(pl){
			slope_lim_pointtype(numprims, pl, dir, idel, jdel, kdel, p2interp, dq, pleft, pright);
		}
	}

}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// GENERAL PRIMITIVE INTERPOLATION CHANGE OF VARIABLES
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


// notice that input order is always:
// 2nd arg: true primitive
// last arg: scaled primitive

// this also allows any fancy remapping, such as characteristic interpolation -- rest of code is setup to allow any remapping as long as you have an inversion
int rescale(int which, int dir, FTYPE *pr, struct of_geom *ptrgeom,FTYPE *p2interp)
{
  FTYPE scale[NPR2INTERP],r,th,X[NDIM],V[NDIM];
  int ii,jj,kk,pl;
  struct of_state q;
  extern int quasivsq_compute(FTYPE *pr, struct of_geom *geom, FTYPE *quasivsq);
  FTYPE quasivsq;
  extern int limit_quasivsq(FTYPE quasivsqnew, struct of_geom *geom, FTYPE *pr);
  FTYPE vcon[NDIM],ucon[NDIM];
  int j;
  FTYPE newpr[NPR];
  FTYPE uconrel[NDIM],ucovrel[NDIM];
  FTYPE vconrel[NDIM];
  extern int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma);
  FTYPE normuconrel,normuconrel_fromui;
  FTYPE gamma;

#if(VARTOINTERPFIELD==PULSARFIELD)
  extern void getconsts(FTYPE *uconmetp, FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconconst);
  extern void undoconsts(FTYPE *uconconst, FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetp);
  FTYPE Bconin[NDIM],Bconout[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
#endif


  ii=ptrgeom->i;
  jj=ptrgeom->j;
  kk=ptrgeom->k;

  coord(ii,jj,kk,ptrgeom->p,X);
  bl_coord(X,V);
  r=V[1]; th=V[2];

#if(VARTOINTERP==PRIMTOINTERP)

  //  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if RESCALEINTERP=%d\n",VARTOINTERP,RESCALEINTERP);
  //  myexit(1);

  // allow multiple types of rescales, and so just identity if not doing this (i.e. PRIMTOINTERP)
  if(which==1){ // rescale before interpolation
    PLOOP(pl) p2interp[pl]=pr[pl];
  }
  else if(which==-1){ // unrescale after interpolation
    PLOOP(pl) pr[pl]=p2interp[pl];
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }




#elif(VARTOINTERP==PRIMTOINTERP_JONRESCALED1)

  if(dir==1){
    // optimized for pole
    if(rescaletype==0){
      scale[RHO]=pow(r,-1.5);
    }
    else if(rescaletype==1){
      scale[RHO]=pow(r,-2.7);
    }
    scale[UU]=scale[RHO]/r;
    scale[U1]=scale[RHO];
    scale[U2]=1.0;
    scale[U3]=1.0/(r*r);
    if(rescaletype==0){
      scale[B1]=scale[U3];
    }
    else if(rescaletype==1){
      scale[B1]=pow(r,-2.4);
    }
    //    if(statpos[2]+jj < 0 || startpos[2]+jj >= totalsize[2]) scale[B1] *= -1. ;
    scale[B2]=scale[B1];
    scale[B3]=scale[B1];

    if(DOENTROPY) scale[ENTROPY]=1.0;
  }
  else if(dir==2){
    scale[RHO]=1.0;
    scale[UU]=1.0;
    scale[U1]=1.0;
    scale[U2]=1.0;
    scale[U3]=1.0;
    scale[B1]=1.0;
    scale[B2]=1.0;
    scale[B3]=1.0;
    if(DOENTROPY) scale[ENTROPY]=1.0;
  }
  else{
    dualfprintf(fail_file,"rescale(): no such direction! dir=%d\n",dir);
    myexit(100);
  }


  if(which==1){ // rescale before interpolation
    PLOOP(pl) p2interp[pl]=pr[pl]/scale[pl];
  }
  else if(which==-1){ // unrescale after interpolation
    PLOOP(pl) pr[pl]=p2interp[pl]*scale[pl];
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==CONSTOINTERP)
  // this doesn't work at all, even if no bug.
  // doesn't work even if setup better guess, as in interpU code.


  if(which==1){ // rescale before interpolation
    MYFUN(get_state(pr, ptrgeom, &q),"interp.c:rescale()", "get_state()", 1);
    MYFUN(primtoU(UDIAG,pr, &q, ptrgeom, p2interp),"interp.c:rescale()", "primtoU()", 1);
  }
  else if(which==-1){ // unrescale after interpolation
    MYFUN(Utoprimgen(OTHERUTOPRIM,UDIAG,p2interp, ptrgeom, pr),"interp.c:rescale()", "Utoprimgen", 1);
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#elif(VARTOINTERP==PRIMTOINTERPLGDEN)


  if(which==1){ // rescale before interpolation
    PLOOP(pl) p2interp[pl]=pr[pl];
    p2interp[RHO]=log(pr[RHO]);
    p2interp[UU]=log(pr[UU]);
  }
  else if(which==-1){ // unrescale after interpolation
    PLOOP(pl) pr[pl]=p2interp[pl];
    pr[RHO]=exp(p2interp[RHO]);
    pr[UU]=exp(p2interp[UU]);
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }

#elif(VARTOINTERP==PRIMTOINTERP_LGDEN_RHOU)
  // unstable!
  // probably because using LGDEN with RHOU.
  // Between 1E-6 and 1, density would be 1E-3 in log.
  // But v=0 and 1, so rho*v~0 and 1, but then resulting v=10^3 !

  if(which==1){ // rescale before interpolation
    PLOOP(pl) p2interp[pl]=pr[pl];
    p2interp[RHO]=log(pr[RHO]);
    p2interp[UU]=log(pr[UU]);

    p2interp[U1]=pr[RHO]*pr[U1];
    p2interp[U2]=pr[RHO]*pr[U2];
    p2interp[U3]=pr[RHO]*pr[U3];
  }
  else if(which==-1){ // unrescale after interpolation
    PLOOP(pl) pr[pl]=p2interp[pl];
    pr[RHO]=exp(p2interp[RHO]);
    pr[UU]=exp(p2interp[UU]);

    pr[U1]=p2interp[U1]/pr[RHO];
    pr[U2]=p2interp[U2]/pr[RHO];
    pr[U3]=p2interp[U3]/pr[RHO];

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_RHOU)
  // kinda works, except for test=7 (peak problem)

  if(which==1){ // rescale before interpolation
    PLOOP(pl) p2interp[pl]=pr[pl];

    p2interp[U1]=pr[RHO]*pr[U1];
    p2interp[U2]=pr[RHO]*pr[U2];
    p2interp[U3]=pr[RHO]*pr[U3];
  }
  else if(which==-1){ // unrescale after interpolation
    PLOOP(pl) pr[pl]=p2interp[pl];

    pr[U1]=p2interp[U1]/pr[RHO];
    pr[U2]=p2interp[U2]/pr[RHO];
    pr[U3]=p2interp[U3]/pr[RHO];

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }




#elif(VARTOINTERP==PRIMTOINTERP_VSQ)



#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif


  if(which==1){ // before interpolation, get quantities to interpolate
    quasivsq_compute(pr, ptrgeom, &quasivsq);
    PLOOP(pl) p2interp[pl]=pr[pl];
    
    p2interp[U1]=pr[RHO]*pr[U1];
    p2interp[U2]=pr[RHO]*pr[U2];
    p2interp[U3]=pr[RHO]*pr[U3];
    
    // max helps limit oscillatory behavior for non-limiter schemes
    p2interp[VSQ]=max(quasivsq,0.0);
    //p2interp[VSQ]=log(quasivsq); // assumes positive
  }
  else  if(which==-1){ // after interpolation
    PLOOP(pl) pr[pl]=p2interp[pl];

    pr[U1]=p2interp[U1]/pr[RHO];
    pr[U2]=p2interp[U2]/pr[RHO];
    pr[U3]=p2interp[U3]/pr[RHO];

    // now rescale velocities to agree with quasivsq
    // max helps limit oscillatory behavior for non-limiter schemes
    limit_quasivsq(max(p2interp[VSQ],0.0),ptrgeom,pr);
    //    limit_quasivsq(exp(p2interp[VSQ]),ptrgeom,pr);
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_3VEL_GAMMA)


#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif


  if(which==1){ // before interpolation, get quantities to interpolate
    pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
    
    // 3-velocity
    SLOOPA(j) vcon[j]=ucon[j]/ucon[TT];

    PLOOP(pl) p2interp[pl]=pr[pl];

    for(pl=U1;pl<=U3;pl++) p2interp[pl]= vcon[pl-U1+1];

    p2interp[VSQ]=ucon[TT];

  }
  else  if(which==-1){ // after interpolation

    PLOOP(pl) pr[pl]=p2interp[pl];

    for(pl=U1;pl<=U3;pl++)  vcon[pl-U1+1] = p2interp[pl]; 

    SLOOPA(j) ucon[j] = vcon[j]*p2interp[VSQ];

    ucon2pr(WHICHVEL,ucon,ptrgeom,pr);

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_RHOV_GAMMA)

#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif


  if(which==1){ // before interpolation, get quantities to interpolate
    pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
    
    // 3-velocity
    SLOOPA(j) vcon[j]=ucon[j]/ucon[TT];

    PLOOP(pl) p2interp[pl]=pr[pl];

    //comment this out if you do not want to interpolate rho * v
    for(pl=U1;pl<=U3;pl++) p2interp[pl]= pr[RHO] * vcon[pl-U1+1];

    p2interp[VSQ]=ucon[TT];

  }
  else  if(which==-1){ // after interpolation

    PLOOP(pl) pr[pl]=p2interp[pl];

    //comment this out if you do notwant to interp. rhov
    for(pl=U1;pl<=U3;pl++)  vcon[pl-U1+1] = p2interp[pl]/pr[RHO]; 

    SLOOPA(j) ucon[j] = vcon[j]*p2interp[VSQ];

    ucon2pr(WHICHVEL,ucon,ptrgeom,pr);

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#elif(VARTOINTERP==PRIMTOINTERP_VELREL4SQ)


#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif

#if(WHICHVEL==VEL3)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if WHICHVEL=%d\n",VARTOINTERP,WHICHVEL);
  myexit(1);
#endif



  if(which==1){ // before interpolation, get quantities to interpolate

    // get relative 4-velocity
    if(WHICHVEL!=VELREL4){
      pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
      ucon2pr(VELREL4,ucon,ptrgeom,newpr);
      uconrel[TT]=0.0;
      SLOOPA(j) uconrel[j]=newpr[UU+j];
    }
    else{
      uconrel[TT]=0.0;
      SLOOPA(j) uconrel[j]=pr[UU+j];
    }
    
    // get |\tilde{u}|
    lower_vec(uconrel,ptrgeom,ucovrel);
    //normuconrel=sqrt(dot(uconrel,ucovrel));
    normuconrel=dot(uconrel,ucovrel);


    PLOOP(pl) p2interp[pl]=pr[pl];
    p2interp[VSQ]=normuconrel;

  }
  else  if(which==-1){ // after interpolation

    // assign over everything, adjust velocity below
    PLOOP(pl) pr[pl]=p2interp[pl];

    // renormalize by relative 4-velocity
    if(WHICHVEL!=VELREL4){
      // get relative 4-velocity from interpolated velocity
      pr2ucon(WHICHVEL,p2interp, ptrgeom ,ucon);
      ucon2pr(VELREL4,ucon,ptrgeom,newpr);
      uconrel[TT]=0.0;
      SLOOPA(j) uconrel[j]=newpr[UU+j];
    }
    else{
      uconrel[TT]=0.0;
      SLOOPA(j) uconrel[j]=p2interp[UU+j];
    }

    // get |\tilde{u}| from interpolated \tilde{u}^i
    lower_vec(uconrel,ptrgeom,ucovrel);
    //    normuconrel_fromui=sqrt(dot(uconrel,ucovrel));
    normuconrel_fromui=dot(uconrel,ucovrel);

    if(WHICHVEL==VEL3){
      //      pr2ucon(WHICHVEL,p2interp,ptrgeom,ucon);
      //      ucon2pr(VELREL4,ucon,ptrgeom,newpr); // now have relative 4-velocity
      //      uconrel[TT]=0.0;
      //      SLOOPA(j) uconrel[j]=newpr[UU+j];
      // not done, but should be possible
      myexit(1);
    }
    else{ // WHICHVEL==VEL4 or WHICHVEL==VELREL4 : acts same on both
      // directly renormalize primitives
      //      SLOOPA(j) pr[UU+j] = p2interp[UU+j]*p2interp[VSQ]/absuconrel_fromui;
      SLOOPA(j) pr[UU+j] = p2interp[UU+j]*fabs(p2interp[VSQ]/normuconrel_fromui);
      // done!
    }





  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_3VELREL_GAMMAREL)


#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif


  if(which==1){ // before interpolation, get quantities to interpolate

    // get relative 4-velocity
    if(WHICHVEL!=VELREL4){
      pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
      ucon2pr(VELREL4,ucon,ptrgeom,newpr);
      uconrel[TT]=0.0;
      SLOOPA(j) uconrel[j]=newpr[UU+j];
    }
    else{
      uconrel[TT]=0.0;
      SLOOPA(j) uconrel[j]=pr[UU+j];
    }

    // get Lorentz factor w.r.t. relative 4-velocity
    gamma_calc_fromuconrel(uconrel,ptrgeom,&gamma);

    PLOOP(pl) p2interp[pl]=pr[pl];

    // interpolate relative 3-velocity
    for(pl=U1;pl<=U3;pl++) p2interp[pl]= uconrel[pl-U1+1]/gamma;

    // interpolate \gamma separately
    p2interp[VSQ]=gamma;

  }
  else  if(which==-1){ // after interpolation

    // assign over everything, adjust velocity below
    PLOOP(pl) pr[pl]=p2interp[pl];

    // get relative 4-velocity from \gamma and relative 3-velocity
    uconrel[TT]=0;
    SLOOPA(j) uconrel[j]=p2interp[UU+j]*p2interp[VSQ];

    // get WHICHVEL velocity
    if(WHICHVEL!=VELREL4){
      pr2ucon(VELREL4,p2interp, ptrgeom ,ucon);
      ucon2pr(WHICHVEL,ucon,ptrgeom,newpr);
      SLOOPA(j) pr[UU+j]=newpr[UU+j];
    }
    else{
      SLOOPA(j) pr[UU+j]=uconrel[j];
    }

    

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_RAMESH1)


  if(which==1){ // rescale before interpolation
    PLOOP(pl) p2interp[pl]=pr[pl];

    // for the fields, do ramesh-like interpolation
    pl=B1; p2interp[pl] = pr[pl] * (sqrt(ptrgeom->gcov[1][1])*pow(r,2.0-nu) );
    pl=B2; p2interp[pl] = pr[pl] * (sqrt(ptrgeom->gcov[2][2])*pow(r,2.0-nu) );
    pl=B3; p2interp[pl] = pr[pl] * (sqrt(ptrgeom->gcov[3][3])*pow(r,2.0-nu) );
  }
  else if(which==-1){ // unrescale after interpolation
    PLOOP(pl) pr[pl]=p2interp[pl];

    // for the fields, do ramesh-like interpolation
    pl=B1; pr[pl] = p2interp[pl] / (sqrt(ptrgeom->gcov[1][1])*pow(r,2.0-nu) );
    pl=B2; pr[pl] = p2interp[pl] / (sqrt(ptrgeom->gcov[2][2])*pow(r,2.0-nu) );
    pl=B3; pr[pl] = p2interp[pl] / (sqrt(ptrgeom->gcov[3][3])*pow(r,2.0-nu) );

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#endif // end over choices for VARTOINTERP


#if(VARTOINTERPFIELD==PULSARFIELD&&0)


  if(dir==1){
    if(which==1){ // rescale before interpolation

      Bconin[0]=0.0;
      Bconin[1]=pr[B1];
      Bconin[2]=pr[B2];
      Bconin[3]=pr[B3];

      dxdxprim(X, V, dxdxp);

      getconsts(Bconin, V, ptrgeom, dxdxp,Bconout);

      p2interp[B1]=Bconout[1];
      p2interp[B2]=Bconout[2];
      p2interp[B3]=Bconout[3];
    }
    else if(which==-1){ // unrescale after interpolation

      Bconin[0]=0.0;
      Bconin[1]=p2interp[B1];
      Bconin[2]=p2interp[B2];
      Bconin[3]=p2interp[B3];

      dxdxprim(X, V, dxdxp);

      undoconsts(Bconin, V, ptrgeom, dxdxp, Bconout);

      pr[B1]=Bconout[1];
      pr[B2]=Bconout[2];
      pr[B3]=Bconout[3];
      
    }
    else{
      dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
      myexit(100);
    }
  }
  // problem when used with dir=2 since near axis end up dividing by 0


#endif // end over choices for VARTOINTERPFIELD



  return(0);
}


