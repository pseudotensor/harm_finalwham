#include "decs.h"




// pi: initial values at t=t0 to compute Ui
// pb: values used to compute flux/source
// pf: solution using flux(pb) from pi's Ui -> Uf

// pi, pb, and pf can all be the same since
// 1) pb used first on a stencil, not modified, to compute fluxes
// 2) pf=pi is assigned by value at each zone
// 3) pf is modified using Utoprim at each zone using pb for sources (to correspond to fluxes which used pb)
//
// So in the end only pf is modified at each zone, so the loop changing p at previous (i,j) location doesn't affect the any new location in (i,j)
int advance(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
						FTYPE ui[][N2M][N3M][NPR],FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
						FTYPE *CUf,FTYPE *Cunew,int timeorder, int numtimeorders, FTYPE *ndt)
{
	//int advance(int stage,
	//	    FTYPE pi[][N2M][N3M][NPR],
	//	    FTYPE ulast[][N2M][N3M][NPR], 
	//	    FTYPE pb[][N2M][N3M][NPR],
	//	    FTYPE *CUf, FTYPE pf[][N2M][N3M][NPR],
	//	    FTYPE *Cunew, FTYPE unew[][N2M][N3M][NPR],
	//	    int timeorder, int numtimeorders,
	//	    FTYPE *ndt)
	//{
	int advance_standard(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
		FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
		FTYPE *CUf,FTYPE *Cunew,int stagenow, int numstages, FTYPE *ndt);
	int advance_finitevolume(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
		FTYPE ui[][N2M][N3M][NPR],FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
		FTYPE *CUf,FTYPE *Cunew,int stagenow, int numstages, FTYPE *ndt);
	//  int advance_standard(int stage, FTYPE pi[][N2M][N3M][NPR], FTYPE ulast[][N2M][N3M][NPR], FTYPE pb[][N2M][N3M][NPR], FTYPE *CUf, FTYPE pf[][N2M][N3M][NPR], FTYPE *Cunew, FTYPE unew[][N2M][N3M][NPR], int timeorder, int numtimeorders, FTYPE *ndt);


	if(DOENOFLUX==ENOFINITEVOLUME){
		MYFUN(advance_finitevolume(stage,pi,pb,pf,ui,uf,ucum,CUf,Cunew,timeorder,numtimeorders,ndt),"advance.c:advance()", "advance_finitevolume()", 1);
	}
	else if((DOENOFLUX==NOENOFLUX)||(DOENOFLUX==ENOFLUXRECON)||(DOENOFLUX==ENOFLUXSPLIT)){
		MYFUN(advance_standard(stage,pi,pb,pf,uf,ucum,CUf,Cunew,timeorder,numtimeorders,ndt),"advance.c:advance()", "advance_standard()", 1);
	}
	else{
		dualfprintf(fail_file,"No such DOENOFLUX=%d\n",DOENOFLUX);
		myexit(1);
	}


	return(0);


}




int advance_standard(int stage,
										 FTYPE pi[][N2M][N3M][NPR],
										 FTYPE pb[][N2M][N3M][NPR],
										 FTYPE pf[][N2M][N3M][NPR],
										 FTYPE uf[][N2M][N3M][NPR],
										 FTYPE ucum[][N2M][N3M][NPR], 
										 FTYPE *CUf, 
										 FTYPE *Cunew, 
										 int timeorder, int numtimeorders,
										 FTYPE *ndt)
{
	int i, j, k, pl, sc;
	FTYPE ndt1, ndt2, ndt3;
	FTYPE Uf[NPR], Ui[NPR], Ub[NPR];
	FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR],dUriemann2[NPR],dUriemann3[NPR],dUcomp[NUMSOURCES][NPR];
	struct of_geom geom;
	struct of_state q;
	FTYPE dUtot;
	FTYPE idx1,idx2;
	SFTYPE dt4diag;
	int finalstep;
	void flux2dUavg(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],FTYPE *dUavg1,FTYPE *dUavg2,FTYPE *dUavg3);
	void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *uf, FTYPE *Uf, FTYPE *ucum);
	void flux2U(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE *dU, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui, FTYPE *uf, FTYPE *Uf, FTYPE *ucum);
	extern int fluxcalc(int stage, FTYPE pr[][N2M][N3M][NPR],
		FTYPE F1[][N2M][N3M][NPR], 
		FTYPE F2[][N2M][N3M][NPR], 
		FTYPE F3[][N2M][N3M][NPR], 
		FTYPE CUf,
		FTYPE *ndt1,
		FTYPE *ndt2,
		FTYPE *ndt3
		);
	extern void diag_source_comp(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt);
	extern void diag_source_all(struct of_geom *ptrgeom, FTYPE *dU,SFTYPE Dt);
	extern int diag_flux(FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],SFTYPE Dt);
	extern int diss_compute(int evolvetype, int inputtype, FTYPE *U, struct of_geom *ptrgeom, FTYPE *pr);
	int asym_compute_1(FTYPE (*prim)[N2M][N3M][NPR]);
	int asym_compute_2(FTYPE (*prim)[N2M][N3M][NPR]);


#if(ASYMDIAGCHECK)
	dualfprintf(fail_file,"BEGINNING steppart=%d nstep=%ld\n",steppart,nstep);


	dualfprintf(fail_file,"pi in advance\n");
	asym_compute_1(pi);

	dualfprintf(fail_file,"pb in advance\n");
	asym_compute_1(pb);
#endif


	// tells diagnostics functions if should be accounting or not
	if(timeorder==numtimeorders-1){
		dt4diag=dt;
		finalstep=1;
	}
	else{
		dt4diag=-1.0;
		finalstep=0;
	}

#if(PRODUCTION==0)
	trifprintf( "#0f");
#endif

	// pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
	MYFUN(fluxcalc(stage,pb,F1,F2,F3,CUf[2],&ndt1,&ndt2,&ndt3),"advance.c:advance_standard()", "fluxcalcall", 1);

	// compute flux diagnostics (accurately using all substeps)
#if(ACCURATEDIAG)
	diag_flux(pb,F1, F2, F3, Cunew[1]*dt); // Cunew[1]*dt is true dt for this flux as added in dUtoU()
#endif


#if(PRODUCTION==0)
	trifprintf( "1f");
#endif
	// from here on, pi/pb/pf are only used a zone at a time rather than on a stencil




	/** now update pi to pf **/
	CZLOOP {



		// initialize uf and ucum if very first time here
		if(timeorder==0) PLOOP(pl) uf[i][j][k][pl]=ucum[i][j][k][pl]=0.0;


		// set geometry for centered zone to be updated
		get_geometry(i, j, k, CENT, &geom);

		// find Ui(pi)
		MYFUN(get_state(pi[i][j][k], &geom, &q),"step_ch.c:advance()", "get_state()", 1);
		MYFUN(primtoU(UEVOLVE,pi[i][j][k], &q, &geom, Ui),"step_ch.c:advance()", "primtoU()", 1);

		// find dU(pb)
		// source() doesn't actually use CUf[2]=dt right now
		//    MYFUN(source(pb[i][j][k], &geom, i, j, k, dUcomp, dU,CUf[2]),"step_ch.c:advance()", "source", 1);
		MYFUN(source(pb[i][j][k], &geom, dUcomp, dUgeom),"step_ch.c:advance()", "source", 1);
		// assumes final dUcomp is nonzero and representative of source term over this timestep

#if(ACCURATEDIAG)
		diag_source_comp(&geom,dUcomp,Cunew[1]*dt);
		diag_source_all(&geom,dUgeom,Cunew[1]*dt);
#else
		diag_source_comp(&geom,dUcomp,dt4diag);
		diag_source_all(&geom,dUgeom,dt4diag);
#endif

		// dUriemann is actually average quantity, but we treat is as a point quantity at the zone center
		flux2dUavg(i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
		PLOOP(pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is one type of avg->point mistake

		// find uf==Uf and additional terms to ucum
		dUtoU(dUgeom, dUriemann, CUf, Cunew, Ui, uf[i][j][k], Uf, ucum[i][j][k]);
		//    flux2U(i,j,k,F1,F2,F3,dU,CUf,Cunew,Ui, uf[i][j][k],Uf,ucum[i][j][k]);


		/////////////
		//
		// Utoprim as initial conditions : can't assume want these to be same in the end, so assign
		//
		// Since final step often sets pointers of pi=pf, in order to use arbitrary guess, must set here once done with pi,pb,pf.
		//
		////////////
		PLOOP(pl) pf[i][j][k][pl] = pb[i][j][k][pl];


		// invert U->p
		if(timeorder==numtimeorders-1){ // last call, so ucum is cooked and ready to eat!
			MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,ucum[i][j][k], &geom, pf[i][j][k]),"step_ch.c:advance()", "Utoprimgen", 1);
#if(DODISS||DODISSVSR)
			// then see what entropy inversion would give
			diss_compute(EVOLVEUTOPRIM,UEVOLVE,ucum[i][j][k],&geom,pf[i][j][k]);
#endif

		}
		else{ // otherwise still iterating on primitives
			MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,Uf, &geom, pf[i][j][k]),"step_ch.c:advance()", "Utoprimgen", 1);
		}



		// immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
		// SUPERGODMARK: Below should differentiate beteween whether want negative densities fixed or not, but right now fixup1zone() does all
		if((STEPOVERNEGU==0)||(STEPOVERNEGRHO==0)||(STEPOVERNEGRHOU==0)||(timeorder==numtimeorders-1)){
			MYFUN(fixup1zone(pf[i][j][k],&geom,finalstep),"fixup.c:fixup()", "fixup1zone()", 1);
		}
#endif
	}// end CZLOOP


#if(ASYMDIAGCHECK)
	dualfprintf(fail_file,"ucum in advance\n");
	asym_compute_2(ucum);

	dualfprintf(fail_file,"ENDING steppart=%d nstep=%ld\n",steppart,nstep);
#endif


#if(FIXUPZONES==FIXUPALLZONES)
	fixup(stage,pf,finalstep);
#endif  




	*ndt = defcon * 1. / (1. / ndt1 + 1. / ndt2 + 1. / ndt3);

#if(PRODUCTION==0)
	trifprintf( "2f");
#endif

	return (0);
}






int advance_finitevolume(int stage,
												 FTYPE pi[][N2M][N3M][NPR],
												 FTYPE pb[][N2M][N3M][NPR],
												 FTYPE pf[][N2M][N3M][NPR],
												 FTYPE ui[][N2M][N3M][NPR],
												 FTYPE uf[][N2M][N3M][NPR],
												 FTYPE ucum[][N2M][N3M][NPR], 
												 FTYPE *CUf, 
												 FTYPE *Cunew, 
												 int timeorder, int numtimeorders,
												 FTYPE *ndt)
{
	int i, j, k, pl, sc;
	FTYPE ndt1, ndt2, ndt3;
	FTYPE Uf[NPR], Ui[NPR], Ub[NPR];
	FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR],dUriemann2[NPR],dUriemann3[NPR],dUcomp[NUMSOURCES][NPR];
	struct of_geom geom;
	struct of_state q;
	FTYPE dUtot;
	FTYPE idx1,idx2;
	SFTYPE dt4diag;
	int finalstep;
	void flux2dUavg(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],FTYPE *dUavg1,FTYPE *dUavg2,FTYPE *dUavg3);
	void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *uf, FTYPE *Uf, FTYPE *ucum);
	void flux2U(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE *dU, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui, FTYPE *uf, FTYPE *Uf, FTYPE *ucum);
	extern int fluxcalc(int stage, FTYPE pr[][N2M][N3M][NPR],
		FTYPE F1[][N2M][N3M][NPR], 
		FTYPE F2[][N2M][N3M][NPR], 
		FTYPE F3[][N2M][N3M][NPR], 
		FTYPE CUf,
		FTYPE *ndt1,
		FTYPE *ndt2,
		FTYPE *ndt3
		);
	extern int avg2cen_interp(int whichquantity, int whichavg2cen, FTYPE (*prims_from_avg_cons)[N2M][N3M][NPR], FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR]);
	extern void diag_source_comp(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt);
	extern void diag_source_all(struct of_geom *ptrgeom, FTYPE *dU,SFTYPE Dt);
	extern int diag_flux(FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],SFTYPE Dt);
	extern int diss_compute(int evolvetype, int inputtype, FTYPE *U, struct of_geom *ptrgeom, FTYPE *pr);

	FTYPE (*utoinvert)[N2M][N3M][NPR];
	static int firsttime=1;
	void debugeno_compute(FTYPE p[][N2M][N3M][NPR],FTYPE U[][N2M][N3M][NPR],FTYPE debugvars[][N2M][N3M][NUMENODEBUGS]);
	FTYPE (*myupoint)[N2M][N3M][NPR];
	int check_point_vs_average(int timeorder, int numtimeorders, int *lpflag, FTYPE *pb, FTYPE *pf, FTYPE *upoint, FTYPE *uavg, struct of_geom *ptrgeom);




	// tells diagnostics functions if should be accounting or not
	if(timeorder==numtimeorders-1){
		dt4diag=dt;
		finalstep=1;
	}
	else{
	  dt4diag=-1.0; // tells diag_source to not consider
		finalstep=0;
	}

#if(PRODUCTION==0)
	trifprintf( "#0f");
#endif

	// pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
	MYFUN(fluxcalc(stage,pb,F1,F2,F3,CUf[2],&ndt1,&ndt2,&ndt3),"advance.c:advance_standard()", "fluxcalcall", 1);

	// compute flux diagnostics (accurately using all substeps)
#if(ACCURATEDIAG)
	diag_flux(pb,F1, F2, F3, Cunew[1]*dt); // Cunew[1]*dt is true dt for this flux as added in dUtoU()
#endif


#if(PRODUCTION==0)
	trifprintf( "1f");
#endif
	// from here on, pi/pb/pf are only used a zone at a time rather than on a stencil



	////////////////
	//
	// Copy last full timestep's ucum into ui
	//
	///////////////
	if(firsttime){
		// assume ui is set correctly from initial conditions
		firsttime=0;
	}
	else{
		if(timeorder==0){
			FULLLOOP PLOOP(pl){
				ui[i][j][k][pl]=ucum[i][j][k][pl];
			}
		}
		// otherwise assume ui didn't change.  Present RK schemes assume this.  Otherwise have to keep track of pf/Uf pairs in RK stepping
	}


	/////////////////////////
	//
	// SOURCE TERM
	//
	////////////////////////
	// GODMARK: other/more special cases?
	FULLLOOP{
		// find dU(pb)
		// only find source term if non-Minkowski and non-Cartesian
		// set geometry for centered zone to be updated
		get_geometry(i, j, k, CENT, &geom);

		// get source term
		MYFUN(source(pb[i][j][k], &geom, dUcomp, dUgeomarray[i][j][k]),"step_ch.c:advance()", "source", 1);

#if(ACCURATEDIAG)
		diag_source_comp(&geom,dUcomp,Cunew[1]*dt);
#else
		// assumes final dUcomp is nonzero and representative of source term over this timestep
		diag_source_comp(&geom,dUcomp,dt4diag);
#endif
	}
	// volume integrate dUgeom
	// c2a_1 c2a_2 c2a_3
	if((avgscheme>3)&&(do_source_integration)){
		avg2cen_interp( ENOSOURCETERM, ENOCENT2AVGTYPE, pb, dUgeomarray, dUgeomarray);
	}





	// update Ui to Uf (ultimately to ucum)
	// expanded loop using expanded range of fluxes so update Uf and ucum in layer of ghost zones so don't have to bound flux or Uf/Ui/ucum
	ZSLOOP( Uconsloop[FIS], Uconsloop[FIE], Uconsloop[FJS], Uconsloop[FJE], Uconsloop[FKS], Uconsloop[FKE] ){

	  // get geometry at center where source is located
		get_geometry(i, j, k, CENT, &geom);

		// get source term (volume integrated)
		PLOOP(pl) dUgeom[pl]=dUgeomarray[i][j][k][pl];
#if(ACCURATEDIAG)
		// do diagnostics for volume integrated source term
		diag_source_all(&geom,dUgeom,Cunew[1]*dt);
#else
		diag_source_all(&geom,dUgeom,dt4diag);
#endif



		/////////////
		//
		// Utoprim as initial conditions : can't assume want these to be same in the end, so assign
		//
		////////////


		// initialize uf and ucum if very first time here since ucum is cumulative (+=)
		if(timeorder==0) PLOOP(pl) uf[i][j][k][pl]=ucum[i][j][k][pl]=0.0;

		// NEED TO DEFINE Ui on other substeps besides the 0th one
		// find Ui(pi)
		PLOOP(pl) Ui[pl]=ui[i][j][k][pl];

		// dUriemann is volume averaged quantity
		flux2dUavg(i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
		PLOOP(pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is entirely consistent with point->averages

		// find uf==Uf and additional terms to ucum
		dUtoU(dUgeom, dUriemann, CUf, Cunew, Ui, uf[i][j][k], Uf, ucum[i][j][k]);
		//    flux2U(i,j,k,F1,F2,F3,dU,CUf,Cunew,Ui, uf[i][j][k],Uf,ucum[i][j][k]);

	}

	/////////////////////////////
	//
	// volume differentiate the conserved quantity
	//
	//////////////////////////////
	if(timeorder==numtimeorders-1){ // last call, so ucum is cooked and ready to eat!
		utoinvert=ucum;
		ubound=utoinvert;
	}
	else{ // otherwise still iterating on primitives
		utoinvert=uf;
		ubound=utoinvert;
	}

	// to debug ENO
//#if(DOENODEBUG)
	//debugeno_compute(pb,utoinvert,enodebugarray); //SASDEBUG -- OLD USE: now assign debugs inside reconstructeno code
//#endif


	// a2c_1 a2c_2 a2c_3
	if((avgscheme>3)&&(do_conserved_integration)){
		myupoint=upoint;
		avg2cen_interp(ENOCONSERVED, ENOAVG2CENTTYPE, pb, utoinvert, myupoint);  //SASMARK:  pb's for shock indicators should be defined on ALL grid, not only on ghost+active.  Maybe should use pi instead because define everywhere?
	}
	else{
		myupoint=utoinvert;
		//    FULLLOOP PLOOP(pl) upoint[i][j][k][pl]=utoinvert[i][j][k][pl];
	}




	CZLOOP {
		// set geometry for centered zone to be updated
		get_geometry(i, j, k, CENT, &geom);

		// setup initial guess for inversion
		// use pb since should be closer to pf
		PLOOP(pl) pf[i][j][k][pl] = pb[i][j][k][pl];

		// invert point U-> point p
		MYFUN(Utoprimgen(EVOLVEUTOPRIM, UEVOLVE, myupoint[i][j][k], &geom, pf[i][j][k]),"step_ch.c:advance()", "Utoprimgen", 1);

		//If using a high order scheme, need to choose whether to trust the point value
		if( (do_conserved_integration)&& ( avgscheme>3 )) {
			MYFUN(check_point_vs_average(timeorder, numtimeorders, pflag[i][j][k],pb[i][j][k],pf[i][j][k],myupoint[i][j][k],utoinvert[i][j][k],&geom),"advance.c:advance_finitevolume()", "check_point_vs_average()", 1);
		}


#if(DODISS||DODISSVSR)
		if(timeorder==numtimeorders-1){
		  // then see what entropy inversion would give
		  diss_compute(EVOLVEUTOPRIM,UEVOLVE,ucum[i][j][k],&geom,pf[i][j][k]);
		}
#endif


		// immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
		// SUPERGODMARK: Below should differentiate beteween whether want negative densities fixed or not, but right now fixup1zone() does all
		if((STEPOVERNEGU==0)||(STEPOVERNEGRHO==0)||(STEPOVERNEGRHOU==0)||(timeorder==numtimeorders-1)){
			MYFUN(fixup1zone(pf[i][j][k], &geom, finalstep),"advance.c:advance_finitevolume()", "fixup1zone()", 1);
		}
#endif
	}// end CZLOOP




#if(FIXUPZONES==FIXUPALLZONES)
	fixup(stage,pf,finalstep);
#endif  




	*ndt = defcon * 1. / (1. / ndt1 + 1. / ndt2 + 1. / ndt3);

#if(PRODUCTION==0)
	trifprintf( "2f");
#endif

	return (0);
}




// check whether point conserved quantity inverted successfully to point primitive.
//   if unsuccessful, then see if want to revert to average conserved quantity and invert that
//   if Uavg->p unsuccessful, then leave as failure
// if Upoint->p is good, then check if p from Upoint is much different than p from Uavg.  If so, limit change

// upoint only needed for diagnostics
int check_point_vs_average(int timeorder, int numtimeorders, int *lpflag, FTYPE *pb, FTYPE *pf, FTYPE *upoint, FTYPE *uavg, struct of_geom *ptrgeom)
{
  FTYPE pavg[NPR];  //atch for temporary storage of primitives obtained from inverting the averaged conserved quantities
  int invert_from_point_flag, invert_from_average_flag;
  FTYPE frac_avg_used;  //this is be used for flux interpolation limiting
  int pl;
  FTYPE limit_prim_correction( FTYPE fractional_difference_threshold, struct of_geom *geom, FTYPE *pin, FTYPE *pout );
  FTYPE fractional_diff( FTYPE a, FTYPE b );
  int is_convergence_failure;



  invert_from_point_flag = lpflag[FLAGUTOPRIMFAIL];


  if( 0 && debugfail >= 1 && (invert_from_point_flag == UTOPRIMFAILUNEG || invert_from_point_flag == UTOPRIMFAILRHONEG) ) {
    dualfprintf( fail_file, "t = %g, nstep = %ld, steppart = %d, i = %d, j = %d, rho = %21.15g, u = %21.15g, fracneg = %21.15g\n",
		 t, realnstep, steppart, ptrgeom->i + startpos[1], ptrgeom->j + startpos[2],
		 pf[RHO], pf[UU], (pf[RHO]>0)?(-pf[UU]/(pf[RHO]+DBL_MIN)):(-pf[RHO]/(pf[UU]+DBL_MIN)) );
  }


  //WHAT IF INTERNAL ENERGY BECOMES SLIGHTLY NEGATIVE?  WE STILL CAN DO THE LIMITING IN PRIM QUANTITIES! -- coorrected but check! -- SUPERSASMARK TODO atch
  if( LIMIT_AC_PRIM_FRAC_CHANGE &&
      (
       invert_from_point_flag == UTOPRIMNOFAIL || //atch added the below to still do the pt. vs. avg. check on primitives if the internal energy goes neg.
       ( (invert_from_point_flag==UTOPRIMFAILUNEG) && (0 != STEPOVERNEGU) ) || //intermediate substep with stepping over u < 0
       ( (invert_from_point_flag==UTOPRIMFAILRHONEG) && (0 != STEPOVERNEGRHO) ) //intermediate substep with stepping over rho < 0
       )
      ) {

    //make a copy of the initial guess so that not to modify the original pb's
    PLOOP(pl) pavg[pl] = pb[pl];
    //invert the average U -> "average" p
    MYFUN(Utoprimgen(EVOLVEUTOPRIM, UEVOLVE, uavg, ptrgeom, pavg),"step_ch.c:advance()", "Utoprimgen", 3);

    invert_from_average_flag = pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL];

    //Inversion from the average value succeeded or has a negative density or internal energy
    if( invert_from_average_flag == UTOPRIMNOFAIL  || invert_from_average_flag == UTOPRIMFAILUNEG ||  invert_from_average_flag == UTOPRIMFAILRHONEG ) {
      //Inversion from both from the point and the average values succeeded
      //checks if the states' gamma factors and densities are different by more than a certain fraction
      //and if different, modify the point values such that they are not further than by MAX_AC_PRIM_FRAC_CHANGE
      //from the average ones
      frac_avg_used = limit_prim_correction(MAX_AC_PRIM_FRAC_CHANGE, ptrgeom, pavg, pf);

      #if( DOENODEBUG )  //atch debug; note that use the location with pl =0 , interporflux = 0, & dir = 1 since limiting the change in prim quantities is not a per-component operation
         if(  frac_avg_used > 0.01 ) {
	   enodebugarray[ptrgeom->i][ptrgeom->j][ptrgeom->k][1-1][0][0][ENODEBUGPARAM_LIMCORR_PRIM]++;
         }
      #endif
      
      if( pf[UU] < 0.0 && timeorder == numtimeorders-1 ) {
        lpflag[FLAGUTOPRIMFAIL] = UTOPRIMFAILU2AVG2;
      }

      //lpflag[FLAGUTOPRIMFAIL] = invert_from_point_flag;  //unneeded since it is alrady = UTOPRIMNOFAIL

    } // end if both point and average did NOT fail
    else {
      //Inversion from the point values succeeded but that from the average failed:
      //retain the point value

      //set the inversion error flag to correspond to the inversion from the average value
      lpflag[FLAGUTOPRIMFAIL] = invert_from_point_flag;
      frac_avg_used = 0.0;  //used point value, i.e. zero fracion of the average value
    }
  }
  else if( INVERTFROMAVERAGEIFFAILED && (invert_from_point_flag!=UTOPRIMNOFAIL) ) {  //failure  //atch correct
    //inversion from the point value failed

    // if last substep -> revert to the average value, else if only negative densities then allow for substep.  If other type of failures, then never allow and revert to Uavg


    // GODMARK:
    // CHECK IF INVERSION FAILED.  IF FAILED, TRY TO USE uavg:
    // invert U->p

    // GODMARK: decide whether best to revert to average or not

    //=1 if it is a non-convergence failure; =0 if it is an occurrence of a negative density
      is_convergence_failure =  invert_from_point_flag!=UTOPRIMFAILUNEG &&
	invert_from_point_flag!=UTOPRIMFAILRHONEG &&
	invert_from_point_flag!=UTOPRIMFAILRHOUNEG;

      if((timeorder==numtimeorders-1 /*&& (1 == is_nonconvergence_failure || FIXUPZONES == FIXUPNOZONES)*/ ) ||   // last substep, then DO invert from average IF no fixups or non-convergence failure (ADT) SASMARKx

	 ( 1 == is_convergence_failure ) || //non-convergence (no solution in primitive quantities) error, so have to fix it up

	 ( (timeorder<numtimeorders-1) && (
					   ((invert_from_point_flag==UTOPRIMFAILUNEG) && (0 == STEPOVERNEGU))||
					   ((invert_from_point_flag==UTOPRIMFAILRHONEG) && (0 == STEPOVERNEGRHO))||
					   ((invert_from_point_flag==UTOPRIMFAILRHOUNEG) && (0 == STEPOVERNEGRHOU))
					   )
	   )  //intermediate substep with no stepping over u < 0, rho<0, or both <0
	 ) {
	if(debugfail >= 1) {
	  dualfprintf( fail_file, "Inversion from the point value failed.  Using the inversion from the average value.\n" );
	}

	//make a copy of the initial guess so that not to modify the original pb's
	PLOOP(pl) pf[pl] = pb[pl];
	//invert the average U -> "average" p
	MYFUN(Utoprimgen(EVOLVEUTOPRIM, UEVOLVE, uavg, ptrgeom, pf),"step_ch.c:advance()", "Utoprimgen", 3);
	//      invert_from_average_flag = lpflag[FLAGUTOPRIMFAIL];


	//Have the results from the inversion from the average value -- copy the result over
	//      PLOOP(pl) pf[pl] = pavg[pl];
	//      lpflag[FLAGUTOPRIMFAIL] = invert_from_average_flag;

	//old code:
	//MYFUN(Utoprimgen(EVOLVEUTOPRIM, UEVOLVE, avg, &geom, pf),"step_ch.c:advance()", "Utoprimgen", 2);

	frac_avg_used = 1.0; //reverted to the average value

      }
      else {
	frac_avg_used = 0.0; //used the point value
      }
  }

  return(0);
}





#define COMPARE_GAMMA 0

//If density or gamma-factors are different by more than fractional_difference_threshold for states pin & pout, 
//if different -- correct pout such that it is not more than fractional_difference_threshold away from pin.
FTYPE limit_prim_correction( FTYPE fractional_difference_threshold, struct of_geom *geom, FTYPE *pin, FTYPE *pout )
{
	FTYPE fractional_diff( FTYPE a, FTYPE b );
	FTYPE gammain = 0.0, gammaout = 0.0;
	FTYPE frac_start, frac_end, frac_diff;
	FTYPE fraction_input_value;
	FTYPE comovingenergyin, comovingenergyout;
	int pl;

#if( COMPARE_GAMMA ) 
	gamma_calc( pin, geom, &gammain );
	gamma_calc( pout, geom, &gammaout );
#endif

	comovingenergyin = pin[RHO] + pin[UU];  //SASMARKx: HD
	comovingenergyout = pout[RHO] + pout[UU];  //SASMARKx: HD
	
#if( COMPARE_GAMMA ) 
	frac_diff = MAX( fractional_diff(gammain, gammaout), 
		fractional_diff( comovingenergyin, comovingenergyout ) );
#else
	frac_diff = fractional_diff( comovingenergyin, comovingenergyout );
#endif

	//fractional difference at which the reduction to the input value starts
	frac_start = 0.5 * fractional_difference_threshold;

	//fractional difference after which only the input value is used
	frac_end = fractional_difference_threshold;

	//the fraction of the input value used in the output; increases linearly from 0 to 1 for fractions going from frac_start to frac_end
	fraction_input_value = MAX( 0., MIN(1., (frac_diff - frac_start)/(frac_end - frac_start) ) );

	if( 0.0 != fraction_input_value ){
		//states are too different: reverted to primitives that correspond to average conserved quantities because trust them more than point values
		dualfprintf( fail_file, "States are too different, using %3d%% of the average values: i = %d, j = %d, k = %d, nstep = %ld, steppart = %d, t = %21.15g\n", 
			(int)(100. * fraction_input_value), geom->i, geom->j, geom->k, nstep, steppart, t );
		if( debugfail >= 2 ){
			dualfprintf( fail_file, "Prim. pt. value (gamma, rho, u): " );
			dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n",  gammaout, pout[RHO], pout[UU] );
			dualfprintf( fail_file, "Prim. avg value (gamma, rho, u): " );
			dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n", gammain, pin[RHO], pin[UU] );
			dualfprintf( fail_file, "Frac. difference(ganna, rho, u): " );
			dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n", 
				fractional_diff(gammain, gammaout),
				fractional_diff(pin[RHO], pout[RHO]),  
				fractional_diff(pin[UU], pout[UU])
				);
		}
	}

	PLOOP(pl) {
		pout[pl] = fraction_input_value * pin[pl] + (1. - fraction_input_value) * pout[pl];
	}

	return( fraction_input_value );
}



//Returns the fractional difference between a & b
FTYPE fractional_diff( FTYPE a, FTYPE b )
{
	FTYPE frac_diff;

	frac_diff = 2. * fabs( a - b ) / ( fabs(a) + fabs(b) + DBL_MIN );

	return( frac_diff );

}

























// get dUavg
void flux2dUavg(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],FTYPE *dU1avg,FTYPE *dU2avg,FTYPE *dU3avg)
{
	FTYPE idx1,idx2,idx3;
	int pl;

#if(VOLUMEDIFF==0)
	idx1=1.0/dx[RR];
	idx2=1.0/dx[TH];
	idx3=1.0/dx[PH];
#else
	idx1=idxvol[i][j][k][RR];
	idx2=idxvol[i][j][k][TH];
	idx3=idxvol[i][j][k][PH];
#endif

	// initialize for simplicity so don't have to do it conditionally on N>1
	PLOOP(pl){
		dU1avg[pl]=0;
		dU2avg[pl]=0;
		dU3avg[pl]=0;
	}


	if(FLUXB==FLUXCD){ // don't use volume reg. since differencing is large
		for(pl=0;pl<B1;pl++){
#if(N1>1)
			dU1avg[pl]=(
				- (F1[ip1][j][k][pl] - F1[i][j][k][pl]) *idx1
				);
#endif
#if(N2>1)
			dU2avg[pl]=(
				- (F2[i][jp1][k][pl] - F2[i][j][k][pl]) *idx2
				);
#endif

#if(N3>1)
			dU3avg[pl]=(
				- (F3[i][j][kp1][pl] - F3[i][j][k][pl]) *idx3
				);
#endif

		}

		// simple version that assumes Fi[Bi] is set to 0 in flux.c for FLUXCD (which it is currently)
		for(pl=B1;pl<=B3;pl++){
#if(N1>1)
			dU1avg[pl]=(
				- (F1[ip1][j][k][pl] - F1[im1][j][k][pl]) *idx1
				);
#endif
#if(N2>1)
			dU2avg[pl]=(
				- (F2[i][jp1][k][pl] - F2[i][jm1][k][pl]) *idx2
				);
#endif
#if(N3>1)
			dU3avg[pl]=(
				- (F3[i][j][kp1][pl] - F3[i][j][km1][pl]) *idx3
				);
#endif
		}

		/*
		// old version
		// FIELDS are special.  The 0's would have come automatically, but spacing is twice the normal size.
		// the double spacing is accounted for in fluxct().
		pl=B1;
		#if(N1>1)
		dU1avg[pl]=(
		0.0
		);
		#endif
		#if(N2>1)
		dU2avg[pl]=(
		- (F2[i][jp1][k][pl] - F2[i][jm1][k][pl]) *idx2
		);
		#endif
		#if(N3>1)
		dU3avg[pl]=(
		- (F3[i][j][kp1][pl] - F3[i][j][km1][pl]) *idx3
		);
		#endif

		pl=B2;
		#if(N1>1)
		dU1avg[pl]=(
		- (F1[ip1][j][k][pl] - F1[im1][j][k][pl]) *idx1
		);
		#endif
		#if(N2>1)
		dU2avg[pl]=(
		0.0
		);
		#endif
		#if(N3>1)
		dU3avg[pl]=(
		- (F3[i][j][kp1][pl] - F3[i][j][km1][pl]) *idx3
		);
		#endif


		pl=B3;
		#if(N1>1)
		dU1avg[pl]=(
		- (F1[ip1][j][k][pl] - F1[im1][j][k][pl]) *idx1
		);
		#endif
		#if(N2>1)
		dU2avg[pl]=(
		- (F2[i][jp1][k][pl] - F2[i][jm1][k][pl]) *idx2
		);
		#endif
		#if(N3>1)
		dU3avg[pl]=(
		0.0
		);
		#endif
		// end old version
		*/


		// rest of variables (if any) are normal
		for(pl=B3+1;pl<NPR;pl++){
#if(N1>1)
			dU1avg[pl]=(
				- (F1[ip1][j][k][pl] - F1[i][j][k][pl]) *idx1
				);
#endif
#if(N2>1)
			dU2avg[pl]=(
				- (F2[i][jp1][k][pl] - F2[i][j][k][pl]) *idx2
				);
#endif
#if(N3>1)
			dU3avg[pl]=(
				- (F3[i][j][kp1][pl] - F3[i][j][k][pl]) *idx3
				);
#endif

		}
	}

	else{


		// other (normal) FLUXB methods
		PLOOP(pl) {


#if(N1>1)
			dU1avg[pl]=(
				- (F1[ip1][j][k][pl] - F1[i][j][k][pl]) *idx1
				);
#endif
#if(N2>1)
			dU2avg[pl]=(
				- (F2[i][jp1][k][pl] - F2[i][j][k][pl]) *idx2
				);
#endif
#if(N3>1)
			dU3avg[pl]=(
				- (F3[i][j][kp1][pl] - F3[i][j][k][pl]) *idx3
				);
#endif
		}


	}




}






// convert point versions of U_i^{n} and dU -> U_i^{n+1} and other versions
void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *uf, FTYPE *Uf, FTYPE *ucum)
{
	int pl;
	// finally assign new Uf and ucum
	// store uf to avoid recomputing U(pf) used later as pb for advance()
	PLOOP(pl) 
		uf[pl]=Uf[pl]   = CUf[0]*Ui[pl] + CUf[1]*uf[pl] + CUf[2]*dt*(dUriemann[pl]+dUgeom[pl]);

	// how much of Ui, dU, and Uf to keep for final solution
	// ultimately ucum is actual solution used to find final pf
	PLOOP(pl) ucum[pl] += Cunew[0]*Ui[pl] + Cunew[1]*dt*(dUriemann[pl]+dUgeom[pl]) + Cunew[2]*Uf[pl];

#if(PRODUCTION==0)
	PLOOP(pl){
	  if(!isfinite(uf[pl])){
	    dualfprintf(fail_file,"dUtoU after: nan found for uf[%d]=%21.15g\n",pl,uf[pl]);
	    dualfprintf(fail_file,"pl=%d Ui=%21.15g dUriemann=%21.15g dugeom=%21.15g\n",pl,Ui[pl],dUriemann[pl],dUgeom[pl]);
	  }
	}
	PLOOP(pl){
	  if(!isfinite(ucum[pl])){
	    dualfprintf(fail_file,"dUtoU after: nan found for ucum[%d]=%21.15g\n",pl,ucum[pl]);
	    dualfprintf(fail_file,"pl=%d Ui=%21.15g dUriemann=%21.15g dugeom=%21.15g\n",pl,Ui[pl],dUriemann[pl],dUgeom[pl]);
	  }
	}
#endif

}






// find global dt.  Useful if needed not during evolution, such as at t=0 or for restarting the run if restarting finished run that has a generally smaller dt than should use (including possibly dt=0)
int set_dt(FTYPE (*prim)[N2M][N3M][NPR], SFTYPE *dt)
{
	struct of_state state;
	struct of_geom geom;
	int i,j,k;
	int dir,ignorecourant;
	FTYPE cmax1,cmin1;
	FTYPE cmax2,cmin2;
	FTYPE cmax3,cmin3;
	SFTYPE dtij;
	FTYPE ndt;

	ndt=BIG;
	//ZLOOP{
	FULLLOOP{ //SUPERSASMARK
		get_geometry(i, j, k, CENT, &geom);

		MYFUN(get_state(prim[i][j][k], &geom, &state),"restart.c:set_dt()", "get_state()", 1);

#if(N1>1)
		dir=1;
		MYFUN(vchar(prim[i][j][k], &state, dir, &geom, &cmax1, &cmin1,&ignorecourant),"restart.c:set_dt()", "vchar() dir=1", 1);
		dtij = cour * dx[dir] / MAX(fabs(cmax1),fabs(cmin1));
		if (dtij < ndt) ndt = dtij;
#endif

#if(N2>1)
		dir=2;
		MYFUN(vchar(prim[i][j][k], &state, dir, &geom, &cmax2, &cmin2,&ignorecourant),"restart.c:set_dt()", "vchar() dir=2", 1);
		dtij = cour * dx[dir] / MAX(fabs(cmax2),fabs(cmin2));
		if (dtij < ndt) ndt = dtij;
#endif

#if(N3>1)
		dir=3;
		MYFUN(vchar(prim[i][j][k], &state, dir, &geom, &cmax3, &cmin3,&ignorecourant),"restart.c:set_dt()", "vchar() dir=3", 1);
		dtij = cour * dx[dir] / MAX(fabs(cmax3),fabs(cmin3));
		if (dtij < ndt) ndt = dtij;
#endif

	}

	// find global minimum value of ndt over all cpus
	mpifmin(&ndt);
	*dt = ndt;
	// don't step beyond end of run
	if (t + *dt > tf) *dt = tf - t;

	return(0);
}




