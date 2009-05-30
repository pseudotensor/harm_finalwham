
#include "decs.h"
#include "reconstructeno.h"  //atch -- to ensure that the declarations are in sync with reconstructeno.h (stupid)

extern weno_weights_t *stencil_weights_array;  //declare the array where the WENO routines put the weights, used for communication with stencil reduction

static weno_weights_t a_weights_array_minimal[NBIGM];  //for storing the minimums of the weights for a one-dimensional array of values; to feed to stencil reduction
static weno_weights_t *weights_array_minimal = &(a_weights_array_minimal[NBIGBND]);

static weno_weights_t weno_weights_backup[NPR2INTERP][NBIGM]; //for storing the weights for all quantities


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// LINE METHODS
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


// flux_interp() is provided p2interp and returns pleft/pright
//
// |=interface
// i=zone center of ith zone
//
// |              |p2interpm(i) if interporflux==ENOFLUXSPLITTYPE
// |              |      p2interp(i)   |
// |         pl(i)|pr(i)    i          |
// |         Fl(i)|Fr(i)    i          |
// |         Ul(i)|Ur(i)    i          |
// |              |pleft(i)   pright(i)|
// |              |F(i)                |
//


void flux_interp(int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interpm)[N2M][N3M][NPR], FTYPE (*p2interpp)[N2M][N3M][NPR], FTYPE (*pleft)[N2M][N3M][NPR], FTYPE (*pright)[N2M][N3M][NPR])
{

	void slope_lim_linetype(int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prims_guess)[N2M][N3M][NPR], FTYPE (*p2interpm)[N2M][N3M][NPR],FTYPE (*p2interpp)[N2M][N3M][NPR], FTYPE (*pleft)[N2M][N3M][NPR], FTYPE (*pright)[N2M][N3M][NPR]);

	int Nvec[NDIM];
	int i, j, k, pl;

	//use ptemparray as a temp array
	FTYPE (*prim_goodlocation)[N2M][N3M][NPR] = ptemparray;


	//////////////////////////////
	//
	// Get primitive to determine shock indicator
	//
	//////////////////////////////

	if( avgscheme == WENO5FLAT || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ) {

		if(USEAVGPRIMITIVEFORWENOFLAT && interporflux!=ENOFLUXSPLITTYPE){
			if(dir==1){
				// need to average primitive in dir=1 direction so in same location
				// GODMARK: loop is too big, but data isn't accessed later in bad region
				LOOPINFP1dir23full PLOOP(pl) prim_goodlocation[i][j][k][pl] = 0.5*(primreal[i][j][k][pl]+primreal[i-1][j][k][pl]);
			}
			else if(dir==2){
				LOOPINFP1dir13full PLOOP(pl) prim_goodlocation[i][j][k][pl] = 0.5*(primreal[i][j][k][pl]+primreal[i][j-1][k][pl]);
			}
			else if(dir==3){
				LOOPINFP1dir12full PLOOP(pl) prim_goodlocation[i][j][k][pl] = 0.5*(primreal[i][j][k][pl]+primreal[i][j][k-1][pl]);
			}
			else{
				dualfprintf(fail_file,"No such dir=%d in flux_interp()\n",dir);
				myexit(100);
			}
		}
		else{
			// flux split is zone centered
			//      prim_goodlocation=pr; // could also loop as above if want to protect pr
			FULLLOOP PLOOP(pl) prim_goodlocation[i][j][k][pl] = primreal[i][j][k][pl];
		}
	}
	else {
		//do not compute the primitive quantities that correspond to conserved ones because they will not be used; the case with NULL should be set up to be handled properly
		prim_goodlocation = NULL;
	}


	//////////////////////////////
	//
	// Perform interpolation
	//
	//////////////////////////////

	slope_lim_linetype(ENOFLUX, interporflux, dir, idel, jdel, kdel, prim_goodlocation, p2interpm, p2interpp, pleft, pright);

}


// convert avg->cent or cent->avg
// whichavg2cen = ENOCENT2AVGTYPE or ENOAVG2CENTTYPE
int avg2cen_interp(int whichquantity, int whichavg2cen, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR])
{
	//Utorpimgen is already declared in global.h

	void slope_lim_linetype(int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prim_goodlocation)[N2M][N3M][NPR], FTYPE (*p2interpm)[N2M][N3M][NPR],FTYPE (*p2interpp)[N2M][N3M][NPR], FTYPE (*pleft)[N2M][N3M][NPR], FTYPE (*pright)[N2M][N3M][NPR]);
	int dir;
	int idel,jdel,kdel;
	int Nvec[NDIM];
	struct of_geom geom;
	int i, j, k, pl;
	FTYPE (*current_in)[N2M][N3M][NPR]; //atch
	FTYPE (*current_out)[N2M][N3M][NPR]; //atch
	int numdirs,dirlist[NDIM],dirlisti;
	int numperms;
	long long itemp,permi;
	int get_dirs_list(int *numdirs, int *dirlist);
	int get_compdimen(int *numdirs, long long *itemp);
	int dirlist3d(long long itemp, int *dirlist);
	int dirlist2d(long long itemp, int *dirlist);
	int dirlist1d(long long itemp, int *dirlist);

	//use ptemparray as a temp array
	FTYPE (*prim_goodlocation)[N2M][N3M][NPR] = ptemparray;


	Nvec[1]=N1;
	Nvec[2]=N2;
	Nvec[3]=N3;


	//////////////////////////////
	//
	// Get primitive to determine shock indicator
	//
	//////////////////////////////

	if( avgscheme == WENO5FLAT || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ) {

		if(USEPRIMITIVEFROMAVGCONSERVED){

			FULLLOOP{
				get_geometry(i, j, k, CENT, &geom);
				// set guess
				PLOOP(pl) prim_goodlocation[i][j][k][pl]=primreal[i][j][k][pl];
				MYFUN(Utoprimgen(EVOLVEUTOPRIM, UEVOLVE, in[i][j][k], &geom, prim_goodlocation[i][j][k]),"interpline.c:avg2cen_interp()", "Utoprimgen", 1);
				// if problem with inversion, then reduce to using primreal
				if(pflag[i][j][k][FLAGUTOPRIMFAIL]) PLOOP(pl) prim_goodlocation[i][j][k][pl]=primreal[i][j][k][pl];
				// pflag will be reset by real inversion routine in advance.c before used elsewhere
			}

		}
		else{
			// then just use primreal, which is offset in time when used on updated average conserved quantity, but spatially in the right place
			// probably very good indicator of the future behavior

			if( primreal == NULL ) {
				//no primitive quantities values supplied -> obtain them by treating input as conserved quantities and inverting them
				dualfprintf( fail_file, "interpline.c: avg2cen_interp: WENO5FLAT requires supplying the primitive quantities\n" );
				myexit( 1 );
			}

			FULLLOOP  PLOOP(pl) prim_goodlocation[i][j][k][pl] = primreal[i][j][k][pl];
			//      prim_goodlocation = primreal; // could just assign pointer
		}
	}
	else {
		//do not compute the primitive quantities that correspond to conserved ones because they will not be used; the case with NULL should be set up to be handled properly
		prim_goodlocation = NULL;
	}


	//////////////////////////////
	//
	// Perform interpolation
	//
	//////////////////////////////

	// get number of dimensions (ignore itemp)
	get_compdimen(&numdirs,&itemp);

	// choose quasistrang or do this process if only 1 direction since no splitting
	if( (A2CDIMENSPLIT==QUASISTRANG)||(numdirs==1) ){
		// Quasi-strang splitting
		// Only applicable for 2D and 3D, where compared to UNSPLIT method this is 2X faster in 2D and 6X faster in 3D

		// get # of dimension and the directions to loop over
		get_dirs_list(&numdirs, dirlist);

		// setup input array so don't have to assume user sends in=out
		current_in = in; 
		current_out = out;

		// loop over dimensions/directions
		for(dirlisti=1;dirlisti<=numdirs;dirlisti++){

			// get direction to integrate
			dir=dirlist[dirlisti];
			//    dualfprintf(fail_file,"steppart=%d nstep=%d : dirlisti=%d dir=%d\n",steppart,nstep,dirlisti,dir);

			// get loop details
			// GODMARK: determine idel,jdel,kdel
			idel = fluxloop[dir][FIDEL];
			jdel = fluxloop[dir][FJDEL];
			kdel = fluxloop[dir][FKDEL];

			// a2c or c2a
			slope_lim_linetype(whichquantity, whichavg2cen, dir, idel, jdel, kdel, prim_goodlocation, current_in, NULL, current_out, NULL);

			// apply deaveraging one after another on the same array (out)
			current_in = current_out;
		}
	}
	else if(A2CDIMENSPLIT==UNSPLIT){
		// UNSPLIT a2c method (2X more expensive in 2D and 6X more expensive in 3D)

		if(numdirs==3){
			// then have 6 permutations to average out
			numperms=6;
		}
		else if(numdirs==2){
			// then have 2 permutations to average out
			numperms=2;
		}
		else numperms=1;


		FULLLOOP PLOOP(pl){
			// need to store input in case in=out (i.e. points to same memory)
			a2cin[i][j][k][pl] = in[i][j][k][pl];
			// inialize output (if in=out, then ok since just above stored in)
			out[i][j][k][pl] = 0.0;
		}


		// loop over permutations
		for(permi=0;permi<numperms;permi++){

			// get permutation
			if(numdirs==3) dirlist3d(permi,dirlist);
			else if(numdirs==2) dirlist2d(permi,dirlist);
			else if(numdirs==1) dirlist1d(permi,dirlist);

			// setup input array so don't have to assume user sends in=out
			current_in = a2cin;
			// output array is the temporary storage for this permutation
			current_out = a2cout;

#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
			//tell to SWENO that should look at other directions when performing reduction
			do_weno_lower_order_fraction = 1;
			FULLLOOP PLOOP(pl) {  //atch zero out the fractions so that initially use WENO5 everywhere
				weno_lower_order_fraction[i][j][k][pl] = 0.0;
			}
#endif

			// loop over dimensions/directions
			for(dirlisti=1;dirlisti<=numdirs;dirlisti++){

				// get direction to integrate
				dir=dirlist[dirlisti];
				//    dualfprintf(fail_file,"steppart=%d nstep=%d : dirlisti=%d dir=%d\n",steppart,nstep,dirlisti,dir);

				// get loop details
				// GODMARK: determine idel,jdel,kdel
				idel = fluxloop[dir][FIDEL];
				jdel = fluxloop[dir][FJDEL];
				kdel = fluxloop[dir][FKDEL];

				// a2c or c2a
				//atch here every next interpolation will use the global value of whether we have reduced while interpolating in 
				//the other direction, and if we did reduce then, then it will reduce now in the same way
				slope_lim_linetype(whichquantity, whichavg2cen, dir, idel, jdel, kdel, prim_goodlocation, current_in, NULL, current_out, NULL);

				// apply de/averaging one after another on the a2cout array
				// and don't overwrite the a2cin array since needed for every permutation
				current_in = current_out;

			}// end loop over dirlisti

#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
			do_weno_lower_order_fraction = 0;  //do not look at the other directions when performing reduction below this point
#endif

			// add this permutation to out (current_out=a2cout)
			FULLLOOP PLOOP(pl) out[i][j][k][pl] += current_out[i][j][k][pl];

		}// end loop over permutations

		// normalize the permutations
		FULLLOOP PLOOP(pl) out[i][j][k][pl] /= ( (FTYPE)numperms ) ;

	}// end UNSPLIT method



	return( 0 );
}



// Quasi-Strang Splitting (alternate order of dimension) of ordering of c2a integration of flux over 2 orthogonal directions
int get_dirs_list(int *numdirs, int *dirlist)
{
	long long itemp;
	int get_compdimen(int *numdirs, long long *itemp);
	int dirlist3d(long long itemp, int *dirlist);
	int dirlist2d(long long itemp, int *dirlist);
	int dirlist1d(long long itemp, int *dirlist);

	get_compdimen(numdirs,&itemp);

	if(*numdirs==3) dirlist3d(itemp,dirlist);
	else if(*numdirs==2) dirlist2d(itemp,dirlist);
	else if(*numdirs==1) dirlist1d(itemp,dirlist);

	return(0);
}



// get number of computational directions and the iterator unique to that dimension for choosing number of permutations
int get_compdimen(int *numdirs, long long *itemp)
{

	if(N1NOT1&&N2NOT1&&N3NOT1){
		*itemp=(nstep*(long long)numstepparts+(long long)steppart)%6;
		*numdirs=3;
	}
	else if( (N1NOT1&&N2NOT1)||(N1NOT1&&N3NOT1)||(N2NOT1&&N3NOT1) ){
		*itemp=(nstep*(long long)numstepparts+(long long)steppart)%2;
		*numdirs=2;
	}
	else{
		*itemp=(nstep*(long long)numstepparts+(long long)steppart)%1;
		*numdirs=1;
	}

	return(0);

}




// sequence of directions to interpolate in order of 1,2,3
int dirlist3d(long long itemp, int *dirlist)
{
	// 3D, number of permutations per call is: 2,2,1,2,2,1
	// normal ZEUS code has 1,2,1,2,1,1, so we are better mixed
	if(itemp==0){
		dirlist[1]=1;
		dirlist[2]=2;
		dirlist[3]=3;
	}
	else if(itemp==1){
		dirlist[1]=3;
		dirlist[2]=1;
		dirlist[3]=2;
	}
	else if(itemp==2){
		dirlist[1]=2;
		dirlist[2]=3;
		dirlist[3]=1;
	}
	else if(itemp==3){
		dirlist[1]=3;
		dirlist[2]=2;
		dirlist[3]=1;
	}
	else if(itemp==4){
		dirlist[1]=1;
		dirlist[2]=3;
		dirlist[3]=2;
	}
	else if(itemp==5){
		dirlist[1]=2;
		dirlist[2]=1;
		dirlist[3]=3;
	}

	return(0);
}


// sequence of directions to interpolate in order of 1,2
int dirlist2d(long long itemp, int *dirlist)
{
	if(N1NOT1&&N2NOT1){
		// 1 2
		if(itemp==0){
			dirlist[1]=1;
			dirlist[2]=2;
		}
		else if(itemp==1){
			dirlist[1]=2;
			dirlist[2]=1;
		}
	}
	else if(N1NOT1&&N3NOT1){
		// 1 3
		if(itemp==0){
			dirlist[1]=1;
			dirlist[2]=3;
		}
		else if(itemp==1){
			dirlist[1]=3;
			dirlist[2]=1;
		}
	}
	else if(N2NOT1&&N3NOT1){
		// 2 3
		if(itemp==0){
			dirlist[1]=2;
			dirlist[2]=3;
		}
		else if(itemp==1){
			dirlist[1]=3;
			dirlist[2]=2;
		}
	}

	return(0);
}



// sequence of directions to interpolate in order of 1
int dirlist1d(long long itemp, int *dirlist)
{

	if(N1NOT1){
		// 1
		dirlist[1]=1;
	}
	else if(N2NOT1){
		// 2
		dirlist[1]=2;
	}
	else if(N3NOT1){
		// 3
		dirlist[1]=3;
	}

	return(0);

}


#define NUMDFS 3

#define NUMMONOINDICATORS 3

// NPR2INTERP always larger than NPR, so can use one memory space for both c2e and others
#define LINETYPEDEFINES   int i,j,k;\
	extern int choose_limiter(int i, int j, int k, int pl);\
	FTYPE a_yin[NPR2INTERP][2][NBIGM];\
	FTYPE a_yprim[NPR2INTERP][2][NBIGM];\
	FTYPE a_yout[NPR2INTERP][2][NBIGM];\
	FTYPE a_df[NPR2INTERP][NUMDFS][NBIGM];\
	FTYPE a_drho[NUMDFS][NBIGM];\
	FTYPE a_dP[NUMDFS][NBIGM];\
	int a_shift[NBIGM];\
	FTYPE a_monoindicator[NPR2INTERP][NUMMONOINDICATORS][NBIGM];\
	int a_minorder[NBIGM];\
	int a_maxorder[NBIGM];\
	FTYPE a_shockindicator[NBIGM];\
	FTYPE (*yin)[2][NBIGM];\
	FTYPE (*yprim)[2][NBIGM];\
	FTYPE (*yout)[2][NBIGM];\
	FTYPE (*df)[NUMDFS][NBIGM];\
	FTYPE (*drho)[NBIGM];\
	FTYPE (*dP)[NBIGM];\
	int *shift;\
	FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM];\
	int *minorder,*maxorder;\
	FTYPE *shockindicator;\
	int reallim;\
	int bs,be,ps,pe;\
	int di,dj,dk;\
	int is,ie,js,je,ks,ke;\
	int preforder, whichreduce;\
	int pl;\
	int recontype;\
	int counter;\
	void set_preforder(int interporflux, int *preforder, int*whichreduce);\
	int get_recon_type(int interporflux);\
	extern void compute_monotonicity_line(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*df)[NBIGM],  FTYPE (*monoindicator)[NBIGM] , FTYPE *yin, FTYPE (*yout)[NBIGM]);\
	void get_df_line_gen_new(int numprims, int interporflux, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (**drhoptr)[NBIGM], FTYPE (**dPptr)[NBIGM]);\
	void set_interp_loop_gen(int dir, int interporflux, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);\
	void get_1d_line(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interpm)[N2M][N3M][NPR],FTYPE (*p2interpp)[N2M][N3M][NPR], FTYPE (*yin)[NBIGM]);\
	int get_shock_indicator(int dir, int interporflux, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yin)[2][NBIGM], FTYPE *shockindicator);\
	void causal_shift_order(int dir, int interporflux, int preforder, int bs, int ps, int pe, int be,  int i, int j, int k,  int idel, int jdel, int kdel, int *shift, int *minorder, int *maxorder);\
	void pass_1d_line(int whichquantity, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM])



// shift pointers to account for boundary zones (similar to as in set_arrays.c)
#define LINETYPESHIFTS   {yin =(FTYPE (*)[2][NBIGM]) (&(a_yin[0][0][NBIGBND]));\
	yprim =(FTYPE (*)[2][NBIGM]) (&(a_yprim[0][0][NBIGBND]));\
	yout=(FTYPE (*)[2][NBIGM]) (&(a_yout[0][0][NBIGBND]));\
	df=(FTYPE (*)[NUMDFS][NBIGM]) (&(a_df[0][0][NBIGBND]));\
	drho=(FTYPE (*)[NBIGM]) (&(a_drho[0][NBIGBND]));\
	dP=(FTYPE (*)[NBIGM]) (&(a_dP[0][NBIGBND]));\
	shift=(int (*)) (&(a_shift[NBIGBND]));\
	monoindicator=(FTYPE (*)[NUMMONOINDICATORS][NBIGM]) (&(a_monoindicator[0][0][NBIGBND]));\
	minorder=(int (*)) (&(a_minorder[NBIGBND]));\
	maxorder=(int (*)) (&(a_maxorder[NBIGBND]));\
	shockindicator=(FTYPE (*)) (&(a_shockindicator[NBIGBND]));}





// very similar to slope_lim_linetype, but different number of quantities to handle
void slope_lim_linetype_c2e(int numprims, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP])
{


	LINETYPEDEFINES;

	void get_1d_line_c2e(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP],FTYPE *yin);
	// below for shock indicator used on quantities with NPR elements
	void assign_eno_result_c2e(int interporflux, int pl, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[NBIGM], FTYPE (*result0)[N2M][N3M][NPR2INTERP], FTYPE (*result1)[N2M][N3M][NPR2INTERP]);


	LINETYPESHIFTS;

	// RESIDUAL NOTES:
	// DEATHMARK
	// yin is of size yin[bf-bs+1] +   yin[bs] is starting point for yin data.  Size:  yout[2][pf-ps+1]  yout[0/1][0] is first data point
	// shift is same size as yout.  It indicates the location of the center of the stencil w.r.t. the point of interest (i).  So a shift if 2 would mean center of stencil is at i+2.
	// primitives of some form are interpolated here and their dimensions are in standard [i][j][k] form.


	/////////////////
	//
	// Define loop over starting positions and range of loop for each starting position
	//
	/////////////////


	set_interp_loop_gen(dir, interporflux, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk, &bs, &ps, &pe, &be);

	///////////////////
	//
	// get reconstruction type (c2e, a2c, c2a)
	//
	////////////////////

	recontype=get_recon_type(interporflux);

	///////////////////
	//
	// set preferred order
	//
	////////////////////

	set_preforder(interporflux, &preforder, &whichreduce);



	////////////////
	//
	// LOOP OVER STARTING POSITIONS: this loop is over starting positions only
	//
	///////////////
	SUPERGENLOOP(i,j,k,is,ie,js,je,ks,ke,di,dj,dk){



		///////////////
		//
		// Figure out shifting of stencil and order of stencil to more closely match with causality.
		//  If no stencil shifting, then set shift and min/max order to default
		//
		///////////////
		causal_shift_order(dir, interporflux, preforder, bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, shift, minorder, maxorder);


		// subloop loops from starting position to end of that line


		///////////////////
		//
		// GET 1D LINES
		//
		////////////////////

		if(whichreduce == WENO_REDUCE_TYPE_PPM || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ){
			if(VARTOINTERP==PRIMTOINTERP){
				// means primreal=p2interp
				// then get all primitives and assign to normal yin
				PLOOP(pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, primreal, NULL, yin[pl]);
				yprim=yin;
			}
			else{
				// means primreal!=p2interp
				// then need to get separate p2interp line
				PLOOPINTERP(pl) get_1d_line_c2e(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, p2interp, yin[pl][0]);

				// get prim line
				if(whichreduce == WENO_REDUCE_TYPE_PPM || SHOCKINDICATOR ){
					// then need to get primitives separately from interpolated quantities
					PLOOP(pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[pl]);
				}
				else if(CONTACTINDICATOR||COMPUTEDRHODP){
					// then only need RHO and UU
					get_1d_line(dir, interporflux, RHO, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[RHO]);
					get_1d_line(dir, interporflux,  UU, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[UU]);
				}
			}
		}
		else{
			// then no ppm reduce or contact indicator, so just get normal line
			///////////////////
			//
			// 1D GET LINE (interpolated quantity is never the same as the prims4shock quantity since not doing c2e here)
			//
			////////////////////
			PLOOPINTERP(pl) get_1d_line_c2e(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, p2interp, yin[pl][0]);
		}


		///////////////////
		//
		// 1D GET SHOCK INDICATOR
		//
		////////////////////
		if(whichreduce == WENO_REDUCE_TYPE_PPM || SHOCKINDICATOR ){
			//use primitive values that correspond to the quantities being interpolated
			get_shock_indicator(dir, interporflux,  bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, yprim, shockindicator);
		}

		///////////////////
		//
		// get df's for contact indicator or in general
		//
		////////////////////
		get_df_line_gen_new(numprims,interporflux,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yprim, yin, df, &drho, &dP);


		///////////////
		//
		// Compute monotonicity indicator and if monotonic or rough then compute interface value
		//
		///////////////

		PLOOPINTERP(pl) compute_monotonicity_line(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], yin[pl][0], yout[pl]);

		///////////////
		//
		// PASS 1D LINE
		//
		///////////////

		PLOOPINTERP(pl) pass_1d_line( ENOPRIMITIVE, ALL_CALC, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], dP, monoindicator[pl], yprim, yin[pl], yout[pl]);


		//////////////////////////////////////////
		//
		// now have 1D line of point data corresponding to left and right interpolated values
		//
		// Assign result back to pleft/pright (or just pleft if one result)
		//
		/////////////////////////////////////////
		PLOOPINTERP(pl) assign_eno_result_c2e(recontype, pl, bs, ps, pe, be, i, j, k, idel, jdel, kdel, yout[pl], pleft, pright);


	} // done with main loop over starting points






}


// linetype assumed to return pleft/pright directly since linetype's are usually higher order
// this function used for primitive interpolation AND for flux interpolation
void slope_lim_linetype(int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interpm)[N2M][N3M][NPR], FTYPE (*p2interpp)[N2M][N3M][NPR], FTYPE (*pleft)[N2M][N3M][NPR], FTYPE (*pright)[N2M][N3M][NPR])
{

	void assign_eno_result(int interporflux, int pl, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[NBIGM], FTYPE (*result0)[N2M][N3M][NPR], FTYPE (*result1)[N2M][N3M][NPR]);
  void minimize_weno5_weights( int preforder, weno_weights_t *weno_weights_current, weno_weights_t *weno_weights_minimal, weno_weights_t *weno_weights_minimal_out );
	void copy_weno5_weights( int preforder, weno_weights_t *weno_weights_source, weno_weights_t *weno_weights_destination );
	void reset_weno5_weights( int preforder, weno_weights_t *weno_weights_to_set );
	void compute_lower_order_fraction_weno5_weights( int preforder, weno_weights_t *weno_weights );
	void rescale_weno5_weights( int preforder, FTYPE rescale_factor, weno_weights_t *weights_array_in, weno_weights_t *weights_array_out );
	extern int apply_bc_line(int doinverse, int iterglobal, int recontype, int bs, int be, FTYPE (*yin)[2][NBIGM],FTYPE (*yout)[2][NBIGM]);


	LINETYPEDEFINES;
	int numprims=NPR;
	int is_copy;
	int mypl,myi;
	int plstart;

	LINETYPESHIFTS;



	// RESIDUAL NOTES:
	// DEATHMARK
	// yin is of size yin[bf-bs+1] +   yin[bs] is starting point for yin data.  Size:  yout[2][pf-ps+1]  yout[0/1][0] is first data point
	// shift is same size as yout.  It indicates the location of the center of the stencil w.r.t. the point of interest (i).  So a shift if 2 would mean center of stencil is at i+2.
	// primitives of some form are interpolated here and their dimensions are in standard [i][j][k] form.


	/////////////////
	//
	// Define loop over starting positions and range of loop for each starting position
	//
	/////////////////

	set_interp_loop_gen(dir, interporflux, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk, &bs, &ps, &pe, &be);

	///////////////////
	//
	// get reconstruction type (c2e, a2c, c2a)
	//
	////////////////////

	recontype=get_recon_type(interporflux);

	///////////////////
	//
	// set preferred order
	//
	////////////////////

	set_preforder(interporflux, &preforder, &whichreduce);




	////////////////
	//
	// LOOP OVER STARTING POSITIONS: this loop is over starting positions only
	//
	///////////////
	SUPERGENLOOP(i,j,k,is,ie,js,je,ks,ke,di,dj,dk){


		///////////////
		//
		// Figure out shifting of stencil and order of stencil to more closely match with causality.
		//  If no stencil shifting, then set shift and min/max order to default
		//
		///////////////
		causal_shift_order(dir, interporflux, preforder, bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, shift, minorder, maxorder);



		///////////////////
		//
		// GET 1D LINES
		//
		////////////////////
		// subloop loops from starting position to end of that line

		if(whichreduce == WENO_REDUCE_TYPE_PPM || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ){
			if((interporflux==ENOINTERPTYPE)&&(VARTOINTERP==PRIMTOINTERP)){
				// means primreal=p2interpm
				// then get all primitives and assign to normal yin
				PLOOP(pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, primreal, NULL, yin[pl]);
				yprim=yin;
			}
			else{
				// means primreal!=p2interpm
				// then need to get separate p2interp line
				PLOOP(pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, p2interpm, p2interpp, yin[pl]);

				// get prim line
				if(whichreduce == WENO_REDUCE_TYPE_PPM || SHOCKINDICATOR ){
					// then need to get primitives separately from interpolated quantities
					PLOOP(pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[pl]);
				}
				else if(CONTACTINDICATOR||COMPUTEDRHODP){
					// then only need RHO and UU
					get_1d_line(dir, interporflux, RHO, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[RHO]);
					get_1d_line(dir, interporflux,  UU, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[UU]);
				}
			}
		}
		else{
			// then no ppm reduce or contact indicator
			///////////////////
			//
			// 1D GET LINE (interpolated quantity is never the same as the prims4shock quantity since not doing c2e here)
			//
			////////////////////
			PLOOP(pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, p2interpm, p2interpp, yin[pl]);
		}



		// user-defined modification of the line -- usually adjusting line data so interpolation can be higher-order
		apply_bc_line(0,iterglobal,recontype,bs,be,yin,NULL);
		
		///////////////////
		//
		// 1D GET SHOCK INDICATOR
		//
		////////////////////

		if(whichreduce == WENO_REDUCE_TYPE_PPM|| SHOCKINDICATOR ){
			//use primitive values that correspond to the quantities being interpolated
			get_shock_indicator(dir, interporflux,  bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, yprim, shockindicator);
		}

		///////////////////
		//
		// get df's for contact indicator or in general
		// GODMARK: not done for flux splitting method (i.e. assume only 1 input always)
		//
		////////////////////
		get_df_line_gen_new(numprims,interporflux,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yprim, yin, df, &drho, &dP);


		///////////////
		//
		// Compute monotonicity indicator and if monotonic or rough then compute interface value
		// GODMARK: not done for flux splitting method (i.e. assume only 1 input always)
		//
		///////////////

		PLOOP(pl) compute_monotonicity_line( recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], yin[pl][0], yout[pl]);


		///////////////
		//
		// PASS 1D LINE
		//
		///////////////

		if( (recontype == CVT_A2C || recontype == CVT_C2A) ){
			if( DO_SPLITA2C == NOSPLITA2C || whichquantity == ENOSOURCETERM ) { //do not do any weights minimization for source term integration
				//unsplit version, compute weights and reconstruction in one call
				PLOOP(pl) 
					pass_1d_line( whichquantity, ALL_CALC, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], dP, monoindicator[pl], yprim, yin[pl], yout[pl]);
				//SUPERSASMARK : insert a2c/c2a limiting code here
			}
			else if(	DO_SPLITA2C != NOSPLITA2C ) {
				//split version, compute weights in the 1st call and do the reconstruction in 2nd call
				reset_weno5_weights( preforder, weights_array_minimal );

				if( DO_SPLITA2C == MINIMIZE_ALL_WEIGHTS ) //a2c on conserved
				{
					PLOOP(pl){
						pass_1d_line( whichquantity, WEIGHT_CALC, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], dP, monoindicator[pl], yprim, yin[pl], yout[pl]);
						minimize_weno5_weights( preforder, stencil_weights_array, weights_array_minimal, weights_array_minimal );
					}
					//set the lower_order_fraction of the weights to match that of the sum of the weno5 weights
					compute_lower_order_fraction_weno5_weights( preforder, weights_array_minimal );

					//put the minimal weights back into the weno storage so that weno uses those weights; 1 stands for copy
					copy_weno5_weights( preforder, weights_array_minimal, stencil_weights_array );
					PLOOP(pl) {
						pass_1d_line( whichquantity, RECON_CALC, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], dP, monoindicator[pl], yprim, yin[pl], yout[pl]);
					}
					//SUPERSASMARK : insert a2c/c2a limiting code here
				}
				else if(  DO_SPLITA2C == ENERGY_CONTROLS_ALL_WEIGHTS ) {
					if( whichquantity == ENOFLUX ) {
						//for fluxes, use common weights computed for the flux of dir-momentum in the direction of dir (because that flux does not vanish even if the velocities are zero)
						plstart = dirglobal + U1 - 1; 
					}
					else if( whichquantity == ENOCONSERVED ){
						//for conserved quantities, use common weights computed for the energy (because in general total energy is unlikely to vanish)
						plstart = UU;
					}
					else {
						dualfprintf( fail_file, "slope_lim_linetype(): Integration of source terms should be handled in an unsplit way.\n" );
						myexit( 1 );
					}

					PLOOP(mypl) {
						pl = (mypl + plstart) % NPR;  //shift the start of the loop to pl = plstart
						pass_1d_line( whichquantity, WEIGHT_CALC, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], dP, monoindicator[pl], yprim, yin[pl], yout[pl]);
						if( pl == plstart ) {
							//store the weights for this quantity (energy or momentum depending on the interpolation type) as common weights;
							//they will be used for comparison to weights computed for other quantities
							minimize_weno5_weights( preforder, stencil_weights_array, weights_array_minimal, weights_array_minimal );
						}
						else {
							//otherwise, the common weights have already been precomputed, so modify the current weights to be smallest of the two:  the energy weights and twice the current weights
							minimize_weno5_weights( preforder, stencil_weights_array, weights_array_minimal, stencil_weights_array );
							//synchronise the lower_order_fraction with the normalization of the resulting weights
							compute_lower_order_fraction_weno5_weights( preforder, stencil_weights_array );
						}
						pass_1d_line( whichquantity, RECON_CALC, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], dP, monoindicator[pl], yprim, yin[pl], yout[pl]);
					}
					//SUPERSASMARK : insert a2c/c2a limiting code here
				}
			}
			else {
				dualfprintf( fail_file, "Unknown reconstruction split type: DO_SPLITA2C = %d\n", DO_SPLITA2C );
				myexit( 1 );
			}
		}
		else{
			//unsplit version, compute weights and reconstruction in one call
			PLOOP(pl) pass_1d_line( whichquantity, ALL_CALC, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], dP, monoindicator[pl], yprim, yin[pl], yout[pl]);
		}



		// user-defined de-modification of the line -- usually adjusting line data so interpolation can be higher-order
		apply_bc_line(1,iterglobal,recontype,bs,be,yin,yout);
		

		//////////////////////////////////////////
		//
		// now have 1D line of point data corresponding to left and right interpolated values
		//
		// Assign result back to pleft/pright (or just pleft if one result)
		//
		/////////////////////////////////////////
		PLOOP(pl) assign_eno_result(recontype, pl, bs, ps, pe, be, i, j, k, idel, jdel, kdel, yout[pl], pleft, pright);


	} // done with main loop over starting points


}


//NOTE: preforder should be set correctly; if it is not, the behaviour will be wrong.  
//      preforder is used because some of the weights at the boundary are not defined and 
//      for those weights one cannot use the order in the summing procedures.
//      Basically, the order that comes with the weights is copied over but is ignored
//      for the purposes of the summing, etc. in these functions

void rescale_weno5_weights( int preforder, FTYPE rescale_factor, weno_weights_t *weights_array_in, weno_weights_t *weights_array_out )
{
	int weight_count, i;
	int order = (preforder + 1) / 2;

	for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
		for( weight_count = 0; weight_count < order; weight_count++ ) {
			weights_array_out[i].weights[weight_count] = rescale_factor * weights_array_in[i].weights[weight_count];  //rescale the weights
		}
		weights_array_out[i].lower_order_fraction = 0.0;
		weights_array_out[i].order = order;
		weights_array_out[i].len = 0;
	}
}



//resets the unoptimized weights to unity; lower_order_fraction to 0.
void reset_weno5_weights( int preforder, weno_weights_t *weno_weights_to_set )
{
	int weight_count, i;
	int order = (preforder + 1) / 2;

	for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
		for( weight_count = 0; weight_count < order; weight_count++ ) {
			weno_weights_to_set[i].weights[weight_count] = 1.0;  //set the weights to unity
		}
		weno_weights_to_set[i].lower_order_fraction = 0.0;
		weno_weights_to_set[i].order = order;
		weno_weights_to_set[i].len = 0;
	}
}

//copies weights from source to destination
void copy_weno5_weights( int preforder, weno_weights_t *weno_weights_source, weno_weights_t *weno_weights_destination )
{
	int weight_count, i;
	int order = (preforder + 1) / 2;

	//simply copy the weights over
	for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
		weno_weights_destination[i] = weno_weights_source[i];  //copy the structures as a whole (as single objects)
	}
}

//computes the sum of the weno5 weights and set the lower order fraction to 1 minus this value
void compute_lower_order_fraction_weno5_weights( int preforder, weno_weights_t *weno_weights )
{
	FTYPE weights_sum;
	int order = (preforder + 1) / 2;
	int weight_count, i;

	//simply copy the weights over
	for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
		//SASMARK order hard-cored
		weights_sum = get_sum_of_elements( order, weno_weights[i].weights );  //copy the structures as a whole (as single objects)
		weno_weights[i].lower_order_fraction = 1.0 - weights_sum;
	}
}

//compute the minimum of each weight in the weno_weights_current and weno_weights_minimal, and put these minimums into weno_weights_minimal_out
//the weno_weights_current's weights are multipulied by (1-lower_order_fraction) before comparison.
void minimize_weno5_weights( int preforder, weno_weights_t *weno_weights_current, weno_weights_t *weno_weights_minimal, weno_weights_t *weno_weights_minimal_out )
{
	int i, weight_count;
	int order = (preforder + 1) / 2;
	
	//keep track of the minimum value of each weight in weno_weights_minimal
	for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
		weno_weights_minimal_out[i].lower_order_fraction = -9;  //lower order fraction is not set
		weno_weights_minimal_out[i].order = weno_weights_current[i].order;
		weno_weights_minimal_out[i].len = weno_weights_current[i].len;
		//only keep track of the maximal lower order fraction but 
		//do not copy the order and len since they do not change from quantity to quantity; however, the weights ratios will change, but we don't use them after this point
		//keep track of the *unoptimized* minimal weights and put them into the destination
		for( weight_count = 0; weight_count < order; weight_count++ ) {
			weno_weights_minimal_out[i].weights[weight_count] = MIN( 
				weno_weights_current[i].weights[weight_count] * ( 1.0 - weno_weights_current[i].lower_order_fraction ), //rescale the weights by the weno5 fraction
				weno_weights_minimal[i].weights[weight_count] );
		}
	}
}

// get type of reconstruction to perform
int get_recon_type(int interporflux)
{

	// c2a stuff
	if(interporflux==ENOFLUXAVG1TYPE|| interporflux==ENOFLUXAVG2TYPE || interporflux==ENOFLUXAVG3TYPE || interporflux==ENOCENT2AVGTYPE) return(CVT_C2A);
	// c2e stuff
	else if(interporflux==ENOFLUXSPLITTYPE || interporflux==ENOINTERPTYPE) return(CVT_C2E);
	// a2c stuff
	else if(interporflux==ENOFLUXRECONTYPE || interporflux==ENOAVG2CENTTYPE) return(CVT_A2C);
	else{
		dualfprintf(fail_file,"No such interporflux=%d\n", interporflux);
		myexit(177);
	}

	return(-100);

}



// sets preferred order based upon scheme and size of stencil
void set_preforder(int interporflux, int *preforder, int*whichreduce)
{

	if(interporflux==ENOINTERPTYPE){
		//      reallim=choose_limiter(i,j,k,pl); // starting point chooses limiter type
		*preforder=interporder[lim]; // get order of scheme

		if( lim == WENO5FLAT ) {  //correct the order for WENO5FLAT: in this case need more points than the order for stencil reduction
			*whichreduce = WENO_REDUCE_TYPE_PPM;
			*preforder -= 2;
		}
		else if( lim == WENO5BND ) {  //correct the order for WENO5FLAT: in this case need more points than the order for stencil reduction
			*whichreduce = WENO_REDUCE_TYPE_DEFAULT;
			*preforder = 5;  //correct the order of the scheme because the number of points passed to it is larger than its order (the order is 5)
		}
		else {
			*whichreduce = WENO_REDUCE_TYPE_DEFAULT;
		}
	}
	else{
		*preforder=interporder[avgscheme];

		if( avgscheme == WENO5FLAT ) {  //correct the order for WENO5FLAT: in this case need more points than the order for stencil reduction
			*whichreduce = WENO_REDUCE_TYPE_PPM;
			*preforder -= 2;
		}
		else if( avgscheme == WENO5BND ) {  //correct the order for WENO5FLAT: in this case need more points than the order for stencil reduction
			*whichreduce = WENO_REDUCE_TYPE_DEFAULT;
			*preforder = 5;  //correct the order of the scheme because the number of points passed to it is larger than its order (the order is 5)
		}
		else {
			*whichreduce = WENO_REDUCE_TYPE_DEFAULT;
		}
	}


}


// just a wrapper
void set_interp_loop_gen(int dir, int interporflux, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be)
{
	void set_interp_loop(int dir, int interporflux, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);
	void set_interp_loop_expanded(int dir, int interporflux, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);

	// straight-forward average or de-average along dir and just one ghost layer (original method)
	if(DOENOFLUX!=ENOFINITEVOLUME){
		// GODMARK: later should convert ENOFLUXRECON method to have expanded ghost layer and ghost+active layer
		set_interp_loop(dir, interporflux, is, ie, js, je, ks, ke, di, dj, dk, bs, ps, pe, be);
	}
	else{
		// now assume active + active/ghost + ghost layers so only bound primitives during calculation
		set_interp_loop_expanded(dir, interporflux, is, ie, js, je, ks, ke, di, dj, dk, bs, ps, pe, be);
	}

}




// wrapper for compute_df_line and can be used by both c2e and other methods
void get_df_line_gen(int numprims, int interporflux, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (**drhoptr)[NBIGM], FTYPE (**dPptr)[NBIGM])
{
	int pl;
	int compute_df_line(int interporflux, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE *yin, FTYPE (*df)[NBIGM]);


#define NUMPRIMLOOP(pl) for(pl=0;pl<numprims;pl++)


	if(CONTACTINDICATOR||COMPUTEDRHODP){
		if((interporflux==ENOINTERPTYPE)&&(VARTOINTERP==PRIMTOINTERP)){
			// then don't need to get separate drho and dP

			// then get all df's
			NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][0], df[pl]);

			// then primitives are same as interpolated quantities
			// assign drho and dP.  This overwrites original pointer that pointed to independent memory for drho and dP
			*drhoptr=df[RHO];
			*dPptr=df[UU];
		}
		else{
			// then need to compute separate drho and dP
			///////////////
			// Compute differentials needed for contactindicator
			///////////////
			compute_df_line(interporflux,whichreduce,preforder, RHO, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yprim[RHO][0], *drhoptr);
			compute_df_line(interporflux,whichreduce,preforder, UU, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yprim[UU][0], *dPptr);

			// then get all p2interp df's
			NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][0], df[pl]);
			//	if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][1], dfp[pl]);
		}
	}
	else{
		// then don't need drho or dP at all, just get p2interp df's
		// then get all p2interp df's
		NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][0], df[pl]);
		//	if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][1], dfp[pl]);
	}
}







// wrapper for compute_df_line and can be used by both c2e and other methods
void get_df_line_gen_new(int numprims, int interporflux, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (**drhoptr)[NBIGM], FTYPE (**dPptr)[NBIGM])
{
	int pl;
	int compute_df_line(int interporflux, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE *yin, FTYPE (*df)[NBIGM]);
	int compute_df_line_new(int numprims, int interporflux, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*yin)[2][NBIGM], FTYPE *df);



	if(CONTACTINDICATOR||COMPUTEDRHODP){
		if((interporflux==ENOINTERPTYPE)&&(VARTOINTERP==PRIMTOINTERP)){
			// then don't need to get separate drho and dP


			NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][0], df[pl]);

			// then get all df's
			compute_df_line_new(numprims, interporflux,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin, (*dPptr)[0]);

			// then primitives are same as interpolated quantities
			// assign drho and dP.  This overwrites original pointer that pointed to independent memory for drho and dP
			*drhoptr=df[RHO];
			//      *dPptr=df[UU];
		}
		else{
			// then need to compute separate drho and dP
			///////////////
			// Compute differentials needed for contactindicator
			///////////////

			// then get all p2interp df's
			NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][0], df[pl]);
			//	if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][1], dfp[pl]);

			compute_df_line(interporflux,whichreduce,preforder, RHO, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yprim[RHO][0], *drhoptr);

			compute_df_line_new(numprims,interporflux,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yprim, (*dPptr)[0]);

		}
	}
	else{
		// then don't need drho or dP at all, just get p2interp df's
		// then get all p2interp df's
		NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][0], df[pl]);
		//	if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pl) compute_df_line(interporflux,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, yin[pl][1], dfp[pl]);
	}
}





// Setup loop over *starting positions* that and define the line of data corresponding to data needed by the given (dir && interporflux) scenario
// i.e. Define loop over starting positions and range of loop for each starting position
void set_interp_loop(int dir, int interporflux, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be)
{

	// determine range of outer loop and range to feed to eno scheme
	if(dir==1){
		*is=*ie=0; // anything so di=1 iterates, next subloop overwrites i actually useds

		if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE)){
			// then centered quantity, and need from -1 .. N and have ijkminmax[1][0/1] outer range of ghost zones
			*bs=ijkminmax[interporflux][1][0];
			*be=ijkminmax[interporflux][1][1];

			*ps=-1;
			*pe=N1;
		}
		else if(interporflux==ENOFLUXRECONTYPE){
			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][1][0]+1;
			*be=ijkminmax[interporflux][1][1];

			*ps=0;
			*pe=N1;
		}

		*js=INFULL2; // 0
		*je=OUTFULL2; // N2-1;

		*ks=INFULL3; //0;
		*ke=OUTFULL3; // N3-1;

		*di=1;
		*dj=1;
		*dk=1;
	}
	else if(dir==2){
		*is=INFULL1; //0;
		*ie=OUTFULL1; //N1-1;


		*js=*je=0; // anything so dj=1 iterates, next subloop overwrites j actually useds

		if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE)){
			// then centered quantity, and need from -1 .. N and have ijkminmax[1][0/1] outer range of ghost zones
			*bs=ijkminmax[interporflux][2][0];
			*be=ijkminmax[interporflux][2][1];

			*ps=-1;
			*pe=N2;
		}
		else if(interporflux==ENOFLUXRECONTYPE){
			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][2][0]+1;
			*be=ijkminmax[interporflux][2][1];

			*ps=0;
			*pe=N2;
		}

		*ks=INFULL3; // 0;
		*ke=OUTFULL3; // N3-1;

		*di=1;
		*dj=1;
		*dk=1;
	}
	else if(dir==3){
		*is=INFULL1; // 0;
		*ie=OUTFULL1; // N1-1;

		*js=INFULL2; // 0;
		*je=OUTFULL2; // N2-1;

		*ks=*ke=0; // anything so dk=1 iterates, next subloop overwrites k actually used


		if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE)){
			// then centered quantity, and need from -1 .. N and have ijkminmax[1][0/1] outer range of ghost zones
			*bs=ijkminmax[interporflux][3][0];
			*be=ijkminmax[interporflux][3][1];

			*ps=-1;
			*pe=N3;
		}
		else if(interporflux==ENOFLUXRECONTYPE){
			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][3][0]+1;
			*be=ijkminmax[interporflux][3][1];

			*ps=0;
			*pe=N3;
		}

		*di=1;
		*dj=1;
		*dk=1;

	}

}





// Setup loop over *starting positions* that and define the line of data corresponding to data needed by the given (dir && interporflux) scenario
// i.e. Define loop over starting positions and range of loop for each starting position
// This function is for any method using the expanded ghost+ghost/active+active layers (presently only DOENOFLUX==ENOFINITEVOLUME, but can be used for ENOFLUXRECON if recode that routine)
void set_interp_loop_expanded(int dir, int interporflux, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be)
{
	int dir_exception[NDIM];

	//  if((interporflux==ENOFLUXAVG1TYPE)||(interporflux==ENOFLUXAVG2TYPE)||(interporflux==ENOFLUXAVG3TYPE) ){

	if(
		((interporflux==ENOFLUXAVG1TYPE)&&(dir==1))||
		((interporflux==ENOFLUXAVG2TYPE)&&(dir==2))||
		((interporflux==ENOFLUXAVG3TYPE)&&(dir==3))
		){
			dualfprintf(fail_file,"No such method with interporflux=%d and dir=%d\n",interporflux,dir);
			myexit(1);
	}

	dir_exception[1] =  (interporflux==ENOFLUXAVG2TYPE) || (interporflux==ENOFLUXAVG3TYPE);
	dir_exception[2] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG3TYPE);
	dir_exception[3] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG2TYPE);

	// determine range of outer loop and range to feed to eno scheme
	if( ( (!dir_exception[1]) && (dir==1) ) || (interporflux==ENOFLUXAVG1TYPE) ){
		*is=*ie=0; // anything so di=1 iterates, next subloop overwrites i actually used

		// input and output at different location

		if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE)){
			// then centered quantity, and need from -1 .. N and have ijkminmax[1][0/1] outer range of ghost zones
			*bs=ijkminmax[interporflux][1][0];
			*be=ijkminmax[interporflux][1][1];

			*ps=Uconsloop[FIS]-1;
			*pe=Uconsloop[FIE]+1;

			*js=INFULL2;
			*je=OUTFULL2;

			*ks=INFULL3;
			*ke=OUTFULL3;
		}
		else if(interporflux==ENOFLUXRECONTYPE){
			// input and output at same location

			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][1][0]+1;
			*be=ijkminmax[interporflux][1][1];

			*ps=Uconsloop[FIS];
			*pe=Uconsloop[FIE]+1;

			*js=INFULL2;
			*je=OUTFULL2;

			*ks=INFULL3;
			*ke=OUTFULL3;

		}
		else if( interporflux==ENOAVG2CENTTYPE ){
			// input and output at same location

			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=Uconsloop[FIS];  
			*be=Uconsloop[FIE];

			*ps=0;
			*pe=N1 - 1;

			*js=INFULL2;
			*je=OUTFULL2;

			*ks=INFULL3;
			*ke=OUTFULL3;


		}
		else if( interporflux==ENOCENT2AVGTYPE ){
			*bs=ijkminmax[interporflux][1][0]+1;
			*be=ijkminmax[interporflux][1][1]-1;

			*ps=Uconsloop[FIS];
			*pe=Uconsloop[FIE];

			*js=INFULL2;
			*je=OUTFULL2;

			*ks=INFULL3;
			*ke=OUTFULL3;

		}
		else if(interporflux==ENOFLUXAVG1TYPE){

			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][1][0]+1;
			*be=ijkminmax[interporflux][1][1]-1;

			*ps=Uconsloop[FIS];
			*pe=Uconsloop[FIE];

			if(dir==2){
				// This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
				*js=Uconsloop[FJS];
				*je=Uconsloop[FJE]+1;
			}
			else{
				// GODMARK: needed?
				*js=Uconsloop[FJS];
				*je=Uconsloop[FJE];
				//	*js=INFULL2;
				//	*je=OUTFULL2;
			}

			if(dir==3){
				// This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
				*ks=Uconsloop[FKS];
				*ke=Uconsloop[FKE]+1;
			}
			else{
				// GODMARK: needed?
				*ks=Uconsloop[FKS];
				*ke=Uconsloop[FKE];
				//*ks=INFULL3;
				//	*ke=OUTFULL3;
			}

		}


		*di=1;
		*dj=1;
		*dk=1;
	}
	else if( ( (!dir_exception[2]) && (dir==2) ) || (interporflux==ENOFLUXAVG2TYPE) ){
		*js=*je=0; // anything so dj=1 iterates, next subloop overwrites j actually used

		if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE)){
			// then centered quantity, and need from -1 .. N and have ijkminmax[1][0/1] outer range of ghost zones
			*bs=ijkminmax[interporflux][2][0];
			*be=ijkminmax[interporflux][2][1];

			*ps=Uconsloop[FJS]-1;
			*pe=Uconsloop[FJE]+1;

			*is=INFULL1;
			*ie=OUTFULL1;

			*ks=INFULL3;
			*ke=OUTFULL3;
		}
		else if(interporflux==ENOFLUXRECONTYPE){
			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][2][0]+1;
			*be=ijkminmax[interporflux][2][1];

			*ps=Uconsloop[FJS];
			*pe=Uconsloop[FJE]+1;

			*is=INFULL1;
			*ie=OUTFULL1;

			*ks=INFULL3;
			*ke=OUTFULL3;
		}
		else if( interporflux==ENOAVG2CENTTYPE ){
			// input and output at same location

			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=Uconsloop[FJS];  //atch correct  SASMARK; ALSO need to set up the same thing for other dimensions-- not sure if this is enough
			*be=Uconsloop[FJE];

			*ps=0;
			*pe=N2 - 1;

			*is=INFULL1;
			*ie=OUTFULL1;

			*ks=INFULL3;
			*ke=OUTFULL3;
		}
		else if(interporflux==ENOCENT2AVGTYPE){
			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][2][0]+1;
			*be=ijkminmax[interporflux][2][1]-1;

			*ps=Uconsloop[FJS];
			*pe=Uconsloop[FJE];

			*is=INFULL1;
			*ie=OUTFULL1;

			*ks=INFULL3;
			*ke=OUTFULL3;
		}
		else if(interporflux==ENOFLUXAVG2TYPE){
			//is bs and be to be corrected for the shock indicator? SASMARK
			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][2][0]+1;
			*be=ijkminmax[interporflux][2][1]-1;

			*ps=Uconsloop[FJS];
			*pe=Uconsloop[FJE];

			if(dir==1){
				// This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
				*is=Uconsloop[FIS];
				*ie=Uconsloop[FIE]+1;
			}
			else{
				*is=Uconsloop[FIS];
				*ie=Uconsloop[FIE];
				// GODMARK: needed?
				//	*is=INFULL1;
				//	*ie=OUTFULL1;
			}

			if(dir==3){
				// This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
				*ks=Uconsloop[FKS];
				*ke=Uconsloop[FKE]+1;
			}
			else{
				// GODMARK: needed?
				*ks=Uconsloop[FKS];
				*ke=Uconsloop[FKE];
				//	*ks=INFULL3;
				//	*ke=OUTFULL3;
			}

		}


		*di=1;
		*dj=1;
		*dk=1;
	}
	else if( ( (!dir_exception[3]) && (dir==3) ) || (interporflux==ENOFLUXAVG3TYPE) ){
		*ks=*ke=0; // anything so dk=1 iterates, next subloop overwrites k actually used


		if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE)){
			// then centered quantity, and need from -1 .. N and have ijkminmax[1][0/1] outer range of ghost zones
			*bs=ijkminmax[interporflux][3][0];
			*be=ijkminmax[interporflux][3][1];

			*ps=Uconsloop[FKS]-1;
			*pe=Uconsloop[FKE]+1;

			*is=INFULL1;
			*ie=OUTFULL1;

			*js=INFULL2;
			*je=OUTFULL2;

		}
		else if(interporflux==ENOFLUXRECONTYPE){
			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][3][0]+1;
			*be=ijkminmax[interporflux][3][1];

			*ps=Uconsloop[FKS];
			*pe=Uconsloop[FKE]+1;

			*is=INFULL1;
			*ie=OUTFULL1;

			*js=INFULL2;
			*je=OUTFULL2;

		}
		else if( interporflux==ENOAVG2CENTTYPE ){
			// input and output at same location

			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=Uconsloop[FKS];  //atch correct  SASMARK; ALSO need to set up the same thing for other dimensions-- not sure if this is enough
			*be=Uconsloop[FKE];

			*ps=0;
			*pe=N3 - 1;

			*is=INFULL1;
			*ie=OUTFULL1;

			*js=INFULL2;
			*je=OUTFULL2;
		}
		else if(interporflux==ENOCENT2AVGTYPE){
			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][3][0]+1;
			*be=ijkminmax[interporflux][3][1]-1;

			*ps=Uconsloop[FKS];
			*pe=Uconsloop[FKE];

			*is=INFULL1;
			*ie=OUTFULL1;

			*js=INFULL2;
			*je=OUTFULL2;

		}
		else if(interporflux==ENOFLUXAVG3TYPE){

			// then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
			*bs=ijkminmax[interporflux][3][0]+1;
			*be=ijkminmax[interporflux][3][1]-1;

			*ps=Uconsloop[FKS];
			*pe=Uconsloop[FKE];

			if(dir==1){
				// This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
				*is=Uconsloop[FIS];
				*ie=Uconsloop[FIE]+1;
			}
			else{
				// GODMARK: needed?
				*is=Uconsloop[FIS];
				*ie=Uconsloop[FIE];
				//	*is=INFULL1;
				//	*ie=OUTFULL1;
			}

			if(dir==2){
				// This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
				*js=Uconsloop[FJS];
				*je=Uconsloop[FJE]+1;
			}
			else{
				*js=Uconsloop[FJS];
				*je=Uconsloop[FJE];
				// GODMARK: needed?
				//	*js=INFULL2;
				//	*je=OUTFULL2;
			}

		}

		*di=1;
		*dj=1;
		*dk=1;

	}

}



// Get 1D line of data so can pass it to ENO scheme (used for c2e only)
void get_1d_line_c2e(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE *yin)
{
	int yiniter;
	int di2,dj2,dk2;
	int i2,j2,k2;
	int is2,ie2,js2,je2,ks2,ke2;
	int dir_exception[NDIM];


	dirglobal=dir;
	iterglobal=dir;
	interporfluxglobal=interporflux;


	// determine range of outer loop and range to feed to eno scheme
	if(dir==1){

		iglobal=0;
		is2=bs;
		ie2=be;

		jglobal=js2=je2=j;

		kglobal=ks2=ke2=k;

		di2=1;
		dj2=1;
		dk2=1;

	}
	else if(dir==2){

		iglobal=is2=ie2=i;

		jglobal=0;
		js2=bs;
		je2=be;

		kglobal=ks2=ke2=k;

		di2=1;
		dj2=1;
		dk2=1;
	}
	else if(dir==3){

		iglobal=is2=ie2=i;

		jglobal=js2=je2=j;

		kglobal=0;
		ks2=bs;
		ke2=be;

		di2=1;
		dj2=1;
		dk2=1;
	}

	// reset to 0 so eno schemes don't care about values there (assume weights set to also 0 there)
	for(yiniter=-NBIGBND;yiniter<NBIG+NBIGBND;yiniter++){
		yin[yiniter] = 0;
	}


	// 1 input, assumed interporflux==ENOINTERPTYPE
	yiniter=bs;
	SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){

		yin[yiniter]=p2interp[i2][j2][k2][pl];
		yiniter++;
	}

	if(interporflux!=ENOINTERPTYPE){
		dualfprintf(fail_file,"get_1d_line_c2e only handles interporflux==ENOINTERPTYPE\n");
		myexit(25);
	}



}




// Get 1D line of data so can pass it to ENO scheme
void get_1d_line(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interpm)[N2M][N3M][NPR],FTYPE (*p2interpp)[N2M][N3M][NPR], FTYPE (*yin)[NBIGM])
{
	int yiniter;
	int di2,dj2,dk2;
	int i2,j2,k2;
	int is2,ie2,js2,je2,ks2,ke2;
	int dir_exception[NDIM];


	if(
		((interporflux==ENOFLUXAVG1TYPE)&&(dir==1))||
		((interporflux==ENOFLUXAVG2TYPE)&&(dir==2))||
		((interporflux==ENOFLUXAVG3TYPE)&&(dir==3))
		){
			dualfprintf(fail_file,"No such method with interporflux=%d and dir=%d\n",interporflux,dir);
			myexit(1);
	}

	dir_exception[1] =  (interporflux==ENOFLUXAVG2TYPE) || (interporflux==ENOFLUXAVG3TYPE);
	dir_exception[2] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG3TYPE);
	dir_exception[3] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG2TYPE);

	dirglobal=dir;
	interporfluxglobal=interporflux;

	// determine range of outer loop and range to feed to eno scheme
	if( ( (!dir_exception[1]) && (dir==1) ) || (interporflux==ENOFLUXAVG1TYPE) ){

		iterglobal=1;

		iglobal=0;
		is2=bs;
		ie2=be;

		jglobal=js2=je2=j;

		kglobal=ks2=ke2=k;

		di2=1;
		dj2=1;
		dk2=1;

	}
	else if( ( (!dir_exception[2]) && (dir==2) ) || (interporflux==ENOFLUXAVG2TYPE) ){

		iterglobal=2;

		iglobal=is2=ie2=i;

		jglobal=0;
		js2=bs;
		je2=be;

		kglobal=ks2=ke2=k;

		di2=1;
		dj2=1;
		dk2=1;
	}
	else if( ( (!dir_exception[3]) && (dir==3) ) || (interporflux==ENOFLUXAVG3TYPE) ){

		iterglobal=3;

		iglobal=is2=ie2=i;

		jglobal=js2=je2=j;

		kglobal=0;
		ks2=bs;
		ke2=be;

		di2=1;
		dj2=1;
		dk2=1;
	}

	// reset to 0 so eno schemes don't care about values there (assume weights set to also 0 there)
	for(yiniter=-NBIGBND;yiniter<NBIG+NBIGBND;yiniter++){
		yin[0][yiniter] = yin[1][yiniter] = 0;
	}


	if( (interporflux==ENOINTERPTYPE)||(interporflux==ENOFLUXRECONTYPE)||(interporflux==ENOFLUXAVG1TYPE)||(interporflux==ENOFLUXAVG2TYPE)||(interporflux==ENOFLUXAVG3TYPE)||(interporflux==ENOAVG2CENTTYPE)||(interporflux==ENOCENT2AVGTYPE) ){// these have only 1 input
		yiniter=bs;
		SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){

			yin[0][yiniter]=p2interpm[i2][j2][k2][pl];
			yiniter++;
		}
	}
	else if(interporflux==ENOFLUXSPLITTYPE){ // this method has 2 inputs
		yiniter=bs;
		SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){

			yin[0][yiniter]=p2interpm[i2][j2][k2][pl];
			if(p2interpp!=NULL)  yin[1][yiniter]=p2interpp[i2][j2][k2][pl];
			yiniter++;
		}
	}





}





// Figure out shifting of stencil and order of stencil to more closely match with causality
// (used for both c2e and other routines)
void causal_shift_order(int dir, int interporflux, int preforder, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, int *shift, int *minorder, int *maxorder)
{
	int i3,j3,k3;
	int is3,ie3,js3,je3,ks3,ke3;
	int di3, dj3, dk3;
	int temporder;
	int superdiv;
	FTYPE	wspeed0l,wspeed0r,wspeed1l,wspeed1r;
	FTYPE wspeed0ll,wspeed1ll;
	int yiniter;
	FTYPE localspeed[2];
	int shifttemp;

	int dir_exception[NDIM];

	if(
		((interporflux==ENOFLUXAVG1TYPE)&&(dir==1))||
		((interporflux==ENOFLUXAVG2TYPE)&&(dir==2))||
		((interporflux==ENOFLUXAVG3TYPE)&&(dir==3))
		){
			dualfprintf(fail_file,"No such method with interporflux=%d and dir=%d\n",interporflux,dir);
			myexit(1);
	}

	dir_exception[1] =  (interporflux==ENOFLUXAVG2TYPE) || (interporflux==ENOFLUXAVG3TYPE);
	dir_exception[2] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG3TYPE);
	dir_exception[3] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG2TYPE);

	if( ( (!dir_exception[1]) && (dir==1) ) || (interporflux==ENOFLUXAVG1TYPE) ){
		is3=ps;
		ie3=pe;

		js3=je3=j;

		ks3=ke3=k;

		di3=1;
		dj3=1;
		dk3=1;
	}
	else if( ( (!dir_exception[2]) && (dir==2) ) || (interporflux==ENOFLUXAVG2TYPE) ){
		is3=ie3=i;

		js3=ps;
		je3=pe;

		ks3=ke3=k;

		di3=1;
		dj3=1;
		dk3=1;
	}
	else if( ( (!dir_exception[3]) && (dir==3) ) || (interporflux==ENOFLUXAVG3TYPE) ){
		is3=ie3=i;

		js3=je3=j;

		ks3=ps;
		ke3=pe;

		di3=1;
		dj3=1;
		dk3=1;
	}


	//#if( (STOREWAVESPEEDS)&& ((VCHARTYPE==GLOBALVCHAR)||(VCHARTYPE==LOCALVCHAR)) ) // this procedure makes no sense with GLOBALVCHAR
#if( (STOREWAVESPEEDS)&& ((VCHARTYPE==LOCALVCHAR)||(VCHARTYPE==VERYLOCALVCHAR) ) )
	// GODMARK: This will use precomputed wave speeds in wspeed[dir][2][i][j][k]
	// wspeed located at cell interface.  Take this into account when forming shifter.
	// therefore wspeed is wave speeds from interface point of view.  For centered quantities should consider average (or max/min) of wspeed.
	//
	// e.g. for  c2e, average wspeed[i] and wspeed[i+1]
	//      for  a2c (as used for ENOFLUXRECONTYPE, not really center, but edge!) then wspeed[i] is correct one
	//      for  a2em/p average wspeed[i] and wspeed[i+1]
	//
	// shift -order/2-1 (e.g. -2 for WENO5) if flow is superRIGHT (when both left/right chars are +)
	// shift order/2-1 (e.g. 2 for WENO5) if flow is superLEFT (when both left/right chars are -)
	//
	// smoothly vary between and feed integer value (which is interpreted correctly for avg2cent vs. avg2edge vs. cent2edge
	yiniter=ps; // note starts at ps not bs since shift only used on points of interest
	SUPERGENLOOP(i3,j3,k3,is3,ie3,js3,je3,ks3,ke3,di3,dj3,dk3){ // only over points of interest

		// first determine if should reduce order because of supersonic divergence
		// first assume maximum preferred order
		temporder=preforder;

		// get standard wavespeed (used by any method)
		wspeed0l=wspeed[dir][0][i3][j3][k3];
		wspeed1l=wspeed[dir][1][i3][j3][k3];			


		if((interporflux==ENOINTERPTYPE)||(interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOFLUXAVG1TYPE)||(interporflux==ENOFLUXAVG2TYPE)||(interporflux==ENOFLUXAVG3TYPE)){ // quantities are at CENT-dir (assumes ENOFLUXAVG?TYPE is orthogonal to dir)


			// get grid-on-right wavespeed
			wspeed0r=wspeed[dir][0][i3+idel][j3+jdel][k3+kdel];
			wspeed1r=wspeed[dir][1][i3+idel][j3+jdel][k3+kdel];


			// check for superfast divergence
			superdiv=0;
			if(SUPERFASTDIVREDUCE&&(wspeed1l<0)&&(wspeed0r>0)){
				superdiv=1;
				temporder=1;
			}

			// get correctly positioned wave speed
			localspeed[0]=min(wspeed0l,wspeed0r);
			localspeed[1]=max(wspeed1l,wspeed1r);

		}
		else if((interporflux==ENOFLUXRECONTYPE)||(interporflux==ENOAVG2CENTTYPE)||(interporflux==ENOCENT2AVGTYPE)){ // quantities are at FACE-dir

			// get grid-on-right wavespeed
			wspeed0r=wspeed[dir][0][i3+idel][j3+jdel][k3+kdel];
			wspeed1r=wspeed[dir][1][i3+idel][j3+jdel][k3+kdel];

			// get grid-on-left wavespeed
			wspeed0ll=wspeed[dir][0][i3-idel][j3-jdel][k3-kdel];
			wspeed1ll=wspeed[dir][1][i3-idel][j3-jdel][k3-kdel];


			// check for superfast divergence
			superdiv=0;
			if(SUPERFASTDIVREDUCE&&(wspeed1ll<0)&&(wspeed0r>0)){
				superdiv=1;
				temporder=1;
			}


			// get correctly positioned wave speed as average of surroundings
			// wave speeds at interface
			localspeed[0]=wspeed0l;
			localspeed[1]=wspeed1l;
		}

		// only shift if not diverging superfast
		if(superdiv==0){
			// linear interpolation between left and right super"sonic" directed shifts back downstream
			shifttemp=-((temporder+1)/2)-(temporder+1)*(FTYPE)(localspeed[0]/(localspeed[1]-localspeed[0]+SMALL));
			if(shifttemp>(temporder+1)/2) shifttemp=(temporder+1)/2;
			if(shifttemp<-((temporder+1)/2)) shifttemp=-((temporder+1)/2);
		}
		else shifttemp=0; // no shift if super divergence

		// now assign results to arrays of minorder,maxorder, and shifts.
		minorder[yiniter]=MIN(MINPREFORDER,temporder);
		maxorder[yiniter]=temporder;
		shift[yiniter]=shifttemp;
		yiniter++;
	}
#else
	// NO SHIFT OR CHANGE OF ORDER
	yiniter=ps; // GODMARK: was bs
	SUPERGENLOOP(i3,j3,k3,is3,ie3,js3,je3,ks3,ke3,di3,dj3,dk3){ // only over points of interest  //atch correct -- was di3, dj3, dk3
		minorder[yiniter]=MIN(preforder,MINPREFORDER); // minimum preferred order
		maxorder[yiniter]=preforder;
		shift[yiniter]=0;// then no shift since don't have wave speeds (not stored)
		yiniter++;
	}
#endif

}




// determine shock indicator (used for both c2e and other routines, only operates on quantities with NPR elements)
int get_shock_indicator(int dir, int interporflux, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yin)[2][NBIGM], FTYPE *shockindicator)
{
	//extern int Utoprimgen(int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr);
	//extern void get_geometry(int ii, int jj, int kk, int pp, struct of_geom *geom);

	int i3,j3,k3;
	int is3,ie3,js3,je3,ks3,ke3;
	int di3, dj3, dk3;
	FTYPE interplistpl[NPR][MAXSPACEORDER];
	FTYPE *ypl[NPR];
	int plpl;
	int startorderi,endorderi;
	int yiniter;
	int l;
	extern FTYPE  Ficalc(int dir, FTYPE **ypl);





	if(((interporflux!=ENOFLUXAVG1TYPE)&&(dir==1))||(interporflux==ENOFLUXAVG1TYPE)){
		is3=ps;
		ie3=pe;

		js3=je3=j;

		ks3=ke3=k;

		di3=1;
		dj3=1;
		dk3=1;
	}
	else if(((interporflux!=ENOFLUXAVG2TYPE)&&(dir==2))||(interporflux==ENOFLUXAVG2TYPE)){
		is3=ie3=i;

		js3=ps;
		je3=pe;

		ks3=ke3=k;

		di3=1;
		dj3=1;
		dk3=1;
	}
	else if(((interporflux!=ENOFLUXAVG3TYPE)&&(dir==3))||(interporflux==ENOFLUXAVG3TYPE)){
		is3=ie3=i;

		js3=je3=j;

		ks3=ps;
		ke3=pe;

		di3=1;
		dj3=1;
		dk3=1;
	}

	yiniter=ps; // note starts at ps not bs since shift only used on points of interest

	startorderi = - (7)/2; // order=7 fixed for shock detector (must make sure MAXBND==5)
	endorderi   = - startorderi;

	// shift pointer
	PLOOP(plpl) ypl[plpl] = interplistpl[plpl] - startorderi;

	SUPERGENLOOP(i3,j3,k3,is3,ie3,js3,je3,ks3,ke3,di3,dj3,dk3){ // only over points of interest

		PLOOP(plpl){
			// get interpolation points, where y[0] is point of interest for which interpolation is found.
			for(l=startorderi;l<=endorderi;l++){
				ypl[plpl][l]=yin[plpl][0][i3*idel+j3*jdel+k3*kdel + l];
			}
		}

		shockindicator[yiniter]=Ficalc(dir,ypl);

		yiniter++;
	}

	return( 0 );

}




// real calculation.  Gets derivatives of input function and passes this to various other function
int compute_df_line(int interporflux, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE *yin, FTYPE (*df)[NBIGM])
{
	int i;


	//  return(0);

	// df
	for(i=bs+1;i<=be;i++) df[0][i] = yin[i]-yin[i-1];

	// d^2f
	for(i=bs+1;i<=be-1;i++) df[1][i] = df[0][i+1]-df[0][i];

	if(CONTACTINDICATOR){
		// only used for contact indicator, so don't compute if not needed
		// d^3f
		for(i=bs+2;i<=be-1;i++) df[2][i] = df[1][i]-df[1][i-1];
	}


#if(NUMDFS>=4)
	// df centered, for Sasha's smoothness indicators
	for(i=bs+1;i<=be-1;i++) df[3][i] = yin[i+1]-yin[i-1];
#endif

#if(NUMDFS>=5)
	// df centered 2 apart, for Sasha's smoothness indicators
	for(i=bs+2;i<=be-2;i++) df[4][i] = yin[i+2]-yin[i-2];
#endif


	return(0);
}


// returns the improved pressure jump indicators
int compute_df_line_new(int numprims, int interporflux, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*yin)[2][NBIGM], FTYPE *df)
{
	int i;
	struct of_geom geom,geoml,geomr;
	int num;
	int pl;
	int di,dj,dk;
	FTYPE myprim[NPR];
	FTYPE mypriml[NPR];
	FTYPE myprimr[NPR];
	//  FTYPE U[NPR];
	struct of_state q,ql,qr;
	void compute_1plusud0(FTYPE *pr,struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0); // plus1ud0=(1+q->ucov[TT])
	FTYPE plus1ud0,plus1ud0l,plus1ud0r;
	FTYPE par,parl,parr;
	FTYPE uparl,uparr,upar;



	di=(iterglobal==1);
	dj=(iterglobal==2);
	dk=(iterglobal==3);

	for(num=bs+1;num<=be-1;num++){  //SASMARK added "-1"
		// GODMARK: metric constant
		get_geometry(iglobal+di*num,jglobal+dj*num,kglobal+dk*num,CENT,&geom);
		get_geometry(iglobal+di*(num-1),jglobal+dj*(num-1),kglobal+dk*(num-1),CENT,&geoml);
		get_geometry(iglobal+di*(num+1),jglobal+dj*(num+1),kglobal+dk*(num+1),CENT,&geomr);

		PLOOP(pl){
			mypriml[pl] = yin[pl][0][num-1];
			myprim[pl] = yin[pl][0][num];
			myprimr[pl] = yin[pl][0][num+1];
		}

		get_state(mypriml,&geoml,&ql);
		get_state(myprim,&geom,&q);
		get_state(myprimr,&geomr,&qr);

		compute_1plusud0(mypriml,&geoml,&ql,&plus1ud0l); // plus1ud0=(1+q->ucov[TT])
		compute_1plusud0(myprim,&geom,&q,&plus1ud0); // plus1ud0=(1+q->ucov[TT])
		compute_1plusud0(myprimr,&geomr,&qr,&plus1ud0r); // plus1ud0=(1+q->ucov[TT])

		//The expression is:
		// -(rho v^2/2 + u) = rho * ucon * (1+ucov)  +  gam * u * (ucon*ucov) + (gam-1)*u

		parl=plus1ud0l*(ql.ucon[TT]);
		par=plus1ud0*(q.ucon[TT]);
		parr=plus1ud0r*(qr.ucon[TT]);

		// GODMARK: non-rel only right
		uparl=gam*mypriml[UU]*(ql.ucon[TT]*ql.ucov[TT])  + (gam-1.0)*mypriml[UU];
		upar =gam*myprim[UU] *(q.ucon[TT]  *q.ucov[TT])  +  (gam-1.0)*myprim[UU];
		uparr=gam*myprimr[UU]*(qr.ucon[TT]*qr.ucov[TT])  + (gam-1.0)*myprimr[UU];

		//normalized difference using the central value -- think may lead to problems
		//df[num] = (fabs(myprim[RHO]*0.5*(parr-parl)) + fabs(0.5*(uparr-uparl)))/(myprim[RHO]*par + upar);

		//normalized difference without using the central value; normalize the difference by the values themselves
		//df[num] = ( fabs(myprim[RHO]*(parr-parl)) + fabs(uparr-uparl) ) / ( fabs(myprim[RHO])*(fabs(parr)+fabs(parl)) + fabs(uparr)+fabs(uparl) + DBL_MIN );

		//divide by one point value, not by the sum of the values themselves
		df[num] = ( fabs(myprim[RHO]*(parr-parl)) + fabs(uparr-uparl) ) / ( fabs(myprim[RHO]*par) + fabs(upar) + DBL_MIN );

		// Original dP/P prescription
		//df[num] = ( fabs(uparr-uparl) ) / ( fabs(upar) + DBL_MIN );



	}

	return(0);
}





//    PLOOPINTERP(pl) pass_1d_line(interporflux, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], yin[pl], yout[pl]);

// Pass 1D line to ENO scheme
void pass_1d_line(int whichquantity, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM])
{
	// These are ENO functions
	//extern int eno_line_c2e(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE *shockindicator, FTYPE *yin, FTYPE *yout_left, FTYPE *yout_right);
	//extern int eno_line_a2c(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE *shockindicator, FTYPE *yin, FTYPE *yout);
	//extern int eno_line_c2a(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE *shockindicator, FTYPE *yin, FTYPE *yout);
	//extern int eno_line_a2em(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE *shockindicator, FTYPE *yin, FTYPE *yout);
	//extern int eno_line_a2ep(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE *shockindicator, FTYPE *yin, FTYPE *yout);

	// pass yin so that yin[0] is first active point

	// DEATHMARK : Sasha, these are your functions
	if(recontype==CVT_C2E){
		// yout filled correctly as 1 quantity interpolated to edges
		//      eno_line_c2e(whichreduce,maxorder,minorder,bs,ps,pe,be,shift,yin[0]+bs,yout);
		// returns left primitive from 0 .. N and right primitive from -1 .. N-1
		// i.e. pleft(0..N) and pright(-1..N-1) as defined below
		// so inclusive on pleft/right need -1..N
		//////////////////////////////////////
		//
		// interpolate primitive using slope
		//
		// |=interface
		// i=zone center of ith zone
		//
		// |         pl(i)|pr(i)    i          |
		// |         Fl(i)|Fr(i)    i          |
		// |         Ul(i)|Ur(i)    i          |
		// |              |pleft(i)   pright(i)|
		// |              |F(i)                |
		//
		//
		//
		//////////////////////////////////////
		//atch: do not shift yin[] -- otherwise need to shift all other indices; spit out two arrays for left and right
		//eno_line_c2e(whichreduce,maxorder,minorder,bs,ps,pe,be,shift,yin[0], yout[0], yout[1]);  

		//GODMARK: pass maxorder and minorder as arrays for each point of interest
		eno_line_c2e( whichquantity, do_weight_or_recon, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift,  shockindicator, df, dP, monoindicator, yprim[UU][0], yin[0], yout[0], yout[1] );  
	}
	else if(recontype==CVT_A2C){
		// yout filled correctly as 1 quantity "deconstructed" to 1 quantity
		//	 eno_line_a2c(whichreduce,maxorder,minorder,bs,ps,pe,be,shift,yin[0],yout[0]);
		// for ENOFLUXRECONTPE: assume previously bounded fluxes, then returns fluxes from 0 .. N (not just N-1) ( i.e. F(0..N) as defined below)
		// for ENOAVG2CENTTYPE: assume previously bounded average U's from 0-order/2 .. N-1+order/2, then return U's from 0 .. N-1
		//////////////////////////////////////
		//
		// interpolate primitive using slope
		//
		// |=interface
		// i=zone center of ith zone
		//
		// |         pl(i)|pr(i)    i          |
		// |         Fl(i)|Fr(i)    i          |
		// |         Ul(i)|Ur(i)    i          |
		// |              |F(i)                |
		//
		//
		//
		//////////////////////////////////////

#if(0) // DEBUG TEST
		// no i2,j2, etc. present anymore
		//      yiniter=bs;
		//      SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){
		//
		//	yout[0][yiniter]=yin[0][yiniter];
		//	yiniter++;
		//      }
#endif
		eno_line_a2c( whichquantity, do_weight_or_recon, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df, dP, monoindicator, yprim[UU][0], yin[0], yout[0] ); 
	}
	// if(interporflux==ENOFLUXSPLITTYPE){
	// yout filled as each left/right prepared quantity is interpolated as an average to the left/right edges
	//	eno_line_a2em(1,preforder,bs,ps,pe,be,minorder,maxorder,shift,yin[0],yout[0]);
	//	eno_line_a2ep(1,preforder,bs,ps,pe,be,minorder,maxorder,shift,yin[1],yout[1]);
	// GODMARK: should return left fluxes from 0 .. N and right fluxes from -1 .. N-1
	// i.e. Fleft(0..N) and Fright(-1..N-1) as defined below
	// so inclusive on Fleft/right need -1..N
	//////////////////////////////////////
	//
	// convert point fluxes to quasi-de-averaged point fluxes at cell interface
	//
	// |=interface
	// i=zone center of ith zone
	//
	// |         pl(i)|pr(i)    i          |
	// |         Fl(i)|Fr(i)    i          |
	// |         Ul(i)|Ur(i)    i          |
	// |              |        F+-(i)      |
	// |              |Fleft(i)   Fright(i)|
	// |              |F(i)                |
	//
	//
	//
	//////////////////////////////////////
	//  }
	else if(recontype==CVT_C2A){
		// yout filled as each left/right prepared quantity is interpolated as an average to the left/right edges
		//	eno_line_c2a(1,preforder,bs,ps,pe,be,minorder,maxorder,shift,yin[0],yout[0]);
		// For F1:
		// INPUT  : Inputs fluxes from j=-N2BND .. N2-1+N2BND inclusive for F1(i,j) for a single i
		// OUTPUT : Returns fluxes from j=0 .. N2-1 for F1(i,j) for a single i (i.e. F1(i,j=0..N2-1)
		// where for F1 interporflux==ENOFLUXAVG1TYPE doesn't apply, ENOFLUXAVG2TYPE involves averaging over second direction and ENOFLUXAVG3TYPE involves averaging over third direction
		// the input data is controlled by the loops above
		//////////////////////////////////////
		//
		// convert point flux to surface integrated flux (this call does only 1 of the NDIMEN-1 directions)
		//
		// |=interface
		// i=zone center of ith zone
		// |              |                    |
		// |            F(i)        i          |
		// |              |                    |
		// |              |                    |
		// |              |                    |
		//
		//
		//
		//////////////////////////////////////
		eno_line_c2a( whichquantity, do_weight_or_recon, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df, dP, monoindicator, yprim[UU][0], yin[0], yout[0] );   //atch correct
	}


}







// Assign result of ENO operation to final array
// for c2e
void assign_eno_result_c2e(int recontype, int pl, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[NBIGM], FTYPE (*result0)[N2M][N3M][NPR2INTERP], FTYPE (*result1)[N2M][N3M][NPR2INTERP])
{
	int l;


	for(l=ps;l<=pe;l++){ // inclusive loop
		result0[i+l*idel][j+l*jdel][k+l*kdel][pl] = yout[0][l];
		result1[i+l*idel][j+l*jdel][k+l*kdel][pl] = yout[1][l];
	}

	if(recontype!=CVT_C2E){
		dualfprintf(fail_file,"assign_eno_result_c2e only handles recontype==CVT_C2E\n");
		myexit(26);
	}

}


// Assign result of ENO operation to final array
// only care about how many outputs
void assign_eno_result(int recontype, int pl, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[NBIGM], FTYPE (*result0)[N2M][N3M][NPR], FTYPE (*result1)[N2M][N3M][NPR])
{
	int l;


	if(recontype==CVT_C2E){// these methods generated 2 output values
		for(l=ps;l<=pe;l++){ // inclusive loop
			result0[i+l*idel][j+l*jdel][k+l*kdel][pl] = yout[0][l];
			result1[i+l*idel][j+l*jdel][k+l*kdel][pl] = yout[1][l];
		}
	}
	else if( (recontype==CVT_C2A)||(recontype==CVT_A2C) ){// these methods generated 1 output value
		for(l=ps;l<=pe;l++){ // inclusive loop
			// result0 can be equal to input p2interpm array since now done with that line and dimensionally split, so never operate in another direction on this quantity
			result0[i+l*idel][j+l*jdel][k+l*kdel][pl] = yout[0][l];
		}
	}

}







