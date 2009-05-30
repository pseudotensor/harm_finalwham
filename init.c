//LIST OF CHANGES TO THIS VERSION OF CODE --- atch
//This version of the code has following corrected:
// global.h -- UTORPIMFAILRETURNTYPE = UTOPRIMRETURNADJUSTED (so that u < 0 can be handled by the inversion and can be stepped over)
// set_dt() is called in initbase.c; all test-specific setting of dt is now done in post_init_specific_init()
// declared static functions to be static in reconstructeno.c
// added init.h to the list of dependencies of init.o

/* 
 *
 * generates initial conditions for a set of 1-D non-relativistic hydro test problems
 * those with positive indices are shock problems; those with indices <= 0 are
 * variations of "Hubble" expanding flows
 *
 */

#include "decs.h"
#include "init.h" //atch
#include "coord.h"
#include "reconstructeno.h"

static int pi2Uavg_specific(FTYPE (*prim)[N2M][N3M][NPR], FTYPE (*Upoint)[N2M][N3M][NPR], FTYPE (*Uavg)[N2M][N3M][NPR]);
static void write_riemannproblem_params_to_file( FTYPE rhol, FTYPE pl, FTYPE ul, FTYPE uly, FTYPE rhor, FTYPE pr, FTYPE ur, FTYPE ury, FTYPE *Bl, FTYPE *Br, FTYPE tf, FTYPE timescalefactor, FTYPE x0, FTYPE a, FTYPE b );

#define REL1DRIEMANNTEST( no ) ( no >= 101 && no <= 149  || no >=1000 && no <=1100)  //whether a test is a relativistic 1d Riemann problem
#define NONREL1DRIEMANNTEST( no ) ((no >= 1 && no <= 7) || (no >=10 && no <= 13) || (no == 99) || (no == 666) || (no >= 51 && no <= 99)) //whether a test is a non-relativistic 1d Riemann problem

#define MAGTEST( no ) (FFDETEST(no) || ( no >=1000 && no <=1100 )) // whether a test uses magnetic field

#define MAGTESTNOPOT( no ) ( (no >= 200 && no <= 208) || no == 1001 || no == 1002 ) // whether a test sets magnetic field (1) or vector potential (0)


#if( TESTNUMBER == 29 )
FILE *test29fp;
FTYPE test29rho[N1];
FTYPE test29u[N1];
FTYPE test29p[N1];
#endif

#if( TESTNUMBER == 153 )
  //bondi problem analytical solution storage (with boundary zones)
	static FTYPE a_rho_bondi[NBIGM];
	static FTYPE a_u_bondi[NBIGM];
	static FTYPE a_uu1_bondi[NBIGM];
	static FTYPE gam_bondi;
	//shift the array index so that it starts at (-NBIGBND)
	static FTYPE *rho_bondi = &a_rho_bondi[NBIGBND];
	static FTYPE *u_bondi = &a_u_bondi[NBIGBND];
	static FTYPE *uu1_bondi = &a_uu1_bondi[NBIGBND];
#endif

int pre_init_specific_init(void)
{
#if( TESTNUMBER == 29 )
	FTYPE r;
	FTYPE uoverrho;
	int i;
#endif

		//assert( MCOORD != CYLMINKMETRIC && MCOORD != CARTMINKMETRIC && MCOORD !=UNIGRAVITY, "pre_init_specific_init: MCOORD != CARTMINKMETRIC: 1-D tests will not be running correctly" );
	
		//default setting for values that can be used outside of init.c (in initibase.c), e.g. grid parameters.
		//Decide if the test is non-relativistic or not and set timescalefactor appropriately
	  if(NONREL1DRIEMANNTEST(TESTNUMBER) || TESTNUMBER < 100){
	  //		if( TESTNUMBER < 100 || TESTNUMBER == 666 ) {
			//non-relativistic test
			coordparams.timescalefactor = 1e10;
		}
		else {
			//relativistic test: no time stretching
			coordparams.timescalefactor = 1.0;
		}

		trifprintf("Setting timescalefactor=%21.15g\n", coordparams.timescalefactor);

#if( DO_WENO_DEBUG )
		//atch debug SASMARK
		debugfp = fopen( "0_wenodebugfile.out", "wb" ); //overwrite the debug file

		assert( debugfp == NULL, "Could not open the 0_wenodebugfile.out for writing\n" );
#endif //DO_WENO_DEBUG

#if( TESTNUMBER == 29 )
		test29fp = fopen( "../e1rpex.out", "rb" ); //open

		assert( test29fp == NULL, "Cannot open e1rpex.out for TESTNUMBER = %d", (int) TESTNUMBER );

		for( i = 0; i < N1; i++ ) {
			fscanf( test29fp, "%lf %lf %lf %lf %lf\n", &r, &test29rho[i], &test29u[i], &test29p[i], &uoverrho );
		}
#endif
    
    return(0);
}


int init_conservatives(FTYPE p[][N2M][N3M][NPR], FTYPE Utemp[][N2M][N3M][NPR], FTYPE U[][N2M][N3M][NPR])
{
  extern int pi2Uavg(FTYPE (*prim)[N2M][N3M][NPR], FTYPE (*Upoint)[N2M][N3M][NPR], FTYPE (*Uavg)[N2M][N3M][NPR]);

	// use ulast as temporary space
  // uinitial is initial volume averaged conserved quantity
  // needed for advance()

	/////////////
	//
	//  Set up analytic expressions for averaged conserved quantities if available
  // use ulast as temporary space
  // uinitial is initial volume averaged conserved quantity
  // needed for advance()
  //  trifprintf("pi2Uavg\n");
  //	is_pi2Uavg = 1;
  pi2Uavg(p, Utemp, U);
  //	is_pi2Uavg = 0;

  //trifprintf("pi2Uavg_specific\n");
  pi2Uavg_specific(p, Utemp, U);

  return(0);

}


int post_init_specific_init(void)
{

	int testno = TESTNUMBER;

  // globally used parameters set by specific initial condition routines, reran for restart as well *after* all other calculations
	extern int set_dt(FTYPE (*prim)[N2M][N3M][NPR], SFTYPE *dt);  //atch update
	
	set_dt( p, &dt);  //atch update: always use dt determined by the courant factor

	//dt = 0.; //SASMARK
	//SAFE = 1.0;
	
	trifprintf( "TESTNUMBER = %d, dt = %f\n", testno, dt );

	//	if( TESTNUMBER >= 101 && TESTNUMBER <= 199 && REMOVERESTMASSFROMUU == 1 ) {  //Relativistic tests	
	//UTOPRIMVERSION = UTOPRIM2DFINAL;
	//	}
	UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;
	//UTOPRIMVERSION = UTOPRIM5D1;


	trifprintf( "Using REMOVERESTMASSFROMUU = %d, ", REMOVERESTMASSFROMUU );

	if( UTOPRIMVERSION == UTOPRIM2DFINAL ) {
		trifprintf( "UTOPRIM2DFINAL\n" );
	}
	else if( UTOPRIMVERSION == UTOPRIMJONNONRELCOMPAT ) {
		trifprintf( "UTOPRIMJONNONRELCOMPAT\n" );
	}
	else if( UTOPRIMVERSION == UTOPRIM5D1 ) {
		trifprintf( "UTOPRIM5D1\n" );
	}
	//	else {
	//		dualfprintf( fail_file, "cannot use REMOVERESTMASSFROMUU != 1 and != 2\n" );
	//		myexit(1);
	//	}

	//assert( REMOVERESTMASSFROMUU != 1, 
	//	"For relativistic test problems should use the UTOPRIM2D inversion which requires REMOVERESTMASSFROMUU = 1, currently: %d\n", 
	//	(int)REMOVERESTMASSFROMUU );
	//
	////////
			



	////////////
	
  return(0);
}

int pi2Uavg_specific(FTYPE (*prim)[N2M][N3M][NPR], FTYPE (*Upoint)[N2M][N3M][NPR], FTYPE (*Uavg)[N2M][N3M][NPR])
{
	FILE *fp;
	char fname[256];
	int ti, tj;
	FTYPE rho, u, v1, v2, rhoa, Ua, p1a, p2a;
	int nscanned;
	long lineno, numused;
	int i, j, k;
  struct of_geom geom;
	int procno;
	int iswithgdet[NPR];
	int pl;
	FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];


	extern int avg2cen_interp( int whichquantity, int whichavg2cen, FTYPE (*prims_from_avg_cons)[N2M][N3M][NPR], FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR]);

	if( TESTNUMBER == 22 ) { //Implosion
		strcpy( fname, "implosion.dat" );
	}
	else if( TESTNUMBER == 23 ) { //Explosion
		strcpy( fname, "explosion.dat" );
	}
	else if( TESTNUMBER == 26 ) { //stationary Gresho
		strcpy( fname, "gresho.dat" );
	}
	else if( TESTNUMBER == 27 ) { //moving Gresho
		strcpy( fname, "mgresho.dat" );
	}
	else {
		return( 0 );
	}

	if( REMOVERESTMASSFROMUU != 2 ) {
		dualfprintf( fail_file, "REMOVERESTMASSFROMUU has to be 2 for pi2Uavg_specific() to operate correctly, currently it is = %d\n", 
			REMOVERESTMASSFROMUU );
		myexit(1);
	}

	trifprintf( "pi2Uavg_specific(): opening %s (proc = %d)... ", fname, myid );

	///////
	//
	// Find out if the conserved quantities are to be multiplied by gdet = sqrt(-g)
	if( WHICHEOM!=WITHGDET ) {
		// 0 if no gdet is used, 1 if with gdet
		iswithgdet[RHO]=(NOGDETRHO>0) ? 0 : 1;
		iswithgdet[UU]=(NOGDETU0>0) ? 0 : 1;
		iswithgdet[U1]=(NOGDETU1>0) ? 0 : 1;
		iswithgdet[U2]=(NOGDETU2>0) ? 0 : 1;
		iswithgdet[U3]=(NOGDETU3>0) ? 0 : 1;
		iswithgdet[B1]=(NOGDETB1>0) ? 0 : 1;
		iswithgdet[B2]=(NOGDETB2>0) ? 0 : 1;
		iswithgdet[B3]=(NOGDETB3>0) ? 0 : 1;
	}
	else {
		PLOOP(pl) iswithgdet[pl] = 1;  //conserved quantities get multiplied by it later
	}
	//
	///////

#if( USEMPI )
	for( procno = 0; procno < numprocs; procno++ ) {
		if( myid == procno ){
#endif
			//////
			// 
			//  Open file for reading
			fp = fopen( fname, "rb" );
			if( NULL != fp ) {
				trifprintf( "done.\npi2Uavg_specific(): reading data from %s... ", fname );

				lineno = 0;
				numused = 0;
				while( (nscanned = fscanf(fp, "%d %d %lg %lg %lg %lg %lg %lg %lg %lg\r\n", &ti, &tj, &rhoa, &Ua, &p1a, &p2a, &rho, &u, &v1, &v2 )) == 10 ) {
					///////
					//
					//  Check if within range and if within, save it to the array
					//  i, j - local grid indices for this processor, ti, tj - global indices
					i = ti - startpos[1];  
					j = tj - startpos[2]; 
					if( i >= INFULL1 && i <= OUTFULL1 && 
						j >= INFULL2 && j <= OUTFULL2 ) {
							for( k = INFULL3; k <= OUTFULL3; k++ ) {
								Uavg[i][j][k][RHO] = rhoa;
								Uavg[i][j][k][UU] = -Ua / (coordparams.timescalefactor * coordparams.timescalefactor);  //change the sign -- according to HARM convention the total energy is negative
								Uavg[i][j][k][U1] = p1a / coordparams.timescalefactor;
								Uavg[i][j][k][U2] = p2a / coordparams.timescalefactor;
								Uavg[i][j][k][U3] = 0.0;
								Uavg[i][j][k][B1] = 0.0;
								Uavg[i][j][k][B2] = 0.0;
								Uavg[i][j][k][B3] = 0.0;

								get_geometry(i, j, k, CENT, &geom);

								PLOOP(pl) {
									Uavg[i][j][k][pl] *= geom.e[pl];  //multiply conserved quantities by the gdet to have the same representation as inside the code
								}
								

								// center
								coord(i, j, k, CENT, X);
								bl_coord(X, V);
								get_geometry(i,j,k,CENT,&geom) ;
								dxdxprim(X, V, dxdxp);

								//convert the conserved momenta (covariant quantities) from orthonormal basis to coordinate one.
								Uavg[i][j][k][U1] *= dxdxp[1][1];
								Uavg[i][j][k][U2] *= dxdxp[2][2];
							}
							numused++;
					}
					lineno++;
				}

				// Close file
				fclose( fp );
				//
				///////
			}
#if(USEMPI)
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif


	if( NULL == fp ) {
		trifprintf( "file %s does not exist; no changes to conserved quantities done.\n", fname );
		return( -2 );
	}

	if( EOF != nscanned && 0 != nscanned ) {
		dualfprintf( fail_file, "Syntax error in file %s at line %ld.", fname, lineno );
		return(-3);
	}

	trifprintf( "done (%ld lines read, %ld of those used).\n", lineno, numused );

    if( TESTNUMBER == 22 || TESTNUMBER == 23 ) {
	//modify the primitive quantities for nonsmooth functions; for smooth functions the
	//primitive quantities are already good

	trifprintf( "pi2Uavg_specific(): de-averaging the conserved quantities... " );

	///////////
	// 
	//  Average conserved quantities set up; now need to de-average them and invert to get consistent primitive
	//  quantities; use the point values of primitive quantities as a guess
	//avg2cen_interp( ENOCONSERVED, ENOAVG2CENTTYPE, p, Uavg, Upoint );

	trifprintf( "done.\npi2Uavg_specific(): inverting _average_ conserved to primitive quantities... " );

	//  Invert the point conserved quantitites to obtain the primitive ones
	//ZSLOOP( Uconsloop[FIS], Uconsloop[FIE], Uconsloop[FJS], Uconsloop[FJE], Uconsloop[FKS], Uconsloop[FKE] ){
	FULLLOOP {
    get_geometry(i, j, k, CENT, &geom);
		MYFUN(Utoprimgen(EVOLVEUTOPRIM, UEVOLVE, Uavg[i][j][k], &geom, p[i][j][k]),"init.c:pi2Uavg_specific()", "Utoprimgen", 1);
	}

	trifprintf( "done.\npi2Uavg_specific(): bounding primitive quantities... " );

	if (bound_prim(STAGEM1,p) >= 1)
      FAILSTATEMENT("initbase.c:init()", "bound_prim()", 1);

	trifprintf( "done (proc = %d).\n", myid );
     }

	return(1);

}


int readin_bondi_solution(int *pN_bondi, FTYPE *pR0_bondi, FTYPE *pRin_bondi, FTYPE *pRout_bondi, FTYPE *pgam_bondi, FTYPE *prho_bondi, FTYPE *puu1_bondi, FTYPE *pu_bondi )
{
	FILE *fp;
	char fname[256];
	int ti, tj;
	FTYPE rho, u, v1, v2, rhoa, Ua, p1a, p2a;
	int nscanned;
	long lineno, numused;
	int i, j, k;
  struct of_geom geom;
	int procno;
	int iswithgdet[NPR];
	int pl;
	FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
	FTYPE Rcr;
	FTYPE uu1cr;
	FTYPE rhocr;
	FTYPE xstart_bondi;
	FTYPE dx_bondi;
	FTYPE istart_bondi;
	FTYPE uu1;


	sprintf( fname, "bondigrid%d.dat", (int)(N1*ncpux1) );
	//strcpy( fname, "bondigrid.dat" );

	trifprintf( "readin_bondi_solution(): opening %s (proc = %d)... ", fname, myid );

#if( USEMPI )
	for( procno = 0; procno < numprocs; procno++ ) {
		if( myid == procno ){
#endif
			//////
			// 
			//  Open file for reading
			fp = fopen( fname, "rb" );
			if( NULL != fp ) {
				trifprintf( "done.\nreadin_bondi_solution(): reading data from %s... ", fname );

				lineno = 0;
				numused = 0;
					
				nscanned = fscanf(fp, "#Gamma = %lg Rcr = %lg uu1cr = %lg rhocr = %lg\r\n", pgam_bondi, &Rcr, &uu1cr, &rhocr );
				lineno++;

				if( nscanned != 4 ) {
					trifprintf( "Could not read from line %ld of %s\n", lineno, fname );
					return( 1 );
				}

				nscanned = fscanf(fp, "#Nx = %d Rin = %lg Rout = %lg R0 = %lg xstart = %lg dx = %lg istart = %lg\r\n", 
															pN_bondi, pRin_bondi, pRout_bondi, pR0_bondi, 
															&xstart_bondi, &dx_bondi, &istart_bondi);
				lineno++;

				if( nscanned != 7 ) {
					trifprintf( "Could not read from line %ld of %s\n", lineno, fname );
					return( 1 );
				}

				nscanned = fscanf(fp, "#%*[^\n]s\r\n" );  //skip the header desciption line

				lineno++;

				while( (nscanned = fscanf(fp, "%d %*lg %*lg %lg %lg %lg\r\n", &ti, &rho, &uu1, &u )) == 4 ) {
					///////
					//
					//  Check if within range and if within, save it to the array
					//  i, j - local grid indices for this processor, ti, tj - global indices
					i = ti - startpos[1];  
					if( i >= INFULL1 && i <= OUTFULL1 ) {
						prho_bondi[i] = rho;
						puu1_bondi[i] = uu1;
						pu_bondi[i] = u;
						numused++;
					}
					lineno++;
				}

				// Close file
				fclose( fp );
				//
				///////
			}
#if(USEMPI)
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif

	if( NULL == fp ) {
		trifprintf( "file %s does not exist.\n", fname );
		return( -2 );
	}
	
	if( numused < OUTFULL1 - INFULL1 + 1 ) {  //did not fill in all the required values (assuming each value of i appeared only once)
		trifprintf( "Not enough values provided in %s: n_required = %d, n_readin = %d.\n", fname, (int)(OUTFULL1 - INFULL1 + 1), numused );
		return( 2 );
	}

	if( EOF != nscanned && 0 != nscanned ) {
		dualfprintf( fail_file, "Syntax error in file %s at line %ld.", fname, lineno );
		return(-3);
	}

	trifprintf( "done (proc = %d, %ld lines read, %ld of those used).\n", myid, lineno, numused );

	return(0);

}
int init_grid(void)
{
  SFTYPE rh;
	int dim;
	int res;
	int n1;
	FTYPE myangle;

  // coord def type
  defcoord = UNIFORMCOORDS;

  R0 = 0.0; //set for compatibility with the rest of the code that wants some value for a 
  Rhor=0.0; //(for calculations of whether the grid is inside or outside the horizon, etc)
  a = 0.0;  

	//defaults for grid size; so for 1D tests only need to set up the limits only for the test direction unless want to change the defaults
	//so, 1d tests only have to define the x-grid size
	Rin_array[1] = 0.;
	Rout_array[1] = 1.;
	Rin_array[2] = 0.;
	Rout_array[2] = 1.;
	Rin_array[3] = 0.;
	Rout_array[3] = 1.;

  //set_coord_parms();  //does nothing for now -- atch
#if( TESTNUMBER == -3 )
  // Size of the computational domain
  Rin_array[DIRX] = 1.e-20;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == -2 )
  // Size of the computational domain
  Rin_array[DIRX] =  1.e-20;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == -1 )
  // Size of the computational domain
  Rin_array[DIRX] = 1.e-20;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 0 )
  // Size of the computational domain
  Rin_array[DIRX] = -1.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 1 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 2 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 3 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 666 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 2.0;
#elif( TESTNUMBER == 4 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 5 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 6 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 7 ) //peak    
  Rin_array[DIRX] = 0.1;
  Rout_array[DIRX] = 0.6;
#elif( TESTNUMBER == 99 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0 * 1.4;
#elif( TESTNUMBER == 8 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 9 ) //entropy    
  Rin_array[DIRX] = -5.;
  Rout_array[DIRX] = 5.;
#elif( TESTNUMBER == 10 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 11 ) //Lax problem from Qiu and Shu, from "On the Construction, Comparison, and Local
                          //Characteristic Decomposition for High-Order Central WENO Schemes"
  Rin_array[DIRX] = -0.5;
  Rout_array[DIRX] = 0.5;
#elif( TESTNUMBER == 12 ) //Lax problem from Qiu and Shu, from "On the Construction, Comparison, and Local
                          //Characteristic Decomposition for High-Order Central WENO Schemes"
  Rin_array[DIRX] = -0.5;
  Rout_array[DIRX] = 0.5;
#elif( TESTNUMBER == 13 )
  Rin_array[DIRX] = 0.0;
  Rout_array[DIRX] = 1.0;
#elif( TESTNUMBER == 14) 
	Rin_array[DIRX] = 0.;
	Rout_array[DIRX] = 1.;
	Rin_array[DIRY] = 0.;
	Rout_array[DIRY] = 1.;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 15) 
	Rin_array[DIRX] = -1.;
	Rout_array[DIRX] = 1.;
	Rin_array[DIRY] = -1.;
	Rout_array[DIRY] = 1.;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER >= 16 && TESTNUMBER <= 21 )  //test case 3 from L&W: S, S, S, S
	Rin_array[DIRX] = 0.;
	Rout_array[DIRX] = 1.;
	Rin_array[DIRY] = 0.;
	Rout_array[DIRY] = 1.;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 22)  //Implosion, L&W
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 0.3;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 0.3;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;  
#elif( TESTNUMBER == 23 )  //Explosion, L&W
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.5;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.5;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 24 )  //Smooth problem, L&W
	Rin_array[DIRX] = -1.0;
	Rout_array[DIRX] = 1.;
	Rin_array[DIRY] =-1.0;
	Rout_array[DIRY] = 1.;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 25 )  //Odd-even decoupling, L&W
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 0.0125;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 26 )  //Stationary vortex, L&W
	Rin_array[DIRX] = -0.5;
	Rout_array[DIRX] = 0.5;
	Rin_array[DIRY] =-0.5;
	Rout_array[DIRY] = 0.5;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 27 )  //Moving vortex, L&W
	Rin_array[DIRX] = -0.5;
	Rout_array[DIRX] = 3.5;
	Rin_array[DIRY] =-0.5;
	Rout_array[DIRY] = 0.5;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 28 )  // RT instability
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 0.25;
	Rin_array[DIRY] =-0.75;
	Rout_array[DIRY] = 0.75;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 29 )  // rarefaction from peak
  Rin_array[DIRX] = 0.1;
  Rout_array[DIRX] = 0.6;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 30 )  //Smooth sound wave problem (2D, arbitrary angle)
	myangle = 26.565 / 180. * M_PI;
	Rin_array[DIRX] = 0;
	Rout_array[DIRX] = 1. / cos(myangle);
	Rin_array[DIRY] = 0;
	Rout_array[DIRY] = 1. / sin(myangle);
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 31 )  //Smooth sound wave problem (1D)
	Rin_array[DIRX] = 0;
	Rout_array[DIRX] = 1;
	Rin_array[DIRY] = 0;
	Rout_array[DIRY] = 1;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 32 )  //Smooth sound wave problem (2:1 ratio)
	Rin_array[DIRX] = 0;
	Rout_array[DIRX] = 2/sqrt(5);
	Rin_array[DIRY] = 0;
	Rout_array[DIRY] = 1./sqrt(5);
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 33 )  //Smooth density wave problem (2:1 ratio)
	Rin_array[DIRX] = 0;
	Rout_array[DIRX] = 2/sqrt(5);
	Rin_array[DIRY] = 0;
	Rout_array[DIRY] = 1/sqrt(5);
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 49 )
  Rin_array[DIRX] = -0.5;
  Rout_array[DIRX] = 0.5;
#elif( TESTNUMBER == 51 ) // Sod's 1D Riemann problem
	Rin_array[DIRX] = -1.5;
	Rout_array[DIRX] = 1.5;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 52 ) // Sod's 1D Riemann problem boosted
	Rin_array[DIRX] = -1.5;
	Rout_array[DIRX] = 25.5;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 101 ) // 1D Riemann problem #1 from RAM paper
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 102 ) // 1D Riemann problem #2 from RAM paper
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 103 ) // 1D Riemann problem #3 from RAM paper
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 104 ) // 1D Riemann problem #4 from RAM paper
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 105 ) // 1D Riemann problem #5 from RAM paper
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 106 ) // 1D Riemann problem #6.1 from RAM paper (Hard test)
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 107 ) // left-going shock from the 1D Riemann problem #3 from RAM paper
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 151 ) // 1D Riemann problem #5 from RAM paper
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 152 ) // Slab jet test from RAM
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 15.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 45.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 2.*M_PI;  
#elif( TESTNUMBER == 153 ) // Bondi problem
	assert( DIRX != 1 || DIRY != 2 || DIRZ != 3, "Axes should be directed along the default directions for the Bondi problem, TESTNUMBER = %d\n", (int)TESTNUMBER );
	
	//sets the grid parameters: R0, Rin, Rout
	res = readin_bondi_solution( &n1, &R0, &Rin, &Rout, &gam_bondi, rho_bondi, uu1_bondi, u_bondi );

	if( 0 != res ) { //the above call did not succeed
		myexit(1);
	}

	if( N1 * ncpux1 != n1 ) { //compare the resolution of the read-in analytic solution to the current total resolution
		trifprintf( "The total resolution of the simulation, %d in the r-direction, does not match that of the analytic solution, %d\n", (int)(N1 * ncpux1), n1 );
		myexit(1);
	}

	Rin_array[DIRX] = Rin;
	Rout_array[DIRX] = Rout;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = M_PI;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 2. * M_PI;

	a = 0.0;
	gam = gam_bondi;

	//choose coordinates
	defcoord = LOGRSINTH;  //logarithmic in r and uniform in theta
	hslope = 1.0;

#elif( TESTNUMBER == 154 ) //Torus problem
  // metric stuff first
  a = 0.95;
  

  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.0;
  Rhor=rhor_calc(0);
  Rout = 20.;
 
  //Rin=setRin(setihor());
  Rin = 0.98 * Rhor;
  
	Rin_array[DIRX] = Rin;      //these used below in this fcn
	Rout_array[DIRX] = Rout;

	//Rin_array[DIRY] =0.0;     //don't need to set these since not used anywhere for defcoord = LOGRSINTH  
	//Rout_array[DIRY] = M_PI;
	//Rin_array[DIRZ] = 0.;
	//Rout_array[DIRZ] = 2. * M_PI;

	// define coordinate type
	defcoord = LOGRSINTH;   //logarithmic in r and non-uniform in theta
  hslope = 0.2;           //nonuniforming parameter of the coord. system: hslope = 1 would be uniform in theta

#elif(TESTNUMBER==200)
  Rin=-0.5;
  Rout=1.5 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif(TESTNUMBER==201)
  Rin=-0.5;
  Rout=1.5 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif(TESTNUMBER==202)
  Rin=-1.5;
  Rout=0.5 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif(TESTNUMBER==203)
  Rin=-0.5;
  Rout=1.5 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif(TESTNUMBER==204)
  Rin=-0.5;
  Rout=1.5 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif(TESTNUMBER==205)
  Rin=-0.5;
  Rout=1.5 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif(TESTNUMBER==206)
  Rin=-0.5;
  Rout=1.5 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif(TESTNUMBER==207)
  Rin=-1.0;
  Rout=2.0 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif(TESTNUMBER==208)
  Rin=-1.0;
  Rout=2.0 ;
  // define coordinate type
  //defcoord = 0;
  set_coord_parms();
  Rin_array[DIRX] = Rin;  
  Rout_array[DIRX] = Rout;
#elif( TESTNUMBER == 667 ) // uniform dentisy distribution in slab jet test with cylindrical bc's at the axis: this should be steady state but its apparently not!
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 10.9;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 12.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 13.0;
#elif( TESTNUMBER == 1001 ) // Mignone mildly relativistic blast wave
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#elif( TESTNUMBER == 1002 ) // Mignone very relativistic blast wave
	Rin_array[DIRX] = 0.0;
	Rout_array[DIRX] = 1.0;
	Rin_array[DIRY] =0.0;
	Rout_array[DIRY] = 1.0;
	Rin_array[DIRZ] = 0.;
	Rout_array[DIRZ] = 1.;
#else
	#error "init.c: init_grid(): grid dimensions should be defined for test number" ## XSTRINGIFY( TESTNUMBER )
#endif
  
  //need to set these in order to properly output them to the dump file
  Rin = Rin_array[DIRX];
  Rout = Rout_array[DIRX];
  
  return(0);
}

int init_global(void)
{
  FTYPE R, v0;	
  int whichbcextrap;
  int dir;


  rescaletype=0;


  prfloorcoef[RHO]=RHOMIN;
  prfloorcoef[UU]=UUMIN;



  SAFE=1.3;
  cour = 0.5;

  INVERTFROMAVERAGEIFFAILED = 1;
  LIMIT_AC_PRIM_FRAC_CHANGE = 1;
  MAX_AC_PRIM_FRAC_CHANGE = 0.1;

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  avgscheme=WENO5BND;
  //avgscheme=DONOR;

  do_transverse_flux_integration = 1;
  do_source_integration = 1;
  do_conserved_integration = 1;

  lim = WENO5BND;
  //lim = WENO3;
  //  lim = PARA;
  //lim = DONOR;
  TIMEORDER=4;
  DOENOFLUX = ENOFINITEVOLUME;
  //DOENOFLUX = NOENOFLUX;
  fluxmethod= HLLFLUX; //LAXFFLUX;
  //fluxmethod= LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  //UTOPRIMVERSION=UTOPRIM5D1;
 	UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;
#elif(EOMTYPE==EOMFFDE)
  avgscheme=DONOR;
  do_transverse_flux_integration = 0;
  do_source_integration = 0;
  do_conserved_integration = 0;

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


  
  if( WHICH_INITVEL == VEL3 ) {
    //extrapolate linearly 3-velocity in BC's
    whichbcextrap = BCEXTRAP_VEL3;
  }
  else {
    //extrapolate linearly WHICHVEL
    whichbcextrap = BCEXTRAP;
  }
  
	if( TESTNUMBER == 105 ) {
		//ultrarelativistic analog of the Noh problem requires handling of large lorentz factors
		GAMMAMAX = 1.0e10;
		GAMMAFAIL = 2.0e10;
	} 
	else {
		//leave to defaults from init_defglobal()
		GAMMAMAX = 100.0;
		GAMMAFAIL = 100.0 * GAMMAMAX;
	}

	RHOMIN=0.0;  //changed to lower values by atch
  UUMIN=0.0;

  ranc(7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  cooling=0;

	//set up defaults for tests	

	//default value of gamma
	gam = 1.4; 

	//default boundary condition types:  outflow on all boundaries
	//need to set prior to setting specific cases so that specific cases overwrite the 
	//defaults if need to 
	for( dir = 1; dir < NDIM; dir++ ) {  //set all boundaries to OUTFLOW by default
				BCtype[BC_TYPE_UP_INDEX(dir)] = OUTFLOW;  
				BCtype[BC_TYPE_DN_INDEX(dir)] = OUTFLOW;
	}

	/* output choices */
  tf = 1e4;
  //SAFE = 1.3;
  
#if( TESTNUMBER == -3 )
  BCtype[X_UP] = whichbcextrap;
  BCtype[X_DN] = ASYMM;

  /* output choices */
  tf = 1e4;
  SAFE = 1.0;
  
  gam = 1.4; // define this in pre_init_specific_init, too
#elif( TESTNUMBER == -2 )
  BCtype[X_UP] = whichbcextrap;
  BCtype[X_DN] = CYLAXIS;

  /* output choices */
  tf = 1e4;
  SAFE = 1.3;
  

  gam = 1.4; // define this in pre_init_specific_init, too

  coordparams.rho0 = 1.;
  coordparams.L0 = 0.;
  coordparams.Omega0 = 1.e-4;
  coordparams.tstart = 0.0; //1 / coordparams.Omega0;

  coordparams.vz0 = 0.0; //2 * coordparams.Omega0 * coordparams.L0;
 
  R=0.01;
  v0 = coordparams.Omega0 * R;
  coordparams.u0 = coordparams.rho0 * v0 * v0 * 0.1 * 100;
#elif( TESTNUMBER == -1 )
  BCtype[X_UP] = whichbcextrap;
  BCtype[X_DN] = CYLAXIS;

  /* output choices */
  tf = 1e4;
  SAFE = 1.3;
  
  gam = 1.4; // define this in pre_init_specific_init, too

  coordparams.rho0 = 1.;
  coordparams.L0 = 0.;
  coordparams.Omega0 = 1.e-4;
  coordparams.tstart = 1 / coordparams.Omega0;

  coordparams.vz0 = 2 * coordparams.Omega0 * coordparams.L0;
 
  R=0.01;
  v0 = coordparams.Omega0 * R;
  coordparams.u0 = coordparams.rho0 * v0 * v0 * 0.1 * 100;
#elif( TESTNUMBER == 0 )
  BCtype[X_UP] = whichbcextrap;
  BCtype[X_DN] = whichbcextrap;

  tf = 1e4;
  //SAFE = 1.0;
  
  gam = 1.4; // define this in pre_init_specific_init, too
#elif( TESTNUMBER == 3 )
  gam = 5./3.;
#elif( TESTNUMBER == 666 )
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = PERIODIC;
  BCtype[Y_DN] = PERIODIC;
#elif( TESTNUMBER == 99 ) //asymmetry test
  BCtype[X_UP] = ASYMM;
  BCtype[X_DN] = ASYMM;
#elif( TESTNUMBER == 8 ) //blast wave
  BCtype[X_UP] = ASYMM;
  BCtype[X_DN] = ASYMM;
#elif( TESTNUMBER == 14 ) //spherical Noh, 1st quadrant
  BCtype[X_UP] = FIXED;
  BCtype[X_DN] = ASYMM;
  BCtype[Y_UP] = FIXED;
  BCtype[Y_DN] = ASYMM;
	gam = 5./3.;
#elif( TESTNUMBER == 15 ) //spherical Noh, all four quadrants
  BCtype[X_UP] = FIXED;
  BCtype[X_DN] = FIXED;
  BCtype[Y_UP] = FIXED;
  BCtype[Y_DN] = FIXED;
	gam = 5./3.;
#elif( TESTNUMBER >= 16 && TESTNUMBER <= 21 )   //2D Riemann problems from L&W
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
	gam = 1.4;
#elif( TESTNUMBER == 22) //Implosion, L&W
  BCtype[X_UP] = ASYMM;
  BCtype[X_DN] = ASYMM;
  BCtype[Y_UP] = ASYMM;
  BCtype[Y_DN] = ASYMM;
	gam = 1.4;
#elif( TESTNUMBER == 23 )  //Explosion, L&W
  BCtype[X_UP] = OUTFLOWNOINFLOW;
  BCtype[X_DN] = ASYMM;
  BCtype[Y_UP] = OUTFLOWNOINFLOW;
  BCtype[Y_DN] = ASYMM;
	gam = 1.4;
#elif( TESTNUMBER == 24 )  //Smooth problem, L&W
  BCtype[X_UP] = PERIODIC;
  BCtype[X_DN] = PERIODIC;
  BCtype[Y_UP] = PERIODIC;
  BCtype[Y_DN] = PERIODIC;
	gam = 1.4;
#elif( TESTNUMBER == 25 )  //Odd-even decoupling
  BCtype[X_UP] = ASYMM;
  BCtype[X_DN] = ASYMM;
  BCtype[Y_UP] = PERIODIC;
  BCtype[Y_DN] = PERIODIC;
	gam = 1.4;
#elif( TESTNUMBER == 26 )  //Stationary vortex, Gresho problem
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
	gam = 1.4;
#elif( TESTNUMBER == 27 )  //Moving vortex, Gresho problem
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
	gam = 1.4;
#elif( TESTNUMBER == 28 )  // RT instability
  BCtype[X_UP] = ASYMM;
	BCtype[X_DN] = ASYMM;
  BCtype[Y_UP] = ASYMM;
  BCtype[Y_DN] = ASYMM;
  gam = 1.4;
#elif( TESTNUMBER == 29 )  // smooth problem 4.1 from L&W
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 1.4;
#elif( TESTNUMBER == 30 )  //Smooth sound wave problem, L&W (2D, arbitrary angle)
  BCtype[X_UP] = PERIODIC;
  BCtype[X_DN] = PERIODIC;
  BCtype[Y_UP] = PERIODIC;
  BCtype[Y_DN] = PERIODIC;
	gam = 1.4;
#elif( TESTNUMBER == 31 )  //Smooth 1d sound wave problem, L&W (1D)
  BCtype[X_UP] = PERIODIC;
  BCtype[X_DN] = PERIODIC;
  BCtype[Y_UP] = PERIODIC;
  BCtype[Y_DN] = PERIODIC;
	gam = 1.4;
#elif( TESTNUMBER == 32 )  //Smooth sound wave problem (2D, 2:1 ratio)
  BCtype[X_UP] = PERIODIC;
  BCtype[X_DN] = PERIODIC;
  BCtype[Y_UP] = PERIODIC;
  BCtype[Y_DN] = PERIODIC;
	gam = 5./3.;
#elif( TESTNUMBER == 33 )  //Smooth density wave problem (2D, 2:1 ratio)
  BCtype[X_UP] = PERIODIC;
  BCtype[X_DN] = PERIODIC;
  BCtype[Y_UP] = PERIODIC;
  BCtype[Y_DN] = PERIODIC;
	gam = 5./3.;
#elif( TESTNUMBER == 49 )
  gam = 5./3.;
  BCtype[X_UP] = PERIODIC;
  BCtype[X_DN] = PERIODIC;
#elif( TESTNUMBER == 51 )  // test #1 from RAM paper
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 52 )  // test #1 from RAM paper boosted
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 101 )  // test #1 from RAM paper
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 102 )  // test #2 from RAM paper
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 103 )  // test #3 from RAM paper
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 4./3.;
#elif( TESTNUMBER == 104 )  // test #4 from RAM paper
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 105 )  // test #5 from RAM paper
  BCtype[X_UP] = ASYMM;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 4./3.;
#elif( TESTNUMBER == 106 )  // test #5 from RAM paper
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 107 )  // single left-going shock from test #3 from RAM paper
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 4./3.;
#elif( TESTNUMBER == 151 )  // shock tube test from RAM paper, fig. 8
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 152 )  // Slab jet from RAM
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = ASYMM;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = JETINJECTION;  //should set up the jet; there will partially be outflow bc, partially fixed, see the init_ds_and_vels()
  gam = 5./3.;
#elif( TESTNUMBER == 153 )  // Bondi flow
  BCtype[X_UP] = BONDIMDOTOUTFLOW;
  BCtype[X_DN] = BONDIINTOUTFLOW;
  BCtype[Y_UP] = POLARAXIS;
  BCtype[Y_DN] = POLARAXIS;
  BCtype[Z_UP] = PERIODIC;
  BCtype[Z_DN] = PERIODIC;
  //gam = gam_bondi;  //set it to the value read in from the file with the analytical solution
	//gamma will be set in init_grid
#elif( TESTNUMBER == 154 )  // Torus problem
  BCtype[X_UP] = FIXEDOUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = POLARAXIS;
  BCtype[Y_DN] = POLARAXIS;
  BCtype[Z_UP] = PERIODIC;
  BCtype[Z_DN] = PERIODIC;
  gam = 4./3.;
#elif( TESTNUMBER == 1001 )  // Mignone mild blast wave
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 1002 )  // Mignone very rel blast wave
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = OUTFLOW;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;
  gam = 5./3.;
#elif( TESTNUMBER == 667 )  // uniform density distribution in cyl.coords -- not stationary?
  BCtype[X_UP] = OUTFLOW;
  BCtype[X_DN] = ASYMM;
  BCtype[Y_UP] = OUTFLOW;
  BCtype[Y_DN] = OUTFLOW;  //should set up the jet; there will partially be outflow bc, partially fixed, see the init_ds_and_vels()
  gam = 5./3.;
//#else
//	#error "Boundary conditions have to be specified for test number TESTNUMBER"
#endif

  return(0);

}

int init_grid_post_set_grid(void)
{

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
  int pl;


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

#if(EOMTYPE==EOMGRMHD)
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
  
  // copy over initial solution as analytic solution
  // NEEDED FOR BOUND in case uses panalytic                                                                                                                                                            
  FULLLOOP{
    PLOOP(pl) panalytic[i][j][k][pl]=p[i][j][k][pl];

    // 0 out these things so dump files are readable by SM                                                                                                                                              
    PLOOP(pl) udump[i][j][k][pl]=0.0;
#if(CALCFARADAYANDCURRENTS)
    DLOOPA(pl) jcon[i][j][k][pl]=0.0;
    for(pl=0;pl<NUMFARADAY;pl++) fcon[i][j][k][pl]=0.0;
#endif
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
  if(MAGTEST( TESTNUMBER )){
    if(!MAGTESTNOPOT( TESTNUMBER )){
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
    }
  }
  else{
    LOOPF{
      p[i][j][k][B1]=p[i][j][k][B2]=p[i][j][k][B3]=0.0;
     }
  }




  // copy over initial solution as analytic solution
  // NEEDED FOR BOUND in case uses panalytic                                                                                                                                                            
  FULLLOOP{
    PLOOP(pl) panalytic[i][j][k][pl]=p[i][j][k][pl];

    // 0 out these things so dump files are readable by SM                                                                                                                                              
    PLOOP(pl) udump[i][j][k][pl]=0.0;
#if(CALCFARADAYANDCURRENTS)
    DLOOPA(pl) jcon[i][j][k][pl]=0.0;
    for(pl=0;pl<NUMFARADAY;pl++) fcon[i][j][k][pl]=0.0;
#endif
  }

  /////////////////////////////
  //
  // Fixup and Bound variables since some primitive quantities may have changed
  // These may be used to define vector potential below
  // Also setup pre_fixup() type quantities
  //
  //
  // BOUND AGAIN IN CASE USING PANALYTIC TO BOUND!
  //
  /////////////////////////////// 
  trifprintf("Fixup and Bound #2\n");

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


  return(0);


}




// assumes normalized density
//GODMARK: this program is not supposed to be changing anything == atch
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  
  //int pl;
  //struct of_geom realgeom,geom;
  //FTYPE pratm[NPR];


  //get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
  //set_atmosphere(0,WHICHVEL,&realgeom,pratm); // set velocity in chosen WHICHVEL frame in any coordinate system

  //if(pr[RHO]<pratm[RHO]){
  //  PLOOP(pl) pr[pl]=pratm[pl];
  //}
  

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}

int bounds_generate( int i, int j, int k, FTYPE *prim )
{
  int init_dsandvels(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *p);
	
	int initreturn;
	int whichvel;
	int whichcoord;
	FTYPE old_t;
	int pl;

	old_t = t;  //save the current value of tiem
  
	t = tsteppartf;  //modify the current value of time to point to the end of the current substep so that the boundaries are OK

	initreturn=init_dsandvels( &whichvel, &whichcoord, i, j, k, prim ); // request densities for all computational centers

	t = old_t;  //reset the global time variable to what it was before

	if( initreturn > 0 ) return( initreturn );

	// transform from whichcoord to MCOORD
	if (bl2met2metp2v(whichvel,whichcoord,prim, i,j,k) >= 1)  //SASMARK: replaced p[i][j][k] with prim (luckily, the transformation appears not to be doing anything for the Noh problem since the domain is [0,1]x[0,1]
		FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);




	return( initreturn );
}

// unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *prim)
{
	extern int get_compdimen(int *numdirs, long long *itemp);
	void 	init_torus( int *whichcoord, int *whichvel, int i, int j, int k, FTYPE *pr );
	int dir, numdirs;
	long long itemp;
	FTYPE X[NDIM],V[NDIM], r, th;
	int ti, tj;
	struct of_geom geom;
	struct of_geom geom_prim;

	FTYPE R, z;
	FTYPE sine, cosine;

	int p_no;
	int dorescaletime = 1;

#if(TESTNUMBER == 30 || TESTNUMBER == 31 || TESTNUMBER == 32 || TESTNUMBER == 33)
	FTYPE delta_ampl;
	FTYPE k_vec_x, k_vec_y, k_vec_len;
	FTYPE myrho, myu, mycs, myv;
	FTYPE delta_rho;
  FTYPE cosa, sina;
#endif

#if( TESTNUMBER < 200 || TESTNUMBER > 250 ) //not a force-free test
  FTYPE rho0 = 1.0, v0p = 1.e-4, v0, u0, den, tstart, tau, vz0, Omega0, L0;
  FTYPE x0 = 0.5, x1 = x0;
  FTYPE rhol, rhor, rhom, ul, ur, um, pl, pr, pm, uly = 0.0, ury = 0.0;
  FTYPE Bxl = 0.,Byl=0.,Bzl=0.,Bxr=0.,Byr=0.,Bzr=0.;
  FTYPE Bl[NDIM], Br[NDIM];
  int sign = 1;


  FTYPE prinit[2][2][NPR];  //array to hold the initial primitives
#endif

#if( FFDETEST(TESTNUMBER) ) //force-free tests
  extern void EBtopr(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern void EBtopr_2(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  void vbtopr(FTYPE *vcon,FTYPE *bcon,struct of_geom *geom, FTYPE *pr);
  void computeKK(FTYPE *pr, struct of_geom *geom, FTYPE *KK);
  void EBvetatopr(FTYPE *Econ, FTYPE *Bcon, FTYPE *veta, struct of_geom *geom, FTYPE *pr);

  FTYPE E[NDIM],B[NDIM];
  FTYPE x0, dx0;
  FTYPE bcon[NDIM],vcon[NDIM],econ[NDIM];
  FTYPE phi0;
  FTYPE KK;
  FTYPE B0;
  FTYPE muf;
#endif


#if( FFDETEST(TESTNUMBER) ) //force-free tests
  gset(0,CARTMINKMETRIC,i,j,k,&geom);
  get_geometry(i, j, k, CENT, &geom_prim); // true coordinate system
#else //not a force-free test
  //initially set the primitive quantities to zero; these quantities are only for 2D Riemann problems
  for( ti = 0; ti <= 1; ti++) {
	  for( tj = 0; tj <= 1; tj++) { 
		  PLOOP( p_no ) prinit[ti][tj][p_no] = 0.0;
	  }
  }

  DIMENLOOP(dir) {  //zero out magnetic field by default
	  Bl[dir] = 0.;
	  Br[dir] = 0.;
  }
#endif

  coord(i, j, k, CENT, X);  //get coordinates of (i, j, k) point in prime (internal) coordinate system, X
  bl_coord(X, V);           //get these coords in the orthonormal basis, V

  r  = R = V[DIRX];  //the coordinate along the chosen direcion of x-axis
  th = z = V[DIRY];

  *whichvel=WHICH_INITVEL;

  PLOOP( p_no ) {  //init the primitives with zero
	  prim[p_no] = 0.0;
  }

#if( TESTNUMBER == -3 )    
  *whichcoord=CARTMINKMETRIC;
  dorescaletime = 0;

  v0 =v0p * 0.01;
  u0 = rho0 * v0 * v0 * 0.1;
 
  // the below realizes Ramesh's diverging flow solution -- on half of the grid with the anti-symmetry condition?
  prim[RHO]= rho0;
  prim[UU] = u0;
  prim[U1] = (v0p * R);
  prim[U2] = 0.0;
  prim[U3] = 0.0;
  return( 0 );
#endif

#if( TESTNUMBER == -2 ) 
	#error "this test is not set up"  //need to set up the contents fo coordparams struct first
  *whichcoord=CYLMINKMETRIC;
  
	rho0 = coordparams.rho0;
  L0 = coordparams.L0;
  Omega0 = coordparams.Omega0;
  tstart = coordparams.tstart;

  vz0 = coordparams.vz0;

  dorescaletime = 0;
  
  u0 = coordparams.u0;
  
  tau = t + tstart;
  den = 1. + Omega0 * tau;
 
  // the below realizes 1d solution in 2d, no rotation
  prim[RHO]= rho0 / pow( den, 2. );
  prim[UU] = u0 / pow( den, 2. * gam );
  prim[U1] = Omega0 * R / den;
  if( 0. == L0 ) 
    prim[U2] = 0.0;
  else
    prim[U2] = (- Omega0 * L0 * pow( R / L0, 2. ) / den + vz0);
  prim[U3] = 0.0; //Omega0 / den;  //contravariant component = physical / R for NO GRAVITY
  return( 0 );
#endif

#if( TESTNUMBER == -1 )    
	#error "this test is not set up"  //need to set up the contents fo coordparams struct first
  *whichcoord=CYLMINKMETRIC;
  
	rho0 = coordparams.rho0;
  L0 = coordparams.L0;
  Omega0 = coordparams.Omega0;
  tstart = coordparams.tstart;

  vz0 = coordparams.vz0;

  dorescaletime = 0;
  
  u0 = coordparams.u0;
  
  tau = t + tstart;
  den = 1. + pow( Omega0 * tau, 2. );
 
  coord(i, j, k, CENT, X);
  bl_coord_2d(X, &R, &z);  //get (R, z) coords of the cell's center

  // the below realizes Ramesh's rotating flow solution:
  prim[RHO]= rho0 / den;
  prim[UU] = u0 / pow( den, gam );
  prim[U1] = Omega0 * Omega0 * R * tau / den;
  if( 0. == L0 ) 
    prim[U2] = 0.0;
  else
    prim[U2] = (- Omega0 * L0 * pow( R / L0, 2. ) / den + vz0);
  prim[U3] = Omega0 / den;  //contravariant component = physical / R for NO GRAVITY
  return( 0 );
#endif
  
#if( TESTNUMBER == 0 )    
  dorescaletime = 0;
  *whichcoord=CARTMINKMETRIC;

  v0 = v0p * 0.01;
  u0 = rho0 * v0 * v0 * 0.1 / 250;
  tf = 1. / v0p; //set the final time
 
  // the below realizes Ramesh's diverging flow solution:
  prim[RHO]= rho0;
  prim[UU] = u0;
  prim[U1] = v0p * R;
  prim[U2] = 0.0;
  prim[U3] = 0.0;
#endif

#if( TESTNUMBER == 1 )    
  *whichcoord=CARTMINKMETRIC;
  rhol =  1.;
  ul = 0.75;
  pl = 1.;

  rhor =  0.125;
  ur = 0.;
  pr = 0.1;

  tf = 0.2;
  x0 = 0.3;
#endif

#if( TESTNUMBER == 2 )    
  *whichcoord=CARTMINKMETRIC;
  rhol =  1.;
  ul = -2.;
  pl = 0.40;

  rhor =  1.   ;
  ur = 2.;
  pr = 0.4;

  tf = 0.15;
  x0 = 0.5;
#endif

#if( TESTNUMBER == 3 )    
  *whichcoord=CARTMINKMETRIC;
  rhol =  1.;
  ul = 1.;
  //ul = 1.+1.0/3.0;
  pl = 1.e-6;
  //pl = 0.3;

  rhor =  1.;
  ur = -1.;
  //ur = -1.+1.0/3.0;
  pr = 1.e-6;
  //  pr = 1.5;

  tf = 1.; //(decrease the resolution and the final time by a factor of 5)
  x0 = 0.5;
  
  gam = 5./3.;
#endif

#if( TESTNUMBER == 666 )    
  *whichcoord=CARTMINKMETRIC;
  // like diagonal case3
  //  rhol =  0.138;
  //  ul = sqrt(1.206*1.206 + 1.206*1.206);
  //  pl = 0.029;
  //  rhor =  1.5;
  //  ur = sqrt(0.0*0.0 + 1.206*1.206);
  //  pr = 1.5;

  // like diagonal case4
  rhol =  1.1;
  ul = sqrt(0.8939*0.8939 + 0.8939*0.8939);
  pl = 1.1;
  rhor =  1.1;
  ur = sqrt(0.0*0.0 + 0.0*0.0);
  pr = 1.1;

  //  tf = 0.3; //(decrease the resolution and the final time by a factor of 5)
  tf = 0.25;
  x0 = 0.5;
  gam=1.4;
  //gam = 5./3.;
#endif

#if( TESTNUMBER == 4 )    
  *whichcoord=CARTMINKMETRIC;
  rhol =  5.99924;
  ul = 19.5975;
  pl = 460.894;

  rhor =  5.99242;
  ur = -6.19633;
  pr = 46.095;

  x0 = 0.4;
  tf = 0.035;
#endif

#if( TESTNUMBER == 5 )    
  *whichcoord=CARTMINKMETRIC;
  rhol =   1.4;
  ul = 0.0;
  pl = 1.0;

  rhor =  1.0;
  ur = 0.0;
  pr = 1.0;

  x0 = 0.5;
  tf = 2.0;
#endif

#if( TESTNUMBER == 6 )    
  *whichcoord=CARTMINKMETRIC;
  rhol =   1.4;
  ul = 0.1;
  pl = 1.0;

  rhor =  1.0;
  ur = 0.1;
  pr = 1.0;

  x0 = 0.5;
  tf = 2.0;
#endif

#if( TESTNUMBER == 7 )    //peak problem
  *whichcoord=CARTMINKMETRIC;
  rhol =   0.1261192;
  ul = 8.9047029;
  pl = 782.92899;

  rhor =  6.591493;
  ur = 2.2654207;
  pr = 3.1544874;

  x0 = 0.5;
  tf = 0.0039;
#endif

#if( TESTNUMBER == 99 )     //asymmetry test
  *whichcoord=CARTMINKMETRIC;
  rhol =   0.125;
  ul = 0;
  pl = 0.14;

  rhom =   1.;
  um = 0;
  pm = 1.;

  rhor =  0.125;
  ur = 0;
  pr = 0.14;

  x0 = 0.4 * 1.4;
  x1 = 0.6 * 1.4;
  tf = 2.5;
#endif

#if( TESTNUMBER == 8 )
  *whichcoord=CARTMINKMETRIC;
  rhol =   1;
  ul = 0;
  pl = 1000;
      
  rhom =   1;
  um = 0;
  pm = 0.01;
	    
  rhor =  1;
  ur = 0;
  pr = 100;
	  
  x0 = 0.1;
  x1 = 0.9;
  tf = 0.038;
#endif
			
			
#if( TESTNUMBER == 9 )     //shock entropy wave
  *whichcoord=CARTMINKMETRIC;
  rhol = 3.857143  ;
  ul = 2.629369; 
  pl = 10.33333;

  rhor =  1. + 0.2 * sin( 5.*R );
  ur = 0.0;
  pr = 1.0;

  x0 = -4.;
  tf = 1.8;
#endif

#if( TESTNUMBER == 10 )    //test 3a
  *whichcoord=CARTMINKMETRIC;
  rhol =  1.;
  ul = -19.59745;
  pl = 1000.;

  rhor =  1.;
  ur = -19.59745;
  pr = 0.01;

  tf = 0.012;
  x0 = 0.8;

  //coordparams.timescalefactor = 1e8;
#endif

#if( TESTNUMBER == 11 )    //Lax test
  *whichcoord=CARTMINKMETRIC;
  rhol =  0.445;
  ul = 0.698;
  pl = 3.528;

  rhor =  0.5;
  ur = 0.0;
  pr = 0.571;

  tf = 0.16;
  x0 = 0.0;

  //coordparams.timescalefactor = 1e8;
#endif

#if( TESTNUMBER == 12 )    //Lax test high res
  *whichcoord=CARTMINKMETRIC;
  rhol =  0.445;
  ul = 0.698;
  pl = 3.528;

  rhor =  0.5;
  ur = 0.0;
  pr = 0.571;

  tf = 0.16;
  x0 = 0.0;

  //coordparams.timescalefactor = 1e8;
#endif

#if( TESTNUMBER == 13 )    //Jon's smoothed version of double rarefaction problem for testing the symmetry issues
  *whichcoord=CARTMINKMETRIC;
  rhol =  1.;
  ul = -2.;
  pl = 0.40;

  rhor =  1.   ;
  ur = 2.;
  pr = 0.4;

  tf = 0.15;
  //  tf = 0.3;
  x0 = 0.5;
#endif

#if( TESTNUMBER == 14 || TESTNUMBER == 15 )    
  *whichcoord=CARTMINKMETRIC;
  tf = 2.; 
  
  gam = 5./3.;

	r = sqrt(V[DIRX]*V[DIRX] + V[DIRY]*V[DIRY]);

	prim[RHO] = 1. + t / r / coordparams.timescalefactor;
	prim[UU] = 1.e-6 * ( 1. + gam * t / r / coordparams.timescalefactor ) / ( gam - 1. );  //pressure equal to 1.e-6, need to divide by (gam-1) to get internal energy value
	prim[UU + DIRX] =  - V[DIRX] / r;  //radially directed veloicity of magnitude one
	prim[UU + DIRY] =  - V[DIRY] / r;

#endif

#if( TESTNUMBER == 16 )   //case 3; 2D Riemann problems from L&W
	prinit[0][0][UU] = 0.3 / (gam - 1);
	prinit[0][0][RHO] = 0.5323;
	prinit[0][0][UU+DIRX] = 1.206;
	prinit[0][0][UU+DIRY] = 0.0;
	
	prinit[0][1][UU] = 0.029 / (gam - 1);
	prinit[0][1][RHO] = 0.138;
	prinit[0][1][UU+DIRX] = 1.206;
	prinit[0][1][UU+DIRY] = 1.206;

	prinit[1][0][UU] = 1.5 / (gam - 1);
	prinit[1][0][RHO] = 1.5;
	prinit[1][0][UU+DIRX] = 0.0;
	prinit[1][0][UU+DIRY] = 0.0;

  prinit[1][1][UU] = 0.3 / (gam - 1);
	prinit[1][1][RHO] = 0.5323;
	prinit[1][1][UU+DIRX] = 0.0;
	prinit[1][1][UU+DIRY] = 1.206;

	tf = 0.3;
#endif 

#if( TESTNUMBER == 17 )   //case 4; 2D Riemann problems from L&W
	prinit[0][0][UU] = 0.35 / (gam - 1);
	prinit[0][0][RHO] = 0.5065;
	prinit[0][0][UU+DIRX] = 0.8939;
	prinit[0][0][UU+DIRY] = 0.0;
	
	prinit[0][1][UU] = 1.1 / (gam - 1);
	prinit[0][1][RHO] = 1.1;
	prinit[0][1][UU+DIRX] = 0.8939;
	prinit[0][1][UU+DIRY] = 0.8939;

	prinit[1][0][UU] = 1.1 / (gam - 1);
	prinit[1][0][RHO] = 1.1;
	prinit[1][0][UU+DIRX] = 0.0;
	prinit[1][0][UU+DIRY] = 0.0;

  prinit[1][1][UU] = 0.35 / (gam - 1);
	prinit[1][1][RHO] = 0.5065;
	prinit[1][1][UU+DIRX] = 0.0;
	prinit[1][1][UU+DIRY] = 0.8939;

	tf = 0.25;
#endif 

#if( TESTNUMBER == 18 )   //case 6; 2D Riemann problems from L&W
	prinit[0][0][UU] = 1.0 / (gam - 1);
	prinit[0][0][RHO] = 2.0;
	prinit[0][0][UU+DIRX] = 0.75;
	prinit[0][0][UU+DIRY] = 0.5;
	
	prinit[0][1][UU] = 1.0 / (gam - 1);
	prinit[0][1][RHO] = 1.0;
	prinit[0][1][UU+DIRX] = -0.75;
	prinit[0][1][UU+DIRY] = 0.5;

	prinit[1][0][UU] = 1.0 / (gam - 1);
	prinit[1][0][RHO] = 1.0;
	prinit[1][0][UU+DIRX] = 0.75;
	prinit[1][0][UU+DIRY] = -0.5;

  prinit[1][1][UU] = 1.0 / (gam - 1);
	prinit[1][1][RHO] = 3.0;
	prinit[1][1][UU+DIRX] = -0.75;
	prinit[1][1][UU+DIRY] = -0.5;

	tf = 0.3;
#endif 

#if( TESTNUMBER == 19 )   //case 12; 2D Riemann problems from L&W
	prinit[0][0][UU] = 1.0 / (gam - 1);
	prinit[0][0][RHO] = 1.0;
	prinit[0][0][UU+DIRX] = 0.7276;
	prinit[0][0][UU+DIRY] = 0.0;
	
	prinit[0][1][UU] = 1.0 / (gam - 1);
	prinit[0][1][RHO] = 0.8;
	prinit[0][1][UU+DIRX] = 0.0;
	prinit[0][1][UU+DIRY] = 0.0;

	prinit[1][0][UU] = 0.4 / (gam - 1);
	prinit[1][0][RHO] = 0.5313;
	prinit[1][0][UU+DIRX] = 0.0;
	prinit[1][0][UU+DIRY] = 0.0;

  prinit[1][1][UU] = 1.0 / (gam - 1);
	prinit[1][1][RHO] = 1.0;
	prinit[1][1][UU+DIRX] = 0.0;
	prinit[1][1][UU+DIRY] = 0.7276;

	tf = 0.25;
#endif 

#if( TESTNUMBER == 20 )   //case 15; 2D Riemann problems from L&W
	prinit[0][0][UU] = 0.4 / (gam - 1);
	prinit[0][0][RHO] = 0.5197;
	prinit[0][0][UU+DIRX] = -0.6259;
	prinit[0][0][UU+DIRY] = -0.3;
	
	prinit[0][1][UU] = 0.4 / (gam - 1);
	prinit[0][1][RHO] = 0.8;
	prinit[0][1][UU+DIRX] = 0.1;
	prinit[0][1][UU+DIRY] = -0.3;

	prinit[1][0][UU] = 1.0 / (gam - 1);
	prinit[1][0][RHO] = 1.0;
	prinit[1][0][UU+DIRX] = 0.1;
	prinit[1][0][UU+DIRY] = -0.3;

  prinit[1][1][UU] = 0.4 / (gam - 1);
	prinit[1][1][RHO] = 0.5313;
	prinit[1][1][UU+DIRX] = 0.1;
	prinit[1][1][UU+DIRY] = 0.4276;

	tf = 0.2;
#endif 

#if( TESTNUMBER == 21 )   //case 17; 2D Riemann problems from L&W
	prinit[0][0][UU] = 1.0 / (gam - 1);
	prinit[0][0][RHO] = 2.0;
	prinit[0][0][UU+DIRX] = 0.0;
	prinit[0][0][UU+DIRY] = -0.3;
	
	prinit[0][1][UU] = 0.4 / (gam - 1);
	prinit[0][1][RHO] = 1.0625;
	prinit[0][1][UU+DIRX] = 0.0;
	prinit[0][1][UU+DIRY] = 0.2145;

	prinit[1][0][UU] = 1.0 / (gam - 1);
	prinit[1][0][RHO] = 1.0;
	prinit[1][0][UU+DIRX] = 0.0;
	prinit[1][0][UU+DIRY] = -0.4;

  prinit[1][1][UU] = 0.4 / (gam - 1);
	prinit[1][1][RHO] = 0.5197;
	prinit[1][1][UU+DIRX] = 0.0;
	prinit[1][1][UU+DIRY] = -1.1259;

	tf = 0.3;
#endif 

#if( TESTNUMBER == 22 )   //Implosion, L&W
  *whichcoord=CARTMINKMETRIC;
	if( fabs(V[DIRX] + V[DIRY]) < 0.15 - 1.e-10 ) {  //interior grid cells
		prim[RHO] = 0.125;
		prim[UU] = 0.14 / (gam - 1);
	}
	else if( fabs(V[DIRX] + V[DIRY]) >= 0.15 - 1.e-10 &&  fabs(V[DIRX] + V[DIRY]) <= 0.15 + 1.e-10 ) {  //boundary grid cells -- take an average
		prim[RHO] = (0.125 + 1.0) / 2.;
		prim[UU] = (0.14 + 1.0)/ 2. / (gam - 1) ;
	}
	else {  //outer grid cells
		prim[RHO] = 1.0;
		prim[UU] = 1.0 / (gam - 1);
	}
	tf = 2.5;
#endif 

#if( TESTNUMBER == 23 )   //Explosion, L&W   //SASMARK need to smooth the interface
  *whichcoord=CARTMINKMETRIC;
	if( V[DIRX] * V[DIRX] + V[DIRY] * V[DIRY] < 0.4 * 0.4 ) {
		prim[RHO] = 1.0;
		prim[UU] = 1.0 / (gam - 1);
	}
	else {
		prim[RHO] = 0.125;
		prim[UU] = 0.1 / (gam - 1);
	}
	tf = 3.2;
#endif 

#if( TESTNUMBER == 24 )   //Smooth periodic problem, L&W
  *whichcoord=CARTMINKMETRIC;
	prim[U1] = 1.;
	prim[U2] = -0.5;
	prim[RHO] = 1. + 0.2 * sin( M_PI * ( V[1] + V[2] - t * ( prim[U1] + prim[U2] ) ) );
	prim[UU] = 1. / (gam - 1);
	tf = 4.;
#endif

#if( TESTNUMBER == 25 )     //blast wave in 2d (odd-even decoupling)
  *whichcoord=CARTMINKMETRIC;
  rhol =   1;
  ul = 0;
  pl = 1000;

  rhom =   1;
  um = 0;
  pm = 0.01;

  rhor =  1;
  ur = 0;
  pr = 100;

  x0 = 0.1;
  x1 = 0.9;
  tf = 0.038;
#endif

#if( TESTNUMBER == 26 || TESTNUMBER == 27 )   //Stationary/Moving Gresho vortex problem, assume density is 1. and gam = 1.4
  *whichcoord=CARTMINKMETRIC;
	r = sqrt( V[DIRX] * V[DIRX] + V[DIRY] * V[DIRY] );
	sine = V[DIRY] / r;
	cosine = V[DIRX] / r;

	prim[RHO] = 1.; //SASMARK guess the value for density, not given in L&W
	if( r < 0.2 ){
		prim[UU] = 5+ 25./2. * r * r;
		prim[UU+DIRX] = - 5. * r * sine;
		prim[UU+DIRY] = 5. * r * cosine;
	}
	else if( r < 0.4 ){ 
		prim[UU] = 9. - 4. * log( 0.2 ) + 25./2.*r*r - 20. * r + 4. * log( r );
		prim[UU+DIRX] = - ( 2. - 5. * r) * sine;
		prim[UU+DIRY] =  (2. - 5. * r) * cosine;
	}
	else {
		prim[UU] = 3. + 4. * log( 2.0 );
		prim[UU+DIRX] = 0.0;
		prim[UU+DIRY] = 0.0;
	}

	prim[UU] /= ( gam - 1. );

	tf = 3;
#endif 

#if( TESTNUMBER == 27 )   //Moving Gresho vortex problem, assume density is 1. and gam = 1.4
	prim[UU+DIRX] += 1.;
#endif

#if( TESTNUMBER == 28 )   // Rayleigh-Taylor instability
  *whichcoord=UNIGRAVITY;

  // not moving
  prim[U1+DIRX-1] = 0.0;
  prim[U1+DIRY-1] = 0.01 * (1 + cos(4.0*M_PI*V[DIRX])) * (1 + cos(4./3.*M_PI*V[DIRY]))/4.;  //corrected: cos(3*pi*y) -> cos(4/3*pi*y) to have only a single mode: research/simu/rt.nb
  //prim[U1+DIRY-1] = 0.01 * (1 + cos(4.0*M_PI*V[DIRX])) * (1 + cos(4./3.*M_PI*V[DIRY]))/4.;   //wrong because not a single mode in the y-dir
  prim[U1+DIRZ-1] = 0.0;

  // mass on top is heavier
  // note in LW2003 they say cos, but should be sin
  /*if( V[DIRY] < 0.5 + 0.01*cos(6.0*M_PI*V[DIRX]) ){  //atch modified sin -> cos; compute only half of the mushroom, use mirror bc's to get the full one
    prim[RHO] = 1.;
  }
  else{
    prim[RHO] = 2.;
  }
  */

	if( V[DIRY] > 0. ){  //from http://www.astro.princeton.edu/~jstone/tests/rt/rt.html
    prim[RHO] = 2.;  //heavy on top of light
	}
  else{
    prim[RHO] = 1.;
  }

	// hydrostatic pressure
	prim[UU] = (2.5 - 0.1 * prim[RHO] * V[DIRY])/( gam - 1. );

  //tf = 8.5;
	tf = 12.75;
#endif 

#if( TESTNUMBER == 29 )   //Rarefaction from the peak problem
  *whichcoord=CARTMINKMETRIC;
	/*rhol =  0.126119;
	ul = 8.904703;
	pl = 782.929016;

	rhor =  0.122060 ;
	ur = 11.944731;
	pr = 747.877502;

	x0 = 0.5;*/
  tf = 0.0029;
	
	prim[RHO] = test29rho[i];
	prim[UU] = test29p[i] / (gam - 1);
	prim[U1] = test29u[i];
	prim[U2] = 0.0;
#endif

#if( TESTNUMBER == 30 )   //Smooth periodic sound wave problem
    *whichcoord=CARTMINKMETRIC;
	delta_ampl = 3.2e-8;
	k_vec_x = 2 * M_PI/Rout_array[DIRX];  //wavevector
	k_vec_y = 2 * M_PI/Rout_array[DIRY] ;
	k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );

	delta_rho = delta_ampl * cos( k_vec_x * V[1] + k_vec_y * V[2] );
	myrho = 1.;
	myu = myrho / (gam * (gam-1));  //so that mycs is unity
 
	mycs = sqrt(gam * (gam-1) * myu / myrho);  //background sound speed

	//applying the perturbations

	prim[RHO] = myrho + delta_rho;
  prim[UU] = myu + gam * myu * delta_rho / myrho;
 	prim[U1] = delta_rho/myrho * mycs * k_vec_x / k_vec_len;
	prim[U2] = delta_rho/myrho * mycs * k_vec_y / k_vec_len;

	tf = 2. * M_PI / (k_vec_len * mycs);
#endif

#if( TESTNUMBER == 31 )   //Smooth periodic sound wave problem
    *whichcoord=CARTMINKMETRIC;
	delta_ampl = 3.2e-8;
	k_vec_x = 2 * M_PI/Rout_array[DIRX];  //wavevector
	k_vec_y = 0;
	k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );

	delta_rho = delta_ampl * cos( k_vec_x * V[1] + k_vec_y * V[2] );
	myrho = 1.;
 	myu = myrho / (gam * (gam-1));  //so that mycs is unity
 
	mycs = sqrt(gam * (gam-1) * myu / myrho);  //background sound speed

	//applying the perturbations

	prim[RHO] = myrho + delta_rho;
  prim[UU] = myu + gam * myu * delta_rho / myrho;
	prim[U1] = delta_rho/myrho * mycs * k_vec_x / k_vec_len;
	prim[U2] = delta_rho/myrho * mycs * k_vec_y / k_vec_len;

	tf = 1./mycs;
#endif

#if( TESTNUMBER == 32 )   //Smooth periodic sound wave problem - Jim Stone way
    *whichcoord=CARTMINKMETRIC;
	delta_ampl = 3.2e-8;
	k_vec_x = 2 * M_PI/Rout_array[DIRX];  //wavevector
	k_vec_y = 2 * M_PI/Rout_array[DIRY] ;
	k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );

	delta_rho = delta_ampl * cos( k_vec_x * V[1] + k_vec_y * V[2] );
	myrho = 1.;
	myu = myrho / (gam * (gam-1));  //so that mycs is unity, Jim Stone says P = 1/\gamma, which is equivalent to what i do
 
	mycs = sqrt(gam * (gam-1) * myu / myrho);  //background sound speed

	//applying the perturbations

	prim[RHO] = myrho + delta_rho;
  prim[UU] = myu + gam * myu * delta_rho / myrho;
 	prim[U1] = delta_rho/myrho * mycs * k_vec_x / k_vec_len;
	prim[U2] = delta_rho/myrho * mycs * k_vec_y / k_vec_len;

	tf = 2. * M_PI / (k_vec_len * mycs);  //so that the wave travels for one wavelength
#endif

#if( TESTNUMBER == 33 )   //Smooth periodic density wave problem - Jim Stone way
    *whichcoord=CARTMINKMETRIC;
	delta_ampl = 3.2e-8;
	k_vec_x = 2 * M_PI/Rout_array[DIRX];  //wavevector, k_x = 2 * k_y
	k_vec_y = 2 * M_PI/Rout_array[DIRY] ;
	k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );

  //cosine and sine of the angle at which wave is propagating (i.e. pi/2 - angle btw. the diagonal and the x-axis)
  cosa = Rout_array[DIRY] / sqrt( Rout_array[DIRX] * Rout_array[DIRX] + Rout_array[DIRY] * Rout_array[DIRY] );
  sina = Rout_array[DIRX] / sqrt( Rout_array[DIRX] * Rout_array[DIRX] + Rout_array[DIRY] * Rout_array[DIRY] );

	delta_rho = delta_ampl * cos( k_vec_x * V[1] + k_vec_y * V[2] );
	myrho = 1.;
	myu = myrho / (gam * (gam-1));  //so that mycs is unity
 
	mycs = sqrt(gam * (gam-1) * myu / myrho);  //background sound speed
  myv = 1.;  //velocity with a magnitude of 1

	//applying the perturbations

	prim[RHO] = myrho + delta_rho;  //perturbations only to density
  prim[UU] = myu; 
 	prim[U1] = myv * cosa;  
	prim[U2] = myv * sina;

	tf = 2. * M_PI / (k_vec_len * myv);  //so that the wave travels for one wavelength
#endif

#if( TESTNUMBER == 49 )     //1D caustics from Ryu et al. 1993
  *whichcoord=CARTMINKMETRIC;

  //this is not used
  rhol = 1.;
  ul = -sin( 2 * M_PI * R ) / ( 2 * M_PI );
  pl = 1.e-10;
  //end of not used

  rhor =  1.;
  ur = -sin( 2 * M_PI * R ) / ( 2 * M_PI );
  pr = 1.e-10;  //1.e-10;

  x0 = -10;  //so that only the right state is used
  tf = 3.;
#endif


#if( TESTNUMBER == 51 )    // Sod's Riemann problem
  *whichcoord=CARTMINKMETRIC;
  rhol =  1.0;
  ul = 0.0;
	uly = 0.0;
  pl = 1.0;

  rhor =  0.2;
  ur = 0.0;
	ury = 0.0;
  pr = 0.01;

  tf = 0.769231;
  x0 = 0.0;
#endif

#if( TESTNUMBER == 52 )    // Sod's Riemann problem boosted
  *whichcoord=CARTMINKMETRIC;
  rhol =  1.0;
  ul = 100 * 0.288675;
	uly = 0.0;
  pl = 1.0;

  rhor =  0.2;
  ur = 100 * 0.288675;
	ury = 0.0;
  pr = 0.01;

  tf = 0.769231;
  x0 = 0.0;
#endif


#if( TESTNUMBER == 101 )    // 1D Riemann problem #1 from RAM paper
  *whichcoord=CARTMINKMETRIC;
  rhol =  10.0;
  ul = 0.0;
	uly = 0.0;
  pl = 13.33;

  rhor =  1.0;
  ur = 0.0;
	ury = 0.0;
  pr = 1.e-8;

  tf = 0.4;
  x0 = 0.5;

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif

#if( TESTNUMBER == 102 )    // 1D Riemann problem #2 from RAM paper
  *whichcoord=CARTMINKMETRIC;
  rhol =  1.0;
  ul = 0.0;
	uly = 0.0;
  pl = 1000;

  rhor =  1.0;
  ur = 0.0;
	ury = 0.0;
  pr = 1.e-2;

  tf = 0.4;
  x0 = 0.5;

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif

#if( TESTNUMBER == 103 )    // 1D Riemann problem #3 from RAM paper
  *whichcoord=CARTMINKMETRIC;
  rhol = 1.0;
  ul = 0.9;
	uly = 0.0;
  pl = 1.0;

  rhor = 1.0;
  ur = 0.0;
	ury = 0.0;
  pr = 10.;

  tf = 0.4;
  x0 = 0.5;

	*whichvel=VEL3; //the above uses 3-velocities

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif

#if( TESTNUMBER == 104 )    // 1D Riemann problem #4 from RAM paper
  *whichcoord=CARTMINKMETRIC;
  rhol = 1.0;
  ul = 0.0;
	uly = 0.0;
  pl = 1000.0;

  rhor = 1.0;
  ur = 0.0;
	ury = 0.99;
  pr = 0.01;

  tf = 0.4;
  x0 = 0.5;

	*whichvel=VEL3; //the above uses 3-velocities

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif

	// stupid version that sets 3-velocity -- which has itself artificial limits
#if( TESTNUMBER == 105 )    // 1D Riemann problem #5 from RAM paper
  *whichcoord=CARTMINKMETRIC;
  rhol = 1.0;
  //  ul = 1.0 - 1.e-11;
  ul = 1.0 - 1.e-10;
  //ul = 1.0 - 1.e-8;
  //ul = (1.0 - 1.e-10) * 70710.675;
	uly = 0.0;
  pl = 0.003 * (gam - 1);

  rhor = 1.0;
  //ur = 1.0 - 1.e-11;
  ur = 1.0 - 1.e-10;
  // ur = 1.0 - 1.e-8;
  //ur =  (1.0 - 1.e-10) * 70710.675;
	ury = 0.0;
  pr = 0.003 * (gam - 1);

  tf = 2.0;
  x0 = 0.5;   //not a Riemann problem, just bouncing off the boundary.  Can make it into a Riemann problem, just like the Noh problem

	*whichvel=VEL3; //the above uses 3-velocities

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif
#if( TESTNUMBER == 1055555 )    // 1D Riemann problem #5 from RAM paper
  *whichcoord=CARTMINKMETRIC;
  rhol = 1.0;
  ul = 1e7;
  uly = 0.0;
  pl = 0.003 * (gam - 1);

  rhor = rhol;
  ur = ul;
  ury = uly;
  pr = pl;

  tf = 2.0;
  x0 = 0.5;   //not a Riemann problem, just bouncing off the boundary.  Can make it into a Riemann problem, just like the Noh problem

  *whichvel=VELREL4; //the above uses relative 4-velocities

  assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif

#if( TESTNUMBER == 106 )    // 1D Riemann problem #5 from RAM paper
  *whichcoord=CARTMINKMETRIC;
  rhol = 1.0;
  ul = 0.0;
        uly = 0.9;
  pl = 1000.0;

  rhor = 1.0;
  ur = 0.0;
        ury = 0.9;
  pr = 1.e-2;

  tf = 0.6;
  x0 = 0.5;   //not a Riemann problem, just bouncing off the boundary.  Can make it into a Riemann problem, just like the Noh problem

	*whichvel=VEL3; //the above uses 3-velocities

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif

#if( TESTNUMBER == 107 )    // 1D Riemann problem #3 from RAM paper
	//Have found that this single shock generates the same wiggles as the full Riemann problem.
	*whichcoord=CARTMINKMETRIC;
	rhol = 1.0;
	ul = 0.9;
	uly = 0.0;
	pl = 1.0;

	rhor = 0.659660744e+01;  //the right state is taken from the exact Riemann solver
	ur = 0.242538591e+00;
	ury = 0.0;
	pr = 0.177916477e+02;

	tf = 0.4;
	x0 = 0.5;

	*whichvel=VEL3; //the above uses 3-velocities

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif

#if( TESTNUMBER == 151 )   //Two-dimensional shock-tube problem, fig. 8, RAM paper
	//   locations on the grid
	// (0,0) | (1,0)
	// ------+------
	// (0,1) | (1,1)

	prinit[1][0][RHO] = 0.1;
	prinit[1][0][UU+DIRX] = 0.0;
	prinit[1][0][UU+DIRY] = 0.0;
	prinit[1][0][UU] = 0.01 / (gam - 1);

	prinit[0][0][RHO] = 0.1;
	prinit[0][0][UU+DIRX] = 0.99;
	prinit[0][0][UU+DIRY] = 0.0;
	prinit[0][0][UU] = 1.0 / (gam - 1);
	
	prinit[0][1][RHO] = 0.5;
	prinit[0][1][UU+DIRX] = 0.0;
	prinit[0][1][UU+DIRY] = 0.0;
	prinit[0][1][UU] = 1.0 / (gam - 1);

	prinit[1][1][RHO] = 0.1;
	prinit[1][1][UU+DIRX] = 0.0;
	prinit[1][1][UU+DIRY] = 0.99;
  prinit[1][1][UU] = 1.0 / (gam - 1);

	*whichvel=VEL3; //the above uses 3-velocities

	tf = 0.4;
#endif 

#if( TESTNUMBER == 152 ) //slab jet  //SUPERSASMARK
	*whichcoord=CYLMINKMETRIC;
	*whichvel = VEL3;
	tf = 100.;
	if( z < 0.0 && R < 1.0 ) {
		prim[RHO] = 0.01;
		prim[UU]  = 0.000170305;
    prim[UU + DIRX] = 0.0;
    prim[UU + DIRY] = 0.99;    //z-direction velocity
    prim[UU + DIRZ] = 0.0;     //phi-direction velocity
	}
	else if( z < 0.0 ) {  //OUTFLOW BC
		//should already be populated with the outflow boundary condition, so do not do anything here
	}
	else {
		//ambient medium
		prim[RHO] = 1.0;
		prim[UU]  = 0.000170305;
    prim[UU + DIRX] = 0.0;
    prim[UU + DIRY] = 0.0;     //z-direction velocity
    prim[UU + DIRZ] = 0.0;     //phi-direction velocity
	}
	tf = 100.;
#endif

#if( TESTNUMBER == 153 ) //bondi problem; 
	//the initial conditions were generated in BLCOORDS but there is no conversion for radial velocity, density, and internal energy
	//so use KSCOORDS to set them so that the internal conversion does not happen
	*whichcoord=KSCOORDS;
	*whichvel = VEL4;
	prim[RHO] = rho_bondi[i];
	prim[UU]  = u_bondi[i];
	prim[UU + DIRX] = uu1_bondi[i];
	prim[UU + DIRY] = 0.0; 
	prim[UU + DIRZ] = 0.0; 
	dorescaletime = 0;
	tf = 1000.;
#endif


#if( TESTNUMBER == 154 ) //Torus problem
	//the initial conditions were generated in BLCOORDS but there is no conversion for radial velocity, density, and internal energy
	//so use KSCOORDS to set them so that the internal conversion does not happen
	init_torus( whichcoord, whichvel, i, j, k, prim );
	dorescaletime = 0;
  tf = 10.;  //final time for the torus problem;
#endif

#if(TESTNUMBER==200) // Fast wave
  tf = 1;
  DTd=tf/10.0;

  //  tf = 1;
  //  DTd=1E-5;
  

  E[1]=0;
  E[2]=0;
  B[3]=0.0;
  B[1]=0.0;
  x0=0.0;
  dx0=0.1;
  if(r-x0<=-dx0) B[2]=1.0;
  else if((r-x0>-dx0)&&(r-x0<dx0)) B[2]=1.0-(0.3/(2.0*dx0))*((r-x0)+dx0);
  else if(r-x0>=dx0) B[2]=0.7;

  muf=1.0;

  E[3]=1.0-muf*B[2];

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention


  EBtopr(E,B,&geom,prim);
  //EBtopr_2(E,B,&geom,prim);

  //  prim[U1]=0.9;

  //dualfprintf(fail_file,"prim[U1]=%21.15g prim[U2]=%21.15g\n",prim[U1],prim[U2]);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;

#endif
#if(TESTNUMBER==201) // comoving Fast wave (NOT a Komissarov test)
  //tf = 1;
  //  DTd=tf/10.0;

  tf = 1;
  DTd=1E-5;
  

  bcon[3]=0.0;
  bcon[1]=0.0;
  x0=0.0;
  dx0=0.1;
  if(r-x0<=-dx0) bcon[2]=1.0;
  else if((r-x0>-dx0)&&(r-x0<dx0)) bcon[2]=1.0-(0.3/(2.0*dx0))*((r-x0)+dx0);
  else if(r-x0>=dx0) bcon[2]=0.7;

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention


  vcon[1]=0.9999;
  vcon[2]=vcon[3]=0;

  vbtopr(vcon,bcon,&geom,prim);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;
#endif
#if(TESTNUMBER==202) // (nondegenerate) Alfven wave
  tf = 2;
  DTd=tf/10.0;

  bcon[1]=bcon[2]=1.0;
  
  x0=0.0;
  if(r-x0<=-0.1) bcon[3]=1.0;
  else if((r-x0>-0.1)&&(r-x0<0.1)) bcon[3]=1.0+3.0/2.0*((r-x0)+0.1);
  else if(r-x0>=0.1) bcon[3]=1.3;

  econ[2]=econ[3]=0.0;

  // can be + or -
  //#define CONSTECON1 1.3
  //  econ[1]=-sqrt(-CONSTECON1+bcon[1]*bcon[1]+bcon[2]*bcon[3]+bcon[3]*bcon[3]);
  //  econ[1]=1-0.5*bcon[3];
#define CONSTECON1 1.0
  econ[1]=sqrt(-CONSTECON1 + bcon[3]*bcon[3]);
  //  econ[1]=0.0;
  //  econ[1]=-bcon[3];

  vcon[1]=-0.5;
  vcon[2]=vcon[3]=0;

  EBvetatopr(econ, bcon, vcon, &geom, prim);
  //  vbtopr(vcon,bcon,&geom,prim);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;
#endif

#if(TESTNUMBER==203) // Degenerate Alfven wave
  tf = 2;
  DTd=tf/10.0;
  bcon[1]=0.0;

  
  x0=0.0;
  if(r-x0<=-0.1) phi0=0.0;
  else if((r-x0>-0.1)&&(r-x0<0.1)) phi0=5.0/2.0*M_PI*((r-x0)+0.1);
  else if(r-x0>=0.1) phi0=M_PI*0.5;

  bcon[2]=2.0*cos(phi0);
  bcon[3]=2.0*sin(phi0);


  vcon[1]=0.5;
  vcon[2]=vcon[3]=0;

  vbtopr(vcon,bcon,&geom,prim);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;

#endif
#if(TESTNUMBER==204) // Three-wave problem
  tf = .75;
  DTd=tf/10.0;

  x0=0.5;
  if(r<x0){
    B[1]=1.0;
    B[2]=1.5;
    B[3]=3.5;
    E[1]=-1.0;
    E[2]=-0.5;
    E[3]=0.5;
  }
  else{
    B[1]=1.0;
    B[2]=2.0;
    B[3]=2.3;
    E[1]=-1.5;
    E[2]=1.3;
    E[3]=-0.5;
  }

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

  EBtopr(E,B,&geom,prim);
  //EBtopr_2(E,B,&geom,prim);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;

#endif

#if(TESTNUMBER==205) // B^2-E^2<0 problem
  tf = .02;
  DTd=tf/10.0;

  x0=0.5;
  if(r<x0){
    B[0]=0.0;
    B[1]=1.0;
    B[2]=1.0;
    B[3]=1.0;
    E[0]=0.0;
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }
  else{
    B[0]=0.0;
    B[1]=1.0;
    B[2]=-1.0;
    B[3]=-1.0;
    E[0]=0.0;
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

  EBtopr(E,B,&geom,prim);
  //EBtopr_2(E,B,&geom,prim);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;

#endif
#if(TESTNUMBER==206) // smoothed B^2-E^2<0 problem
  tf = .02;
  DTd=tf/10.0;

  x0=0.5;
  if(r-x0<-0.1){
    B[1]=1.0;
    B[2]=1.0;
    B[3]=1.0;
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }
  else if(r-x0>0.1){
    B[1]=1.0;
    B[2]=-1.0;
    B[3]=-1.0;
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }
  else if((r-x0>=-0.1)&&(r-x0<=0.1)){
    B[1]=1.0;
    B[2]=1.0+(r-x0+0.1)*(-2.0/0.2);
    B[3]=1.0+(r-x0+0.1)*(-2.0/0.2);
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention


  EBtopr(E,B,&geom,prim);
  //EBtopr_2(E,B,&geom,prim);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;

#endif

#if(TESTNUMBER==207) // Komissarov 2004 C3.1 Alfven wave
  prim[RHO]=prim[UU]=0;
  prim[U1]=prim[U2]=prim[U3]=0.0;
  prim[B2]=prim[B3]=0;
  prim[B1]=0;
  // PARA generates crap on left side, but wave doesn't move
  // MC does very well
  // no obvious difference between HLL and LAXF
  // Athena1/2 ok
  tf = 2.0;
  DTd=tf/10.0;

  //  B[1]=B[2]=E[3]=E[2]=0;
  B[1]=B[2]=E[3]=1.0;
  E[2]=0;


  if(r<=0.5){
    B[3]=1.0;
  }
  else if(r>=0.2+0.5){
    B[3]=1.3;
  }
  else{
    B[3]=1.0+0.15*(1.0+sin(5.0*M_PI*(r-0.1-0.5)));
  }
  E[1]=-B[3];

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

  EBtopr(E,B,&geom,prim);
  //EBtopr_2(E,B,&geom,prim);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;

#endif
#if(TESTNUMBER==208) // Komissarov 2004 C3.2 Current Sheet
  tf = 1.0;
  DTd=tf/10.0;

  E[1]=E[2]=E[3]=0.0;
  B[3]=0.0;
  B[1]=1.0;

  //B0=0.5; // fine
  B0 = 2.0;

  if(r<=0.5){
    B[2]=B0;
  }
  else{
    B[2]=-B0;
  }


  EBtopr(E,B,&geom,prim);
  //EBtopr_2(E,B,&geom,prim);

  computeKK(prim,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC;

#endif

#if( TESTNUMBER == 1001 )    // Mignone mild blast wave
  *whichcoord=CARTMINKMETRIC;
  rhol = 1.0;
  ul = 0.0;
  uly = 0.0;
  pl = 30.0;
  Bxl=5.0;
  Byl=6.0;
  Bzl=6.0;

  rhor = 1.0;
  ur = 0.0;
  ury = 0.0;
  pr = 1.0;
  Bxr=5.0;
  Byr=0.7;
  Bzr=0.7;


	//assign magnetic field
	Bl[DIRX] = Bxl;
	Bl[DIRY] = Byl;
	Bl[DIRZ] = Bzl;

	Br[DIRX] = Bxr;
	Br[DIRY] = Byr;
	Br[DIRZ] = Bzr;


  tf = 0.4;
  x0 = 0.5;   //not a Riemann problem, just bouncing off the boundary.  Can make it into a Riemann problem, just like the Noh problem

	*whichvel=VEL3; //the above uses 3-velocities

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif

#if( TESTNUMBER == 1002 )    // Mignone strong blast wave
  *whichcoord=CARTMINKMETRIC;
  rhol = 1.0;
  ul = 0.0;
  uly = 0.0;
  pl = 1.0E3;
  Bxl=10.0;
  Byl=7.0;
  Bzl=7.0;

  rhor = 1.0;
  ur = 0.0;
  ury = 0.0;
  pr = 0.1;
  Bxr=10.0;
  Byr=0.7;
  Bzr=0.7;

	//assign magnetic field
	Bl[DIRX] = Bxl;
	Bl[DIRY] = Byl;
	Bl[DIRZ] = Bzl;

	Br[DIRX] = Bxr;
	Br[DIRY] = Byr;
	Br[DIRZ] = Bzr;

  tf = 0.4;
  x0 = 0.5;   //not a Riemann problem, just bouncing off the boundary.  Can make it into a Riemann problem, just like the Noh problem

	*whichvel=VEL3; //the above uses 3-velocities

	assert( coordparams.timescalefactor != 1.0, "For test no. %d timescalefactor should be 1.0 (currently = %g)\n", (int)TESTNUMBER, coordparams.timescalefactor );
#endif



#if( TESTNUMBER == 667 ) //// uniform density distribution in cyl.coords -- not stationary?
	*whichcoord=CYLMINKMETRIC;
	*whichvel = VEL3;
	tf = 100.;

	//ambient medium
	prim[RHO] = 1.0;
	prim[UU]  = 0.000170305;
  prim[UU + DIRX] = 0.0;
  prim[UU + DIRY] = 0.0;     //z-direction velocity
  prim[UU + DIRZ] = 0.0;     //phi-direction velocity
#endif


	//TEST numbering conventions:
	// < 100 & 666 -- non-relativistic tests
	// 101 -- 149 -- 1d relativistic tests
	// 151 -- 199 -- 2d relativistic tests
#if( NONREL1DRIEMANNTEST(TESTNUMBER) || REL1DRIEMANNTEST(TESTNUMBER) || TESTNUMBER == 8 ||TESTNUMBER == 9 || TESTNUMBER == 25 || TESTNUMBER == 99 || TESTNUMBER == 49 )  //1D test problems or 2D blast wave (#25) or play around with 1D problem
	//only for 1d Riemann problems
	//set up the problems to be directed along the chosen x-direction
  if( V[DIRX] < x0 ) {
    prim[RHO] = rhol;
    prim[UU] = u_rho0_p(rhol,pl);
    prim[U1 + DIRX - 1] = ul;
    prim[UU + DIRY] = uly;
    prim[B1 + DIRX - 1] = Bxl;
    prim[B1 + DIRY - 1] = Byl;
    prim[B1 + DIRZ - 1] = Bzl;
    //    prim[U1 + DIRY - 1] = 2.0;// for TESTNUMBER == 5 w/ tangential velocity
  }
  #if( TESTNUMBER == 8 || TESTNUMBER == 25 || TESTNUMBER == 99 )  //the prolems requiring three initial states
  else if( V[DIRX] < x1 ) {
    prim[RHO] = rhom;
    prim[UU] = u_rho0_p(rhom,pm);
    prim[UU + DIRX] = um;
  }
  #endif	  
  else {
    prim[RHO] = rhor;
    prim[UU] = u_rho0_p(rhor,pr);
    prim[UU + DIRX] = ur;
    prim[UU + DIRY] = ury;
    prim[B1 + DIRX - 1] = Bxr;
    prim[B1 + DIRY - 1] = Byr;
    prim[B1 + DIRZ - 1] = Bzr;
    //    prim[U1 + DIRY - 1] = -2.0; // for TESTNUMBER == 5 w/ tangential velocity
  }

#elif( (TESTNUMBER >= 16 && TESTNUMBER <= 21) || (TESTNUMBER == 151) )  //2D test problems
	*whichcoord=CARTMINKMETRIC;
	//get the index of the current quadrant: 
	//   
	// (0,0) | (1,0)
	// ------+------
	// (0,1) | (1,1)
	//
	PLOOP( p_no ) prim[p_no] = prinit[V[DIRX] > x0][V[DIRY] < x1][p_no];
#endif


#if( NONREL1DRIEMANNTEST(TESTNUMBER) || REL1DRIEMANNTEST(TESTNUMBER) )  //1d riemann problems
  if( 0 == i && 0 == j && 0 ==k ) {
      write_riemannproblem_params_to_file( rhol, pl, ul, uly, rhor, pr, ur, ury, Bl, Br, tf, coordparams.timescalefactor, x0, Rin_array[DIRX], Rout_array[DIRX] );
    }
#endif
	//rescale the time to make the problem non-relativistic (this stretches the time and lowers the velocities and internal energies)
	//all test problems
  if( dorescaletime ) {
    prim[U1] /= coordparams.timescalefactor;
		prim[U2] /= coordparams.timescalefactor;
		prim[U3] /= coordparams.timescalefactor;
		prim[UU] /= ( coordparams.timescalefactor * coordparams.timescalefactor );
		prim[B1] /= coordparams.timescalefactor;  //SASMARK: shoud fields scale linearly with timescalefactor so that ratio of energy densities, B^2 / u, is invariant?
		prim[B2] /= coordparams.timescalefactor;
		prim[B3] /= coordparams.timescalefactor;


    tf *= coordparams.timescalefactor;
  }


  //tf /= 10.;

  DTd = tf / 10;
  DTdebug = DTd;

  get_compdimen(&numdirs,&itemp);
  if(numdirs>1) DTi = tf/100.0; 
  else DTi = tf; //output no images
  DTavg = tf;
  DTener = tf/100.0;
	DTr = 30000; //restart period in steps

	//print out the final time only once
	//  if( 0 == i && 0 == j && 0 == k ) trifprintf( "tf = %g\n", tf );
  
//  PLOOP(p_no) {
//	prim[p_no]*=(1.+1.e-13*(ranc(0)-0.5));
//  }

  return(0);

}

//initializes the torus problem according to the harm paper and init.fishmon.c
void 	init_torus( int *whichcoord, int *whichvel, int i, int j, int k, FTYPE *pr )
{
  FTYPE sth, cth;
  FTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom realgeom;
 

  /* for disk interior */
  SFTYPE l, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;
  SFTYPE rmax, lfish_calc(SFTYPE rmax);
  SFTYPE rh;
  //  FTYPE pratm[NPR];
	FTYPE rin;
	FTYPE pr_atmosphere[NPR];

  rin = 3.7;
  l = 3.85;        //= u^t u_\phi "specific angular momentum"; from this u^\phi is obtained
  kappa = 1.e-3 ;  // p = (\Gamma - 1) \kappa \rho^\Gamma
  

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

		RHOMIN = 1.0e-4;  //atmosphere prefactor values
		UUMIN = 1.0e-6;

	  RHOMINLIMIT=1.e-20;  //minimum actual values of density and internal energy {e.g., density_floor = MAX[ RHOMIN * r^{-1.5}, RHOMINLIMIT ] }
		UUMINLIMIT =1.e-20;

		set_atmosphere(-1,WHICHVEL,&realgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
  }
  /* region inside magnetized torus; u^i is calculated in
     Boyer-Lindquist coordinates, as per Fishbone & Moncrief, so it
     needs to be transformed at the end */
  else {
    hm1 = exp(lnh) - 1.;  //= h - 1, where h is specific enthalpy: h = (p + \rho)/\rho = 1 + p/\rho
    rho = pow(hm1 * (gam - 1.) / (kappa * gam), 1. / (gam - 1.));  //obtained from combining the above and below eqs.: \rho = ( (h-1)/(\kappa(\Gamma-1)) )^( 1/(\Gamma-1) )
    u = kappa * pow(rho, gam) / (gam - 1.);  // u = p / (\Gamma - 1) = \kappa \rho^\Gamma/(\Gamma - 1)
    
		ur = 0.;  //purely toroidal motion
    uh = 0.;
    
    /* calculate u^phi */
    expm2chi = SS * SS * DD / (AA * AA * sth * sth);
    up1 = sqrt((-1. + sqrt(1. + 4. * l * l * expm2chi)) / 2.);
    up = 2. * a * r * sqrt(1. + up1 * up1) / sqrt(AA * SS * DD) +
      sqrt(SS / AA) * up1 / sth;
    

		/////////
		//
		// make sure torus densities do not fall below the atmospheric level
		//
		//set_density_floors(&realgeom,pr,prfloor);
    get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
		set_atmosphere( -1, VELREL4, &realgeom, pr_atmosphere );

		if( rho < pr_atmosphere[RHO] )  rho = pr_atmosphere[RHO];
		if(  u  < pr_atmosphere[UU]  )   u  = pr_atmosphere[UU];
		//
		/////
   
    pr[RHO] = rho ;
    pr[UU] = u;
    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = up;

    *whichvel=VEL4;
    *whichcoord=BLCOORDS;
  }

	pr[B1] = 0.0;
	pr[B2] = 0.0;
	pr[B3] = 0.0;
}



void write_riemannproblem_params_to_file( FTYPE rhol, FTYPE pl, FTYPE ul, FTYPE uly, FTYPE rhor, FTYPE pr, FTYPE ur, FTYPE ury, FTYPE *Bl, FTYPE *Br, 
																				 FTYPE tf, FTYPE timescalefactor, FTYPE x0, FTYPE a, FTYPE b ) 
{
	FILE *fp;

	FTYPE prim_l[NPR], prim_r[NPR];

	prim_l[RHO] = rhol;
	prim_l[U1] = ul;
	prim_l[U2] = uly;
	prim_l[U3] = 0;
	prim_l[UU] = pl / (gam - 1);
	prim_l[B1] = Bl[DIRX];
	prim_l[B2] = Bl[DIRY];
	prim_l[B3] = Bl[DIRZ];

	prim_r[RHO] = rhor;
	prim_r[U1] = ur;
	prim_r[U2] = ury;
	prim_r[U3] = 0;
	prim_r[UU] = pr / (gam - 1);
	prim_r[B1] = Br[DIRX];
	prim_r[B2] = Br[DIRY];
	prim_r[B3] = Br[DIRZ];

	if( REL1DRIEMANNTEST(TESTNUMBER) == 1 ) {
		//First write grid parameters for the Riemann Solver -- need to feed them into it separately from stdin ( e.g. ./riemannrmhd <ParamsInput.txt &> RiemannOutput.txt )
		fp = fopen( "ParamsInput.txt", "wt" );

		if( fp == NULL ) {
			dualfprintf( fail_file, "Could not open output file for writing Relativistic Riemann problem parameters file, %s\n", "ParamsInput.txt" );
			return;
		}

		//Output the test info: internal testnumber inside of Riemann Solver (always 0), desired resolution of the computed solution, etc.
#if( TESTNUMBER == 1001 )
		fprintf( fp, "%-21d ! Test number\n", ((int)9) );
#elif( TESTNUMBER == 1002 )
		fprintf( fp, "%-21d ! Test number\n", ((int)10) );
#else
		fprintf( fp, "%-21d ! Test number\n", ((int)0) );
#endif
		fprintf( fp, "%-21.15g ! Left boundary\n", (a - x0) );
		fprintf( fp, "%-21.15g ! Right boundary\n", (b - x0) );
		fprintf( fp, "%-21.15g ! Final time\n", (tf * timescalefactor) );
		fprintf( fp, "%-21d! Number of grid cell points\n", 10 * (int)(ncpux1 * N1 * (DIRX == 1) + ncpux2 * N2 * (DIRX == 2) + ncpux3 * N3 * (DIRX == 3)) );
		fprintf( fp, "%-21.15g ! Accuracy (default = 1e-10)\n", (double) 1e-10 );
		fprintf( fp, "%-21d ! Verbose ( 0 = least info )\n", ((int)0)  );


		fclose( fp );

		//Write the actual Riemann problem set up
		fp = fopen( "RInput.txt", "wt" );

		if( fp == NULL ) {
			dualfprintf( fail_file, "Could not open output file for writing Relativistic Riemann problem input data file, %s\n", "RInput.txt" );
			return;
		}

		//timescalefactor should be zero for Rel. test problems, but we can do this anyways for consistency
		prim_l[U1] /= timescalefactor;
		prim_l[U2] /= timescalefactor;
		prim_l[U3] /= timescalefactor;
		prim_l[UU] /= ( timescalefactor * timescalefactor );
		prim_l[B1] /= timescalefactor; 
		prim_l[B2] /= timescalefactor;
		prim_l[B3] /= timescalefactor;

		prim_r[U1] /= timescalefactor;
		prim_r[U2] /= timescalefactor;
		prim_r[U3] /= timescalefactor;
		prim_r[UU] /= ( timescalefactor * timescalefactor );
		prim_r[B1] /= timescalefactor; 
		prim_r[B2] /= timescalefactor;
		prim_r[B3] /= timescalefactor;

		fprintf( fp, "Bx                   = 0.0D0\n" );
		fprintf( fp, "gamma                = %-15.10E\n", gam );
		fprintf( fp, "rho     (LEFT STATE) = %-15.10E\n", prim_l[RHO] );
		fprintf( fp, "Pgas    (LEFT STATE) = %-15.10E\n", prim_l[UU] * (gam - 1) );
		fprintf( fp, "vx      (LEFT STATE) = %-15.10E\n", prim_l[U1] );
		fprintf( fp, "vy      (LEFT STATE) = %-15.10E\n", prim_l[U2] );
		fprintf( fp, "vz      (LEFT STATE) = %-15.10E\n", prim_l[U3] );
		fprintf( fp, "By      (LEFT STATE) = %-15.10E\n", prim_l[B2] );
		fprintf( fp, "Bz      (LEFT STATE) = %-15.10E\n", prim_l[B3] );
		fprintf( fp, "rho     (RIGHT STATE)= %-15.10E\n", prim_r[RHO] );
		fprintf( fp, "Pgas    (RIGHT STATE)= %-15.10E\n", prim_r[UU] * (gam - 1) );
		fprintf( fp, "vx      (RIGHT STATE)= %-15.10E\n", prim_r[U1] );
		fprintf( fp, "vy      (RIGHT STATE)= %-15.10E\n", prim_r[U2] );
		fprintf( fp, "vz      (RIGHT STATE)= %-15.10E\n", prim_r[U3] );
		fprintf( fp, "By      (RIGHT STATE)= %-15.10E\n", prim_r[B2] );
		fprintf( fp, "Bz      (RIGHT STATE)= %-15.10E\n", prim_r[B3] );

		fclose( fp );
	}
	else if( NONREL1DRIEMANNTEST(TESTNUMBER) == 1 ) {
		fp = fopen( "e1rpex.ini", "wt" );  //weird name Toro riemann solver require the .ini file to have...

		if( fp == NULL ) {
			dualfprintf( fail_file, "Could not open output file for writing Riemann problem parameters\n" );
			return;
		}

		fprintf( fp, "%-15.10g\t! DOMLEN : Domain length, TEST #%d\n", b - a, TESTNUMBER );
		fprintf( fp, "%-15.10g\t! DIAPH1 : Position of diaphragm\n", x0 - a );
		fprintf( fp, "%-15d\t! CELLS  : Number of cells in evaluating exact solution\n", (int) (NXTEST * 5) );  //resolution of the analytic solution should be larger than that of the numerical one
		fprintf( fp, "%-15.10g\t! GAMMA  : Ratio of specific heats\n", gam );
		fprintf( fp, "%-15.10g\t! TIMEOU : Output time\n", tf );
		fprintf( fp, "%-15.10g\t! DL     : Initial density  on left  section of tube\n", prim_l[RHO] );
		fprintf( fp, "%-15.10g\t! UL     : Initial velocity on left  section of tube\t\n", prim_l[U1] );
		fprintf( fp, "%-15.10g\t! PL     : Initial pressure on left  section of tube\n", prim_l[UU] * (gam - 1) );
		fprintf( fp, "%-15.10g\t! DR     : Initial density  on right section of tube\n", prim_r[RHO] );
		fprintf( fp, "%-15.10g\t! UR     : Initial velocity on right section of tube\n", prim_r[U1] );
		fprintf( fp, "%-15.10g\t! PR     : Initial pressure on right section of tube\n", prim_r[UU] * (gam - 1) );
		fprintf( fp, "%-15.10g\t! PSCALE : Normalising factor for pressure and energy\n", (double)1.0 );

		fclose( fp );

	}
	else {
		dualfprintf( fail_file, "Unknown test number for a Riemann problem: %d\n", (int) TESTNUMBER );
	}

	return;
}




// assumes normal field in pr
int init_vpot(int l, int i, int j, int k, FTYPE *A)
{

	//no magnetic field -- atch
  *A=0;
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
	//do not normalize densities

  return(0);
}


// assumes normal field definition
int normalize_field(FTYPE p[][N2M][N3M][NPR])
{
  //do not normalize field -- atch
  return(0);
}

FTYPE interpn( int order, FTYPE x_eval,  FTYPE x1, FTYPE f1, FTYPE x2, FTYPE f2, FTYPE x3, FTYPE f3, FTYPE x4, FTYPE f4, FTYPE x5, FTYPE f5, FTYPE x6, FTYPE f6 ) 
//Arguments:
//x_eval -- abscissa of a point where an interpolated value is to be evaluated
//x1, x2, x3 -- set of abscissas; should all be different; can come in any order
//f1, f2, f3 -- set of function values at the above abscissas, i.e. f# = f( x# ), # = 1..3

//Description:
//Constructs a second-order polynomial interpolating the set of points {x#, f#}, #=1..3 and
//evaluates it at the point x.
//Interpolation is performed using the standard Lagrange method.
{
	int i, j; //iterator through the abscissa points
#define max_interp_order  6 //should be equal to the number of points input to the function; degree of interpolation is only limited by the number of arguments to the function --
	//so increase that number if require a higher order; also, modify the definitions of x[] and f[] to include the added arguments

	FTYPE f_eval = 0.0; //interpolated value to be returned
	FTYPE x[] = { x1, x2, x3, x4, x5, x6};
	FTYPE f[] = { f1, f2, f3, f4, f5, f6};
	FTYPE x_a, x_b; //temporary variables
	FTYPE delta_f_eval;

	if( order < 1 || order > max_interp_order) {
		trifprintf( "interpn: requested order of interpolation non-positive or too high.\n" );
		myexit( 1 );
	}

	//test if the supplied x-values are valid: there should be no identical values
#if(DO_ASSERTS)
	for( i = 1; i < order; i++ ) {
		for( j = 0; j < i; j++ ) {
			if( x[i] - x[j] == 0.0 ) {
				fprintf( stderr, "interpn: abscissas of values to be interpolated are to be all different\n" );
				return 0.;
			}
		}
	}
#endif

	for( i = 0; i < order; i++ ) {
		//construct a polynomial that would be equal to f[i] at x[i] and equal to 0 at other x-points ( x[(i + j) % order], j = 1..order )
		x_a = x[i];  
		delta_f_eval = f[i];

		for( j = 1; j < order; j++ ) {
			x_b = x[(i + j) % order];
			delta_f_eval *= (x_eval - x_b) / (x_a - x_b);
		}

		//add up the value of this polynomial at x = x_eval to the final result
		f_eval += delta_f_eval;

	}

	return f_eval;
}

