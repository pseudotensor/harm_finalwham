#include "global.h"
#include "mpidefs.h"


/*

to add a new variable:

0) global.h: setup macro to turn on/off memory (a_name) below.

1) Add a variable called a_name that is the *definition* of the memory space

2) Lower down in this file, define a pointer of the same type

3) set_arrays.c: add pointer shifting

4) set_arrays.c: down lower, assign to 0 or valueinit

5) set_array.c: use global.h macro to avoid assignments to memory is turned off

6) Use it!

*/



FTYPE a_p[N1M][N2M][N3M][NPR];	/* space for primitive vars */

FTYPE a_panalytic[N1M][N2M][N3M][NPR];       /* space for primitive vars */
FTYPE a_omegafanalytic[N1M][N2M][N3M][NPR];


// emf has extra zone on upper end since corner quantity and some components exist that are needed for cell centered quantities
FTYPE a_emf[COMPDIM][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3];	/* space for emf */
FTYPE a_vconemf[N1M][N2M][N3M][NDIM-1];	/* used for Athena EMFs */

FTYPE a_pleft[N1M][N2M][N3M][NPR2INTERP]; /* for parabolic interpolation */
FTYPE a_pright[N1M][N2M][N3M][NPR2INTERP]; /* for parabolic interpolation */

// for higher order RK time stepping integrations
FTYPE a_ulast[N1M][N2M][N3M][NPR]; 
FTYPE a_unew[N1M][N2M][N3M][NPR];

#if(DOENOFLUXMEMORY)
FTYPE a_uinitial[N1M][N2M][N3M][NPR];
FTYPE a_upoint[N1M][N2M][N3M][NPR];
FTYPE a_dUgeomarray[N1M][N2M][N3M][NPR];
#endif


#if(STOREWAVESPEEDS)
FTYPE a_wspeed[COMPDIM][2][N1M][N2M][N3M]; // wave speeds (left/right)
#endif

#if(DODQMEMORY)
#if(N1>1)
FTYPE a_dq1[N1M][N2M][N3M][NPR2INTERP];	/* slopes */
#endif
#if(N2>1)
FTYPE a_dq2[N1M][N2M][N3M][NPR2INTERP];	/* slopes */
#endif
#if(N3>1)
FTYPE a_dq3[N1M][N2M][N3M][NPR2INTERP];	/* slopes */
#endif
#endif

#if(N1>1)
FTYPE a_F1[N1M][N2M][N3M][NPR];	/* fluxes */
#endif
#if(N2>1)
FTYPE a_F2[N1M][N2M][N3M][NPR];	/* fluxes */
#endif
#if(N3>1)
FTYPE a_F3[N1M][N2M][N3M][NPR];	/* fluxes */
#endif

#if( LIMIT_FLUXC2A_PRIM_CHANGE )
FTYPE a_fluxvectemp[N1M][N2M][N3M][NPR];	/* termporary storage for flux */ //atch
#endif

// unsplit 3D flux_interp storage or memory for a2c unsplit method
#if( (FLUXDIMENSPLIT==UNSPLIT)&&(N1NOT1*N2NOT1*N3NOT1) || ((A2CDIMENSPLIT==UNSPLIT)&&( (N1NOT1*N2NOT1*N3NOT1)||(N1NOT1&&N2NOT1)||(N1NOT1&&N3NOT1)||(N2NOT1&&N3NOT1) ) ) )
FTYPE a_Fa[N1M][N2M][N3M][NPR];	/* fluxes */
FTYPE a_Fb[N1M][N2M][N3M][NPR];	/* fluxes */
#endif

#if(N1>1)
//FTYPE a_F1CT[N1M][N2M][N3M][NPR];	/* fluxes */
#endif
#if(N2>1)
//FTYPE a_F2CT[N1M][N2M][N3M][NPR];	/* fluxes */
#endif
#if(N3>1)
//FTYPE a_F3CT[N1M][N2M][N3M][NPR];	/* fluxes */
#endif

FTYPE a_pk[MAXDTSTAGES][N1M][N2M][N3M][NPR];	/* next-step primitives */

FTYPE a_prc[N1M][N2M][N3M][NPR2INTERP];	/* rescaled primitives, also used for temporary storage in fixup_utoprim() */

// arbitrary temporary storage array
FTYPE a_ptemparray[N1M][N2M][N3M][NPR];

#if(CALCFARADAYANDCURRENTS)
// current density stuff
FTYPE a_cfaraday[N1M][N2M][N3M][NUMCURRENTSLOTS][3];
FTYPE a_fcon[N1M][N2M][N3M][NUMFARADAY];
FTYPE a_jcon[N1M][N2M][N3M][NDIM];
#endif

/////////////////
//
// AVERAGES
//
////////////////
#if(DOAVG)
// time average stuff
FTYPE a_normalvarstavg[N1M][N2M][N3M][NUMNORMDUMP];
FTYPE a_anormalvarstavg[N1M][N2M][N3M][NUMNORMDUMP];

#if(CALCFARADAYANDCURRENTS)
FTYPE a_jcontavg[N1M][N2M][N3M][NDIM];
FTYPE a_jcovtavg[N1M][N2M][N3M][NDIM];
FTYPE a_ajcontavg[N1M][N2M][N3M][NDIM];
FTYPE a_ajcovtavg[N1M][N2M][N3M][NDIM];
#endif

FTYPE a_massfluxtavg[N1M][N2M][N3M][NDIM];
FTYPE a_amassfluxtavg[N1M][N2M][N3M][NDIM];

FTYPE a_othertavg[N1M][N2M][N3M][NUMOTHER];
FTYPE a_aothertavg[N1M][N2M][N3M][NUMOTHER];

#if(CALCFARADAYANDCURRENTS)
FTYPE a_fcontavg[N1M][N2M][N3M][NUMFARADAY];
FTYPE a_fcovtavg[N1M][N2M][N3M][NUMFARADAY];
FTYPE a_afcontavg[N1M][N2M][N3M][NUMFARADAY];
FTYPE a_afcovtavg[N1M][N2M][N3M][NUMFARADAY];
#endif

FTYPE a_tudtavg[N1M][N2M][N3M][NUMSTRESSTERMS];
FTYPE a_atudtavg[N1M][N2M][N3M][NUMSTRESSTERMS];
#endif


///////////////////////////
//
// COUNTERS and failure checks
//
////////////////////////////

int a_pflag[N1M][N2M][N3M][NUMPFLAGS];	/* space for flag */


#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
	FTYPE a_weno_lower_order_fraction[N1M][N2M][N3M][NPR];	/* space for lower order indicators */  //atch
	FTYPE do_weno_lower_order_fraction;
#endif 

#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )
	FTYPE a_weno_prim_lower_order_fraction[NDIM][N1M][N2M][N3M];	/* space for lower order fraction of primitive quantities */  //atch
#endif 


#if(DOENODEBUG)
// 3: 1,2,3 directions
// 5: c2e on P, a2c on U, c2a for U, c2a for Fperp1, c2a for Fperp2
CTYPE a_enodebugarray[N1M][N2M][N3M][3-1][NUMINTERPTYPES][NPR-4][NUMENODEBUGS];// space for debugging eno
#endif

/* for debug */
#if(DODEBUG)
CTYPE a_failfloorcount[N1M][N2M][N3M][NUMTSCALES][NUMFAILFLOORFLAGS]; // number of failures and floor adjustments for each zone
#endif


#if(DODISS)
FTYPE a_dissfunpos[N1M][N2M][N3M][NUMDISSFUNPOS]; // storage for dissipation function
#endif

///////////////////////////////
//
// grid functions (+1 size larger so can have geometry at upper corners -- say for vector potential or whatever)
//
///////////////////////////////
// dummy global variables
FTYPE gengcov[NDIM][NDIM]; // gengcov/gengcon used as pointers for get_geometry() in gset() in metric.c .  This allows native coords to be different than requested coords.
FTYPE gengcon[NDIM][NDIM];

FTYPE gengcovpert[NDIM];

#if(MCOORD!=CARTMINKMETRIC)

// grid functions that exist at multiple locations within a cell
FTYPE a_gcon[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
FTYPE a_gcov[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
FTYPE a_gcovpert[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM];
FTYPE a_gdet[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG];
FTYPE a_eomfunc[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG];
//FTYPE a_dxdxp[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];

// rest of grid functions always at center
FTYPE a_conn[N1M][N2M][N3M][NDIM][NDIM][NDIM];
FTYPE a_conn2[N1M][N2M][N3M][NDIM];
#if(VOLUMEDIFF)
FTYPE a_idxvol[N1M][N2M][N3M][NDIM]; // volume regularization 1/dx for each direction
#endif


#else

// grid functions that exist at multiple locations within a cell
FTYPE a_gcon[1][1][1][NPG][NDIM][NDIM];
FTYPE a_gcov[1][1][1][NPG][NDIM][NDIM];
FTYPE a_gcovpert[1][1][1][NPG][NDIM];
FTYPE a_gdet[1][1][1][NPG];
FTYPE a_eomfunc[1][1][1][NPG];
//FTYPE a_dxdxp[1][1][1][NPG][NDIM][NDIM];

// rest of grid functions always at center
FTYPE a_conn[1][1][1][NDIM][NDIM][NDIM];
FTYPE a_conn2[1][1][1][NDIM];

#if(VOLUMEDIFF)
FTYPE a_idxvol[1][1][1][NDIM]; // volume regularization 1/dx for each direction
#endif


#endif









/////////////////////////////////////////////////////////////////
//
//
// (possibly shifted) pointers to memory  real space
//
//
/////////////////////////////////////////////////////////////////




int (*pflag)[N2M][N3M][NUMPFLAGS];
CTYPE (*enodebugarray)[N2M][N3M][3-1][NUMINTERPTYPES][NPR-4][NUMENODEBUGS];
CTYPE (*failfloorcount)[N2M][N3M][NUMTSCALES][NUMFAILFLOORFLAGS];

FTYPE (*dissfunpos)[N2M][N3M][NUMDISSFUNPOS];


FTYPE (*p)[N2M][N3M][NPR];
FTYPE (*panalytic)[N2M][N3M][NPR];
FTYPE (*omegafanalytic)[N2M][N3M][NPR];


FTYPE (*pdump)[N2M][N3M][NPR];


FTYPE (*emf)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3];
FTYPE (*vconemf)[N2M][N3M][NDIM-1];

FTYPE (*wspeed)[2][N1M][N2M][N3M]; // wave speeds

FTYPE (*ubound)[N2M][N3M][NPR];
FTYPE (*udump)[N2M][N3M][NPR];

FTYPE (*unew)[N2M][N3M][NPR];
FTYPE (*ulast)[N2M][N3M][NPR];
FTYPE (*uinitial)[N2M][N3M][NPR];
FTYPE (*upoint)[N2M][N3M][NPR];
FTYPE (*dUgeomarray)[N2M][N3M][NPR];



FTYPE (*pleft)[N2M][N3M][NPR2INTERP];
FTYPE (*pright)[N2M][N3M][NPR2INTERP];

FTYPE (*dq1)[N2M][N3M][NPR2INTERP];
FTYPE (*dq2)[N2M][N3M][NPR2INTERP];
FTYPE (*dq3)[N2M][N3M][NPR2INTERP];

FTYPE (*prc)[N2M][N3M][NPR2INTERP];

FTYPE (*ptemparray)[N2M][N3M][NPR];


FTYPE (*F1)[N2M][N3M][NPR];
FTYPE (*F2)[N2M][N3M][NPR];
FTYPE (*F3)[N2M][N3M][NPR];

#if( LIMIT_FLUXC2A_PRIM_CHANGE )
FTYPE (*fluxvectemp)[N2M][N3M][NPR];	/* termporary storage for flux */ //atch
#endif

FTYPE (*Fa)[N2M][N3M][NPR];
FTYPE (*Fb)[N2M][N3M][NPR];

FTYPE (*a2cin)[N2M][N3M][NPR];
FTYPE (*a2cout)[N2M][N3M][NPR];

//FTYPE (*F1CT)[N2M][N3M][NPR];
//FTYPE (*F2CT)[N2M][N3M][NPR];
//FTYPE (*F3CT)[N2M][N3M][NPR];
FTYPE (*pk)[N1M][N2M][N3M][NPR];





///////////////////////////////
//
// grid functions (+1 size larger so can have geometry at upper corners -- say for vector potential or whatever)
//
///////////////////////////////
FTYPE (*gcon)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
FTYPE (*gcov)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
FTYPE (*gcovpert)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM];
FTYPE (*gdet)[N2M+SHIFT2][N3M+SHIFT3][NPG];
FTYPE (*eomfunc)[N2M+SHIFT2][N3M+SHIFT3][NPG];
//FTYPE (*dxdxp)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM];

FTYPE (*conn)[N2M][N3M][NDIM][NDIM][NDIM];
FTYPE (*conn2)[N2M][N3M][NDIM];
FTYPE (*idxvol)[N2M][N3M][NDIM];





// current density stuff
FTYPE (*cfaraday)[N2M][N3M][NUMCURRENTSLOTS][3];
FTYPE (*fcon)[N2M][N3M][NUMFARADAY];
FTYPE (*jcon)[N2M][N3M][NDIM];

// time average stuff
FTYPE (*normalvarstavg)[N2M][N3M][NUMNORMDUMP];
FTYPE (*anormalvarstavg)[N2M][N3M][NUMNORMDUMP];

FTYPE (*jcontavg)[N2M][N3M][NDIM];
FTYPE (*jcovtavg)[N2M][N3M][NDIM];
FTYPE (*ajcontavg)[N2M][N3M][NDIM];
FTYPE (*ajcovtavg)[N2M][N3M][NDIM];

FTYPE (*massfluxtavg)[N2M][N3M][NDIM];
FTYPE (*amassfluxtavg)[N2M][N3M][NDIM];

FTYPE (*othertavg)[N2M][N3M][NUMOTHER];
FTYPE (*aothertavg)[N2M][N3M][NUMOTHER];

FTYPE (*fcontavg)[N2M][N3M][NUMFARADAY];
FTYPE (*fcovtavg)[N2M][N3M][NUMFARADAY];
FTYPE (*afcontavg)[N2M][N3M][NUMFARADAY];
FTYPE (*afcovtavg)[N2M][N3M][NUMFARADAY];

FTYPE (*tudtavg)[N2M][N3M][NUMSTRESSTERMS];
FTYPE (*atudtavg)[N2M][N3M][NUMSTRESSTERMS];

// for image, see image.c
FTYPE (*pimage)[N2M][N3M][NPR];


//#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
	FTYPE (*weno_lower_order_fraction)[N2M][N3M][NPR]; //atch
//#endif

//#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )
	FTYPE (*weno_prim_lower_order_fraction)[N1M][N2M][N3M]; //atch
//#endif


//#if(DOENO)
FTYPE (*uenotmp0)[N2M][N3M][NPR]; /* for ENO reconstruction of U or dU*/
FTYPE (*uenotmp1)[N2M][N3M][NPR]; /* for ENO reconstruction of U or dU*/
FTYPE (*uenotmp2)[N2M][N3M][NPR]; /* for ENO reconstruction of U or dU*/
//#endif










/////////////////////////////////////////////////////////////////
//
//
// GLOBAL VARIABLES that are not for every point in space
//
//
/////////////////////////////////////////////////////////////////


SFTYPE *lumvsr,*lumvsr_tot;

SFTYPE *dissvsr,*dissvsr_tot;


/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
FTYPE a;
FTYPE gam;
FTYPE Bpole,Omegastar; 

/* numerical parameters */
int defcoord;
FTYPE Rin, R0, Rout, hslope, Zin, Zout;
FTYPE Rin_array[NDIM], Rout_array[NDIM];  //atch -- arrays for a more general way of handling the grid dimensions
FTYPE Risco,Rhor;
FTYPE cour;
FTYPE dV, dVF, dx[NDIM], startx[NDIM];
SFTYPE dt,t,tf,tstepparti,tsteppartf;
FTYPE rcurr, hcurr;

//int istart, istop, jstart, jstop;
#if(SIMULBCCALC!=-1) 
	int isc,iec,jsc,jec;
	int isf1,ief1,jsf1,jef1,ksf1,kef1;
	int isf2,ief2,jsf2,jef2,ksf2,kef2;
	int isf3,ief3,jsf3,jef3,ksf3,kef3;
	int ise,iee,jse,jee;
	int isf1ct,ief1ct,jsf1ct,jef1ct;// GODMARK: other stage type requires more
	int isf2ct,ief2ct,jsf2ct,jef2ct;
	int isf3ct,ief3ct,jsf3ct,jef3ct;
	int isdq,iedq,jsdq,jedq;
	int ispdq,iepdq,jspdq,jepdq;
#endif

FTYPE mydminarg1, mydminarg2;
long nstep;

int steppart,numstepparts;


/* output parameters */
SFTYPE DTd;
SFTYPE DTavg;
SFTYPE DTener;
SFTYPE DTi;
SFTYPE DTdebug;
long DTr;
long dump_cnt;
long avg_cnt;
long debug_cnt;
long image_cnt;
long rdump_cnt;
long fieldline_cnt; // assumed to keep track with images (as in diag.c), so no need to include in restart()
int nstroke;

/* global flags */
int failed;
int lim,fluxmethod,FLUXB,UTOPRIMVERSION,TIMEORDER,DOENOFLUX,avgscheme,do_transverse_flux_integration,do_conserved_integration,do_source_integration;
FTYPE defcon;

/* diagnostics */
// don't track this separately in other regions except global region
SFTYPE frdot[N1][NPR];
SFTYPE pdottermsjet2[COMPDIM*2][NUMFLUXTERMS][NPR];
CTYPE failfloorcountlocal[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet
CTYPE failfloorcountlocal_tot[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet

// general stuff for ener.out file for regions to completely track, including terms within flux
int dofluxreg[NUMENERREGIONS][COMPDIM*2];
int enerposreg[NUMENERREGIONS][COMPDIM*2];
// these quantities contain diagnostics
// all these require writing to restart file
// other _tot quantities appear in dump_ener.c that don't need to be written to restart file since easily computed from existing data.
SFTYPE fladdreg[NUMENERREGIONS][NPR];
SFTYPE fladdreg_tot[NUMENERREGIONS][NPR];
SFTYPE fladdtermsreg[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
SFTYPE fladdtermsreg_tot[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
SFTYPE Ureg_init[NUMENERREGIONS][NPR];
SFTYPE Ureg_init_tot[NUMENERREGIONS][NPR];
SFTYPE pcumreg[NUMENERREGIONS][COMPDIM*2][NPR];
SFTYPE pcumreg_tot[NUMENERREGIONS][COMPDIM*2][NPR];
SFTYPE pdotreg[NUMENERREGIONS][COMPDIM*2][NPR];
SFTYPE pdottermsreg[NUMENERREGIONS][COMPDIM*2][NUMFLUXTERMS][NPR];
SFTYPE sourceaddreg[NUMENERREGIONS][NPR];
SFTYPE sourceaddreg_tot[NUMENERREGIONS][NPR];
SFTYPE sourceaddtermsreg[NUMENERREGIONS][NUMSOURCES][NPR];
SFTYPE sourceaddtermsreg_tot[NUMENERREGIONS][NUMSOURCES][NPR];
SFTYPE dissreg[NUMENERREGIONS][1];
SFTYPE dissreg_tot[NUMENERREGIONS][1];

// used for each region, related to global quantities
// _tot quantities here are global since used in restart.
int *doflux;
int *enerpos;
SFTYPE *fladd;
SFTYPE *fladd_tot;
SFTYPE (*fladdterms)[NPR];
SFTYPE (*fladdterms_tot)[NPR];
SFTYPE *U_init;
SFTYPE *U_init_tot;
SFTYPE (*pcum)[NPR];
SFTYPE (*pcum_tot)[NPR];
SFTYPE (*pdot)[NPR];
SFTYPE (*pdotterms)[NUMFLUXTERMS][NPR];
SFTYPE *sourceadd;
SFTYPE *sourceadd_tot;
SFTYPE (*sourceaddterms)[NPR];
SFTYPE (*sourceaddterms_tot)[NPR];
SFTYPE *diss;
SFTYPE *diss_tot;

// end changes after ...

/* current local position */
int icurr, jcurr, kcurr, pcurr, ihere, jhere, phere;

/* Jon's addition */
int horizoni;
long realnstep;
int partialstep;
int mpicombine;
int mpicombinetype;
int truempicombinetype;
int halftimep;
int whichrestart;
int appendold;
int whocalleducon;
// global flags
long restartsteps[2];
int binaryoutput,sortedoutput;
long steptofaildump,steptofailmap;
int ifail,jfail,kfail,dofailmap,dofaildump,restartonfail;
// IC
FTYPE h_over_r;
// BC
FTYPE h_over_r_jet;
int BCtype[COMPDIM*2];
int rescaletype;
int cooling;
int DOENERDIAG,DOGDUMPDIAG,DORDUMPDIAG,DODUMPDIAG,DOAVGDIAG, DOIMAGEDIAG,DOAREAMAPDIAG;
int GAMMIEDUMP,GAMMIEIMAGE,GAMMIEENER,DODIAGS,RESTARTMODE,WHICHFILE,POSDEFMETRIC,DOENODEBUGEVERYSUBSTEP,DODIAGEVERYSUBSTEP,
		INVERTFROMAVERAGEIFFAILED,LIMIT_AC_PRIM_FRAC_CHANGE; //atch
FTYPE RHOMIN,UUMIN,RHOMINLIMIT,UUMINLIMIT,MAX_AC_PRIM_FRAC_CHANGE; //atch
FTYPE prMAX[NPR];
FTYPE BSQORHOLIMIT,BSQOULIMIT,GAMMAMAX,GAMMADAMP,GAMMAFAIL;
FTYPE SAFE;
int debugfail;
FTYPE uttdiscr; // for check_pr for now
int jonchecks;
int dnumcolumns[NUMDUMPTYPES];
struct blink * blinkptr0[NUMDUMPTYPES];
struct blink * cpulinkptr0[NUMDUMPTYPES];
int DOCOLSPLIT[NUMDUMPTYPES];
int docolsplit; // global var for now
int nextcol;

/* physical consts */
FTYPE msun,lsun,G,H,C,mn,me,kb,arad,sigmasb,sigmamat,mevocsq,ergPmev;
FTYPE mb,M,Mdot,Mdotc;
FTYPE Lunit,Tunit,rhounit,Munit,mdotunit,energyunit,edotunit,Pressureunit,Tempunit,Bunit,massunitPmsun;
FTYPE ledd,leddcode;

int NUMBUFFERS;

int ijkminmax[NUMINTERPTYPES][NDIM][2]; // 0=in/1=out (or similar)
// as for BCtype
int ijkminmaxud[NUMINTERPTYPES][COMPDIM*2];
int fluxloop[NDIM][NUMFLUXLOOPNUMBERS];
int Uconsloop[NUMFLUXLOOPNUMBERS];
int interporder[NUMINTERPS];


// ENO DEBUG GLOBAL VARS
int dirglobal,iglobal,jglobal,kglobal,iterglobal,interporfluxglobal;

FTYPE prfloorcoef[NPR];

int numbercpu[ 3+1 ];

FTYPE globalinv[NUMGLOBALINV];
char globalinvtext[NUMGLOBALINV][10];

// Ramesh stuff
FTYPE nu,ss,ucrit,Ttpow,jetalpha;

int lntries;
FTYPE lerrx;

// EOS related functions
FTYPE (*ptr_pressure_rho0_u)(FTYPE rho0, FTYPE u);
FTYPE (*ptr_compute_u_from_entropy)(FTYPE rho0, FTYPE entropy);
FTYPE (*ptr_u_rho0_p)(FTYPE rho0, FTYPE p);
FTYPE (*ptr_dpdu_rho0_u)(FTYPE rho0, FTYPE u);
FTYPE (*ptr_dpdrho0_rho0_u)(FTYPE rho0, FTYPE u);
FTYPE (*ptr_cs2_compute)(FTYPE rho0, FTYPE u);
FTYPE (*ptr_compute_dSdrho)(FTYPE rho0, FTYPE u);
FTYPE (*ptr_compute_dSdu)(FTYPE rho0, FTYPE u);
FTYPE (*ptr_compute_entropy)(FTYPE rho0, FTYPE u);
FTYPE (*ptr_pressure_wmrho0) (FTYPE rho0, FTYPE wmrho0);
FTYPE (*ptr_compute_idwmrho0dp) (FTYPE rho0, FTYPE wmrho0);
FTYPE (*ptr_compute_idrho0dp) (FTYPE rho0, FTYPE wmrho0);
