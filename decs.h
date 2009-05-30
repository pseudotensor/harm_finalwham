#include "global.h"
#include "mpidecs.h"


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



extern FTYPE a_p[N1M][N2M][N3M][NPR];	/* space for primitive vars */

extern FTYPE a_panalytic[N1M][N2M][N3M][NPR];       /* space for primitive vars */
extern FTYPE a_omegafanalytic[N1M][N2M][N3M][NPR];


// emf has extra zone on upper end since corner quantity and some components exist that are needed for cell centered quantities
extern FTYPE a_emf[COMPDIM][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3];	/* space for emf */
extern FTYPE a_vconemf[N1M][N2M][N3M][NDIM-1];	/* used for Athena EMFs */

extern FTYPE a_pleft[N1M][N2M][N3M][NPR2INTERP]; /* for parabolic interpolation */
extern FTYPE a_pright[N1M][N2M][N3M][NPR2INTERP]; /* for parabolic interpolation */

// for higher order RK time stepping integrations
extern FTYPE a_ulast[N1M][N2M][N3M][NPR]; 
extern FTYPE a_unew[N1M][N2M][N3M][NPR];

#if(DOENOFLUXMEMORY)
extern FTYPE a_uinitial[N1M][N2M][N3M][NPR];
extern FTYPE a_upoint[N1M][N2M][N3M][NPR];
extern FTYPE a_dUgeomarray[N1M][N2M][N3M][NPR];
#endif


#if(STOREWAVESPEEDS)
extern FTYPE a_wspeed[COMPDIM][2][N1M][N2M][N3M]; // wave speeds (left/right)
#endif

#if(DODQMEMORY)
#if(N1>1)
extern FTYPE a_dq1[N1M][N2M][N3M][NPR2INTERP];	/* slopes */
#endif
#if(N2>1)
extern FTYPE a_dq2[N1M][N2M][N3M][NPR2INTERP];	/* slopes */
#endif
#if(N3>1)
extern FTYPE a_dq3[N1M][N2M][N3M][NPR2INTERP];	/* slopes */
#endif
#endif

#if(N1>1)
extern FTYPE a_F1[N1M][N2M][N3M][NPR];	/* fluxes */
#endif
#if(N2>1)
extern FTYPE a_F2[N1M][N2M][N3M][NPR];	/* fluxes */
#endif
#if(N3>1)
extern FTYPE a_F3[N1M][N2M][N3M][NPR];	/* fluxes */
#endif

#if( LIMIT_FLUXC2A_PRIM_CHANGE )
extern FTYPE a_fluxvectemp[N1M][N2M][N3M][NPR];	/* termporary storage for flux */ //atch
#endif

// unsplit 3D flux_interp storage or memory for a2c unsplit method
#if( (FLUXDIMENSPLIT==UNSPLIT)&&(N1NOT1*N2NOT1*N3NOT1) || ((A2CDIMENSPLIT==UNSPLIT)&&( (N1NOT1*N2NOT1*N3NOT1)||(N1NOT1&&N2NOT1)||(N1NOT1&&N3NOT1)||(N2NOT1&&N3NOT1) ) ) )
extern FTYPE a_Fa[N1M][N2M][N3M][NPR];	/* fluxes */
extern FTYPE a_Fb[N1M][N2M][N3M][NPR];	/* fluxes */
#endif

#if(N1>1)
//extern FTYPE a_F1CT[N1M][N2M][N3M][NPR];	/* fluxes */
#endif
#if(N2>1)
//extern FTYPE a_F2CT[N1M][N2M][N3M][NPR];	/* fluxes */
#endif
#if(N3>1)
//extern FTYPE a_F3CT[N1M][N2M][N3M][NPR];	/* fluxes */
#endif

extern FTYPE a_pk[MAXDTSTAGES][N1M][N2M][N3M][NPR];	/* next-step primitives */

extern FTYPE a_prc[N1M][N2M][N3M][NPR2INTERP];	/* rescaled primitives, also used for temporary storage in fixup_utoprim() */

// arbitrary temporary storage array
extern FTYPE a_ptemparray[N1M][N2M][N3M][NPR];

#if(CALCFARADAYANDCURRENTS)
// current density stuff
extern FTYPE a_cfaraday[N1M][N2M][N3M][NUMCURRENTSLOTS][3];
extern FTYPE a_fcon[N1M][N2M][N3M][NUMFARADAY];
extern FTYPE a_jcon[N1M][N2M][N3M][NDIM];
#endif

/////////////////
//
// AVERAGES
//
////////////////
#if(DOAVG)
// time average stuff
extern FTYPE a_normalvarstavg[N1M][N2M][N3M][NUMNORMDUMP];
extern FTYPE a_anormalvarstavg[N1M][N2M][N3M][NUMNORMDUMP];

#if(CALCFARADAYANDCURRENTS)
extern FTYPE a_jcontavg[N1M][N2M][N3M][NDIM];
extern FTYPE a_jcovtavg[N1M][N2M][N3M][NDIM];
extern FTYPE a_ajcontavg[N1M][N2M][N3M][NDIM];
extern FTYPE a_ajcovtavg[N1M][N2M][N3M][NDIM];
#endif

extern FTYPE a_massfluxtavg[N1M][N2M][N3M][NDIM];
extern FTYPE a_amassfluxtavg[N1M][N2M][N3M][NDIM];

extern FTYPE a_othertavg[N1M][N2M][N3M][NUMOTHER];
extern FTYPE a_aothertavg[N1M][N2M][N3M][NUMOTHER];

#if(CALCFARADAYANDCURRENTS)
extern FTYPE a_fcontavg[N1M][N2M][N3M][NUMFARADAY];
extern FTYPE a_fcovtavg[N1M][N2M][N3M][NUMFARADAY];
extern FTYPE a_afcontavg[N1M][N2M][N3M][NUMFARADAY];
extern FTYPE a_afcovtavg[N1M][N2M][N3M][NUMFARADAY];
#endif

extern FTYPE a_tudtavg[N1M][N2M][N3M][NUMSTRESSTERMS];
extern FTYPE a_atudtavg[N1M][N2M][N3M][NUMSTRESSTERMS];
#endif


///////////////////////////
//
// COUNTERS and failure checks
//
////////////////////////////

extern int a_pflag[N1M][N2M][N3M][NUMPFLAGS];	/* space for flag */


#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
	extern FTYPE a_weno_lower_order_fraction[N1M][N2M][N3M][NPR];	/* space for lower order indicators */  //atch
	extern FTYPE do_weno_lower_order_fraction;
#endif 

#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )
	extern FTYPE a_weno_prim_lower_order_fraction[NDIM][N1M][N2M][N3M];	/* space for lower order fraction of primitive quantities */  //atch
#endif 


#if(DOENODEBUG)
// 3: 1,2,3 directions
// 5: c2e on P, a2c on U, c2a for U, c2a for Fperp1, c2a for Fperp2
extern CTYPE a_enodebugarray[N1M][N2M][N3M][3-1][NUMINTERPTYPES][NPR-4][NUMENODEBUGS];// space for debugging eno
#endif

/* for debug */
#if(DODEBUG)
extern CTYPE a_failfloorcount[N1M][N2M][N3M][NUMTSCALES][NUMFAILFLOORFLAGS]; // number of failures and floor adjustments for each zone
#endif


#if(DODISS)
extern FTYPE a_dissfunpos[N1M][N2M][N3M][NUMDISSFUNPOS]; // storage for dissipation function
#endif

///////////////////////////////
//
// grid functions (+1 size larger so can have geometry at upper corners -- say for vector potential or whatever)
//
///////////////////////////////
// dummy global variables
extern FTYPE gengcov[NDIM][NDIM]; // gengcov/gengcon used as pointers for get_geometry() in gset() in metric.c .  This allows native coords to be different than requested coords.
extern FTYPE gengcon[NDIM][NDIM];

extern FTYPE gengcovpert[NDIM];

#if(MCOORD!=CARTMINKMETRIC)

// grid functions that exist at multiple locations within a cell
extern FTYPE a_gcon[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
extern FTYPE a_gcov[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
extern FTYPE a_gcovpert[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM];
extern FTYPE a_gdet[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG];
extern FTYPE a_eomfunc[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG];
//extern FTYPE a_dxdxp[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];

// rest of grid functions always at center
extern FTYPE a_conn[N1M][N2M][N3M][NDIM][NDIM][NDIM];
extern FTYPE a_conn2[N1M][N2M][N3M][NDIM];
#if(VOLUMEDIFF)
extern FTYPE a_idxvol[N1M][N2M][N3M][NDIM]; // volume regularization 1/dx for each direction
#endif


#else

// grid functions that exist at multiple locations within a cell
extern FTYPE a_gcon[1][1][1][NPG][NDIM][NDIM];
extern FTYPE a_gcov[1][1][1][NPG][NDIM][NDIM];
extern FTYPE a_gcovpert[1][1][1][NPG][NDIM];
extern FTYPE a_gdet[1][1][1][NPG];
extern FTYPE a_eomfunc[1][1][1][NPG];
//extern FTYPE a_dxdxp[1][1][1][NPG][NDIM][NDIM];

// rest of grid functions always at center
extern FTYPE a_conn[1][1][1][NDIM][NDIM][NDIM];
extern FTYPE a_conn2[1][1][1][NDIM];

#if(VOLUMEDIFF)
extern FTYPE a_idxvol[1][1][1][NDIM]; // volume regularization 1/dx for each direction
#endif


#endif









/////////////////////////////////////////////////////////////////
//
//
// (possibly shifted) pointers to memory  real space
//
//
/////////////////////////////////////////////////////////////////




extern int (*pflag)[N2M][N3M][NUMPFLAGS];
extern CTYPE (*enodebugarray)[N2M][N3M][3-1][NUMINTERPTYPES][NPR-4][NUMENODEBUGS];
extern CTYPE (*failfloorcount)[N2M][N3M][NUMTSCALES][NUMFAILFLOORFLAGS];

extern FTYPE (*dissfunpos)[N2M][N3M][NUMDISSFUNPOS];


extern FTYPE (*p)[N2M][N3M][NPR];
extern FTYPE (*panalytic)[N2M][N3M][NPR];
extern FTYPE (*omegafanalytic)[N2M][N3M][NPR];


extern FTYPE (*pdump)[N2M][N3M][NPR];


extern FTYPE (*emf)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3];
extern FTYPE (*vconemf)[N2M][N3M][NDIM-1];

extern FTYPE (*wspeed)[2][N1M][N2M][N3M]; // wave speeds

extern FTYPE (*ubound)[N2M][N3M][NPR];
extern FTYPE (*udump)[N2M][N3M][NPR];

extern FTYPE (*unew)[N2M][N3M][NPR];
extern FTYPE (*ulast)[N2M][N3M][NPR];
extern FTYPE (*uinitial)[N2M][N3M][NPR];
extern FTYPE (*upoint)[N2M][N3M][NPR];
extern FTYPE (*dUgeomarray)[N2M][N3M][NPR];



extern FTYPE (*pleft)[N2M][N3M][NPR2INTERP];
extern FTYPE (*pright)[N2M][N3M][NPR2INTERP];

extern FTYPE (*dq1)[N2M][N3M][NPR2INTERP];
extern FTYPE (*dq2)[N2M][N3M][NPR2INTERP];
extern FTYPE (*dq3)[N2M][N3M][NPR2INTERP];

extern FTYPE (*prc)[N2M][N3M][NPR2INTERP];

extern FTYPE (*ptemparray)[N2M][N3M][NPR];


extern FTYPE (*F1)[N2M][N3M][NPR];
extern FTYPE (*F2)[N2M][N3M][NPR];
extern FTYPE (*F3)[N2M][N3M][NPR];

#if( LIMIT_FLUXC2A_PRIM_CHANGE )
extern FTYPE (*fluxvectemp)[N2M][N3M][NPR];	/* termporary storage for flux */ //atch
#endif

extern FTYPE (*Fa)[N2M][N3M][NPR];
extern FTYPE (*Fb)[N2M][N3M][NPR];

extern FTYPE (*a2cin)[N2M][N3M][NPR];
extern FTYPE (*a2cout)[N2M][N3M][NPR];

//extern FTYPE (*F1CT)[N2M][N3M][NPR];
//extern FTYPE (*F2CT)[N2M][N3M][NPR];
//extern FTYPE (*F3CT)[N2M][N3M][NPR];
extern FTYPE (*pk)[N1M][N2M][N3M][NPR];





///////////////////////////////
//
// grid functions (+1 size larger so can have geometry at upper corners -- say for vector potential or whatever)
//
///////////////////////////////
extern FTYPE (*gcon)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
extern FTYPE (*gcov)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
extern FTYPE (*gcovpert)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM];
extern FTYPE (*gdet)[N2M+SHIFT2][N3M+SHIFT3][NPG];
extern FTYPE (*eomfunc)[N2M+SHIFT2][N3M+SHIFT3][NPG];
//extern FTYPE (*dxdxp)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM];

extern FTYPE (*conn)[N2M][N3M][NDIM][NDIM][NDIM];
extern FTYPE (*conn2)[N2M][N3M][NDIM];
extern FTYPE (*idxvol)[N2M][N3M][NDIM];





// current density stuff
extern FTYPE (*cfaraday)[N2M][N3M][NUMCURRENTSLOTS][3];
extern FTYPE (*fcon)[N2M][N3M][NUMFARADAY];
extern FTYPE (*jcon)[N2M][N3M][NDIM];

// time average stuff
extern FTYPE (*normalvarstavg)[N2M][N3M][NUMNORMDUMP];
extern FTYPE (*anormalvarstavg)[N2M][N3M][NUMNORMDUMP];

extern FTYPE (*jcontavg)[N2M][N3M][NDIM];
extern FTYPE (*jcovtavg)[N2M][N3M][NDIM];
extern FTYPE (*ajcontavg)[N2M][N3M][NDIM];
extern FTYPE (*ajcovtavg)[N2M][N3M][NDIM];

extern FTYPE (*massfluxtavg)[N2M][N3M][NDIM];
extern FTYPE (*amassfluxtavg)[N2M][N3M][NDIM];

extern FTYPE (*othertavg)[N2M][N3M][NUMOTHER];
extern FTYPE (*aothertavg)[N2M][N3M][NUMOTHER];

extern FTYPE (*fcontavg)[N2M][N3M][NUMFARADAY];
extern FTYPE (*fcovtavg)[N2M][N3M][NUMFARADAY];
extern FTYPE (*afcontavg)[N2M][N3M][NUMFARADAY];
extern FTYPE (*afcovtavg)[N2M][N3M][NUMFARADAY];

extern FTYPE (*tudtavg)[N2M][N3M][NUMSTRESSTERMS];
extern FTYPE (*atudtavg)[N2M][N3M][NUMSTRESSTERMS];

// for image, see image.c
extern FTYPE (*pimage)[N2M][N3M][NPR];


//#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
	extern FTYPE (*weno_lower_order_fraction)[N2M][N3M][NPR]; //atch
//#endif

//#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )
	extern FTYPE (*weno_prim_lower_order_fraction)[N1M][N2M][N3M]; //atch
//#endif


//#if(DOENO)
extern FTYPE (*uenotmp0)[N2M][N3M][NPR]; /* for ENO reconstruction of U or dU*/
extern FTYPE (*uenotmp1)[N2M][N3M][NPR]; /* for ENO reconstruction of U or dU*/
extern FTYPE (*uenotmp2)[N2M][N3M][NPR]; /* for ENO reconstruction of U or dU*/
//#endif










/////////////////////////////////////////////////////////////////
//
//
// GLOBAL VARIABLES that are not for every poextern int in space
//
//
/////////////////////////////////////////////////////////////////


extern SFTYPE *lumvsr,*lumvsr_tot;

extern SFTYPE *dissvsr,*dissvsr_tot;


/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
extern FTYPE a;
extern FTYPE gam;
extern FTYPE Bpole,Omegastar; 

/* numerical parameters */
extern int defcoord;
extern FTYPE Rin, R0, Rout, hslope, Zin, Zout;
extern FTYPE Rin_array[NDIM], Rout_array[NDIM];  //atch -- arrays for a more general way of handling the grid dimensions
extern FTYPE Risco,Rhor;
extern FTYPE cour;
extern FTYPE dV, dVF, dx[NDIM], startx[NDIM];
extern SFTYPE dt,t,tf,tstepparti,tsteppartf;
extern FTYPE rcurr, hcurr;

//extern int istart, istop, jstart, jstop;
#if(SIMULBCCALC!=-1) 
	extern int isc,iec,jsc,jec;
	extern int isf1,ief1,jsf1,jef1,ksf1,kef1;
	extern int isf2,ief2,jsf2,jef2,ksf2,kef2;
	extern int isf3,ief3,jsf3,jef3,ksf3,kef3;
	extern int ise,iee,jse,jee;
	extern int isf1ct,ief1ct,jsf1ct,jef1ct;// GODMARK: other stage type requires more
	extern int isf2ct,ief2ct,jsf2ct,jef2ct;
	extern int isf3ct,ief3ct,jsf3ct,jef3ct;
	extern int isdq,iedq,jsdq,jedq;
	extern int ispdq,iepdq,jspdq,jepdq;
#endif

extern FTYPE mydminarg1, mydminarg2;
extern long nstep;

extern int steppart,numstepparts;


/* output parameters */
extern SFTYPE DTd;
extern SFTYPE DTavg;
extern SFTYPE DTener;
extern SFTYPE DTi;
extern SFTYPE DTdebug;
extern long DTr;
extern long dump_cnt;
extern long avg_cnt;
extern long debug_cnt;
extern long image_cnt;
extern long rdump_cnt;
extern long fieldline_cnt; // assumed to keep track with images (as in diag.c), so no need to include in restart()
extern int nstroke;

/* global flags */
extern int failed;
extern int lim,fluxmethod,FLUXB,UTOPRIMVERSION,TIMEORDER,DOENOFLUX,avgscheme,do_transverse_flux_integration,do_conserved_integration,do_source_integration;
extern FTYPE defcon;

/* diagnostics */
// don't track this separately in other regions except global region
extern SFTYPE frdot[N1][NPR];
extern SFTYPE pdottermsjet2[COMPDIM*2][NUMFLUXTERMS][NPR];
extern CTYPE failfloorcountlocal[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet
extern CTYPE failfloorcountlocal_tot[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet

// general stuff for ener.out file for regions to completely track, including terms within flux
extern int dofluxreg[NUMENERREGIONS][COMPDIM*2];
extern int enerposreg[NUMENERREGIONS][COMPDIM*2];
// these quantities contain diagnostics
// all these require writing to restart file
// other _tot quantities appear in dump_ener.c that don't need to be written to restart file since easily computed from existing data.
extern SFTYPE fladdreg[NUMENERREGIONS][NPR];
extern SFTYPE fladdreg_tot[NUMENERREGIONS][NPR];
extern SFTYPE fladdtermsreg[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
extern SFTYPE fladdtermsreg_tot[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
extern SFTYPE Ureg_init[NUMENERREGIONS][NPR];
extern SFTYPE Ureg_init_tot[NUMENERREGIONS][NPR];
extern SFTYPE pcumreg[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pcumreg_tot[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pdotreg[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pdottermsreg[NUMENERREGIONS][COMPDIM*2][NUMFLUXTERMS][NPR];
extern SFTYPE sourceaddreg[NUMENERREGIONS][NPR];
extern SFTYPE sourceaddreg_tot[NUMENERREGIONS][NPR];
extern SFTYPE sourceaddtermsreg[NUMENERREGIONS][NUMSOURCES][NPR];
extern SFTYPE sourceaddtermsreg_tot[NUMENERREGIONS][NUMSOURCES][NPR];
extern SFTYPE dissreg[NUMENERREGIONS][1];
extern SFTYPE dissreg_tot[NUMENERREGIONS][1];

// used for each region, related to global quantities
// _tot quantities here are global since used in restart.
extern int *doflux;
extern int *enerpos;
extern SFTYPE *fladd;
extern SFTYPE *fladd_tot;
extern SFTYPE (*fladdterms)[NPR];
extern SFTYPE (*fladdterms_tot)[NPR];
extern SFTYPE *U_init;
extern SFTYPE *U_init_tot;
extern SFTYPE (*pcum)[NPR];
extern SFTYPE (*pcum_tot)[NPR];
extern SFTYPE (*pdot)[NPR];
extern SFTYPE (*pdotterms)[NUMFLUXTERMS][NPR];
extern SFTYPE *sourceadd;
extern SFTYPE *sourceadd_tot;
extern SFTYPE (*sourceaddterms)[NPR];
extern SFTYPE (*sourceaddterms_tot)[NPR];
extern SFTYPE *diss;
extern SFTYPE *diss_tot;

// end changes after ...

/* current local position */
extern int icurr, jcurr, kcurr, pcurr, ihere, jhere, phere;

/* Jon's addition */
extern int horizoni;
extern long realnstep;
extern int partialstep;
extern int mpicombine;
extern int mpicombinetype;
extern int truempicombinetype;
extern int halftimep;
extern int whichrestart;
extern int appendold;
extern int whocalleducon;
// global flags
extern long restartsteps[2];
extern int binaryoutput,sortedoutput;
extern long steptofaildump,steptofailmap;
extern int ifail,jfail,kfail,dofailmap,dofaildump,restartonfail;
// IC
extern FTYPE h_over_r;
// BC
extern FTYPE h_over_r_jet;
extern int BCtype[COMPDIM*2];
extern int rescaletype;
extern int cooling;
extern int DOENERDIAG,DOGDUMPDIAG,DORDUMPDIAG,DODUMPDIAG,DOAVGDIAG, DOIMAGEDIAG,DOAREAMAPDIAG;
extern int GAMMIEDUMP,GAMMIEIMAGE,GAMMIEENER,DODIAGS,RESTARTMODE,WHICHFILE,POSDEFMETRIC,DOENODEBUGEVERYSUBSTEP,DODIAGEVERYSUBSTEP,
		INVERTFROMAVERAGEIFFAILED,LIMIT_AC_PRIM_FRAC_CHANGE; //atch
extern FTYPE RHOMIN,UUMIN,RHOMINLIMIT,UUMINLIMIT,MAX_AC_PRIM_FRAC_CHANGE; //atch
extern FTYPE prMAX[NPR];
extern FTYPE BSQORHOLIMIT,BSQOULIMIT,GAMMAMAX,GAMMADAMP,GAMMAFAIL;
extern FTYPE SAFE;
extern int debugfail;
extern FTYPE uttdiscr; // for check_pr for now
extern int jonchecks;
extern int dnumcolumns[NUMDUMPTYPES];
extern struct blink * blinkptr0[NUMDUMPTYPES];
extern struct blink * cpulinkptr0[NUMDUMPTYPES];
extern int DOCOLSPLIT[NUMDUMPTYPES];
extern int docolsplit; // global var for now
extern int nextcol;

/* physical consts */
extern FTYPE msun,lsun,G,H,C,mn,me,kb,arad,sigmasb,sigmamat,mevocsq,ergPmev;
extern FTYPE mb,M,Mdot,Mdotc;
extern FTYPE Lunit,Tunit,rhounit,Munit,mdotunit,energyunit,edotunit,Pressureunit,Tempunit,Bunit,massunitPmsun;
extern FTYPE ledd,leddcode;

extern int NUMBUFFERS;

extern int ijkminmax[NUMINTERPTYPES][NDIM][2]; // 0=in/1=out (or similar)
// as for BCtype
extern int ijkminmaxud[NUMINTERPTYPES][COMPDIM*2];
extern int fluxloop[NDIM][NUMFLUXLOOPNUMBERS];
extern int Uconsloop[NUMFLUXLOOPNUMBERS];
extern int interporder[NUMINTERPS];


// ENO DEBUG GLOBAL VARS
extern int dirglobal,iglobal,jglobal,kglobal,iterglobal,interporfluxglobal;

extern FTYPE prfloorcoef[NPR];

extern int numbercpu[ 3+1 ];

extern FTYPE globalinv[NUMGLOBALINV];
extern char globalinvtext[NUMGLOBALINV][10];

// Ramesh stuff
extern FTYPE nu,ss,ucrit,Ttpow,jetalpha;

extern int lntries;
extern FTYPE lerrx;

// EOS related functions
extern FTYPE (*ptr_pressure_rho0_u)(FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_compute_u_from_entropy)(FTYPE rho0, FTYPE entropy);
extern FTYPE (*ptr_u_rho0_p)(FTYPE rho0, FTYPE p);
extern FTYPE (*ptr_dpdu_rho0_u)(FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_dpdrho0_rho0_u)(FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_cs2_compute)(FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_compute_dSdrho)(FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_compute_dSdu)(FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_compute_entropy)(FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_pressure_wmrho0) (FTYPE rho0, FTYPE wmrho0);
extern FTYPE (*ptr_compute_idwmrho0dp) (FTYPE rho0, FTYPE wmrho0);
extern FTYPE (*ptr_compute_idrho0dp) (FTYPE rho0, FTYPE wmrho0);
