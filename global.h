#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#ifndef WIN32
#include <unistd.h>
#endif
#include <float.h>

#ifndef GLOBAL_H
#define GLOBAL_H



#include "metric.h"
#include "coord.h"




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Various physics and model setup parameters that are macros either for performance reasons or since no need to change them at runtime.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////
//
// nmenomics associated with definit.h
//
//
//
////////////////////////
#define NUMJETS 2
#define INNERJET 0
#define OUTERJET 1


#define QUASISTRANG 0
#define UNSPLIT 1

#define GLOBALENERREGION 0
#define INNERJETREGION 1
#define OUTERJETREGION 2

// for VCHARTYPE
#define VERYLOCALVCHAR 0
#define LOCALVCHAR 1
#define GLOBALVCHAR 2



#define MAXCPUS 1500
#define MAXFILENAME 200

// for WHICHVEL
// 0: 4-velocity (leads to ambiguous u^t +- discr part)
// 1: 3-velocity (unambiguous u^t but interpolation is not constrained to be a good 3-velocity)
// 2: relative 4-velocity (unambiguous u^t and any interpolation gives good value)
#define VEL4 0
#define VEL3 1
#define VELREL4 2

// for WHICHEOM
#define WITHGDET 0
#define WITHNOGDET 1
#define WITHSINSQ 2

// This stays naturally, simply consistent with how code evolves conserved quantities.
// for WHICHEOM
#define NOGDETRHO 0
#define NOGDETU0 0
#define NOGDETU1 0
#define NOGDETU2 0
#define NOGDETU3 0
#define NOGDETB1 0
#define NOGDETB2 0
#define NOGDETB3 0
#define NOGDETENTROPY 0





// for RELTYPE
#define RELEOM 0
#define NONRELEOM 1 // NOT FINISHED // NOT RIGHT
// whether relativistic or nonrelativistic EOMs (speed of light limitation)

// for EOMTYPE
// 0 = GRMHD
// 1 = FF(D)E force-free electrodynamics
// for force-free, must turn off:
// ok now, but effectively setup already the below 2 lines implicitly
// global.h : FIXUPAFTERINIT, FIXUPAFTERRESTART,CHECKSOLUTION,LIMADJUST,FLUXADJUST
// global.h FIXUPZONES->FIXUPNOZONES
#define EOMGRMHD 0
#define EOMFFDE 1
#define EOMCOLDGRMHD 2




// mnenomics
// for DOENTROPY
#define DONOENTROPY  0
#define DOEVOLVECOMPAREENTROPY 1
#define DOEVOLVEDIRECTENTROPY 2


// for WHICHENTROPYEVOLVE
#define EVOLVENOENTROPY 0
#define EVOLVESIMPLEENTROPY 1 // should be used with DOENTROPY==DOEVOLVECOMPAREENTOPY
#define EVOLVEFULLENTROPY 2 // should only be used with DOENTROPY==DOEVOLVEDIRECTENTROPY or DOENTROPY==DOEVOLVECOMPAREENTOPY



// defines how Utoprimgen is used
#define EVOLVEUTOPRIM 0
#define OTHERUTOPRIM 1


// defines data return types for primtoU() and primtoflux()
#define UEVOLVE 0
#define UDIAG 1
#define UNOTHING 2
#define UENTROPY 3 // implicit UNOTHING but cons UU is overwritten by cons entropy (see primtoflux() in phys.c as used by utoprim() in utoprim.c)




// for LIMADJUST
// 0: use fixed limiter
// 1: use limiter based upon b^2/rho
// 2: use limiter based upon b^2/u
// 3: use limiter based upon both b^2/rho or b^2/u
#define LIMITERFIXED 0
#define LIMITERBSQORHO 1
#define LIMITERBSQOU 2
#define LIMITERBSQORHOANDU 3

// for FLUXADJUST
#define FLUXFIXED 0 // (see get_bsqflags() in fixup.c)
#define FLUXBSQORHO 1
#define FLUXBSQOU 2
#define FLUXBSQORHOANDU 3


// for UTOPRIMFAILRETURNTYPE  --  controls the behaviour of inversion: does allow the return of solutions with negative densities, etc.
#define UTOPRIMRETURNNOTADJUSTED 0
#define UTOPRIMRETURNADJUSTED 1



// for UTOPRIMADJUST  -- controls the behaviour of fixups:  UTOPRIMAVG means fix it up, UTOPRIMSTATIC means do not do it
// 0=just use static solution
// 1=use average surrounding solution, and if no good surrounding solution use the normal observer velocity with static densities
#define UTOPRIMSTATIC 0
#define UTOPRIMAVG 1



// used to choose which method interp.c uses
#define NUMINTERPTYPES 9

#define NONENOINTERPTYPE 0
#define ENOINTERPTYPE 1 // ce2
#define ENOFLUXRECONTYPE 2
#define ENOFLUXSPLITTYPE 3
#define ENOAVG2CENTTYPE 4
#define ENOCENT2AVGTYPE 5
#define ENOFLUXAVG1TYPE 6
#define ENOFLUXAVG2TYPE 7
#define ENOFLUXAVG3TYPE 8


//quantities to interp
#define ENOSOURCETERM 0
#define ENOCONSERVED  1
#define ENOPRIMITIVE  2
#define ENOFLUX       3


#define NOENOFLUX 0
#define ENOFLUXRECON 1
#define ENOFLUXSPLIT 2
#define ENOFINITEVOLUME 3
// 0: no ENO flux reconstruction
// 1: reconstruct F for finite difference rep. of U
// 2 : flux splitting (not done yet)
// 3: reconstruct dU for finite volume rep. of U


#define NUMDISSFUNPOS 2


// see interp_loop_set() in initbase.c
#define NUMFLUXLOOPNUMBERS 10
#define FIDEL 0
#define FJDEL 1
#define FKDEL 2
#define FFACE 3
#define FIS 4
#define FIE 5
#define FJS 6
#define FJE 7
#define FKS 8
#define FKE 9


// number of inversion quantities to report when inversion fails if CHECKONINVERSION = 1
#define NUMGLOBALINV 13

// for WHICHCURRENTCALC
// 0: original time is on edge and spatial on edge, but spatials are different locations.  old time.
// 1: all centered in space and all time, present time (best)
// 2: like 0, but spatially centered (i.e. old time)
#define CURRENTCALC0 0
#define CURRENTCALC1 1
#define CURRENTCALC2 2



#define CURRENTPRECALCTYPES 5

#define CURTYPET 0
#define CURTYPEX 1
#define CURTYPEY 2
#define CURTYPEZ 3
#define CURTYPEFARADAY 4


// whether and which type of fixups to be used
#define FIXUP1ZONE 0
#define FIXUPALLZONES 1
#define FIXUPNOZONES 2



/* mnemonics for flux method (Riemann solver) */
// ordered from most diffusive to least diffusive, so can back track
// 0 should be reasonable most diffusive
#define LAXFFLUX 0
#define HLLFLUX 1
#define FORCEFLUX 2
#define MUSTAFLUX 3
#define HLLLAXF1FLUX 4


// DIVB constraint method
#define FLUXCTHLL 0
#define FLUXCTTOTH 1
#define FLUXCD 2
#define ATHENA1 3
#define ATHENA2 4
/* these are different ways of calculating the EMFs */
//#define FLUXB FLUXCTTOTH
// 0: HLL
// 1: FLUXCT TOTH version (toth 2000 eq. 25)
// 2: FLUXCD TOTH version (toth 2000 eq. 31)
// 3: Athena type eq 39
// 4: Athena type eq 48


//#define UTOPRIMVERSION 6
// 0: original gammie 5D method
#define UTOPRIM5D1 0
// 1: ldz method
#define UTOPRIMLDZ 1
// 2: SCN 2D method
#define UTOPRIM2D 2
// 3: SCN 1D method
#define UTOPRIM1D 3
// 4: SCN 1D OPTIMIZED method -- not sure if identical to 3 otherwise
#define UTOPRIM1DOPT 4
// 5: SCN 1D final and optimized
#define UTOPRIM1DFINAL 5
// 6: SCN 2D final and optimized and recommended by Scott
#define UTOPRIM2DFINAL 6
// 7: SCN 5D final -- bit less accurate compared to 1D and 2D
#define UTOPRIM5D2 7
// 8: Jon 1D/2D final version -- can handle non-rel problems
#define UTOPRIMJONNONRELCOMPAT 8
// 100: use 5D, but compare with ldz in runtime
#define UTOPRIMCOMPARE 100


/* mnemonics for slope limiter */
// ordered from most diffusive to least diffusive, so can start high and go down if needed
// 0 should be reasonble most diffusive, highest should be least diffusive
#define NUMINTERPS 18

#define DONOR	0
#define VANL	1
#define MINM	2
#define MC      3
#define PARA    4
#define PARAFLAT 5
#define CSSLOPE      6 // not tested/compared against others
#define WENO3 7
#define WENO4 8
#define WENO5  9
#define WENO6  10
#define WENO7  11
#define WENO8  12
#define WENO9  13
#define ENO3 14
#define ENO5 15
#define WENO5FLAT 16
#define WENO5BND 17

// negative versions for testing only
#define NLIM    -1 // no limiter
#define NLIMCENT    -2 // no limiter
#define NLIMUP    -3 // no limiter
#define NLIMDOWN    -4 // no limiter


// see orders_set() in initbase.c
#define MAXSPACEORDER 11 // maximum number of points in stencil
//#define MAXSPACESHIFT ((MAXSPACEORDER-1)/2) // center point for symmetric stencil


#define MAXTIMEORDER 4

//#define TIMEORDER 3
// order of algorithm in time from 1 to 4.
// TIMEORDER: 1 : single step (Euler method -- error term is 2nd order for smooth flows)
// TIMEORDER: 2 : 2 steps in halfs (midpoint method -- error term is 3rd order for smooth flows)
// TIMEORDER: 3 : 4 steps (classic RK3 method -- error term is 4th order for smooth flows)
// TIMEORDER: 4 : 4 steps (classic RK4 method -- error term is 5th order for smooth flows)

//////////////////////////////////
//
// which variable to interpolate
//
/////////////////////////////////

#define PRIMTOINTERP -1
#define PRIMTOINTERP_JONRESCALED1 0
#define CONSTOINTERP 1
#define PRIMTOINTERPLGDEN 2
#define PRIMTOINTERP_LGDEN_RHOU 3
#define PRIMTOINTERP_RHOU 4
#define PRIMTOINTERP_VSQ 5
#define PRIMTOINTERP_3VEL_GAMMA 6
#define PRIMTOINTERP_RHOV_GAMMA 7
#define PRIMTOINTERP_VELREL4SQ 8
#define PRIMTOINTERP_3VELREL_GAMMAREL 9
#define PRIMTOINTERP_RAMESH1 10


#define WENO_REDUCE_TYPE_DEFAULT 0 
#define WENO_REDUCE_TYPE_PPM 1

// definition of minmod operator
#define MINMOD(a,b) ( ((a)*(b)<=0) ? 0.0 : MINMODB(a,b) )
#define MINMODB(a,b) ( (fabs(a)<fabs(b)) ? (a) : (b) )
//#define MINMOD3( x, y, z )   ( 0.25 * (sign(x) + sign(y)) * (sign(x) + sign(z)) * MIN( MIN(fabs(x), fabs(y)), fabs(z)) )    

///////////////////////////////
//
// parabolic interpolation stuff
//
////////////////////////////////

#define PARA1 0 // old
#define PARA2 1 // works
#define PARA3 2 // broken
#define PARA4 3 // latest


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Default global and user choices for various code options that can change for each run without significant modifications
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// default values
#include "definit.h"

// user specific values
#include "init.h"





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// MNEMONICS or other things that should rarely change, or things that on depend on the above items (e.g. if statements, loops, etc.)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//equals to unity if the interpolation of gamma is performed and if requested to use prim. reduction
#define STORE_GAMMA_PRIM_REDUCTION_FRACTION (WENO_USE_PRIM_REDUCTION && (VARTOINTERP == PRIMTOINTERP_3VEL_GAMMA || VARTOINTERP == PRIMTOINTERP_RHOV_GAMMA || VARTOINTERP == PRIMTOINTERP_3VELREL_GAMMAREL) )



//
// below for para2() and para3()
//

#define PMC 100
// PMC: originally used in PPM paper -- Xiaoyue says fails on nonlinear tests since has no limiter -- but they used it, right?

#define MINM_STEEPENER 101

#if(WHICHPARA==PARA2)
#define PARA2LIM VANL
// jon tests show for PARA2:
// PARA2LIM = MC completely fails
// PARA2LIM = VANL best ok
// PARA2LIM = PMC kinda ok
#elif(WHICHPARA==PARA3)
// MC recommended by Xiaoyue
#define PARA2LIM MC
#elif(WHICHPARA==PARA4)
// MC recommended by Xiaoyue
#define PARA2LIM MC
//#define PARA2LIM MINM
#endif




#if(DOJETDIAG)
#define NUMENERREGIONS 3
#else
#define NUMENERREGIONS 1
#endif

// whether to check if rho<=0 or u<=0 when restarting.
#if( (EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD) )
#define CHECKRHONEGZERORESTART 1
#else
#define CHECKRHONEGZERORESTART 0
#endif


#if((WHICHCURRENTCALC==0)||(WHICHCURRENTCALC==2))
#define NUMCURRENTSLOTS 5
#elif(WHICHCURRENTCALC==1)
#define NUMCURRENTSLOTS 6
#endif





//////////////////////////////////////
//
// PURE mnemonics
//
///////////////////////////////////
#define NUMBOUNDTYPES 2
#define BOUNDPRIMTYPE 0
#define BOUNDFLUXTYPE 1


// -------------> r
// |	     3    
// |	    1-0   
// |	     2    
// v	        
// theta      
// and likewise for 4,5 (4=out,5=in)

#define X1UP	0
#define X1DN	1
#define X2UP	2
#define X2DN	3
#define X3UP	4
#define X3DN	5

// long double constants
# define M_El           2.7182818284590452353602874713526625L  /* e */
# define M_LOG2El       1.4426950408889634073599246810018922L  /* log_2 e */
# define M_LOG10El      0.4342944819032518276511289189166051L  /* log_10 e */
# define M_LN2l         0.6931471805599453094172321214581766L  /* log_e 2 */
# define M_LN10l        2.3025850929940456840179914546843642L  /* log_e 10 */

#ifdef WIN32
# define M_PI           3.1415926535897932384626433832795029  /* pi */
#endif

# define M_PIl          3.1415926535897932384626433832795029L  /* pi */
# define M_PI_2l        1.5707963267948966192313216916397514L  /* pi/2 */
# define M_PI_4l        0.7853981633974483096156608458198757L  /* pi/4 */
# define M_1_PIl        0.3183098861837906715377675267450287L  /* 1/pi */
# define M_2_PIl        0.6366197723675813430755350534900574L  /* 2/pi */
# define M_2_SQRTPIl    1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
# define M_SQRT2l       1.4142135623730950488016887242096981L  /* sqrt(2) */
# define M_SQRT1_2l     0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
# define SIXTH          0.1666666666666666666666666666666666L  /* 1/6 */
# define FOURTHIRD      1.3333333333333333333333333333333333L  /* 4/3 */
# define THIRD          0.3333333333333333333333333333333333L  /* 1/3 */
# define ONE            1.0000000000000000000000000000000000L
# define PTFIVE         0.5L
# define TWO            2.0L
# define ONEPT25        1.25L
# define THREE          3.0L
# define SIX            6.0L
# define EIGHT          8.0L


// setup for various boundary situations
// so doesn't produce differences in irrelevant directions, whether boundary zones or not
// mac(?) macros are for use in definitions of other macros since macro with args needs to be directly a function of what's hardcoded, not some "replacement" since nothing is replaced, not to be used inside code for USING a macro.

#if(N1>1)
#define im1 (i-1)
#define im1mac(i) (i-1)
#define ip1 (i+1)
#define ip1mac(i) (i+1)

#define irefshift (2*N1-i-1)

#else

#define im1 (i)
#define im1mac(i) (i)
#define ip1 (i)
#define ip1mac(i) (i)

#define irefshift (i)

#endif

#if(N2>1)
#define jm1 (j-1)
#define jm1mac(j) (j-1)
#define jp1 (j+1)
#define jp1mac(j) (j+1)
#define jrefshift (2*N2-j-1)
#else
#define jm1 (j)
#define jm1mac(j) (j)
#define jp1 (j)
#define jp1mac(j) (j)
#define jrefshift (j)
#endif

#if(N3>1)
#define km1 (k-1)
#define km1mac(k) (k-1)
#define kp1 (k+1)
#define kp1mac(k) (k+1)
#define krefshift (2*N3-k-1)
#else
#define km1 (k)
#define km1mac(k) (k)
#define kp1 (k)
#define kp1mac(k) (k)
#define krefshift (k)
#endif


// used to tell if N>1 or not (can just ask directly)
#define N1NOT1 ((N1>1) ? 1 : 0)
#define N2NOT1 ((N2>1) ? 1 : 0)
#define N3NOT1 ((N3>1) ? 1 : 0)

// GODMARK: looks like I can set just above and rest are set for me for any case

// GODMARK: check these new conditions

// restrict loops only over relevant domain in reduced dimension case



// 3 maximum boundary zones needed if doing Parabolic interpolation
// maximum number of boundary zones needed for all calculations

// have to set manually if going to set DOENOFLUX at runtime.



// x1
#define SHIFT1 N1NOT1
#define N1BND MAXBND*N1NOT1

#define INFULL1 (-N1BND)
#define INFULLP11 (-N1BND+SHIFT1)
#define OUTFULL1 (N1-1+N1BND) 
#define OUTFULLM11 (N1-1+N1BND-SHIFT1) 
#define OUTFULLP11 (N1-1+N1BND+SHIFT1)

#define INHALF1 (-N1BND/2)
#define OUTHALF1 (N1-1+N1BND/2)
#define INP11 (-N1BND+SHIFT1)
#define OUTP11 (N1-1+N1BND-SHIFT1)
#define INM1 -SHIFT1
#define OUTM1 N1-1+SHIFT1

// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO1 (-N1BND)
#define INBOUNDHI1 (-1)

#define OUTBOUNDLO1 (N1)
#define OUTBOUNDHI1 (N1+N1BND-1)

//#define INFACEBOUNDLO1 (-N1BND) // (-N1BND+1) // GODMARK: large domain used for easy checking of fluxes after bound_flux().
//#define INFACEBOUNDHI1 (-1+SHIFT1)

// up to -1 since 0 is actually defined with original primitives
#define INFACEBOUNDLO1 (-N1BND+1)
#define INFACEBOUNDHI1 (-1)


// from N1+1 since N1 is actually defined with original primitives
#define OUTFACEBOUNDLO1 (N1+1)
#define OUTFACEBOUNDHI1 (N1+N1BND-1)

// x2
#define SHIFT2 N2NOT1
#define N2BND MAXBND*N2NOT1

#define INFULL2 (-N2BND)
#define INFULLP12 (-N2BND+SHIFT2)
#define OUTFULL2 (N2-1+N2BND)
#define OUTFULLM12 (N2-1+N2BND-SHIFT2) 
#define OUTFULLP12 (N2-1+N2BND+SHIFT2)

#define INHALF2 (-N2BND/2)
#define OUTHALF2 (N2-1+N2BND/2)
#define INP12 (-N2BND+SHIFT2)
#define OUTP12 (N2-1+N2BND-SHIFT2)
#define INM2 -SHIFT2
#define OUTM2 N2-1+SHIFT2

// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO2 (-N2BND)
#define INBOUNDHI2 (-1)

#define OUTBOUNDLO2 (N2)
#define OUTBOUNDHI2 (N2+N2BND-1)

//#define INFACEBOUNDLO2 (-N2BND)
//#define INFACEBOUNDHI2 (-1+SHIFT2)

#define INFACEBOUNDLO2 (-N2BND+1)
#define INFACEBOUNDHI2 (-1)


#define OUTFACEBOUNDLO2 (N2+1)
#define OUTFACEBOUNDHI2 (N2+N2BND-1)


// x3
#define SHIFT3 N3NOT1
#define N3BND MAXBND*N3NOT1

#define INFULL3 (-N3BND)
#define INFULLP13 (-N3BND+SHIFT3)
#define OUTFULL3 (N3-1+N3BND)
#define OUTFULLM13 (N3-1+N3BND-SHIFT3)
#define OUTFULLP13 (N3-1+N3BND+SHIFT3)

#define INHALF3 (-N3BND/2)
#define OUTHALF3 (N3-1+N3BND/2)
#define INP13 (-N3BND+SHIFT3)
#define OUTP13 (N3-1+N3BND-SHIFT3)
#define INM3 -SHIFT3
#define OUTM3 N3-1+SHIFT3

// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO3 (-N3BND)
#define INBOUNDHI3 (-1)

#define OUTBOUNDLO3 (N3)
#define OUTBOUNDHI3 (N3+N3BND-1)

//#define INFACEBOUNDLO3 (-N3BND)
//#define INFACEBOUNDHI3 (-1+SHIFT3)

#define INFACEBOUNDLO3 (-N3BND+1)
#define INFACEBOUNDHI3 (-1)


#define OUTFACEBOUNDLO3 (N3+1)
#define OUTFACEBOUNDHI3 (N3+N3BND-1)




/* allocated memory uses this for active zones 0-N1-1 and bc beyond that */
#define N1M (N1+N1BND*2)
#define N2M (N2+N2BND*2)
#define N3M (N3+N3BND*2)


/* NBIG is bigger of N1 and N2 and N3 */
#define NBIG1 ((N1>N2) ? N1 : N2)
#define NBIG  ((NBIG1>N3) ? NBIG1 : N3)

#define NBIGBND1 ((N1BND>N2BND) ? N1BND : N2BND)
#define NBIGBND  ((NBIGBND1>N3BND) ? NBIGBND1 : N3BND)

// N?OFF and N?NOT1 are a bit redundant
//#define N1OFF (((N1BND>0)&&(N1>1)) ? 1 : 0)
//#define N2OFF (((N2BND>0)&&(N2>1)) ? 1 : 0)
//#define N3OFF (((N3BND>0)&&(N3>1)) ? 1 : 0)



/* NBIGM is bigger of N1M and N2M and N3M */
// not currently used
#define NBIG1M ((N1M>N2M) ? N1M : N2M)
#define NBIGM  ((NBIG1M>N3M) ? NBIG1M : N3M)

// surface areas of sides WITH boundary zones
// reduces to 0 if that dimension not used
#define SURFA1 (N2M*N3M*N1NOT1)
#define SURFA2 (N1M*N3M*N2NOT1)
#define SURFA3 (N1M*N2M*N3NOT1)

// maximal surface of boundary exchange
// notice that this works in any number of dimensions and any N1,N2,N3
// that is, it reduces correctly when a dimension is degenerate (N=1)
#define NBIGS1M ((SURFA1>SURFA2) ? SURFA1 : SURFA2)
#define NBIGSM ((NBIGS1M>SURFA3) ? NBIGS1M : SURFA3)

// GODMARK: Could make a volume that is not NBIGBND*NBIGSM but may be smaller?
// used in init_mpi.c for workbc and workbc_int




// for wavespeeds.c and fluxcompute.c
#define NUMCS 2
#define CMIN 0
#define CMAX 1

#define MINMAX(q,a,b) ( ((q)==CMIN) ? MIN(a,b) : MAX(a,b) )



#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define SIGN(a) ( ((a) <0.) ? -1. : 1. )

#define PROGRADERISCO 0
#define RETROGRADERISCO 1

#define NUMTSCALES 4
// number of times scales to watch failure rates at
#define ALLTS 0 // full cumulative
#define ENERTS 1 // cumulative each dump_ener (over all grid)
#define IMAGETS 2 // cumulative each image dump (full grid)
#define DEBUGTS 3 // debug dump time scale (full grid)

// dump.c's fieldlinedump()
#define NUMFIELDLINEQUANTITIES 11
// rho, u, -hu_t, -T^t_t/U0, u^t, v1,v2,v3,B1,B2,B3

// see failfloorcount counter
#define NUMFAILFLOORFLAGS 9
//  mnemonics
#define COUNTUTOPRIMFAILCONV 0 // if failed to converge
#define COUNTFLOORACT 1 // if floor activated
#define COUNTLIMITGAMMAACT 2 // if Gamma limiter activated
#define COUNTINFLOWACT 3 // if inflow check activated
#define COUNTUTOPRIMFAILRHONEG 4
#define COUNTUTOPRIMFAILUNEG 5
#define COUNTUTOPRIMFAILRHOUNEG 6
#define COUNTGAMMAPERC 7 // see fixup_checksolution()
#define COUNTUPERC 8 // see fixup_checksolution()

// failure codes for utoprim failures
#define UTOPRIMFAILFIXED -1
#define UTOPRIMNOFAIL 0
#define UTOPRIMFAILCONVRET   101
#define UTOPRIMFAILCONV 1
#define UTOPRIMFAILCONVW     2
#define UTOPRIMFAILCONVUTSQ  3
#define UTOPRIMFAILCONVGUESSUTSQ 4
#define UTOPRIMFAILCONVUTSQVERYBAD 5
#define UTOPRIMFAILCONVBADINVERTCOMPARE 6
#define UTOPRIMFAILNANGUESS 7
#define UTOPRIMFAILNANRESULT 8
#define UTOPRIMFAILRHONEG 9
#define UTOPRIMFAILUNEG 10
#define UTOPRIMFAILRHOUNEG 11
#define UTOPRIMFAILGAMMAPERC 12
#define UTOPRIMFAILUPERC 13
#define UTOPRIMFAILU2AVG1 14
#define UTOPRIMFAILU2AVG2 15

/* failure modes */
#define FAIL_UTOPRIM_NEG	1
#define FAILSTR01 "UTOPRIM_NEG"
#define FAIL_UTOPRIM_TEST	2
#define FAILSTR02 "UTOPRIM_TEST"
#define FAIL_VCHAR_DISCR	3
#define FAILSTR03 "VCHAR_DISCR"
#define FAIL_COEFF_NEG		4
#define FAILSTR04 "COEFF_NEG"
#define FAIL_COEFF_SUP		5
#define FAILSTR05 "COEFF_SUP"
#define FAIL_UTCALC_DISCR	6
#define FAILSTR06 "UTCALC_DISCR"
#define FAIL_LDZ	        7
#define FAILSTR07 "FAIL_LDZ"
#define FAIL_BCFIX	        8
#define FAILSTR08 "FAIL_BCFIX"
#define FAIL_VSQ_NEG	        9
#define FAILSTR09 "FAIL_VSQ_NEG"


#if((DODISS||DODISSVSR)&&(DOENTROPY==DONOENTROPY))
#error Turn on entropy evolution if want dissipation
#endif




/* size of global arrays */
#if(DOENTROPY==DONOENTROPY) // normal total energy equation
#define NPR	8		/* number of primitive variables */
#define NPRDUMP NPR
#define NPRINVERT 5
#define NPRBOUND NPR
#define NFLUXBOUND NPR
#else
#define NPR	9		/* number of primitive variables */
#define NPRDUMP 8
#define NPRINVERT 5
#define NPRBOUND 8 // don't care about entropy for primitive bounding since entropy just simple function of primitives and don't use entropy in ghost zones
#define NFLUXBOUND NPR // must be equal to NPR for BOUNDFLUXTYPE
#endif

// number of interpolated quantities (independent from actual list of primitives)
#if(DOEXTRAINTERP)
#define NPR2INTERP (NPR+1)
#else
#define NPR2INTERP (NPR)
#endif

#define NMAXBOUND ((NPRBOUND>NFLUXBOUND) ? NPRBOUND : NFLUXBOUND)

#define NDIM	4		/* number of total dimensions.  Never
				   changes */

// cent,face1,face2,face3,corn



// flag failures/problems for correction/check in fixup
#define NUMPFLAGS (5)
// the below needs to be bounded since one CPU doesn't know if the other failed, and neighbor failure determines nature of how failure is treated
// also, bounded values at real boundaries need to identify if bad copy
#define FLAGUTOPRIMFAIL 0 // changes behavior of fixup()
// the below flags are done after bound_prim, and can be determined at any time, so just come after bound.
#define FLAGREALLIM 1 // value of limiter to be used
#define FLAGBSQORHO 2 // set when B^2/RHO > BSQORHOLIMIT ; currently changes  behavior of slope_lim
#define FLAGBSQOU 3 // set when B^2/u > BSQOULIMIT
#define FLAGREALFLUX 4 // type of flux to use


#define NUMSOURCES 7
// number of source terms.  Currently includes: 1) geometry, 2) radiative cooling, and 3-7) 3=total neutrino cooling (+ the 4 processes)
#define GEOMSOURCE 0
#define RADSOURCE 1
#define NEUTRINOSOURCE 2 // total of all neutrino processes
#define NEUTRINOANN 3
#define NEUTRINOPLASMA 4
#define NEUTRINOECAP 5
#define NEUTRINOBREM 6



// max number of terms in stress tensor (for term-level flux diagnostic)
#define NUMFLUXTERMS (7)


#define ORDERDEBUG 3
//#define NUMENODEBUGS (NPR*(3+ORDERDEBUG*3 + 1 + 2 + ORDERDEBUG*2 + 1))
// p p_l p_r  order*3 per point reduce
// NPR*3 order*3*NPR   NPR
// Uavg Upoint  order*2 per point reduce
// NPR*2 order*2*NPR NPR

//#define NUMENODEBUGS 21
// see email

// short switches:
//1)       SMONO (0,1)
//2)       WENO5 (0,1)
//3)       WENO3 (0,1)
//4)       -> dP/P (0,1)
//5)       limit c2e/c2a/a2c correction (0,1) through checking the change of the quantity being interpolated
//6)       limit c2e/c2a/a2c correction (0,1) through checking the change of primitives

#define NUMENODEBUGS 6


// allow for pk[0] and pk[1]
#define MAXDTSTAGES 2
// maximum number of allowed temporal integration stages



/* mnemonics for primitive vars; conserved vars */
#define RHO	0
#define UU	1
#define U1	2
#define U2	3
#define U3	4
#define B1	5
#define B2	6
#define B3	7
#define ENTROPY 8

#if(DOENTROPY!=DONOENTROPY)
#define VSQ (ENTROPY+1)
#else
#define VSQ (B3+1)
#endif



/* mnemonics for dimensional indices */
#define TT	0
#define RR	1
#define TH	2
#define PH	3


/* mnemonics for centering of grid functions */
// GODMARK: is there a way to pick and choose the dimension and number of grid positions?
/* number of positions on grid for grid 
				   functions */
#define NPG 8
#define FACE1	0
#define FACE2	1
#define FACE3   2
#define CORN1	3 // corner in 2-3 plane
#define CORN2	4 // corner in 1-3 plane
#define CORN3	5 // corner in 1-2 plane
#define CORNT   6 // true corner: full 3D corner (only required for 3D)
#define CENT	7

/* mnemonics for diagnostic calls */
#define INIT_OUT	0
#define DUMP_OUT	1
#define IMAGE_OUT	1
#define LOG_OUT		1
#define FINAL_OUT	2

// size of certain dumped tavg quantities

// was 29 =>
#define NUMNORMDUMP (NPR+1+4*4+6) // number of "normal" dump variables
// for above see diag.c and set_varstavg()
#define NUMFARADAY 6
#define NUMOTHER 1
#define NUMSTRESSTERMS (NUMFLUXTERMS*NDIM*NDIM)

/** GLOBAL ARRAY SECTION **/






// size of data type used for all floats
#define FLOATTYPE 0
#define DOUBLETYPE 1
#define LONGDOUBLETYPE 2
#define LONGLONGINTTYPE 3


#if(REALTYPE>SENSITIVE)
god=deathadflkjasdflkjasdlfkja242424
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON (1.19209290e-07F)
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON (2.2204460492503131e-16)
#endif

#ifndef LDBL_EPSILON
#define LDBL_EPSILON (1.08420217248550443401e-19L)
#endif


// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define NUMEPSILON FLT_EPSILON
#define FTYPE float
#elif(REALTYPE==DOUBLETYPE)
#define NUMEPSILON DBL_EPSILON
#define FTYPE double
#elif(REALTYPE==LONGDOUBLETYPE)
#define NUMEPSILON LDBL_EPSILON
#define FTYPE long double
#endif

#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define SFTYPE float
#elif(SENSITIVE==DOUBLETYPE)
#define SFTYPE double
#elif(SENSITIVE==LONGDOUBLETYPE)
#define SFTYPE long double
#endif


// used for numerical differencing
#define NUMSQRTEPSILON (sqrt(NUMEPSILON))



#if(COUNTTYPE==DOUBLETYPE)
#define CTYPE double
#elif(COUNTTYPE==LONGLONGINTTYPE)
#define CTYPE long long int
#endif


// GODMARK: NUMENERVAR outdated?
#define NUMENERVAR (6+NPR+NPR+3)


/* numerical convenience */
#define BIG (1.e+50)
#define SMALL	(1.e-50)

#define SLEPSILON (1.e-6)


/* size of step in numerical derivative evaluations */
#define HSTEP	1.e-5


#define SURFACETOTAL 0
#define VOLUMETOTAL 1

#define CONSTYPE 0
#define SURFACETYPE 1
#define CUMULATIVETYPE 2
#define CONSJETINNERTYPE 3
#define CONSJETOUTERTYPE 4
#define CUMULATIVETYPE2 5
#define CONSTYPE2 6

#define WRITEHEAD 0
#define READHEAD 1

#define TIMESERIESAREAMAP 0
#define FINALTDUMPAREAMAP 1

#define WRITEFILE 0
#define READFILE 1

#define ENERFNAME "ener.out"
#define GENERFNAME "gener.out"

#define NUMDUMPTYPES 10

#define IMAGECOL 0
#define RDUMPCOL 1
#define DUMPCOL 2
#define GDUMPCOL 3
#define AVGCOL 4
#define AVG2COL 5 // used when needing AVG2COL to avoid too large a file size for avgdump
#define DEBUGCOL 6
#define FIELDLINECOL 7
#define ENODEBUGCOL 8
#define DISSDUMPCOL 9

#if(SENSITIVE==LONGDOUBLETYPE)
// assume sensitive>=realtype in precision
#if(REALTYPE==LONGDOUBLETYPE) // was FLOATTYPE==REALTYPE and SENS=DOUBLETYPE
#define HEADERONEIN "%Lf"
#define HEADER2IN "%Lf %Lf"
#define HEADER3IN "%Lf %Lf %Lf"
#define HEADER4IN "%Lf %Lf %Lf %Lf"
#define HEADER5IN "%Lf %Lf %Lf %Lf %Lf"
#define HEADER6IN "%Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER7IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER8IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER9IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define RESTARTHEADER "%d %d %d "\
	      "%Lf %Lf %ld %Lf %Lf %Lf "\
	      "%Lf %Lf %Lf %Lf %ld %Lf %ld %ld %ld %ld %ld "\
	      "%Lf %d %d %d %d %d %d "\
	      "%Lf %Lf %Lf %Lf %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#elif(REALTYPE==DOUBLETYPE)
#define HEADERONEIN "%lf"
#define HEADER2IN "%lf %lf"
#define HEADER3IN "%lf %lf %lf"
#define HEADER4IN "%lf %lf %lf %lf"
#define HEADER5IN "%lf %lf %lf %lf %lf"
#define HEADER6IN "%lf %lf %lf %lf %lf %lf"
#define HEADER7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADER8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADER "%d %d %d "\
	      "%Lf %Lf %ld %lf %lf %lf "\
	      "%Lf %Lf %Lf %Lf %ld %Lf %ld %ld %ld %ld %ld "\
	      "%Lf %d %d %d %d %d %d "\
	      "%lf %lf %lf %lf %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "\
	      "%lf %lf %ld %f %f %f "\
	      "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "\
	      "%lf %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==DOUBLETYPE)
// assume sensitive>=realtype in precision
#if(REALTYPE==DOUBLETYPE)
#define HEADERONEIN "%lf"
#define HEADER2IN "%lf %lf"
#define HEADER3IN "%lf %lf %lf"
#define HEADER4IN "%lf %lf %lf %lf"
#define HEADER5IN "%lf %lf %lf %lf %lf"
#define HEADER6IN "%lf %lf %lf %lf %lf %lf"
#define HEADER7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADER8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADER "%d %d %d "\
	      "%lf %lf %ld %lf %lf %lf "\
	      "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "\
	      "%lf %d %d %d %d %d %d "\
	      "%lf %lf %lf %lf %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "\
	      "%lf %lf %ld %f %f %f "\
	      "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "\
	      "%lf %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==FLOATTYPE)

#if(REALTYPE==DOUBLETYPE)
#define RESTARTHEADER "" // dumb, so crash on compile
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "\
	      "%f %f %ld %f %f %f "\
	      "%f %f %f %f %ld %f %ld %ld %ld %ld %ld "\
	      "%f %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#endif


// 2
// 6
// 10
// 3
// 5
// 4
// SUM=30
// 23+10=33
// total=63
#define WRITERESTARTHEADER "%d %d %d " \
		 "%21.15g %21.15g %ld %21.15g %21.15g %21.15g " \
		 "%21.15g %21.15g %21.15g %21.15g %ld %21.15g %ld %ld %ld %ld %ld" \
		 "%21.15g %d %d %d %d %d %d " \
		 "%21.15g %21.15g %21.15g %21.15g %d " \
                 "%d %d %d %d %d %d " \
                 "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "

#define HEADERONEOUT "%21.15g "









//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// LOOPS
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// these loops used for general purposes
#define LOOPF3 for(k=INFULL3;k<=OUTFULL3;k++)
#define LOOPF2 for(j=INFULL2;j<=OUTFULL2;j++)
#define LOOPF1 for(i=INFULL1;i<=OUTFULL1;i++)

// boundary loops for CENT quantities
#define LOOPBOUND1IN for(i=INBOUNDLO1;i<=INBOUNDHI1;i++)
#define LOOPBOUND1OUT for(i=OUTBOUNDLO1;i<=OUTBOUNDHI1;i++)

#define LOOPBOUND2IN for(j=INBOUNDLO2;j<=INBOUNDHI2;j++)
#define LOOPBOUND2OUT for(j=OUTBOUNDLO2;j<=OUTBOUNDHI2;j++)

#define LOOPBOUND3IN for(k=INBOUNDLO3;k<=INBOUNDHI3;k++)
#define LOOPBOUND3OUT for(k=OUTBOUNDLO3;k<=OUTBOUNDHI3;k++)

// boundary loops for FACE quantities
#define LOOPFACEBOUND1IN for(i=INFACEBOUNDLO1;i<=INFACEBOUNDHI1;i++)
#define LOOPFACEBOUND1OUT for(i=OUTFACEBOUNDLO1;i<=OUTFACEBOUNDHI1;i++)

#define LOOPFACEBOUND2IN for(j=INFACEBOUNDLO2;j<=INFACEBOUNDHI2;j++)
#define LOOPFACEBOUND2OUT for(j=OUTFACEBOUNDLO2;j<=OUTFACEBOUNDHI2;j++)

#define LOOPFACEBOUND3IN for(k=INFACEBOUNDLO3;k<=INFACEBOUNDHI3;k++)
#define LOOPFACEBOUND3OUT for(k=OUTFACEBOUNDLO3;k<=OUTFACEBOUNDHI3;k++)


// full loop + 1 on outer edge for emf or corner quantities
#define LOOPFP13 for(k=INFULL3;k<=OUTFULLP13;k++)
#define LOOPFP12 for(j=INFULL2;j<=OUTFULLP12;j++)
#define LOOPFP11 for(i=INFULL1;i<=OUTFULLP11;i++)


// full loop + 1 (shift away from boundary) on inner edge for comptuing emf or corner quantities
#define LOOPINFP13 for(k=INFULLP13;k<=OUTFULL3;k++)
#define LOOPINFP12 for(j=INFULLP12;j<=OUTFULL2;j++)
#define LOOPINFP11 for(i=INFULLP11;i<=OUTFULL1;i++)


#define LOOPOUTFM13 for(k=INFULL3;k<=OUTFULLM13;k++)
#define LOOPOUTFM12 for(j=INFULL2;j<=OUTFULLM12;j++)
#define LOOPOUTFM11 for(i=INFULL1;i<=OUTFULLM11;i++)

#define LOOPH3 for(k=INHALF3;k<=OUTHALF3;k++)
#define LOOPH2 for(j=INHALF2;j<=OUTHALF2;j++)
#define LOOPH1 for(i=INHALF1;i<=OUTHALF1;i++)

#define LOOPP13 for(k=INP13;k<=OUTP13;k++)
#define LOOPP12 for(j=INP12;j<=OUTP12;j++)
#define LOOPP11 for(i=INP11;i<=OUTP11;i++)

#define LOOPN3 for(k=0;k<=N3-1;k++)
#define LOOPN2 for(j=0;j<=N2-1;j++)
#define LOOPN1 for(i=0;i<=N1-1;i++)

#define LOOPFMHP3 for(k=INFULL3;k<=OUTHALF3;k++)
#define LOOPFMHP2 for(j=INFULL2;j<=OUTHALF2;j++)
#define LOOPFMHP1 for(i=INFULL1;i<=OUTHALF1;i++)

#define LOOPHMFP3 for(k=INHALF3;k<=OUTFULL3;k++)
#define LOOPHMFP2 for(j=INHALF2;j<=OUTFULL2;j++)
#define LOOPHMFP1 for(i=INHALF1;i<=OUTFULL1;i++)

#define LOOPHP3 for(k=0;k<=OUTHALF3;k++)
#define LOOPHP2 for(j=0;j<=OUTHALF2;j++)
#define LOOPHP1 for(i=0;i<=OUTHALF1;i++)

// below used for initialization and such, not a computational issue
#define LOOPF LOOPF3 LOOPF2 LOOPF1
#define LOOPH LOOPH3 LOOPH2 LOOPH1
#define LOOPP1 LOOPP13 LOOPP12 LOOPP11
#define LOOP LOOPN3 LOOPN2 LOOPN1
#define LOOPFMHP LOOPFMHP3 LOOPFMHP2 LOOPFMHP1
#define LOOPHMFP LOOPHMFP3 LOOPHMFP2 LOOPHMFP1
#define LOOPHP LOOPHP3 LOOPHP2 LOOPHP1

#define LOOPINT3 for(k=intix3;k<intox3;k++)
#define LOOPINT2 for(j=intix2;j<intox2;j++)
#define LOOPINT1 for(i=intix1;i<intox1;i++)



#define LOOPFC LOOPF
#define LOOPHC LOOPH
#define LOOPFMHPC LOOPFMHP
#define LOOPHMFPC LOOPHMFP
#define LOOPHPC LOOPHP


#define LOOPC3 LOOPN3
#define LOOPC2 LOOPN2
#define LOOPC1 LOOPN1

#define LOOPC LOOPC3 LOOPC2 LOOPC1

// general loop, but assumes i,j,k used
#define ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) \
	for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)\
	for(k=kstart;k<=kstop;k++)

// general loop for any indicies
#define GENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop) \
        for((i)=(istart);(i)<=(istop);(i)++)\
        for((j)=(jstart);(j)<=(jstop);(j)++)\
        for((k)=(kstart);(k)<=(kstop);(k)++)

// general loop for any indicies
#define SUPERGENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk) \
        for((i)=(istart);(i)<=(istop);(i)+=(di))\
        for((j)=(jstart);(j)<=(jstop);(j)+=(dj))\
        for((k)=(kstart);(k)<=(kstop);(k)+=(dk))

/* loop over all active zones */
//#define ZLOOP for(i=0;i<=N1-1;i++)for(j=0;j<=N2-1;j++)for(k=0;k<=N3-1;k++)
#define ZLOOP ZSLOOP(0,N1-1,0,N2-1,0,N3-1)


//#define FULLLOOP ZSLOOP(-N1BND, N1 -1 + N1BND, -N2BND, N2 -1 + N2BND, -N3BND, N3 -1 + N3BND)
#define FULLLOOP LOOPF1 LOOPF2 LOOPF3

// computing emf for FLUXCT
#define LOOPINFP1 LOOPINFP11 LOOPINFP12 LOOPINFP13

#define LOOPINFP1dir1full LOOPF1 LOOPINFP12 LOOPINFP13

#define LOOPINFP1dir2full LOOPINFP11 LOOPF2 LOOPINFP13

#define LOOPINFP1dir3full LOOPINFP11 LOOPINFP12 LOOPF3

#define LOOPINFP1dir23full LOOPINFP11 LOOPF2 LOOPF3

#define LOOPINFP1dir13full LOOPF1 LOOPINFP12 LOOPF3

#define LOOPINFP1dir12full LOOPF1 LOOPF2 LOOPINFP13


// computing emf for FLUXCD
#define LOOPOUTFM1 LOOPOUTFM11 LOOPOUTFM12 LOOPOUTFM13

#define LOOPOUTFM1dir1full LOOPF1 LOOPOUTFM12 LOOPOUTFM13

#define LOOPOUTFM1dir2full LOOPOUTFM11 LOOPF2 LOOPOUTFM13

#define LOOPOUTFM1dir3full LOOPOUTFM11 LOOPOUTFM12 LOOPF3





// larger loop than full for cornered quantities such as emf defined on corners that need to be initialized for boundary condition reasons
#define FULLLOOPP1 LOOPFP11 LOOPFP12 LOOPFP13

//#define WSPEEDLOOP ZSLOOP(-SHIFT1,N1-1+2*SHIFT1,-SHIFT2,N2-1+2*SHIFT2,-SHIFT3,N3-1+2*SHIFT3)


//#define PLUSLOOP ZSLOOP(-1, N1, -1, N2, -1, N3)

// below same as FULLLOOP if NBND=2
//#define PLUSPLUSLOOP ZSLOOP(-2, N1+1, -2, N2+1, -2, N3+1)




// divb loop
//#define LOOPDIVB LOOPP11 LOOPP12 LOOPP13
// boundary zones may not require divb=0 since proxy for flux
#define LOOPDIVB LOOPC1 LOOPC2 LOOPC3
//ZSLOOP(-N1BND+1,N1+1,-1,N2+1)


/* want dump output to be ordered in radius first!! */
#define DUMPLOOP(istart,istop,jstart,jstop,kstart,kstop) \
	for(k=kstart;k<=kstop;k++)\
	for(j=jstart;j<=jstop;j++)\
	for(i=istart;i<=istop;i++)

#if(FULLOUTPUT==0)
#define EXTRADUMP1 0
#define EXTRADUMP2 0
#define EXTRADUMP3 0
#else
#define EXTRADUMP1T FULLOUTPUT*N1NOT1
#define EXTRADUMP2T FULLOUTPUT*N2NOT1
#define EXTRADUMP3T FULLOUTPUT*N3NOT1

#define EXTRADUMP1 ((EXTRADUMP1T>N1BND) ? N1BND : EXTRADUMP1T)
#define EXTRADUMP2 ((EXTRADUMP2T>N2BND) ? N2BND : EXTRADUMP2T)
#define EXTRADUMP3 ((EXTRADUMP3T>N3BND) ? N3BND : EXTRADUMP3T)

#endif

#if(FULLOUTPUT==0)
#define DUMPGENLOOP DUMPLOOP(0,N1-1,0,N2-1,0,N3-1)
#else
#define DUMPGENLOOP DUMPLOOP(-EXTRADUMP1,N1-1+EXTRADUMP1,-EXTRADUMP2,N2-1+EXTRADUMP2,-EXTRADUMP3,N3-1+EXTRADUMP3)
#endif


// defines whether within the enerregion
#define WITHINENERREGION(i,j,k) (i>=enerpos[X1DN])&&(i<=enerpos[X1UP])&&(j>=enerpos[X2DN])&&(j<=enerpos[X2UP])&&(k>=enerpos[X3DN])&&(k<=enerpos[X3UP]) 


#define NUMIMAGEPARMS 2

#define ORIGIN 0
#define LMIN 1



  //#define IMAGELOOP(istart,istop,jstart,jstop,kstart,kstop) \
  //	for(k=kstart;k<=kstop;k++)\
  //	for(j=jstart;j<=jstop;j++)\
  //	for(i=istart;i<=istop;i++)

//#define OLDIMAGELOOP for(j=N2-1;j>=0;j--) for(i=0;i<N1;i++)	// nasty 
								// to
								// deal 
								// with

  // check for existence in bad form using:
// grep "PLOOP" *.c | grep --invert-match "PLOOP("
/* loop over all Primitive variables */
#define PLOOP(pl) for(pl=0;pl<NPR;pl++)
/* loop over all dumped Primitive variables */
#define PDUMPLOOP(pd) for(pd=0;pd<NPRDUMP;pd++)
/* loop over all inversion Primitive variables */
#define PINVERTLOOP(pi) for(pi=0;pi<NPRINVERT;pi++)
/* loop over all bounding Primitive variables */
#define PBOUNDLOOP(pb) for(pb=0;pb<NPRBOUND;pb++)
/* loop over all center to edge variables */
#define PLOOPINTERP(pl) for(pl=0;pl<NPR2INTERP;pl++)


/* loop over all Dimensions; second rank loop */
#define DLOOP(j,k) for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA(j) for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP(j,k) for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA(j) for(j=1;j<NDIM;j++)
/* loop over all for j and Space for k; second rank loop */
#define DSLOOP(j,k) for(j=0;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* space-space */
#define SSLOOP(j,k) for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all for k and Space for j; second rank loop */
#define SDLOOP(j,k) for(j=1;j<NDIM;j++)for(k=0;k<NDIM;k++)

// loop over directions
#define DIRLOOP(dir) for(dir=0;dir<COMPDIM*2;dir++)

#define DIMENLOOP(dir) for(dir=1;dir<=COMPDIM;dir++)

// loop over fail flag in boundary code
#define FLOOP(ff) for(ff=FLAGUTOPRIMFAIL;ff<=FLAGUTOPRIMFAIL;ff++)

// loop over jet regions
#define JETLOOP(jetio) for(jetio=0;jetio<NUMJETS;jetio++)

// loop over ener/flux regions
#define ENERREGIONLOOP(enerregion) for(enerregion=0;enerregion<NUMENERREGIONS;enerregion++)

// loop over fair/floor types
#define FLOORLOOP(floor) for(floor=0;floor<NUMFAILFLOORFLAGS;floor++)

// loop over debug time scales
#define TSCALELOOP(tscale) for(tscale=0;tscale<NUMTSCALES;tscale++)


// loop over sources
#define SCLOOP(sc) for(sc=0;sc<NUMSOURCES;sc++)

// loop over fluxterms
#define FLLOOP(fl) for(fl=0;fl<NUMFLUXTERMS;fl++)


// loop over pflag flags
#define PFLAGLOOP(pf) for(pf=0;pf<NUMPFLAGS;pf++)

// for USEMPI&&USEROMIO==1
#define ROMIOCOLLOOP(romiocoliter) for(romiocoliter=0;romiocoliter<romiocloopend;romiocoliter++)

#define BUFFERINIT nextbuf=0
#define COLINIT nextcol=0

// for mpicombie==0
#define COLLOOP(coliter) for(coliter=0;coliter<numfiles;coliter++)


#define DTSTAGELOOP(dtstage) for(dtstage=0;dtstage<MAXDTSTAGES;dtstage++)

#define INTERPLOOP(interpi) for(interpi=0;interpi<NUMINTERPTYPES;interpi++)

#define ENODEBUGLOOP(enodebugi) for(enodebugi=0;enodebugi<NUMENODEBUGS;enodebugi++)


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Various macros
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if(MCOORD!=CARTMINKMETRIC)

// e.g. MYGDET(i,j,k,CENT)
#define MYGDET(i,j,k,p) (gdet[i][j][k][p])

#else

#define MYGDET(i,j,k,p) (gdet[0][0][0][p])

#endif



// for FLUXCT ala Toth 2000

// below 4 macros used for vpot2field() in initbase.c
#define FgCORN(F,i,j,k) (F[i][j][k])
// AVG_x functions that average in x direction
// typically operates on corner quantities to get center quantities
#define AVGCORN_1(F,i,j,k) (0.5*(FgCORN(F,ip1mac(i),j,k) + FgCORN(F,i,j,k) ))
#define AVGCORN_2(F,i,j,k) (0.5*(FgCORN(F,i,jp1mac(j),k) + FgCORN(F,i,j,k) ))
#define AVGCORN_3(F,i,j,k) (0.5*(FgCORN(F,i,j,kp1mac(k)) + FgCORN(F,i,j,k) ))

// below macros are similar to those farther below for divb, but no gdet stuff
#define FgN(F,i,j,k,pl) (F[i][j][k][pl])

// AVG_x functions that average in x direction
// typically operates on centered quantities to get edge/corner quantities
#define AVGN_1(F,i,j,k,pl) (0.5*(FgN(F,i,j,k,pl) + FgN(F,im1mac(i),j,k,pl) ))
#define AVGN_2(F,i,j,k,pl) (0.5*(FgN(F,i,j,k,pl) + FgN(F,i,jm1mac(j),k,pl) ))
#define AVGN_3(F,i,j,k,pl) (0.5*(FgN(F,i,j,k,pl) + FgN(F,i,j,km1mac(k),pl) ))

// just average of 4 quantities that eventually lives on (0.5,0,0) (i.e. CORN1)
#define AVGN_for1(F,i,j,k,pl) (0.5*(AVGN_2(F,i,j,k,pl)+AVGN_2(F,i,j,km1mac(k),pl))) // 2 and 3 directions averaged
// lives on (0,0.5,0) (i.e. CORN2)
#define AVGN_for2(F,i,j,k,pl) (0.5*(AVGN_1(F,i,j,k,pl)+AVGN_1(F,i,j,km1mac(k),pl))) // 1 and 3 directions averaged
// lives on (0,0,0.5) (i.e. CORN3)
#define AVGN_for3(F,i,j,k,pl) (0.5*(AVGN_2(F,i,j,k,pl)+AVGN_2(F,im1mac(i),j,k,pl))) // 1 and 2 directions averaged



// below macros used for divb definition
#define Fg(F,i,j,k,pl) (F[i][j][k][pl]*MYGDET(i,j,k,CENT))

// AVG_x functions that average in x direction
// typically operates on centered quantities to get edge/corner quantities
#define AVG_1(F,i,j,k,pl) (0.5*(Fg(F,i,j,k,pl) + Fg(F,im1mac(i),j,k,pl) ))
#define AVG_2(F,i,j,k,pl) (0.5*(Fg(F,i,j,k,pl) + Fg(F,i,jm1mac(j),k,pl) ))
#define AVG_3(F,i,j,k,pl) (0.5*(Fg(F,i,j,k,pl) + Fg(F,i,j,km1mac(k),pl) ))

// just average of 4 quantities that eventually lives on (0.5,0,0) (i.e. CORN1)
#define AVG_for1(F,i,j,k,pl) (0.5*(AVG_2(F,i,j,k,pl)+AVG_2(F,i,j,km1mac(k),pl))) // 2 and 3 directions averaged
// lives on (0,0.5,0) (i.e. CORN2)
#define AVG_for2(F,i,j,k,pl) (0.5*(AVG_1(F,i,j,k,pl)+AVG_1(F,i,j,km1mac(k),pl))) // 1 and 3 directions averaged
// lives on (0,0,0.5) (i.e. CORN3)
#define AVG_for3(F,i,j,k,pl) (0.5*(AVG_2(F,i,j,k,pl)+AVG_2(F,im1mac(i),j,k,pl))) // 1 and 2 directions averaged


// #define DIVBCONDITION(p,i,j)
// if((i>=-1)&&(j>=-1)&&(startpos[2]+j!=0)&&(startpos[2]+j!=N2TOT))
//#define DIVBCONDITION(p,i,j) if((startpos[1]+i>0)&&(startpos[2]+j>0)&&(startpos[1]+i<totalsize[1])&&(startpos[2]+j<totalsize[2]))
//#define DIVBCONDITION(p,i,j) if((startpos[1]+i!=0)&&(startpos[1]+i!=totalsize[1])&&(startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2]))

// only needed for polar axis condition
//GODMARK
//#define DIVBCONDITION(p,i,j,k) if((startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2]))

// GODMARK: No, need for all boundaries
//#define DIVBCONDITION(p,i,j,k) if(((startpos[1]+i!=0)&&(startpos[1]+i!=totalsize[1]) || (!N1NOT1))&&((startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2]) || (!N2NOT1))&&((startpos[3]+k!=0)&&(startpos[3]+k!=totalsize[3]) || (!N3NOT1)))

// only consider within bounds or ignore condition if no such dimension
#define DIVBCDIR1 ((i >= -N1BND+1 && i <= N1 - 1 + N1BND - 1) || (!N1NOT1))
#define DIVBCDIR2 ((j >= -N2BND+1 && j <= N2 - 1 + N2BND - 1) || (!N2NOT1))
#define DIVBCDIR3 ((k >= -N3BND+1 && k <= N3 - 1 + N3BND - 1) || (!N3NOT1))

#define DIVBCONDITION(p,i,j,k) if(DIVBCDIR1&&DIVBCDIR2&&DIVBCDIR3)

//#define DIVBCONDITION(p,i,j,k) if(1)




//(startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2])&&(startpos[3]+k!=0)&&(startpos[3]+k!=totalsize[3])

// the below lines results in quantity at TRUE CORNER
#define DIVBDIFFFLUXCTx(p,i,j,k) (AVG_for1(p,i,j,k,B1)-AVG_for1(p,im1mac(i),j,k,B1))

#define DIVBDIFFFLUXCTy(p,i,j,k) (AVG_for2(p,i,j,k,B2)-AVG_for2(p,i,jm1mac(j),k,B2))

#define DIVBDIFFFLUXCTz(p,i,j,k) (AVG_for3(p,i,j,k,B3)-AVG_for3(p,i,j,km1mac(k),B3))

#define DIVBNORMFLUXCTx(p,i,j,k) (AVG_for1(p,i,j,k,B1)+AVG_for1(p,im1mac(i),j,k,B1))

#define DIVBNORMFLUXCTy(p,i,j,k) (AVG_for2(p,i,j,k,B2)+AVG_for2(p,i,jm1mac(j),k,B2))

#define DIVBNORMFLUXCTz(p,i,j,k) (AVG_for3(p,i,j,k,B3)+AVG_for3(p,i,j,km1mac(k),B3))

#define DIVBNORMFLUXCT(p,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(THIRD*MYGDET(i,j,k,CORNT)*fabs(DIVBNORMFLUXCTx(p,i,j,k)+DIVBNORMFLUXCTy(p,i,j,k)+DIVBNORMFLUXCTz(p,i,j,k)) +SMALL))

#define DIVBFLUXCT(p,i,j,k) ((\
                             DIVBDIFFFLUXCTx(p,i,j,k)/dx[1] + DIVBDIFFFLUXCTy(p,i,j,k)/dx[2] + DIVBDIFFFLUXCTz(p,i,j,k)/dx[3]\
                             )*(DIVBNORMFLUXCT(p,i,j,k)))



// FLUXCD DIVB

// FLUXCD divb lives at CENT
#define DIVBNORMFLUXCD(p,i,j,k) (MAX(MAX(dx[1],dx[2]),dx[3])/(THIRD*MYGDET(i,j,k,CENT)*fabs(\
                                   Fg(p,ip1mac(i),j,k,B1) + Fg(p,im1mac(i),j,k,B1)\
                                  +Fg(p,i,jp1mac(j),k,B2) + Fg(p,i,jm1mac(j),k,B2)\
                                  +Fg(p,i,j,kp1mac(k),B3) + Fg(p,i,j,km1mac(k),B3)\
				 )+SMALL))

#define DIVBFLUXCD(p,i,j,k)  (0.5*(\
                                  (Fg(p,ip1mac(i),j,k,B1) - Fg(p,im1mac(i),j,k,B1))/dx[1]\
                                 +(Fg(p,i,jp1mac(j),k,B2) - Fg(p,i,jm1mac(j),k,B2))/dx[2]\
                                 +(Fg(p,i,j,kp1mac(k),B3) - Fg(p,i,j,km1mac(k),B3))/dx[3]\
                                )*DIVBNORMFLUXCD(p,i,j,k))

// poles defined as divb=0, can't divide due to singularity (could use
// volume regularization)
#define SETFDIVBFLUXCT(divb,p,i,j,k) {DIVBCONDITION(p,i,j,k){ divb = fabs(DIVBFLUXCT(p,i,j,k)) ;} else divb = 0.;}

#define SETFDIVBFLUXCD(divb,p,i,j,k) {DIVBCONDITION(p,i,j,k){ divb = fabs(DIVBFLUXCD(p,i,j,k)) ;} else divb = 0.;}

// for now // GODMARK
#define SETFDIVB(divb,p,i,j,k) SETFDIVBFLUXCT(divb,p,i,j,k)




#define MYDMIN(a,b) (mydminarg1=(a),mydminarg2=(b),(mydminarg1) < (mydminarg2) ?\
        (mydminarg1) : (mydminarg2))

#define delta(i,j) ((i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

#define mink(I,J) (I != J ? (0.) : (I == 0 ? (-1.) : (1.)))

#define pfixupeach(pr,i,j,k,which,min) {if(pr[which]<min){ fladd[which]+=dV*MYGDET(i,j,k,CENT)*(min-pr[which]); pr[which]=min;}}

#define pfixup(pr,i,j,k) {pfixupeach(pr,i,j,k,RHO,RHOMIN); pfixupeach(pr,i,j,k,UU,UUMIN); }

// #define FAILSTATEMENT(file,function,number) {fprintf(fail_file,"%s
// %d-%s(): failure\n",file,number,function); fflush(fail_file);
// fprintf(fail_file,"rho[i][j][k]: %21.15g uu[i][j][k]: %21.15g rho2[i][j][k]:
// %21.15g uu2[i][j][k]: %21.15g i: %d j: %d pl:
// %d\n",p[i][j][k][RHO],p[i][j][k][UU],ph[i][j][k][RHO],ph[i][j][k][UU],i,j,pl);
// return(1);}

#define FAILSTATEMENT(file,function,number) {if(debugfail>=1){ dualfprintf(fail_file,"%s %d-%s(): failure\n",file,number,function); dualfprintf(fail_file,"i: %d j: %d k: %d p: %d\n",icurr,jcurr,kcurr,pcurr);} return(1);}

#define FAILSTATEMENTVOID(file,function,number) {if(debugfail>=1){ dualfprintf(fail_file,"%s %d-%s(): failure\n",file,number,function); dualfprintf(fail_file,"i: %d j: %d k: %d p: %d\n",icurr,jcurr,kcurr,pcurr);} }

#if(JONCHECKS2)
#define MYFUN(fun,one,two,three) if(fun>=1){ FAILSTATEMENT(one,two,three);}
#else
#define MYFUN(fun,one,two,three) {fun;}
#endif



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// structure definitions
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct blink {
  int num;
  struct blink * np;
  // only used by cpu=0
  int cpu; // which cpu
  int i,j,k,col; // starting values for cpu=0
  int ri,rj,rk,rcol; // reference values for first cpu in sequence of nodes for a single buffer
  int end;
};



// structure declarations
/* set global variables that indicate current local metric, etc. */
struct of_geom {
  //FTYPE gcon[NDIM][NDIM];
  //  FTYPE gcov[NDIM][NDIM];
  // bit faster since not all values always used
  FTYPE (*gcov)[NDIM];
  FTYPE (*gcon)[NDIM];
  FTYPE g;
  FTYPE *gcovpert;
  FTYPE e[NPR]; // eomfunc
  int i,j,k,p;
};


struct of_state {
  FTYPE rho;
  FTYPE ie;
  FTYPE ucon[NDIM];
  FTYPE ucov[NDIM];
  FTYPE bcon[NDIM];
  FTYPE bcov[NDIM];
};



// now that all hashes have been defined, get mpi header
#include "mympi.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// function declarations
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern int main(int argc, char *argv[]);
extern int error_check(int wherefrom);
extern int find_horizon(void);

// initialize DUMP stuff
extern int init_dumps(void);
int setuplinklist(int numcolumns,int which);
extern struct blink * addlink(struct blink * clinkptr);

// ENER file stuff
extern int dump_ener(int doener, int dordump, int call_code);

extern int diag(int call_code);

extern void frdotout(void);

extern void makedirs(void);

extern void appendener(FILE* ener_file,SFTYPE pdot_tot[][NPR],SFTYPE*fladd_tot,SFTYPE*sourceadd_tot);

extern void divbmaxavg(FTYPE p[][N2M][N3M][NPR],SFTYPE*ptrdivbmax,SFTYPE*ptrdivbavg);
extern void gettotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[]);
extern void gettotali(int numvars, int* vars[],int*sizes,int*vars_tot[]);
extern int constotal(int enerregion, SFTYPE *vars_tot);
extern int integrate(SFTYPE * var,SFTYPE *var_tot,int type, int enerregion);

extern int counttotal(int tscale, int enerregion, CTYPE *vars_tot, int num);
extern int integratel(CTYPE * var,CTYPE *var_tot,int type, int tscale, int enerregion);


// DUMP file stuff
extern int isenoughfreespace(unsigned long long need);

extern int dump_gen(int readwrite, long dump_cnt, int bintxt, int whichdump,MPI_Datatype datatype, char *fileprefix, char *fileformat, char *filesuffix, int (*headerfun) (int bintxt, FILE*headerptr),int (*content) (int i, int j, int k, MPI_Datatype datatype, void*setbuf));

extern int dump(long dump_cnt);
extern int dump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int dump_header(int bintxt, FILE *headerptr);

extern int avgdump(long avg_cnt);
extern int avg_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int avgdump2(long avg_cnt);
extern int avg2_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int debugdump(long debug_cnt);
extern int debug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int gdump(void);
extern int gdump_content(int i, int j, int k, MPI_Datatype datatype, void *writebuf);

extern int fieldlinedump(long fieldline_cnt);
extern int fieldline_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int dissdump(long dump_cnt);
extern int dissdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);


extern int image_dump(long image_cnt);
extern int imagedefs(int whichk, int scale, int limits, int vartype);
extern int image(long dump_cnt, int whichk, int scale, int limits, int vartype);
extern int image_header(int bintxt, FILE *headerptr);
extern int image_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern void prminmaxsum(FTYPE p[][N2M][N3M][NPR], int start,int nmemb, FTYPE *max, FTYPE*min,FTYPE*sum);

extern int restart_init(int which);

extern int restart_read(long which);
extern int read_restart_header(int bintxt, FILE* headerptr);
extern int restart_read_defs(void);
extern int rdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int restart_write(long dump_cnt);
extern int write_restart_header(int bintxt, FILE* headerptr);
extern int rdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern void myfopen(char*fname, char*fmt, char*message, FILE ** fileptr);
extern void myfclose(FILE ** fileptr,char*message);

extern void myset(MPI_Datatype datatype, void *ptr, int start, int nmemb, void*writebuf);
extern void myget(MPI_Datatype datatype, void *ptr, int start, int nmemb, void*writebuf);

extern void myfwrite(int bintxt, MPI_Datatype datatype, void *ptr, int start, int nmemb, int i, int j, int k, FILE**stream,void*writebuf);

extern void myfread(int bintxt, MPI_Datatype datatype, void *ptr, int start, int nmemb, int i, int j, int k, FILE**stream,void*writebuf);



// initialize stuff
// specific to init.c's and used in initbase.c and init.c, so leave global
extern int post_init_specific_init(void);
extern int pre_init_specific_init(void);
extern int init_grid(void);
extern int init_global(void);
extern int vpot2field(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE p[][N2M][N3M][NPR]);
// called in restart.c and initbase.c
extern void set_grid(void);


// some physics

extern int sourcephysics(FTYPE *ph, struct of_geom *geom,
		       struct of_state *q,FTYPE (*dUcomp)[NPR]);

extern void postdt(void);
extern int primtoU(int returntype, FTYPE *p, struct of_state *q, struct of_geom *geom,
		   FTYPE *U);

extern int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
extern int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
extern int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
#if(RELTYPE==RELEOM)

#if(WHICHVEL==VEL4)
#define ucon_calc ucon_calc_4vel
#define dudp_calc dudp_calc_gen
#elif(WHICHVEL==VEL3)
#define ucon_calc ucon_calc_3vel
#define dudp_calc dudp_calc_3vel
#elif(WHICHVEL==VELREL4)
#define ucon_calc ucon_calc_rel4vel
#define dudp_calc dudp_calc_gen

#elif(RELTYPE==NONRELEOM) // not really right
#define ucon_calc ucon_calc_nonrel
#define dudp_calc dudp_calc_nonrel
#endif

#endif

extern int ucon_calcother(FTYPE *pr, FTYPE *ucon);
extern void ucon_precalc(FTYPE *ucon, FTYPE *AA, FTYPE *BB,
			 FTYPE *CC, FTYPE *discr);



extern FTYPE ranc(int seed);

// fixup stuff

extern int check_pr(FTYPE *pr, FTYPE *prmodel, struct of_geom *geom, int modelpos,int finalstep);
extern int ucon_fix(FTYPE disc, FTYPE AA, FTYPE BB, FTYPE CC,
		    FTYPE *ucon);

/* // dudp stuff */

/* extern void dutdui_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void duiduj_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui); */
/* extern void dbiduj_calc(FTYPE *dbtdui, FTYPE *dutdui, FTYPE *ucon, */
/* 			FTYPE *b, FTYPE dbiduj[][NDIM]); */
/* extern void db2dui_calc(FTYPE dbiduj[][NDIM], FTYPE *b, */
/* 			FTYPE *db2dui); */
/* extern void duudud_calc(FTYPE *ucon, FTYPE duudud[][NDIM]); */

/* extern void dbsqdui_calc(FTYPE dbiduj[][NDIM], FTYPE *b, */
/* 			 FTYPE *dbsqdui); */
/* extern void dgdvi_calc(FTYPE *pr,FTYPE *dgdvi); */
/* extern void duidvj_calc(FTYPE *dgdv,FTYPE duidvj[][NDIM]); */
/* extern void dudduu_calc(FTYPE*dutdui, FTYPE dudduu[][NDIM]); */
/* extern void dbdiduj_calc(FTYPE dbiduj[][NDIM],FTYPE dbdiduj[][NDIM]); */
/* extern void ducon_dv3_calc(struct of_state *q,FTYPE ducon_dv[][NDIM]); */
extern int sp_stress_calc(FTYPE *pr, FTYPE tens_matt[][NDIM],
			  FTYPE tens_em[][NDIM], FTYPE *b,
			  FTYPE *ucon);



// log file stuff
extern void myfprintf(FILE* fileptr, char *format, ...);
extern void dualfprintf(FILE* fileptr,char *format, ...);
extern void logsfprintf(char *format, ...);
extern void trifprintf(char *format, ...);

// boundary stuff
extern int bound_prim(int boundstage, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_pflag(int boundstage, int primbase[][N2M][N3M][NUMPFLAGS]);
extern int inflow_check_4vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_3vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_rel4vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);


// transform stuff
extern int bl2met2metp2v(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk);
extern int metp2met2bl(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk);
extern int pr2ucon(int whichvel, FTYPE *pr, struct of_geom *geom, FTYPE*ucon);
extern int coordtrans(int whichcoordin, int whichcoordout, int ii, int jj, int kk, FTYPE*ucon);
extern void bltoks(int ii, int jj, int kk, FTYPE*ucon);
extern void kstobl(int ii, int jj, int kk, FTYPE*ucon);
extern void mettometp(int ii, int jj, int kk, FTYPE*ucon);
extern void metptomet(int ii, int jj, int kk, FTYPE*ucon);
extern void ucon2pr(int whichvel, FTYPE *ucon, struct of_geom *geom, FTYPE *pr);
extern int vcon2pr(int whichvel, FTYPE *vcon, struct of_geom *geom, FTYPE *pr);



// metric stuff
extern void gset(int getprim, int whichcoord, int i, int j, int k, struct of_geom *geom);
extern FTYPE gdet_func(FTYPE gcov[][NDIM]);
//extern FTYPE bl_gdet_func(FTYPE r, FTYPE th);
//extern void bl_gcov_func(FTYPE r, FTYPE th, FTYPE gcov[][NDIM]);
//extern void bl_gcon_func(FTYPE r, FTYPE th, FTYPE gcon[][NDIM]);
extern void conn_func(int whichcoord, FTYPE *X, struct of_geom *geom,
		      FTYPE lconn[][NDIM][NDIM],FTYPE *conn2);
extern void mks_unitheta_idxvol_func(int i, int j, int k, FTYPE *idxvol);

//extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM]);
//extern void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
extern void matrix_inverse(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);

// coordinate stuff
extern void set_coord_parms(void);
extern void write_coord_parms(void);
extern void read_coord_parms(void);
extern void coord(int i, int j, int k, int loc, FTYPE *X);
extern void bl_coord(FTYPE *X, FTYPE *V);
extern int setihor(void);
extern FTYPE setRin(int ihor);



// eos stuff
extern int pickeos_eomtype(int whicheom);

extern FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
extern FTYPE u_rho0_p(FTYPE rho0, FTYPE p);
extern FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);
extern FTYPE dpdu_rho0_u(FTYPE rho0, FTYPE u);
extern FTYPE dpdrho0_rho0_u(FTYPE rho0, FTYPE u);
extern FTYPE cs2_compute(FTYPE rho0, FTYPE u);
extern FTYPE compute_entropy(FTYPE rho0, FTYPE u);
extern FTYPE compute_u_from_entropy(FTYPE rho0, FTYPE entropy);

// eos stuff used by inversion
extern FTYPE pressure_wmrho0_grmhd(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idwmrho0dp_grmhd(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idrho0dp_grmhd(FTYPE rho0, FTYPE wmrho0);

extern FTYPE pressure_wmrho0_coldgrmhd(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idwmrho0dp_coldgrmhd(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idrho0dp_coldgrmhd(FTYPE rho0, FTYPE wmrho0);

extern FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idwmrho0dp(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idrho0dp(FTYPE rho0, FTYPE wmrho0);



// physics stuff
extern FTYPE contract(FTYPE *vcon, FTYPE *wcon);

extern int bsq_calc(FTYPE *pr, struct of_geom *geom, FTYPE *b2);
extern void b_calc(FTYPE *pr, FTYPE *ucon, FTYPE *b);


extern int gamma_calc(FTYPE *pr, struct of_geom *geom,FTYPE *gamma);


extern int dudp_calc_gen(int whichcons,FTYPE *pr, struct of_state *q,
		     struct of_geom *geom, FTYPE **alpha);

extern int dudp_calc_3vel(int whichcons,FTYPE *pr, struct of_state *q,
		     struct of_geom *geom, FTYPE **alpha);



extern int Utoprimgen(int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr);
extern int Utoprimloop(FTYPE unew[][N2M][N3M][NPR],FTYPE pf[][N2M][N3M][NPR]);
extern int primtoUloop(FTYPE pi[][N2M][N3M][NPR],FTYPE unew[][N2M][N3M][NPR]);

extern int Utoprim(int entropyeom, FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);
extern int Utoprim_ldz(FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);
extern int Utoprim_1d(FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);
extern int Utoprim_1d_opt(FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);
extern int Utoprim_2d(FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);
extern int Utoprim_1d_final(FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);
extern int Utoprim_2d_final(FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);
//extern int Utoprim_2d_final_nonrelcompat_inputnorestmass(FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);  //wrong function name, corrected by atch, see below
extern int Utoprim_jon_nonrelcompat_inputnorestmass(int whicheom, FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);
extern int Utoprim_5d2_final(FTYPE *U, struct of_geom *geom, int *lpflag, FTYPE *pr);




extern void tetr_func(FTYPE tetr_cov[][NDIM], FTYPE tetr_con[][NDIM]);
extern void get_geometry(int i, int j, int k, int loc, struct of_geom *geom);
extern int get_state(FTYPE *pr, struct of_geom *geom,
		     struct of_state *q);
extern int primtoflux(int returntype, FTYPE *pa, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *fl);
extern void mks_source_conn(FTYPE *ph, struct of_geom *ptrgeom,
		     struct of_state *q,FTYPE *dU);
extern int source(FTYPE *pa, struct of_geom *geom,
		  FTYPE (*Uacomp)[NPR], FTYPE *Ua);

extern FTYPE taper_func(FTYPE R,FTYPE rin) ;
extern FTYPE rhor_calc(int which);
extern FTYPE rmso_calc(int which) ;
extern FTYPE uphi_isco_calc(int which,FTYPE r);

extern int set_atmosphere(int whichcond, int whichvel, struct of_geom *geom, FTYPE *pr);
extern int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *scaler);
extern int fixup(int stage,FTYPE (*var)[N2M][N3M][NPR],int finalstep);
extern int fixup1zone(FTYPE *pr,struct of_geom *ptrlgeom, int finalstep);
extern int diag_fixup(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int finalstep,int whocalled);

extern int superdebug(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled);

extern int limit_gamma(FTYPE gammamax, FTYPE*pr,struct of_geom *geom, int finalstep);

extern int fixup_checksolution(int stage, FTYPE (*pv)[N2M][N3M][NPR],int finalstep);
extern int fixup_utoprim(int stage, FTYPE (*pv)[N2M][N3M][NPR],FTYPE (*pbackup)[N2M][N3M][NPR], int finalstep);
extern int post_fixup(int stage,FTYPE (*pv)[N2M][N3M][NPR],FTYPE (*pbackup)[N2M][N3M][NPR],int finalstep);
extern int pre_fixup(int stage, FTYPE (*pv)[N2M][N3M][NPR]);
extern int get_bsqflags(int stage, FTYPE (*pv)[N2M][N3M][NPR]);


extern void tet_func(FTYPE metr[][NDIM], FTYPE tetr[][NDIM]);
extern int dsyev_(char *jobz, char *uplo, int *n, FTYPE *a, int *lda,
		  FTYPE *w, FTYPE *work, int *lwork, int *iwork,
		  int *liwork, int *info);
extern int fail(int fail_type);
extern void setfailresponse(int restartonfail);
extern int setflux(void);
extern int setjetflux(void);
extern void setrestart(int*appendold);


extern int vchar(FTYPE *pr, struct of_state *q, int dir,
		 struct of_geom *geom, FTYPE *cmax, FTYPE *cmin,int *ignorecourant);
extern FTYPE chk_disp(FTYPE v);
extern void make_co_to_comov(FTYPE *ucon, FTYPE ecov[][NDIM],
			     FTYPE econ[][NDIM]);
extern void transform(FTYPE *vec, FTYPE t[][NDIM]);
extern void coeff_set(FTYPE rho, FTYPE u);
extern void transform(FTYPE *ucon, FTYPE t[][NDIM]);

extern void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);

extern int area_map(int call_code, int type, int size, int i, int j, int k, FTYPE prim[][N2M][N3M][NPR]);
extern void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov,
		      FTYPE *bcon);

extern FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
extern void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE faraday[][NDIM]);
extern void current_precalc(int which, struct of_geom *geom, struct of_state *q, SFTYPE Dt,FTYPE faraday[][3]);
extern void init_varstavg(void);
extern void final_varstavg(FTYPE IDT);
extern int set_varstavg(FTYPE tfrac);
extern void current_calc(FTYPE cfaraday[][N2M][N3M][NUMCURRENTSLOTS][3]);
extern int current_doprecalc(int which, FTYPE p[][N2M][N3M][NPR]);
extern int average_calc(int doavg);

/////////////////////////////////////
//
// NR STUFF
//
/////////////////////////////////////
extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
extern void lubksb(FTYPE **a, int n, int *indx, FTYPE *d);
//extern FTYPE zbrent(FTYPE (*func) (FTYPE), FTYPE v1, FTYPE v2,
//		     FTYPE tol);


/* NR routines from nrutil.h */
extern int *ivector(long nl, long nh);
extern void free_ivector(int *v, long nl, long nh);
extern FTYPE *dvector(long nl, long nh);
extern void free_dvector(FTYPE *v, long nl, long nh);
extern FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch);
extern void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl,
			 long nch);
extern FTYPE ***dtensor(long nrl, long nrh, long ncl, long nch,
			 long ndl, long ndh);
extern void free_dtensor(FTYPE ***t, long nrl, long nrh, long ncl,
			 long nch, long ndl, long ndh);
extern void nrerror(char error_text[]);

//////////////////////////////
//
// specialty functions
//
//////////////////////////////
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs, FTYPE *Urs,
			FTYPE *Edot);
extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf, FTYPE r,
			  FTYPE rs, FTYPE urs);
extern void timestep(FTYPE ndtr, FTYPE ndth);
extern FTYPE dtset(FTYPE ndtr, FTYPE ndth);

extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf,
			  FTYPE r, FTYPE rs, FTYPE urs);
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs,
			FTYPE *Urs, FTYPE *Edot);
extern FTYPE edot_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedr_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edr2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edur2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edrdur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);

extern void lower_vec(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void lowerf(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void raise_vec(FTYPE *v1, struct of_geom *geom, FTYPE *v2);
extern void gaussj(FTYPE **tmp, int n, FTYPE **b, int m);
extern void set_points(void);
// extern FTYPE delta(int j, int k) ;
// extern FTYPE mink(int j, int k) ;
extern void make_tetr(FTYPE *ucon, FTYPE econ[][NDIM]);
extern void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);


extern FTYPE sign(FTYPE a);
extern FTYPE sign2(FTYPE a);

#ifndef WIN32
extern FTYPE max(FTYPE a, FTYPE b);

extern FTYPE min(FTYPE a, FTYPE b);
#endif
///////////////////////////////////
//
// SUPERLONGDOUBLE declarations
//
///////////////////////////////////

#if(SUPERLONGDOUBLE)
#include "mconf.h"
extern long double fabsl ( long double );
extern long double sqrtl ( long double );
extern long double cbrtl ( long double );
extern long double expl ( long double );
extern long double logl ( long double );
extern long double tanl ( long double );
extern long double atanl ( long double );
extern long double sinl ( long double );
extern long double asinl ( long double );
extern long double cosl ( long double );
extern long double acosl ( long double );
extern long double powl ( long double, long double );
extern long double tanhl ( long double );
extern long double atanhl ( long double );
extern long double sinhl ( long double );
extern long double asinhl ( long double );
extern long double coshl ( long double );
extern long double acoshl ( long double );
extern long double exp2l ( long double );
extern long double log2l ( long double );
extern long double exp10l ( long double );
extern long double log10l ( long double );
extern long double gammal ( long double );
extern long double lgaml ( long double );
extern long double jnl ( int, long double );
extern long double ynl ( int, long double );
extern long double ndtrl ( long double );
extern long double ndtril ( long double );
extern long double stdtrl ( int, long double );
extern long double stdtril ( int, long double );
extern long double ellpel ( long double );
extern long double ellpkl ( long double );
long double lgammal(long double);
extern int isfinitel ( long double );
#define finite(arg) isfinitel(arg)
extern int merror;
#else



#include <math.h>

#ifdef WIN32
#define finite(arg) _finite(arg)
#define isfinite(arg) _finite(arg)
#else
//#define isfinite(arg) finite(arg)  //atch -- on mako, in force-free it would complain about multiply-defined __finite() if not include this line
#endif

#endif


#ifdef WIN32
#define copysign( a, b ) ( fabs(a) * sign(b) ) 
#endif

#if( DO_WENO_DEBUG )
///debug atch SASMARK
FILE *debugfp;
#endif //DO_WENO_DEBUG

//defines of the interpolation direction unit vector -- can be used outside of ENODEBUG stuff, so put it here
#define ITERGLOBALDEF {di = (iterglobal==1); dj = (iterglobal==2); dk = (iterglobal==3);}

#if(DOENODEBUG) //atch 
#define ENODEBUGPARAM_SMONO 0
#define ENODEBUGPARAM_WENO5 1
#define ENODEBUGPARAM_WENO3 2
#define ENODEBUGPARAM_dPP 3
#define ENODEBUGPARAM_LIMCORR 4
#define ENODEBUGPARAM_LIMCORR_PRIM 5
#endif 


#endif// endif for #ifndef GLOBAL_H



