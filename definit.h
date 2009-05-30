// whether to check symmetry along diagonal (for test==151 in Sasha's tests)
#define ASYMDIAGCHECK 0

// whether to use accurate (bit more expensive) diagnostics.  Only necessary for timeorder>2
#define ACCURATEDIAG 1

#define ACCURATESINCOS 0 // whether to use sin/cos exactly symmetrized around pi/2

// whether to try cold inversion if hot fails
#define HOT2COLD 0


// *NUMBER* OF DIMENSIONS FOR COMPUTATION
// Can choose 3 and still do 1D optimized problems
// Choosing 2 or even 1 reduces some memory load for those things that are accessed by using 1,2,3 instead of arbitrarily accessed
// not sure if 1 works.
#define COMPDIM 3


// size of global arrays
// ALSO CONTROLS WHETHER MEMORY IS ALLOCATED FOR SPECIAL 1D/2D/3D directional arrays
// e.g. if N3==1, then second dimension is neglected in all calculations

#define N1	64		// number of zones in 1-direction
#define N2	64       	// number of zones in 2-direction
#define N3      1               // number of zones in 3-direction

//#if(DOENOFLUX==ENOFINITEVOLUME)
#define MAXBND 5
//#else
//#define MAXBND 3
//#endif


// choose metric
#define MCOORD KSCOORDS


#define DO_VORTICITY_IMAGE 0

// define your user type here
// (normal non-sensitive or performance critical datatypes)
#define REALTYPE DOUBLETYPE
 // (non-perf critical or sensitive data types) 
#define SENSITIVE DOUBLETYPE
// WE ASSUME SENSITIVE>=REALTYPE !

// counter (integer) stuff where counts can exceed integer (2 billion)
#define COUNTTYPE DOUBLETYPE // can't make long long int work, so use double
//#define COUNTTYPE LONGLONGINTTYPE // can't make long long int work, so use double

#define PRODUCTION 0
// 0: full images, dumps, etc., few more key failure stderr/fail_file messages
// 1: only log density image since too many images (takes alot of time), no utoprim failure messages -- assume debug.out and debug???? will have info if needed


// 0: normal computational zones outputted on diagnostics
// else, FULLOUTPUT # of extra boundary zones on each side (if allowed by dimensionality)
// If # of requested boundary zones is larger than real #, then real # used

// ONLY CAN BE USED WITH USEMPI==0
#define FULLOUTPUT 0


#define MAILWHENDONE 1
#define MAILFROMREMOTE 0
#define REMOTEHOST "bh.astro.uiuc.edu"
#define EMAILADDRESS "jmckinney@cfa.harvard.edu"
#define EMAILMESSAGE "Done with GRMHD run DEFAULT"



// whether doing super long doubles
#define SUPERLONGDOUBLE 0


#define PERFTEST 0
// 0: don't perform performance test
// 1: do

#define DOAVG 0
// 0: don't do time average dumps, so don't allocate memory
// 1: do

// GODMARK: only setup for full Pi theta grid, not Pi/2 or any other theta grid.
// e.g. will fail for defcoord=3 Pi/2 grid.  Causes code to crash
#define DOJETDIAG 1
// 0: don't do jet diagnostics (ener, not flener that is always done)
// 1: do

#define DOAVG2 0 // only make 1 if DOAVG 1 above
// 0: don't split AVG file
// 1: do

#define DODEBUG 1
// 0: don't output debug dump file or ener file(ener is based on dump counts)
// 1: do

#define DOENODEBUG 0
// whether to do ENO debug

#define WENO_USE_PRIM_REDUCTION 1
// 0: do not perform the c2a limiting on the fluxes
// 1: perform the limiting of c2a flux reconstruction based on the limiting of the a2c primitive correction:
//    if the a2c reconstruction of the conserved quantities leads to very different primitives, then no a2c reconstruction is done; in this case
//    no reconstruction is done on the fluxes either.


#define WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS 1
// 0: do not perform additional reduction
// 1: perform additional reduction:
//    for NDIM > 1, subsequent 1d a2c reconstructions reduce to lower order (=no reconstruction) if any of the previous reconstructions reduced

#define LIMIT_FLUXC2A_PRIM_CHANGE 0
// 0: do not limit primtive correction due to c2a done on fluxes
// 1: limit c2a correction
 
#define DO_WENO_DEBUG 0
// whether to debug WENO

#define DOSUPERDEBUG 0
// 0: don't output super debug output
// 1: do


#define DODISS 0

// see diag_source()
#define DOLUMVSR 0

// see diag_source()
#define DODISSVSR 0

#define DOFIELDLINE 1
// 0: don't output energy@infinity and field line stuff
// 1: do

// whether to use Roe-averaged state in determining diffusive term in flux
#define ROEAVERAGEDWAVESPEED 0
// whether to use Athena's version of Roe averaging and estimating wave speeds
#define ATHENAROE 0

// whether to store wave speeds over whole grid before computing flux
// useful to avoid extra calculations if doing "local" LAXF or "global" LAXF.
// default HARM was using VERY local LAXF (only wavespeeds from primitives interpolated to the edge).
#define STOREWAVESPEEDS 0

// whether to use stored wave speeds for flux calculation (allows one to store wave speeds for interp.c but use true VERYLOCALVCHAR that is vchar's estimated from boundary as in standard HARM -- rather than maximum from center zones as done by STORED version of VERYLOCALVCHAR)
// silly choice would be 0 if VCHARTYPE=GLOBALVCHAR since interp.c doesn't use the stored wave speeds if this is the case.  So shouldn't store in this case since no point.
#define USESTOREDSPEEDSFORFLUX 0

#define VCHARTYPE VERYLOCALVCHAR

// whether to check on inversion and report problem
#define CHECKONINVERSION 1

#define PRECISEINVERSION 1
// whether we use ultimately precise inversion or "workable" inversion precision (see utoprim's for how defined)

#define WHICHVEL VELREL4
// which velocity to compute in (init can be anything (see init.c's for transformations)

#define WHICHEOM WITHGDET
//#define WHICHEOM WITHNOGDET

// if WHICHEOM==WITHNOGDET, then below determines which EOMs get what geometric prefactor.  Notice (as described in phys.c's source_conn() ) that geometry issue applies AFTER additions/subtractions of EOMs (as done by REMOVERESTMASSFROMUU).
#define REMOVERESTMASSFROMUU 1
// whether to subtract rest-mass from energy equation for divT=0 equation of motion
// 0: use MHD stress tensor with explicit rest-mass included
// 1: as 0, but subtract out rest-mass from conservation and flux terms for evolution
// 2: use MHD stress tensor withOUT rest-mass
// this changes mhd_calc() in phys.c and assumes rest of code uses mhd stress-energy tensor without restmass also!!
// DIAGNOSTICS also without rest-mass for UU terms.

#define RELTYPE RELEOM

#define EOMTYPE EOMGRMHD
//#define EOMTYPE EOMCOLDGRMHD

// whether to try other methods for the inversion if primary choices fails
// created because utoprim_2d_final fails for large b^2/rho when other methods (even utoprim_2d) do not fail.
#define UTOPRIMTRYAGAIN 1

// whether to coevolve the entropy to check for shock+reconnection dissipation
#define DOENTROPY DONOENTROPY // normal total energy equation
//#define DOENTROPY DOEVOLVECOMPAREENTROPY // coevolve and compare
//#define DOENTROPY DOEVOLVEDIRECTENTROPY // directly evolve entropy equation instead of total energy equation


// generically which type of entropy evolution to do, except not used when doing DOENTROPY==DOEVOLVEDIRECTENTROPY
#define WHICHENTROPYEVOLVE EVOLVESIMPLEENTROPY

// whether to call fixup() after initialization
#define FIXUPAFTERINIT 1

// whether to call fixup() after restart
#define FIXUPAFTERRESTART 1

// checks if solution is reasonble and fails a point if not
// only checks if b^2/\rho>>BSQORHOMAX, etc.
#define CHECKSOLUTION 1



// factor by which zone-zone values can be different for \gamma and internal energy, used of CHECKSOLUTION==1
#define GAMMAPERCDIFFMAX 2.0
#define UPERCDIFFMAX 10.0


// whether to interpolate extra quantity (e.g. v^2) and use to renormalize velocities after they are interpolated
#define DOEXTRAINTERP 0
// must also set RESCALEINTERP=1

// whether to control the limiter
#define LIMADJUST LIMITERFIXED
//#define LIMADJUST LIMITERBSQORHOANDU

// whether to limit all variables or just hydro variables
#define HYDROLIMADJUSTONLY 1

// determine how flux is computed per zone
#define FLUXADJUST FLUXFIXED
//#define FLUXADJUST FLUXBSQORHOANDU

// whether to change flux calculation for all variables or just hydro variables
#define HYDROFLUXADJUSTONLY 1

// whether to allow negative internal energy in substep
// UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED should be set in global.h
#define STEPOVERNEGU 1
// seems to make things worse (more failures, but also worse failures)

#define STEPOVERNEGRHO 1

#define STEPOVERNEGURHO 1




//#define UTOPRIMADJUST UTOPRIMSTATIC
#define UTOPRIMADJUST UTOPRIMAVG

// 0: return non-updated quantities
// 1: if u<0 or rho<0, then return updated quantities, else return non-updated quantities
#define UTOPRIMFAILRETURNTYPE UTOPRIMRETURNADJUSTED

#define COORDSINGFIX 1
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to

#define COORDSINGFIXCYL 0 //whether perform the same fix for CYLMINKMETRIC's z axis

//#define SINGSMALL (1E-3)
#define SINGSMALL (1E-20)
// Hawley uses 0.06283 (0.02Pi)

#define VOLUMEDIFF 0
// whether to use volume regularization or not

#define MINDT 1.e-20 // minimum dt

#define JONCHECKS 1 // for vel=3 extra checks (standard things)
#define JONCHECKS2 1 // check_pr() crazy thing for vel3 (crazy u^t check and fix -- 2D only)


#define FLOORDIAGS 1 // whether to compute floor diagnostics


#define ANALYTICCONNECTION 0 // whether to use analytic connection
// only applies to certain metric's and coordinates, see set_grid.c

#define ANALYTICSOURCE 0 // whether to use analytic source function
// very slow with gcc, some extend pgcc, not a bit problem with icc
// only applies to certain metric's and coordaintes, see phys.c

// whether to (with ENO method) avoid boundary zones if doing outflow on that boundary
#define OUTFLOWAVOIDBC 0

#define FLUXDIMENSPLIT QUASISTRANG
#define A2CDIMENSPLIT QUASISTRANG

// whether to have dq's : only used by second order methods
#define DODQMEMORY 1

#define DOENOFLUXMEMORY 1 // whether to allocate memory for doing ENO version of FINITEVOLUME scheme

// whether to adjust interpolation based upon boundary conditions
// defunct.  Not used anymore by ENO scheme and shouldn't be activated for other schemes.
// LEAVE 0!
#define BOUNDARYINTERPADJUST 0

// whether to compute \dot{F_r} as in diag.c
//disable the analysis if N1 <= 1 because leads to segfault if F1 == NULL (as is in the case of N1 == 1)
#define COMPUTEFRDOT 0

// whether to compute faraday and currents and output them to dumps
#define CALCFARADAYANDCURRENTS 1

#define WHICHCURRENTCALC CURRENTCALC1


// whether want faraday at t=0 for dump
#define FARADAYT0 1

// whether want partial currents for t=0 dump
#define CURRENTST0 1

// boundary condition mnemonics
// can be reset/added to by user init.h
#define OUTFLOW	0
#define SYMM	1
#define ASYMM	2
#define FIXED	3
#define POLARAXIS 4
#define FIXEDOUTFLOW 5 // means fixed inflow but allows outflow -- basically outflow if no inflow, but if inflow then set values to fixed quantities
#define NSSURFACE 6 // whatever in bounds.c for NS surface
#define PERIODIC 7 // periodic boundary conditions
#define OUTFLOWNOINFLOW 8  //copy if velocity directed outward; force velocity to zero and copy if directed inward (into the grid from the boundary)
#define RAMESHOUTFLOW 9 //OUTFLOW quantities according the Ramesh's power-law solution
#define BCEXTRAP 100  //atch: extrapolation into the boundary with high order (can be specified)
#define CYLAXIS 101   //atch: cylindrical axis BC
#define BCEXTRAP_VEL3 102  //atch: the same as BCEXTRAP but extrapolates 3-velocity rather than 4-velocity; important for Hubble flow
#define JETINJECTION 103 //Fixed boundary condition for jet injection surrounded by an outflow condition
#define BCU1EXTRAPOTHERFIXED 104
#define BCEXTRAPCONSTRAINED 105
#define RESCALEOUTFLOW 106
#define RESCALEFIXEDALLOUTFLOWU1 107
#define FIXED_RESCALEOUTFLOWU1 108
#define BONDIMDOTOUTFLOW 109
#define BONDIINTOUTFLOW 110

#define EVOLVECHECKS 0
// whether to check boundary conditions and limit gamma during advance()

// whether and which type of fixups to be used
#define FIXUPZONES FIXUP1ZONE

/** algorithmic choices **/

#define HLLBOUNDARY 0
// use HLL on boundary fluxes


#define FIXUPFLUX 0
// fix up the flux using fix_flux() in fixup.c
// caution...
// GODMARK: Not sure how metric and how emf works around axes.  Not unlike \Omega and such things that are "symmetric" while v^\phi might be antisymmetric?


#define ZEROOUTFLOWFLUX 0
// 0: don't do anything special, let boundary conditions control fux
// 1: zero the x1-velocity at the outflow boundary if inflow in flux routine so interpolation is overridden
// seems to cause problems at outer edge in internal energy

// seems to cause problems with EOMFFDE at boudaries.
// GODMARK



// assumes inner and outer theta boundaries are at 0 and Pi.
#define ZEROPOLEFLUX 0
// 0: don't do anything special
// 1: zero polar theta flux since should be 0

// seems to cause problems at 0 and pi boundaries with FFDE, or even just one boundary???
// GODMARK

// REASON: Cannot just set F=0.  The flux of the kinetic energy still has the pressure term even if the velocity into the wall is 0!



// whether to rescale interpolation
#define RESCALEINTERP 0
// 0: don't rescale
// 1: do rescale

#define BDIRCONT 1
// 0: don't
// 1: make field along flux direction continuous

// whether to use hyperbolic function for cmin and cmax instead of sharp max(0,l,r)
#define HYPERHLL 0

// whether to shift stencil inside horizon
#define HORIZONSUPERFAST 0



//////////////////////////////////
//
// which variable to interpolate
//
/////////////////////////////////


//#define VARTOINTERP PRIMTOINTERP_VSQ
#define VARTOINTERP PRIMTOINTERP



/////////////////////////////////////////////////////////
//
// STUFF FOR INTERPLINE.C
//
/////////////////////////////////////////////////////////



// whether to use appropriately centered primitive
// GODMARK: averaging is diffusive, so diffuses shock -- may make indicator weak
#define USEAVGPRIMITIVEFORWENOFLAT 1


// whether to invert average conserved quantity to get a primitive to use as shock indicator
// then correctly positioned in time, but more diffused value
// may be a bad idea since inversion is expensive, may find no solution, or may lead to negative internal energy in shock.
#define USEPRIMITIVEFROMAVGCONSERVED 0


// whether to use contact discontinuity indicator to steepen contacts in interpline.c
#define CONTACTINDICATOR 0

// send Sasha dP
#define COMPUTEDRHODP 1


// whether to reduce to 1-point stencil in case of superfast divergence at that point
#define SUPERFASTDIVREDUCE 0 //atch

// minimum preferred order
#define MINPREFORDER 3

// whether to compute shock indicator
#define SHOCKINDICATOR 1

// interp general stuff
#define WHICHPARA PARA4


#define BONDI_BOUNDARY_SET_PL_PR 0  //do not analytically set p_l & p_r at the outer boundary for the Bondi problem


#define NOSPECIALFIELD 0
#define PULSARFIELD 1

#define VARTOINTERPFIELD NOSPECIALFIELD

