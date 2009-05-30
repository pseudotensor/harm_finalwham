
/* 
 *
 * invert U (conserved variables) to obtain
 * p (primitive variables).  
 * NB: pr must contain initial guess
 *
 */

#include "decs.h"

static struct of_geom *ptrlgeom;
static FTYPE U_target[NPR];
static FTYPE pr0[NPR];
static int numnormterms;
static FTYPE normf;
static int whichconsglobal,primtoUcons;
/* pr *MUST* contain initial guess */

// todo:
// 1) choose tolx/tolf more intelligently (perhaps scales with
// difference between maximum and minimum density)

static int *glpflag; // global pflag for local file

/* pr *MUST* contain initial guess */
int Utoprim(int whichcons, FTYPE *U, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr)
{
  FTYPE priter[NPR];
  FTYPE prguess[NPR];
  FTYPE tolx, tolf,tolxallowed, tolfallowed, tolxreport, tolfreport;
  int ntrial, k, mnewtfailed,mintrial;
  FTYPE Ustart[NPR];
  FTYPE dU[NPR];
  struct of_state q;
  int i,j;
  int nutoprims,countutoprims;
  FTYPE frac,SUPERFRAC;
  int faildebug1(FTYPE *pr0,struct of_geom *ptrgeom);
  int faildebug2(FTYPE *pr0,struct of_geom *ptrgeom);
  extern int invertentropyflux_calc(FTYPE entropyflux,int dir, struct of_state *q, FTYPE *pr);
  void fixUtarget(int which,FTYPE *Ui,FTYPE *Uf);
  extern int mnewt(int whethertoreport,int ntrail, int mintrial, FTYPE *p, int n, int startp, FTYPE tolx,
		 FTYPE tolf, FTYPE tolxallowed,
		 FTYPE tolfallowed, FTYPE tolxreport,
		 FTYPE tolfreport);
  int usrfun(FTYPE *pr, int n, FTYPE *beta, FTYPE **alpha,FTYPE*normf);
  int Utoprim_nonrel(FTYPE *U, struct of_geom *ptrgeom, FTYPE *pr);
  FTYPE prnonrel[NPR];
  FTYPE prrel[NPR];
  int pl;
  int flagifnonrel;
  int whethertoreport;


#if(EOMTYPE==EOMGRMHD)
#define INVERTNPR (5) // always same
#define STARTPINVERT (RHO) // where to start in primitive space
#elif(EOMTYPE==EOMCOLDGRMHD)
#define INVERTNPR (4) // always same
#define STARTPINVERT (RHO) // where to start in primitive space
  // need to skip UU if using this (SUPERGODMARK -- not done yet)
#elif(EOMTYPE==EOMFFDE)
#define INVERTNPR (3) // always same
#define STARTPINVERT (U1) // where to start in primitive space
#endif



  // assign global int pointer to lpflag pointer
  glpflag=lpflag;
  *glpflag= UTOPRIMNOFAIL; // set as no failure at first

  ptrlgeom=ptrgeom;
  whichconsglobal=whichcons;


  // assume normally do not want special messages
  whethertoreport=0;


  /////////////
  //
  // check (disabled since should allow negative density so RK can recover)
  //
  ////////////
  if(0&&( (EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD) )){ // don't check for now, just fail later
    if (U[RHO] < 0.) {
      if (fail(FAIL_UTOPRIM_NEG) >= 1)
	return (1);
    }
  }
  
  ////////////////
  //
  // evolve B before (also better than initial guess)
  /* solution is known for B, so best guess is solution (from point of view of other variables) */
  //
  ///////////////
  for (k = B1; k <= B3; k++) pr[k] = U[k];




  //////////////////////
  //
  // rest determines rho,u,u^i for GRMHD or u^i for FFDE
  //


#if(PRECISEINVERSION)
  ntrial = 200;
  mintrial = 3;
  tolx = tolf = 1.e-11; // 1E-11 is best mnewt() can do with DODAMP=2
  tolxallowed = tolfallowed = 1.e-8 ;
  tolxreport=1.e4*tolx;
  tolfreport=1.e4*tolf;
  
#else
  // choice
  ntrial = 20;
  mintrial=2;
  tolx = 1.e-10;
  tolf = 1.e-10;
  tolxallowed=tolfallowed=1.e-4;
  tolxreport=1.e3*tolx;
  tolfreport=1.e3*tolf;
#endif



  /////////////
  //
  // Setup inversion
  //
  ////////////


  // store initial guess and U as target
  PLOOP(k) pr0[k]=priter[k]=pr[k];


#if(EOMTYPE==EOMGRMHD)
  PLOOP(k) U_target[k] = U[k];

#define NORMMETHOD (1)
// 0=use Utarget[RHO]
// 1=use generaleralized mean so matrix stays well conditioned

  // assume normally want special messages
  whethertoreport=1;

  if((whichcons==EVOLVENOENTROPY)||(whichcons==EVOLVESIMPLEENTROPY)){
    primtoUcons=UNOTHING;
    // normal total energy inversion
    // if evolvesimpleentropy, then just invert in utoprim() after mnewt succeeds
  }
  else if(whichcons==EVOLVEFULLENTROPY){
    primtoUcons=UENTROPY; // implicit UNOTHING with entropy cons. quant. assigned to UU cons. quant.
    // overwrite total energy EOM with entropy EOM
    U_target[UU]=U_target[ENTROPY];
    // assumes dudp_calc replaces total energy equation with entropy equation
  }


  // assume don't want special messages if just using entropy evolution for comparison and not primary evolution
  if(whichcons==EVOLVEFULLENTROPY && DOENTROPY!=DOEVOLVEDIRECTENTROPY){
    whethertoreport=0;
    // then simplify inversion a bit
    ntrial = 20;
    mintrial=2;
    tolx = 1.e-6;
    tolf = 1.e-6;
    tolxallowed=tolfallowed=1.e-4;
    tolxreport=1.e3*tolx;
    tolfreport=1.e3*tolf;
  }


#elif(EOMTYPE==EOMCOLDGRMHD)

  // SUPERGODMARK: Not done yet

#define NORMMETHOD (1)  // no choice
#define WHICHFFDEINVERSION 0 // choice (see below)

  primtoUcons=UNOTHING;

  fixUtarget(WHICHFFDEINVERSION,U,U_target);

  // zero internal energy density so can use standard dU/dp
  priter[UU]=0;

#elif(EOMTYPE==EOMFFDE)

#define NORMMETHOD (1)  // no choice
#define WHICHFFDEINVERSION 0 // choice (see below)

  primtoUcons=UNOTHING;

  fixUtarget(WHICHFFDEINVERSION,U,U_target);

  // zero densities so can use standard dU/dp
  priter[RHO]=priter[UU]=0;

#endif


  /////////////
  //
  // checks if initial pr gives good u^t and determines U(pr)
  //
  ////////////
  if (get_state(pr, ptrgeom, &q) >= 1){
    //    dualfprintf(fail_file,"bad guess: i=%d j=%d k=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
    //    PLOOP(k) dualfprintf(fail_file,"bad guess: pr=%21.15g\n",pr[k]);
    // then guess bad, adjust guess
    set_atmosphere(0,WHICHVEL,ptrgeom,prguess);
    for(k=U1;k<=U3;k++) pr[k]=prguess[k]; // only reset velocities
    // try again
    if (get_state(pr, ptrgeom, &q) >= 1)
      FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
    //    PLOOP(k) dualfprintf(fail_file,"good guess: pr=%21.15g\n",pr[k]);
  }
  if (primtoU(primtoUcons,pr, &q, ptrgeom, Ustart) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 1);
  // now we have U(pr)



  ////////////////////
  //
  // if u^t=1, then can do non-relativistic inversion
  //
  ////////////////////
#define DONONRELINV 1


#if(DONONRELINV&&( (EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD) ) )
  flagifnonrel=0;
  if(fabs(q.ucon[TT]-1.0)<NUMEPSILON*5.0){
    flagifnonrel=1;

    if(Utoprim_nonrel(U_target,ptrgeom,pr)>0){
      *glpflag= UTOPRIMFAILRHONEG;
      flagifnonrel=0; //then can't trust non-rel inversion
    }
    else{
    
      // report if failure.  Only negative densities can occur with non-relativistic inversion
      if((pr[RHO]<=0.)&&(pr[UU]>=0.)) *glpflag= UTOPRIMFAILRHONEG;
      if((pr[RHO]>=0.)&&(pr[UU]<=0.)) *glpflag= UTOPRIMFAILUNEG;
      if((pr[RHO]<=0.)&&(pr[UU]<=0.)) *glpflag= UTOPRIMFAILRHOUNEG;

      if(*glpflag!=UTOPRIMNOFAIL){
	//	PLOOP(k) dualfprintf(fail_file,"U[%d]=%21.15g pr0[%d]=%21.15g pr[%d]=%21.15g\n",k,U_target[k],k,pr0[k],k,pr[k]);
      }
      PLOOP(pl) prnonrel[pl]=pr[pl];
      
      /////////////
      //
      // checks if initial pr gives good u^t and determines U(pr)
      //
      ////////////
      if (get_state(pr, ptrgeom, &q) >= 1)
	FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
      if(fabs(q.ucon[TT]-1.0)>NUMEPSILON*5.0) flagifnonrel=0;
    }
  }
  else{
    flagifnonrel=0;
    PLOOP(pl) prnonrel[pl]=-1.0;
  }
#else
  flagifnonrel=0;
  PLOOP(pl) prnonrel[pl]=-1.0;
#endif

  //  dualfprintf(fail_file,"uu0diff1=%21.15g flagifnonrel=%d\n",fabs(q.ucon[TT]-1.0),flagifnonrel);




  if(flagifnonrel==0){
    //////////////////
    //
    // setup mnewt
    //
    //////////////////
    mnewtfailed=1; // assume failed
    nutoprims=1; // >1 usually, but 1 means directly to static if failed (with one chance for SUPERFRAC)
    countutoprims=nutoprims;
    frac=1.0/((FTYPE)nutoprims); // so that by nutoprims adjustments we are back at U_target=Ustart
    SUPERFRAC=0.0001;
    PLOOP(k) dU[k]=U_target[k]-Ustart[k];


    //  PLOOP(k) dualfprintf(fail_file,"IN: i=%d j=%d k=%d :: pr0[%d]=%21.15g U[%d]=%21.15g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,k,pr0[k],k,U_target[k]);
    //////////////////
    //
    // mnewt (loop)
    //
    //////////////////
    while(mnewtfailed){
      if(mnewtfailed){
	if(countutoprims==-1){
	  // forced to use static solution. Probably better ways than a
	  // linear approx, but better than nothing.  Must improve flux,
	  // or lower courant factor to otherwise improve this situation
	  break;
	}
	if(countutoprims==0){
	  // force to be nearly identical instead of exactly identical, last hope!
	  PLOOP(k) U_target[k]=Ustart[k]+dU[k]*SUPERFRAC;
	}
	else{
	  // Only failed because U_target is trying to send u^t to imaginary values, or densities to negative values
	  // try backup up towards original solution in conservative variable space
	  PLOOP(k) U_target[k]=Ustart[k]+dU[k]*(FTYPE)countutoprims*frac;
	}
      }

      // mnewt is set to only take 5 entries and deal with k=0-4 in the inversion.
      //
      // notice that priter is used.  Assumes setup Utarget to be consistent with priter and INVERTNPR
      //
      mnewtfailed = mnewt(whethertoreport,ntrial, mintrial, priter - 1, INVERTNPR , STARTPINVERT, tolx, tolf, tolxallowed, tolfallowed, tolxreport, tolfreport);

      if(mnewtfailed) countutoprims--;
      if(nutoprims==1){
	break; // if failed, static immediately if failure since apparently goes that way anyways
	// if good, then break cause we are done
      }

    }

    /////////////////////
    //
    // convert the priter->pr
    //
    /////////////////////
    PLOOP(k) prrel[k]=priter[k];


  }
  else{
    mnewtfailed=0; // don't trigger debug stuff for mnewt if not using it
    PLOOP(k) prrel[k]=-1;
  }



  /////////////////////
  //
  // ASSIGN ANSWER from mnewt or non-rel inversion
  //
  /////////////////////

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)



#if(0)
  //  PLOOP(pl)  dualfprintf(fail_file,"%ld %d %d inversions prnonrel[%d]=%21.15g prrel[%d]=%21.15g\n",nstep,steppart,ptrgeom->i,pl,prnonrel[pl],pl,prrel[pl]);
  //  if(flagifnonrel)  PLOOP(pl)  if(fabs(prnonrel[pl]-prrel[pl])>0.0) dualfprintf(fail_file,"%ld %d %d inversions prnonrel[%d]=%21.15g prrel[%d]=%21.15g diff=%21.15g\n",nstep,steppart,ptrgeom->i,pl,prnonrel[pl],pl,prrel[pl],prnonrel[pl]-prrel[pl]);
  if(flagifnonrel)  PLOOP(pl)  if(fabs(prnonrel[pl]-prrel[pl])/fabs(fabs(prnonrel[pl])+fabs(prrel[pl])+SMALL)>0.0) dualfprintf(fail_file,"%ld %d %d inversions prnonrel[%d]=%21.15g prrel[%d]=%21.15g diff=%21.15g\n",nstep,steppart,ptrgeom->i,pl,prnonrel[pl],pl,prrel[pl],prnonrel[pl]-prrel[pl]);
#endif

  // default is to use nonrel if was good solution
  if(flagifnonrel)  PLOOP(pl) pr[pl]=prnonrel[pl];
  else PLOOP(pl) pr[pl]=prrel[pl];

#elif(EOMTYPE==EOMFFDE)
  // v^i
  // don't overwrite rest
  pr[U1]=prrel[U1];
  pr[U2]=prrel[U2];
  pr[U3]=prrel[U3];
#endif


  if(whichcons==EVOLVEFULLENTROPY){
    // in case user wants internal energy from entropy evolution in pr[ENTROPY]
    pr[ENTROPY]=pr[UU];
  }



  //////////////////
  //
  // check if mnewt failed or got negative densities
  //
  //////////////////
  if(mnewtfailed){



    //    PLOOP(k) dualfprintf(fail_file,"U_target[%d]=%21.15g\n",k,U_target[k]);



    // reset to initial pr since new pr is bad
    PLOOP(k) pr[k]=pr0[k];

    // report failure
    *glpflag= UTOPRIMFAILCONV;

    // then old p should be new p, do nothing since initial guess must be final solution
    if(whethertoreport && debugfail>=1) dualfprintf(fail_file,"static solution at t=%21.15g step=%ld i=%d j=%d k=%d wtf=%d\n",t,realnstep,ptrgeom->i,ptrgeom->j,ptrgeom->k,debugfail);

    // no need to continue unless want debug info, then pflag is messed up
    return(0); 
    // remove above return if want to get debug below
  }// otherwise got good solution
  else{
    // check densities for positivity
    if( ( (EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD))&&((pr[RHO]<0.0)||(pr[UU]<0.0))){
      // it doesn't seem reasonable to just "fix" negative densities since they could be arbitrarily negative and the other primitive variables depend on the negativity of the density which is unphysical.  So we fail if negative.  Could fail for a lower threshold on density than the floor to avoid those numbers as well.
      if((pr[RHO]<=0.)&&(pr[UU]>=0.)) *glpflag= UTOPRIMFAILRHONEG;
      if((pr[RHO]>=0.)&&(pr[UU]<=0.)) *glpflag= UTOPRIMFAILUNEG;
      if((pr[RHO]<=0.)&&(pr[UU]<=0.)) *glpflag= UTOPRIMFAILRHOUNEG;

      // mostly we hit pr[UU]<0, or maybe never pr[RHO]<0 and pr[UU]<0 alot
      if(debugfail>=2) dualfprintf(fail_file,"utoprim found negative density: t=%21.15g step=%ld i=%d j=%d k=%d %21.15g %21.15g\n",t,realnstep,ptrgeom->i,ptrgeom->j,ptrgeom->k,pr[RHO],pr[UU]);
      return(0);
    }
  }



  ////////////////////
  //
  // make sure solution is valid
  //
  ////////////////////
#if(WHICHVEL==VEL3)
  // this is just a check, fixup() will soon do both in step_ch.c

  // currently need to check this since last iteration may lead to bad pr
  if(jonchecks) fixup1zone(pr,ptrgeom,0); // actually this is done in a general fixup call soon after this point in step_ch.c
  // check and see if u^t is actually a good solution now we think we have the solution? or fix it if u^t is bad
  if(jonchecks) if(check_pr(pr,pr0,ptrgeom,0)>=1)
    FAILSTATEMENT("utoprim.c:Utoprim()", "mnewt check:check_pr()", 1);
#endif



  ////////////////////
  //
  // determine entropy inversion if not doing full entropy evolution
  //
  ///////////////////

  if(whichcons==EVOLVESIMPLEENTROPY){
    // if evolvesimpleentropy, then just invert in utoprim() after mnewt succeeds
    if (get_state(pr, ptrgeom, &q) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 2);
    invertentropyflux_calc(U_target[ENTROPY],0,&q,pr);// doesn't use pr[UU], but uses rest of updated primitives
    
  }// otherwise no entropy evolution or direct full entropy evolution




  /////////////////////////////
  //
  // debug stuff
  // don't fail anymore even if failed-static solution
  //
  /////////////////////////////
  if ((debugfail>=2)&&(mnewtfailed >= 1)) {

    //    faildebug1(pr0,ptrgeom);
    faildebug2(pr0,ptrgeom);

    if (fail(FAIL_UTOPRIM_TEST) >= 1)
      return (1);

    
    if(!jonchecks){
      if (fail(FAIL_UTOPRIM_TEST) >= 1)
	return (1);
    }
  }



  ////////////////////
  //
  // otherwise just safely fail or succeed
  //
  ////////////////////
  *glpflag = UTOPRIMNOFAIL;
  return (0);


}


// Nonrelativistic inversion
// not setup for gravity (i.e. assumes g_{tt}=-1)
// probably not right for magnetic fields

/*

Conserved quantities:


rest-mass : U[RHO] = \rho

REMOVERESTMASSFROMUU==2
-Energy   : U[UU]  = -\rho v^2/2 - u  = -(\rho v^2/2 + u)

Momentums : U[UU+j] = \rho v_j  j=1..3





*/
int Utoprim_nonrel(FTYPE *U, struct of_geom *ptrgeom, FTYPE *pr)
{
  int j,k;
  FTYPE vd[NDIM];
  FTYPE vsq;
  FTYPE rhovsq,rhosqvsq;
  int pl;
  //  FTYPE bt;
  FTYPE Bsq;//,bsq;

  if(U[RHO]<=0.0){
    dualfprintf(fail_file,"density was <=0.0\n");
    //    myexit(0);
    return(1);
  }

  // field B^i
  //  pr[B1]=U[B1];
  //  pr[B2]=U[B2];
  //  pr[B3]=U[B3];

  // \rho
  pr[RHO]=U[RHO];

  // v_1 = (\rho v_1)/\rho
  vd[1]=U[U1]/U[RHO];
  // v_2 = (\rho v_2)/\rho
  vd[2]=U[U2]/U[RHO];
  // v_3 = (\rho v_3)/\rho
  vd[3]=U[U3]/U[RHO];

  // b^t = B^i u_i
  //  bt=0.0;
  //  SLOOPA(j) bt+=pr[B1+k-1]*vd[k];

  // get v^i, which are primitives
  pr[U1]=pr[U2]=pr[U3]=0.0;
  SLOOP(j,k) pr[UU+j]+=vd[k]*ptrgeom->gcon[j][k];

#if(0)
  // this method leads to catastrophic cancellation of some kind
  
  // get square of velocity: (v)^2
  vsq=0.0;
  SLOOP(j,k) vsq+=pr[UU+j]*vd[k];

  rhovsq=U[RHO]*vsq;

#else

  // this method required to avoid catastrophic cancellation due to how $\rho$ is used
  // that changes the machine error
  rhosqvsq=0.0;
  SLOOP(j,k) rhosqvsq+=U[UU+k]*U[UU+j]*ptrgeom->gcon[j][k];
  rhovsq=rhosqvsq/U[RHO];


#endif

  // (B)^2
  Bsq=0.0;
  SLOOP(j,k) Bsq+=pr[B1+j-1]*pr[B1+k-1]*ptrgeom->gcov[j][k];

  //  dualfprintf(fail_file,"rho=%21.15g v1=%21.15g vsq=%21.15g\n",pr[RHO],pr[U1],vsq);
  //  SLOOP(j,k) dualfprintf(fail_file,"%d %d : %21.15g\n",j,k,ptrgeom->gcov[j][k]);

  //  bsq=Bsq+bt*bt*(vsq+2.0*(1.0+ptrgeom->gcov[0][0])); // machine precision accurate for Minkowski

#if(REMOVERESTMASSFROMUU==2)

  // u = (\rho v^2/2+u+B^2/2) - rho*v^2/2 - Bsq/2 - \rho \Phi
  pr[UU]=-U[UU] - 0.5*rhovsq - Bsq*0.5 + 0.5*U[RHO]*ptrgeom->gcovpert[TT];


#elif( (REMOVERESTMASSFROMUU==1)||REMOVERESTMASSFROMUU==0)// stupid idea, but let user do it

  // u = (\rho+\rho v^2/2+u+B^2/2) - rho - rho*v^2/2 - Bsq/2 - \rho \Phi
  pr[UU]=-U[UU] - U[RHO] - 0.5*rhovsq - Bsq*0.5 + 0.5*U[RHO]*ptrgeom->gcovpert[TT];

#endif

  PLOOP(pl) if(!finite(pr[pl])){
    dualfprintf(fail_file,"Bad primitive pr[%d]=%21.15g\n",pl,pr[pl]);
    return(1);
  }

  //    PLOOP(pl) dualfprintf(fail_file,"U[%d]=%21.15g\n",pl,U[pl]);

    //  PLOOP(pl)  dualfprintf(fail_file,"%ld %d %d nonrel inversion pr[%d]=%21.15g\n",nstep,steppart,ptrgeom->i,pl,pr[pl]);

  return(0);// can't fail
}

// convert between normal U and iterative U
void fixUtarget(int which,FTYPE *Ui,FTYPE *Uf)
{


  // U1,U2,U3 for U_target are just effective storage locations
  if(which==0){
    // energy and angular momentum conserving inversion #1
    // dumps truncation error into U2
    Uf[0] = Ui[UU];
    Uf[1] = Ui[U1];
    Uf[2] = Ui[U3];
  }
  else if(which==1){
    // energy and angular momentum conserving inversion #2
    // dumps truncation error into U1
    Uf[0] = Ui[UU];
    Uf[1] = Ui[U2];
    Uf[2] = Ui[U3];
  }
  else if(which==2){
    // as with first analytic inversion
    Uf[0] = Ui[U1];
    Uf[1] = Ui[U2];
    Uf[2] = Ui[U3];
  }
  
}




/* auxiliary function required by mnewt */
int usrfun(FTYPE *prguess, int n, FTYPE *beta, FTYPE **alpha,FTYPE*norm)
{
  static FTYPE **alpha5;
  FTYPE prstart[NPR];
  FTYPE *pr;
  static FTYPE U[NPR],U_curr[NPR];
  struct of_state q;
  int i = 0, j = 0, k = 0;
  int failreturn=0;
  void removerow(int row, FTYPE **alpha5,FTYPE**alpha);
  static int firstc=1;

  if (firstc) {
    firstc = 0;
    alpha5 = dmatrix(1, 5, 1, 5);
  }


  pr=prguess+1;

  // Doing this constrained iteration may lead to non-convergence, but then we adjust U if that happens until we converge.

  // store old pr
  PLOOP(k) prstart[k]=pr[k];

  // these checks seem to cause problems with convergence of any kind // GODMARK
  // we need to make sure densities are ok while iterating.
#if(WHICHVEL==VEL3)
  if(0&&jonchecks) fixup1zone(pr,ptrlgeom,0);

  // once pr lead to imaginary u^t, state is invalid and U cannot be found.
  // Thus, need to adjust pr while iterating to keep pr in line
  // essentially we are doing a constrained damped Newton's method.
  if(0&&jonchecks) failreturn=check_pr(pr,pr0,ptrlgeom,-100);
  if(0&&jonchecks) if(failreturn>=1)
    FAILSTATEMENT("utoprim.c:usrfun()", "mnewt check:check_pr()", 2);

  // just check
  if(jonchecks) failreturn=check_pr(pr,pr0,ptrlgeom,-2);
#endif

  if(!failreturn){

    // if didn't fail, we are guarenteed a good state and U.
    
    // problem here is that now pr is different and doesn't help
    // converge to Utarget, maybe kills possibility of convergence which
    // would later be just fixed at the end
    // so should probably adjust Utarget somehow to compensate for adjustment in pr.  Use dudp to fix U_target along the way
    
    /* normalize error = beta to \rho u^t */
    failreturn=get_state(pr, ptrlgeom, &q);
    if(failreturn>=1)
      FAILSTATEMENT("utoprim.c:usrfun()", "get_state()", 1);

    failreturn=primtoU(primtoUcons,pr, &q, ptrlgeom, U);
    if(failreturn>=1)
      FAILSTATEMENT("utoprim.c:usrfun()", "primtoU()", 1);

#if(EOMTYPE==EOMGRMHD)
    PLOOP(k) U_curr[k]=U[k];
#elif(EOMTYPE==EOMFFDE)
    // convert from normal U to iterative U
    fixUtarget(WHICHFFDEINVERSION,U,U_curr);
#endif


    failreturn=dudp_calc(whichconsglobal,pr, &q, ptrlgeom, alpha5);
    if(failreturn>=1)
      FAILSTATEMENT("utoprim.c:usrfun()", "dudp_calc()", 1);


    // convert from normal alpha to iterative alpha
#if(EOMTYPE==EOMGRMHD)
    for (j = 0; j < INVERTNPR ; j++){
      for (k = 0; k < INVERTNPR ; k++){
	alpha[j+1][k+1]=alpha5[j+1][k+1];
      }
    }    
#elif(EOMTYPE==EOMFFDE)

#if(WHICHFFDEINVERSION==0)
    removerow(U2,alpha5,alpha);
#elif(WHICHFFDEINVERSION==1)
    removerow(U1,alpha5,alpha);
#elif(WHICHFFDEINVERSION==2)
    removerow(UU,alpha5,alpha);
#endif

#endif

    /////////////////////
    //
    // normalize matrix
    //
    ////////////////////
#if(NORMMETHOD==0)
    *norm=1.0/U_target[RHO];
#elif(NORMMETHOD==1)
    // more general normalization
    *norm=0.0;
    numnormterms=0;
    for (j = 0; j < INVERTNPR ; j++){
      for (k = 0; k < INVERTNPR ; k++){
	// GODMARK: correction:
	if(fabs(alpha[j + 1][k + 1])>NUMEPSILON){
	// old version (bug)
	//	if(alpha[j + 1][k + 1]>NUMEPSILON){
	  *norm+=fabs(alpha[j + 1][k + 1]);
	  numnormterms++;
	}
      }
    }
    *norm=(FTYPE)(numnormterms)/(*norm); // (i.e. inverse of average)

#endif
    
    for (j = 0; j < INVERTNPR ; j++)
      for (k = 0; k < INVERTNPR ; k++)
	alpha[j + 1][k + 1] *= (*norm);

    
    // below assumes alpha at U_curr isn't too different from at U_target
    // doesn't seem to work
    /*
      for (j = 0; j < INVERTNPR ; j++){
      for (k = 0; k < INVERTNPR ; k++){
      // should really use average alpha from initial pr and new pr, but try just new pr's dudp
      U_target[j]+=U_target[RHO]*alpha[j+1][k+1]*(pr[k]-prstart[k]);
      }
      }
    */   

    // determine normalized error
    for (k = 0; k < INVERTNPR ; k++)
      beta[k + 1] = (U_curr[k] - U_target[k]) *(*norm);
    
  }
  else{
    // error is huge if failed, a marker for failure
    for (k = 0; k < INVERTNPR ; k++)
      beta[k + 1] = 1.E30;
    
  }
  
  return (0);
}


// this function remove one U row and makes alpha[1][1] start at alpha5[1+{3 U's}][1+{U1..U3}]
void removerow(int row, FTYPE **alpha5, FTYPE**alpha)
{
  // alpha[i][j]=dU^i/dp^j
  int j,k;
  int jj,kk;


  
  // FFDE: j=[UU...U3] inclusive
  for(j = UU; j < UU+INVERTNPR+1 ; j++){ // +1 assumes killing 1 row that otherwise exists

    // FFDE: k=[U1..U3] inclusive
    for(k = STARTPINVERT; k < STARTPINVERT+INVERTNPR ; k++){
      
      if(j<row){
	jj=(j+1)-U1;
	kk=k-STARTPINVERT;
      }
      else if(j==row){
	continue;
      }
      else if(j>row){
	jj=(j+1)-(U1+1);
	kk=k-STARTPINVERT;
      }
      alpha[jj+1][kk+1]=alpha5[j+1][k+1];
    }
  }
  


}








///////////////////////////////////////////////////////////////////
//
// debug stuff
//
///////////////////////////////////////////////////////////////////



int faildebug1(FTYPE *pr0,struct of_geom *ptrgeom)
{
  FTYPE tolx, tolf;
  int ntrial, k, mnewtfailed;
  FTYPE Ustart[NPR];
  FTYPE dU[NPR],dUcomp[NUMSOURCES][NPR];
  FTYPE bsq;
  FTYPE mhd[NDIM][NDIM];
  FTYPE ucon[NDIM];
  FTYPE flux[NDIM][NPR];
  FTYPE **alpha;
  struct of_state q;
  int i,j;
  int myii,myjj,mykk;

  if (get_state(pr0, ptrgeom, &q) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 2);
  if (primtoU(UNOTHING,pr0, &q, ptrgeom, Ustart) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 2);
  
  mhd_calc(pr0,0, ptrgeom, &q, mhd[0]);
  mhd_calc(pr0,1, ptrgeom, &q, mhd[1]);
  mhd_calc(pr0,2, ptrgeom, &q, mhd[2]);
  mhd_calc(pr0,3, ptrgeom, &q, mhd[3]);
  primtoflux(UNOTHING,pr0, &q, 1,ptrgeom, flux[0]);
  primtoflux(UNOTHING,pr0, &q, 2,ptrgeom, flux[1]);
  primtoflux(UNOTHING,pr0, &q, 3,ptrgeom, flux[2]);
  source(pr0, ptrgeom, dUcomp,dU);
  bsq_calc(pr0, ptrgeom, &bsq);
  
  
  alpha=dmatrix(1, 5, 1, 5);

  if (dudp_calc(EVOLVENOENTROPY,pr0, &q, ptrgeom, alpha) >= 1)
    FAILSTATEMENT("utoprim.c:usrfun()", "dudp_calc()", 1);
  // more general normalization
  // copied from usrfun
#if(NORMMETHOD==0)
  normf=1.0/U_target[RHO];
#elif(NORMMETHOD==1)
  normf=0.0;
  numnormterms=0;
  for (j = 0; j < INVERTNPR ; j++){
    for (k = 0; k < INVERTNPR ; k++){
      if(alpha[j + 1][k + 1]>NUMEPSILON){
	normf+=fabs(alpha[j + 1][k + 1]);
	numnormterms++;
      }
    }
  }
  normf=(FTYPE)(numnormterms)/(normf); // (i.e. inverse of average)
#endif    
  for (j = 0; j < INVERTNPR; j++)
    for (k = 0; k < INVERTNPR; k++)
      alpha[j + 1][k + 1] *= (normf);
  
  dualfprintf(fail_file,"myTud={");
  for(i=0;i<NDIM;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NDIM;j++){
      dualfprintf(fail_file, "%21.15g``20 ",mhd[i][j]);
      if(j<NDIM-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<NDIM-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  
  
  dualfprintf(fail_file,"myF={");
  for(i=0;i<=2;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NPR;j++){
      dualfprintf(fail_file, "%21.15g``20 ",flux[i][j]);
      if(j<NPR-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<2)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  
  
  dualfprintf(fail_file,"myucon={%21.15g, %21.15g, %21.15g, %21.15g}\n",q.ucon[0],q.ucon[1],q.ucon[2],q.ucon[3]);
  dualfprintf(fail_file,"myucov={%21.15g, %21.15g, %21.15g, %21.15g}\n",q.ucov[0],q.ucov[1],q.ucov[2],q.ucov[3]);
  dualfprintf(fail_file,"mybcon={%21.15g, %21.15g, %21.15g, %21.15g}\n",q.bcon[0],q.bcon[1],q.bcon[2],q.bcon[3]);
  dualfprintf(fail_file,"mybcov={%21.15g, %21.15g, %21.15g, %21.15g}\n",q.bcov[0],q.bcov[1],q.bcov[2],q.bcov[3]);
  dualfprintf(fail_file,"myg=%21.15g\n",ptrgeom->g);
  
  //    PLOOP(k) dualfprintf(fail_file,"pr0[%d]=%21.15g Ustart[%d]=%21.15g U_target[%d]=%21.15g\n",k,pr0[k],k,Ustart[k],k,U_target[k]);
  dualfprintf(fail_file,"pr0={");
  PLOOP(k){
    dualfprintf(fail_file,"%21.15g``20 ",pr0[k]);
    if(k<NPR-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  dualfprintf(fail_file,"myU0={");
  PLOOP(k){
    dualfprintf(fail_file,"%21.15g``20 ",Ustart[k]);
    if(k<NPR-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  dualfprintf(fail_file,"mydU={");
  PLOOP(k){
    dualfprintf(fail_file,"%21.15g``20 ",dU[k]);
    if(k<NPR-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  dualfprintf(fail_file,"Utarget={");
  PLOOP(k){
    dualfprintf(fail_file,"%21.15g``20 ",U_target[k]);
    if(k<NPR-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  dualfprintf(fail_file,"normf=%21.15g``20\n",normf);
  dualfprintf(fail_file,"mybsq=%21.15g``20\n",bsq);
  
  dualfprintf(fail_file,"alpha={");
  for(i=1;i<=5;i++){
    dualfprintf(fail_file,"{");
    for(j=1;j<=5;j++){
      dualfprintf(fail_file, "%21.15g``20 ",alpha[i][j]);
      if(j<5)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<5)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  

#if(MCOORD!=CARTMINKMETRIC)
  myii=icurr;
  myjj=jcurr;
  mykk=kcurr;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif

  
  dualfprintf(fail_file,"myconn={");
  for(i=0;i<NDIM;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NDIM;j++){
      dualfprintf(fail_file,"{");
      for(k=0;k<NDIM;k++){
	dualfprintf(fail_file, "%21.15g``20 ",conn[myii][myjj][mykk][i][j][k]);
	if(k<NDIM-1)       dualfprintf(fail_file,", ");
      }
      dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
      if(j<NDIM-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<NDIM-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  
  dualfprintf(fail_file,"mygcon={");
  for(i=0;i<NDIM;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NDIM;j++){
      dualfprintf(fail_file, "%21.15g``20 ",gcon[myii][myjj][mykk][ptrgeom->p][i][j]);
      if(j<NDIM-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<NDIM-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  
  dualfprintf(fail_file,"mygcov={");
  for(i=0;i<NDIM;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NDIM;j++){
      dualfprintf(fail_file, "%21.15g``20 ",gcov[myii][myjj][mykk][ptrgeom->p][i][j]);
      if(j<NDIM-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<NDIM-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");

  return(0);
}




// same as faildebug1 but with column formatted output, instead of copy/paste into mathematica format
int faildebug2(FTYPE *pr0,struct of_geom *ptrgeom)
{
  FTYPE tolx, tolf;
  int ntrial, k, mnewtfailed;
  FTYPE Ustart[NPR];
  FTYPE dU[NPR],dUcomp[NUMSOURCES][NPR];
  FTYPE bsq;
  FTYPE mhd[NDIM][NDIM];
  FTYPE ucon[NDIM];
  FTYPE flux[NDIM][NPR];
  FTYPE **alpha;
  struct of_state q;
  int i,j;
  FILE *out;
  FTYPE X[NDIM],V[NDIM];
  int myii,myjj,mykk;


  coord(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,X);
  bl_coord(X,V);


  if (get_state(pr0, ptrgeom, &q) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 2);
  if (primtoU(UNOTHING,pr0, &q, ptrgeom, Ustart) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 2);
  
  mhd_calc(pr0,0, ptrgeom, &q, mhd[0]);
  mhd_calc(pr0,1, ptrgeom, &q, mhd[1]);
  mhd_calc(pr0,2, ptrgeom, &q, mhd[2]);
  mhd_calc(pr0,3, ptrgeom, &q, mhd[3]);
  primtoflux(UNOTHING,pr0, &q, 1,ptrgeom, flux[0]);
  primtoflux(UNOTHING,pr0, &q, 2,ptrgeom, flux[1]);
  primtoflux(UNOTHING,pr0, &q, 3,ptrgeom, flux[2]);
  source(pr0, ptrgeom, dUcomp, dU);
  bsq_calc(pr0, ptrgeom, &bsq);
  
  
  alpha=dmatrix(1, 5, 1, 5);

  if (dudp_calc(EVOLVENOENTROPY,pr0, &q, ptrgeom, alpha) >= 1)
    FAILSTATEMENT("utoprim.c:usrfun()", "dudp_calc()", 1);
  // more general normalization
  // copied from usrfun
#if(NORMMETHOD==0)
  normf=1.0/U_target[RHO];
#elif(NORMMETHOD==1)
  normf=0.0;
  numnormterms=0;
  for (j = 0; j < INVERTNPR; j++){
    for (k = 0; k < INVERTNPR; k++){
      if(alpha[j + 1][k + 1]>NUMEPSILON){
	normf+=fabs(alpha[j + 1][k + 1]);
	numnormterms++;
      }
    }
  }
  normf=(FTYPE)(numnormterms)/(normf); // (i.e. inverse of average)
#endif    
  for (j = 0; j < INVERTNPR; j++)
    for (k = 0; k < INVERTNPR; k++)
      alpha[j + 1][k + 1] *= (normf);



  out=fopen("myconsts.txt","wt");  
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",X[1],X[2],X[3],V[1],V[2],V[3],gam,a,hslope,R0);
  fclose(out);

  
  out=fopen("myTud.txt","wt");
  for(i=0;i<NDIM;i++){
    for(j=0;j<NDIM;j++){
      fprintf(out,"%21.15g ",mhd[i][j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);
  
  out=fopen("myF.txt","wt");
  for(i=0;i<=2;i++){ // mathematica rows
    for(j=0;j<NPR;j++){ // mathematica columns
      fprintf(out, "%21.15g ",flux[i][j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);
  
  out=fopen("myub.txt","wt");  
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",q.ucon[0],q.ucon[1],q.ucon[2],q.ucon[3]);
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",q.ucov[0],q.ucov[1],q.ucov[2],q.ucov[3]);
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",q.bcon[0],q.bcon[1],q.bcon[2],q.bcon[3]);
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",q.bcov[0],q.bcov[1],q.bcov[2],q.bcov[3]);
  fclose(out);

  out=fopen("mygnbsq.txt","wt");  
  fprintf(out,"%21.15g %21.15g %21.15g\n",ptrgeom->g,normf,bsq);
  fclose(out);
  
  out=fopen("mypU.txt","wt");  
  PLOOP(k) fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",pr0[k],Ustart[k],U_target[k],dU[k]);
  fclose(out);

  
  out=fopen("alpha.txt","wt");  
  for(i=1;i<=5;i++){
    for(j=1;j<=5;j++){
      fprintf(out, "%21.15g ",alpha[i][j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);

#if(MCOORD!=CARTMINKMETRIC)
  myii=icurr;
  myjj=jcurr;
  mykk=kcurr;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif
  
   
  out=fopen("myconn.txt","wt");  
  for(i=0;i<NDIM;i++){
    for(j=0;j<NDIM;j++){
      for(k=0;k<NDIM;k++){
	fprintf(out, "%21.15g ",conn[myii][myjj][mykk][i][j][k]);
      }
      fprintf(out,"\n");
    }
  }
  fclose(out);
  
  out=fopen("mygcon.txt","wt");  

  for(i=0;i<NDIM;i++){
    for(j=0;j<NDIM;j++){
      fprintf(out, "%21.15g ",gcon[myii][myjj][mykk][ptrgeom->p][i][j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);

  out=fopen("mygcov.txt","wt");  

  for(i=0;i<NDIM;i++){
    for(j=0;j<NDIM;j++){
      fprintf(out, "%21.15g ",gcov[myii][myjj][mykk][ptrgeom->p][i][j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);

  return(0);
}
