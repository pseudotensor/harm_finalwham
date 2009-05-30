
#include "u2p_defs.h"

// whether to turn on jon's modificatoins to Scott's general Newton's method
#define JONGENNEWTSTUFF 1

#define CRAZYNEWCHECK 0

// actual number of dimensions in method (choose method = 1D or 2D).  Just a high-level choice.  Core bits of code don't use this
// 1: 1-D E'(W')  (really has 3 roots, but seemingly always positive-most root found AND correct root)
// 2: 2-D E'(W',v^2) and P^2(W',v^2) (unsure how many roots)
// 3: 1-D E'(W') Jon's method without intermediate v^2 and using \tilde{u}^2 instead of v^2 -- avoids catastrophic cancellation (still has 3 roots)
#define WHICHHOTINVERTER 3

// whether to do final calculation of primitives from v^2 or \tilde{u}^2
// 0: from \tilde{u}^2 -- avoids catastrophic cancellation at high Lorentz factors
// 1: from v^2 -- requires computing \gamma from 1/sqrt(1-v^2) which has problems when v^2\sim 1
#define PRIMFROMVSQ 0


// 0: original 1-D but with E' (W and v^2, but can give nan's when Wp<0 slightly -- maybe fixed with new normalization?) (apparently also not giving unique solution)
// 1: Jon's simple 1-D E' (no use of v^2, but obtains an inconsistent v^i[P^2,E'])  (apparently also not giving unique solution)
//    Uses same as WHICHHOTINVERTER==3
// 2: simple 1-D P^2 (no use of v^2, obtains a consistent v^i[P^2]) (but 4 roots and never sure which one is right -- not always positive root.  If D<0, then can be negative root)
// 3: 1-D P^\alpha error on all 3 components of P^i so there is a unique answer (NOT DONE YET)
#define WHICHCOLDINVERTER 2


#define METHODTYPE  3
#define NEWT_DIM 2 // actually maximum # of dimensions
#define USE_LINE_SEARCH -2

#define MAX_NEWT_RETRIES 0    /*Max. # of retries of N-R procedure, while increasing guess for W by *10 after each time*/

/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
static FTYPE Bsq,QdotBsq,QdotB,Qtsq,Qdotn,Qdotnp,D ;
static FTYPE wglobal;
static int geomtype; // global choice of EOM (EOMTYPE)
static int *glpflag; // global pflag for local file

FTYPE Qcovorig[NDIM],Qconorig[NDIM],Qdotnorig,Qsqorig,Qtsqorig;
FTYPE Qtsqnew;


// pointer to the correct validate_x function
static void (*ptr_validate_x)(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM] );


// Declarations: 
static FTYPE vsq_calc(FTYPE W);


#include "u2p_util.c"
#include "utoprim_jon_eos.c"


/******************************************************************
  
   Driver for new prim. var. solver.  The driver just translates
   between the two sets of definitions for U and P. 

******************************************************************/

int Utoprim_jon_nonrelcompat_inputnorestmass(int eomtype, FTYPE U[NPR], struct of_geom *ptrgeom,  int *lpflag,  FTYPE prim[NPR])
{

  static int Utoprim_new_body(FTYPE U[], struct of_geom *ptrgeom,  FTYPE prim[]);
  FTYPE U_tmp[NPR], U_tmp2[NPR], prim_tmp[NPR];
  int i, j, k,ret; 
  int pl,pl2;
  FTYPE gcov[NDIM][NDIM], gcon[NDIM][NDIM], alpha;
  const int ltrace = 0;



#if(!OPTIMIZED)
  // make sure nan and inf's are not in initial guess or U
  PLOOP(pl){
    if(!isfinite(prim[pl])){
      dualfprintf(fail_file,"Guess for p is NAN\n");
      PLOOP(pl2) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl2,prim[pl2]);
      *glpflag=  UTOPRIMFAILNANGUESS;
      myexit(876);
    }
    if(!isfinite(U[pl])){
      dualfprintf(fail_file,"Value of U is NAN\n");
      PLOOP(pl2) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl2,U[pl2]);
      *glpflag=  UTOPRIMFAILNANGUESS;
      myexit(877);
    }
  }
#endif


  lntries=0;


  // assign global int pointer to lpflag pointer
  glpflag=lpflag;
  geomtype=eomtype;


  // EOS functions used during inversion
  pickeos_eomtype(eomtype);


#if( WHICHVEL != VELREL4 ) 
  fprintf(stderr,"Utoprim_2d() Not implemented for WHICHVEL = %d \n", WHICHVEL );
  return(1); 
#endif

  /* First update the primitive B-fields */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  /* Set the geometry variables: */
  alpha = 1.0/sqrt(-ptrgeom->gcon[0][0]);
 
 
  /* Calculate the transform the CONSERVATIVE variables into the new system */
  for( i = RHO; i <= BCON3; i++ ) {
    U_tmp[i] = alpha * U[i] ;
  }
 
#if(CRAZYDEBUG)
  if(ptrgeom->i==0 && ptrgeom->j==31 && nstep==9 && steppart==2){
    PLOOP(k) dualfprintf(fail_file,"U_tmp[%d]=%21.15g\n",k,U_tmp[k]);
  }
#endif

 
  /* Calculate the transform the PRIMITIVE variables into the new system */
  for( i = 0; i < BCON1; i++ ) {
    prim_tmp[i] = prim[i];
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    prim_tmp[i] = alpha*prim[i];
  }

  ret = Utoprim_new_body(U_tmp, ptrgeom, prim_tmp);

  /* Check conservative variable transformation: */
#if(!OPTIMIZED)
  if( ltrace ) {
    for(i = 0; i < NDIM; i++ ) {
      for(j = 0; j < NDIM; j++ ) {
	gcov[i][j] = ptrgeom->gcov[i][j];
	gcon[i][j] = ptrgeom->gcon[i][j];
	dualfprintf(fail_file,"gcov,gcon %d %d = %21.15g %21.15g \n", i, j, gcov[i][j], gcon[i][j]);
      }
    }
    dualfprintf(fail_file,"gdet = %21.15g \n", ptrgeom->g);
    //    primtoU_g( prim_tmp, gcov, gcon, U_tmp2 ); 
    //    for( i = 0; i < NPR; i++ ) {
    //      dualfprintf(fail_file, "Utoprim_1d(): Utmp1[%d] = %21.15g , Utmp2[%d] = %21.15g , dUtmp[%d] = %21.15g \n", 
    //	       i, U_tmp[i], i, U_tmp2[i], i, fabs( (U_tmp[i]-U_tmp2[i]) / ( (U_tmp2[i]!=0.) ? U_tmp2[i] : 1. ) )  ); 
    //    }
  }
#endif

  /* Transform new primitive variables back : */ 
  if( ret == 0 ) {
    for( i = 0; i < BCON1; i++ ) {
      prim[i] = prim_tmp[i];
    }
  }


  // revert EOS for EOMTYPE
  pickeos_eomtype(EOMTYPE);

  return( 0 ) ;

}


/**********************************************************************************

attempt an inversion from U to prim using the initial guess
prim.

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu  |
             \  alpha B^i        /



             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \  alpha B^i   /


*glpflag:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence (occurrence of "nan" or "+/-inf" ;
	     
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
	     3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: utsq<0 or utsq > UTSQ_TOO_BIG   with new  W;
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(FTYPE U[NPR], struct of_geom *ptrgeom,  FTYPE prim[NPR])
{
  int compute_setup_quantities(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom,  FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov);
  FTYPE Qconp[NDIM],Qcovp[NDIM];
  FTYPE Qtcon[4],Bcon[4],Bcov[4];
  FTYPE primother[NPR];
  FTYPE W_last,W,Wp_last,Wp0,Wp,Wp2,Wp3,Wpother;
  int i,j, retval0, retval, retval2, retval3 ;
  int retvalother,retvalother2;
   
  static int find_root_2D_gen(FTYPE x0, FTYPE *xnew) ;
  static int find_root_1D_gen(FTYPE x0, FTYPE *xnew) ;
  static int find_root_1D_gen_scn(FTYPE x0, FTYPE *xnew);
  static int find_root_1D_gen_Eprime(FTYPE x0, FTYPE *xnew);
  static int find_root_1D_gen_Psq(FTYPE x0, FTYPE *xnew);
  static int find_root_3D_gen_Palpha(FTYPE x0, FTYPE *xnew);


  const int ltrace = 0;
  const int ltrace2 = 0;

  int forcefree_inversion(struct of_geom *ptrgeom, FTYPE *Qtcon, FTYPE Bsq, FTYPE *Bcon, FTYPE *Bcov, FTYPE Qtsq, FTYPE *U, FTYPE *prim);
  int set_guess_Wp(FTYPE *prim, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last, FTYPE *wglobal);
  int check_Wp(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp_last, FTYPE Wp, int retval);
  int Wp2prim(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, int retval);
  static void check_on_inversion(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, int retval);



  /* TEMPORARY */
  /*
//    primtoU_g(prim,gcov,gcon,U) ;
  */

#if(!OPTIMIZED)
  if( ltrace ) {
    for(i=0;i<8;i++) dualfprintf(fail_file,"%d %21.15g %21.15g\n",i,prim[i],U[i]) ;
    
  }
#endif
  // Assume ok initially:
  *glpflag= UTOPRIMNOFAIL;



  // compute quantities necessary for inversion in a certain frame  
  compute_setup_quantities(prim, U, ptrgeom, Qtcon, Bcon, Bcov);


#if(0)
  if(ptrgeom->i==1 && ptrgeom->j==18){
    DLOOPA(i) dualfprintf(fail_file,"Qtcon[%d]=%21.15g Bcon[%d]=%21.15g Bcov[%d]=%21.15g\n",i,Qtcon[i],i,Bcon[i],i,Bcov[i]);
    dualfprintf(fail_file,"Bsq=%21.15g QdotB=%21.15g QdotBsq=%21.15g Qtsq=%21.15g Qdotnp=%21.15g D=%21.15g\n",Bsq,QdotB,QdotBsq,Qtsq,Qdotnp,D);
  }
#endif

  ////////////////////////////////
  //
  // CHECK WHICH EOM TYPE and proceed accordingly
  //
  ////////////////////////////////

  if(geomtype==EOMFFDE){
    // then perform the analytic inversion with Lorentz factor limit
    if(forcefree_inversion(ptrgeom, Qtcon, Bsq, Bcon, Bcov, Qtsq, U, prim)) return(1);
    else return(0) ;
  }




  // SETUP ITERATIVE METHODS (good for GRMHD or cold GRMHD)
  set_guess_Wp(prim, ptrgeom, &W_last, &Wp_last,&wglobal);






  ////////////////////////////////////////
  ///////////////////////////////////////
  //
  //  Obtain Wp
  //
  ////////////////////////////////////// 
  ////////////////////////////////////// 





  /////////////////////////////////////////////////////////////////////////
  //
  //
  // HOT GRMHD INVERTERS
  //
  /////////////////////////////////////////////////////////////////////////
  if(geomtype==EOMGRMHD){

    // METHOD specific:
#if(WHICHHOTINVERTER==3)
    retval = find_root_1D_gen_Eprime(Wp_last, &Wp);
#elif(WHICHHOTINVERTER==2)
    retval = find_root_2D_gen(Wp_last, &Wp);
#elif(WHICHHOTINVERTER==1)
    retval = find_root_1D_gen(Wp_last, &Wp);
#endif


#if(0)
    retval = find_root_2D_gen(Wp_last, &Wp);
    retvalother = find_root_1D_gen(Wp_last, &Wpother);
    retvalother2=find_root_1D_gen_scn(W_last, &W);

    //  Wp=Wpother;
    if(retval!=retvalother || (fabs(Wp-Wpother)/(Wp+Wpother)>1E-13)){
      dualfprintf(fail_file,"2D/1D: i=%d j=%d :: %d %d %d :: Wp_last=%21.15g  Wp=%21.15g Wpother=%21.15g\n",ptrgeom->i,ptrgeom->j,retval,retvalother,retvalother2,Wp_last,Wp,Wpother);
      dualfprintf(fail_file,"scn W_last=%g W=%g D=%g  Wp_last=%g Wp=%g\n",W_last,W,D,W_last-D,W-D);
      dualfprintf(fail_file,"rho0=%g u=%g p=%g gamma=%g gammasq=%g w=%g utsq=%g\n",rho0,u,p,gamma,gammasq,w,utsq);
    }
#endif // end comparing SCN and Jon method


    check_on_inversion(prim, U, ptrgeom, Wp, Qtcon, Bcon, Bcov,retval);


    //  Check if solution was found
    if(check_Wp(prim, U, ptrgeom, Wp_last, Wp, retval)>0) return(retval);


    // find solution
    if(Wp2prim(prim, U, ptrgeom, Wp, Qtcon, Bcon, Bcov, retval)>0) return(retval);


  } // end if hot GRMHD





  /////////////////////////////////////////////////////////////////////////
  //
  //
  // COLD GRMHD INVERTERS
  //
  /////////////////////////////////////////////////////////////////////////
  if(geomtype==EOMCOLDGRMHD){



#if(WHICHCOLDINVERTER==0)
    retval0 = find_root_1D_gen(Wp_last, &Wp0);

    check_on_inversion(prim, U, ptrgeom, Wp, Qtcon, Bcon, Bcov,retval);

    //  Check if solution was found
    check_Wp(prim, U, ptrgeom, Wp_last, Wp0, retval0);

    //  if(*glpflag!=0){
    //    dualfprintf(fail_file,"E'[W,v^2] equation gave bad answer: Wp=%21.15g Qdotnp=%21.15g Qdotn=%21.15g D=%21.15g\n",Wp0,Qdotnp,Qdotn,D);
    //    *glpflag=0; // let E' try again
    //  }
    //  else{
    //    dualfprintf(fail_file,"E'[W',v^2] equation gave good answer: Wp=%21.15g  Qdotnp=%21.15g\n",Wp0,Qdotnp);
    //  }

    retval=retval0;
    Wp=Wp0;
#endif



#if(WHICHCOLDINVERTER==1)

    // Do inversion using E'[W'] equation
    retval = find_root_1D_gen_Eprime(Wp_last, &Wp);

    check_on_inversion(prim, U, ptrgeom, Wp, Qtcon, Bcon, Bcov,retval);

    //  Check if solution was found
    check_Wp(prim, U, ptrgeom, Wp_last, Wp, retval);


    //  if(*glpflag!=0){
    //    dualfprintf(fail_file,"E'[W'] equation gave bad answer: Wp=%21.15g  Bsq=%21.15g Qdotnp=%21.15g Qdotn=%21.15g D=%21.15g\n",Wp,Bsq,Qdotnp,Qdotn,D);

    //    for(i=0;i<4;i++) Qcovp[i] = U[QCOV0+i] ;
    //    raise_g(Qcovp,ptrgeom->gcon,Qconp) ;  

    //    for(i=0;i<4;i++) dualfprintf(fail_file,"%d Qcovp=%21.15g Qconp=%21.15g\n",i,Qcovp[i],Qconp[i]);


    //    *glpflag=0; // let momentum try again
    //  }
    //  else{
    //    dualfprintf(fail_file,"Wp=%21.15g Bsq=%21.15g QdotB=%21.15g QdotBsq=%21.15g Qtsq=%21.15g Qdotnp=%21.15g D=%21.15g\n",Wp,Bsq,QdotB,QdotBsq,Qtsq,Qdotnp,D);
    //    dualfprintf(fail_file,"E'[W'] equation gave good answer: Wp=%21.15g Qdotnp=%21.15g\n",Wp,Qdotnp);
    //  }

  

#endif


#if(WHICHCOLDINVERTER==2)
    // Do inversion using P^2[W'] equation

#if(CRAZYDEBUG)
    if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"Wp_last=%21.15g\n",Wp_last);
    }
#endif

    retval2 = find_root_1D_gen_Psq(Wp_last, &Wp2);

    check_on_inversion(prim, U, ptrgeom, Wp, Qtcon, Bcon, Bcov,retval);

    //  Check if solution was found
    check_Wp(prim, U, ptrgeom, Wp_last, Wp2, retval2);


    //  if(*glpflag!=0){
    //    dualfprintf(fail_file,"P^2[W'] equation gave bad answer: Wp=%21.15g Qdotnp=%21.15g Qdotn=%21.15g\n",Wp2,Qdotnp,Qdotn);
    //    myexit(0);
    //  }
    //  else{
    //    dualfprintf(fail_file,"P^2[W'] equation gave good answer: Wp=%21.15g Qdotnp=%21.15g\n",Wp2,Qdotnp);
    //  }


    //  Wp = (Wp>Wp2) ? Wp2 : Wp;
    retval=retval2;
    Wp=Wp2;

#endif

#if(WHICHCOLDINVERTER==3)
    // Do inversion using P^\alpha[W'] equation

    retval3 = find_root_3D_gen_Palpha(Wp_last, &Wp3);

    check_on_inversion(prim, U, ptrgeom, Wp, Qtcon, Bcon, Bcov,retval);

    //  Check if solution was found
    check_Wp(prim, U, ptrgeom, Wp_last, Wp3, retval3);


    //  if(*glpflag!=0){
    //    dualfprintf(fail_file,"P^2[W'] equation gave bad answer: Wp=%21.15g Qdotnp=%21.15g Qdotn=%21.15g\n",Wp2,Qdotnp,Qdotn);
    //    myexit(0);
    //  }
    //  else{
    //    dualfprintf(fail_file,"P^2[W'] equation gave good answer: Wp=%21.15g Qdotnp=%21.15g\n",Wp2,Qdotnp);
    //  }


    //  Wp = (Wp>Wp3) ? Wp3 : Wp;
    retval=retval3;
    Wp=Wp3;

#endif


    // find solution
    if(Wp2prim(prim, U, ptrgeom, Wp, Qtcon, Bcon, Bcov, retval)>0) return(retval);
  


  } // end if geomtype==EOMCOLDGRMHD






  /* done! */
  return(retval) ;

}




//returns global quanities Bsq,QdotBsq,Qtsq,Qdotn,Qdotnp,D
// returns other local quantities
int compute_setup_quantities(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov) // returns global quantities: D, Bsq, Qtcon
{
  int i,j;
  FTYPE Qcov[4],Qcon[4],Qcovp[4],Qconp[4],ncov[4],ncon[4],Qsq,Qtcov[4];
  static void raise_g(FTYPE vcov[], FTYPE gcon[][4], FTYPE vcon[]);
  static void lower_g(FTYPE vcon[], FTYPE gcov[][4], FTYPE vcov[]);
  static void ncov_calc(FTYPE gcon[][4],FTYPE ncov[]) ;



  //  for(i=0;i<NPR;i++) dualfprintf(fail_file,"U[%d]=%g\n",i,U[i]);
  //  exit(0);

  D = U[RHO] ;


  // get B^i
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  // get B_\mu
  lower_g(Bcon,ptrgeom->gcov,Bcov) ;

  // get B^2
  Bsq = 0. ;
  for(i=1;i<4;i++) Bsq += Bcon[i]*Bcov[i] ;

  // get \eta^\mu and \eta_\mu
  ncov_calc(ptrgeom->gcon,ncov) ;
  raise_g(ncov,ptrgeom->gcon,ncon);

  // -P_\alpha (P_t has D term inside it)
  // Q_\nu = \alpha \{T^t_\nu\}  = -T.\eta
  
  // Q'_\nu = Q_\nu + \delta^t_\nu D
  for(i=0;i<4;i++) Qcovp[i] = U[QCOV0+i] ;
  raise_g(Qcovp,ptrgeom->gcon,Qconp) ;  



  // P^2, which can be written without use of T^t_t
  //  Qsq = 0. ;
  //for(i=0;i<4;i++) Qsq += Qcovp[i]*Qconp[i] ;

  // P^2 = Qtilde^2 can be written without T^t_t
  // Notice U !=Qcon, but only a problem for time term that is not used
  // T^t_t not really needed (cancels in the end)
  //  DLOOP(j) Qtcov[j]=0.0;
  //  DSLOOP(j,i) Qtcov[j] = U[QCOV0+i]*(delta(i,j)+ncon[i]*ncov[j]) ;
  //  raise_g(Qtcov,ptrgeom->gcon,Qtcon) ;

  // P^\alpha = Qtcon^\alpha
  DLOOPA(j) Qtcon[j]=0.0;
  SSLOOP(j,i) Qtcon[j] += U[QCOV0+i]*(ptrgeom->gcon[i][j] + ncon[i]*ncon[j]) ;
  //  for(i=0;i<4;i++){
  //    for(j=1;j<4;j++){
  //      Qtcon[i] = U[QCOV0+j]*(ptrgeom->gcon[j][i]+ncon[j]*ncon[i]);
  //    }
  //  }
  lower_g(Qtcon,ptrgeom->gcov,Qtcov);
  Qtsq = 0. ;
  for(i=0;i<4;i++) Qtsq += Qtcov[i]*Qtcon[i] ;

#if(0)
  DLOOPA(j) dualfprintf(fail_file,"Qtcon[%d]=%21.15g Qtcov[%d]=%21.15g\n",j,Qtcon[j],j,Qtcov[j]);
#endif


  // Qcon/Qcov  = -\eta.T = \alpha T^t_\nu :: where T^t_t *includes* rest-mass
  // avoid using this Q if possible to avoid large errors when v or u is small
  // put back in the D
  for(i=0;i<4;i++){
    if(i==0) Qcov[i]=Qcovp[i]-D;
    else Qcov[i]=Qcovp[i];
  }
  raise_g(Qcov,ptrgeom->gcon,Qcon) ;  


  


  //  dualfprintf(fail_file,"D=%g\n",D);
  //  for(i=0;i<4;i++) dualfprintf(fail_file,"Qcon[%d]=%g Qcov[%d]=%g\n",i,Qcon[i],i,Qcov[i]);
  //  exit(0);



  // P.B (P_t not used because B^t=0, so ok to use Qcov)
  QdotB = 0. ;
  for(i=1;i<4;i++) QdotB += Qcov[i]*Bcon[i] ;
  QdotBsq = QdotB*QdotB ;

  // Energy
  // Qdotnp=-E'=-E+D automatically
  // GODMARK: need to store that metric term to have gravity taken into account in non-rel limit
  Qdotnp = Qconp[0]*ncov[0] + D*(1.0-1.0/fabs(ncov[TT])) ; // -Qdotn-W = -Qdotnp-Wp
  Qdotn  =  Qcon[0]*ncov[0] ;


#if(0)
  // used to check on old vs. new method
  Qtsqnew=Qtsq;
  //  dualfprintf(fail_file,"1: Qtsq=%21.15g\n",Qtsq);
  for(i=0;i<4;i++) Qcovorig[i] = U[QCOV0+i] ;
  Qcovorig[TT] -=D;
  raise_g(Qcovorig,ptrgeom->gcon,Qconorig) ;
  Qdotnorig = Qconorig[0]*ncov[0] ;
  Qsqorig = 0. ;
  for(i=0;i<4;i++) Qsqorig += Qcovorig[i]*Qconorig[i] ;
  Qtsqorig = Qsqorig + Qdotnorig*Qdotnorig ;
  //  dualfprintf(fail_file,"2: Qtsq=%21.15g\n",Qtsq);
  Qtsq=Qtsqorig;
#endif


  // P^2 (since D already removed from E, then cancelling terms are small)
  // so no problem with catastrophic cancellation since involves energy term
  //  Qtsq = Qsq + Qdotnp*Qdotnp ; // P^2 in eta frame


  //  dualfprintf(fail_file,"Qdotn=%g Qtsq=%g\n",Qdotn,Qtsq);
  //exit(0);


  // debug Qsq
#if(CRAZYDEBUG)
  if(ptrgeom->i==0 && ptrgeom->j==31 && nstep==9 && steppart==2){
    Qsq = 0. ;
    for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;
    for(i=0;i<4;i++) dualfprintf(fail_file,"Qcon[%d]=%21.15g Qcov[%d]=%21.15g Qsq=%21.15g Qtsq=%21.15g\n",i,Qcon[i],i,Qcov[i],Qsq,Qtsq);
  }
#endif


  return(0);

}



////////////////////////////////
//
// SETUP ITERATIVE METHODS (good for GRMHD or cold GRMHD)
//
////////////////////////////////
int set_guess_Wp(FTYPE *prim, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last, FTYPE *wglobal)
{
  FTYPE u,p;
  FTYPE utsq;
  FTYPE gammasq,gamma,rho0,w;
  FTYPE bsq;
  int verify_Wlast(FTYPE u, FTYPE p, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last);
  int i,j;
  static void check_utsq_fail(FTYPE Wp);
  static int coldgrhd(FTYPE Qtsq, FTYPE D, FTYPE *Wp);
  FTYPE Wpcoldhd;


  /* calculate W from last timestep and use 
     for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += ptrgeom->gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  if( (utsq < 0.) && (fabs(utsq) < MAXNEGUTSQ) ) { 
    utsq = 0.0;
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    if( debugfail>=2 ) dualfprintf(fail_file,"Utoprim_new(): utsq < 0 in utoprim_1d attempt, utsq = %21.15g \n", utsq) ;
    *glpflag= UTOPRIMFAILCONVGUESSUTSQ; // guess failure actually
    dualfprintf(fail_file,"utsq=%21.15g\n",utsq);

    // check on v^2>1 failure
    check_utsq_fail(-1E30);

    return(1) ;
  }

  gammasq = 1. + utsq ;
  gamma  = sqrt(gammasq);
	
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = D / gamma ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0,u) ;
  w = rho0 + u + p ;

  // need b^2 to normalize W since error in W limited by b^2
  bsq_calc(prim, ptrgeom, &bsq);


  // W'=W-D (D removed from W)
  *Wp_last = D*utsq/(1.0+gamma) + (u+p)*gammasq;
  *W_last = *Wp_last + D;


  //  *wglobal=w+bsq; // used to validate normalized version of W

  if(JONGENNEWTSTUFF&&((w+bsq)*gammasq>=GAMMASQCHECKRESID)){
    if(coldgrhd(Qtsq, D, &Wpcoldhd)) Wpcoldhd=*Wp_last; // if failed to get coldGRHD Wp, then just use Wp_last
    
    //  *wglobal=MAX(MAX(Wpcoldhd,*Wp_last),w+bsq); // used to validate normalized version of W -- accurate for cold GRHD
    *wglobal=1E-5*MAX(MAX(Wpcoldhd,*Wp_last),w+bsq); // used to validate normalized version of W -- accurate for cold GRHD
    //*wglobal=w+bsq;
    
    *Wp_last=MAX(MAX(2.0*Wpcoldhd,*Wp_last),w+bsq); // use cold HD to upgrade Wp if needed
  }
  else{
    *wglobal=w+bsq;
  }


  //  if((nstep==4)&&(ptrgeom->i==0)&&(ptrgeom->j==47)){
  // dualfprintf(fail_file,"utsq=%g p=%g w=%g W_last=%g Wp_last=%g\n",utsq,p,w,*W_last,*Wp_last);
    //    exit(0);
  //  }

#if(CRAZYDEBUG)
  if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"pre verify\n");
    dualfprintf(fail_file,"rho0=%21.15g u=%21.15g p=%21.15g\n",rho0,u,p);
    
    dualfprintf(fail_file,"W_last=%21.15g Wp_last=%21.15g gammasq=%21.15g w=%21.15g\n",*W_last,*Wp_last,gammasq,w);
  }
#endif


  verify_Wlast(u, p, ptrgeom, W_last, Wp_last); // W_last and Wp_last already pointers

#if(CRAZYDEBUG)
  if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"post verify\n");
    
    dualfprintf(fail_file,"W_last=%21.15g Wp_last=%21.15g gammasq=%21.15g w=%21.15g\n",*W_last,*Wp_last,gammasq,w);
  }
#endif




  return(0);


}


///////////////////////////////////////
//
// Good for GRMHD or cold GRMHD
//
// Make sure that W is large enough so that v^2 < 1 : 
//
// GODMARK: apparently this is necessary for 1D method, otherwise it can't find solution (2D method seems to be ok)
//
////////////////////////////////////// 
int verify_Wlast(FTYPE u, FTYPE p, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last)
{
  int i_increase;
  static int vsqgtr1(FTYPE W,FTYPE Bsq,FTYPE QdotBsq, FTYPE Qtsq);



  if(*W_last<-D || *Wp_last<0){
    *W_last = -D;
    *Wp_last = 0;
    if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"W_last=%21.15g Wp_last=%21.15g\n",*W_last,*Wp_last);
    }
  }


  if(1){
    i_increase = 0;
    while( vsqgtr1(*W_last,Bsq,QdotBsq,Qtsq) && (i_increase < 10) ) {
    
      //    dualfprintf(fail_file,"bumped up %d %d %g %g\n",ptrgeom->i,ptrgeom->j,*W_last,*Wp_last);


      *Wp_last*= 10.;
      *W_last = *Wp_last+D;

      //    dualfprintf(fail_file,"bumped up %d %d\n",ptrgeom->i,ptrgeom->j);


      i_increase++;
#if(!OPTIMIZED)
      dualfprintf(fail_file,"badval :  W = %21.15g, i_increase = %d \n", *W_last, i_increase); 
#endif
    }
#if(!OPTIMIZED)
    if( i_increase >= 10 ) { 
      dualfprintf(fail_file,"i_increase is too large, i_increase = %d , W = %21.15g \n", i_increase, *W_last);
    }
#endif
	  

#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file,"u = %21.15g,  p = %21.15g,  Bsq = %21.15g, Qtsq = %21.15g \n",u,p,Bsq,Qtsq);
      dualfprintf(fail_file,"Bcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcon[0],Bcon[1],Bcon[2],Bcon[3]);
      dualfprintf(fail_file,"Bcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcov[0],Bcov[1],Bcov[2],Bcov[3]);
      dualfprintf(fail_file,"Qcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcon[0],Qcon[1],Qcon[2],Qcon[3]);
      dualfprintf(fail_file,"Qcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcov[0],Qcov[1],Qcov[2],Qcov[3]);
      dualfprintf(fail_file,"call find_root\n") ; 	
    
    }
#endif
  }


  return(0);

}




///////////////////////////////////////
//
//  Check if solution was found
//
////////////////////////////////////// 
int check_Wp(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp_last, FTYPE Wp, int retval)
{
  int i;
  FTYPE Wtest;

	  


  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (Wp == FAIL_VAL) ) {
    if( debugfail>=2 ) {
      dualfprintf(fail_file, "Failed to find a prim. var. solution!! %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g \n",Wp_last,Bsq,QdotBsq,Qdotn,D,Qtsq);
      dualfprintf(fail_file, "Utoprim_new_body(): bad newt failure, nstep,steppart :: t,i,j, p[0-7], U[0-7] = %ld %d :: %21.15g %d %d ", nstep,steppart,t, ptrgeom->i, ptrgeom->j );  
      for( i = 0 ; i < NPR; i++ ) {
	dualfprintf(fail_file, "%21.15g ", prim[i]);
      }
      for( i = 0 ; i < NPR; i++ ) {
	dualfprintf(fail_file, "%21.15g ", U[i]);
      }
      dualfprintf(fail_file, "\n");
    }      
    *glpflag= retval*100+1;// related to UTOPRIMFAILCONVRET
    dualfprintf(fail_file,"Got here Wp_last=%21.15g Wp=%21.15g retval=%d\n",Wp_last,Wp,retval);
    return(retval);
  }
  else{
    Wtest=Wp/wglobal; // normalize to old densities

    // GODMARK
    //      if(Wtest<=0. && fabs(Wtest)>1E-4){
    if(Wtest<=0.0){
      Wp = SMALL;
      return(retval) ;
    }
    else if(Wtest <= 0. || Wtest > GAMMASQ_TOO_BIG) {



      //if(Wtest <= 0.) {
      if( debugfail>=2 ) {
		dualfprintf(fail_file,"Wtest failure %21.15g \n",Wtest) ;
	dualfprintf(fail_file, "Utoprim_new_body(): Wtest<0 or Wtest=toobig failure, t,i,j, p[0-7], U[0-7] = %21.15g %d %d ", t, ptrgeom->i, ptrgeom->j );  
	for( i = 0 ; i < NPR; i++ ) {
	  dualfprintf(fail_file, "%21.15g ", prim[i]);
	}
	for( i = 0 ; i < NPR; i++ ) {
	  dualfprintf(fail_file, "%21.15g ", U[i]);
	}
	dualfprintf(fail_file, "\n");
      }      

      *glpflag=  UTOPRIMFAILCONVW;
      return(retval) ;
    }
  }


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"(Wp,Wp_last,Bsq,Qtsq,QdotB,gammasq,Qdotn) %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
	    Wp,Wp_last,
	    Bsq,Qtsq,QdotB,gammasq,Qdotn) ;
    dualfprintf(fail_file,"done find_root\n") ;	
  }
  if( ltrace2 ) {
    dualfprintf(fail_file, "\n <--------- %21.15g %21.15g %21.15g %21.15g %21.15g  \n", Bsq,QdotBsq,Qdotn,D,Qtsq);
    
  }
#endif

  /*
  if( nstep == 3 ) { 
    if( ((ptrgeom->i==89) && (ptrgeom->j==88||ptrgeom->j==111)) || (ptrgeom->i==92&&(ptrgeom->j==86||ptrgeom->j==113)) 
	|| (ptrgeom->i==92&&(ptrgeom->j==86||ptrgeom->j==113) ) || (ptrgeom->i==94&&(ptrgeom->j==84||ptrgeom->j==115) ) 
	|| (ptrgeom->i==105&&(ptrgeom->j==84||ptrgeom->j==115) ) || (ptrgeom->i==107&&(ptrgeom->j==86||ptrgeom->j==113) ) 
	|| (ptrgeom->i==110&&(ptrgeom->j==88||ptrgeom->j==111) ) ) {
      dualfprintf(fail_file, "\n <--------- %21.15g %21.15g %21.15g %21.15g %21.15g  \n", Bsq,QdotBsq,Qdotn,D,Qtsq);
    }
  }
  */

  return(0);


}


// check on v^2>1 failure
static void check_utsq_fail(FTYPE Wp)
{

#if(0)
  if(steppart==1 && nstep==32 && icurr==0 &&  jcurr==19){
    dualfprintf(fail_file,"Qtsq=%21.15g Bsq=%21.15g D=%21.15g QdotB=%21.15g Wp=%21.15g\n",Qtsq,Bsq,D,QdotB,Wp);
    myexit(12245);
  }
#endif

}


static void check_on_inversion(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, int retval)
{
  int checki;

#if(CHECKONINVERSION)
  // then store quantities required to check up on inversion in mathematica
  checki=0;
  globalinv[checki++]=Qdotnp;
  globalinv[checki++]=Qtsq;
  //  globalinv[checki++]=Qtsqorig;
  globalinv[checki++]=Qtsqnew;
  globalinv[checki++]=Bsq;
  globalinv[checki++]=D;
  globalinv[checki++]=QdotB;
  globalinv[checki++]=Wp;
  globalinv[checki++]=Qtcon[1];
  globalinv[checki++]=Qtcon[2];
  globalinv[checki++]=Qtcon[3];
  globalinv[checki++]=Bcon[1];
  globalinv[checki++]=Bcon[2];
  globalinv[checki++]=Bcon[3];
  if(checki>NUMGLOBALINV){
    dualfprintf(fail_file,"Not enough memory for globalinv: checki=%d NUMGLOBALINV=%d\n",checki,NUMGLOBALINV);
  }
#endif
}



int Wp2prim(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, int retval)
{
  FTYPE W,vsq;
  FTYPE utsq;
  static FTYPE utsq_calc(FTYPE W);
  int i;
  FTYPE rho0;
  FTYPE gtmp,gamma,gammasq,w,wmrho0,p,u,tmpdiff;
  FTYPE wmrho0_compute_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma,FTYPE gammasq);
  FTYPE wmrho0_compute_utsq(FTYPE Wp, FTYPE D, FTYPE utsq, FTYPE gamma,FTYPE gammasq);
  static void check_utsq_fail(FTYPE Wp);
  static int validate_vsq(FTYPE vsqold,FTYPE *vsqnew);
  int retvsq;
  FTYPE pr0[NPR];
  int pl,pl2;

  //  if(icurr==1 && jcurr==10 && nstep==137 && steppart==3){
    //Wp = -0.513243850397847;
    //    Wp = -0.4043038500443759;
    //Wp = -0.0010560020881852734;
    //Wp =  0.001055708263718976;
    //dualfprintf(fail_file,"Testing\n");
    //}


#if(!OPTIMIZED)
  if(!isfinite(Wp)){
    PLOOP(pl) prim[pl]=pr0[pl]; // since might contaminate solution
    *glpflag=  UTOPRIMFAILCONVW;
    return(1);
  }
#endif





  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  // set backup of primitive
  PLOOP(pl) pr0[pl]=prim[pl];

  

  //  if(icurr==0 && jcurr==31 && nstep==51 && steppart==2){
    //    Wp=-0.21569833959001; // big error
    //nstep=51 stepart=2 :: i=0 j=31 :: lntries=15 lerrx= 9.64208793876686e-20
    //fdiff[1]=     1.10722611308861 ::     -1133.39047976655        57.67252740565


    //    Wp=-0.20494570836483;
    //nstep=51 stepart=2 :: i=0 j=31 :: lntries=15 lerrx= 9.64208793876686e-20
    //fdiff[1]=    0.903231963163783 ::     -1133.39047976655     -57.6261715957891

    //    Wp=-0.0008728880019961;


    //    Wp=0.0011219944118860;
    //nstep=51 stepart=2 :: i=0 j=31 :: lntries=15 lerrx= 9.64208793876686e-20
    //fdiff[3]=     15.1572386799252 ::      -106.20485235917      121.208463552046
  //}



  W=Wp+D;



#if(PRIMFROMVSQ==1)
  // Calculate vsq
  vsq = vsq_calc(W) ;

  retvsq=validate_vsq(vsq,&vsq);


  if( retvsq==1 || retvsq==2) {
      //(vsq >= 1.0) || (vsq<0.0) ) {
    if( debugfail>=2 ) { 
      dualfprintf(fail_file,"vsq failure:  vsq = %21.15g , W = %21.15g \n",vsq, W) ;
      dualfprintf(fail_file, "Utoprim_new_body(): utsq==bad failure, t,i,j, p[0-7], U[0-7] = %21.15g %d %d ", t, ptrgeom->i, ptrgeom->j );  
      for( i = 0 ; i < NPR; i++ ) {
	dualfprintf(fail_file, "%21.15g ", prim[i]);
      }
      for( i = 0 ; i < NPR; i++ ) {
	dualfprintf(fail_file, "%21.15g ", U[i]);
      }
      dualfprintf(fail_file, "\n");
    }      

    // check on v^2>1 failure
    check_utsq_fail(Wp);

    *glpflag=  UTOPRIMFAILCONVUTSQVERYBAD;
    return(retval) ;
  }
  else if(retvsq==3){
    // if( vsq > VSQ_TOO_BIG ) {

    // check on v^2>VSQ_TOO_BIG failure (not really a failure of convergence, but of confidence in the solution)
    check_utsq_fail(Wp);

    *glpflag=  UTOPRIMFAILCONVUTSQ;
    return(retval) ;
  }

  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  gammasq=gamma*gamma;
  // can compute utsq or gamma from Wp alone instead of using vsq
  utsq = gammasq*vsq;

  // GODMARK: can compute u+p directly from utsq, Wp, and D.
  wmrho0=wmrho0_compute_vsq(Wp,D,vsq,gamma,gammasq);


#elif(PRIMFROMVSQ==0)

  utsq=utsq_calc(W);

  gammasq=1.0+utsq;
  gamma=sqrt(gammasq);

  wmrho0=wmrho0_compute_utsq(Wp,D,utsq,gamma,gammasq);


#endif

  rho0 = D / gamma;

  //  w = W * (1. - vsq) ;
  // not used, but compute correctly anyways
  // w = \rho_0 + u + p
  w = rho0+wmrho0;

  p = pressure_wmrho0(rho0, wmrho0) ;
  //  u = w - (rho0 + p) ;
  u = wmrho0-p;

//  if((nstep==4)&&(ptrgeom->i==0)&&(ptrgeom->j==47)){
//  dualfprintf(fail_file,"Wp=%g W=%g D=%g vsq=%g gamma=%g wmrho0=%g p=%g u=%g\n",Wp,W,D,vsq,gamma,wmrho0,p,u);
    //    exit(0);
  //  }


  // GODMARK
  // fix checks for cold case
  if( (rho0 <= 0.) 
      || ((u < 0.)&&(geomtype==EOMGRMHD))
      ) { 
    if( debugfail>=2 ) {
      tmpdiff = w - rho0;
      dualfprintf(fail_file,
		  "rho or uu < 0 failure: rho,w,(w-rho),p,u  = %21.15g %21.15g %21.15g %21.15g %21.15g \n",
		  rho0,w,tmpdiff,p,u) ;
      dualfprintf(fail_file,
		  "rho or uu < 0 failure: gamma,utsq = %21.15g %21.15g  \n",  gamma, utsq) ;
    }
    if((rho0<=0.)&&(u>=0.)) *glpflag=  UTOPRIMFAILRHONEG;
    if((rho0>0.)&&(u<0.)) *glpflag= UTOPRIMFAILUNEG;
    if((rho0<=0.)&&(u<0.)) *glpflag= UTOPRIMFAILRHOUNEG;
    if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNNOTADJUSTED) return(retval) ; // else let assign -- used to check how bad failure is.
  }



  prim[RHO] = rho0 ;



  if(geomtype==EOMGRMHD){
    prim[UU] = u ;
  }
  else if(geomtype==EOMCOLDGRMHD){
    prim[UU] = 0;
    //  dualfprintf(fail_file,"u=%21.15g\n",u);
  }




  //  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn; // Qtcon already defined in good way

  // GODMARK: gamma should be computed from W and conserved quantities directly as \gamma[W] (see notes)
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;


  //  dualfprintf(fail_file,"%d %d : %g %g %g\n",ptrgeom->i,ptrgeom->j,vsq,W,Wp);
  //if(fabs(prim[U1])>1E-30)
  //  for(i=1;i<4;i++) dualfprintf(fail_file,"u[%d]= %g %g\n",i,Qtcon[i],prim[UTCON1+i-1]);
  //  dualfprintf(fail_file,"retval=%d\n",retval);


  // make sure nan and inf's don't pass through
  PLOOP(pl){
    if(!isfinite(prim[pl])){
      prim[pl]=pr0[pl]; // since might contaminate solution
      if(debugfail>=2)  PLOOP(pl2) dualfprintf(fail_file,"prim is nan: prim[%d]=%21.15g -> pr0[%d]=%21.15g\n",pl2,prim[pl2],pl2,pr0[pl2]);
      *glpflag=  UTOPRIMFAILNANRESULT;
    }
  }


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file," rho final = %21.15g ,  u final = %21.15g \n", rho0, u);
  }
#endif



  return(0);


}






// appears to give qualitatively similar, but still qualitatively different results than original phys.ffde.c code.  Gives more oscillatory results for torus force-free field.
int forcefree_inversion(struct of_geom *ptrgeom, FTYPE *Qtcon, FTYPE Bsq, FTYPE *Bcon, FTYPE *Bcov, FTYPE Qtsq, FTYPE *U, FTYPE *prim)
{
  // force-free variables
  FTYPE Qtconclean[NDIM],Qtcovclean[NDIM],Qtsqclean;
  FTYPE alphasq,gammamax,vcon[NDIM],vcov[NDIM],QtdotB,vconff[NDIM],vcovff[NDIM];
  static void lower_g(FTYPE vcon[], FTYPE gcov[][4], FTYPE vcov[]);
  int i;
  FTYPE vsq,gamma;
  FTYPE realBsq;

  // note that doing this means B.B!=Bsq
  Bsq+=SMALL; // in case Bsq=0 identically


  alphasq=GAMMAMAX*GAMMAMAX/(GAMMAMAX*GAMMAMAX-1.0);


  prim[RHO]=prim[UU]=0.0;


  // first 3-velocity
  vcon[TT]=0.0; // Qtcon[TT]=0
  for(i=1;i<4;i++) vcon[i] = Qtcon[i]/Bsq ;

  // first project out parallel velocity from Qtcon
  // note that \tilde{Q}.B = Q.B
  QtdotB=0.0;
  for(i=1;i<4;i++) QtdotB += Qtcon[i]*Bcov[i];

  Qtconclean[TT]=0.0;
  for(i=1;i<4;i++) Qtconclean[i] = Qtcon[i] - QtdotB*Bcon[i]/Bsq;

  // define new cleaned 3-velocity
  vconff[TT]=0.0;
  for(i=1;i<4;i++) vconff[i] = Qtconclean[i]/Bsq;

  // v_\alpha cleaned
  lower_g(vconff,ptrgeom->gcov,vcovff) ;

  // v^2 cleaned
  vsq=0.0;
  for(i=1;i<4;i++) vsq+=vconff[i]*vcovff[i]; // vcon[TT]=0

  //  dualfprintf(fail_file,"vsq=%21.15g Bsq=%21.15g alphasq=%21.15g Qtsq=%21.15g\n",vsq,Bsq,alphasq,Qtsq);

  // now limit the 3-velocity to be physical (no change in direction, so QtdotB=0 will be preserved)
  if( (vsq<0.0) || (vsq>1.0/alphasq) ){
    // then need to limit Lorentz factor
    // dissipation is expected to occur without affecting force-free equations

    // define cleaned Qtsq
    lower_g(Qtconclean,ptrgeom->gcov,Qtcovclean) ;
    Qtsqclean=0.0;
    for(i=1;i<4;i++) Qtsqclean+=Qtconclean[i]*Qtcovclean[i];
 
    // define new vconff (overwrites above) that is cleaned AND has dissipation
    // v^\alpha = P^\alpha/(B^2+\sqrt{\alpha^2 P^2})
    vconff[TT]=0.0;
    //    for(i=1;i<4;i++) vconff[i] = Qtconclean[i]/(Bsq+sqrt(alphasq*Qtsqclean)) ;
    //realBsq=sqrt(Bsq+sqrt(alphasq*Qtsqclean));
    // Note that E^2 = P^2/B^2, and then use same limiting procedure as in phys.ffde.c
    realBsq=sqrt(Qtsqclean*alphasq);
    for(i=1;i<4;i++) vconff[i] = Qtconclean[i]/realBsq ;

    // v_\alpha cleaned and w/ dissipation
    lower_g(vconff,ptrgeom->gcov,vcovff) ;

    // v^2 cleaned and w/ dissipation
    vsq=0.0;
    for(i=1;i<4;i++) vsq+=vconff[i]*vcovff[i]; // vcon[TT]=0

  }
 
  // notice that vsq<1 enforced now, so should be ok
  if( (vsq<-VSQ_TOO_BIG) || (vsq>1.0/alphasq + 100.0*NUMEPSILON) ){
    dualfprintf(fail_file,"Failed to limit vcon\n");
    if(vsq<-VSQ_TOO_BIG) vsq=-VSQ_TOO_BIG;
    else return(1);
  }

  // gamma cleaned and w/ dissipation if necessary
  // clearly limits accuracy of gamma -- unavoidable when using B^\alpha instead of b^\alpha as fundamental variable
  gamma=1.0/sqrt(1.0-vsq);

  // set relative 4-velocity
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma*vconff[i];
  
  // set field (not really used)
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  //  PLOOP(i) dualfprintf(fail_file,"vsq=%21.15g prim[%d]=%21.15g\n",vsq,i,prim[i]);

  return(0);

}


// returns v^2>=1
// does NOT depend on EOS
// W^3 (W+2B^2) - S^2 (2W+B^2) \le W^2 (P^2-B^4)

// GODMARK: vsqgtr1 based upon energy equation information?  What about P^2 equation inversion?
static int vsqgtr1(FTYPE W,FTYPE Bsq,FTYPE QdotBsq, FTYPE Qtsq)
{

  return(W*W*W * ( W + 2.*Bsq ) - QdotBsq*(2.*W + Bsq)  <= W*W*(Qtsq-Bsq*Bsq) );

}






//////////////////////////////////////////////
//
// GENERAL FUNCTIONS USED BY INVERSION
//
//////////////////////////////////////////////


// returns \tilde{u}^2 = \gamma^2 v^2
static FTYPE utsq_calc(FTYPE W)
{
  FTYPE utsqtop,utsqbottom,utsq;
  FTYPE S,SoW;

  S = QdotB;
  SoW = S/W;

  utsqtop = Qtsq + (Bsq+2.0*W)*SoW*SoW;

  utsqbottom = (Bsq*Bsq-Qtsq) + 2.0*Bsq*W + W*W - SoW*SoW*(Bsq+2.0*W);

  utsq = utsqtop/utsqbottom;

  //  if(icurr>90 && CRAZYNEWCHECK){
  if(CRAZYNEWCHECK&&0){
    dualfprintf(fail_file,"utsqtop=%21.15g utsqbottom=%21.15g Qtsq=%21.15g\n",utsqtop,utsqbottom,Qtsq);
  }

  return(utsq);

}

// returns \gamma^2
static FTYPE gammasq_calc_W(FTYPE W)
{
  FTYPE gammatop,gammabottom,gammasq;
  FTYPE S,SoW;

  S = QdotB;
  SoW = S/W;

  gammatop = Bsq+W;

  gammabottom = W*W+2.0*Bsq*W + (Bsq*Bsq-Qtsq)-2.0*SoW*SoW*W - Bsq*SoW*SoW;

  gammasq = gammatop*gammatop/gammabottom;

  return(gammasq);


}

// returns \gamma
static FTYPE gamma_calc_W(FTYPE W)
{
  static FTYPE gammasq_calc_W(FTYPE W);
  FTYPE gammasq,gamma;

  gammasq = gammasq_calc_W(W);

  // should I check sign of gammasq?  GODMARK

  gamma=sqrt(gammasq);

  return(gamma);

}




/* evaluate v^2 (spatial, normalized velocity) from W = \gamma^2 w */
// does NOT depend on EOS
static FTYPE vsq_calc(FTYPE W)
{
	FTYPE Wsq,Xsq,Ssq;
	FTYPE SoW;
	
	Wsq = W*W ;
	Xsq = (Bsq + W) * (Bsq + W);
	//	Ssq = QdotBsq / Bsq;

	SoW = QdotB/W;
	
	//return(  Ssq * ( 1./Wsq - 1./Xsq )  +  Qtsq / Xsq  ); 
	//	return(  ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq) );

	return(  (SoW*SoW*(2.0*W+Bsq)  +  Qtsq)/Xsq );

}

/* evaluate dv^2/dW */
// does NOT depend on EOS
static FTYPE dvsq_dW(FTYPE W)
{
	FTYPE W3,X3,Ssq,Wsq,X;
	
	X = Bsq + W;
	Wsq = W*W;
	W3 = Wsq*W ;
	X3 = X*X*X;

	//	if(fabs(Bsq)==0.0) Ssq=0.0;
	//	else Ssq = QdotBsq / Bsq;

	//return( -2.*( Ssq * ( 1./W3 - 1./X3 )  +  Qtsq / X3 ) ); 
	//	return( -2.*( W3*Qtsq + QdotBsq * ( 3*W*X + Bsq*Bsq ) ) / ( W3 * X3 )   );

	return( -2.*( Qtsq/X3  +  QdotBsq * (3.0*W*X + Bsq*Bsq) / ( W3 * X3 )  )  );

	//	return( -2.*( Qtsq/X3  +  QdotBsq/Bsq * (1.0/W3 - 1.0/X3)  )  ); // WRONG

}


// re-look at this GODMARK
// Solution to the cold GRHD equation P^2 = ...[Wp]
static int coldgrhd(FTYPE Qtsq, FTYPE D, FTYPE *Wp)
{
  FTYPE det;

  det = D*D+Qtsq;

  if(det<0.0){
    //    dualfprintf(fail_file,"coldgrhd: no solution: det=%21.15g\n",det);
    return(1);
  }
  else{
    if(D>0.0){
      *Wp = Qtsq/(D+sqrt(det));
    }
    else{
      *Wp = Qtsq/(D-sqrt(det));
    }
    return(0);
  }

  return(0);
}


// resultP = 0 = -P^2 + P^2 [ Wp ] for cold GRMHD

// old version with poles
static FTYPE Psq_Wp_old(FTYPE Wp)
{
  FTYPE result;
  FTYPE S,W,W2,DoW,SoW;
  FTYPE Y;

  W = Wp+D;
  Y = Bsq+2.0*W;
  W2 = W*W;
  S = QdotB;
  DoW = D/W;
  SoW = S/W;
 
#if(CRAZYDEBUG)
 if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"Qtsq=%21.15g Bsq=%21.15g D=%21.15g S=%21.15g\n",Qtsq,Bsq,D,S);
  }
#endif

  //  result = -Qtsq + (Bsq+2.0*D)*Bsq + 2.0*(Bsq+D)*Wp + Wp*Wp - (Bsq*DoW*DoW + SoW*SoW)*(Bsq+2.0*W);

  result = -Qtsq + Wp*(D+W)*(1.0 + Bsq*Y/W2) - Y*SoW*SoW;

  return(result);

}


// old version with poles
// result = d(resultP)/dWp = + dP^2[Wp]/dWp
static FTYPE dPsqdWp_Wp_old(FTYPE Wp)
{
  FTYPE result;
  FTYPE S,W,DoW,SoW;

  W = Wp+D;
  S = QdotB;
  DoW = D/W;
  SoW = S/W;

  result = 2.0*(Bsq+W)*(Bsq/W*DoW*DoW + SoW*SoW/W + 1.0);

  return(result);

}

// new version without poles
// returns actually W^2 P^2
static FTYPE Psq_Wp(FTYPE Wp)
{
  FTYPE result;
  FTYPE S,W,W2,DoW,SoW;
  FTYPE Y;

  W = Wp+D;
  Y = Bsq+2.0*W;
  W2 = W*W;
  S = QdotB;

#if(CRAZYDEBUG)
  //  if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
  if(icurr==0 && jcurr==31 && nstep==51 && steppart==2){
    dualfprintf(fail_file,"Qtsq=%21.15g Bsq=%21.15g D=%21.15g S=%21.15g Wp=%21.15g\n",Qtsq,Bsq,D,S,Wp);
  }
#endif

  result = -Qtsq*W2 + Wp*(D+W)*(W2 + Bsq*Y) - Y*S*S;

  return(result);

}


// new version without poles
// result = d(W^2 P^2)/dWp = 
static FTYPE dPsqdWp_Wp(FTYPE Wp)
{
  FTYPE result;
  FTYPE S,W,DoW,SoW;

  W = Wp+D;
  S = QdotB;

  result = 2.0*(D*( (Bsq+D)*(Bsq+D) - Qtsq) - (S*S) + ((Bsq + D)*(Bsq + 5.0*D) - Qtsq)*Wp + 3.0*(Bsq + 2.0*D)*(Wp*Wp) + 2.0*(Wp*Wp*Wp));

  return(result);

}










// 0 = -E' + E'[Wp] for hot GRMHD (can be used for cold GRMHD too)
// while this has poles near bad roots, the E\propto Wp for ultra relativistic case, so easy to find root
static FTYPE Eprime_Wp(FTYPE Wp)
{
  FTYPE result;
  FTYPE X,X2,W;
  FTYPE utsq;
  static FTYPE utsq_calc(FTYPE W);
  static FTYPE pressure_Wp_utsq(FTYPE Wp, FTYPE D, FTYPE utsq);
  FTYPE p_tmp;


  W = Wp+D;
  X = Bsq + W + SMALL;
  X2 = X*X;
 
  utsq=utsq_calc(W); // check for precision problems
  if(utsq<-0.1) utsq=0.0; // really if this happens probably not going to find solution

  p_tmp=pressure_Wp_utsq(Wp,D,utsq);
 
  result = Qdotnp + Wp - p_tmp + 0.5*Bsq + 0.5*(Bsq*Qtsq-QdotBsq)/X2;

  return(result);

}


// dE'/dWp [ Wp ]
// It is ok that there exists a term (1-dpdW) since even with ideal gas EOS in non-rel limit (1-dpdWp) still order unity
// as above, while poles near bad roots, linear dependence makes easier to find roots.
static FTYPE dEprimedWp(FTYPE Wp)
{
  FTYPE result;
  FTYPE X,X3,W;
  FTYPE utsq;
  static FTYPE utsq_calc(FTYPE W);
  static FTYPE dpdWp_calc_utsq(FTYPE Wp, FTYPE D, FTYPE utsq);
  static FTYPE dpdvsq_calc2_Wp(FTYPE Wp, FTYPE D, FTYPE utsq);
  FTYPE dvsq,dp1,dp2,dpdWp;


  W = Wp+D;
  X = Bsq + W + SMALL;
  X3 = X*X*X;

  utsq=utsq_calc(W); // check for precision problems
  if(utsq<-0.1) utsq=0.0; // really if this happens probably not going to find solution

  dvsq = dvsq_dW( W ); // no precision problem
  dp1 = dpdWp_calc_utsq(Wp,D,utsq);
  dp2 = dpdvsq_calc2_Wp( Wp, D, utsq );
  dpdWp = dp1  + dp2*dvsq; // dp/dW = dp/dWp

  result = (1.0 - dpdWp) - (Bsq*Qtsq - QdotBsq)/X3;

  return(result);

}


// this version has no poles
// assumes X2!=0, in which case previous old version gives inf anyways
// 0 = (B^2+W)^2 (-E' + E'[Wp]) for hot GRMHD (can be used for cold GRMHD too)
static FTYPE Eprime_Wp_new1(FTYPE Wp)
{
  FTYPE result;
  FTYPE X,X2,W;
  FTYPE utsq;
  static FTYPE utsq_calc(FTYPE W);
  static FTYPE pressure_Wp_utsq(FTYPE Wp, FTYPE D, FTYPE utsq);
  FTYPE p_tmp;


  W = Wp+D;
  X = Bsq + W + SMALL;
  X2 = X*X;
 
  utsq=utsq_calc(W); // check for precision problems
  if(utsq<-0.1) utsq=0.0;

  p_tmp=pressure_Wp_utsq(Wp,D,utsq);
 
  result = X2*(Qdotnp + Wp - p_tmp + 0.5*Bsq) + 0.5*(Bsq*Qtsq-QdotBsq);

  //  if(icurr>90 && CRAZYNEWCHECK){
  if(CRAZYNEWCHECK&&0){
    dualfprintf(fail_file,"Eprime: Wp=%21.15g W=%21.15g X2=%21.15g p_tmp=%21.15g utsq=%21.15g D=%21.15g result=%21.15g\n",Wp,W,X2,p_tmp,utsq,D,result);
  }

  return(result);

}


// this version has no poles
// d((B^2+W)^2 E')/dWp [ Wp ]
// It is ok that there exists a term (1-dpdW) since even with ideal gas EOS in non-rel limit (1-dpdWp) still order unity
static FTYPE dEprimedWp_new1(FTYPE Wp)
{
  FTYPE result;
  FTYPE X,X3,W;
  FTYPE utsq;
  static FTYPE utsq_calc(FTYPE W);
  static FTYPE dpdWp_calc_utsq(FTYPE Wp, FTYPE D, FTYPE utsq);
  static FTYPE dpdvsq_calc2_Wp(FTYPE Wp, FTYPE D, FTYPE utsq);
  FTYPE dvsq,dp1,dp2,dpdWp;
  static FTYPE pressure_Wp_utsq(FTYPE Wp, FTYPE D, FTYPE utsq);
  FTYPE p_tmp;


  W = Wp+D;
  X = Bsq + W + SMALL;
  X3 = X*X*X;

  utsq=utsq_calc(W); // check for precision problems
  if(utsq<-0.1) utsq=0.0;

  dvsq = dvsq_dW( W ); // no precision problem
  dp1 = dpdWp_calc_utsq(Wp,D,utsq);
  dp2 = dpdvsq_calc2_Wp( Wp, D, utsq );
  dpdWp = dp1  + dp2*dvsq; // dp/dW = dp/dWp

  // just copied from Eprime_Wp() to avoid recomputing utsq_calc()
  p_tmp=pressure_Wp_utsq(Wp,D,utsq); 

  result = X*(2.0*(Qdotnp + Wp - p_tmp + 0.5*Bsq) + (1.0 - dpdWp)*X );

  //  if(icurr>90 && CRAZYNEWCHECK){
  if(CRAZYNEWCHECK&&0){
    dualfprintf(fail_file,"dEprime: Wp=%21.15g W=%21.15g X3=%21.15g p_tmp=%21.15g utsq=%21.15g dp1=%21.15g dp2=%21.15g dpdWp=%21.15g result=%21.15g\n",Wp,W,X3,p_tmp,utsq,dp1,dp2,dpdWp,result);
  }

  return(result);

}




/********************************************************************

  x1_of_x0(): Only Used For 2D Inversion
           
    -- calculates x1 from x0 depending on value of vartype;

                          x0   x1
             vartype = 2 : Wp, gamma
             vartype = 3 : Wp, vsq
             vartype = 4 : Wp, utsq
    
    -- asumes x0 is already physical

*********************************************************************/

static FTYPE x1_of_x0(FTYPE x0 ) 
{
  FTYPE x1,vsq;
  //  FTYPE dv = 1.e-15;
  FTYPE W;
  static int validate_vsq(FTYPE vsqold,FTYPE *vsqnew);



  W = x0+D;
  
#if(CHANGEDTOOLDER)
  vsq = fabs(vsq_calc(W)) ; // guaranteed to be positive (>=0)
#else
  vsq = vsq_calc(W) ; // ok if not positive
#endif

  validate_vsq(vsq,&vsq);

  //  return( ( vsq > 1. ) ? VSQ_TOO_BIG : vsq   ); 

  return(vsq);

}

/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0,1] have physical values, based upon their definitions:

// added wglobal to properly normalize -- global variable set before find_root_2d()

// avoided truncating bottom value of x[0] and x[1] since then can't converge properly to slightly unphysical values or temporarily reach out into unphysical space and come back in.

// avoid truncating x[1] for v^2<0 since all equations are ok for v^2<0 and can help convergence to nearly physical or even physical solutions

// must avoid x[0] too large and v^2>1, however, since that produces nan's.
    
*********************************************************************/


// Wpold is previous Newton iteration's value
// *Wpnew is updated Newton iterations value
static void validate_Wp(FTYPE Wpold, FTYPE *Wpnew)
{
  static int coldgrhd(FTYPE Qtsq, FTYPE D, FTYPE *Wp);
  FTYPE Wp;
  FTYPE dv = NUMEPSILON;

  /* Always take the absolute value of Wpnew[0] and check to see if it's too big:  */ 
  //  Wpnew[0] = fabs(Wpold);
  if(Wpnew[0]<-D && D>0.0){
    //    Wpnew[0] = 10.0*Qtsq/(2.0*D+Wp0);
    //    Wpnew[0]=10.0*D+dv;

    // reset above HD solution
    coldgrhd(Qtsq, D, &Wp);
    Wpnew[0] = Wp*10.0;

#if(CRAZYDEBUG)
    dualfprintf(fail_file,"Wpnew[0] changed from %21.15g to %21.15g (D=%21.15g Qtsq=%21.15g)\n",Wpold,Wpnew[0],D,Qtsq);
#endif    

  }

  // revert to previous Wp if "new" Wp very bad
  // can't go back to previous Wp otherwise won't go anywhere
  // Wpnew[0]/global = gamma^2
  Wpnew[0] = (fabs(Wpnew[0]/wglobal) > GAMMASQ_TOO_BIG) ?  0.5*(Wpold+Wpnew[0]) : Wpnew[0];


}


// retvsq=0 : no problem
// retvsq=1 : >=1.0
// retvsq=2 : < too negative
// retvsq=3 : >VSQ_TOO_BIG (forced failure even if converged)
static int validate_vsq(FTYPE vsqold,FTYPE *vsqnew)
{
  int retval;


  retval=0; // assume no problem

  if(vsqnew[0]>1.0){
    retval=1;
    vsqnew[0]=VSQ_TOO_BIG;
    return(retval);
  }


  if(vsqnew[0]>VSQ_TOO_BIG){
    retval=3;
    vsqnew[0]=VSQ_TOO_BIG;
    return(retval);
  }

#if(VSQNEGCHECK==1)
  if( (vsqnew[0] < 0.) && (fabs(vsqnew[0]) < MAXNEGVSQ) ) { 
    vsqnew[0] = 0.0;
    retval=2;
    return(retval);
  }
#elif(VSQNEGCHECK==2)
  if(vsqnew[0] <=-VSQ_TOO_BIG) { 
    vsqnew[0] = -VSQ_TOO_BIG;
    retval=2;
    return(retval);
  }
#elif(VSQNEGCHECK==3)
  if(vsqnew[0]<0.0) {
    vsqnew[0]=0.0;
    retval=2;
    return(retval);
  }
#endif


  return(retval);

}


static void validate_x_2d(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM] ) 
{
  static void validate_Wp(FTYPE Wpold, FTYPE *Wpnew);
  static int validate_vsq(FTYPE vsqold,FTYPE *vsqnew);


  validate_Wp(x0[0], &x[0]);

  validate_vsq(x0[1], &x[1]);


}


// 1D methods using W'
static void validate_x_1d_old(FTYPE x[1], FTYPE x0[1])
{
  FTYPE dv = NUMEPSILON;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */
  //  x[0] = fabs(x[0]);
  if(x[0]<-D) x[0]=D+dv;
  // x[0]/global = gamma^2
  x[0] = (x[0]/wglobal > GAMMASQ_TOO_BIG) ?  x0[0] : x[0];

}



// 1D methods using W'
static void validate_x_1d(FTYPE x[1], FTYPE x0[1])
{
  static void validate_Wp(FTYPE Wpold, FTYPE *Wpnew);

  validate_Wp(x0[0], &x[0]);

}


static void pick_validate_x(int eomtype)
{

  if(
     (WHICHHOTINVERTER==1 && eomtype==EOMGRMHD) ||
     (WHICHHOTINVERTER==3 && eomtype==EOMGRMHD) ||
     ( (WHICHCOLDINVERTER==0 || WHICHCOLDINVERTER==1) && eomtype==EOMCOLDGRMHD) ||
     (WHICHCOLDINVERTER==2 && eomtype==EOMCOLDGRMHD)
     ){
    ptr_validate_x = &validate_x_1d;
  }
  else if(WHICHHOTINVERTER==2 && eomtype==EOMGRMHD){
    ptr_validate_x = &validate_x_2d;
  }
  else {
    dualfprintf(fail_file,"pick_validate_x: Shouldn't reach here\n");
    myexit(999);
  }


}

static void validate_x(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM] ) 
{
  (*ptr_validate_x)(x,x0);

}



//////////////////////////////////////////////////////
//
//
//////////////  STUFF FOR 1D METHOD
//
//
//////////////////////////////////////////////////////




/**************************************************** 
*****************************************************/

static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE vsq,vsqfake,gamma,gamma_sq,rho0,w,W,Wp,drdW,dpdW,Wsq,p_tmp, dvsq,dp1,dp2,a_term,ap_term ;
  FTYPE X,X2,X3;
  static FTYPE pressure_Wp_vsq(FTYPE Wp, FTYPE D, FTYPE vsq) ;
  static FTYPE dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq);
  static FTYPE dpdvsq_calc_Wp(FTYPE Wp, FTYPE D, FTYPE vsq);
  static FTYPE dvsq_dW(FTYPE W);
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;
  static int validate_vsq(FTYPE vsqold,FTYPE *vsqnew);
  

  Wp = x[0];
  W = Wp+D;
  Wsq = W*W;
  X = Bsq+W;
  X2 = X*X;
  X3 = X2*X;

  //dualfprintf(fail_file,"call utsq in residual\n") ;
  vsq = vsq_calc(W) ;
  //dualfprintf(fail_file,"done w/ utsq in residual\n") ;

  // GODMARK DANGEROUS 
  //  vsqfake = ( vsq < -1.e-10 ) ?  dv : fabs(vsq) ;
  //  vsqfake = ( vsq > 1. ) ?  (1.-dv) : vsq ;
 
  validate_vsq(vsq,&vsqfake);

#if(0)
  // no longer use below methods 
#if(OPTIMIZED)
  vsqfake = ( vsq < -MAXNEGVSQ ) ?  0.0 : fabs(vsq) ;
  vsqfake = ( vsq > 1. ) ?  VSQ_TOO_BIG : vsq ;
#elif(0)
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsqfake=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<1.0+MAXNEGVSQ) vsqfake=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif
#endif

  //rho0 = D * sqrt(1. - vsq) ;
  //dualfprintf(fail_file,"rho0,D,gamma: %21.15g %21.15g %21.15g\n",rho0,D,gamma) ;
  //  w = W * (1. - vsq);

  p_tmp = pressure_Wp_vsq(Wp,D,vsqfake); // needs to be correct



  // Same as previous residual, however, jacobian is calculated using full differentiation w.r.t. W  
  //
  // dp/dW' = dp/dW + dP/dv^2 dv^2/dW


  dvsq = dvsq_dW( W ); // no precision problem
  dp1 = dpdWp_calc_vsq( Wp, D, vsq ); // vsq can be unphysical
  dp2 = dpdvsq_calc_Wp( Wp, D, vsq ); // vsq must be physical
  dpdW = dp1  + dp2*dvsq; // dp/dW = dp/dWp


#if(0) // Scott

  resid[0] = 
    + Wp 
    + 0.5 * Bsq * ( 1. + vsq ) // any vsq
    - 0.5*QdotBsq/Wsq
    + Qdotnp
    - p_tmp;

  jac[0][0] = drdW = 1. - dpdW + QdotBsq/(Wsq*W) + 0.5*Bsq*dvsq;

#elif(0)
//  resid[0] = 
//    + W*Wsq
//    + 0.5 * Wsq * Bsq * ( 1. + vsq )
//    - 0.5*QdotBsq
//    - Qdotn*Wsq
//    - p_tmp*Wsq;
//
//  jac[0][0] = drdW = W * ( 3.*W  +  Bsq * (1. + vsq + 0.5*W*dvsq)  - 2.*Qdotn - 2.*p_tmp - W*dpdW );

#elif(1)  // JON:

  resid[0] = Qdotnp + Wp - p_tmp + 0.5*Bsq + (Bsq*Qtsq - QdotBsq)/X2;

  jac[0][0] = drdW = 1. - dpdW + (Bsq*Qtsq - QdotBsq)/X3*(-2.0);

#endif


  dx[0] = -resid[0]/drdW;

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

#if(!OPTIMIZED)
  if( ltrace ) {
    a_term = 0.5*Bsq*(1.+vsq) + Qdotnp - p_tmp;
    ap_term = dvsq * ( 0.5*Bsq - dp2 ) - dp1;

    dualfprintf(fail_file,"func_1d_orig(): x = %21.15g, dx = %21.15g, resid = %21.15g, drdW = %21.15g \n", x[0], dx[0], resid[0], drdW);
    dualfprintf(fail_file,"func_1d_orig(): dvsq = %21.15g,  dp = %21.15g , dp1 = %21.15g , dp2 = %21.15g \n", dvsq, dpdW,dp1,dp2);
    dualfprintf(fail_file,"func_1d_orig(): W = %21.15g , vsq = %21.15g , a = %21.15g,  a_prime = %21.15g  \n", x[0], vsq, a_term, ap_term);
    
  }
#endif

}





/**************************************************** 
  Routine for line searching for the 1D method . 
*****************************************************/

static FTYPE res_sq_1d_orig(FTYPE x[])
{
  FTYPE vsq,W,Wp,Wsq,p_tmp,resid[1];
  FTYPE vsqfake;
  FTYPE X,X2;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;
  static int validate_vsq(FTYPE vsqold,FTYPE *vsqnew);


  Wp = x[0];
  W = Wp+D;
  Wsq = W*W;
  X = Bsq+W;
  X2 = X*X;

  vsq = vsq_calc(W) ;


  validate_vsq(vsq,&vsqfake);

#if(0)
  // GODMARK DANGEROUS 
#if(OPTIMIZED)
  // scn's version
  //  vsqfake = ( vsq < -dv ) ?  dv : fabs(vsq) ;
  //vsqfake = ( vsq > 1. ) ?  (1.-dv) : vsq ;
  // jon's version
  vsqfake = ( vsq < 0.0 ) ?  0.0 : vsq ;
  vsqfake = ( vsq > VSQ_TOO_BIG ) ? VSQ_TOO_BIG  : vsq ; // GODMARK: must let iterate if vsq>1.  Can't just force vsq or else Newton's method won't work.
#else
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsqfake=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<1.0+MAXNEGVSQ) vsqfake=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif
#endif

  p_tmp = pressure_Wp_vsq(Wp,D,vsqfake); // vsq<=1 for pressure



#if(0) // Scott
  resid[0] = 
    + Wp 
    + 0.5 * Bsq * ( 1. + vsq ) // residual must use unphysical vsq
    - 0.5*QdotBsq/Wsq
    + Qdotnp
    - p_tmp;

#elif(1)  // JON:

  resid[0] = Qdotnp + Wp - p_tmp + 0.5*Bsq + (Bsq*Qtsq - QdotBsq)/X2;

#endif

  return(  0.5*resid[0]*resid[0] );


}






/////////////////////////////////
//
// SCN versions of 1D method
//
/////////////////////////////////

static void func_1d_orig_scn(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE vsq,vsqfake,gamma,gamma_sq,rho0,w,W,drdW,dpdW,Wsq,p_tmp, dvsq,dp1,dp2,a_term,ap_term ;
  static FTYPE pressure_W_vsq(FTYPE W, FTYPE D, FTYPE vsq) ;
  static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE D, FTYPE vsq);
  static FTYPE dpdvsq_calc(FTYPE W, FTYPE D, FTYPE vsq);
  static FTYPE dvsq_dW(FTYPE W);
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;


  W = x[0];
  Wsq = W*W;

  //dualfprintf(fail_file,"call utsq in residual\n") ;
  vsq = vsq_calc(W) ;
  //dualfprintf(fail_file,"done w/ utsq in residual\n") ;

  // GODMARK DANGEROUS 
  //  vsq = ( vsq < -1.e-10 ) ?  dv : fabs(vsq) ;
  //  vsq = ( vsq > 1. ) ?  (1.-dv) : vsq ;
  
#if(OPTIMIZED)
  vsqfake = ( vsq < -MAXNEGVSQ ) ?  0.0 : fabs(vsq) ;
  vsqfake = ( vsq > 1. ) ?  VSQ_TOO_BIG : vsq ;
#else
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsqfake=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<1.0+MAXNEGVSQ) vsqfake=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif

  //rho0 = D * sqrt(1. - vsq) ;
  //dualfprintf(fail_file,"rho0,D,gamma: %21.15g %21.15g %21.15g\n",rho0,D,gamma) ;
  //  w = W * (1. - vsq);

  p_tmp = pressure_W_vsq_scn(W,D,vsqfake);



  // Same as previous residual, however, jacobian is calculated using full differentiation w.r.t. W  

  dvsq = dvsq_dW( W );
  dp1 = dpdW_calc_vsq_scn( W, D, vsq ); // use real vsq
  dp2 = dpdvsq_calc_scn( W, D, vsqfake );
  dpdW = dp1  + dp2*dvsq;

  resid[0] = 
    + W 
    + 0.5 * Bsq * ( 1. + vsq ) // use real vsq
    - 0.5*QdotBsq/Wsq
    + Qdotn
    - p_tmp;

  jac[0][0] = drdW = 1. - dpdW + QdotBsq/(Wsq*W) + 0.5*Bsq*dvsq;

//  resid[0] = 
//    + W*Wsq
//    + 0.5 * Wsq * Bsq * ( 1. + vsq )
//    - 0.5*QdotBsq
//    - Qdotn*Wsq
//    - p_tmp*Wsq;
//
//  jac[0][0] = drdW = W * ( 3.*W  +  Bsq * (1. + vsq + 0.5*W*dvsq)  - 2.*Qdotn - 2.*p_tmp - W*dpdW );

  dx[0] = -resid[0]/drdW;

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

#if(!OPTIMIZED)
  if( ltrace ) {
    a_term = 0.5*Bsq*(1.+vsq) + Qdotn - p_tmp;
    ap_term = dvsq * ( 0.5*Bsq - dp2 ) - dp1;

    dualfprintf(fail_file,"func_1d_orig(): x = %21.15g, dx = %21.15g, resid = %21.15g, drdW = %21.15g \n", x[0], dx[0], resid[0], drdW);
    dualfprintf(fail_file,"func_1d_orig(): dvsq = %21.15g,  dp = %21.15g , dp1 = %21.15g , dp2 = %21.15g \n", dvsq, dpdW,dp1,dp2);
    dualfprintf(fail_file,"func_1d_orig(): W = %21.15g , vsq = %21.15g , a = %21.15g,  a_prime = %21.15g  \n", x[0], vsq, a_term, ap_term);
    
  }
#endif

}

/**************************************************** 
  Routine for line searching for the 1D method . 
*****************************************************/

static FTYPE res_sq_1d_orig_scn(FTYPE x[])
{
  FTYPE vsq,vsqfake,W,Wsq,p_tmp,resid[1];
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;


  W = x[0];
  Wsq = W*W;

  vsq = vsq_calc(W) ;
  
  // GODMARK DANGEROUS 
#if(OPTIMIZED)
  // scn's version
  //  vsq = ( vsq < -dv ) ?  dv : fabs(vsq) ;
  //vsq = ( vsq > 1. ) ?  (1.-dv) : vsq ;
  // jon's version
  vsqfake = ( vsq < 0.0 ) ?  0.0 : vsq ;
  vsqfake = ( vsq > VSQ_TOO_BIG ) ? VSQ_TOO_BIG  : vsq ;
#else
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsqfake=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<1.0+MAXNEGVSQ) vsqfake=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif

  p_tmp = pressure_W_vsq_scn(W,D,vsqfake);


  resid[0] = 
    + W 
    + 0.5 * Bsq * ( 1. + vsq )
    - 0.5*QdotBsq/Wsq
    + Qdotn
    - p_tmp;

  return(  0.5*resid[0]*resid[0] );


}








static void func_Eprime(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE drdW;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;
  static FTYPE Eprime_Wp(FTYPE Wp);
  FTYPE Wp;


  Wp = x[0];

  resid[0] = Eprime_Wp(Wp);
  jac[0][0] = drdW = dEprimedWp(Wp);


  dx[0] = -resid[0]/drdW;
  //  if((fabs(dx[0])/(Wp+D)>0.1){ // don't trust such a big jump
  // dx[0] = sign(dx[0])*0.1;
  // }
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);


}


static FTYPE res_sq_Eprime(FTYPE x[])
{
  FTYPE Wp;
  FTYPE resid;
  static FTYPE Eprime_Wp(FTYPE Wp);

  Wp = x[0];
  resid = Eprime_Wp(Wp);

  return(  0.5*resid*resid );


}




static void func_Psq(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE drdW,Wp;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;
  static FTYPE Psq_Wp(FTYPE Wp);
  static FTYPE dPsqdWp_Wp(FTYPE Wp);


  Wp = x[0];


  resid[0] = Psq_Wp(Wp);
  jac[0][0] = drdW = dPsqdWp_Wp(Wp);


  dx[0] = -resid[0]/drdW;
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);


}


static FTYPE res_sq_Psq(FTYPE x[])
{
  FTYPE Wp;
  FTYPE resid;
  static FTYPE Psq_Wp(FTYPE Wp);

  Wp = x[0];
  resid = Psq_Wp(Wp);

  return(  0.5*resid*resid );


}




/////////////////////////////////////////////////////
//
//////////////// STUFF FOR 2D METHOD
//
/////////////////////////////////////////////////////




/*********************************************************
returns residual  of vsq-1
**********************************************************/

static void func_vsq(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{

  
  FTYPE  Wp,W, vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  const int ltrace = 0;


  Wp = x[0];
  W = Wp+D;
  vsq = x[1];
  
  Wsq = W*W;

  p_tmp  = pressure_Wp_vsq( Wp, D, vsq );
  dPdW   = dpdWp_calc_vsq( Wp, D, vsq );
  dPdvsq = dpdvsq_calc_Wp( Wp, D, vsq );


#if(CRAZYDEBUG)
  if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"p_tmp=%21.15g dPdW=%21.15g dPdvsq=%21.15g\n",p_tmp,dPdW,dPdvsq);
  }
#endif


  /* These expressions were calculated using Mathematica */
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  dx[0] = (-Bsq/2. + dPdvsq)*(Qtsq - vsq*((Bsq+W)*(Bsq+W)) + 
			      (QdotBsq*(Bsq + 2*W))/Wsq) + 
    ((Bsq+W)*(Bsq+W))*(-Qdotnp - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) -Wp + p_tmp);

  dx[1] = -((-1 + dPdW - QdotBsq/(Wsq*W))*
	    (Qtsq - vsq*((Bsq+W)*(Bsq+W)) + (QdotBsq*(Bsq + 2*W))/Wsq)) - 
    2*(vsq + QdotBsq/(Wsq*W))*(Bsq + W)*
    (-Qdotnp - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - Wp + p_tmp);

  detJ = (Bsq + W)*((-1 + dPdW - QdotBsq/(Wsq*W))*(Bsq + W) + 
		     ((Bsq - 2*dPdvsq)*(QdotBsq + vsq*(Wsq*W)))/(Wsq*W));
  
  dx[0] /= -(detJ) ;
  dx[1] /= -(detJ) ;

  jac[0][0] = -2*(vsq + QdotBsq/(Wsq*W))*(Bsq + W);
  jac[0][1] = -(Bsq + W)*(Bsq + W);
  jac[1][0] = -1 - QdotBsq/(Wsq*W) + dPdW;
  jac[1][1] = -Bsq/2. + dPdvsq;
  resid[0]  = Qtsq - vsq*(Bsq + W)*(Bsq + W) + (QdotBsq*(Bsq + 2*W))/Wsq;
  resid[1]  = -Qdotnp - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - Wp + p_tmp;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"func_vsq(): x[0] = %21.15g , x[1] = %21.15g , dx[0] = %21.15g , dx[1] = %21.15g , f = %21.15g \n", 
	    x[0], x[1],dx[0], dx[1], *f);
    
  }
#endif

}


/**********************************************************/

static void func_vsq2(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{

  
  FTYPE  W, Wp,vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  FTYPE t11;
  FTYPE t16;
  FTYPE t18;
  FTYPE t2;
  FTYPE t21;
  FTYPE t23;
  FTYPE t24;
  FTYPE t25;
  FTYPE t3;
  FTYPE t35;
  FTYPE t36;
  FTYPE t4;
  FTYPE t40;
  FTYPE t9;

  const int ltrace = 0;


  Wp = x[0];
  W = Wp+D;
  vsq = x[1];
  
  Wsq = W*W;
  
  p_tmp  = pressure_Wp_vsq( Wp, D, vsq );
  dPdW   = dpdWp_calc_vsq( Wp, D, vsq );
  dPdvsq = dpdvsq_calc_Wp( Wp, D, vsq );

#if(CRAZYDEBUG)
    if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"W=%21.15g D=%21.15g vsq=%21.15g\n",W,D,vsq);
      dualfprintf(fail_file,"p_tmp=%21.15g dPdW=%21.15g dPdvsq=%21.15g\n",p_tmp,dPdW,dPdvsq);
    }
#endif

#if(0)
  dualfprintf(fail_file,"p_tmp=%21.15g dPdW=%21.15g dPdvsq=%21.15g\n",p_tmp,dPdW,dPdvsq);
#endif

  /* These expressions were calculated using Mathematica, but made into efficient code using Maple */
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  t2 = -0.5*Bsq+dPdvsq;
  t3 = Bsq+W;
  t4 = t3*t3;
  t9 = 1/Wsq;
  t11 = Qtsq-vsq*t4+QdotBsq*(Bsq+2.0*W)*t9;
  t16 = QdotBsq*t9;
  t18 = -Qdotnp-0.5*Bsq*(1.0+vsq)+0.5*t16-Wp+p_tmp;
  //t18 = -Qdotn-0.5*Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t21 = 1/t3;
  t23 = 1/W;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(Bsq-2.0*dPdvsq)*(QdotBsq+vsq*Wsq*W)*t9*t23;
  t36 = 1/t35;
  dx[0] = -(t2*t11+t4*t18)*t21*t36;

#if(0)
  dualfprintf(fail_file,"Qdotn=%21.15g Qdotnp=%21.15g Wp=%21.15g p_tmp=%21.15g\n",Qdotn,Qdotnp,Wp,p_tmp);
  dualfprintf(fail_file,"t2=%21.15g t3=%21.15g t4=%21.15g t9=%21.15g t11=%21.15g t16=%21.15g t18=%21.15g t21=%21.15g t23=%21.15g t24=%21.15g t25=%21.15g t35=%21.15g t36=%21.15g dx[0]=%21.15g\n",t2,t3,t4,t9,t11,t16,t18,t21,t23,t24,t25,t35,t36,dx[0]);
#endif

  t40 = (vsq+t24)*t3;
  dx[1] = -(-t25*t11-2.0*t40*t18)*t21*t36;
  detJ = t3*t35;
  jac[0][0] = -2.0*t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

#if(CRAZYDEBUG)
  if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"Qdotn=%21.15g QdotBsq=%21.15g Qtsq=%21.15g Bsq=%21.15g vsq=%21.15g W=%21.15g p_tmp=%21.15g -Qdotnp-Wp=%21.15g\n",Qdotn,QdotBsq,Qtsq,Bsq,vsq,W,p_tmp,-Qdotnp-Wp);
    
    dualfprintf(fail_file,"t2=%21.15g\nt3=%21.15g\nt4=%21.15g\nt9=%21.15g\nt11=%21.15g\nt16=%21.15g\nt18=%21.15g\nt21=%21.15g\nt23=%21.15g\nt24=%21.15g\nt25=%21.15g\nt35=%21.15g\nt36=%21.15g\ndx[0]=%21.15g\nt40=%21.15g\ndx[1]=%21.15g\ndetJ=%21.15g\njac[0][0]=%21.15g\njac[0][1]=%21.15g\njac[1][0]=%21.15g\njac[1][1]=%21.15g\nresid[0]=%21.15g\nresid[1]=%21.15g\n",t2,t3,t4,t9,t11,t16,t18,t21,t23,t24,t25,t35,t36,dx[0],t40,dx[1],detJ,jac[0][0],jac[0][1],jac[1][0],jac[1][1],resid[0],resid[1]);
    
    dualfprintf(fail_file,"*df=%21.15g *f=%21.15g Wp=%21.15g W=%21.15g\n",*df,*f,Wp,W);
  }
#endif


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"func_vsq2(): x[0] = %21.15g , x[1] = %21.15g , dx[0] = %21.15g , dx[1] = %21.15g , f = %21.15g \n", 
	    x[0], x[1],dx[0], dx[1], *f);
    
  }
#endif

}


/*********************************************************
**********************************************************/
static FTYPE res_sq_vsq(FTYPE x[])
{

  
  FTYPE  W, Wp, vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  FTYPE  resid[2];
  const int ltrace = 0;


  Wp = x[0];
  W = Wp+D;
  vsq = x[1];
  
  Wsq = W*W;
  
  p_tmp  = pressure_Wp_vsq( Wp, D, vsq );
  //  dPdW   = dpdW_calc_vsq( Wp, D, vsq );
  //  dPdvsq = dpdvsq_calc_Wp( Wp, D, vsq );

  /* These expressions were calculated using Mathematica */
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  resid[0]  = Qtsq - vsq*(Bsq + W)*(Bsq + W) + (QdotBsq*(Bsq + 2*W))/Wsq;
  resid[1]  = -Qdotnp - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - Wp + p_tmp;


  tmp3 = 0.5 * ( resid[0]*resid[0] + resid[1]*resid[1] ) ;
  
#if(!OPTIMIZED)
  if( ltrace  ) {
    dualfprintf(fail_file,"res_sq_vsq(): W,vsq,resid0,resid1,f =   %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  \n",
	    W,vsq,resid[0],resid[1],tmp3);
    dualfprintf(fail_file,"res_sq_vsq(): Qtsq, Bsq, QdotBsq, Qdotn =   %21.15g  %21.15g  %21.15g  %21.15g  \n",
	    Qtsq, Bsq, QdotBsq, Qdotn);
  }
#endif
    
  return( tmp3 );
}


/*********************************************************
**********************************************************/
static FTYPE res_sq_vsq2(FTYPE x[])
{

  FTYPE  W, Wp, vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  FTYPE  resid[2];
  FTYPE t3, t4 ;
  FTYPE t9 ;
  FTYPE t11;
  FTYPE t16;
  FTYPE t18;

  const int ltrace = 0;


  Wp = x[0];
  W = Wp+D;
  vsq = x[1];
  
  Wsq = W*W;
  
  p_tmp  = pressure_Wp_vsq( Wp, D, vsq );

  /* These expressions were calculated using Mathematica and made into efficient code using Maple*/
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  t3 = Bsq+W;
  t4 = t3*t3;//
  t9 = 1/Wsq;  //
  t11 = Qtsq-vsq*t4+QdotBsq*(Bsq+2.0*W)*t9; //
  t16 = QdotBsq*t9;//
  t18 = -Qdotnp-0.5*Bsq*(1.0+vsq)+0.5*t16-Wp+p_tmp;//
  resid[0] = t11;
  resid[1] = t18;

  tmp3 = 0.5 * ( resid[0]*resid[0] + resid[1]*resid[1] ) ;
  
#if(!OPTIMIZED)
  if( ltrace  ) {
    dualfprintf(fail_file,"res_sq_vsq(): W,vsq,resid0,resid1,f =   %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  \n",
	    W,vsq,resid[0],resid[1],tmp3);
    dualfprintf(fail_file,"res_sq_vsq(): Qtsq, Bsq, QdotBsq, Qdotn =   %21.15g  %21.15g  %21.15g  %21.15g  \n",
	    Qtsq, Bsq, QdotBsq, Qdotn);
  }
#endif

#if(CRAZYDEBUG)
  if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"res: W=%21.15g Wsq=%21.15g p_tmp=%21.15g t3=%21.15g t4=%21.15g t9=%21.15g t11=%21.15g t16=%21.15g t18=%21.15g resid[0]=%21.15g resid[1]=%21.15g tmp3=%21.15g\n",W,Wsq,p_tmp,t3,t4,t9,t11,t16,t18,resid[0],resid[1],tmp3);
  }
#endif

    
  return( tmp3 );
}










/**************************************************** 


GENERAL PROCEDURES BELOW (except validate_x() and must choose function and residual)


*****************************************************/


static int find_root_1D_gen(FTYPE x0, FTYPE *xnew)
{
  FTYPE x_1d[1];
  int retval;
  static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
				     void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
						    FTYPE *, FTYPE *, int), 
				     FTYPE (*res_func) (FTYPE []) );
  static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static FTYPE res_sq_1d_orig( FTYPE [] );


  x_1d[0] = x0;

  if( (retval=general_newton_raphson( x_1d, 1, USE_LINE_SEARCH, func_1d_orig, res_sq_1d_orig ) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}



static int find_root_1D_gen_scn(FTYPE x0, FTYPE *xnew)
{
  FTYPE x_1d[1];
  int retval;
  static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
				     void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
						    FTYPE *, FTYPE *, int), 
				     FTYPE (*res_func) (FTYPE []) );
  static void func_1d_orig_scn(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static FTYPE res_sq_1d_orig_scn( FTYPE [] );


  x_1d[0] = x0;

  if( (retval=general_newton_raphson( x_1d, 1, USE_LINE_SEARCH, func_1d_orig_scn, res_sq_1d_orig_scn ) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}



static int find_root_1D_gen_Eprime(FTYPE x0, FTYPE *xnew)
{
  FTYPE x_1d[1];
  int retval;
  static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
				     void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
						    FTYPE *, FTYPE *, int), 
				     FTYPE (*res_func) (FTYPE []) );
  static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static FTYPE res_sq_1d_orig( FTYPE [] );


  x_1d[0] = x0;

  if( (retval=general_newton_raphson( x_1d, 1, USE_LINE_SEARCH, func_Eprime, res_sq_Eprime ) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}

static int find_root_1D_gen_Psq(FTYPE x0, FTYPE *xnew)
{
  FTYPE x_1d[1];
  int retval;
  static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
				     void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
						    FTYPE *, FTYPE *, int), 
				     FTYPE (*res_func) (FTYPE []) );
  static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static FTYPE res_sq_1d_orig( FTYPE [] );


  x_1d[0] = x0;

  if( (retval=general_newton_raphson( x_1d, 1, USE_LINE_SEARCH, func_Psq, res_sq_Psq ) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}




// below not done yet!
static int find_root_3D_gen_Palpha(FTYPE x0, FTYPE *xnew)
{
  FTYPE x_3d[3];
  int retval;
  static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
				     void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
						    FTYPE *, FTYPE *, int), 
				     FTYPE (*res_func) (FTYPE []) );
  static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static FTYPE res_sq_1d_orig( FTYPE [] );


  x_3d[0] = x0;
  x_3d[1] = x0;
  x_3d[2] = x0;

  if( (retval=general_newton_raphson( x_3d, 3, USE_LINE_SEARCH, func_Psq, res_sq_Psq ) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_3d[0] = %21.15g x_3d[1] = %21.15g x_3d[2] = %21.15g \n", x_3d[0], x_3d[1], x_3d[2] );  
#endif
  }

  *xnew = x_3d[0];
  *xnew = x_3d[1];
  *xnew = x_3d[2];

  return(retval);
}





/********************************************************************

   find_root_2D_gen(): 
  
         -- performs a 2D Newton-Raphson method for solving the \tilde{Q}^2
             and Q.n equations for two variables in general;
              (usually W and something else like utsq, vsq, or gamma)

         residual vector       = { Qtsq eq. , Q.n eq. }
         indep. var. vector, x = { W, vsq/utsq/gamma/? }

       vartype = 2 : W, gamma
       vartype = 3 : W, vsq
       vartype = 4 : W, utsq
    
     -- Uses  general_newton_raphson();

*********************************************************************/

static int find_root_2D_gen(FTYPE x0, FTYPE *xnew) 
{
  
  FTYPE x[2], x_orig[2];
  int ntries = 0;
  int retval, n = 2;
  int it;
  static  void func_vsq(   FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static  void func_vsq2(   FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static  FTYPE  res_sq_vsq( FTYPE [] );
  static  FTYPE  res_sq_vsq2( FTYPE [] );
  static FTYPE x1_of_x0(FTYPE x0 ) ;
  static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
				     void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
						    FTYPE *, FTYPE *, int), 
				     FTYPE (*res_func) (FTYPE []) );

  const int ltrace = 0;

  /* Set presets: */
  x_orig[0] = x[0] = fabs(x0) ;
  x_orig[1] = x[1] = x1_of_x0( x0 ) ;

#if(CRAZYDEBUG)
  if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"x[0]=%21.15g W=%21.15g x[1]=%21.15g\n",x[0],x[0]+D,x[1]);
  }
#endif

#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file, "find_root_2D_gen():  x[0] = %21.15g , x[1] = %21.15g \n", x[0], x[1]);
    
  }
#endif

#if(0)
  for(it=0;it<n;it++) dualfprintf(fail_file,"x[%d]=%21.15g xnew[%d]=%21.15g\n",it,x[it],it,xnew[it]);
#endif

  retval = general_newton_raphson( x, n, USE_LINE_SEARCH,  func_vsq2, res_sq_vsq2 ) ;
  ntries++;

  while( (retval==1) && (ntries <= MAX_NEWT_RETRIES ) ) {
    //       x[0] = x_orig[0] * (1. + ( (1.*rand())/(1.*RAND_MAX) - 0.5 ) );  
    x[0] = x_orig[0];
    for( it=0; it<ntries; it++)  x[0] *= 10.0;
    x[1] = x1_of_x0( x[0] ) ;
    //       retval = general_newton_raphson( x, n, USE_LINE_SEARCH,  func_vsq, res_sq_vsq ) ;
    retval = general_newton_raphson( x, n, USE_LINE_SEARCH,  func_vsq2, res_sq_vsq2 ) ;

#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file, "find_root_2D_gen():  ntries, x[0,1] = %4i  %21.15g   %21.15g  \n", ntries, x[0], x[1]);
      
    }
#endif
    ntries++;
  }

#if( MAX_NEWT_RETRIES > 0 )
  if( (ntries > MAX_NEWT_RETRIES) &&  (retval==1)  ) {
    if( debugfail >= 2 ) { 
      dualfprintf(fail_file, "find_root_2D_gen():  Bad exit value from general_newton_raphson() !! \n");
      dualfprintf(fail_file, "find_root_2D_gen():  ntries = %d , x[0] = %21.15g , x[1] = %21.15g  \n", ntries, x[0], x[1]); 
    }
  }
#endif

  *xnew = x[0];
  return(retval);

}




/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
			    void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int), 
			    FTYPE (*res_func) (FTYPE []) )
{
  FTYPE f, f_old, df, df_old, dx[NEWT_DIM], dx_old[NEWT_DIM], x_old[NEWT_DIM], x_older[NEWT_DIM], resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, errx_old, errx_oldest, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE randtmp, tmp;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;
  FTYPE resid_norm, resid_check, grad_check;

  FTYPE res_func_val, res_func_old, res_func_new;
  FTYPE dn[NEWT_DIM], del_f[NEWT_DIM];

  // Jon's modification
  FTYPE residold[NEWT_DIM],xold[NEWT_DIM],dxold[NEWT_DIM],DAMPFACTOR[NEWT_DIM];
  static int newt_errorcheck(FTYPE errx, FTYPE x0, FTYPE wglobal);
  static int newt_extracheck(FTYPE errx, FTYPE x0, FTYPE wglobal);


  static void my_lnsrch(int, FTYPE [], FTYPE, FTYPE [], FTYPE [], FTYPE [], FTYPE *, 
			  FTYPE, FTYPE, int *, FTYPE (*res_func) (FTYPE []));

  static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  ) ;

  int   keep_iterating, i_increase, retval2,retval = 0;
  const int ltrace  = 0;
  const int ltrace2 = 1;
  int it;
  static void pick_validate_x(int eomtype);


  // pick version of validate_x depending upon scheme
  pick_validate_x(geomtype);


  retval = 0;


  errx = 1. ; 
  errx_old = 2.;
  df = df_old = f = f_old = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_older[id]=x_old[id] = x_orig[id] = x[id] ; // save initial guess as orig and old
  vsq_old = vsq = W = W_old = 0.;

  for( id = 0; id < n ; id++) DAMPFACTOR[id]=1.0; // Jon's mod: start out fast
  for( id = 0; id < n ; id++) dx_old[id]=1E30*x[id];// Assume bad error at first
  n_iter = 0;




  //////////////////////////////////////////////////////////
  //
  /* Start the Newton-Raphson iterations : */
  //
  /////////////////////////////////////////////////////////
  keep_iterating = 1;
  while( keep_iterating ) { 
     nstroke++;
     lntries++;


#if(CRAZYDEBUG)
  if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
    for(it=0;it<n;it++) dualfprintf(fail_file,"before funcd: x[%d]=%21.15g dx[%d]=%2.15g\n",it,x[it],it,dx[it]);
  }
#endif      

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */



#if(CRAZYDEBUG)
  if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
    for(it=0;it<n;it++) dualfprintf(fail_file,"after funcd: x[%d]=%21.15g dx[%d]=%2.15g\n",it,x[it],it,dx[it]);
  }
#endif      


#if(!OPTIMIZED)
    /*  Check for bad untrapped divergences : */
    if( (finite(f)==0) ||  (finite(df)==0) ) {
      if( debugfail >= 2 ) { 
	dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
	dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
      }
      return(1);
    }
#endif


#if(!OPTIMIZED)
    /* Randomly rescale Newton step to break out of iteration cycles: */
    if( ((n_iter+1) % CYCLE_BREAK_PERIOD) == 0 ) {
      randtmp = ( (1.*rand())/(1.*RAND_MAX) );
      for( id = 0; id < n ; id++) dx[id] *= randtmp;
      //	for( id = 0; id < n ; id++) dx[id] *= ( (1.*rand())/(1.*RAND_MAX) );
    }
#endif

    /* Save old values before calculating the new: */
    errx_oldest = errx_old;
    errx_old = errx;
    lerrx=errx;
    errx = 0.;
    f_old = f;
    for( id = 0; id < n ; id++) {
      x_older[id]=x_old[id];  //x_older contains last step's W used by above funcd() to get resid and dx
      x_old[id] = x[id] ; // x_old contains just used W that was used to compute funcd() to get present resid and dx
    }

    /* Make the newton step: */
    if( do_line_search == 1 ) { 

      /* Compare the residual to its initial value */ 
      if( n_iter == 0 ) { 
	resid_norm = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_norm += fabs(resid[id]);
	}
	resid_norm /= 1.0*n ;
	if( resid_norm == 0.0 ) resid_norm = 1.0;
      }

      for( id = 0; id < n ; id++) {
	tmp = 0.;
	for( jd = 0; jd < n ; jd++) {
	  tmp += jac[jd][id] * resid[jd];
	}
	del_f[id] = tmp;
      }
      for( id = 0; id < n ; id++) {
	dn[id] = dx[id];
      }

      my_lnsrch(n, x_old-1, f_old, del_f-1, dn-1, x-1, &f, TOL_LINE_STEP, SCALEMAX, &retval, res_func);

      /* dx is needed for errx calculation below: */
      for( id = 0; id < n ; id++) {
	dx[id] = x[id] - x_old[id];
      }

#if(!OPTIMIZED)
      if( ltrace ) { 
	res_func_val = res_func(x);
	res_func_old = res_func(x_old);
	dualfprintf(fail_file,"gnr(): f_old, f, res_func_old, res_func_val = %21.15g  %21.15g  %21.15g  %21.15g  \n",
		f_old, f, res_func_old, res_func_val );
	dualfprintf(fail_file,"gnr(): x_old = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",x_old[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): x     = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",x[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): dn    = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",dn[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): del_f = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",del_f[id]);
	}
	dualfprintf(fail_file,"\n ");
      }
#endif

      /* Check to see if line search problem is because the residual vector is already small enough */
      if( retval == 1 ) {
	resid_check = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_check += fabs(resid[id]);
	}
	resid_check /= 1.0*n;
	
	if( resid_check <= resid_norm * NEWT_FUNC_TOL ) {
	  retval = 0;
	}
	if( ltrace && retval ) { 
	  dualfprintf(fail_file,"general_newton_raphson():  retval, resid_check = %4i  %21.15g \n",retval, resid_check);
	}	  
      }
      /* If initial Newton step is bad, then try again without line searching: */
      if( (retval == 2) && (USE_LINE_SEARCH == do_line_search) ) { 
	if(debugfail>=2) dualfprintf(fail_file,"bad first step\n");
#if(!OPTIMIZED)	  
	if( ltrace ) { 
	  dualfprintf(fail_file,"gnr(): bad first step: retval, f_old, f  = %4i  %21.15g  %21.15g  \n",retval,f_old,f);
	  dualfprintf(fail_file,"gnr: doing recursive call, retval, errx = %4i  %21.15g \n", retval, errx );
	}
#endif
	retval = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
	for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
	return( retval );
      }

      /* Check to see if it is trapped in a local minimum, i.e. gradient is too small */ 
      if( retval == 1 ) { 
	grad_check = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_check = (x[id] == 0.) ? 1.0 : fabs(x[id]) ;
	  grad_check  +=  del_f[id] * resid_check ;
	}
	resid_check = (f == 0.) ? 1.0 : fabs(f) ;
	grad_check /= resid_check;
	
	/* Then we've most likely found a solution: */
	if( grad_check > GRADMIN ) { 
	  retval = -1;
	}
	else if( ltrace ) { 
	  dualfprintf(fail_file,"general_newton_raphson():  retval, grad_check = %4i  %21.15g \n",retval, grad_check);
	}
      }
    }
    else {
      /* don't use line search : */

#if(JONGENNEWTSTUFF)
      if(n_iter==0){
	for(id=0;id<n;id++) {
	  xold[id]=x[id]; // original W
	  dxold[id]=dx[id]; // original dW
	  residold[id]=resid[id]; // original residual
	}
      }
      else{ // only meaningful to consider goodness of W if compute residual from present W, not previous W.
	for(id=0;id<n;id++){
	  if(sign(residold[id])!=sign(resid[id])){
	    //  then passed over root -- bad in very relativistic case
	    if(fabs(xold[0])>GAMMASQCHECKRESID || fabs(x_old[0])>GAMMASQCHECKRESID){
	      // then ultrarelativistic and need to be more careful
	      DAMPFACTOR[id]=0.5*DAMPFACTOR[id];
	      x[id]=x_older[id]; // since current W is what gave bad residual and already saved it into x_old[0]
	      dx[id]=dx_old[id]; // new dx is probably messed up (didn't yet save dx[0], so use dx_old[0]
	      //	      dualfprintf(fail_file,"resetW: id=%d :: %21.15g %21.15g\n",id,x[id],dx[id]);
	    }
	  }
	}
      }
#endif
      
      for( id = 0; id < n ; id++) {
	x[id] += DAMPFACTOR[id]*dx[id]  ;
      }
#if(CRAZYDEBUG)
    if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"no line search: x[0]=%21.15g dx[0]=%21.15g\n",x[0],dx[0]);
    }
#endif

    } /* End of "to do line search or not to do line search..." */


    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific: 

#if( NEWCONVERGE == 1 )



#if(JONGENNEWTSTUFF==0)
    // SUPERGODMARK: is wglobal always the right normalization for the error?

    //    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);
    //    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/wglobal);
    
    errx  = (wglobal==0.) ?  fabs(dx[0]) : fabs(dx[0]/MAX(fabs(wglobal),fabs(x[0])));

#elif(JONGENNEWTSTUFF==1)
    // makes no sense to compute error from present dx for present x.  The error in x[0] is unknown.  errx here will be error in xold
    // actually error of previous x[0]
    // in this way, final funcd() call is actually an error check only, not an update to x[0]!  So final solution returned should be x_old
    errx  = (wglobal==0.) ?  fabs(dx_old[0]) : fabs(dx_old[0]/MAX(fabs(wglobal),fabs(x[0])));
    //    errx  = (wglobal==0.) ?  fabs(dx_old[0]) : fabs(dx_old[0]/fabs(wglobal));
#endif




#if(CRAZYDEBUG)
    if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"x[0]=%21.15g dx[0]=%21.15g\n",x[0],dx[0]);
    }
#endif


    /* For the old criterion, look at errors in each indep. variable(s) (except for 5D) : */
#else
    for( id = 0; id < n ; id++) {
      errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]);
    }
    errx /= 1.*n;
#endif


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    // METHOD specific:

    
    validate_x( x, x_old ) ;


    /****************************************/
    /* Check to see if we're in a infinite loop with error function: */
    /****************************************/
#if( CHECK_FOR_STALL )
    if( ( (errx_old == errx) || (errx_oldest == errx) ) && (errx <= MIN_NEWT_TOL) )  errx = -errx;
#endif 

    /****************************************/
    /* If there's a problem with line search, then stop iterating: */
    /****************************************/
    if( (retval == 1) || (retval == -1) ) errx = -errx;


#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file," general_newton_raphson(): niter,f_old,f,errx_old,errx = %4i  %21.15g  %21.15g  %21.15g  %21.15g\n",  
	      n_iter,f_old,f,errx_old,errx );
      dualfprintf(fail_file,"gnr(): x_old = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",x_old[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): x     = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",x[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): dx     = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",dx[id]);
      }
      dualfprintf(fail_file,"\n ");
      
    }
#endif

    /****************************************/
    /* Prepare for the next iteration, set the "old" variables: */
    /****************************************/
    for( id = 0; id < n ; id++)  dx_old[id] = dx[id] ;
    f_old  = f;
    df_old = df;


    /****************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations before stopping */
    /****************************************/
   

    
    if(newt_errorcheck(errx, x[0],  wglobal) && (doing_extra == 0) && (newt_extracheck(errx,x[0],wglobal) > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( (newt_errorcheck(errx, x[0],  wglobal)&&(doing_extra == 0)) || (i_extra > newt_extracheck(errx,x[0],wglobal)) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

#if(CRAZYDEBUG||1)
    //i=0 j=63 part=2 step=9 :: n_iter=6
    //  if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
    //    if(icurr==0 && jcurr==63 && nstep==9 && steppart==2){
    //    if(n_iter>20){
    //if(icurr>90 && CRAZYNEWCHECK){
    if(CRAZYNEWCHECK){
    //    if(icurr==0 && jcurr==63 && nstep==12 && steppart==0){
      dualfprintf(fail_file,"i=%d j=%d part=%d step=%d :: n_iter=%d :: errx=%21.15g minerr=%21.15g :: x[0]=%21.15g dx[0]=%21.15g wglobal=%21.15g\n",icurr,jcurr,steppart,nstep,n_iter,errx,MIN_NEWT_TOL,x[0],dx[0],wglobal);
    }
    //  }
#endif

  }   // END of while(keep_iterating)

  
    /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) ||  (finite(df)==0) || (finite(x[0])==0) || (finite(x[1])==0)) {
    if(1||debugfail>=2) dualfprintf(fail_file,"naninf\n");
#if(!OPTIMIZED)
    if( debugfail >= 2 ) { 
      dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
      dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
    }
#endif
    return(1);
  }


  ///////////////////////////
  //
  /// CHECK ERROR
  //
  /////////////////////////

  // assign final result
#if(JONGENNEWTSTUFF==1)
    for(id=0;id<n;id++) {
      x[id]=x_old[id]; // actually evaluated x
    }
#endif


  if(newt_errorcheck(errx, x[0],  wglobal) ){
#if(DOHISTOGRAM)
    bin_newt_data( errx, n_iter, 2, 0 );
#endif
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   2   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra, errx); 
      
    }
    //    dualfprintf(fail_file,"gnr retval5 = %4i \n", 0); 
#endif


    return(0);
  }

#if(0) // stupid code
  // check if ultrarel, the relax error tolerance
  if(
     (fabs(errx) >= NEWT_TOL) && // larger than nairvely expected
     ( (fabs(wglobal)>fabs(GAMMASQCHECKRESID)) || (fabs(x[0])>fabs(GAMMASQCHECKRESID)) )&& // but ultrarelativistic, so relax requirements
      //      fabs(errx*errx)<= MINERROREXPECTED*MINERROREXPECTED*fabs(x[0]))
     (fabs(errx)<= MINERROREXPECTED*fabs(x[0]))
     ){


    return(0);
  }
#endif


  if( (fabs(errx) <= MIN_NEWT_TOL) && newt_errorcheck(errx, x[0],  wglobal) ){
#if(DOHISTOGRAM)
    bin_newt_data( errx, n_iter, 1, 0 );
#endif
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   1   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      
    }
    if(ltrace) {
      dualfprintf(fail_file,"general_newton_raphson(): found minimal solution \n");
      
    }
    //    dualfprintf(fail_file,"gnr retval4 = %4i \n", 0); 
#endif


#if(JONGENNEWTSTUFF==1)
    for(id=0;id<n;id++) {
      x[id]=x_old[id]; // actually evaluated x
    }
#endif

    return(0);
  }



  if( fabs(errx) > MIN_NEWT_TOL){
    if( (do_line_search != USE_LINE_SEARCH) || (USE_LINE_SEARCH < 0) ) { 
#if(DOHISTOGRAM)
      bin_newt_data( errx, n_iter, 0, 0 );
#endif

      if(debugfail>=2) dualfprintf(fail_file,"fabs(errx)=%21.15g > MIN_NEWT_TOL=%21.15g n=%d\n",fabs(errx),MIN_NEWT_TOL,n);
#if(!OPTIMIZED)

      if(ltrace2) {
	dualfprintf(fail_file," totalcount = %d   0   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      }
      if(ltrace) {
	dualfprintf(fail_file,"general_newton_raphson():  did not find solution \n");
	if( retval == -1 ) {
	  dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x = ");
	  for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x[id]);
	  dualfprintf(fail_file,"\n");
	  dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x_old = ");
	  for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x_old[id]);
	  dualfprintf(fail_file,"\n");
	}
	
      }
      //      dualfprintf(fail_file,"gnr retval2 = %4i \n", 1); 
#endif
      return(1);
    } 
    else {
      /* If bad return and we tried line searching, try it without before giving up: */
      //      dualfprintf(fail_file,"gnr: doing recursive call, do_line_search, retval, errx = %4i  %4i  %21.15g \n", do_line_search, retval, errx );
      //      
      retval2 = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
      for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
      //      dualfprintf(fail_file,"gnr retval3 = %4i \n", retval2); 

      if(debugfail>=2) dualfprintf(fail_file,"fabs(errx) > MIN_NEWT_TOL\n");


      return( retval2 );
    }
  }




#if(!OPTIMIZED)
  dualfprintf(fail_file,"gnr retval6 = %4i \n", 0);
#endif
  return(0);

}



static int newt_errorcheck(FTYPE errx, FTYPE x0, FTYPE wglobal)
{

  int ultrarel=0;

  ultrarel=( (fabs(wglobal)>fabs(GAMMASQCHECKRESID)) || (fabs(x0)>fabs(GAMMASQCHECKRESID)) );

  if(!ultrarel){
    return(fabs(errx) <= NEWT_TOL);
  }
  else{
    return(fabs(errx) <= NEWT_TOL_ULTRAREL);
  }
}

static int newt_extracheck(FTYPE errx, FTYPE x0, FTYPE wglobal)
{

  int ultrarel=0;

  ultrarel=( (fabs(wglobal)>fabs(GAMMASQCHECKRESID)) || (fabs(x0)>fabs(GAMMASQCHECKRESID)) );

  if(!ultrarel){
    return(EXTRA_NEWT_ITER);
  }
  else{
    return(EXTRA_NEWT_ITER_ULTRAREL);
  }
}


#define ALF 1.0e-4

static void my_lnsrch( int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], 
		FTYPE *f, FTYPE TOLX, FTYPE stpmax, int *check, FTYPE (*func) (FTYPE []) )
{
  int i;
  FTYPE a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  int icount=0;
  int bad_step;
  FTYPE bad_step_factor = 2.0;
  
  const int ltrace = 0;
	
  
  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for (i=1;i<=n;i++) {
    //    temp=fabs(p[i])/max(fabs(xold[i]),1.0);
    temp= (xold[i] == 0.) ? fabs(p[i]) :  fabs(p[i]/xold[i]);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;

#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"my_lnsrch(): sum, slope, test, alamin =   %21.15g  %21.15g  %21.15g  %21.15g \n",sum,slope,test,alamin); 
  }
#endif

  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];

    // METHOD specific:
    validate_x( (x+1), (xold+1) );

    *f=(*func)(x+1);
      
    bad_step = 0;

    if( finite(*f)==0 ) { 
      bad_step = 1;
    }


    //      if( bad_step ) alam /= bad_step_factor;
    //      if (alam < alamin) bad_step = 0;

    if( bad_step ) { 
      *check = 2;
#if(!OPTIMIZED)
      dualfprintf(fail_file,"my_lnsrch(): bad_step = 1,  f = %21.15g , alam = %21.15g \n", *f, alam); 
      for( i = 1; i <= n; i++ ) { 
	dualfprintf(fail_file,"my_lnsrch(): (xold[%d], x[%d], p[%d]) = %21.15g , %21.15g , %21.15g \n", i,i,i, xold[i],x[i],p[i]);
      }
#endif
      return;
    }
      
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
#if(!OPTIMIZED)
      if( ltrace ) { 
	dualfprintf(fail_file,"my_lnsrch(): alam < alamin: alam, alamin = %21.15g  %21.15g \n", alam,alamin); 
      }
#endif
      return;
    } 
    else if (*f <= fold+ALF*alam*slope) {
#if(!OPTIMIZED)
      if( ltrace ) { 
	dualfprintf(fail_file,"my_lnsrch(): good exit:  alam, alamin, f, fold = %21.15g  %21.15g %21.15g  %21.15g \n", alam,alamin, *f, fold); 
      }
#endif
      return;
    }
    else {
      if (alam == 1.0) {
	tmplam = -slope/(2.0*(*f-fold-slope));
#if(!OPTIMIZED)
	if( ltrace ) {
	  dualfprintf(fail_file,"my_lnsrch(): setting tmplam!!    tmplam, alam =  %21.15g  %21.15g !!\n", tmplam, alam); 
	}
#endif
      }
      else {
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold2-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc<0.0) {
#if(!OPTIMIZED)
	    if( disc < -1.e-10 ) {
	      if( ltrace ) dualfprintf(fail_file,"my_lnsrch(): Big Roundoff problem:  disc = %21.15g \n", disc);
	    }
#endif	      
	    disc = 0.;
	  }
	  else tmplam=(-b+sqrt(disc))/(3.0*a);
	}
	if (tmplam>0.5*alam)
	  tmplam=0.5*alam;
#if(!OPTIMIZED)
	if( ltrace ) {
	  dualfprintf(fail_file,"my_lnsrch(): rhs1, rhs2, a, b, tmplam, alam =  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g !!\n",
		  rhs1, rhs2, a, b, tmplam, alam );
	}
#endif
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=max(tmplam,0.1*alam);
#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file,"my_lnsrch(): icount, alam, alam2, tmplam  =   %4i  %21.15g  %21.15g  %21.15g \n",
	      icount, alam, alam2, tmplam);
    }
#endif    
    icount++;
  }
}
#undef ALF


#define N_CONV_TYPES 3
#define N_NITER_BINS (MAX_NEWT_ITER + 3)
#define NBINS 200

/**************************************************** 

  Used to gather statistics for the Newton solver during a disk run:

  -- bins are inclusive on their minimum bounds, and exlusive on their max.


      lerrx_[min,max] = the min./max. bounds of the errx historgram
      d_errx          = errx histogram bin width
      n_beyond_range  = number of solutions with errx max. beyond range speciified by lerrx_[min,max]
      n_conv_types    = number of types of ways newton procedure can exit;

      print_now       = 1 if you want to just print the histograms without counting a solution;

*****************************************************/

static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  ) 
{

  
  /* General variables */
  int i, j;
  static int first_call   = 1;
  static long int n_bin_calls  = 0L;
  const  long int n_output_freq = 128*128*100;
  FTYPE lerrx;

  /* Variables for the errx histogram */ 
  static const FTYPE lerrx_min = -23.;
  static const FTYPE lerrx_max =   4.;
  static FTYPE d_errx;
  static long int n_errx[N_CONV_TYPES][NBINS];
  static FTYPE xbin[NBINS];
  static long int n_beyond_range = 0L;   /* Number of points that lie out of the bounds of our errx histogram */
  int ibin;

  /* Variables for the histogram of the number of newton iterations : */
  static long int n_niters[N_CONV_TYPES][N_NITER_BINS];
  

  /* Clear arrays, set constants : */
  if( first_call ) {
    d_errx    = ((lerrx_max - lerrx_min)/(1.*NBINS));
    
    for( i = 0; i < N_CONV_TYPES; i++ ) { 
      for( j = 0; j < NBINS; j++ ) { 
	n_errx[i][j] = 0;
      }
      for( j = 0; j < N_NITER_BINS; j++ ) { 
	n_niters[i][j] = 0;
      }
    }

    for( j = 0; j < NBINS; j++ ) { 
      xbin[j] = lerrx_min + (1.*j + 0.5)*d_errx;
    }
      
    first_call = 0;
  }


  if( print_now != 1 ) {

    /* Check validity of arguments : */
    errx = fabs(errx) ;
    lerrx = log10(errx + 1.0e-2*pow(10.,lerrx_min));

    if( (niters < 0) || (niters >= N_NITER_BINS) ) {
      dualfprintf(fail_file,"bin_newt_data(): bad value for niters = %d \n", niters );
      fflush(stdout);
      return;
    }
    
    /* Determine ibin */
    if( lerrx < lerrx_min ) {
      ibin = 0;
      n_beyond_range++ ;
    }
    else if( lerrx >= lerrx_max ) {
      ibin = NBINS - 1;
      n_beyond_range++ ;
    }
    else {
      ibin = (int) ( ( lerrx - lerrx_min ) / d_errx );
    }
		  
    /* Tally this solution  */
    n_errx[ conv_type][ibin]++;
    n_niters[conv_type][niters]++;

  }


  /* Print out the histograms periodically or when asked to : */
  if( print_now  ||  ( (n_bin_calls % n_output_freq) == 0 ) ) {

    dualfprintf(log_file,"t = %21.15g ,  n_beyond_range = %ld , n_bin_calls = %ld \n", t, n_beyond_range, n_bin_calls);

    /* ERRX */
    dualfprintf(log_file,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");
    dualfprintf(log_file,"                         x");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      dualfprintf(log_file,"            N%d",j);
    }
    dualfprintf(log_file,"\n");
    
    for( i = 0; i < NBINS; i++ ) { 
      dualfprintf(log_file,"%21.15g ",xbin[i]);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
	dualfprintf(log_file,"%13ld ", n_errx[j][i]);
      }
      dualfprintf(log_file,"\n");
    }    
    dualfprintf(log_file,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");


    /* NITER */

    dualfprintf(log_file,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");
    dualfprintf(log_file,"        niter");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      dualfprintf(log_file,"            N%d",j);
    }
    dualfprintf(log_file,"\n");
    
    for( i = 0; i < N_NITER_BINS; i++ ) { 
      dualfprintf(log_file,"%13d ", i);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
	dualfprintf(log_file,"%13ld ", n_niters[j][i]);
      }
      dualfprintf(log_file,"\n");
    }    
    dualfprintf(log_file,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");

  }
    
    
  if( print_now != 1 )  n_bin_calls++;


  return;

}
 
 
#undef N_CONV_TYPES 
#undef N_NITER_BINS 
#undef NBINS 

#undef CRAZYDEBUG

