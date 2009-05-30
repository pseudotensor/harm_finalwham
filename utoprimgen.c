


#include "decs.h"



// fractional error above which inversion is reported on as having made a signficant error
#if(PRODUCTION)
#define CHECKONINVFRAC (1E-1)
#else
#define CHECKONINVFRAC (1E-2)
#endif

// whether to fail if check on inversion fails
#define FAILIFBADCHECK 1

// fraction upon if greater will treat as failure if FAILIFBADCHECK is triggered
#define CHECKONINVFRACFAIL (1E-1)



int Utoprimgen(int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr)
{
  // debug
  int i, j, k;
  FTYPE Ugeomfree[NPR],Ugeomfree0[NPR];
  FTYPE pr0[NPR];
  FTYPE prother[NPR];
  FTYPE Upr[NPR],Uprother[NPR];
  int whichentropy;
  struct of_state q;
  FTYPE Uold[NPR],Unew[NPR];
  FTYPE fdiff[NPR];
  int otherfail;
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  extern int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  extern int Utoprim_coldgrmhd(FTYPE *U, struct of_geom *geom, FTYPE *pr, int *positivityproblem);
  int positivityproblem;
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr);
  int Utoprimgen_compare(int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr);
  int Utoprimgen_tryagain(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, int *lpflag, FTYPE *pr0, FTYPE *pr);
  int lpflag,hotpflag,coldpflag;
  FTYPE prhot[NPR],prcold[NPR];
  int trycoldinversion(int hotpflag, FTYPE *pr0, FTYPE *pr, FTYPE *Ugeomfree, FTYPE *Ugeomfree0, struct of_geom *ptrgeom);
  int badinversion,badinversionfail;


  //  if(U[ENTROPY]<0.0){
  //    dualfprintf(fail_file,"BADENTROPY: i=%d j=%d nstep=%ld steppart=%d : bad U[ENTROPY]=%21.15g\n",ptrgeom->i,ptrgeom->j,nstep,steppart,U[ENTROPY]);
  //  }

  // Notice that reconstruct "standard" geometry-free conserved quantities by first modifying the geometry prefactor and THEN dealing with rest-mass or other subtractions and additions.
  // This is consistent with primtoflux() and how primtoflux() is used in source_conn().
  UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);
  // e.g. if REMOVERESTMASSFROMUU==1 and inputtype==UEVOLVE, then Ugeomfree has energy T^t_t that includes rest-mass

  // backup
  PLOOP(k){
    Uold[k]=Ugeomfree0[k]=Ugeomfree[k];
    pr0[k]=pr[k];
  }

  


  ///////////////////////////////////////////////////
  //
  ///////////// Setup U[UU] with or without rest-mass for different inversions and setting of REMOVERESTMASSFROMUU
  //
  // check if inversion is chosen consistent with REMOVERESTMASSFROMUU
  // At this point, Ugeomfree[UU] in "standard form with rest-mass" for REMOVERESTMASSFROMUU=0,1.  If ==2, then no rest-mass
  //
  ///////////////////////////////////////////////////

#if(REMOVERESTMASSFROMUU==2)
  // 5D1 handles REMOVERESTMASSFROMUU internally
  // NONRELCOMPAT only accepts UU without rest-mass, so good to go
  if( UTOPRIMVERSION!=UTOPRIM5D1 && (UTOPRIMVERSION!=UTOPRIMJONNONRELCOMPAT) ){
    //    dualfprintf(fail_file,"This code does not handle REMOVERESTMASSFROMUU==2: other Utoprim's besides 5D1 and 2DFINAL\n");
    //    myexit(1);
    // then change conserved quantity so can work
    Ugeomfree0[UU]-=Ugeomfree0[RHO]; // "add in rest-mass" so in correct form
    Ugeomfree[UU]-=Ugeomfree[RHO]; // "add in rest-mass" so in correct form
  }
#elif(REMOVERESTMASSFROMUU==0 || REMOVERESTMASSFROMUU==1)
  if(UTOPRIMVERSION==UTOPRIMJONNONRELCOMPAT){
    // then change conserved quantity so can work for UTOPRIMJONNONRELCOMPAT
    Ugeomfree0[UU]+=Ugeomfree0[RHO]; // "subtract out rest-mass" so in correct form for NONRELCOMPAT
    Ugeomfree[UU]+=Ugeomfree[RHO]; // "subtract out rest-mass" so in correct form for NONRELCOMPAT
  }
#endif






  ////////////////////////
  //
  // DO INVERSION

  // assume solution is good unless flagged as bad
  pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL]=UTOPRIMNOFAIL;


  if(EOMTYPE==EOMGRMHD){

    ///////////////////////////////////////////////////
    //
    ///////////// HOT GRMHD
    //
    ///////////////////////////////////////////////////

    if(DOENTROPY!=DOEVOLVEDIRECTENTROPY){
      // then do energy equation version
    
      if(UTOPRIMVERSION!=UTOPRIMCOMPARE) Utoprimgen_pick(UTOPRIMVERSION, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr);
      else Utoprimgen_compare(EVOLVENOENTROPY,Ugeomfree,ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr);

      // try other methods (assumes all methods can handle WHICHVEL, etc. used for primary model)
      // right now, all can handle WHICHVEL==VELREL4 and energy equation evolution and REMOVERESTMASSFROMUU=0,1
#if((WHICHVEL==VELREL4)&&(REMOVERESTMASSFROMUU<=1)&&(UTOPRIMTRYAGAIN))
      Utoprimgen_tryagain(EVOLVENOENTROPY, Ugeomfree0, Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr0, pr);
#endif
    }
    else{ // direct entropy evolution

      // do inversion for entropy version of EOMs
      // only one inversion is setup to handle this
      PLOOP(k) prother[k]=pr0[k];
      whichentropy=EVOLVEFULLENTROPY;

      ////////////////////////
      // get entropy evolution (don't use failure -- otherfail)
      MYFUN(Utoprim(whichentropy,Uold, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr),"step_ch.c:Utoprimgen()", "Utoprim", 1);

    } // end if doing direct entropy evolution


    ////////////////////
    //
    //  If hot GRMHD failed or gets suspicious solution, revert to cold GRMHD if solution is cold
    //
    ///////////////////

#if(HOT2COLD)
    hotpflag=pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL];

    if(hotpflag){
      trycoldinversion(hotpflag, pr0, pr, Ugeomfree, Ugeomfree0, ptrgeom);
    }// end if hotpflag
#endif



  }
  else if(EOMTYPE==EOMCOLDGRMHD){

    ///////////////////////////////////////////////////
    //
    ///////////// COLDGRMHD
    //
    ///////////////////////////////////////////////////


    if(1){
      // Jon's inversion
      Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr);
    }

    if(0){ // not working yet
      MYFUN(Utoprim_coldgrmhd(Ugeomfree, ptrgeom, pr,&positivityproblem),"step_ch.c:Utoprimgen()", "Utoprim_coldgrmhd", 1);
      //Utoprim_coldgrmhd(Ugeomfree, ptrgeom, pr,&positivityproblem);
    }



  }
  else if(EOMTYPE==EOMFFDE){

    ///////////////////////////////////////////////////
    //
    ///////////// FORCE FREE
    //
    ///////////////////////////////////////////////////


    // GODMARK: inversions lead to different behavior!

    
    if(0){ // Jon's old inversion
      MYFUN(Utoprim_ffde(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_ffde", 1);
    }
    //    else{ // Jon's new inversion

    if(1){
      // Jon's inversion
      Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr);
    }

    // compare
    if(0){

      MYFUN(Utoprim_ffde(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_ffde", 1);

      PLOOP(k) Ugeomfree[k]=Ugeomfree0[k]; // make sure good conserved quantity
      
      Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], prother);

      // now get conserved quantity from pr to check
      // find U(p)
      MYFUN(get_state(pr, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
      MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Upr),"step_ch.c:advance()", "primtoU()", 1);

      MYFUN(get_state(prother, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
      MYFUN(primtoU(UNOTHING,prother, &q, ptrgeom, Uprother),"step_ch.c:advance()", "primtoU()", 1);


      // now compare

      //PLOOP(k)
      dualfprintf(fail_file,"nstep=%ld stepart=%d i=%d j=%d\n",nstep,steppart,ptrgeom->i,ptrgeom->j);
      for(k=U1;k<B3;k++){
	/*
	  if(
	  ((fabs(Ugeomfree[k]-Upr[k])/(fabs(Ugeomfree[k])+fabs(Upr[k])+SMALL)) > 1E-12)||
	  ((fabs(Ugeomfree[k]-Uprother[k])/(fabs(Ugeomfree[k])+fabs(Uprother[k])+SMALL)) > 1E-12)
	  ){
	  dualfprintf(fail_file,"DIFF: %21.15g %21.15g :: k=%d Ugeomfree=%21.15g prold=%21.15g prnew=%21.15g Upr=%21.15g Uprother=%21.15g\n",(fabs(Ugeomfree[k]-Upr[k])/(fabs(Ugeomfree[k])+fabs(Upr[k])+SMALL)),(fabs(Ugeomfree[k]-Uprother[k])/(fabs(Ugeomfree[k])+fabs(Uprother[k])+SMALL)),k,Ugeomfree[k],pr[k],prother[k],Upr[k],Uprother[k]);
	  }
	*/
	if(
	   ((fabs(pr[k]-prother[k])/(fabs(pr[k])+fabs(prother[k])+SMALL)) > 1E-12)&&( (fabs(pr[k])>1E-15)||(fabs(prother[k])>1E-15) )
	   ){
	  dualfprintf(fail_file,"DIFF: %21.15g :: k=%d Ugeomfree=%21.15g prold=%21.15g prnew=%21.15g Upr=%21.15g Uprother=%21.15g\n",(fabs(pr[k]-prother[k])/(fabs(pr[k])+fabs(prother[k])+SMALL)),k,Ugeomfree[k],pr[k],prother[k],Upr[k],Uprother[k]);
	}
      }
    }

      //    }
    //Utoprim_ffde(Ugeomfree, ptrgeom, pr);
  }



  ///////////////////////////////////////////////////
  //
  ///////////// Report failure
  //
  // complain if pflag set
  // only output if failed
  //
  ///////////////////////////////////////////////////

  // for now only report if not just negative density failure
  lpflag=pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL];

  if( lpflag==UTOPRIMFAILRHONEG || lpflag==UTOPRIMFAILUNEG || lpflag==UTOPRIMFAILRHOUNEG ){
    // then don't report info for now GODMARK
  }
  else if(lpflag &&(debugfail>=1)){
    dualfprintf(fail_file, "Failed to find a prim. var. solution!! t=%21.15g steppart=%d nstep=%ld i=%d j=%d k=%d : fail=%d : errx=%21.15g\n",t,steppart,realnstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag,lerrx);
  }






  ///////////////////////////////////////////////////
  //
  ///////////// LOTS OF DEBUG STUFF
  //
  ///////////////////////////////////////////////////





#define CRAZYDEBUG 0
#if(CRAZYDEBUG) // begin crazy debug stuff
  if(1|| lntries>3){
    if(nstep==9 && steppart==2 && ptrgeom->i==0 && ptrgeom->j==31){
    //    if(nstep==0 && steppart==0 && ptrgeom->i==4 && ptrgeom->j==36){
    //  if(nstep==0 && steppart==0 && ptrgeom->i == 5 && ptrgeom->j == 31 ){
    dualfprintf(fail_file,"nstep=%ld stepart=%d :: i=%d j=%d :: lntries=%d\n",nstep,steppart,ptrgeom->i,ptrgeom->j,lntries);
    
    //    PLOOP(k) dualfprintf(fail_file,"Ugeomfree[%d]=%21.15g pr[%d]=%21.15g\n",k,Ugeomfree[k],k,pr[k]);
    PLOOP(k) dualfprintf(fail_file,"Uoldgeomfree[%d]=%21.15g Uold[%d]=%21.15g pr[%d]=%21.15g (pr0[%d]=%21.15g)\n",k,Uold[k],k,Uold[k]*ptrgeom->g,k,pr[k],k,pr0[k]);

    dualfprintf(fail_file,"g=%21.15g\n",ptrgeom->g);

    DLOOP(j,k) dualfprintf(fail_file,"gcon=%21.15g gcov=%21.15g\n",ptrgeom->gcov[j][k],ptrgeom->gcon[j][k]);
    //    myexit(777);
    //  }
    }
  }


  // test inversion for accuracy in conservative space
  //  if(nstep==0 && steppart==0 && ptrgeom->i == 4 && ptrgeom->j == 36 ){
  if(nstep==9 && steppart==2 && ptrgeom->i==0 && ptrgeom->j==31){

#if(0)
pr[0]= 7.22714301361038e-06 ;
pr[1]=-2.55753797775927e-07 ;
pr[2]= 0.000892168434512681 ;
pr[3]=  -0.0349621334882251 ;
pr[4]=                   -0 ;
pr[5]= -0.00101880685475446 ;
pr[6]=   0.0399035382308458 ;
pr[7]=                    0 ;
#elif(0)
pr[0]= 7.46157819677347e-06;
pr[1]=  7.3407120498644e-06;
pr[2]= 0.000367464781319951;
pr[3]=  -0.0144105143008605;
pr[4]=                    0;
pr[5]= -0.00101880685475446;
pr[6]=   0.0399035382308458;
 pr[7]=                    0;
#elif(0)
pr[0]=  7.5124289176258e-06 ;
pr[1]= 1.33752037209996e-08 ;
pr[2]= 1.33529579432262e-07 ;
pr[3]=-5.23276639757274e-06 ;
pr[4]=                    0 ;
pr[5]= -0.00101880685475444 ;
pr[6]=   0.0399035382308459 ;
 pr[7]=                    0 ;

#endif

  MYFUN(get_state(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
  MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included

  for(k=0;k<4;k++){
    dualfprintf(fail_file,"q.ucon[%d]=%21.15g q.ucov[%d]=%21.15g q.bcon[%d]=%21.15g q.bcov[%d]=%21.15g\n",k,q.ucon[k],k,q.ucov[k],k,q.bcon[k],k,q.bcov[k]);
  }

  PLOOP(k){
    Unew[k]*=ptrgeom->g;
    Uold[k]*=ptrgeom->g;
    fdiff[k] = fabs(Unew[k]-Uold[k])/(fabs(Unew[k]+Uold[k])+1E-30);
    //    if(fdiff[k]>1E-10){
      if((k>=RHO)&&(k<=B3)&&((fabs(Uold[k])>1E-20)||(fabs(Unew[k])>1E-20))){
	dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %21.15g %21.15g\n",k,fdiff[k],Uold[k],Unew[k]);
      }
      //    }
  }

  dualfprintf(fail_file,"dt=%21.15g\n",dt);


  myexit(124);
  }
#endif// end crazy debug stuff








  //////////////////////////////////////
  //
  // check up on solution
  //
  ////////////////////////////////////////
#if(CHECKONINVERSION)

  if(1|| !lpflag){ // only report if not failure
    MYFUN(get_state(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
    MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included


    badinversion=0; // assume not bad
    badinversionfail=0;
    PLOOP(k){
      // no point checking if inversion doesn't handle or is inconsistent with conservation of that quantity
      if(EOMTYPE==EOMFFDE && (k==RHO || k==UU) ) continue;
      if(EOMTYPE==EOMCOLDGRMHD  && (k==UU) ) continue;
      // leave geometry out of it
      //      Unew[k]*=ptrgeom->g;
      //      Uold[k]*=ptrgeom->g;
      fdiff[k] = fabs(Unew[k]-Uold[k])/(fabs(Unew[k])+fabs(Uold[k])+1E-30);
      if(fdiff[k]>CHECKONINVFRAC){
	if((k>=RHO)&&(k<=B3)&&((fabs(Uold[k])>1E-20)||(fabs(Unew[k])>1E-20))){
	  badinversion++;
	  dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %21.15g %21.15g\n",k,fdiff[k],Uold[k],Unew[k]);
	}
      }
      if(fdiff[k]>CHECKONINVFRACFAIL){
	if((k>=RHO)&&(k<=B3)&&((fabs(Uold[k])>1E-20)||(fabs(Unew[k])>1E-20))){
	  badinversionfail++;
	}
      }
    }
    
    if(FAILIFBADCHECK && badinversionfail){
      lpflag=UTOPRIMFAILCONVBADINVERTCOMPARE;
    }
    if(badinversion){
      dualfprintf(fail_file,"Bad inversion:\n");
      dualfprintf(fail_file,"t=%21.15g nstep=%ld stepart=%d :: i=%d j=%d :: lntries=%d lerrx=%21.15g\n",t,nstep,steppart,ptrgeom->i,ptrgeom->j,lntries,lerrx);
      PLOOP(k) dualfprintf(fail_file,"Uoldgeomfree[%d]=%21.15g Uold[%d]=%21.15g pr[%d]=%21.15g (pr0[%d]=%21.15g)\n",k,Uold[k],k,Uold[k]*ptrgeom->g,k,pr[k],k,pr0[k]);
      
      dualfprintf(fail_file,"g=%21.15g\n",ptrgeom->g);
      
      DLOOP(j,k) dualfprintf(fail_file,"gcon=%21.15g gcov=%21.15g\n",ptrgeom->gcov[j][k],ptrgeom->gcon[j][k]);

      for(k=0;k<4;k++){
	dualfprintf(fail_file,"q.ucon[%d]=%21.15g q.ucov[%d]=%21.15g q.bcon[%d]=%21.15g q.bcov[%d]=%21.15g\n",k,q.ucon[k],k,q.ucov[k],k,q.bcon[k],k,q.bcov[k]);
      }

      // only really need the below quantities to check on inversion in mathematica
      for(j=0;j<NUMGLOBALINV;j++) dualfprintf(fail_file,"%sstr=\"%21.15g\";\n",globalinvtext[j],globalinv[j]);
    }

  }


#endif





#undef CRAZYDEBUG


     
  return(0);
}


// use cold grmhd if hot grmhd gives u<0 since then cold is valid
#define USECOLDIFHOTUNEG 1
// use cold grmhd if hot grmhd leads to convergence or other "no solution" type failure
#define USECOLDIFHOTFAILCONV 1

// try cold inversion of hot one fails
int trycoldinversion(int hotpflag, FTYPE *pr0, FTYPE *pr, FTYPE *Ugeomfree, FTYPE *Ugeomfree0, struct of_geom *ptrgeom)
{
  int k;
  FTYPE prhot[NPR],prcold[NPR];
  int coldpflag;


  //  dualfprintf(fail_file,"Got here in trycoldinversion\n");

  if( hotpflag==UTOPRIMFAILRHONEG || hotpflag==UTOPRIMFAILRHOUNEG || (USECOLDIFHOTUNEG==0 && hotpflag==UTOPRIMFAILUNEG) ){
    // then maybe not so bad failure
    
  }
  else if( (USECOLDIFHOTUNEG==1 && hotpflag==UTOPRIMFAILUNEG) || (USECOLDIFHOTFAILCONV==1 && hotpflag!=0) ){
    
    
    
    // then bad failure, so try to use cold grmhd
    // restore backup in case previous inversion changed things
    PLOOP(k){
      prhot[k]=pr[k];
      Ugeomfree[k]=Ugeomfree0[k];
      pr[k]=pr0[k];
    }
    
    
    // get cold inversion
    MYFUN(Utoprim_jon_nonrelcompat_inputnorestmass(EOMCOLDGRMHD,Ugeomfree, ptrgeom, &coldpflag, prcold),"step_ch.c:Utoprimgen()", "Utoprim_2d_final_nonrelcompat_inputnorestmass", 1);
    
    
    ///////////////////////////////////
    //
    // check if cold solution is good
    if(coldpflag==UTOPRIMNOFAIL){
      // then keep cold solution
      PLOOP(k) pr[k]=prcold[k];


      // but set internal energy to previous value (should really evolve with entropy equation, but if negligible and no strong shocks, then ok )
      // GODMARK: another ok way is to set special failure that only averages internal energy!  Then can evolve at least -- in some diffusive way
      
      if(debugfail>=2) dualfprintf(fail_file,"Tried cold and good! hotpflag=%d coldpflag=%d\n",hotpflag,coldpflag);

      ///////////////////////////////
      //
      // decide how to use cold inversion solution    
      if(hotpflag==UTOPRIMFAILUNEG){
	// since cold approximation is very good, then use cold solution and just set internal energy to 0
	// if internal energy is actually small, then just set it to 0
	// works for Hubble flow!
	pr[UU]=0.0;
      }
      else{
	//////////////
	//  if internal energy is not negligible or unknown, then should average or evolve!
	pr[UU]=pr0[UU];
	// then very bad failure, so try cold inversion and average internal energy for now
	pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL]=UTOPRIMFAILU2AVG1; // assume only internal energy needs correcting by averaging
	
      }// end if (else) trying cold inversion

    }
    else{
      // then both hot and cold are bad, so keep hot
      // GODMARK: Could try going to force-free and "failing" the parallel velocity so it gets averaged like internal energy in cold case!
      // keep hotpflag and keep hot solution
      PLOOP(k) pr[k]=prhot[k];
      if(debugfail>=2) dualfprintf(fail_file,"Tried cold and bad! hotpflag=%d coldpflag=%d\n",hotpflag,coldpflag);
    }
        
  }// if bad failure of some kind

  return(0);
}





// Used for dissipation calculation only
int Utoprimdiss(int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, int *otherfail)
{
  // debug
  int i, j, k;
  FTYPE Ugeomfree[NPR],Ugeomfree0[NPR];
  FTYPE pr0[NPR];
  FTYPE prother[NPR];
  FTYPE Upr[NPR],Uprother[NPR];
  int whichentropy;
  struct of_state q;
  FTYPE Uold[NPR],Unew[NPR];
  FTYPE fdiff[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);


  // Notice that reconstruct "standard" geometry-free conserved quantities by first modifying the geometry prefactor and THEN dealing with rest-mass or other subtractions and additions.
  // This is consistent with primtoflux() and how primtoflux() is used in source_conn().
  UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);

  // backup
  PLOOP(k){
    Uold[k]=Ugeomfree0[k]=Ugeomfree[k];
    pr0[k]=pr[k];
  }

  ////////////////////////
  //
  // DO INVERSION

  // assume solution is good unless flagged as bad
  *otherfail=UTOPRIMNOFAIL;

  // do separate inversion for entropy version of EOMs
  // only one inversion is setup to handle this
  whichentropy=WHICHENTROPYEVOLVE;

  ////////////////////////
  // get entropy evolution (don't use failure -- otherfail)
  // this inversion knows about global settings for DOENTROPY and increases tolerance for doing comparison
  MYFUN(Utoprim(whichentropy, Uold, ptrgeom, otherfail, pr),"step_ch.c:Utoprimgen()", "Utoprim", 1);

  
  // now get conserved quantity from pr to check
  // find U(p)
  //  MYFUN(get_state(pr, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
  //  MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1);
  
  
  //  if(0&& *otherfail){
  //    PLOOP(k){
  //  fdiff[k] = fabs(Unew[k]-Uold[k])/(fabs(Unew[k]+Uold[k])+1E-30);
  //      if(fdiff[k]>1E-10){
  //	if((k>=U1)&&(k<=B3)&&((fabs(Uold[k])>1E-20)||(fabs(Unew[k])>1E-20))){
  //	  dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %g %g\n",k,fdiff[k],Uold[k],Unew[k]);
  //	}
  //      }
  //	}
  //      }

     
  return(0);
}





// old case study when BDIRCONT=0
// 2 9 :: 0 31 0 :: lpflag1=101 lpflag2=6 errx1=                    2 errx2= 2.05255876165841e-16
// 200 vs. 9 iteration

/* New case study (happens that BDIRCONT=1)

A) 

utoprimdiff: 2 7 :: 0 9 0 :: 1 ::  -2.27439265121865e-18    1.20924610877995e-10    1.88083520360745e-08  1.88083520383489e-08 :: 6 6 ::  1.19802964119872e-19  4.52705760710304e-17
utoprimdiff: 2 7 :: 0 47 0 :: 1 ::   1.95422412332388e-18    1.04013649723147e-10    1.87881506756607e-08  1.87881506737064e-08 :: 5 5 ::  1.19803248202617e-19  2.99507613474951e-17
utoprimdiff: 2 7 :: 0 63 0 :: 1 ::  -9.77586863333545e-19    1.36537803916274e-10    7.15982559623566e-09  7.15982559721325e-09 :: 14 200 ::  2.91804985413138e-15  1.62522245203594e-14
2 7 :: 0 63 0 :: lpflag1=0 lpflag2=0 errx1= 2.91804985413138e-15 errx2= 1.62522245203594e-14 lntries1=14 lntries2=200


*/


//////////////////////////////
//
// COMPARISON of 2 (currently only 2) methods
int Utoprimgen_compare(int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr)
{
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr);
  int pl;
  FTYPE Uo1[NPR], Uo2[NPR];
  FTYPE ptest1[NPR],ptest2[NPR];
  FTYPE Uo[NPR], po[NPR];
  int test;
  int lpflag1,lpflag2;
  int lntries1,lntries2;
  FTYPE lerrx1,lerrx2;


  // backup
  PLOOP(pl){
    ptest1[pl]=ptest2[pl]=pr[pl];
    Uo1[pl]=Uo2[pl]=Ugeomfree[pl];
  }

  // currently comparing 5D1 and LDZ
  //  Utoprimgen_pick(UTOPRIM5D1, EVOLVENOENTROPY, Uo1, ptrgeom, lpflag, ptest1);
  //  Utoprimgen_pick(UTOPRIMLDZ, EVOLVENOENTROPY, Uo2, ptrgeom, lpflag, ptest2);

  // remove rest-mass
  Uo1[UU]+=Uo1[RHO];
  Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Uo1, ptrgeom, &lpflag1, ptest1);
  lntries1=lntries; lntries=0; lerrx1=lerrx;
  Utoprimgen_pick(UTOPRIM2DFINAL, EVOLVENOENTROPY, Uo2, ptrgeom, &lpflag2, ptest2);
  lntries2=lntries; lerrx2=lerrx;

#define ERRORCONST (1E-10)

  PLOOP(pl){
    if(ptest1[pl]!=0.0) test=fabs(ptest1[pl]-ptest2[pl])/ptest1[pl]>ERRORCONST;
    else test=(ptest1[pl]-ptest2[pl])>ERRORCONST;
    if(test){
      dualfprintf(fail_file,"utoprimdiff: %d %ld :: %d %d %d :: %d ::  %21.15g   %21.15g   %21.15g %21.15g :: %d %d :: %21.15g %21.15g\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,pl,ptest1[pl]-ptest2[pl],(ptest1[pl]!=0.0) ? fabs(ptest1[pl]-ptest2[pl])/ptest1[pl] : ptest1[pl]-ptest2[pl],ptest1[pl],ptest2[pl],lntries1,lntries2,lerrx1,lerrx2);
    }
  }
  if(lpflag1!=0 || lpflag2!=0){
    dualfprintf(fail_file,"%d %ld :: %d %d %d :: lpflag1=%d lpflag2=%d errx1=%21.15g errx2=%21.15g\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag1,lpflag2,lerrx1,lerrx2);
  }
  if(lntries1>10 || lntries2>10){
    dualfprintf(fail_file,"%d %ld :: %d %d %d :: lpflag1=%d lpflag2=%d errx1=%21.15g errx2=%21.15g lntries1=%d lntries2=%d\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag1,lpflag2,lerrx1,lerrx2,lntries1,lntries2);
  }

  
  // always do (use old utoprim)
  if(lpflag1==0) PLOOP(pl){
    pr[pl]=ptest1[pl];
    *lpflag=lpflag1;
  }
  else{
    PLOOP(pl) pr[pl]=ptest2[pl];
    *lpflag=lpflag2;
  }

  

  return(0);
}




// just picks the algorithm to invert
int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr)
{
  if(which==UTOPRIM5D1){      MYFUN(Utoprim(parameter,Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim", 1);}
  else if(which==UTOPRIMLDZ){ MYFUN(Utoprim_ldz(Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_ldz", 1);}
  else if(which==UTOPRIM2D){ MYFUN(Utoprim_2d(Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_2d", 1);}
  else if(which==UTOPRIM1D){ MYFUN(Utoprim_1d(Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_1d", 1);}
  else if(which==UTOPRIM1DOPT){ MYFUN(Utoprim_1d_opt(Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_1d_opt", 1);}
  else if(which==UTOPRIM1DFINAL){ MYFUN(Utoprim_1d_final(Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_1d_final", 1);}
  else if(which==UTOPRIM2DFINAL){ MYFUN(Utoprim_2d_final(Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_2d_final", 1);}
  else if(which==UTOPRIMJONNONRELCOMPAT){ MYFUN(Utoprim_jon_nonrelcompat_inputnorestmass(EOMTYPE,Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_2d_final_nonrelcompat_inputnorestmass", 1);}
  else if(which==UTOPRIM5D2){ MYFUN(Utoprim_5d2_final(Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_5d2_final", 1);}
  
  return(0);
}




// Utoprimgen try again.  Can be whatever sequence or number of inversions
int Utoprimgen_tryagain(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, int *lpflag, FTYPE *pr0, FTYPE *pr)
{
  int Utoprimgen_tryagain_substep(int which, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr0, FTYPE *pr);

  Utoprimgen_tryagain_substep(UTOPRIMJONNONRELCOMPAT, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM5D1, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);
  // LDZ as normally used generates nan's
  //  Utoprimgen_tryagain_substep(UTOPRIMLDZ, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM1DFINAL, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM2DFINAL, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM1DOPT, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM5D2, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM2D, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM1D, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr);

  // don't restore in the end, leave to whatever state inversion leaves it in.

  return(0);

}




int Utoprimgen_tryagain_substep(int which, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr0, FTYPE *pr)
{
  int k;

  // if really a bad failure that don't want / can't handle, then try again
  if(
     (*lpflag!=0)
     &&(!((*lpflag==UTOPRIMFAILUNEG)&&(STEPOVERNEGU)))
     &&(!((*lpflag==UTOPRIMFAILRHONEG)&&(STEPOVERNEGRHO)))
     &&(!((*lpflag==UTOPRIMFAILRHOUNEG)&&(STEPOVERNEGRHOU)))
     ){
    // restore
    PLOOP(k){
      Ugeomfree[k]=Ugeomfree0[k];
      pr[k]=pr0[k];
    }
    // try again
    MYFUN(Utoprimgen_pick(which,parameter,Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen_tryagain_substep()", "Utoprimgen_pick", 1);
  }

  return(0);
}




// there may be something wrong with this function -- didn't work in TIMEORDER==4, had to do standard method
// could have just been that I wasn't bounding after using this
int Utoprimloop(FTYPE unew[][N2M][N3M][NPR],FTYPE pf[][N2M][N3M][NPR])
{
  struct of_geom geom;
  int i,j,k;

  ZLOOP{
    get_geometry(i, j, k, CENT, &geom);
    // invert True U->p
    MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,unew[i][j][k], &geom, pf[i][j][k]),"step_ch.c:advance()", "Utoprimgen", 1);
  }
  return(0);
}


int primtoUloop(FTYPE pi[][N2M][N3M][NPR],FTYPE unew[][N2M][N3M][NPR])
{
  struct of_geom geom;
  struct of_state q;
  int i,j,k;

  ZLOOP{
    MYFUN(get_state(pi[i][j][k], &geom, &q),"step_ch.c:primtoUloop()", "get_state()", 1);
    get_geometry(i, j, k, CENT, &geom);
    // forward calculate U(p)
    MYFUN(primtoU(UEVOLVE,pi[i][j][k], &q, &geom, unew[i][j][k]),"step_ch.c:primtoUloop()", "primtoU()", 1);
  }
  return(0);
}



// filter out velocity along field line
void filterffde(int i, int j, int k, FTYPE *pr)
{
  int pl;
  struct of_geom geom;
  FTYPE prout[NPR],prin[NPR];
  FTYPE U[NPR];
  struct of_state q;
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)  ;


  // now filter velocity
  get_geometry(i,j,k,CENT,&geom);
  get_state(pr,&geom,&q);
  primtoU(UNOTHING,pr,&q,&geom,U);

  //  Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion
  Utoprimgen(EVOLVEUTOPRIM,UNOTHING,U,&geom,prout);

  PLOOP(pl) pr[pl]=prout[pl];
  // kill densities
  pr[RHO]=pr[UU]=0.0;



}
