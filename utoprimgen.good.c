


#include "decs.h"



int Utoprimgen(int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr)
{
  // debug
  int i, j, k;
  FTYPE Ugeomfree[NPR],Ugeomfree0[NPR];
  FTYPE pr0[NPR];
  FTYPE prother[NPR];
  int whichentropy;
  struct of_state q;
  FTYPE Uold[NPR],Unew[NPR];
  FTYPE fdiff[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  extern int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  extern int Utoprim_coldgrmhd(FTYPE *U, struct of_geom *geom, FTYPE *pr, int *positivityproblem);
  int positivityproblem;
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr);
  int Utoprimgen_compare(int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr);
  int Utoprimgen_tryagain(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, int *lpflag, FTYPE *pr0, FTYPE *pr);



  // Notice that reconstruct "standard" geometry-free conserved quantities by first modifying the geometry prefactor and THEN dealing with rest-mass or other subtractions and additions.
  // This is consistent with primtoflux() and how primtoflux() is used in source_conn().
  UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);

  // backup
  PLOOP(k){
    Uold[k]=Ugeomfree0[k]=Ugeomfree[k];
    pr0[k]=pr[k];
  }


  // check if inversion is chosen consistent with REMOVERESTMASSFROMUU
  // At this point, Ugeomfree[UU] in "standard form with rest-mass" for REMOVERESTMASSFROMUU=0,1.  If ==2, then no rest-mass
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
    // then change conserved quantity so can work for UTOPRIM2DFINAL
    Ugeomfree0[UU]+=Ugeomfree0[RHO]; // "subtract out rest-mass" so in correct form for NONRELCOMPAT
    Ugeomfree[UU]+=Ugeomfree[RHO]; // "subtract out rest-mass" so in correct form for NONRELCOMPAT
  }
#endif



  ////////////////////////
  //
  // DO INVERSION

  if(EOMTYPE==EOMGRMHD){

    if(DOENTROPY!=DOEVOLVEDIRECTENTROPY){
    
      if(UTOPRIMVERSION!=UTOPRIMCOMPARE) Utoprimgen_pick(UTOPRIMVERSION, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr);
      else Utoprimgen_compare(EVOLVENOENTROPY,Ugeomfree,ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr);

      // try other methods (assumes all methods can handle WHICHVEL, etc. used for primary model)
      // right now, all can handle WHICHVEL==VELREL4 and energy equation evolution and REMOVERESTMASSFROMUU=0,1
#if((WHICHVEL==VELREL4)&&(REMOVERESTMASSFROMUU<=1)&&(UTOPRIMTRYAGAIN))
      Utoprimgen_tryagain(EVOLVENOENTROPY, Ugeomfree0, Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr0, pr);
#endif
    }


    if(DOENTROPY!=DONOENTROPY){
      // also do separate inversion for entropy version of EOMs
      // only one inversion is setup to handle this
      PLOOP(k) prother[k]=pr[k]; // in any case this is a good guess
      if(DOENTROPY==DOEVOLVECOMPAREENTROPY) whichentropy=WHICHENTROPYEVOLVE;
      else if(DOENTROPY==DOEVOLVEDIRECTENTROPY) whichentropy=EVOLVEFULLENTROPY;

      MYFUN(Utoprim(whichentropy,Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], prother),"step_ch.c:Utoprimgen()", "Utoprim", 1);
      // result now contains internal energy (prother[UU,ENTROPY]) as derived from entropy evolution
      // notice that all other quantities are could also be different (if doentropy==evolvefullentropy), hence the prother variable for temporary storage.

      // at this point, ie version of entropy
      if(DOENTROPY==DOEVOLVEDIRECTENTROPY) PLOOP(k) pr[k]=prother[k]; // full evolution
      else if(DOENTROPY==DOEVOLVECOMPAREENTROPY){
	// just overwrite entropy primitive, leave rest same as from full energy equation
	pr[ENTROPY]=prother[ENTROPY];
	// now compare pr[UU] and pr[ENTROPY] with some kind of diagnostic?
	if(evolvetype==EVOLVEUTOPRIM){
	  // then during evolution and pr[UU]-pr[ENTROPY] is relevant to physical calculation
	  // store difference

	  // use difference to compute something
	}
      }

    }

  }
  else if(EOMTYPE==EOMFFDE){
    // GODMARK: inversions lead to different behavior!
    pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL]=0.0; // how is -100 getting in here?
    if(0){ // Jon's old inversion
      MYFUN(Utoprim_ffde(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_ffde", 1);
    }
    else{ // Jon's new inversion
      Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL], pr);
    }
    //Utoprim_ffde(Ugeomfree, ptrgeom, pr);
  }
  else if(EOMTYPE==EOMCOLDGRMHD){
    MYFUN(Utoprim_coldgrmhd(Ugeomfree, ptrgeom, pr,&positivityproblem),"step_ch.c:Utoprimgen()", "Utoprim_coldgrmhd", 1);
    //Utoprim_coldgrmhd(Ugeomfree, ptrgeom, pr,&positivityproblem);
  }

  // complain if pflag set
  // only output if failed
  if(pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL]&&(debugfail>=1)) 
		dualfprintf(fail_file, "Failed to find a prim. var. solution!! t=%21.15g steppart=%d nstep=%ld i=%d j=%d k=%d : fail=%d\n",t,steppart,realnstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL]);


  // test inversion for accuracy in conservative space
#if(0)
  MYFUN(get_state(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
  MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1);

  PLOOP(k){
    fdiff[k] = fabs(Unew[k]-Uold[k])/(fabs(Unew[k]+Uold[k])+1E-30);
    if(fdiff[k]>1E-10){
      if((k>=U1)&&(k<=B3)&&((fabs(Uold[k])>1E-20)||(fabs(Unew[k])>1E-20))){
	dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %g %g\n",k,fdiff[k],Uold[k],Unew[k]);
      }
    }
  }
#endif

     
  return(0);
}




//////////////////////////////
//
// COMPARISON of 2 (currently only 2) methods
int Utoprimgen_compare(int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr)
{
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, int *lpflag, FTYPE *pr);
  int k;
  FTYPE Uo1[NPR], Uo2[NPR];
  FTYPE ptest1[NPR],ptest2[NPR];
  FTYPE Uo[NPR], po[NPR];
  int test;


  // backup
  PLOOP(k){
    ptest1[k]=ptest2[k]=pr[k];
    Uo1[k]=Uo2[k]=Ugeomfree[k];
  }

  // currently comparing 5D1 and LDZ
  Utoprimgen_pick(UTOPRIM5D1, EVOLVENOENTROPY, Uo1, ptrgeom, lpflag, ptest1);
  Utoprimgen_pick(UTOPRIMLDZ, EVOLVENOENTROPY, Uo2, ptrgeom, lpflag, ptest2);

  PLOOP(k){
    if(ptest1[k]!=0.0) test=fabs(ptest1[k]-ptest2[k])/ptest1[k]>1E-11;
    else test=(ptest1[k]-ptest2[k])>1E-11;
    if(test){
      dualfprintf(fail_file,"utoprimdiff: %d %d %d  %21.15g   %21.15g   %21.15g %21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,k,ptest1[k]-ptest2[k],(ptest1[k]!=0.0) ? fabs(ptest1[k]-ptest2[k])/ptest1[k] : ptest1[k]-ptest2[k],ptest1[k],ptest2[k]);
    }
  }
  // always do (use old utoprim)
  PLOOP(k) pr[k]=ptest1[k];

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
  else if(which==UTOPRIMJONNONRELCOMPAT){ MYFUN(Utoprim_jon_nonrelcompat_inputnorestmass(Ugeomfree, ptrgeom, lpflag, pr),"step_ch.c:Utoprimgen()", "Utoprim_2d_final_nonrelcompat_inputnorestmass", 1);}
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
