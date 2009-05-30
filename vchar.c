
#include "decs.h"

/* 
 * calculate components of magnetosonic velocity 
 * corresponding to primitive variables p 
 *
 * cfg 7-10-01
 * 
 */


#define BOUNDARYSMOOTHER 0
// whether to use larger vmin/vmax near boundaries to help with large changes in primitive quantities that Riemann solver cannot handle
// this also can control courant condition to be consistent, unless set IGNORECOURANT 1
#define SMOOTHFACTOR 2.0
// factor by which to make vmin/vmax larger.
#define IGNORECOURANT 0
// whether to ignore vchar result in courant condition.  Otherwise dt may be controlled by boundary smoother.  Ignoring may lead to numerical instabilities.
// this leads to numerical instabilities!

// whether to check if phase velocity is >c and limit it to c.
// not quite right.  really phase or group velocity can have independent other components.  Other than k direction, can move with different velocity, so not a good comparison
#define CHECKSOL 0

// whether to use group velocity rather than phase velocity
// mutually exclusive with CHECKSOL
// actually doesn't work when va2=1.  ctsq->\infty, should avoid somehow
#define USEGROUPVEL 0

//#define GTOL 1.e-8
// Charles believed that near polar axis one should set vmax=vmin=0, but not right since makes flux completely non-diffusive (leads to numerical instabilities)
#define TRUEFAST 0
// TRUEFAST=0 : use overall larger phase velocity
// TRUEFAST=1 : use true phase velocity (has problems with catastrophic cancellation)

int vchar(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin,int *ignorecourant)
{
  FTYPE vgmax,vgmin,ctsq;
  FTYPE bsq, WW, EF, va2, cs2, cms2, rho, u, P;
  int simplefast(int dir, struct of_geom *geom,struct of_state *q, FTYPE cms2,FTYPE *vmin, FTYPE *vmax);
  int realfast(int dir,FTYPE EF,FTYPE cs2,FTYPE cms2,FTYPE va2,FTYPE *ucon,FTYPE *bcon,FTYPE (*gcon)[NDIM],FTYPE *vmin,FTYPE *vmax);
  int boundary_smoother(struct of_geom *ptrgeom, FTYPE *vmax, FTYPE *vmin, int *ignorecourant);
  extern int limitv3(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *v);
  void groupvel(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin, FTYPE ctsq, FTYPE *vgmax, FTYPE *vgmin);

  int pl;

  


  /* for regions along poles */
  /*
  if (geom->g < GTOL) {
    *vmax = 0.;
    *vmin = 0.;
    return (0);
  }
  */

  /* find fast magnetosonic speed */
  bsq = dot(q->bcon, q->bcov);
  rho = pr[RHO];
  u = pr[UU];
  if(u<0) u=0; // cold fix
  if(rho<0) rho=0; // neg. density fix

  P = pressure_rho0_u(rho,u);
  

  WW = rho + u + P;
  if(WW>0.0){
    cs2 =  cs2_compute(rho,u);
  }
  else WW=cs2=0.0;

  // make sure sound speed is well-defined (needed in case fundamental variables are temporarily unphysical)

  if(cs2>=1.0) cs2=1.0-1e-9; // some upper limit resolvable by machine precision
  else if(cs2<=0.0) cs2=0.0; // lower limit


  //  dualfprintf(fail_file,"rho=%g u=%g p=%g WW=%g cs2=%g\n",rho,u,P,WW,cs2);
    

  EF = bsq + WW;
  if(EF>0.0) va2 = bsq / EF;
  else EF=va2=0.0;

  cms2 = cs2 + va2 - cs2 * va2;	/* and there it is... */

  //atch SASMARK debug print
  //if( nstep == 1 && steppart == 0 && 95 == geom->i && 0 == geom->j && 0 == geom-> k ) {
  //  PLOOP(pl) dualfprintf( fail_file, "%21.15g\n", pr[pl] );
  //}

  /* check on it! */
  if(cms2<0.0){
    
    if (cms2/(cs2+va2) > -NUMEPSILON*100.0) cms2=0.0;
    else{
      if (fail(FAIL_COEFF_NEG) >= 1){
	dualfprintf(fail_file,"dir=%d :\n bsq=%21.15g\n rho=%21.15g\n u=%21.15g\n WW=%21.15g\n EF=%21.15g\n va2=%21.15g\n cs2=%21.15g\n cms2=%21.15g\n",dir,bsq,rho,u,WW,EF,va2,cs2,cms2);
	return (1);
      }
    }
  }

  if (cms2 > 1.) {
    if (cms2 < 1.0+NUMEPSILON*100.0) cms2=1.0;
    else if (rho<0.0) cms2=1.0;
    else{
      dualfprintf(fail_file,"dir=%d :\n bsq=%21.15g\n rho=%21.15g\n u=%21.15g\n WW=%21.15g\n EF=%21.15g\n va2=%21.15g\n cs2=%21.15g\n cms2=%21.15g\n",dir,bsq,rho,u,WW,EF,va2,cs2,cms2);
      if (fail(FAIL_COEFF_SUP) >= 1) return (1);
    }
  }
  
  if(TRUEFAST==1){
    if(realfast(dir,EF,cs2,cms2,va2,q->ucon,q->bcon,geom->gcon,vmin,vmax)>=1) return(1);
  }
  else if(TRUEFAST==0){
    if(simplefast(dir,geom, q, cms2,vmin,vmax)>=1) return(1);

    
    
  }
  
#if(0)
  if(realfast(dir,EF,cs2,cms2,va2,q->ucon,q->bcon,geom->gcon,vmin,vmax)>=1) return(1);
  fprintf(stderr,"%d %d\n",geom->i,geom->j);
  fprintf(stderr,"rf: vmax=%21.15g vmin=%21.15g\n",*vmax,*vmin);
  if(simplefast(dir,geom, q, cms2,vmin,vmax)>=1) return(1);
  fprintf(stderr,"sf: vmax=%21.15g vmin=%21.15g\n",*vmax,*vmin);
#endif

#if(BOUNDARYSMOOTHER)
  boundary_smoother(geom,vmax,vmin,ignorecourant);
#else
  *ignorecourant=0;
#endif

  //  if(cms2==1.0){// isfinite
  //  dualfprintf(fail_file,"1 cms2=%21.15g vmin=%21.15g vmax=%21.15g\n",cms2,*vmin,*vmax);
  // }
  //if(!isfinite(*vmin)){
  //  dualfprintf(fail_file,"2 cms2=%21.15g vmin=%21.15g vmax=%21.15g\n",cms2,*vmin,*vmax);    
  // }
  //if(!isfinite(*vmax)){
  //  dualfprintf(fail_file,"3 cms2=%21.15g vmin=%21.15g vmax=%21.15g\n",cms2,*vmin,*vmax);    
  // }


#if(USEGROUPVEL)
  // ct^2 = (va^2 + cs^2(1-va^2) ) / ( 1 - (va^2+cs^2(1-va^2)))
  // only really correct with simplified dispersion relation
  ctsq=(va2+cs2*(1.0-va2))/(1.0-(va2+cs2*(1.0-va2)));
  groupvel(pr,q,dir,geom,vmax,vmin,ctsq,&vgmax,&vgmin);
  *vmax=vgmax;
  *vmin=vgmin;
#endif


  // check to make sure phase velocity is not greater than speed of light
#if(CHECKSOL)
  limitv3(pr,q,dir,geom,vmax);
  limitv3(pr,q,dir,geom,vmin);
#endif


  return (0);
}


// method didn't work
int boundary_smoother(struct of_geom *geom, FTYPE *vmax, FTYPE *vmin, int *ignorecourant)
{


  if(
     // x3-surfaces
     // 1-zones
     (((geom->p==CENT)||(geom->p==FACE1)||(geom->p==FACE2)||(geom->p==CORN3))&&((startpos[3]+geom->k<=0)||(startpos[3]+geom->k>=totalsize[3]-1))) ||
     // 2-zones
     (((geom->p==FACE3)||(geom->p==CORN1)||(geom->p==CORN2))&&((startpos[3]+geom->k<=1)||(startpos[3]+geom->k>=totalsize[3]-1))) ||

     // x2-surfaces
     // 1-zones
     (((geom->p==CENT)||(geom->p==FACE1)||(geom->p==FACE3)||(geom->p==CORN2))&&((startpos[2]+geom->j<=0)||(startpos[2]+geom->j>=totalsize[2]-1))) ||
     // 2-zones, but for POLARAXIS, absolute face2 on boundary always leads to 0 flux, so only difference can be next zone value.
     (((geom->p==FACE2)||(geom->p==CORN3)||(geom->p==CORN1))&&((startpos[2]+geom->j<=1)||(startpos[2]+geom->j>=totalsize[2]-1))) ||

     // x1-surfaces
     // 1-zones
     (((geom->p==CENT)||(geom->p==FACE2)||(geom->p==FACE3)||(geom->p==CORN1))&&((startpos[1]+geom->i<=0)||(startpos[1]+geom->i>=totalsize[1]-1))) ||
     // 2-zones
     (((geom->p==FACE1)||(geom->p==CORN3)||(geom->p==CORN2))&&((startpos[1]+geom->i<=1)||(startpos[1]+geom->i>=totalsize[1]-1)))
     ){
    *ignorecourant=IGNORECOURANT;
    *vmax *=SMOOTHFACTOR;
    *vmin *=SMOOTHFACTOR;
  }
  else{
    // don't normally ignore vmin/vmax in Courant condition
    *ignorecourant=0;
  }

  return(0);
}


#define USESASHAREWRITE 1

int simplefast(int dir,struct of_geom *geom, struct of_state *q, FTYPE cms2,FTYPE *vmin, FTYPE *vmax)
{
  FTYPE discr;
  FTYPE Acov[NDIM], Bcov[NDIM], Acon[NDIM], Bcon[NDIM];
  FTYPE Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;
  FTYPE vm,vp;
  int j;

  DLOOPA(j) Acov[j] = 0.;
  Acov[dir] = 1.;
  raise_vec(Acov, geom, Acon);

  DLOOPA(j) Bcov[j] = 0.;
  Bcov[TT] = 1.;
  raise_vec(Bcov, geom, Bcon);

  /* now require that speed of wave measured by observer q->ucon is
     cms2 */
  Asq = dot(Acon, Acov);
  Bsq = dot(Bcon, Bcov);
  Au = dot(Acov, q->ucon);
  Bu = dot(Bcov, q->ucon);
  AB = dot(Acon, Bcov);
  Au2 = Au * Au;
  Bu2 = Bu * Bu;
  AuBu = Au * Bu;

  //OBSOLETE -- START
  // JCM: grouped cms2 with 1.0 so good in ultrarelativistic regime.
  // was getting A=B=0 and giving vm=vp=inf, when really wasn't A=B=0!  That is, code was failing when cms=1 and A=B=0
  A = Bu2 * (1.0 - cms2) - Bsq * cms2;
  //  A = Bu2 - (Bsq + Bu2) * cms2;
  B = 2. * (AuBu * (1.0-cms2)  - AB*cms2);
  //B = 2. * (AuBu - (AB + AuBu) * cms2);
  C = Au2*(1.0-cms2) - Asq * cms2;
  //  C = Au2 - (Asq + Au2) * cms2;

#if(USESASHAREWRITE==0)
  // order unity normalized already
  discr = B * B - 4 * A * C;
  //OBSOLETE -- END
#else
  //REWRITTEN without catastrophic cancellation
  discr = 4.0 * cms2 *
    ( 
     (AB * AB - Asq * Bsq) * cms2
     + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (cms2 - 1.0)
     );
#endif

  // DEBUG
  //  if(A<=0.0){// isfinite
  // dualfprintf(fail_file,"A=%21.15g B=%21.15g C=%21.15g discr=%21.15g\n",A,B,C,discr);
  // }



  if((discr<0.0)&&(discr> -1E-10)) discr=0.0;  //atch: negative but not too much -- allow fractional error of 1e-10
  else if(discr<0.0){  //if discriminant is too negative
    dualfprintf(fail_file,"simplefast discr=%21.15g, nstep = %ld, steppart == %d, i = %d\n",
	discr, nstep, steppart, geom->i );

    dualfprintf(fail_file, "\n\t %21.15g %21.15g %21.15g %21.15g %21.15g\n", A, B, C, discr, cms2);
    dualfprintf(fail_file, "\n\t q->ucon: %21.15g %21.15g %21.15g %21.15g\n", q->ucon[0],
	    q->ucon[1], q->ucon[2], q->ucon[3]);
    dualfprintf(fail_file, "\n\t q->bcon: %21.15g %21.15g %21.15g %21.15g\n", q->bcon[0],
	    q->bcon[1], q->bcon[2], q->bcon[3]);
    dualfprintf(fail_file, "\n\t Acon: %21.15g %21.15g %21.15g %21.15g\n", Acon[0], Acon[1],
	    Acon[2], Acon[3]);
    dualfprintf(fail_file, "\n\t Bcon: %21.15g %21.15g %21.15g %21.15g\n", Bcon[0], Bcon[1],
	    Bcon[2], Bcon[3]);
    if (fail(FAIL_VCHAR_DISCR) >= 1)
      return (1);
    discr = 0.;
  }

  discr = sqrt(discr);

  // order not obvious
  vm = -(-B + discr) / (2. * A);
  vp = -(-B - discr) / (2. * A);


  // isfinite check -- non-debug check
  if( (!isfinite(vm)) || (!isfinite(vp)) ){
    dualfprintf(fail_file,"vm=%21.15g vp=%21.15g discr=%21.15g A=%21.15g B=%21.15g C=%21.15g\n",vm,vp,discr,A,B,C);
    myexit(10001); // hard failure so look back and see what's wrong with code
  }

  // order
  if (vp < vm) {
    *vmax = vm;
    *vmin = vp;
  }
  else{
    *vmax = vp;
    *vmin = vm;
  }
  return(0);
}


// See grmhd-sonicpoint.nb and grmhd-quartic?.nb
// related: see grmhdwaves2.nb
int realfast(int dir,FTYPE EF,FTYPE cs2,FTYPE cms2,FTYPE va2,FTYPE *ucon,FTYPE *bcon,FTYPE (*gcon)[NDIM],FTYPE *vmin,FTYPE *vmax)
{
  FTYPE va02,vax2,va0x2,uux;
  FTYPE oEF;
  FTYPE gn300,gn3xx,gn30x;
  FTYPE uu0;
  FTYPE AA,BB,CC,DD,EE;
  FTYPE uuxsq,uux3,uux4,uu0sq,uu03,uu04;
  FTYPE quarticsol(FTYPE sign1, FTYPE sign2,FTYPE AAA,FTYPE BBB,FTYPE CCC,FTYPE DDD,FTYPE EEE);
  int i,j;

#if(0)
  // DEBUG
  if((icurr==32)&&(jcurr==40)&&(kcurr==32)){
    fprintf(stderr,"%d %21.15g %21.15g %21.15g %21.15g\n",dir,EF,cs2,cms2,va2);
    for(i=0;i<NDIM;i++) fprintf(stderr,"%21.15g ",ucon[i]); fprintf(stderr,"\n");
    for(i=0;i<NDIM;i++) fprintf(stderr,"%21.15g ",bcon[i]); fprintf(stderr,"\n");
    for(i=0;i<NDIM;i++){ for(j=0;j<NDIM;j++) fprintf(stderr,"%21.15g ",gcon[i][j]); fprintf(stderr,"\n");}
  }
#endif

  oEF=1.0/EF;
  gn300=gcon[0][0];
  va02=bcon[0]*bcon[0]*oEF;
  uu0=ucon[0];

  vax2=bcon[dir]*bcon[dir]*oEF;
  va0x2=bcon[0]*bcon[dir]*oEF;
  uux=ucon[dir];
  gn3xx=gcon[dir][dir];
  gn30x=gcon[0][dir];

  uuxsq=uux*uux;
  uux3=uuxsq*uux;
  uux4=uuxsq*uuxsq;
  uu0sq=uu0*uu0;
  uu03=uu0sq*uu0;
  uu04=uu0sq*uu0sq;

  // quartic coefficients
  EE=uux4 - cms2*uuxsq*(gn3xx + uuxsq) + cs2*gn3xx*vax2;
  DD=2.*(-2.*uu0*uux3 + cms2*uux*(gn3xx*uu0 + uux*(gn30x + 2.*uu0*uux)) - cs2*gn3xx*va0x2 - cs2*gn30x*vax2);
  CC=6.*uu0sq*uuxsq - cms2*(gn3xx*uu0sq + uux*(4.*gn30x*uu0 + gn300*uux + 6.*uu0sq*uux)) + cs2*(4.*gn30x*va0x2 + gn3xx*va02 + gn300*vax2);
  BB=2.*(-2.*uu03*uux + cms2*uu0*(gn300*uux + uu0*(gn30x + 2.*uu0*uux)) - cs2*gn300*va0x2 - cs2*gn30x*va02);
  AA=uu04 - cms2*uu0sq*(gn300 + uu0sq) + cs2*gn300*va02;

#if(0)
  if((icurr==32)&&(jcurr==40)&&(kcurr==32)){
    fprintf(stderr,"%21.15g %21.15g %21.15g %21.15g %21.15g\n",AA,BB,CC,DD,EE);
  }
#endif

#if(0)
  // test
  AA=4.36134649461842;
  BB=0.953704334932897;
  CC=0.00501787191044863;
  DD=-0.0036931214344023;
  EE=3.97410207658954e-05;
#endif

  *vmax=quarticsol(1.0,1.0,AA,BB,CC,DD,EE);
  *vmin=quarticsol(-1.0,-1.0,AA,BB,CC,DD,EE);

#if(0)
  if((icurr==32)&&(jcurr==40)&&(kcurr==32)){
    fprintf(stderr,"insi vmax=%21.15g vmin=%21.15g\n",*vmax,*vmin);
    *vmax=quarticsol(-1.0,1.0,AA,BB,CC,DD,EE);
    *vmin=quarticsol(1.0,-1.0,AA,BB,CC,DD,EE);
    fprintf(stderr,"slow vmax=%21.15g vmin=%21.15g\n",*vmax,*vmin);
    exit(0);
  }
#endif
  
  return(0);

}

// one should reorder phase speeds for sign1!=sign2 since results in ambiguous order.  Fast terms are unambiguously ordered always
FTYPE quarticsol(FTYPE sign1, FTYPE sign2,FTYPE AAA,FTYPE BBB,FTYPE CCC,FTYPE DDD,FTYPE EEE)
{
  FTYPE unit1,unit2,unit3,newradius,newphi,unit4;
  FTYPE part1,part2,part25,part3;
  FTYPE AAAsq,AAAcu,BBBsq,BBBcu,CCCsq,CCCcu,DDDsq,unit1sq,unit2sq,unit2cu;
  FTYPE disc;

  AAAsq=AAA*AAA;
  AAAcu=AAAsq*AAA;
  BBBsq=BBB*BBB;
  BBBcu=BBBsq*BBB;
  CCCsq=CCC*CCC;
  CCCcu=CCCsq*CCC;
  DDDsq=DDD*DDD;

  unit1=2.*CCCcu - 9.*BBB*CCC*DDD + 27.*AAA*DDDsq + 27.*BBBsq*EEE - 72.*AAA*CCC*EEE;
  // unit2 is always positive or 0
  unit2=CCCsq - 3.*BBB*DDD + 12.*AAA*EEE;

  unit2sq=unit2*unit2;
  unit2cu=unit2sq*unit2;
  unit1sq=unit1*unit1;

  // unit3 is always negative or 0, but may be positive due to machine precision error
  unit3=-4.*unit2cu + unit1sq;
  if((unit3>0.0)&&(unit3<1E-3)) unit3=0.0;
  else if(unit3>1E-3){
    dualfprintf(fail_file,"unit3=%21.15g\n",unit3);
    exit(0);
  }

  newradius=sqrt(unit1sq+fabs(-unit3));
  // ATAN2(y,x) where tan(phi)=y/x
  newphi=atan2(sqrt(-unit3),unit1);
  // unit4 is always positive
  unit4=(1./(3.*AAA))*pow(2.,(2./3.))*pow(newradius,(1./3.))*cos(newphi/3.);

  // part1 can be anything
  part1 = -BBB/(4.*AAA);

  // part2 >0 always
  disc=BBBsq/(4.*AAAsq) - 2./3.*CCC/AAA + unit4;

  // disc already normalized
  if(disc<0.0){
    if(disc>-NUMEPSILON*100.0){
      part2=0.0;
      part3=0.0;
    }
    else{
      dualfprintf(fail_file,"part2 disc=%21.15g\n",disc);
      myexit(1);
    }
  }
  else{
    part2 = sqrt(disc);

    // part25 can be anything
    part25 = (-BBBcu/AAAcu + 4.*BBB*CCC/AAAsq - 8.*DDD/AAA)/(4.*part2);

    // part3>0 always
    disc=BBBsq/(2.*AAAsq) - 4./3.*CCC/AAA - unit4 +sign1*part25;

    //disc already normalized
    if(disc<0.0){
      if(disc>-NUMEPSILON*100.0) disc=0.0;
      else{
      dualfprintf(fail_file,"part3 disc=%21.15g\n",disc);
      myexit(1);
      }
    }
    part3 = sqrt(disc);
  }


  return(part1 +sign1* 0.5*part2 +sign2* 0.5*part3);
}

// group velocities along a particular direction

// assumes vmax and vmin are filled with left and right going phase velocities
// w^2 = ct^2 k^2 (little k)

// e.g. for simplified dispersion relation:
// ct^2 = (va^2 + cs^2(1-va^2) ) / ( 1 - (va^2+cs^2(1-va^2)))

// the velocity of the wave packet can have completely different velocity components than fluid velocity -- even for k along one direction.

void groupvel(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin, FTYPE ctsq, FTYPE *vgmax, FTYPE *vgmin)
{
  FTYPE uu[NDIM],vg[NDIM];
  FTYPE kd[NDIM],ku[NDIM],omega;
  FTYPE uu0,vu[NDIM];
  FTYPE vp,vdiff;
  FTYPE git,gtt,gii;
  FTYPE ctrat;
  int j;
  FTYPE vgd[NDIM];
  FTYPE isone1,isone2;



#if(0)
  uu0=q->ucon[TT];
  SLOOPA(j) vu[j]=q->ucon[j]/uu0;

  ctrat=ctsq/(uu0*uu0);
  
  git=geom->gcon[dir][TT];
  gtt=geom->gcon[TT][TT]; 
  gii=geom->gcon[dir][dir];


  // omega=k.u
  //#define numerator (vu[dir]*vdiff-ctrat*(git*uu0*vdiff+gii)) // vg^i/(u^t)^2
  //#define denominator (vdiff-ctrat*(uu0*vdiff*gtt+git)) // vg^t/(u^t)^2

  // omega=-k.u
  //#define numerator (vu[dir]*vdiff-ctrat*(git*uu0*vdiff-gii)) // vg^i/(u^t)^2
  //#define denominator (vdiff-ctrat*(uu0*vdiff*gtt-git)) // vg^t/(u^t)^2

  // correct omega=-k.u and k=ke(-vp,e)
  //#define numerator (vu[dir]*vdiff-ctrat*(git*vp-gii)) // vg^i/(u^t)^2
  //#define denominator (vdiff-ctrat*(gtt*vp-git)) // vg^t/(u^t)^2

  vp=*vmax;
  vdiff=(vu[dir]-vp);
  *vgmax=numerator/denominator;

  vp=*vmin;
  vdiff=(vu[dir]-vp);
  *vgmin=numerator/denominator;
#endif

  DLOOPA(j) uu[j]=q->ucon[j];

  DLOOPA(j) kd[j]=0.0;
  kd[TT]=-*vmax;
  kd[dir]=1.0;
  raise_vec(kd,geom,ku);
  omega=0.0; DLOOPA(j) omega+=-kd[j]*uu[j];
  DLOOPA(j) vg[j]=(-omega*uu[j]-ctsq*ku[j]);
  *vgmax=vg[dir]/vg[TT];

#if(0)
  // check if really 4-velocity
  lower_vec(vg,geom,vgd);
  isone1=0.0;
  DLOOPA(j) isone1+=vg[j]*vgd[j];
  isone1/=(omega*omega*(1.0+ctsq));

  // check 3-velocity
  dualfprintf(fail_file,"vmax: vfx=%21.15g vgx=%21.15g : vfy=%21.15g vgy=%21.15g : vfz=%21.15g vgz=%21.15g\n",q->ucon[1]/q->ucon[0],vg[1]/vg[0],q->ucon[2]/q->ucon[0],vg[2]/vg[0],q->ucon[3]/q->ucon[0],vg[3]/vg[0]);
#endif

  DLOOPA(j) kd[j]=0.0;
  kd[TT]=-*vmin;
  kd[dir]=1.0;
  raise_vec(kd,geom,ku);
  omega=0.0; DLOOPA(j) omega+=-kd[j]*uu[j];
  DLOOPA(j) vg[j]=(-omega*uu[j]-ctsq*ku[j]);
  *vgmin=vg[dir]/vg[TT];

#if(0)
  // check if really 4-velocity
  lower_vec(vg,geom,vgd);
  isone2=0.0;
  DLOOPA(j) isone2+=vg[j]*vgd[j];
  isone2/=(omega*omega*(1.0+ctsq));

  // check if really 4-velocity
  dualfprintf(fail_file,"%d %d : isone: %21.15g %21.15g\n",geom->i,geom->j,isone1,isone2);
  // checks out.


  // check 3-velocity
  dualfprintf(fail_file,"vmin: vfx=%21.15g vgx=%21.15g : vfy=%21.15g vgy=%21.15g : vfz=%21.15g vgz=%21.15g\n",q->ucon[1]/q->ucon[0],vg[1]/vg[0],q->ucon[2]/q->ucon[0],vg[2]/vg[0],q->ucon[3]/q->ucon[0],vg[3]/vg[0]);

#endif

}
