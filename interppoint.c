
#include "decs.h"




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// POINTS METHODS
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// whether to send all pl's for access by interpolator (such as for steepening and flattening)



// interpolation with loop over POINTS
void slope_lim_pointtype(int numprims, int pl, int dir, int idel, int jdel, int kdel, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP])
{
  int i,j,k;
  void slope_lim_point(int reallim, int pl, int startorderi, int endorderi, FTYPE*y, FTYPE *dq,FTYPE *left,FTYPE *right);
  void slope_lim_point_allpl(int dir, int reallim, int startorderi, int endorderi, FTYPE **y, FTYPE *dq,FTYPE *left,FTYPE *right);
  extern int choose_limiter(int i, int j, int k, int pl);
  void set_interppoint_loop(int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
  void set_interppoint_loop_expanded(int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
  int l;
  FTYPE interplist[MAXSPACEORDER];
  FTYPE interplistpl[NPR2INTERP][MAXSPACEORDER];
  FTYPE *yin;
  int ijkshift;
  FTYPE *y;
  FTYPE *ypl[NPR2INTERP];
  int startorderi,endorderi;
  int reallim;
  int is,ie,js,je,ks,ke,di,dj,dk;
  int plpl;



  // make sure dq's exist if going to use them
  if( (DODQMEMORY==0)&&(LIMADJUST==LIMITERFIXED) ){
    dualfprintf(fail_file,"Must have dq's when using MC or lower second order methods\n");
    myexit(17);
  }


  if(DOENOFLUX!=ENOFINITEVOLUME){
    set_interppoint_loop(dir, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);
  }
  else{
    set_interppoint_loop_expanded(dir, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);
  }


  ZSLOOP( is, ie, js, je, ks, ke ) {

    reallim=choose_limiter(i,j,k,pl);
    //    if((abs(i-N1/2)<3)) reallim=DONOR;
    //    reallim=DONOR;


#if(0&&BOUNDARYINTERPADJUST) // GODMARK: not setup to do anything yet

    // setup shift
    ijkshift = i*idel+j*jdel+k*kdel;

    // get starting point for stencil, assumes quasi-symmetric (including even/odd size of stencil)
    startorderi= ijkshift - interporder[reallim]/2;

    // check if interpolation reaches beyond valid range on "left"
    if(startoderi < ijkminmax[NONENOINTERPTYPE][dir][0]) startorderi = ijkminmax[NONENOINTERPTYPE][dir][0];

    // let endorder go as far as possible
    endorderi = startorderi + interporder[reallim] - 1;

    // check if interpolation reaches beyond valid range on "right"
    if(endorderi > ijkminmax[NONENOINTERPTYPE][dir][1]) endorderi = ijkminmax[NONENOINTERPTYPE][dir][1];

    // shift so that i,j,k  = 0
    startorderi -= ijkshift;
    endorderi   -= ijkshift;

#else

    // get starting point for stencil, assumed quasi-symmetric (including even/odd size of stencil)
    startorderi = - interporder[reallim]/2;
    endorderi   = - startorderi;

#endif

    ////////////////////
    //
    // For limiters that use specific information about the equations and need all variables, use slope_lim_point_allpl()
    //
    ///////////////////
    if(reallim==PARAFLAT){

      startorderi = - (interporder[reallim])/2;
      endorderi   = - startorderi;
      
      PLOOPINTERP(plpl){
	ypl[plpl] = interplistpl[plpl] - startorderi;
	
	// get interpolation points, where y[0] is point of interest for which interpolation is found.
	for(l=startorderi;l<=endorderi;l++){
	  ypl[plpl][l]=p2interp[i + l*idel][j + l*jdel][k + l*kdel][plpl];
	}
      }
      slope_lim_point_allpl(dir,reallim,startorderi,endorderi,ypl,
			    &dq[i][j][k][0],&pleft[i][j][k][0],&pright[i][j][k][0]
			    );
      
    }
    else{
      ////////////////////
      //
      // For limiters that are general, use slope_lim_point()
      //
      ///////////////////

      // shift for easy use and clarity of loop and easy of use by slope_lim_point()
      y = interplist - startorderi;
      
      // get interpolation points, where y[0] is point of interest for which interpolation is found.
      for(l=startorderi;l<=endorderi;l++){
	y[l]=p2interp[i + l*idel][j + l*jdel][k + l*kdel][pl];
      }
      
      slope_lim_point(reallim,pl,startorderi,endorderi,y,
		      &dq[i][j][k][pl],&pleft[i][j][k][pl],&pright[i][j][k][pl]
		      );
    }

  }


  if(reallim==PARAFLAT){ // assumes last chosen reallim is constant // SUPERGODMARK
    pl=NPR2INTERP; // finish the loop outside this function
    if(LIMADJUST!=LIMITERFIXED){
      dualfprintf(fail_file,"Can't use PARAFLAT with LIMADJUST!=LIMITERFIXED\n");
      myexit(1);
    }
  }

}





// Setup loop over region of points
// very similar to set_interp_loop() in interpline.c
void set_interppoint_loop(int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{


  if(dir==1){

    *is=-1;
    *ie=N1;

    *js=INFULL2; // 0
    *je=OUTFULL2; // N2-1;

    *ks=INFULL3; //0;
    *ke=OUTFULL3; // N3-1;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==2){
    *is=INFULL1; //0;
    *ie=OUTFULL1; //N1-1;

    *js=-1;
    *je=N2;

    *ks=INFULL3; // 0;
    *ke=OUTFULL3; // N3-1;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==3){
    *is=INFULL1; // 0;
    *ie=OUTFULL1; // N1-1;

    *js=INFULL2; // 0;
    *je=OUTFULL2; // N2-1;

    *ks=-1;
    *ke=N3;

    *di=1;
    *dj=1;
    *dk=1;

  }
}





// Setup loop over region of points for finite volume method (or any method that uses extended ghost+active grid)
void set_interppoint_loop_expanded(int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{

  if(dir==1){
    *is=Uconsloop[FIS]-1;
    *ie=Uconsloop[FIE]+1;

    *js=INFULL2;
    *je=OUTFULL2;
    
    *ks=INFULL3;
    *ke=OUTFULL3;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==2){

    *is=INFULL1;
    *ie=OUTFULL1;
    
    *js=Uconsloop[FJS]-1;
    *je=Uconsloop[FJE]+1;
    
    *ks=INFULL3;
    *ke=OUTFULL3;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==3){

    *is=INFULL1;
    *ie=OUTFULL1;
    
    *js=INFULL2;
    *je=OUTFULL2;
    
    *ks=Uconsloop[FKS]-1;
    *ke=Uconsloop[FKE]+1;

    *di=1;
    *dj=1;
    *dk=1;

  }

}












// interpolate from a center point to a left/right interface
void slope_lim_point(int reallim, int pl, int startorderi, int endorderi, FTYPE*y, FTYPE *dq,FTYPE *left,FTYPE *right)
{
  void para(FTYPE *y, FTYPE *lout, FTYPE *rout);
  void para2(FTYPE *y, FTYPE *lout, FTYPE *rout);
  void para3(FTYPE *y, FTYPE *lout, FTYPE *rout);
  void para4(int pl, FTYPE *y, FTYPE *lout, FTYPE *rout);
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  void csslope(int pl, FTYPE *y, FTYPE *dq);
  FTYPE mydq = 0.0; //to avoid copying of nan's


  // parabolic reconstruction
  if (reallim == PARA){
    // para
#if(WHICHPARA==PARA1)
    para(y,left,right);
#elif(WHICHPARA==PARA2)
    para2(y,left,right);
#elif(WHICHPARA==PARA3)
    para3(y,left,right);
#elif(WHICHPARA==PARA4)
    para4(pl,y,left,right);
#endif
  }
  else if(reallim == CSSLOPE){
    csslope(pl, y, &mydq);
  }
  else{
    slope_lim_3points(reallim, y[-1], y[0], y[1], &mydq);
  }


  // need to set left/right when not PARA but adjusting limiter for each point so consistent with fluxcalc()
  if(LIMADJUST!=LIMITERFIXED){
    *left =y[0] - 0.5* mydq;
    *right=y[0] + 0.5* mydq;
  }
  else *dq=mydq; // only set to dq if necessary since DODQMEMORY may be 0

}

// interpolate from a center point to a left/right interface
void slope_lim_point_allpl(int dir,int reallim, int startorderi, int endorderi, FTYPE **y, FTYPE *dq,FTYPE *left,FTYPE *right)
{
  void parapl(int dir, FTYPE **y, FTYPE *lout, FTYPE *rout);

  parapl(dir,y,left,right);

}




void slope_lim_3points_old(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq)
{
  FTYPE Dqm, Dqp, Dqc, s;

  if (reallim == MC) {
    Dqm = 2.0 * (yc - yl);
    Dqp = 2.0 * (yr - yc);
    Dqc = 0.5 * (yr - yl);
    s = Dqm * Dqp;
    if (s <= 0.)  *dq= 0.;
    else{
      if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
	*dq= (Dqm);
      else if (fabs(Dqp) < fabs(Dqc))
	*dq= (Dqp);
      else
	*dq= (Dqc);
    }
  }
  /* van leer slope limiter */
  else if (reallim == VANL) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    s = Dqm * Dqp;
    if (s <= 0.)
      *dq= 0.;
    else
      *dq= (2.0 * s / (Dqm + Dqp));
  }
  /* minmod slope limiter (crude but robust) */
  else if (reallim == MINM) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    s = Dqm * Dqp;
    if (s <= 0.) *dq= 0.;
    else{
      if (fabs(Dqm) < fabs(Dqp)) *dq= Dqm;
      else *dq= Dqp;
    }
  }
  else if (reallim == NLIM) {
    Dqc = 0.5 * (yr - yl);
    *dq= (Dqc);
  }
  else if (reallim == DONOR) {
    *dq=(0.0);
  }
  else {
    dualfprintf(fail_file, "unknown slope limiter: %d\n",reallim);
    myexit(10);
  }



}


void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq)
{
  FTYPE Dqm, Dqp, Dqc, s;
  FTYPE theta;

  if (reallim == MC) { // monotonized central (Woodward) (Barth-Jespersen)
    Dqm = 2.0 * (yc - yl);
    Dqp = 2.0 * (yr - yc);
    Dqc = 0.5 * (yr - yl);
    *dq=MINMOD(Dqc,MINMOD(Dqm,Dqp));
  }
  /* van leer slope limiter */
  else if (reallim == VANL) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    s = Dqm * Dqp;
    if (s <= 0.)
      *dq= 0.;
    else
      *dq= (2.0 * s / (Dqm + Dqp));
  }
  /* minmod slope limiter (crude but robust) */
  /*
  else if (reallim == MINM) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    *dq = MINMOD(Dqm,Dqp);
  }
  */
  else if (reallim == MINM) {
    theta = 1.0;
    Dqm = theta * (yc - yl);
    Dqp = theta * (yr - yc);
    Dqc = 0.5 * (yr - yl);
    //*dq = MINMOD3(Dqm,Dqc,Dqp);    
    *dq = MINMOD(MINMOD(Dqm,Dqc),Dqp);
    }
  else if (reallim == NLIMCENT) { // Centered slope (Fromm)
    *dq = 0.5 * (yr - yl);
  }
  else if (reallim == NLIMUP) { // Upwind slope (Beam-Warming)
    *dq = (yc - yl);
  }
  else if (reallim == NLIMDOWN) { // Downwind slope (Lax-Wendroff)
    *dq = (yr - yc);
  }
  else if (reallim == DONOR) { // no slope
    *dq=(0.0);
  }
  else {
    dualfprintf(fail_file, "unknown slope limiter: %d\n",reallim);
    myexit(10);
  }



}


// see Mignone & Bodo (2005) astro-ph/0506414 equations 27-31
// see Colella (1985), Saltzman (1994)
// second order with 4th order steepeners
void csslope(int pl, FTYPE *y, FTYPE *dq)
{
  FTYPE s, sm, sp;
  FTYPE Dql, Dqlm, Dqlp;
  FTYPE Dqp, Dqpp;
  FTYPE Dqm, Dqmm;
  FTYPE Dqc, Dqcm, Dqcp;
  FTYPE Dqbp,Dqbm;
  FTYPE alpha;


  Dqp=(y[1]-y[0]);
  Dqpp=(y[2]-y[1]);

  Dqm=(y[0]-y[-1]);
  Dqmm=(y[-1]-y[-2]);

  Dqc=0.5*(y[1]-y[-1]);
  Dqcm=0.5*(y[0]-y[-2]);
  Dqcp=0.5*(y[2]-y[0]);

  // Dqmm  Dqm Dqp  Dqpp
  //    Dqcm Dqc Dqcp
  //     sm   q   sp
  //    Dqlm Dql Dqlp
  //    Dqbm     Dqbp

  s=0.5*(sign(Dqp)+sign(Dqm));
  sp=0.5*(sign(Dqpp)+sign(Dqp));
  sm=0.5*(sign(Dqm)+sign(Dqmm));

  // MB05 use alpha=2 for 1-D and alpha=2,1.25,1 for rho,v,p (respectively) for 2D.
  // alpha=[1-2].  alpha=2 more compressive, alpha=1 less compressive.
  if(pl==RHO) alpha=2.0;
  else if((pl>=U1)&&(pl<=U3)) alpha=1.25;
  else if(pl==UU) alpha=1.0;
  else if((pl>=B1)&&(pl<=B3)) alpha=1.25;
  else alpha=2.0;

  Dql=alpha*min(fabs(Dqp),fabs(Dqm));
  Dqlp=alpha*min(fabs(Dqpp),fabs(Dqp));
  Dqlm=alpha*min(fabs(Dqm),fabs(Dqmm));

  Dqbp=sp*min(Dqlp,fabs(Dqcp));
  Dqbm=sm*min(Dqlm,fabs(Dqcm));

  *dq=s*min(fabs(FOURTHIRD*Dqc-SIXTH*(Dqbp+Dqbm)),Dql);


}



/*
 * parabolic interpolation subroutin  
 * ref. Colella && Woodward's paper
 * Colella, P., & Woodward, P. R. 1984, J. Comput. Phys., 54, 174-201
 *
 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 * 
 */

// lout/rout is left and right sides of cell
// note how used in step_ch.c to get appropriate interface value

// given by Xiaoyue Guan to Scott Noble on Nov 9, 2004, given to me Jan 7, 2005
void para(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(aDqc,min(aDqm,aDqp))*sign(Dqc);
  }

  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;

  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }

  if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


  *lout=l;   //a_L,j
  *rout=r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}



// given by Xiaoyue Guan on Jan 9, 2005
void para2(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == PMC)
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MC)
    dq[mm] =Dqc;
#endif
  }
  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;

  /*
    l=max(min(y[0],y[-1]),l);
    l=min(max(y[0],y[-1]),l);
    r=max(min(y[0],y[1]),r);
    r=min(max(y[0],y[1]),r);
  */

  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }

  else if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


  *lout=l;   //a_L,j
  *rout=r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}




// 3rd para from Xiaoyue that she bundled with a new TVD-optimal RK3
// given on 02/17/2005
void para3(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);

    if (s <=0.) dq[mm]=0.;
    else dq[mm] = -aDqvanl*sign(Dqc);
    //else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MC)
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=-min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MINM)
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = -aDqm*sign(Dqc);
    else dq[mm]=-aDqp*sign(Dqc);
#elif(PARA2LIM == NLIM) //w/o slope limiter
    //if(s<=0.) dq[mm] = 0.; // DONOR
    dq[mm] = Dqc;
#endif
  }

  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;


  l=max(min(y[0],y[-1]),l);
  l=min(max(y[0],y[-1]),l);
  r=max(min(y[0],y[1]),r);
  r=min(max(y[0],y[1]),r);


  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));

  /*
    if (qa <=0. ) {
    l=y[0];
    r=y[0];
    }

    else if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
    else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


    *lout=l;   //a_L,j
    *rout=r;
    */

  if (qa <=0. ) {
    *lout=y[0];
    *rout=y[0];
  }  
  else {
    *lout = l;
    *rout = r;
  }

  //2.0 at top/bottom of a steep gradient 
  if (qd*(qd-qe)<0.0) *lout=3.0*y[0]-2.0*r;
  else *lout = l;

  if (qd*(qd+qe)<0.0) *rout=3.0*y[0]-2.0*l;
  else *rout = r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}






#define JONPARAREDUCE 0  //atch adjust // jon adjust

// Xiaoyue given on 03/25/05
// she realized sign error in 1st der's in para3()
// noted Matt's paper astro-ph/0503420 suggested CW1.6 uses 1/8 rather than 1/6
// I noted Matt uses MC for field variables and PPM+ for hydro variables
void para4(int pl, FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;


#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);

    if (s <=0.) dq[mm]=0.;
    //else dq[mm] = aDqvanl*sign(Dqc);
    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);

#elif(PARA2LIM == MC)

#if(0)
    // Jon's version
    dq[mm]=MINMOD(Dqc,MINMOD(Dqm,Dqp));
#else
    // Xioyue's version
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]= min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#endif



#elif(PARA2LIM == MINM_STEEPENER)

    // Xioyue's version (steepeneed version of MINM)
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);


#elif(PARA2LIM == MINM) // no steepener, normal MINM

#if(0)
    // Jon's version
    dq[mm] = MINMOD(0.5*Dqm,0.5*Dqp); // no steepening    
#elif(1)
    // Jon's steep version
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);
#elif(0)
    // Xioyue's version
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = 0.5*aDqm*sign(Dqc);
    else dq[mm]=0.5*aDqp*sign(Dqc);
#endif

#elif(PARA2LIM == NLIM) //w/o slope limiter

    dq[mm] = Dqc;
#endif
  }

#if(JONPARAREDUCE)
  //  if(pl==U1){
  if(pl!=RHO){
    if(
       (fabs(dq[-1]-dq[0])/(fabs(dq[-1])+fabs(dq[0])+SMALL)>0.1)||
       (fabs(dq[1]-dq[0])/(fabs(dq[1])+fabs(dq[0])+SMALL)>0.1)
       ){
      slope_lim_3points(MINM, y[-1], y[0], y[1], dq);
      *lout =y[0] - 0.5* (*dq);
      *rout=y[0] + 0.5* (*dq);
      return;
    }
  }

#endif

  /* CW1.6 */

  // modified as per Matt's paper
  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/8.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/8.0;


  l=max(min(y[0],y[-1]),l);
  l=min(max(y[0],y[-1]),l);
  r=max(min(y[0],y[1]),r);
  r=min(max(y[0],y[1]),r);


  // modified as per Matt's paper
  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }
  else{

    if (qd*(qd-qe)<0.0) 
      l=3.0*y[0]-2.0*r;


    if (qd*(qd+qe)<0.0) 
      r=3.0*y[0]-2.0*l;
  }


  *lout=l;   //a_L,j
  *rout=r;

  //  *dqleft=dq[-1];
  //  *dqcenter=dq[0];
  //  *dqright=dq[1];

}



FTYPE ftilde_orig(int dir, int shift, FTYPE **ypl ) 
{
  FTYPE Ftilde;
  FTYPE *P;
  FTYPE Sp;

  P = ypl[UU] + shift;

  Sp = fabs(P[1] - P[-1]) / (fabs(P[2] - P[-2]) + SMALL);
  Ftilde = max( 0, min( 1, 10. * (Sp - 0.6) ) );
  //dualfprintf( fail_file, "nstep = %ld, steppart = %d, icurr = %d, Ftilde1 = %21.15g\n", nstep, steppart, icurr, Ftilde );
  Ftilde *= (FTYPE)(fabs(P[1] - P[-1]) / max( min(fabs(P[1]), fabs(P[-1])), SMALL ) >= 1./3. );
  //dualfprintf( fail_file, "nstep = %ld, steppart = %d, icurr = %d, Ftilde2 = %21.15g\n", nstep, steppart, icurr, Ftilde );
  Ftilde *= (FTYPE)( ypl[UU+dir][1] - ypl[UU+dir][-1] < 0 );
  //dualfprintf( fail_file, "nstep = %ld, steppart = %d, icurr = %d, Ftilde3 = %21.15g\n", nstep, steppart, icurr, Ftilde );
  
  return( Ftilde );
}

FTYPE ftilde(int dir, int shift, FTYPE **ypl ) 
{
  FTYPE Ftilde;
  FTYPE *P;
  FTYPE Sp;

  P = ypl[UU] + shift;

  Sp = fabs(P[1] - P[-1]) / (fabs(P[2] - P[-2]) + SMALL);
  Ftilde = max( 0, min( 1, 10. * (Sp - 0.65) ) );
  //dualfprintf( fail_file, "nstep = %ld, steppart = %d, icurr = %d, Ftilde1 = %21.15g\n", nstep, steppart, icurr, Ftilde );
  Ftilde *= (FTYPE)(fabs(P[1] - P[-1]) / max( min(fabs(P[1]), fabs(P[-1])), SMALL ) >= 1./3. );
  //dualfprintf( fail_file, "nstep = %ld, steppart = %d, icurr = %d, Ftilde2 = %21.15g\n", nstep, steppart, icurr, Ftilde );
  Ftilde *= (FTYPE)( ypl[UU+dir][1] - ypl[UU+dir][-1] < 0 );
  //dualfprintf( fail_file, "nstep = %ld, steppart = %d, icurr = %d, Ftilde3 = %21.15g\n", nstep, steppart, icurr, Ftilde );
  
  return( Ftilde );
}


FTYPE  Ficalc(int dir, FTYPE **ypl)
{
  FTYPE ftilde( int dir, int shift, FTYPE **ypl );
  FTYPE *P;
  FTYPE *rho;
  int signdP;
  FTYPE Fi;

  P=ypl[UU];
  rho=ypl[RHO];

  signdP = (P[1] - P[-1] > 0) * 2 - 1;
  Fi = max( ftilde(dir, 0, ypl), ftilde(dir, -signdP, ypl) );

  return(Fi);
}



#define DOPPMREDUCE 1 //atch adjust
#define DOPPMCONTACTSTEEP 0 // doesn't seem to work

void parapl(int dir, FTYPE **ypl, FTYPE *loutpl, FTYPE *routpl)
{
  void para4(int pl, FTYPE *y, FTYPE *lout, FTYPE *rout);
  FTYPE etaicalc(FTYPE *rho, FTYPE *P, FTYPE *rhold, FTYPE *rhord);
  FTYPE  Ficalc(int dir, FTYPE **ypl);
  FTYPE *y;
  int pl;
  FTYPE *P;
  FTYPE *rho;
  FTYPE Fi;
  FTYPE etai;
  FTYPE rhold, rhord;



#if( DOPPMREDUCE )
  Fi = Ficalc(dir,ypl);
#else
  Fi = 0.0;
#endif


#if(DOPPMCONTACTSTEEP)
  P=ypl[UU];
  rho=ypl[RHO];
  // get contact indicator
  etai=etaicalc(rho,P,&rhold,&rhord);
#else
  etai = rhold = rhord = 0.0;
#endif

  //dualfprintf( fail_file, "nstep = %ld, steppart = %d, icurr = %d, Fi = %21.15g\n", nstep, steppart, icurr, Fi );

  PLOOPINTERP(pl){

    y=ypl[pl];
    
    para4(pl,y,&loutpl[pl],&routpl[pl]);


#if(DOPPMCONTACTSTEEP)
    if(pl==RHO){
      // assign density value
      loutpl[pl] = loutpl[pl] * ( 1.0 - etai ) + rhold*etai;
      routpl[pl] = routpl[pl] * ( 1.0 - etai ) + rhord*etai;
    }
#endif


#if( DOPPMREDUCE )
    loutpl[pl] = Fi * ypl[pl][0] + ( 1.0 - Fi ) * loutpl[pl];
    routpl[pl] = Fi * ypl[pl][0] + ( 1.0 - Fi ) * routpl[pl];
#endif
    
  }



}


// PPM steepener parameter, where etai=1 is steep and etai=0 is normal not steep
FTYPE etaicalc(FTYPE *rho, FTYPE *P, FTYPE *rhold, FTYPE *rhord)
{
  FTYPE dqlmono,dqrmono;
  FTYPE delta2rhol,delta2rhor;
  FTYPE etatilde;
  FTYPE *dqlcr;
  FTYPE a_dqlcr[3];
  FTYPE etai;
  int ii;


  dqlcr=&a_dqlcr[1];


  // not sure if those para4 dqlcr's are the right one's to use, so just compute them
#if(1)
  ii=-1;
  dqlcr[ii]=0.5*(rho[ii+1]-rho[ii-1]);

  ii=0;
  dqlcr[ii]=0.5*(rho[ii+1]-rho[ii-1]);

  ii=1;
  dqlcr[ii]=0.5*(rho[ii+1]-rho[ii-1]);
#endif

  // equation 26 in FLASH
  ii=-1;
  if( (rho[ii+1]-rho[ii])*(rho[ii]-rho[ii-1])>0.0){
    dqlmono=min(min(fabs(dqlcr[ii]),2.0*fabs(rho[ii]-rho[ii-1])),2.0*fabs(rho[ii]-rho[ii+1]))*copysign(dqlcr[ii],1.0);
  }
  else dqlmono=0.0;

  // equation 26 in FLASH
  ii=1;
  if( (rho[ii+1]-rho[ii])*(rho[ii]-rho[ii-1])>0.0){
    dqrmono=min(min(fabs(dqlcr[ii]),2.0*fabs(rho[ii]-rho[ii-1])),2.0*fabs(rho[ii]-rho[ii+1]))*copysign(dqlcr[ii],1.0);
  }
  else dqrmono=0.0;
      
  // equation 33 and 34 in FLASH
  *rhold=rho[ii-1]+0.5*dqlmono;
  *rhord=rho[ii+1]-0.5*dqrmono;

  // equation 35 in FLASH
  ii=-1;
  delta2rhol=SIXTH*( (rho[ii+1]-rho[ii]) - (rho[ii] - rho[ii-1]) );

  // equation 35 in FLASH
  ii=1;
  delta2rhor=SIXTH*( (rho[ii+1]-rho[ii]) - (rho[ii] - rho[ii-1]) );

  // equation 36 in FLASH
  // sign of denominator is important
#if(0)
  if(rho[1]-rho[-1]!=0.0){
    etatilde=-(delta2rhor-delta2rhol)/(rho[1]-rho[-1]);
  }
  else etatilde=0.0;
#else
  etatilde=-0.5*(delta2rhor-delta2rhol)/(fabs(rho[1]-rho[-1])+SMALL);
#endif

  // equation 37 in FLASH
  if(fabs(rho[1]-rho[-1])/max(min(rho[1],rho[-1]),SMALL)<0.01) etatilde=0.0;

  // equation 38 in FLASH
  if(delta2rhol*delta2rhor>0.0) etatilde=0.0;

  // equation 39 in FLASH
  if( (P[1]-P[-1])/max(min(P[1],P[-1]),SMALL)> 0.1*(rho[1]-rho[-1])/max(min(rho[1],rho[-1]),SMALL)) etatilde=0.0;

  // equation 40 in FLASH
  etai=max(0.0,min(20*(etatilde-0.05),1.0));

  return(etai);

}

