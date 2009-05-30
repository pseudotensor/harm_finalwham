#include "decs.h"



// this file includes metric dependent terms, including for initial
// condition routines for IC coords.

// crucially need to setup analytic form of gcov.  All rest can be done numerically or analytically if wanted.

// obtain gcov in primcoords of whichcoord type metric/coords
void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  void set_gcov_cylminkmetric    (FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void set_gcov_spcminkmetric    (FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void set_gcov_cartminkmetric   (FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void set_gcov_unigravity   (FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void set_gcov_htmetric         (FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void set_gcov_htmetric_accurate(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void set_gcov_ksmetric         (FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void set_gcov_blmetric         (FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void gcov2gcovprim(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert, FTYPE gcovprim[][NDIM], FTYPE *gcovpertprim);

  FTYPE V[NDIM];
  int j,k;


  // before here, X set by i,j or chosen
  bl_coord(X, V);

 
  if(whichcoord>=0){
    if(whichcoord==BLCOORDS){
      set_gcov_blmetric(V, gcov, gcovpert);
    }
    else if(whichcoord==KSCOORDS){
      set_gcov_ksmetric(V, gcov, gcovpert);
    }
    else if(whichcoord==HTMETRIC){
      set_gcov_htmetric(V, gcov, gcovpert);
    }
    else if(whichcoord==HTMETRICACCURATE){
      set_gcov_htmetric_accurate(V, gcov, gcovpert);
    }
    else if(whichcoord==CARTMINKMETRIC){
      set_gcov_cartminkmetric(V, gcov, gcovpert);
    }
    else if(whichcoord==UNIGRAVITY){
      set_gcov_unigravity(V, gcov, gcovpert);
    }
    else if(whichcoord==CYLMINKMETRIC){
      set_gcov_cylminkmetric(V, gcov, gcovpert);
    }
    else if(whichcoord==SPCMINKMETRIC){
      set_gcov_spcminkmetric(V, gcov, gcovpert);
    }
    else{
      dualfprintf(fail_file,"gcov_func(): no such whichcoord=%d\n",whichcoord);
      myexit(1);
    }
  }
  else{
    dualfprintf(fail_file,"your request makes no sense (i.e. can't get prim gcov from prim gcov\n");
    myexit(1);
  }

  //  DLOOP(j,k) { fprintf(stderr,"1gcov[%d][%d]=%21.15g\n",j,k,gcov[j][k]); fflush(stderr);}

  // whether to convert to prim coords
  if(getprim==1){
    // all the above are analytic, so have to convert to prim coords.
    gcov2gcovprim(X, V, gcov,gcovpert, gcov, gcovpert);
  }
  //  DLOOP(j,k) { fprintf(stderr,"2gcov[%d][%d]=%21.15g\n",j,k,gcov[j][k]); fflush(stderr);}

}


// obtain prim gcon in primcoords of whichcoord type metric/coords
void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
{
  void set_gcon_blmetric(FTYPE *V, FTYPE gcon[][NDIM]);
  void set_gcon_ksmetric(FTYPE *V, FTYPE gcon[][NDIM]);
  void gcon2gconprim(FTYPE *X, FTYPE *V, FTYPE gcon[][NDIM],FTYPE gconprim[][NDIM]);
  void matrix_inverse(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
  int j,k;
  FTYPE V[NDIM];


  if(whichcoord>=0){
    if(whichcoord==BLCOORDS){
      bl_coord(X, V);
      set_gcon_blmetric(V, gcon);
      if(getprim) gcon2gconprim(X, V, gcon,gcon);
    }
    else if(whichcoord==KSCOORDS){
      bl_coord(X, V);
      set_gcon_ksmetric(V,gcon);
      if(getprim) gcon2gconprim(X, V,gcon,gcon);
    }
    else if(whichcoord==HTMETRIC){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(gcov,gcon);
    }
    else if(whichcoord==HTMETRICACCURATE){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(gcov,gcon);
    }
    else if(whichcoord==CARTMINKMETRIC){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(gcov,gcon);
    }
    else if(whichcoord==UNIGRAVITY){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(gcov,gcon);
    }
    else if(whichcoord==CYLMINKMETRIC){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(gcov,gcon);
    }
    else if(whichcoord==SPCMINKMETRIC){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(gcov,gcon);
    }
    else{
      dualfprintf(fail_file,"gcon_func(): no such whichcoord=%d\n",whichcoord);
      myexit(1);
    }
  }
  else{
    dualfprintf(fail_file,"your request makes no sense (i.e. can't get prim gcon from prim gcon\n");
    myexit(1);
  }

}




// delta is simply how big the differencing is, should be small, but not so small to lead to errors due to erros in the metric itself (i.e. keep larger than machine precision)
//#define DELTA (NUMEPSILON*1000.0)
//#define DELTA 1.e-5
//#define DELTA 1.e-5
//#define DELTA (pow(NUMEPSILON,1.0/3.0))
//#define DELTA NUMSQRTEPSILON // as in NR's fdjac() -- smaller isn't always better
// how to generically set this?  Too high, even slightly (10^{-10} for long doubles) and connection is screwed)

// Avery mentions that long double trig. functions only return double precision answer.  see ~/research/utils/triglongdouble.c

#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
#define CONNDELTA 1E-5 // default -- seems to work pretty good generally to reduce max error
//#define CONNDELTA 5.4E-5 // default -- seems to work pretty good
//#define CONNDELTA 4.6E-5 // min of error for a specific case, but apparently not generally good
#elif(REALTYPE==LONGDOUBLETYPE)
//#define CONNDELTA 7.17E-6 // based on min of error for specific case
#define CONNDELTA 1E-6 // based on min of error for specific case
// polar region likes 6.5E-8 (min of error for specific case)
#endif

// connection not simply transformed -- so compute directly from final metric (primcoords)
void conn_func(int whichcoord, FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void set_conn_cylminkmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_spcminkmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_cartminkmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_unigravity(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_htmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_htmetric_accurate(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_ksmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_blmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);


  if(whichcoord==BLCOORDS){
    set_conn_blmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==KSCOORDS){
    set_conn_ksmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==HTMETRIC){
    set_conn_htmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==HTMETRICACCURATE){
    set_conn_htmetric_accurate(X,geom,conn,conn2);
  }
  else if(whichcoord==CARTMINKMETRIC){
    set_conn_cartminkmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==UNIGRAVITY){
    set_conn_unigravity(X,geom,conn,conn2);
  }
  else if(whichcoord==CYLMINKMETRIC){
    set_conn_cylminkmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==SPCMINKMETRIC){
    set_conn_spcminkmetric(X,geom,conn,conn2);
  }
  else{
    dualfprintf(fail_file,"conn_func(): no such whichcoord=%d\n",whichcoord);
    myexit(1);
  }
}



// obtain eomfunc (f_{(\nu)} factor: see how used in connection below and get_geometry())
// assumes X set before eomfun_func()
void eomfunc_func(int getprim, int whichcoord, FTYPE *X, FTYPE *eomfunc)
{

  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE gcovpertcoord[NDIM];
  FTYPE r,th;
  int j,k;
  FTYPE V[NDIM];
  void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);



  if(WHICHEOM==WITHGDET){
    gcov_func(getprim, whichcoord,X,gcovmcoord,gcovpertcoord); // actually returns primcoords version of whichcoord
    *eomfunc=gdet_func(gcovmcoord);
  }
  else if(WHICHEOM==WITHNOGDET){
    *eomfunc=1.0;
  }
  else if(WHICHEOM==WITHSINSQ){ // obviously coordinate dependent
    bl_coord(X, V);

    r=V[1];
    th=V[2];
    *eomfunc=sin(th)*sin(th);
  }
  else{
    dualfprintf(fail_file,"eomfunc_func(): no such WHICHEOM=%d\n",WHICHEOM);
    myexit(1);
  }
}



///////////////////////////
//
// g_{\mu\nu} (always analytic -- returns some coordinate system)
//
///////////////////////////



// needs M, J, and Q
// M: total M+dM mass of star
// J: total angular momentum
// Q: mass quadrapole moment
// external space-time of slowly rotating star, accurate to second order in $\Omega_{\star}$
// c=G=1
// NS mass sheds at v=0.8c (eq 28), allowed approx up to v=0.1c
// WD mass sheds at v=0.01c, allowed approx up to v=0.001c
// equally restrictive as measured by ratio of rotational energy to gravitational energy
void set_gcov_htmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE P2,Q12, Q22;
  FTYPE z,OO,RO,SS,AA,alphasq,betaphi;
  FTYPE SSpert,OOpert,AApert;
  FTYPE M,J,Q;
  FTYPE r,th;
  FTYPE ftemp;


  r=V[1];
  th=V[2];


#define FLATSPACE 0
#define NOSPIN 1
#define FULLHT 2

#define HTMETRICTYPE FLATSPACE

  // here G=c=1 and M is either 0 or 1 usually, such that
  // J=a M^2 , Q~J^2/M =a^2 M

  if(HTMETRICTYPE==FULLHT){
    M=1.0;
    J=a; // I\Omega = J
    Q=a*a; // approximately, perhaps, like Kerr geometry, if hard enough EOS
  }
  else if(HTMETRICTYPE==NOSPIN){
    M=1.0;
    J=0; // no spin of metric
    Q=0;
  }
  else if(HTMETRICTYPE==FLATSPACE){
    M=1E-100; // no mass of star, just flat space (GJ model)
    J=0;
    Q=0;
  }
  // surface radius is not constant radius if Q!=0
  // if Q=0, then surface radius is 3M?  depends on EOS.

  // let's set M=1, such that J is up to M^2 and Q is up to M^3 such that radius is in units of GM/c^2

  // For Sun: M/r<M/R~2E-6, J/r^2<J/R^2~1E-12 (J/M^2=4E-24), Q/r^3<Q/R^3<~1E-10 (Q/R^3=8E-28) and surface roughly spherical

  // Fastest pulsar PSR 1937+21 \tau=1.56ms , 640times/sec
  // Deepto Chakrabarty of MIT
  // Vela: PSR 0833-45 \tau=89.3ms
  // strongest: PSR 0329+54 \tau=0.715s
  // peak: 760times/sec (limit by GR?), 2-3X is theoretical limit but would mass shed.
  // formation spin rate: 30/sec


  // see Hartle & Thorne 1968 or Kim Lee Lee Lee 2005

  // Legendre polynomial
  P2=0.5*(3.0*cos(th)*cos(th)-1.0);

  z=r/M-1.0;

  // Mathematica's LegendreQ (associated Legendre of 2nd kind) is such that:
  //   Q12=Im[Legendre[2,1,z]] and Q22=-Re[Legendre[2,2,z]]
  Q12=sqrt(z*z-1.0)*( (3.0*z*z-2.0)/(z*z-1.0)-3.0/2.0*z*log((z+1)/(z-1)) );

  Q22=(3.0/2.0*(z*z-1.0)*log((z+1)/(z-1))-(3.0*z*z*z-5.0*z)/(z*z-1.0) );



  // metric stuff
  OOpert=-2.0*M/r+2.0*J*J/(r*r*r*r);
  OO=1.0+OOpert;

  RO=1.0+2.0*(J*J/(M*r*r*r)*(1.0+M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*Q22)*P2;
  
  SSpert=-2.0*(J*J/(M*r*r*r)*(1.0-5.0*M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*Q22)*P2;
  SS=1.0+SSpert;

  AApert=+2.0*(-J*J/(M*r*r*r)*(1.0+2.0*M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*( (2.0*M/(sqrt(r*r*(1.0-2.0*M/r))))*Q12 - Q22  ))*P2;
  AA=1.0+AApert;
  
  alphasq=OO*RO;

  betaphi=-2.0*J/(r*r*r); // note that omega_{zamo}=-betaphi

  gcov[PH][PH] = r*r*AA*sin(th)*sin(th) ;
  gcovpert[PH] = gcov[PH][PH]-1.0;
  
  gcov[TT][TT] = -(alphasq-betaphi*betaphi*gcov[PH][PH]) ;
  ftemp = (1.0-alphasq);
  // g_{tt} + 1
  gcovpert[TT] = ftemp + betaphi*betaphi*gcov[PH][PH] ;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = betaphi*gcov[PH][PH] ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = SS/OO ;
  // SS/OO = (SSpert+1.0)/OO = SSpert/OO + 1.0/OO , so SS/OO-1.0 = SSpert/OO + (1.0/OO-1.0)
  ftemp=(1.0/OO - 1.0);
  gcovpert[RR] = SSpert/OO + ftemp ; // might help with radial gravity terms

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = r*r*AA ;
  // r^2 AA - 1  = r^2(AApert+1.0) -1 = r^2 AApert + r^2 -1, so not a big issue with catastrophic cancellation in non-rel limit
  // probably have to define background as being r^2 to have gravity work at large radii (GODMARK)
  gcovpert[TH] = gcov[TH][TH] - 1.0;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;

  
								     



}

#undef HTMETRICTYPE

// from Berti, White, Maniopoulou, and Bruni (2005)
void set_gcov_htmetric_accurate(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE M,J,Q;
  FTYPE u,p,A1,A2,W,F1,F2,H1,H2,L,G1;
  FTYPE r,th;
  FTYPE ftemp;


  r=V[1];
  th=V[2];



#define FLATSPACE 0
#define NOSPIN 1
#define FULLHT 2

  //#define HTMETRICTYPE FULLHT
  //#define HTMETRICTYPE FLATSPACE
#define HTMETRICTYPE NOSPIN

  // here G=c=1 and M is either 0 or 1 usually, such that
  // J=a M^2 , Q~J^2/M =a^2 M


  // these J,Q are actually dimensionless

  if(HTMETRICTYPE==FULLHT){
    M=1.0;
    J=a; // I\Omega = J (J=trueJ/M^2)
    Q=a*a; // approximately, perhaps, like Kerr geometry, if hard enough EOS
    // Q=Qtrue/M^3
  }
  else if(HTMETRICTYPE==NOSPIN){
    M=1.0;
    J=0; // no spin of metric
    Q=0;
  }
  else if(HTMETRICTYPE==FLATSPACE){
    M=1E-100; // no mass of star, just flat space (GJ model)
    J=0;
    Q=0;
  }

  L=(80.0*pow(M,6)+8.0*pow(M,4)*r*r+10.0*M*M*M*r*r*r+20.0*M*M*pow(r,4)-45.0*M*pow(r,5)+15.0*pow(r,6));

  u=cos(th);
  p=1.0/(8.0*M*pow(r,4)*(r-2.0*M)); // typo in paper for parenthesis


  A1=(15.0*r*(r-2.0*M)*(1.0-3.0*u*u))/(16.0*M*M)*log(r/(r-2.0*M));
  A2=(15.0*(r*r-2.0*M*M)*(3.0*u*u-1.0))/(16.0*M*M)*log(r/(r-2.0*M));
  
  W=(r-M)*(16.0*pow(M,5)+8.0*pow(M,4)*r-10*M*M*r*r*r-30.0*M*pow(r,4)+15.0*pow(r,5))+u*u*(48.0*pow(M,6)-8.0*pow(M,5)*r-24.0*pow(M,4)*r*r-30.0*M*M*M*r*r*r-60.0*M*M*pow(r,4)+135.0*M*pow(r,5)-45.0*pow(r,6));

  F1=-p*W+A1;
  F2=5.0*r*r*r*p*(3.0*u*u-1.0)*(r-M)*(2.0*M*M+6.0*M*r-3.0*r*r)-A1;

  H1=A2+(1.0/(8.0*M*r*r*r*r))*(1.0-3.0*u*u) * (16.0*pow(M,5)+8.0*pow(M,4)*r-10.0*M*M*r*r*r+15.0*M*pow(r,4)+15.0*pow(r,5)) ;
  H2=-A2 + (1.0/(8.0*M*r))*5.0*(1.0-3.0*u*u)*(2.0*M*M-3.0*M*r-3.0*r*r);

  G1=p*((L-72.0*M*M*M*M*M*r)-3.0*u*u*(L-56.0*M*M*M*M*M*r))-A1;
  




  ftemp  = J*J*F1-Q*F2;
  gcov[TT][TT] = -(1.0-2.0*M/r)*(1.0+ftemp);
  // g_{tt} + 1
  //  gcovpert[TT] =  -(J*J*F1-Q*F2) + (2.0*M/r)*(1.0+J*J*F1-Q*F2);
  gcovpert[TT] =  -ftemp +(1.0+ftemp)*(2.0*M/r);



  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = (2.0*J*M*M/r)*sin(th)*sin(th);
    
  gcov[RR][TT] = 0.0 ;

  ftemp=J*J*G1+Q*F2;
  gcov[RR][RR] = (1.0/(1.0-2.0*M/r))*(1.0+ftemp);
  // g_{rr} - 1
  gcovpert[RR] = (ftemp+2.0*M/r)/(1.0-2.0*M/r); // should help with radial gravity in non-rel regime

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = r*r*(1.0+J*J*H1-Q*H2) ;
  gcovpert[TH] = gcov[TH][TH]-1.0;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;
  gcov[PH][PH] = gcov[TH][TH]*sin(th)*sin(th);
  gcovpert[PH] = gcov[PH][PH]-1.0;

  
								     



}
#undef HTMETRICTYPE



// Cartesian minkowski
// (t,x,y,z)
void set_gcov_cartminkmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  
  gcov[TT][TT] = -1.0;
  gcovpert[TT] = 0.0;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = 0.0 ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = 1.0 ;
  gcovpert[RR] = 0.0;

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = 1.0 ;
  gcovpert[TH] = 0.0;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;
  gcov[PH][PH] = 1.0 ;
  gcovpert[PH] = 0.0;

}


// (t,x,y,z)
void set_gcov_unigravity(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{

  FTYPE phi;
  FTYPE x,y,z;
  FTYPE FORCEMAG;

  x=V[1];
  y=V[2];
  z=V[3];

	//this is the hard-cored value for test 28, right? SASMARK
  FORCEMAG=0.1/(coordparams.timescalefactor*coordparams.timescalefactor); // strength of force, scales like v^2

  // uniform force field in +y direction (i.e. dphi/dy = FORCEMAG, so Force = -FORCEMAG)
  phi = FORCEMAG*y;

  gcov[TT][TT] = -1.0-2.0*phi;
  gcovpert[TT] = -2.0*phi;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = 0.0 ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = 1.0-2.0*phi ;
  gcovpert[RR] = -2.0*phi;

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = 1.0-2.0*phi ;
  gcovpert[TH] = -2.0*phi;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;

  gcov[PH][PH] = 1.0-2.0*phi ;
  gcovpert[PH] = -2.0*phi;

  //  dualfprintf(fail_file,"tsf=%21.15g :: %21.15g %21.15g :: %21.15g %21.15g %21.15g %21.15g\n",coordparams.timescalefactor,FORCEMAG,phi,gcovpert[TT],gcovpert[RR],gcovpert[TH],gcovpert[PH]);



}



// (t,R,z,\phi)
void set_gcov_cylminkmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{

  FTYPE r,Mass,RSTAR;
  FTYPE phi;
  FTYPE R,z;



  R=V[1];
  z=V[2];




  Mass=0.0;
  RSTAR=0.1; // could use dxdxp[RR][RR]*dx[RR] to get size
  r=sqrt(R*R+z*z);

  if(r<RSTAR){
    phi = (Mass/(2.0*RSTAR))*((r/RSTAR)*(r/RSTAR)-3.0);
  }
  else phi = -Mass/r;

  
  gcov[TT][TT] = -1.0-2.0*phi;
  gcovpert[TT] = -2.0*phi;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = 0.0 ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = 1.0-2.0*phi*R*R/(r*r); 
  gcovpert[RR] = -2.0*phi*R*R/(r*r); 

  gcov[RR][TH] = -2.0*phi*R*z/(r*r) ; // order v^4
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = 1.0 -2.0*phi*z*z/(r*r) ;
  gcovpert[TH] =-2.0*phi*z*z/(r*r) ; // order v^4

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;
  gcov[PH][PH] = R*R;
  gcovpert[PH] = gcov[PH][PH]-1.0; // doesn't matter that -1.0 subtracted

}






// (t,r,\theta,\phi)
void set_gcov_spcminkmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE r,th;

  r=V[1];
  th=V[2];
  
  
  gcov[TT][TT] = -1.0;
  gcovpert[TT] = 0.0;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = 0.0 ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = 1.0; 
  gcovpert[RR] = 0.0;

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = r*r ;
  gcovpert[TH] = gcov[TH][TH]-1.0;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;

  gcov[PH][PH] = (r*sin(th))*(r*sin(th));
  gcovpert[PH] = gcov[PH][PH]-1.0;

}



// (~t,r,\theta,~\phi)
void set_gcov_ksmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE sth, cth, s2, a2, r2, r3, DD, mu;
  FTYPE rho2;
  FTYPE r,th;
  FTYPE mysin(FTYPE th);
  FTYPE mycos(FTYPE th);



  r=V[1];
  th=V[2];


  sth=mysin(th);
  cth=mycos(th);


  s2 = sth * sth;
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;
  rho2 = r * r + a * a * cth * cth;


#define ks_gcov00 (-1. + 2.*r/rho2)
#define ks_gcov01 (2.*r/rho2)
#define ks_gcov02 (0)
#define ks_gcov03 (-2.*a*r*s2/rho2)
#define ks_gcov10 (ks_gcov01)
#define ks_gcov11 (1. + 2.*r/rho2)
#define ks_gcov12 (0)
#define ks_gcov13 (-a*s2*(1. + 2.*r/rho2))
#define ks_gcov20 (0)
#define ks_gcov21 (0)
#define ks_gcov22 (rho2)
#define ks_gcov23 (0)
#define ks_gcov30 (ks_gcov03)
#define ks_gcov31 (ks_gcov13)
#define ks_gcov32 (0)
#define ks_gcov33 (s2*(rho2 + a*a*s2*(1. + 2.*r/rho2)))


  gcov[TT][TT] = ks_gcov00 ;
  gcovpert[TT] = 2.*r/rho2;


  gcov[TT][RR] = ks_gcov01 ;
  gcov[TT][TH] = ks_gcov02 ;
  gcov[TT][PH] = ks_gcov03 ;
    
  gcov[RR][TT] = ks_gcov10 ;

  gcov[RR][RR] = ks_gcov11 ;
  gcovpert[RR] = 2.*r/rho2 ;

  gcov[RR][TH] = ks_gcov12 ;
  gcov[RR][PH] = ks_gcov13 ;
    
  gcov[TH][TT] = ks_gcov20 ;
  gcov[TH][RR] = ks_gcov21 ;

  gcov[TH][TH] = ks_gcov22 ;
  gcovpert[TH] = gcov[TH][TH]-1.0;

  gcov[TH][PH] = ks_gcov23 ;
    
  gcov[PH][TT] = ks_gcov30 ;
  gcov[PH][RR] = ks_gcov31 ;
  gcov[PH][TH] = ks_gcov32 ;

  gcov[PH][PH] = ks_gcov33 ;
  gcovpert[PH] = gcov[PH][PH]-1.0;

}


// (t,r,\theta,\phi)
void set_gcov_blmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE sth, cth, s2, a2, r2, r3, DD, mu;
  FTYPE mupert,DDpert;
  FTYPE r,th;

  r=V[1];
  th=V[2];




  cth = cos(th);
  sth=sin(th);


  s2 = sth * sth;
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DDpert =- 2. / r + a2 / r2; 
  DD = 1. + DDpert;
  mupert=+ a2 * cth * cth / r2;
  mu = 1. + mupert;



#define bl_gcov00 (-(1. - 2./(r*mu)))
#define bl_gcov01 (0)
#define bl_gcov02 (0)
#define bl_gcov03 (-2.*a*s2/(r*mu))
#define bl_gcov10 (0)
#define bl_gcov11 (mu/DD)
#define bl_gcov12 (0)
#define bl_gcov13 (0)
#define bl_gcov20 (0)
#define bl_gcov21 (0)
#define bl_gcov22 (r2*mu)
#define bl_gcov23 (0)
#define bl_gcov30 (bl_gcov03)
#define bl_gcov31 (0)
#define bl_gcov32 (0)
#define bl_gcov33 (r2*sth*sth*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu)))

  gcov[TT][TT] = bl_gcov00 ;
  gcovpert[TT] = 2./(r*mu);

  gcov[TT][RR] = bl_gcov01 ;
  gcov[TT][TH] = bl_gcov02 ;
  gcov[TT][PH] = bl_gcov03 ;
    
  gcov[RR][TT] = bl_gcov10 ;

  gcov[RR][RR] = bl_gcov11 ;
  // mu/DD = (1.0+mupert)/(1.0+DDpert) = 
  gcovpert[RR] = (DDpert - mupert) / (1.0 + DDpert);

  gcov[RR][TH] = bl_gcov12 ;
  gcov[RR][PH] = bl_gcov13 ;
    
  gcov[TH][TT] = bl_gcov20 ;
  gcov[TH][RR] = bl_gcov21 ;

  gcov[TH][TH] = bl_gcov22 ;
  gcovpert[TH] = gcov[TH][TH]-1.0;

  gcov[TH][PH] = bl_gcov23 ;
    
  gcov[PH][TT] = bl_gcov30 ;
  gcov[PH][RR] = bl_gcov31 ;
  gcov[PH][TH] = bl_gcov32 ;

  gcov[PH][PH] = bl_gcov33 ;
  gcovpert[PH] = gcov[PH][PH]-1.0;

}


///////////////////////////
//
// g^{\mu\nu} (analytic or numerical)
//
///////////////////////////


// (~t,r,\theta,~\phi)
void set_gcon_ksmetric(FTYPE *V, FTYPE gcon[][NDIM])
{
  FTYPE sth, cth, s2, a2, r2, r3, DD, mu;
  FTYPE rho2;
  FTYPE r,th;



  r=V[1];
  th=V[2];

  cth = cos(th);
  sth=sin(th);


  s2 = sth * sth;
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;
  rho2 = r * r + a * a * cth * cth;


#define ks_gcon00 (-(1.+2.*r/rho2))
#define ks_gcon01 (2.*r/rho2)
#define ks_gcon02 (0)
#define ks_gcon03 (0)
#define ks_gcon10 (ks_gcon01)
#define ks_gcon11 ((r*(r-2.)+a*a)/rho2)
#define ks_gcon12 (0)
#define ks_gcon13 (a/rho2)
#define ks_gcon20 (ks_gcon02)
#define ks_gcon21 (ks_gcon12)
#define ks_gcon22 (1./rho2)
#define ks_gcon23 (0)
#define ks_gcon30 (ks_gcon03)
#define ks_gcon31 (ks_gcon13)
#define ks_gcon32 (ks_gcon23)
#define ks_gcon33 (1./(rho2*s2))


  gcon[TT][TT] = ks_gcon00 ;
  gcon[TT][RR] = ks_gcon01 ;
  gcon[TT][TH] = ks_gcon02 ;
  gcon[TT][PH] = ks_gcon03 ;
    
  gcon[RR][TT] = ks_gcon10 ;
  gcon[RR][RR] = ks_gcon11 ;
  gcon[RR][TH] = ks_gcon12 ;
  gcon[RR][PH] = ks_gcon13 ;
    
  gcon[TH][TT] = ks_gcon20 ;
  gcon[TH][RR] = ks_gcon21 ;
  gcon[TH][TH] = ks_gcon22 ;
  gcon[TH][PH] = ks_gcon23 ;
    
  gcon[PH][TT] = ks_gcon30 ;
  gcon[PH][RR] = ks_gcon31 ;
  gcon[PH][TH] = ks_gcon32 ;
  gcon[PH][PH] = ks_gcon33 ;

}

// (t,r,\theta,\phi)
void set_gcon_blmetric(FTYPE *V, FTYPE gcon[][NDIM])
{
  FTYPE sth, cth, s2, a2, r2, r3, DD, mu;
  FTYPE r,th;



  r=V[1];
  th=V[2];

  cth = cos(th);
  sth=sin(th);


  s2 = sth * sth;
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;


#define bl_gcon00 (-1. - 2.*(1. + a2/r2)/(r*DD*mu))
#define bl_gcon01 (0)
#define bl_gcon02 (0)
#define bl_gcon03 (-2.*a/(r3*DD*mu))
#define bl_gcon10 (0)
#define bl_gcon11 (DD/mu)
#define bl_gcon12 (0)
#define bl_gcon13 (0)
#define bl_gcon20 (0)
#define bl_gcon21 (0)
#define bl_gcon22 (1./(r2*mu))
#define bl_gcon23 (0)
#define bl_gcon30 (bl_gcon03)
#define bl_gcon31 (0)
#define bl_gcon32 (0)
#define bl_gcon33 ((1. - 2./(r*mu))/(r2*sth*sth*DD))


  gcon[TT][TT] = bl_gcon00 ;
  gcon[TT][RR] = bl_gcon01 ;
  gcon[TT][TH] = bl_gcon02 ;
  gcon[TT][PH] = bl_gcon03 ;
    
  gcon[RR][TT] = bl_gcon10 ;
  gcon[RR][RR] = bl_gcon11 ;
  gcon[RR][TH] = bl_gcon12 ;
  gcon[RR][PH] = bl_gcon13 ;
    
  gcon[TH][TT] = bl_gcon20 ;
  gcon[TH][RR] = bl_gcon21 ;
  gcon[TH][TH] = bl_gcon22 ;
  gcon[TH][PH] = bl_gcon23 ;
    
  gcon[PH][TT] = bl_gcon30 ;
  gcon[PH][RR] = bl_gcon31 ;
  gcon[PH][TH] = bl_gcon32 ;
  gcon[PH][PH] = bl_gcon33 ;

}

///////////////////////////
//
// CONNECTIONS (analytic or numerical)
//
///////////////////////////

void set_conn_cylminkmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  int i,j,k;
  FTYPE gcovmid[NDIM][NDIM];
  FTYPE gcovpertmid[NDIM];
  FTYPE gdetmid;

	FTYPE dxdxp[NDIM][NDIM]; //atch
	FTYPE V[NDIM];           //atch


	if(defcoord==UNIFORMCOORDS){ // uniform grid  SUPERSASMARK
		// could directly use gdet in global memory
		// only works for X1=R and X2=z
		gcov_func(1,CYLMINKMETRIC,X, gcovmid,gcovpertmid);
		gdetmid=gdet_func(gcovmid);

		bl_coord( X, V );  //actually, dxdxprim() does not use V or X since metric is uniform

	  // see transforms.c and mettometp() and see gcov2gcovprim()  //atch
		dxdxprim(X, V, dxdxp);                                       //atch

		for (k = 0; k < NDIM; k++) conn2[k]= 0.0;
		//conn2[RR]=-1.0/gdetmid;  //wrong as well... shouldn't it be 0 since there is no 2nd connection when WITHGDET?  SUPERSASMARK
 

		for (i = 0; i < NDIM; i++)
			for (j = 0; j < NDIM; j++)
				for (k = 0; k < NDIM; k++) {
					conn[i][j][k] = 0.;
				}
				conn[PH][RR][PH]=1.0 / X[1]; //1.0/gdetmid; //apparently, wrong because should not care about the 2nd dimension  SUPERSASMARK
				conn[PH][PH][RR]=1.0 / X[1]; //1.0/gdetmid; //apparently, wrong because should not care about the 2nd dimension  SUPERSASMARK
				conn[RR][PH][PH]= - X[1] * dxdxp[3][3] * dxdxp[3][3]; //-gdetmid;
	}
	else{
		conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
	}


  
}

// only works for X1=R and X2=z
void set_conn_cartminkmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  int i,j,k;

  if(defcoord==UNIFORMCOORDS){// uniform grid
    for (k = 0; k < NDIM; k++) conn2[k]= 0.0;
    
    for (i = 0; i < NDIM; i++)
      for (j = 0; j < NDIM; j++)
	for (k = 0; k < NDIM; k++) {
	  conn[i][j][k] = 0.;
	}
  }
  else{
    conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  }

}

void set_conn_unigravity(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}


void set_conn_spcminkmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}


void set_conn_htmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}

void set_conn_htmetric_accurate(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}

void set_conn_ksmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void mks_conn_func(FTYPE *X, struct of_geom *geom,
		      FTYPE lconn[][NDIM][NDIM],FTYPE *conn2);

  //  FTYPE DELTA; // for debug below

  // the analytic form can be used with WHICHEOM=WITHNOGDET
  // determine for which we have analytic expressions
  // this currently competes with mks_source_conn in phys.c
  if((!ANALYTICCONNECTION)||(defcoord!=LOGRSINTH)) conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  else mks_conn_func(X,geom,conn,conn2);


  /*
  // some debug stuff
  //  if((geom->i==62)&&(geom->j==32)){ // where c002 is bad
  if(0&&(geom->i==32)&&(geom->j==0)){ // where c023 is bad
    for(DELTA=1E-15;DELTA<1E5;DELTA*=1.01){
      conn_func_numerical1(DELTA,X, geom, conn, conn2);
      dualfprintf(fail_file,"%30.20Lg %30.20Lg\n",DELTA,conn[0][2][3]); fflush(fail_file);
      //      dualfprintf(fail_file,"%21.15g %21.15g\n",DELTA,conn[0][2][3]); fflush(fail_file);
    }
    exit(0);
  }
  else{
    conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  }
  //mks_conn_func(X,geom,conn,conn2);
  */

}


void set_conn_blmetric(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}














// ///////////////////////////////////////////////////////////
// 
// below are independent of user choice of metric/coords/grid
// 
// ///////////////////////////////////////////////////////////



// find the con/cov forms of the chosen metric
void gset(int getprim, int whichcoord, int i, int j, int k, struct of_geom *ptrgeom)
{
  FTYPE X[NDIM];
  struct of_geom tempgeom;
  extern void assign_eomfunc(struct of_geom *geom, FTYPE eomfuncgen);
  void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);

  FTYPE eomfuncgen;

  ptrgeom->i=i;
  ptrgeom->j=j;
  ptrgeom->k=k;
  ptrgeom->p=CENT;
  icurr=i;
  jcurr=j;
  jcurr=k;
  pcurr=CENT;
  
  if(whichcoord>=0){
    coord(i, j, k, CENT, X);
    gcov_func(getprim,whichcoord,X,gengcov,gengcovpert);
    ptrgeom->g=gdet_func(gengcov); // must come after gcov_func() above
    gcon_func(getprim,whichcoord,X,gengcov,gengcon); // must come after gcov_func() above
    
    eomfunc_func(getprim, whichcoord,X,&eomfuncgen);
    assign_eomfunc(ptrgeom,eomfuncgen); // must come after assigning ptrgeom->g above
    
    ptrgeom->gcov=gengcov; // pointer
    ptrgeom->gcon=gengcon; // pointer
    ptrgeom->gcovpert=gengcovpert; // pointer
  }
  else if(whichcoord==PRIMECOORDS){ // special case
    get_geometry(i,j,k,CENT,ptrgeom);
  }
  else{
    dualfprintf(fail_file,"gset(): no such whichcoord=%d\n",whichcoord);
    myexit(1);
  }

}


void gcov2gcovprim(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert, FTYPE gcovprim[][NDIM], FTYPE *gcovpertprim)
//void gcov2gcovprim(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM],FTYPE gcovprim[][NDIM])
{
  int j, k, l, m;
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE tmp[NDIM][NDIM];
  FTYPE ftemp1,ftemp2;
  int q;

  // now take term by term:
  // g_{u v} = \vec{e_{\mu}}\cdot\vec{e_{\nu}} 
  //           * (dx/dx')_{mu} * (dx/dx')_{\nu} =
  //          \vec{e'_{\mu}}\cdot\vec{e'_{\nu}} 

  // dx/dx' where '=prim coords (i.e. nonuni coords)
  dxdxprim(X, V, dxdxp);

  DLOOP(j,k) tmp[j][k] = 0.;
  DLOOP(j,k) for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
    // g_{mup nup} = g_{mu nu} T^mu_mup T^nu_nup
    // where T^mu_mup == dx^mu[BL]/dx^mup[KSP uni grid]
    tmp[j][k] += gcov[l][m] * dxdxp[l][j] * dxdxp[m][k];
  }
  // use tmp since gcon might be same address as gconprim
  DLOOP(j,k) gcovprim[j][k] = tmp[j][k];

  // perturbed terms
  // SUPERGODMARK GODMARK: This only works for non-rel gravity if dxdxp is diagonal.  Not sure what to do in general
  DLOOPA(q){
    //    dualfprintf(fail_file,"gcovpert[%d]=%21.15g\n",q,gcovpert[q]);


    // get q-q term
    // -1 + g_{q'q'} for q!=TT, else  SOMECONSTANT + g_{q'q'} for q=TT

    // for spatial parts, only care about deviations from constant in getting connection coefficients
    // avoids catastrophic cancellation (GODMARK)
    if( (q!=TT)&&(defcoord==UNIFORMCOORDS)) ftemp1=0.0;
    // then unsure how to handle in general, so just leave alone (SUPERGODMARK GODMARK)
    // problem is with machine error in dxdxp leading to apparently large connection coefficients due to that machine error
    // also don't know ahead of time what constant to choose for prim coordinate quantities
    // and don't know how to handle variations in dxdxp, since then no constant to subtract off, but grid accelerations apparently wash out non-rel gravity accelerations?
    // except for time term, which is solid
    else ftemp1=(-1.0 + dxdxp[q][q] * dxdxp[q][q]);
    if(q==TT) ftemp1 *=-1.0;
    gcovpertprim[q]  = (gcovpert[q] * dxdxp[q][q] * dxdxp[q][q]) + ftemp1;

    // now add 15 other terms
    ftemp2 = 0.;
    for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
      if((l!=q)&&(m!=q)) ftemp2+= gcov[l][m] * dxdxp[l][q] * dxdxp[m][q];
    }
    // add other 15 terms to answer for total of 16 terms
    gcovpertprim[q]+=ftemp2;

    
    //    dualfprintf(fail_file,"dxdxp[%d][%d]=%21.15g\n",q,q,dxdxp[q][q]);
    //    dualfprintf(fail_file,"ftemp1[%d]=%21.15g ftemp2[%d]=%21.15g gcovpertprim[%d]=%21.15g\n",q,ftemp1,q,ftemp2,q,gcovpertprim[q]);
  }

  // override for now
  gcovpertprim[TT]=gcovprim[TT][TT]-mink(TT,TT);
  SLOOPA(q){
    gcovpertprim[q]=gcovprim[q][q]-mink(q,q);
  }

 

}

void gcon2gconprim(FTYPE *X, FTYPE *V, FTYPE gcon[][NDIM],FTYPE gconprim[][NDIM])
{
  int j, k, l, m;
  FTYPE dxdxp[NDIM][NDIM],idxdxp[NDIM][NDIM];
  FTYPE tmp[NDIM][NDIM];

  // see transforms.c and mettometp() and see gcov2gcovprim()
  dxdxprim(X, V, dxdxp);
  matrix_inverse(dxdxp,idxdxp);
  
  DLOOP(j,k) tmp[j][k] = 0.;
  DLOOP(j,k) for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
    tmp[j][k] += idxdxp[j][l] * idxdxp[k][m] * gcon[l][m] ;
  }
  // use tmp since gcon might be same address as gconprim
  DLOOP(j,k) gconprim[j][k] = tmp[j][k];

}


// find determinant in general of a metric
/* assumes gcov has been set first; returns determinant */

// determinant not simply transformed from analytic function -- so no analytic form possible yet
FTYPE gdet_func(FTYPE gcov[][NDIM])
{
  static int firstc = 1;
  static FTYPE **tmp;
  FTYPE d;
  int j, k, indx[NDIM];

  if (firstc) {
    tmp = dmatrix(1, NDIM, 1, NDIM);
    firstc = 0;
  }


  DLOOP(j,k) tmp[j + 1][k + 1] = gcov[j][k];
  if(ludcmp(tmp, NDIM, indx - 1, &d)>=1){
    dualfprintf(fail_file,"ludcmp failure\n");
    myexit(1);
  }
  // below from 1..NDIM due to ludcmp requiring 1..N
  for (j = 1; j <= NDIM; j++)
    d *= tmp[j][j];

  if(d>=0){
    dualfprintf(fail_file,"Metric has bad signature: d=%21.15g\n",d);
    DLOOP(j,k) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",j,k,gcov[j][k]);
    myexit(1);
  }
  else{  return (sqrt(-d)); }

  return(-1); // shouldn't ever get here
}

/* invert gcov to get gcon */
// can be used to invert any 2nd rank tensor (symmetric or not)
// actually returns the inverse transpose, so if
// gcov=T^j_k then out pops (iT)^k_j such that T^j_k (iT)^k_l = \delta^j_l
void matrix_inverse(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
{
  static int firstc = 1;
  int j, k;
  static FTYPE **tmp;

  if (firstc) {
    tmp = dmatrix(1, NDIM, 1, NDIM);
    firstc = 0;
  }

  DLOOP(j,k) tmp[j + 1][k + 1] = gcov[j][k];
  gaussj(tmp, NDIM, NULL, 0);
  // assign but also transpose (shouldn't do in general, confusing)
  //DLOOP(j,k) gcon[j][k] = tmp[k + 1][j + 1];
  DLOOP(j,k) gcon[j][k] = tmp[j + 1][k + 1];
}





/* 
   this gives the connection coefficient \Gamma^{i}_{j,k} =
   conn[..][i][j][k] where i,j,k = {0,1,2,3} corresponds to {t,r,theta,phi} 
 */

/*
  we also compute the 2nd connection:
  -d/dj(ln(\detg))
*/


#define GAMMIEDERIVATIVE 0
#define NUMREC 1


// GODMARK: had problems with large run (jetnewnoenv,jetnew on sauron)
//#define DERTYPE NUMREC
#define DERTYPE GAMMIEDERIVATIVE

/* NOTE: parameter hides global variable */
void conn_func_numerical1(FTYPE DELTA, FTYPE *X, struct of_geom *geom,
													FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
	int i, j, k, l;
	FTYPE tmp[NDIM][NDIM][NDIM];
	FTYPE Xh[NDIM], Xl[NDIM];
	FTYPE gh[NDIM][NDIM];
	FTYPE gl[NDIM][NDIM];
	FTYPE lngdeth,lngdetl;
	void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
	FTYPE dfridr(FTYPE (*func)(FTYPE*,int,int), FTYPE *X,int ii, int jj, int kk);
	FTYPE gcov_func_mcoord(FTYPE* X, int i, int j);
	FTYPE lngdet_func_mcoord(FTYPE* X, int i, int j);
	FTYPE gcovperth[NDIM],gcovpertl[NDIM];

	//	dualfprintf(fail_file,"myi=%d myj=%d myk=%d\n",geom->i,geom->j,geom->k);

	// gabc_{ijk}=dg_{ij}/dx^k
	// gammie derivative
	if(DERTYPE==GAMMIEDERIVATIVE){
		for (k = 0; k < NDIM; k++) {

			for (l = 0; l < NDIM; l++){
				Xh[l] = X[l];
				Xl[l] = X[l];
			}
			Xh[k] += DELTA;
			Xl[k] -= DELTA;

			if(WHICHEOM!=WITHGDET){
				lngdeth=lngdet_func_mcoord(Xh,0,0); // doesn't use 0,0
				lngdetl=lngdet_func_mcoord(Xl,0,0); // doesn't use 0,0
				conn2[k]= (lngdeth - lngdetl) / (Xh[k] - Xl[k]);

			}
			else{
				conn2[k]=0.0; // no 2nd connection then
			}

			gcov_func(1,MCOORD,Xh, gh,gcovperth);
			gcov_func(1,MCOORD,Xl, gl,gcovpertl);
			for (i = 0; i < NDIM; i++) for (j = 0; j < NDIM; j++){
				// d(1+g_{tt}) -> dg_{tt}, so can use 1+g_{tt} for accurate non-relativistic gravity
			  //if(i==j) conn[i][j][k] = (gcovperth[i] - gcovpertl[i]) / (Xh[k] - Xl[k]);
			  // else 
			  conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
			}
		}
	}
	else if(DERTYPE==NUMREC){
		for (k = 0; k < NDIM; k++) {

			if(WHICHEOM!=WITHGDET){
				conn2[k]= dfridr(lngdet_func_mcoord,X,0,0,k); // 0,0 not used
			}
			else conn2[k]=0.0; // then no 2nd connection

			for (i = 0; i < NDIM; i++) for (j = 0; j < NDIM; j++){
				conn[i][j][k] = dfridr(gcov_func_mcoord,X,i,j,k);
			}
		}
	}


	/* now rearrange to find \Gamma_{ijk}=1/2*(gabc_{jik}+gabc_{kij}-gabc_{kji}) */
	for (i = 0; i < NDIM; i++)
		for (j = 0; j < NDIM; j++)
			for (k = 0; k < NDIM; k++)
				tmp[i][j][k] =
				0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);

	/* finally, raise first index */
	for (i = 0; i < NDIM; i++)
		for (j = 0; j < NDIM; j++)
			for (k = 0; k < NDIM; k++) {
				conn[i][j][k] = 0.;
				for (l = 0; l < NDIM; l++)
					conn[i][j][k] += geom->gcon[i][l] * tmp[l][j][k];
			}

			/* done! */
}
#undef GAMMIEDERIVATIVE
#undef NUMREC
#undef DERTYPE
#undef CONNDELTA



// returns MCOORD gcov value for i,j element
// excessive to compute other elements, but ok for now
FTYPE gcov_func_mcoord(FTYPE* X, int i, int j)
{
  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE gcovpert[NDIM];
  void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);

  gcov_func(1,MCOORD,X,gcovmcoord,gcovpert);
  //  if(i==j) return(gcovpert[i]); // for accurate non-rel gravity (works since d(1+g_{tt}) = dg_{tt})
  //  else
     return(gcovmcoord[i][j]);
}



/*

Based upon EOM:

(f_{(\nu)} T^t_\nu)_{,t} 
= -(f_{(\nu)}T^j_\nu)_{,j}
+  f_{(\nu)} (T^\lambda_\nu[\ln(f_{(\nu)}/\detg)]_{,\lambda}
+T^\mu_\lambda \Gamma^\lambda_{\mu\nu}
+\ln(f_{(\nu)})_{,t} T^t_\nu
)
*/

// returns MCOORD  value for log(gdet).  Doesn't use i,j (these are not grid locations)
FTYPE lngdet_func_mcoord(FTYPE* X, int i, int j)
{
  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE gcovpertcoord[NDIM];
  FTYPE toreturn;
  FTYPE eomfunc;
  void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);

  gcov_func(1,MCOORD,X,gcovmcoord,gcovpertcoord);
  eomfunc_func(1,MCOORD,X,&eomfunc);

  toreturn=log(eomfunc/gdet_func(gcovmcoord));

  return(toreturn);
}






// jon's version of NR's dfridr modified to accept more general, needed, function
#if(1) // uniformly low values although not always lower than original version
#define MAXITER 50
#define NRANSI
#define CON 1.1
#define CON2 (CON*CON)
//#define NTAB 130 // number of function evaluations is 2XNTAB
#define NTAB 30 // number of function evaluations is 2XNTAB
#define SAFE 2.0
//#define NRHMAX 1 // maximum size where metric changes substantially
#define NRHMAX 1.e-1
#define TRYTOL (pow(NUMEPSILON,2.0/3.0)) // attempted tolerance
#define OKTOL 1e-1 // error must be below this
#define MINH (pow(NUMEPSILON,2.0/3.0)) // minimum h
//#define HSTARTCHANGE 10.0
#define HSTARTCHANGE 2.0
#endif

#if(0) // original version (gets pretty damn low for many, but not all, derivatives -- some 1E-6 rather than 1E-13)
#define MAXITER 128
#define NRANSI
#define CON 1.3
#define CON2 (CON*CON)
#define NTAB 10 // number of function evaluations is 2XNTAB
#define SAFE 2.0
//#define NRHMAX 1 // maximum size where metric changes substantially
#define NRHMAX 1E-2
#define TRYTOL 1E-10 // attempted tolerance
#define OKTOL 1e-5 // error must be below this
#define MINH 1E-15 // minimum h
#define HSTARTCHANGE 10.0
#endif

// whether to turn on extensive recent debugging
#define DEBUGDF 0

FTYPE dfridr(FTYPE (*func)(FTYPE*,int,int), FTYPE *X,int ii, int jj, int kk)
{
  int i,j,k;
  FTYPE errt,fac,hh,**a,ans;
  FTYPE dX[NDIM],Xh[NDIM],Xl[NDIM];
  FTYPE h,err;
  FTYPE hstart;
  FTYPE newdx,temp;
  int shit;
  int iter;
  FTYPE errlist[MAXITER];
  FTYPE hhlist[MAXITER];
  FTYPE anslist[MAXITER];
  FTYPE minerror,minerrorhstart;
  FTYPE hstartfrac;
  int miniter;
  int iterdebug;
  FTYPE minans;


  hstartfrac=HSTARTCHANGE; // starting fractional drop in hstart
	
  // allocate memory
  a=dmatrix(1,NTAB,1,NTAB);

  //    hstart=NRHMAX;
  //  hstart=MAX(dx[kk],NRHMAX);
  //  hstart=MIN(NRHMAX,hstart*10.0);
  hstart=1.0;



  miniter=0; // didn't find minimum is assumed, which means first is minimum!
  iter=0;
  minerror=BIG;
  minans=BIG;

  while(1){

#if(DEBUGDF)
    dualfprintf(fail_file,"iter=%d\n",iter);
#endif


    h=hstart;
    if (h == 0.0) nrerror("h must be nonzero in dfridr.");

    hh=h;
    // HARM STUFF
    for(k=0;k<NDIM;k++) dX[k]=0.0; // other components will remains 0 for this function
    dX[kk]=hh;
    for(k=0;k<NDIM;k++){
      //      Xl[k]=X[k]-dX[k];
      temp=X[k]-dX[k];
      newdx=temp-X[k];
      Xl[k]=X[k]+newdx;
    }
    for(k=0;k<NDIM;k++){
      //      Xh[k]=X[k]+dX[k];
      temp=X[k]+dX[k];
      newdx=temp-X[k];
      Xh[k]=X[k]+newdx;
    }
    // end HARM STUFF
    //    for(k=0;k<NDIM;k++) dualfprintf(fail_file,"k=%d %21.15g %21.15g\n",k,X[k],dX[k]);
    a[1][1]=((*func)(Xh,ii,jj)-(*func)(Xl,ii,jj))/(2.0*hh);
    err=BIG;

    for (i=2;i<=NTAB;i++) {
      hh /= CON;
      // HARM STUFF
      dX[kk]=hh;
      for(k=0;k<NDIM;k++) Xl[k]=X[k]-dX[k];
      for(k=0;k<NDIM;k++) Xh[k]=X[k]+dX[k];
      //      for(k=0;k<NDIM;k++) dualfprintf(fail_file,"i=%d k=%d %21.15g %21.15g\n",i,k,X[k],dX[k]);
      // end HARM STUFF
      a[1][i]=((*func)(Xh,ii,jj)-(*func)(Xl,ii,jj))/(2.0*hh);
      fac=CON2;

#if(DEBUGDF)
	// debug
      //	if(ii==1 && jj==1 && kk==1){
      //	  if( icurr==29 && jcurr==-11){
	    for(shit=0;shit<4;shit++) dualfprintf(fail_file,"hh=%21.15g Xl[%d]=%21.15g Xh[%d]=%21.15g\n",hh,shit,Xl[shit],shit,Xh[shit]);
	    //	  }
    //	}
#endif


      for (j=2;j<=i;j++) {
	a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
	fac=CON2*fac;
	//normalized error
	//	errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]))/((*func)(X,ii,jj)+SMALL);
	//	errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
	// normalized error
	errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]))/(fabs(a[j][i])+SMALL);

#if(DEBUGDF)
	// debug
	//	if(ii==3 && jj==3 && kk==2){
	//	if(ii==2 && jj==2 && kk==2){
	//	if(ii==1 && jj==1 && kk==1){
	//	  if( icurr==29 && jcurr==-11){
	    //	  if( icurr==40 && ((jcurr==54)||(jcurr==9)) ){
	    dualfprintf(fail_file,"a11=%21.15g a[%d][%d]=%21.15g errt=%21.15g fac=%21.15g\n",a[1][1],i,j,a[i][j],errt,fac);
	    //	  }
	    //	}
#endif

	if (errt <= err) {
	  err=errt;
	  ans=a[j][i];
	}
      }
      // normalized error
      //      if (fabs((a[i][i]-a[i-1][i-1])/( (*func)(X,ii,jj)+SMALL)) >= SAFE*(err)) break;
      if (fabs((a[i][i]-a[i-1][i-1])) >= SAFE*(err)) break;
      if (fabs((a[i][i]-a[i-1][i-1]))/(fabs(ans)+SMALL) >= SAFE*(err)) break;
    }
	  
    // now check error is not crazy with the starting large h, decrease h if crazy until not crazy and get good error
    errlist[iter]=err; // store error
    //    hhlist[iter]=hstart;
    hhlist[iter]=hh; // current hh to focus around
    anslist[iter]=ans;

    if(err>TRYTOL){
      if(hstart<MINH){
	if(err<OKTOL) break;
	else{

	  if(miniter>=0){
	    if(minerror<OKTOL){
	      ans=minans;
	      break;	      // then accept as answer
	    }
	  }
	  else{ // then must fail
	    dualfprintf(fail_file,"never found error below %21.15g: err=%21.15g : ii=%d jj=%d kk=%d\n",OKTOL,err,ii,jj,kk);
	    dualfprintf(fail_file,"gi=%d gj=%d gk=%d\n",icurr,jcurr,kcurr);
	    dualfprintf(fail_file,"miniter=%d errlist[miniter]=%21.15g hhlist[miniter]=%21.15g\n",miniter,errlist[miniter],hhlist[miniter]);	    
	    dualfprintf(fail_file,"minerror=%21.15g minans=%21.15g\n",minerror,minans);
	  
	    for(iterdebug=0;iterdebug<=iter;iterdebug++){
	      dualfprintf(fail_file,"h[%d]=%21.15g err[%d]=%21.15g ans[%d]=%21.15g\n",iterdebug,hhlist[iterdebug],iterdebug,errlist[iterdebug],iterdebug,anslist[iterdebug]);
	    }
	    myexit(66);
	  }
	}
      }
      else{ // keep going since not below MINH
	if(errlist[iter]<minerror){
	  // store min error event
	  minerror=errlist[iter];
	  minerrorhstart=hhlist[iter];
	  miniter=iter;
	  minans=ans;
#if(DEBUGDF)
	  dualfprintf(fail_file,"minerr=%21.15g minhstart=%21.15g miniter=%d minans=%21.15g\n",minerror,minerrorhstart,miniter,minans);
#endif
	}
	if(iter>0){
	  // check on error
	  if(0&& (miniter>=1) && (errlist[iter]>1000.0*minerror)){
	    // then something bad is happening ... probably skipped over minimum error region
	    hstartfrac=(hstartfrac-1.0)*0.5+1.0; // slow down
	    //	    hstart=hhlist[miniter]*(1.0+NTAB*0.5/hstartfrac); // back up a bit
	    hstart=minerrorhstart*10.0;
	    iter=-1; // start over
	    //	    miniter=0;
#if(DEBUGDF||1)
	    dualfprintf(fail_file,"hstart=%21.15g hstartfrac=%21.15g\n",hstart,hstartfrac);
#endif
	  }
	  else{
	    // keep dropping h
	    hstart/=hstartfrac;
#if(DEBUGDF)
	    dualfprintf(fail_file,"hstartnew=%21.15g\n",hstart);
#endif
	  }
	}
      } // end of keep going part
    }// end if err > TRYTOL
    else break; // then error<TRYTOL!  GOOD!
    iter++;
    if(iter>=MAXITER){
      if(err<OKTOL) break;
      else{ // then maybe problems
	if(minerror<OKTOL){
	  ans=minans;
	  break;	      // then accept as answer
	}
	else{
	  // then must fail
	  dualfprintf(fail_file,"iter>MAXITER: never found error below %21.15g: err=%21.15g : ii=%d jj=%d kk=%d\n",OKTOL,err,ii,jj,kk);
	  dualfprintf(fail_file,"gi=%d gj=%d gk=%d\n",icurr,jcurr,kcurr);
	  dualfprintf(fail_file,"miniter=%d errlist[miniter]=%21.15g hhlist[miniter]=%21.15g\n",miniter,errlist[miniter],hhlist[miniter]);
	  dualfprintf(fail_file,"minerror=%21.15g minans=%21.15g\n",minerror,minans);
	  for(iterdebug=0;iterdebug<iter;iterdebug++){
	    dualfprintf(fail_file,"h[%d]=%21.15g err[%d]=%21.15g ans[%d]=%21.15g\n",iterdebug,hhlist[iterdebug],iterdebug,errlist[iterdebug],iterdebug,anslist[iterdebug]);
	  }
	  myexit(67);
	}
      }      
    }
  }
  // done
  free_dmatrix(a,1,NTAB,1,NTAB);
  //  dualfprintf(fail_file,"hfinal=%21.15g errfinal=%21.15g\n",hh,err);
  //fflush(fail_file);
  return ans;

       

}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE
#undef NRANSI
// (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. 






/* 
   FTYPE delta(int i, int j) { if(i == j) return(1.) ; else return(0.) 
   ; } */

/* Minkowski metric; signature +2 */
/* 
   FTYPE mink(int i, int j) { if(i == j) { if(i == 0) return(-1.) ;
   else return(1.) ; } else return(0.) ; } */




//////////////////////////////
//
// below very specific to defcoord==LOGRSINTH and MCOORD=KSCOORDS:
// gives: analytic connection and source
//
/////////////////////////////////


// below 2 mysin/mycos preserve symmetry for arg=-M_PI/2 .. 2M_PI by restricting arg to 0..PI/2
FTYPE mysin(FTYPE th)
{
  FTYPE arg;

#if(ACCURATESINCOS)
  if(th<0.0){
    arg=-th;
    return(-sin(arg));
  } 
  if(th<=0.5*M_PI){
    // base
    arg=th;
    return(sin(arg));
  }
  else if(th>0.5*M_PI && th<=M_PI){
    arg=M_PI-th;
    return(sin(arg));
  }
  else if(th>M_PI && th <=1.5*M_PI){
    arg=-(M_PI-th);
    return(-sin(arg));
  }
  else if(th>1.5*M_PI && th <=2.0*M_PI){
    arg=2.0*M_PI-th;
    return(-sin(arg));
  }
  else if(th>2.0*M_PI && th <=2.5*M_PI){
    arg=-(2.0*M_PI-th);
    return(sin(arg));
  }
  else{
    arg=-(2.0*M_PI-th);
    return(sin(arg));
  }
  //  if(arg<0.0 || arg>M_PI*0.5){
  //    dualfprintf(fail_file,"Problem with mysin\n");
  //  }

#else
  return(sin(th));
#endif

}


FTYPE mycos(FTYPE th)
{
  FTYPE arg;


#if(ACCURATESINCOS)
  if(th<0.0){
    arg=-th;
    return(cos(arg));
  } 
  if(th<=0.5*M_PI){
    // base
    arg=th;
    return(cos(arg));
  }
  else if(th>0.5*M_PI && th<=M_PI){
    arg=M_PI-th;
    return(-cos(arg));
  }
  else if(th>M_PI && th <=1.5*M_PI){
    arg=-(M_PI-th);
    return(-cos(arg));
  }
  else if(th>1.5*M_PI && th <=2.0*M_PI){
    arg=2.0*M_PI-th;
    return(cos(arg));
  }
  else if(th>2.0*M_PI && th <=2.5*M_PI){
    arg=-(2.0*M_PI-th);
    return(-cos(arg));
  }
  else{
    arg=-(2.0*M_PI-th);
    return(-cos(arg));
  }
  //  if(arg<0.0 || arg>M_PI*0.5){
  //    dualfprintf(fail_file,"Problem with mycos\n");
  //  }
#else
  return(cos(th));
#endif

}



//#ifdef WIN32
FTYPE cot(FTYPE arg)
{
  return(1.0/tan(arg));
}
//#endif 

FTYPE csc(FTYPE arg)
{
  return(1.0/sin(arg));
}

// jon's MKS connection (and conn2)
// only applies to defcoord==LOGRSINTH
void mks_conn_func(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  int i, j, k, l;
  FTYPE V[NDIM];
  FTYPE r,th,sigma,dxdxptrue[NDIM][NDIM];
#ifdef WIN32
  FTYPE cot(FTYPE arg);
#endif
  FTYPE dxdxp[NDIM];

  // get bl coordinates
  bl_coord(X,V);
  r=V[1];
  th=V[2];


  // the connection

  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }

  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  // set aux vars
  dxdxprim(X,V,dxdxptrue);
  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal
  sigma=r*r+a*a*cos(th)*cos(th);


  conn[0][0][0]=(-2.*r*sigma + 4.*pow(r,3.))*pow(sigma,-3.);
conn[0][0][1]=dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
conn[0][0][2]=-1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th);
conn[0][0][3]=-2.*a*r*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[0][1][0]=dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
conn[0][1][1]=2.*(r + sigma)*pow(dxdxp[1],2.)*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.);
conn[0][1][2]=-1.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
    sin(2.*th);
conn[0][1][3]=dxdxp[1]*a*(2.*r + sigma)*(sigma - 2.*pow(r,2.))*
    pow(sigma,-3.)*pow(sin(th),2.);
conn[0][2][0]=-1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th);
conn[0][2][1]=-1.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
    sin(2.*th);
conn[0][2][2]=-2.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-1.);
conn[0][2][3]=2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-2.)*
    pow(sin(th),3.);
conn[0][3][0]=-2.*a*r*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[0][3][1]=dxdxp[1]*a*(2.*r + sigma)*(sigma - 2.*pow(r,2.))*
    pow(sigma,-3.)*pow(sin(th),2.);
conn[0][3][2]=2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-2.)*
    pow(sin(th),3.);
conn[0][3][3]=2.*r*pow(sigma,-3.)*pow(sin(th),2.)*
    (-1.*r*pow(sigma,2.) + pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*
       pow(sin(th),2.));
conn[1][0][0]=pow(dxdxp[1],-1.)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.)*(-2.*r + sigma + pow(a,2.)*pow(sin(th),2.));
conn[1][0][1]=0.5*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (sigma - 2.*pow(r,2.))*pow(sigma,-3.);
conn[1][0][2]=0.;
conn[1][0][3]=0.5*a*pow(dxdxp[1],-1.)*
    (4.*r - 2.*sigma - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*pow(sin(th),2.);
conn[1][1][0]=0.5*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (sigma - 2.*pow(r,2.))*pow(sigma,-3.);
conn[1][1][1]=pow(sigma,-3.)*
    (-1.*dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.)) + 
      pow(sigma,3.) + dxdxp[1]*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*
       pow(sin(th),2.));
conn[1][1][2]=-1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)
   ;
conn[1][1][3]=0.5*a*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
      cos(2.*th)*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.)) + 
      2.*r*((-2. + sigma)*sigma + 4.*pow(r,2.)))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[1][2][0]=0.;
conn[1][2][1]=-1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)
   ;
conn[1][2][2]=-1.*r*pow(dxdxp[1],-1.)*pow(dxdxp[2],2.)*
    pow(sigma,-1.)*(-2.*r + sigma + pow(a,2.)*pow(sin(th),2.));
conn[1][2][3]=0.;
conn[1][3][0]=0.5*a*pow(dxdxp[1],-1.)*
    (4.*r - 2.*sigma - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*pow(sin(th),2.);
conn[1][3][1]=0.5*a*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
      cos(2.*th)*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.)) + 
      2.*r*((-2. + sigma)*sigma + 4.*pow(r,2.)))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[1][3][2]=0.;
conn[1][3][3]=-1.*pow(dxdxp[1],-1.)*pow(sigma,-3.)*pow(sin(th),2.)*
    (-2.*r + sigma + pow(a,2.)*pow(sin(th),2.))*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
conn[2][0][0]=-1.*r*pow(dxdxp[2],-1.)*pow(a,2.)*pow(sigma,-3.)*
    sin(2.*th);
conn[2][0][1]=-1.*dxdxp[1]*r*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
conn[2][0][2]=0.;
conn[2][0][3]=2.*a*r*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (sigma + pow(a,2.)*pow(sin(th),2.))*sin(th);
conn[2][1][0]=-1.*dxdxp[1]*r*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
conn[2][1][1]=-1.*r*pow(dxdxp[1],2.)*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
conn[2][1][2]=dxdxp[1]*r*pow(sigma,-1.);
conn[2][1][3]=dxdxp[1]*a*pow(dxdxp[2],-1.)*pow(sigma,-3.)*sin(th)*
    (sigma*(2.*r + sigma)*cos(th) + r*pow(a,2.)*sin(th)*sin(2.*th));
conn[2][2][0]=0.;
conn[2][2][1]=dxdxp[1]*r*pow(sigma,-1.);
conn[2][2][2]=4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*
     pow(M_PI,2.) - 1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)\
    ;
conn[2][2][3]=0.;
conn[2][3][0]=2.*a*r*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (sigma + pow(a,2.)*pow(sin(th),2.))*sin(th);
conn[2][3][1]=dxdxp[1]*a*pow(dxdxp[2],-1.)*pow(sigma,-3.)*sin(th)*
    (sigma*(2.*r + sigma)*cos(th) + r*pow(a,2.)*sin(th)*sin(2.*th));
conn[2][3][2]=0.;
conn[2][3][3]=-1.*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (pow(sigma,3.) + sigma*(4.*r + sigma)*pow(a,2.)*pow(sin(th),2.) + 
      2.*r*pow(a,4.)*pow(sin(th),4.))*sin(th);
conn[3][0][0]=a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.);
conn[3][0][1]=dxdxp[1]*a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)
   ;
conn[3][0][2]=-2.*dxdxp[2]*a*r*cot(th)*pow(sigma,-2.);
conn[3][0][3]=-1.*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[3][1][0]=dxdxp[1]*a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)
   ;
conn[3][1][1]=a*pow(dxdxp[1],2.)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
conn[3][1][2]=-1.*dxdxp[1]*dxdxp[2]*a*(2.*r + sigma)*cot(th)*
    pow(sigma,-2.);
conn[3][1][3]=dxdxp[1]*pow(sigma,-3.)*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
conn[3][2][0]=-2.*dxdxp[2]*a*r*cot(th)*pow(sigma,-2.);
conn[3][2][1]=-1.*dxdxp[1]*dxdxp[2]*a*(2.*r + sigma)*cot(th)*
    pow(sigma,-2.);
conn[3][2][2]=-1.*a*r*pow(dxdxp[2],2.)*pow(sigma,-1.);
conn[3][2][3]=dxdxp[2]*
    (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th));
conn[3][3][0]=-1.*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[3][3][1]=dxdxp[1]*pow(sigma,-3.)*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
conn[3][3][2]=dxdxp[2]*
    (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th));
conn[3][3][3]=pow(sigma,-3.)*
    (-1.*a*r*pow(sigma,2.)*pow(sin(th),2.) + 
      pow(a,3.)*(-1.*sigma + 2.*pow(r,2.))*pow(sin(th),4.));
conn2[0]=0.;
conn2[1]=-1.*pow(sigma,-1.)*(2.*dxdxp[1]*r + pow(r,2.) + 
      pow(a,2.)*pow(cos(th),2.));
conn2[2]=-1.*dxdxp[2]*cot(th) + 
    4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
    dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th);
conn2[3]=0.;


}



/*********************************************************************************************
   Scott's s MKS connection that can be used for any transformation between r,th <-> X1,X2 :
*********************************************************************************************/
#define M (1.)
void mks_conn_func_general(FTYPE *X, struct of_geom *geom, FTYPE conn[][NDIM][NDIM] )
{
  int i, j, k, l;
  FTYPE V[NDIM];
  FTYPE r,th,sigma,dxdxp[NDIM][NDIM],dxdxp_dxp[NDIM][NDIM][NDIM];

  double t1,   t10,   t102,   t1024,   t1035,   t1037,   t104,   t11,   t114,   t116,   t119;
  double  t12,   t121,   t123,   t126,   t129,   t130,   t132,   t14,   t148,   t149,   t15,   t152;
  double t154,   t156,   t157,   t158,   t159,   t161,   t169,   t17,   t171,   t172,   t175,   t177;
  double t185,   t2,   t203,   t204,   t208,   t209,   t21,   t210,   t212,   t214,   t22,   t221;
  double   t222,   t224,   t227,   t23,   t236,   t24,   t240,   t241,   t242,   t243,   t245,   t246;
  double    t247,   t248,   t25,   t250,   t251,   t258,   t26,   t260,   t261,   t263,   t264,   t271;
  double   t273,   t275,   t276,   t278,   t28,   t280,   t281,   t283,   t284,   t285,   t286,   t288;
  double    t289,   t29,   t297,   t299,   t3,   t30,   t300,   t303,   t305,   t306,   t308,   t309;
  double    t31,   t310,   t313,   t314,   t320,   t325,   t327,   t328,   t329,   t330,   t333,   t336;
  double    t338,   t34,   t340,   t342,   t344,   t346,   t35,   t358,   t361,   t363,   t366,   t367;
  double    t368,   t370,   t372,   t375,   t38,   t380,   t381,   t384,   t385,   t387,   t39,   t399;
  double    t4,   t40,   t400,   t402,   t404,   t405,   t406,   t408,   t409,   t41,   t411,   t412;
  double    t418,   t42,   t421,   t425,   t428,   t431,   t432,   t433,   t434,   t437,   t440,   t442;
  double    t448,   t451,   t453,   t454,   t459,   t462,   t467,   t469,   t480,   t481,   t486,   t487;
  double    t488,   t491,   t492,   t498,   t501,   t504,   t507,   t508,   t510,   t512,   t52,   t521;
  double    t528,   t53,   t530,   t553,   t556,   t56,   t57,   t588,   t60,   t607,   t627,   t628;
  double    t63,   t630,   t631,   t632,   t634,   t636,   t637,   t64,   t651,   t652,   t654,   t656;
  double    t657,   t659,   t661,   t662,   t670,   t673,   t675,   t677,   t686,   t689,   t7,   t712;
  double    t74,   t748,   t75,   t78,   t793,   t794,   t795,   t799,   t8,   t800,   t801,   t803;
  double    t806,  t807,   t813,   t816,   t822,   t83,   t831,   t84,   t845,   t86,   t89,   t891; 
  double    t90,   t91,   t916,   t917,   t920,   t924,   t928,   t940,   t946,   t968,   t97, t970,t991;

  // get bl coordinates
  bl_coord(X,V);
  r=V[1];
  th=V[2];

  // the connection

  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }
  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  // set aux vars
  dxdxprim(X,V,dxdxp);
  //  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal
  //  dx_dxp_calc(X,r,th,dx_dxp);
  //dx_dxp_dxp_calc(X,r,th,dx_dxp_dxp);



  // GODMARK
  // need to set second derivative analytically
  for(i=0;i<NDIM;i++)  for(j=0;j<NDIM;j++)  for(k=0;k<NDIM;k++){
    dxdxp_dxp[i][j][k]=0.0;
  }


  //  t1 = rf(X1,X2);
  t1 = r;
  //  t2 = thf(X1,X2);
  t2 = th;
  t3 = cos(t2);
  t4 = a*t3;
  t7 = (-t1+t4)*(t1+t4);
  t8 = M*M;
  t10 = t1*t1;
  t11 = a*a;
  t12 = t3*t3;
  t14 = t10+t11*t12;
  t15 = t14*t14;
  t17 = 1/t15/t14;
  conn[0][0][0] = -2.0*t7*t8*t1*t17;
  t21 = t10*t10;
  //  t22 = diff(rf(X1,X2),X1);
  t22 = dxdxp[RR][1];
  t23 = t21*t22;
  t24 = t10*t1;
  t25 = M*t24;
  t26 = t25*t22;
  t28 = t24*t3;
  t29 = sin(t2);
  //  t30 = diff(thf(X1,X2),X1);
  t30 = dxdxp[TH][1];
  t31 = t29*t30;
  t34 = M*t1;
  t35 = t22*t12;
  t38 = t12*t12;
  t39 = t38*t22;
  t40 = t29*t1;
  t41 = t12*t3;
  t42 = t41*t30;
  conn[0][0][1] = -M*(-t23-2.0*t26+(2.0*t28*t31+2.0*t34*t35+(t39+2.0*t40*t42)*t11)*t11)*t17;
  //  t52 = diff(rf(X1,X2),X2);
  t52 = dxdxp[RR][2];
  t53 = t21*t52;
  //  t56 = diff(thf(X1,X2),X2);
  t56 = dxdxp[TH][2];
  t57 = t29*t56;
  t60 = t52*t12;
  t63 = t38*t52;
  t64 = t41*t1;
  conn[0][0][2] = -M*(-t53-2.0*t25*t52+(2.0*t28*t57+2.0*t34*t60+(t63+2.0*t64*t57)*t11)*t11)*t17;
  t74 = -1.0+t3;
  t75 = 1.0+t3;
  t78 = t7*t74*t75;
  conn[0][0][3] = -2.0*t78*a*t1*t8*t17;
  conn[0][1][0] = conn[0][0][1];
  t83 = t30*t30;
  t84 = t21*t10;
  t86 = t22*t22;
  t89 = t24*t22;
  t90 = t3*t29;
  t91 = t90*t30;
  t97 = t86*t12;
  t102 = t22*t29;
  t104 = t64*t102;
  conn[0][1][1] = -2.0*(t83*t84-t21*t86-t25*t86+(2.0*t89*t91+2.0*t83*t21*t12+
						 t34*t97+(t83*t10*t38+t38*t86+2.0*t104*t30)*t11)*t11)*M*t17;
  t114 = t22*t52;
  t116 = t30*t56;
  t119 = t24*t52;
  t121 = t114*t12;
  t123 = t21*t12;
  t126 = t90*t56;
  t129 = t52*t29;
  t130 = t129*t30;
  t132 = t10*t38;
  conn[0][1][2] = -2.0*(-t25*t114+t116*t84-t23*t52+(t119*t91+t34*t121+2.0*
						    t116*t123+t89*t126
						    +(t39*t52+t64*t130+t116*t132+t104*t56)*t11)*t11)*M*t17;
  t148 = 2.0*t28*t30;
  t149 = t102*t12;
  t152 = t41*t24;
  t154 = 2.0*t152*t30;
  t156 = 2.0*t64*t30;
  t157 = t39*t29;
  t158 = t30*t1;
  t159 = t38*t3;
  t161 = 2.0*t158*t159;
  t169 = a*t29*t17;
  conn[0][1][3] = -(2.0*t25*t102+t23*t29+(-t148-2.0*t34*t149+t154+(-t156-t157+t161)*t11)*t11)*M*t169;
  conn[0][2][0] = conn[0][0][2];
  conn[0][2][1] = conn[0][1][2];
  t171 = t52*t52;
  t172 = t171*M;
  t175 = t56*t56;
  t177 = t1*t12;
  t185 = t52*t41;
  conn[0][2][2] = -2.0*(-t172*t24-t171*t21+t175*t84+(t172*t177+2.0*t119*t126+
						     2.0*t175*t21*t12
						     +(t171*t38+2.0*t185*t40*t56+t175*t10*t38)*t11)*t11)*M*t17;
  t203 = 2.0*t152*t56;
  t204 = t129*t12;
  t208 = 2.0*t28*t56;
  t209 = t63*t29;
  t210 = t56*t1;
  t212 = 2.0*t210*t159;
  t214 = 2.0*t64*t56;
  conn[0][2][3] = (-2.0*t25*t129-t53*t29+(-t203+2.0*t34*t204+t208+(t209-t212+t214)*t11)*t11)*M*t169;
  conn[0][3][0] = conn[0][0][3];
  conn[0][3][1] = conn[0][1][3];
  conn[0][3][2] = conn[0][2][3];
  t221 = t21*t1;
  t222 = t24*t12;
  t224 = t10*t12;
  t227 = t1*t38;
  t236 = (-t221+(-2.0*t222+(-t224+t10)*M+(-t227+(t38-t12)*M)*t11)*t11)*t74*t75;
  conn[0][3][3] = -2.0*t236*t34*t17;
  t240 = t56*M;
  t241 = t240*t24;
  t242 = 2.0*t241;
  t243 = t56*t21;
  t245 = 2.0*t240*t177;
  t246 = t56*t10;
  t247 = t246*t12;
  t248 = t1*t3;
  t250 = 2.0*t248*t129;
  t251 = t56*t12;
  t258 = t30*t52;
  t260 = 1/(-t22*t56+t258);
  t261 = t260*t17;
  conn[1][0][0] = -(-t242+t243+(t245-t247+t250+t246-t251*t11)*t11)*M*t261;
  t263 = M*t22;
  t264 = t56*t38;
  t271 = (-t242+(t245-t247+t250+t246+(-t251+t264)*t11)*t11)*t260*t17;
  conn[1][0][1] = -t263*t271;
  t273 = M*t52;
  conn[1][0][2] = -t273*t271;
  t275 = t24*t29;
  t276 = t240*t275;
  t278 = t119*t3;
  t280 = t243*t29;
  t281 = t40*t12;
  t283 = 2.0*t240*t281;
  t284 = t29*t12;
  t285 = t246*t284;
  t286 = t52*t3;
  t288 = 2.0*t286*t1;
  t289 = t57*t10;
  t297 = a*t29*t260*t17;
  conn[1][0][3] = (-2.0*t276+2.0*t278+t280+(t283-t285+t288+t289-t56*t11*t284)*t11)*M*t297;
  conn[1][1][0] = conn[1][0][1];
  //  t299 = diff(diff(thf(X1,X2),X1),X1);
  t299 = dxdxp_dxp[TH][1][1];
  t300 = t52*t299;
  t303 = t258*t221*t22;
  t305 = t56*t84;
  //  t306 = diff(diff(rf(X1,X2),X1),X1);
  t306 = dxdxp_dxp[RR][1][1] ;
  t308 = t83*t56;
  t309 = t21*t24;
  t310 = t308*t309;
  t313 = 2.0*t308*t84;
  t314 = t56*t24;
  t320 = t30*t22;
  t325 = t258*t89*t12;
  t327 = t308*t221;
  t328 = t52*t83;
  t329 = t90*t21;
  t330 = t328*t329;
  t333 = t306*t12;
  t336 = t221*t12;
  t338 = 2.0*t308*t336;
  t340 = 4.0*t308*t123;
  t342 = t248*t29;
  t344 = 2.0*t52*t86*t342;
  t346 = t56*t86;
  t358 = t306*t38;
  t361 = t1*t22;
  t363 = t258*t361*t38;
  t366 = 2.0*t308*t222;
  t367 = t24*t38;
  t368 = t308*t367;
  t370 = t41*t29*t10;
  t372 = 2.0*t328*t370;
  t375 = 2.0*t308*t132;
  t380 = t159*t29;
  t381 = t380*t56;
  t384 = t308*t227;
  t385 = t38*t12;
  t387 = t328*t380;
  conn[1][1][1] = -(-t300*t84-2.0*t303+t305*t306-t310+(-t243*t86+t313-2.0*t314*t86*M)*M
		    +(-2.0*t320*t3*t280-4.0*t325-t327+t330-3.0*t300*t123+3.0*t243*t333
		      -t338+(t340+t344-t246*t97+t346*t10+2.0*t210*t97*M)*M
		      +(-4.0*t320*t41*t289-3.0*t300*t132+3.0*t246*t358-2.0*t363-t366-t368+t372
			+(-t346*t12+t375+2.0*t346*t38)*M
			+(-2.0*t320*t381-t384-t300*t385+t387+t56*t306*t385)*t11)*t11)*t11)*t260*t17;
  //  t399 = diff(diff(thf(X1,X2),X1),X2);
  t399 = dxdxp_dxp[TH][1][2]; 
  t400 = t52*t399;
  //  t402 = diff(diff(rf(X1,X2),X1),X2);
  t402 = dxdxp_dxp[RR][1][2]; 
  t404 = t30*t175;
  t405 = t404*t309;
  t406 = t171*t30;
  t408 = t56*t221;
  t409 = t114*t408;
  t411 = 2.0*t404*t84;
  t412 = t52*t56;
  t418 = t402*t12;
  t421 = t404*t221;
  t425 = t114*t314*t12;
  t428 = 2.0*t404*t336;
  t431 = t22*t175;
  t432 = t431*t329;
  t433 = t10*t22;
  t434 = t433*t12;
  t437 = 4.0*t404*t123;
  t440 = 2.0*t171*t22*t342;
  t442 = t210*M;
  t448 = 2.0*t370*t431;
  t451 = t114*t210*t38;
  t453 = 2.0*t404*t222;
  t454 = t402*t38;
  t459 = t404*t367;
  t462 = 2.0*t404*t132;
  t467 = t404*t227;
  t469 = t431*t380;
  conn[1][1][2] = (t400*t84-t305*t402+t405+t406*t221+t409+(-t411+t412*t23+2.0*t114*t241)*M
		   +(-3.0*t243*t418+t421+3.0*t400*t123+2.0*t425+t428+2.0*t406*t222+
		     t432+(t412*t434-t437-t440-t412*t433-2.0*t121*t442)*M
		     +(t448+t406*t227+t451+t453-3.0*t246*t454+3.0*t400*t132+t459
		       +(t114*t251-t462-2.0*t114*t264)*M+(t467+t400*t385+t469
							  -t56*t402*t385)*t11)*t11)*t11)*t260*t17;
  t480 = t286*t21;
  t481 = t57*t221;
  t486 = 2.0*t185*t10;
  t487 = t57*t222;
  t488 = 2.0*t487;
  t491 = t57*t227;
  t492 = t52*t159;
  t498 = (-t491+t492+(-t57*t12+t57*t38)*M)*t11;
  t501 = t480-t481+(2.0*t278-2.0*t276)*M+(t486-t488+(-t285+t289+t288+t283)*M+t498)*t11;
  conn[1][1][3] = t501*t22*t297;
  conn[1][2][0] = conn[1][0][2];
  conn[1][2][1] = conn[1][1][2];
  t504 = t171*t56;
  //  t507 = diff(diff(thf(X1,X2),X2),X2);
  t507 = dxdxp_dxp[TH][2][2];
  t508 = t52*t507;
  //  t510 = diff(diff(rf(X1,X2),X2),X2);
  t510 = dxdxp_dxp[RR][2][2];
  t512 = t175*t56;
  t521 = t512*t221;
  t528 = t52*t175;
  t530 = t510*t12;
  t553 = t510*t38;
  t556 = t512*t24;
  conn[1][2][2] = -(-2.0*t504*t221-t508*t84+t305*t510-t512*t309
		    +(2.0*t512*t84-t504*t21-2.0*t504*t25)*M
		    +(-2.0*t521*t12-3.0*t508*t123-4.0*t504*t222-t521-t528*t329+3.0*t243*t530
		      +(-t504*t224+t504*t10+2.0*t342*t171*t52+4.0*t512*t21*t12+2.0*t12*t171*t442)*M
		      +(-3.0*t508*t132-2.0*t504*t227-2.0*t528*t370+3.0*t246*t553-2.0*t556*t12-t556*t38
			+(2.0*t512*t10*t38-t504*t12+2.0*t504*t38)*M
			+(-t528*t380-t512*t1*t38-t508*t385+t56*t510*t385)*t11)*t11)*t11)*t260*t17;
  conn[1][2][3] = t501*t52*t297;
  conn[1][3][0] = conn[1][0][3];
  conn[1][3][1] = conn[1][1][3];
  conn[1][3][2] = conn[1][2][3];
  t588 = t84*t29;
  t607 = t29*t38;
  conn[1][3][3] = -(-t56*t309*t29+t286*t84+2.0*t240*t588
		    +(-t481-2.0*t408*t284+t480+2.0*t185*t21+(-4.0*t185*t24+t280+3.0*t243*t284+4.0*t278
							     +(-2.0*t314*t29+2.0*t487)*M)*M
		      +(t492*t10+t486-t314*t607-t488+(-2.0*t492*t1+t289+t288-2.0*t285+3.0*t246*t607
						      +(-2.0*t491+2.0*t210*t284)*M)*M+t498)*t11)*t11)*t29*t261;
  t627 = t30*t21;
  t628 = t30*M;
  t630 = 2.0*t628*t24;
  t631 = t30*t10;
  t632 = t631*t12;
  t634 = 2.0*t628*t177;
  t636 = 2.0*t361*t90;
  t637 = t30*t12;
  conn[2][0][0] = -(-t627+t630+(t632-t634-t631-t636+t637*t11)*t11)*M*t261;
  t651 = (-t630+(t634+t636+t631-t632+(-t637+t30*t38)*t11)*t11)*t260*t17;
  conn[2][0][1] = t263*t651;
  conn[2][0][2] = t273*t651;
  t652 = t628*t275;
  t654 = t89*t3;
  t656 = t627*t29;
  t657 = t631*t284;
  t659 = 2.0*t628*t281;
  t661 = 2.0*t361*t3;
  t662 = t31*t10;
  conn[2][0][3] = (2.0*t652-2.0*t654-t656+(t657-t659-t661-t662+t30*t11*t284)*t11)*M*t297;
  conn[2][1][0] = conn[2][0][1];
  t670 = t30*t86;
  t673 = t30*t84;
  t675 = t22*t299;
  t677 = t83*t30;
  t686 = t677*t221;
  t689 = t83*t22;
  t712 = t677*t24;
  conn[2][1][1] = -(2.0*t670*t221-t673*t306+t675*t84+t677*t309+(-2.0*t677*t84+t627*t86+2.0*t670*t25)*M
		    +(2.0*t686*t12+t689*t329-3.0*t627*t333+t686+4.0*t670*t222+3.0*t675*t123
		      +(-2.0*t86*t22*t1*t90+t631*t97-t670*t10-4.0*t677*t21*t12-2.0*t637*t86*t34)*M
		      +(2.0*t712*t12+t712*t38+2.0*t670*t227+3.0*t675*t132-3.0*t631*t358+2.0*t689*t370
			+(t670*t12-2.0*t677*t10*t38-2.0*t670*t38)*M
			+(t675*t385-t30*t306*t385+t689*t380+t677*t1*t38)*t11)*t11)*t11)*t260*t17;
  t748 = t22*t399;
  conn[2][1][2] = -(t303+t310-t673*t402+t748*t84+t346*t221
		    +(-t313+t258*t23+2.0*t258*t26)*M
		    +(t327+t330+2.0*t325-3.0*t627*t418+2.0*t346*t222+3.0*t748*t123+t338
		      +(-t340-t258*t433-t344+t258*t434-2.0*t60*t30*t361*M)*M
		      +(t346*t227+t372+t363+t366+t368-3.0*t631*t454+3.0*t748*t132
			+(t258*t35-t375-2.0*t258*t39)*M+(t387+t748*t385+t384
							 -t30*t402*t385)*t11)*t11)*t11)*t260*t17;
  t793 = t31*t221;
  t794 = t3*t22;
  t795 = t794*t21;
  t799 = t31*t222;
  t800 = 2.0*t799;
  t801 = t22*t41;
  t803 = 2.0*t801*t10;
  t806 = t31*t227;
  t807 = t22*t159;
  t813 = (-t806+t807+(t31*t38-t31*t12)*M)*t11;
  t816 = -t793+t795+(2.0*t654-2.0*t652)*M+(-t800+t803+(t661-t657+t662+t659)*M+t813)*t11;
  conn[2][1][3] = -t816*t22*t297;
  conn[2][2][0] = conn[2][0][2];
  conn[2][2][1] = conn[2][1][2];
  t822 = t22*t507;
  t831 = t3*t56;
  t845 = t41*t56;
  conn[2][2][2] = -(-t673*t510+2.0*t409+t405+t822*t84+(t406*t21-t411+2.0*t406*t25)*M
		    +(-3.0*t627*t530+t421-t432+2.0*t130*t831*t21+t428+3.0*t822*t123+4.0*t425
		      +(t406*t224-t406*t10-t440-t437-2.0*t406*t177*M)*M
		      +(4.0*t130*t845*t10+2.0*t451+3.0*t822*t132+t453+t459-t448-3.0*t631*t553
			+(t406*t12-t462-2.0*t406*t38)*M+(2.0*t258*t381+t822*t385-t30*t510*t385
							 +t467-t469)*t11)*t11)*t11)*t260*t17;
  conn[2][2][3] = -t816*t52*t297;
  conn[2][3][0] = conn[2][0][3];
  conn[2][3][1] = conn[2][1][3];
  conn[2][3][2] = conn[2][2][3];
  t891 = t30*t24;
  conn[2][3][3] = (t794*t84+2.0*t628*t588-t30*t309*t29
		   +(t795-t793-2.0*t30*t221*t284+2.0*t801*t21
		     +(3.0*t627*t284+t656-4.0*t801*t24+4.0*t654+(-2.0*t891*t29+2.0*t799)*M)*M
		     +(t807*t10+t803-t891*t607-t800+(3.0*t631*t607-2.0*t807*t1+t661-2.0*t657+t662
						     +(2.0*t158*t284-2.0*t806)*M)*M+t813)*t11)*t11)*t29*t261;
  t917 = a*M;
  conn[3][0][0] = -t7*t917*t17;
  t920 = t102*t10;
  t924 = 1/t29;
  t916 = t924*t17;
  conn[3][0][1] = -t917*(-t920+t148+(t156+t149)*t11)*t916;
  t928 = t129*t10;
  conn[3][0][2] = -t917*(t208-t928+(t214+t204)*t11)*t916;
  conn[3][0][3] = -t78*M*t11*t17;
  conn[3][1][0] = conn[3][0][1];
  t940 = t83*t29;
  t946 = t86*t29;
  t968 = t916;
  conn[3][1][1] = -(t940*t221+2.0*t794*t627+(4.0*t794*t891-t946*t10)*M
		    +(4.0*t801*t631+2.0*t940*t222+(4.0*t361*t42+t946*t12)*M+(t940*t227+2.0*t807*t30)*t11)
		    *t11)*a*t968;
  t970 = t29*t221;
  t991 = t52*t1;
  conn[3][1][2] = -(t116*t970+t286*t627+t794*t243
		    +(2.0*t286*t891-t102*t52*t10+2.0*t794*t314)*M
		    +(2.0*t185*t631+2.0*t801*t246+2.0*t116*t275*t12+(2.0*t361*t845+2.0*t991*t42+t102*t60)*M
		      +(t116*t40*t38+t492*t30+t807*t56)*t11)*t11)*a*t968;
  t1024 = t38*t41;
  conn[3][1][3] = (t3*t30*t84+t970*t22+(3.0*t42*t21+2.0*t275*t35
					+(t102*t224+t148-t154-t920)*M
					+(t40*t39+3.0*t159*t30*t10+(t156-t161+t149-t157)*M
					  +t1024*t30*t11)*t11)*t11)*t924*t17;
  conn[3][2][0] = conn[3][0][2];
  conn[3][2][1] = conn[3][1][2];
  t1035 = t175*t29;
  t1037 = t171*t29;
  conn[3][2][2] = -(2.0*t286*t243+t1035*t221+(-t1037*t10+4.0*t286*t314)*M
		    +(4.0*t185*t246+2.0*t1035*t222+(t1037*t12+4.0*t991*t845)*M
		      +(t1035*t227+2.0*t492*t56)*t11)*t11)*a*t968;
  conn[3][2][3] = -(-t970*t52-t831*t84+(-2.0*t275*t60-3.0*t845*t21
					+(t203-t208-t129*t224+t928)*M
					+(-t40*t63-3.0*t159*t56*t10+(t209+t212-t204-t214)*M
					  -t1024*t56*t11)*t11)*t11)*t924*t17;
  conn[3][3][0] = conn[3][0][3];
  conn[3][3][1] = conn[3][1][3];
  conn[3][3][2] = conn[3][2][3];
  conn[3][3][3] = -t236*a*t17;

  return;

}

#undef M






// jon's MKS source in mid-simplified form (used to UPDATE the source)
// only for ideal EOS (gam)
void mks_source_conn(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE *dU)
{
  int ii,jj,kk;
  int i=0, j=0, k=0, l=0;
  FTYPE r,th,X[NDIM],V[NDIM],sigma,dxdxptrue[NDIM][NDIM];
  //  FTYPE cot(FTYPE arg),csc(FTYPE arg);
  FTYPE b[NDIM],u[NDIM],bsq,en,rho;
  FTYPE dxdxp[NDIM];

  ii=ptrgeom->i;
  jj=ptrgeom->j;
  kk=ptrgeom->k;


  bsq = dot(q->bcon, q->bcov);
  u[TT]=q->ucon[TT];
  u[RR]=q->ucon[RR];
  u[TH]=q->ucon[TH];
  u[PH]=q->ucon[PH];

  b[TT]=q->bcon[TT];
  b[RR]=q->bcon[RR];
  b[TH]=q->bcon[TH];
  b[PH]=q->bcon[PH];

  rho=pr[RHO];
  en=pr[UU];

  coord(ptrgeom->i, ptrgeom->j, ptrgeom->k,ptrgeom->p, X);
  // get bl coordinates
  bl_coord(X,V);
  r=V[1];
  th=V[2];


  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }
  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif
  // set aux vars
  dxdxprim(X,V,dxdxptrue);
  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal


  sigma=r*r+a*a*cos(th)*cos(th);



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU0==1)){
  // see grmhd-fullsource-simplify.nb

    dU[UU]+=pow(sigma,-1.)*(-1.*(-1. - 2.*dxdxp[1]*r*pow(sigma,-1.))*
       (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
               2.*b[PH]*a*pow(sin(th),2.))) + 
         u[RR]*(bsq + en*gam + rho)*
          (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
               u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
      (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
               2.*b[PH]*a*pow(sin(th),2.))) - 
         1.*u[TH]*(bsq + en*gam + rho)*
          (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
               u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*
       (-1.*dxdxp[2]*cot(th) + 
         4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
         dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)));

  }
  // else nothing to add then


  if((WHICHEOM==WITHNOGDET)&&(NOGDETU1==1)){

// see grmhd-ksksp-mhd.nb and the other source*simplify*.nb files

dU[U1]+=0.0625*(16.*dxdxp[1]*r*
       (0.5*bsq + en*(-1. + gam) - 
         1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
         (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
       pow(sigma,-1.) + 2.*(-1. - 2.*dxdxp[1]*r*pow(sigma,-1.))*
       (4.*bsq + 8.*en*(-1. + gam) - 
         1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
            8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
             cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
            1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
            8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
          pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
          (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
            4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
            4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
            u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
            4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.)) - 
      1.*dxdxp[1]*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
       (-1.*b[TT]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
            8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
             cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
            1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
            8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
         u[TT]*(bsq + en*gam + rho)*
          (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
            4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
            4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
            u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
            4.*u[PH]*a*pow(r,2.)))*pow(sigma,-4.)*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)) + 
      8.*a*pow(dxdxp[1],2.)*(-1.*b[RR]*
          (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
            b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
               b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
            b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
         u[RR]*(bsq + en*gam + rho)*
          (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
               dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
               u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
            u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-4.)*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.) + 
      8.*dxdxp[1]*a*(-1.*b[TT]*
          (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
            b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
               b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
            b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
         u[TT]*(bsq + en*gam + rho)*
          (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
               dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
               u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
            u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-4.)*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.) + 
      dxdxp[1]*a*(-1.*b[PH]*
          (16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 8.*b[PH]*a*r + 
            4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*cos(2.*th) + 
            4.*b[RR]*dxdxp[1]*pow(a,2.) - 1.*b[PH]*pow(a,3.) + 
            b[PH]*cos(4.*th)*pow(a,3.) + 8.*b[RR]*dxdxp[1]*pow(r,2.) - 
            4.*b[PH]*a*pow(r,2.)) + 
         u[PH]*(bsq + en*gam + rho)*
          (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
            4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
            4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
            u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
            4.*u[PH]*a*pow(r,2.)))*pow(sigma,-4.)*
       (pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
         2.*r*(pow(1.,2.)*(-2.*sigma + 4.*pow(r,2.)) + pow(sigma,2.)) + 
         cos(2.*th)*pow(a,2.)*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)))*
       pow(sin(th),2.) + 8.*dxdxp[1]*pow(sigma,-3.)*
       (bsq + 2.*en*(-1. + gam) - 
         1.*b[PH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
            b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
               b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
            b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.))*pow(sigma,-1.)*
          pow(sin(th),2.) + u[PH]*(bsq + en*gam + rho)*
          (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
               dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
               u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
            u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.))*pow(sigma,-1.)*
          pow(sin(th),2.))*(r*pow(sigma,2.) + 
         pow(a,2.)*(-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)
         ) + 2.*pow(sigma,-3.)*(4.*bsq + 8.*en*(-1. + gam) - 
         1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
            8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
             cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
            1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
            8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
          pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
          (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
            4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
            4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
            u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
            4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.))*
       (-1.*dxdxp[1]*(2. + r)*pow(r,3.) + pow(sigma,3.) + 
         dxdxp[1]*pow(a,2.)*(2.*r*pow(cos(th),2.) + 
            pow(a,2.)*pow(cos(th),4.) + 
            (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))) + 
      16.*dxdxp[1]*a*pow(sigma,-4.)*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)*
       (b[PH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
               2.*b[PH]*a*pow(sin(th),2.))) - 
         1.*u[PH]*(bsq + en*gam + rho)*
          (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
               u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) - 
      32.*pow(dxdxp[1],2.)*pow(sigma,-4.)*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(1. + r) + pow(a,2.)*pow(cos(th),2.))*
       (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
               2.*b[PH]*a*pow(sin(th),2.))) + 
         u[RR]*(bsq + en*gam + rho)*
          (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
               u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
      16.*dxdxp[1]*pow(sigma,-3.)*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (0.5*bsq + en*(-1. + gam) + 
         b[TT]*pow(sigma,-1.)*
          (b[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
               2.*b[PH]*a*pow(sin(th),2.))) - 
         1.*u[TT]*(bsq + en*gam + rho)*pow(sigma,-1.)*
          (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
               u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) - 
      2.*dxdxp[1]*dxdxp[2]*cos(th)*pow(a,2.)*
       (-1.*b[TH]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
            8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
             cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
            1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
            8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
         u[TH]*(bsq + en*gam + rho)*
          (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
            4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
            4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
            u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
            4.*u[PH]*a*pow(r,2.)))*pow(sigma,-2.)*sin(th) - 
      8.*dxdxp[1]*dxdxp[2]*a*cos(th)*
       (-1.*b[TH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
            b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
               b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
            b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
         u[TH]*(bsq + en*gam + rho)*
          (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
               dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
               u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
            u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*sin(th) + 
      16.*dxdxp[1]*dxdxp[2]*r*
       (b[TH]*b[TT] - 1.*u[TH]*u[TT]*(bsq + en*gam + rho))*pow(a,2.)*
       pow(sigma,-2.)*sin(2.*th) + 
      16.*dxdxp[2]*r*(b[RR]*b[TH] - 
         1.*u[RR]*u[TH]*(bsq + en*gam + rho))*pow(dxdxp[1],2.)*pow(a,2.)*
       pow(sigma,-2.)*sin(2.*th) - 
      16.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
       (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
               2.*b[PH]*a*pow(sin(th),2.))) - 
         1.*u[TH]*(bsq + en*gam + rho)*
          (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
            r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
               u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) + 
      2.*dxdxp[1]*(-1.*b[TH]*
          (16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 8.*b[PH]*a*r + 
            4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*cos(2.*th) + 
            4.*b[RR]*dxdxp[1]*pow(a,2.) - 1.*b[PH]*pow(a,3.) + 
            b[PH]*cos(4.*th)*pow(a,3.) + 8.*b[RR]*dxdxp[1]*pow(r,2.) - 
            4.*b[PH]*a*pow(r,2.)) + 
         u[TH]*(bsq + en*gam + rho)*
          (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
            4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
            4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
            u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
            4.*u[PH]*a*pow(r,2.)))*pow(sigma,-1.)*
       (-1.*dxdxp[2]*cot(th) + 
         4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
         dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)) - 
      16.*dxdxp[1]*dxdxp[2]*a*
       (b[PH]*b[TH] - 1.*u[PH]*u[TH]*(bsq + en*gam + rho))*
       pow(sigma,-2.)*sin(th)*(r*(2. + r)*sigma*cos(th) + 
         sigma*pow(a,2.)*pow(cos(th),3.) + r*pow(a,2.)*sin(th)*sin(2.*th)));

  }
  else{

dU[U1]+=0.5*dxdxp[1]*pow(sigma,-5.)*
    (r*pow(sigma,4.)*(bsq - 2.*en + 2.*en*gam - 
         2.*sigma*pow(dxdxp[2],2.)*pow(b[TH],2.) + 
         (bsq + en*gam + rho)*pow(dxdxp[2],2.)*
          (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[TH],2.)) + 
      a*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (-1.*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
          (r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.)) - 2.*a*
          ((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            0.5*dxdxp[1]*(bsq + en*gam + rho)*u[TT]*u[RR]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         4.*a*r*sigma*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
            (bsq + en*gam + rho)*pow(u[TT],2.)))*pow(sin(th),2.) + 
      dxdxp[1]*a*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (-4.*a*r*((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*
          pow(dxdxp[1],-1.) + 
         sigma*(r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.))*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           - 1.*a*pow(dxdxp[1],-1.)*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          ((bsq + 2.*en*(-1. + gam))*((-2. + r)*r + pow(a,2.)) - 
            2.*sigma*pow(dxdxp[1],2.)*pow(b[RR],2.) + 
            (bsq + en*gam + rho)*pow(dxdxp[1],2.)*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[RR],2.)))*
       pow(sin(th),2.) - 1.*sigma*
       (4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*pow(sigma,-1.)*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
         2.*r*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
            (bsq + en*gam + rho)*pow(u[TT],2.)) - 
         1.*a*(-1.*b[TT]*b[PH] + (bsq + en*gam + rho)*u[TT]*u[PH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)) + 
      (4.*a*r*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH]) - 
         1.*a*(a*(bsq + 2.*en*(-1. + gam)) - 
            1.*dxdxp[1]*b[RR]*b[PH]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)) + 
            dxdxp[1]*(bsq + en*gam + rho)*u[RR]*u[PH]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
         0.5*(r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.))*((bsq + 2.*en*(-1. + gam))*pow(csc(th),2.) - 
            2.*sigma*pow(b[PH],2.) + 2.*(bsq + en*gam + rho)*sigma*pow(u[PH],2.))
         )*pow(sin(th),2.)*(r*pow(sigma,2.) + 
         pow(a,2.)*(-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)
         ) + 2.*a*sigma*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)*
       (r*(a*(bsq + 2.*en*(-1. + gam)) - 2.*dxdxp[1]*sigma*b[RR]*b[PH] + 
            2.*dxdxp[1]*(bsq + en*gam + rho)*sigma*u[RR]*u[PH])*pow(sigma,-1.)\
          - 1.*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
          ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) - 
         2.*a*r*(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             pow(csc(th),2.) - 1.*pow(b[PH],2.) + 
            (bsq + en*gam + rho)*pow(u[PH],2.))*pow(sin(th),2.)) + 
      a*sigma*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
         2.*r*(pow(1.,2.)*(-2.*sigma + 4.*pow(r,2.)) + pow(sigma,2.)) + 
         cos(2.*th)*pow(a,2.)*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)))*
       pow(sin(th),2.)*(2.*r*(-1.*b[TT]*b[PH] + 
            (bsq + en*gam + rho)*u[TT]*u[PH]) + 
         dxdxp[1]*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         1.*a*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          (0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*pow(csc(th),2.) - 
            1.*pow(b[PH],2.) + (bsq + en*gam + rho)*pow(u[PH],2.))*
          pow(sin(th),2.)) + sigma*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       ((bsq + 2.*en*(-1. + gam))*sigma + 
         2.*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.))*pow(b[TT],2.) - 
         1.*(bsq + en*gam + rho)*
          (2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.))*pow(u[TT],2.) + 
         4.*r*b[TT]*(-1.*dxdxp[1]*b[RR] + a*b[PH]*pow(sin(th),2.)) + 
         4.*r*(bsq + en*gam + rho)*u[TT]*
          (dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.))) + 
      2.*sigma*(2.*r*(-1.*b[TT]*b[RR] + (bsq + en*gam + rho)*u[TT]*u[RR] + 
            (bsq + 2.*en*(-1. + gam))*r*pow(dxdxp[1],-1.)*pow(sigma,-1.)) + 
         dxdxp[1]*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          (0.5*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-2.)*
             ((-2. + r)*r + pow(a,2.))*pow(sigma,-1.) - 1.*pow(b[RR],2.) + 
            (bsq + en*gam + rho)*pow(u[RR],2.)) - 
         1.*a*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*
       (-1.*dxdxp[1]*(2. + r)*pow(r,3.) + pow(sigma,3.) + 
         dxdxp[1]*pow(a,2.)*(2.*r*pow(cos(th),2.) + 
            pow(a,2.)*pow(cos(th),4.) + 
            (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))) + 
      4.*dxdxp[1]*sigma*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(1. + r) + pow(a,2.)*pow(cos(th),2.))*
       (b[TT]*b[RR]*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.)) - 
         2.*dxdxp[1]*r*pow(b[RR],2.) + 2.*a*r*b[RR]*b[PH]*pow(sin(th),2.) + 
         0.5*(bsq + en*gam + rho)*u[RR]*
          (-1.*u[TT]*(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.)) + 
            4.*r*(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.)))) - 
      1.*dxdxp[2]*a*cos(th)*pow(sigma,2.)*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (4.*a*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH]) - 
         1.*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*
          (r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.)) + 2.*dxdxp[1]*a*
          (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)))*sin(th) - 
      2.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,3.)*
       (2.*r*(-1.*b[TT]*b[TH] + (bsq + en*gam + rho)*u[TT]*u[TH]) + 
         dxdxp[1]*(-1.*b[RR]*b[TH] + (bsq + en*gam + rho)*u[RR]*u[TH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         1.*a*(-1.*b[TH]*b[PH] + (bsq + en*gam + rho)*u[TH]*u[PH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*sin(th) + 
      2.*dxdxp[2]*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
       pow(a,2.)*pow(sigma,3.)*sin(2.*th) + 
      2.*dxdxp[1]*dxdxp[2]*r*
       (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*pow(a,2.)*pow(sigma,3.)*
       sin(2.*th) - 2.*dxdxp[2]*r*pow(a,2.)*pow(sigma,2.)*
       (-2.*dxdxp[1]*r*(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH]) - 
         1.*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
          ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) + 
         2.*a*r*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sin(th),2.)
         )*sin(2.*th) - 2.*dxdxp[2]*a*
       (b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sigma,3.)*sin(th)*
       (r*(2. + r)*sigma*cos(th) + sigma*pow(a,2.)*pow(cos(th),3.) + 
         r*pow(a,2.)*sin(th)*sin(2.*th)));
  }



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU2==1)){

dU[U2]+=dxdxp[1]*r*(-1.*b[RR]*b[TH] + 
       u[RR]*u[TH]*(bsq + en*gam + rho))*pow(dxdxp[2],2.) - 
    1.*(-1.*b[RR]*b[TH] + u[RR]*u[TH]*(bsq + en*gam + rho))*
     (2.*dxdxp[1]*r + sigma)*pow(dxdxp[2],2.) - 
    0.125*r*pow(dxdxp[2],2.)*((-2. + r)*r + pow(a,2.))*
     (-1.*b[TH]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
          8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
           cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
          1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
          8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
       u[TH]*(bsq + en*gam + rho)*
        (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
          4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
          4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
          u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
          4.*u[PH]*a*pow(r,2.)))*pow(sigma,-2.) - 
    0.5*a*r*pow(dxdxp[2],2.)*(-1.*b[TH]*
        (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
          b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
          cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
             b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
          1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
          2.*b[PH]*pow(r,4.)) + 
       u[TH]*(bsq + en*gam + rho)*
        (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
          u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
          cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
             u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
          1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
          2.*u[PH]*pow(r,4.)))*pow(sigma,-2.)*pow(sin(th),2.) - 
    2.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-2.)*
     (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
          r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
             2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[TH]*(bsq + en*gam + rho)*
        (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
          r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
             u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
    2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-3.)*
     (b[PH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
          r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
             2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[PH]*(bsq + en*gam + rho)*
        (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
          r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
             u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*pow(sin(th),3.) - 
    1.*dxdxp[2]*a*r*cos(th)*(-1.*b[TT]*
        (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
          b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
          cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
             b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
          1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
          2.*b[PH]*pow(r,4.)) + 
       u[TT]*(bsq + en*gam + rho)*
        (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
          u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
          cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
             u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
          1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
          2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*sin(th) - 
    0.125*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*
     (4.*bsq + 8.*en*(-1. + gam) - 
       1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
          8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
           cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
          1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
          8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
        pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
        (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
          4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
          4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
          u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
          4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.))*sin(th) - 
    0.5*dxdxp[1]*dxdxp[2]*a*cos(th)*
     (-1.*b[RR]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
          b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
          cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
             b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
          1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
          2.*b[PH]*pow(r,4.)) + 
       u[RR]*(bsq + en*gam + rho)*
        (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
          u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
          cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
             u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
          1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
          2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*
     (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*sin(th) + 
    (0.5*bsq + en*(-1. + gam) - 1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
       (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
     (4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) - 
       1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)) + 
    dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
     (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
          r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
             2.*b[PH]*a*pow(sin(th),2.))) + 
       u[RR]*(bsq + en*gam + rho)*
        (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
          r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
             u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) - 
    1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
     (0.5*bsq + en*(-1. + gam) + 
       b[TT]*pow(sigma,-1.)*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
          r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
             2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[TT]*(bsq + en*gam + rho)*pow(sigma,-1.)*
        (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
          r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
             u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) + 
    0.5*dxdxp[2]*(bsq + 2.*en*(-1. + gam) - 
       1.*b[PH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
          b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
          cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
             b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
          1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
          2.*b[PH]*pow(r,4.))*pow(sigma,-1.)*pow(sin(th),2.) + 
       u[PH]*(bsq + en*gam + rho)*
        (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
          u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
          cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
             u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
          1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
          2.*u[PH]*pow(r,4.))*pow(sigma,-1.)*pow(sin(th),2.))*
     (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th)) + 
    (0.5*bsq + en*(-1. + gam) - 1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
       (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
     (-1.*dxdxp[2]*cot(th) + 4.*(-1.*M_PI*X[2] + th)*
        pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
       dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th));
  }
  else{


dU[U2]+=0.5*(-2.*dxdxp[1]*r*
       (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*pow(dxdxp[2],2.) - 
      1.*a*r*pow(dxdxp[2],2.)*pow(sigma,-2.)*
       (4.*a*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH]) - 
         1.*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*
          (r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.)) + 2.*dxdxp[1]*a*
          (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)))*pow(sin(th),2.) - 
      4.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-2.)*
       (-2.*dxdxp[1]*r*(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH]) - 
         1.*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
          ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) + 
         2.*a*r*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sin(th),2.)
         ) - 2.*r*pow(dxdxp[2],2.)*((-2. + r)*r + pow(a,2.))*pow(sigma,-2.)*
       (2.*r*(-1.*b[TT]*b[TH] + (bsq + en*gam + rho)*u[TT]*u[TH]) + 
         dxdxp[1]*(-1.*b[RR]*b[TH] + (bsq + en*gam + rho)*u[RR]*u[TH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         1.*a*(-1.*b[TH]*b[PH] + (bsq + en*gam + rho)*u[TH]*u[PH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)) + 
      4.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-3.)*
       (r*(a*(bsq + 2.*en*(-1. + gam)) - 2.*dxdxp[1]*sigma*b[RR]*b[PH] + 
            2.*dxdxp[1]*(bsq + en*gam + rho)*sigma*u[RR]*u[PH])*pow(sigma,-1.)\
          - 1.*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
          ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) - 
         2.*a*r*(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             pow(csc(th),2.) - 1.*pow(b[PH],2.) + 
            (bsq + en*gam + rho)*pow(u[PH],2.))*pow(sin(th),2.))*pow(sin(th),3.)
        - 2.*dxdxp[2]*a*r*cos(th)*pow(sigma,-4.)*
       (-1.*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
          (r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.)) - 2.*a*
          ((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            0.5*dxdxp[1]*(bsq + en*gam + rho)*u[TT]*u[RR]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         4.*a*r*sigma*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
            (bsq + en*gam + rho)*pow(u[TT],2.)))*sin(th) - 
      1.*dxdxp[1]*dxdxp[2]*a*cos(th)*pow(sigma,-4.)*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (-4.*a*r*((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*
          pow(dxdxp[1],-1.) + 
         sigma*(r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.))*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           - 1.*a*pow(dxdxp[1],-1.)*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          ((bsq + 2.*en*(-1. + gam))*((-2. + r)*r + pow(a,2.)) - 
            2.*sigma*pow(dxdxp[1],2.)*pow(b[RR],2.) + 
            (bsq + en*gam + rho)*pow(dxdxp[1],2.)*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[RR],2.)))*
       sin(th) - 2.*dxdxp[1]*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-2.)*
       (2.*r*(-1.*b[TT]*b[RR] + (bsq + en*gam + rho)*u[TT]*u[RR] + 
            (bsq + 2.*en*(-1. + gam))*r*pow(dxdxp[1],-1.)*pow(sigma,-1.)) + 
         dxdxp[1]*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          (0.5*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-2.)*
             ((-2. + r)*r + pow(a,2.))*pow(sigma,-1.) - 1.*pow(b[RR],2.) + 
            (bsq + en*gam + rho)*pow(u[RR],2.)) - 
         1.*a*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*sin(th)\
       + (bsq - 2.*en + 2.*en*gam - 2.*sigma*pow(dxdxp[2],2.)*pow(b[TH],2.) + 
         (bsq + en*gam + rho)*pow(dxdxp[2],2.)*
          (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[TH],2.))*
       (4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) - 
         1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)) - 
      1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
       ((bsq + 2.*en*(-1. + gam))*sigma + 
         2.*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.))*pow(b[TT],2.) - 
         1.*(bsq + en*gam + rho)*
          (2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.))*pow(u[TT],2.) + 
         4.*r*b[TT]*(-1.*dxdxp[1]*b[RR] + a*b[PH]*pow(sin(th),2.)) + 
         4.*r*(bsq + en*gam + rho)*u[TT]*
          (dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.)))*sin(2.*th) - 
      2.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
       (b[TT]*b[RR]*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.)) - 
         2.*dxdxp[1]*r*pow(b[RR],2.) + 2.*a*r*b[RR]*b[PH]*pow(sin(th),2.) + 
         0.5*(bsq + en*gam + rho)*u[RR]*
          (-1.*u[TT]*(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.)) + 
            4.*r*(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.))))*sin(2.*th) + 
      dxdxp[2]*pow(sigma,-2.)*
       (4.*a*r*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH]) - 
         1.*a*(a*(bsq + 2.*en*(-1. + gam)) - 
            1.*dxdxp[1]*b[RR]*b[PH]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)) + 
            dxdxp[1]*(bsq + en*gam + rho)*u[RR]*u[PH]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
         0.5*(r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.))*((bsq + 2.*en*(-1. + gam))*pow(csc(th),2.) - 
            2.*sigma*pow(b[PH],2.) + 2.*(bsq + en*gam + rho)*sigma*pow(u[PH],2.))
         )*pow(sin(th),2.)*(cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th)))\
    ;
  }



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU3==1)){

dU[U3]+=0.5*pow(sigma,-1.)*pow(sin(th),2.)*
    ((-1.*b[RR]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
            b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
               b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
            b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
         u[RR]*(bsq + en*gam + rho)*
          (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
               dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
               u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
            u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*
       (-1. - 2.*dxdxp[1]*r*pow(sigma,-1.)) + 
      (-1.*b[TH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
            b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
               b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
            b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
         u[TH]*(bsq + en*gam + rho)*
          (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
               dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*
             (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
               u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
            u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*
       (-1.*dxdxp[2]*cot(th) + 
         4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
         dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)));

  }
  // else nothing to add then


}



// GODMARK: remove? Worthless idea? 2D only also.
// jon's volume diff for uni theta grid and exp r grid (defcoord==LOGRUNITH)
void mks_unitheta_idxvol_func(int i, int j, int k, FTYPE *idxvol)
{

  /*
  int k, l;
  FTYPE r,th,ph;
  FTYPE r1[2],th1[2],r2[2],th2[2];
  FTYPE X0[NDIM],X1[2][NDIM],X2[2][NDIM];
  //  FTYPE cot(FTYPE arg);


  coord(i, j, CENT, X0);
  coord(i, j, FACE1, X1[0]);
  coord(i+1, j, FACE1, X1[1]);
  coord(i, j, FACE2, X2[0]);
  coord(i, j+1, FACE2, X2[1]);

  // get bl coordinates
  bl_coord(X0,&r,&th);
  bl_coord(X1[0],&r1[0],&th1[0]);
  bl_coord(X1[1],&r1[1],&th1[1]);
  bl_coord(X2[0],&r2[0],&th2[0]);
  bl_coord(X2[1],&r2[1],&th2[1]);

  */




  /* comment out for now until adjust everything

  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }
  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif
  */

#define IDXR(a,R0,r,th,rl,rh) ((pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))/((rh - 1.*rl)*(2.*R0 + rh + rl) + (pow(a,2) + 2.*pow(R0,2))*(log(-1.*R0 + rh) - 1.*log(-1.*R0 + rl)) + pow(a,2)*cos(2.*th)*(log(-1.*R0 + rh) - 1.*log(-1.*R0 + rl))))

#define IDXTH(a,R0,r,th,thl,thh) ((-3.*M_PI*(pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))*sin(th))/(cos(thh)*(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*thh)) - 1.*cos(thl)*(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*thl))))


  //#define IDXTH(a,R0,r,th,thl,thh) ((3.*(pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))*sin(th))/((-1.*cos(thh) + cos(thl))*(2.*(pow(a,2) + 3.*pow(r,2)) + pow(a,2)*(cos(2.*thh) + 2.*cos(thh)*cos(thl) + cos(2.*thl)))))

  // fullth integrated -- not used currently
#define FIDXR(a,R0,r,th) ((pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))/(4.*R0*(-1.*R0 + r) + pow(-1.*R0 + r,2) + (pow(a,2) + 2.*pow(R0,2) + pow(a,2)*cos(2.*th))*log(-1.*R0 + r)))

#define FIDXTH(a,R0,r,th) ((-3.*pow(a,2)*sin(2.*th) - 6.*pow(r,2)*tan(th))/(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*th)))


  /*

  idxvol[TT]=1.0; // really 1/dt, but changes in time
  // comment out non-volume regularization
  //  idxvol[RR]=1.0/dx[1];
  //idxvol[TH]=1.0/dx[2];
  idxvol[RR]=IDXR(a,R0,r,th,r1[0],r1[1]);
  idxvol[TH]=IDXTH(a,R0,r,th,th2[0],th2[1]);
  idxvol[PH]=1.0/dx[3];

  fprintf(fail_file,"%d %d %21.15g %21.15g\n",i,j,idxvol[RR]*dx[1],idxvol[TH]*dx[2]);
  
  */




}
