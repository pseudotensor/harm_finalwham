#include "decs.h"

////////////////////////////
//
// CHOOSE EOS
//
////////////////////////////

#define IDEALGAS 0
#define MIGNONE 1

#define WHICHEOS IDEALGAS
//#define WHICHEOS MIGNONE





////////////////////////////
//
// IMPLEMENT EOS
//
// u_rho0_p : used by initial conditions
//
// pressure_rho0_u : used by inversion for initial guess and by rest of code to set pressure as functions of rho0 and u
// 
// dpdu_rho0_u  dpdrho0_rho0_u : used by sources or other such derivatives
//
// cs2_compute : used by vchar.c for characteristics
//
// pressure_wmrho0 : used by inversion
//
// compute_idwmrho0dp : used by 2D inversion
//
//
// compute_entropy and compute_u_from_entropy : used by entropy evolution and inversion
//
////////////////////////////

#define OLDCALC 0


#if(WHICHEOS==IDEALGAS)

// IDEAL GAS EOS

// P = (\GAMMA -1) u
// h = 1+\Gamma_r \Theta  ; \Theta=p/\rho_0  \Gamma_r = (\GAMMA/(\GAMMA-1))

#define GAMMA (gam)

// 1/\Gamma_r
#define GAMMAM1 (GAMMA-1.0)
#define IGAMMAR (GAMMAM1/GAMMA)


// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
#if(OLDCALC)
  return(GAMMAM1*u) ;
#else
  return((GAMMA - 1.)*u) ;
#endif

}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_grmhd(FTYPE rho0, FTYPE p)
{
	return(p/GAMMAM1) ;
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
	return(GAMMAM1) ;
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
	return(0.0) ;
}

// sound speed squared (for vchar.c)
FTYPE cs2_compute_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE pressure;
  FTYPE h;
  FTYPE cs2;

  pressure = pressure_rho0_u(rho0,u);
  h=rho0+u+pressure;

  cs2=GAMMA*pressure/h;

  return(cs2);

}

// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
  FTYPE pressure,indexn,entropy;


  pressure=pressure_rho0_u(rho0,u);
  indexn=1.0/GAMMAM1;

  if(rho0<0.0) rho0=SMALL;
  if(pressure<0.0) pressure=SMALL;
  
  entropy=rho0*log(pow(pressure,indexn)/pow(rho0,indexn+1.0));

  return(entropy);

}

// u(rho0,S)
FTYPE compute_u_from_entropy_grmhd(FTYPE rho0, FTYPE entropy)
{
  FTYPE rho,ie,pressure;
  FTYPE indexn;
  FTYPE u;

  indexn=1.0/GAMMAM1;

  if(rho0<0.0) rho0=SMALL;
  
  // entropy version of ie
  u=pow(pow(rho0,indexn+1.0)*exp(entropy/rho0),1.0/indexn)/GAMMAM1;

  return(u);

}

// used for dudp_calc
FTYPE compute_dSdrho_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE entropy;
  FTYPE dSdrho;
  FTYPE compute_entropy(FTYPE rho0, FTYPE u);

  entropy=compute_entropy(rho0,u);

  // ideal gas
  indexn=1.0/GAMMAM1;

  dSdrho=entropy/rho0-(indexn+1.0);

  return(dSdrho);

}


// used for dudp_calc
FTYPE compute_dSdu_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE dSdu;

  // ideal gas
  indexn=1.0/GAMMAM1;

  dSdu=indexn*rho0/u;

  return(dSdu);

}



// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  return(IGAMMAR*wmrho0) ;
}

// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  return(GAMMAM1/GAMMA);
}


// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);
}




#elif(WHICHEOS==MIGNONE)

// TM EOS


// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE pressure;

  pressure = u*(2.0*rho0+u)/(3.0*(rho0+u));

  //  dualfprintf(fail_file,"rho0=%g u=%g pressure=%g\n",rho0,u,pressure);

  return(pressure);
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_grmhd(FTYPE rho0, FTYPE p)
{
  return( 1.5*(p + 3.0*p*p/(2.0*rho0+sqrt(9.0*p*p+4.0*rho0*rho0))) );
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE dpdu;

  dpdu = 1.0/3.0*(1.0 + rho0*rho0/( (rho0+u)*(rho0+u)));

  return(dpdu);
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE dpdrho0;

  dpdrho0 = u*u/(3.0*(rho0+u)*(rho0+u));

  return(dpdrho0) ;
}


// sound speed squared (for vchar.c)
FTYPE cs2_compute_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE pressure;
  FTYPE h;
  FTYPE cs2;

  pressure = pressure_rho0_u(rho0,u);
  h=rho0+u+pressure; // not specific h

  cs2=pressure*(5.0*h - 8.0*pressure) / (3.0*h*(h-pressure));


  //  dualfprintf(fail_file,"cs2=%21.15g pressure=%21.15g\n",cs2,pressure);

  return(cs2);

}

// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
  FTYPE entropy;

  entropy=0.0; // GODMARK: not set yet

  return(entropy);

}

// u(rho0,S)
FTYPE compute_u_from_entropy_grmhd(FTYPE rho0, FTYPE entropy)
{
  FTYPE u;

  u=0.0; // GODMARK: not set yet

  return(u);

}


// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE Q,delta,delta2;
  FTYPE pressure;

  Q=wmrho0/rho0;
  delta=9.0/25.0*wmrho0*(2.0+Q);
  delta2=delta/rho0;

  pressure=(5.0/8.0)*(wmrho0 - delta/(1.0+sqrt(1.0+delta2)));

  return(pressure);
}


// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp_grmhd_old(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE Q,delta,delta2;
  FTYPE ddeltadwmrho0,idwmrho0dp;

  Q=wmrho0/rho0;
  delta=9.0/25.0*wmrho0*(2.0+Q);
  delta2=delta/rho0;
  
  ddeltadwmrho0=18.0/25.0*(1.0+Q);

  idwmrho0dp = 5.0/16.0*(2.0-ddeltadwmrho0/sqrt(1.0+delta2));

  return(idwmrho0dp);

}

// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
  FTYPE idwmrho0dp;
  FTYPE p;

  p = pressure_wmrho0(rho0, wmrho0);

  idwmrho0dp = (2.0*wmrho0 + 2.0*rho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);

  return(idwmrho0dp);

}

// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
  FTYPE idrho0dp;
  FTYPE p;

  p = pressure_wmrho0(rho0, wmrho0);

  idrho0dp = (2.0*wmrho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);

  return(idrho0dp);
}




#endif // end if MIGNONE EOS




//////////////////////////////////////////////////////
//
// COLD EOS or COLD EOMTYPE
//
// then explicitly set pressure/ie to 0 so can keep some parts of code similar
//
//////////////////////////////////////////////////////

// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_coldgrmhd(FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_coldgrmhd(FTYPE rho0, FTYPE p)
{
  return(0.0) ;
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_coldgrmhd(FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_coldgrmhd(FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

// sound speed squared (for vchar.c)
FTYPE cs2_compute_coldgrmhd(FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}


// used for dudp_calc
FTYPE compute_dSdrho_coldgrmhd(FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}


// used for dudp_calc
FTYPE compute_dSdu_coldgrmhd(FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}


// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_coldgrmhd(FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

// u(rho0,S)
FTYPE compute_u_from_entropy_coldgrmhd(FTYPE rho0, FTYPE entropy)
{
  return(0.0) ;
}


// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_coldgrmhd(FTYPE rho0, FTYPE wmrho0)
{
  return(0.0) ;
}


// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp_coldgrmhd(FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);
}


// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_coldgrmhd(FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);
}








//////////////////////////////////////////////////
//
// wrappers for any EOS and any EOMTYPE
//
//////////////////////////////////////////////////

// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u)
{
  return( (*ptr_pressure_rho0_u)(rho0,u) );
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p(FTYPE rho0, FTYPE p)
{
  return( (*ptr_u_rho0_p)(rho0,p) );
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u(FTYPE rho0, FTYPE u)
{
  return( (*ptr_dpdu_rho0_u)(rho0,u) );
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u(FTYPE rho0, FTYPE u)
{
  return( (*ptr_dpdrho0_rho0_u)(rho0,u) );
}

// sound speed squared (for vchar.c)
FTYPE cs2_compute(FTYPE rho0, FTYPE u)
{
  return( (*ptr_cs2_compute)(rho0,u) );
}


// used for dudp_calc
FTYPE compute_dSdrho(FTYPE rho0, FTYPE u)
{
  return( (*ptr_compute_dSdrho)(rho0,u) );
}


// used for dudp_calc
FTYPE compute_dSdu(FTYPE rho0, FTYPE u)
{
  return( (*ptr_compute_dSdu)(rho0,u) );
}


// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy(FTYPE rho0, FTYPE u)
{
  return( (*ptr_compute_entropy)(rho0,u) );
}

// u(rho0,S)
FTYPE compute_u_from_entropy(FTYPE rho0, FTYPE entropy)
{
  return( (*ptr_compute_u_from_entropy)(rho0,entropy) );
}


// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0)
{
  return( (*ptr_pressure_wmrho0)(rho0,wmrho0) );
}


// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp(FTYPE rho0, FTYPE wmrho0)
{
  return( (*ptr_compute_idwmrho0dp)(rho0,wmrho0) );
}


// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp(FTYPE rho0, FTYPE wmrho0)
{
  return( (*ptr_compute_idrho0dp) (rho0,wmrho0) );
}







//////////////////////////////////////////
//
// ANY EOS
//

// old function
// p(rho0,w)
FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w)
{
  FTYPE wmrho0=w-rho0;

#if(OLDCALC)
  return((GAMMA-1.)*(w - rho0)/GAMMA) ;
#else
  return(pressure_wmrho0(rho0,wmrho0)) ;
#endif
}


// pick EOMTYPE for EOS
int pickeos_eomtype(int whicheom)
{

  // EOS functions used during inversion and other places
  if(whicheom==EOMGRMHD){
    ptr_pressure_rho0_u = &pressure_rho0_u_grmhd;
    ptr_compute_u_from_entropy = &compute_u_from_entropy_grmhd;
    ptr_u_rho0_p = &u_rho0_p_grmhd;
    ptr_dpdu_rho0_u = &dpdu_rho0_u_grmhd;
    ptr_dpdrho0_rho0_u = &dpdrho0_rho0_u_grmhd;
    ptr_cs2_compute = &cs2_compute_grmhd;
    ptr_compute_dSdrho = &compute_dSdrho_grmhd;
    ptr_compute_dSdu = &compute_dSdu_grmhd;
    ptr_compute_entropy = &compute_entropy_grmhd;

    ptr_pressure_wmrho0=&pressure_wmrho0_grmhd;
    ptr_compute_idwmrho0dp=&compute_idwmrho0dp_grmhd;
    ptr_compute_idrho0dp=&compute_idrho0dp_grmhd;
  }
  else{
    ptr_pressure_rho0_u = &pressure_rho0_u_coldgrmhd;
    ptr_compute_u_from_entropy = &compute_u_from_entropy_coldgrmhd;
    ptr_u_rho0_p = &u_rho0_p_coldgrmhd;
    ptr_dpdu_rho0_u = &dpdu_rho0_u_coldgrmhd;
    ptr_dpdrho0_rho0_u = &dpdrho0_rho0_u_coldgrmhd;
    ptr_cs2_compute = &cs2_compute_coldgrmhd;
    ptr_compute_dSdrho = &compute_dSdrho_coldgrmhd;
    ptr_compute_dSdu = &compute_dSdu_coldgrmhd;
    ptr_compute_entropy = &compute_entropy_coldgrmhd;

    ptr_pressure_wmrho0=&pressure_wmrho0_coldgrmhd;
    ptr_compute_idwmrho0dp=&compute_idwmrho0dp_coldgrmhd;
    ptr_compute_idrho0dp=&compute_idrho0dp_coldgrmhd;
  }

  return(0);

}
