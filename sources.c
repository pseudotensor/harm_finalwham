#include "decs.h"


int sourcephysics(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE (*dUcomp)[NPR])
{
  int coolfunc_thindisk(FTYPE h_over_r, FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE (*dUcomp)[NPR]);
  int coolfunc_neutrino(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE (*dUcomp)[NPR]);

  

  //  SCLOOP(sc) PLOOP(k) dUcomptmp[sc][k] = 0.; // use when duplicating source type (sc) in dUcomp[sc][k]
  // currently all are unique so can overlap

 // cooling function returns energy density per unit time in coordinate frame, typically computed in comoving frame inside the cooling function
  if(cooling==1){
    return(coolfunc_thindisk(h_over_r, pr, ptrgeom, q,dUcomp));
  }
  else if(cooling==2){
    return(coolfunc_neutrino(pr, ptrgeom, q,dUcomp));
  }

  // random physics
  //misc_source(ph, geom, &q, dU, Dt) ;

  return(0);
}

/* JON: here's the cooling function, w/ two parameters */

#define THETACOOL       (h_over_r)	/* should be same as h_over_r */
#define TAUCOOL         (1.0)	        /* cooling time in number of rotational times : really TAUCOOL=2*M_PI would be 1 rotational time */
#define NOCOOLTHFACT    (3.0)           /* this times h_over_r and no more cooling there*/
#define COOLTAPER1(th)   (exp(-pow((th)-M_PI*0.5,2.0)/(2.0*pow(NOCOOLTHFACT*h_over_r,2.0))))
#define SLOPE2 (1.0/(M_PI*0.5-NOCOOLTHFACT*h_over_r))
#define COOLTAPER2(th)   (((th)<M_PI*0.5-NOCOOLTHFACT*h_over_r) ? (SLOPE2*(th)) : ( ((th)>M_PI*0.5+NOCOOLTHFACT*h_over_r) ? (-SLOPE2*((th)-M_PI)) : 1.0 ) )
#define SLOPE3 (1.0/(NOCOOLTHFACT*h_over_r))
#define WIDTHTAPER (NOCOOLTHFACT*h_over_r)
#define TAPERPOS1 (M_PI*0.5-NOCOOLTHFACT*h_over_r-WIDTHTAPER)
#define TAPERPOS2 (M_PI*0.5-NOCOOLTHFACT*h_over_r)
#define TAPERPOS3 (M_PI*0.5+NOCOOLTHFACT*h_over_r)
#define TAPERPOS4 (M_PI*0.5+NOCOOLTHFACT*h_over_r+WIDTHTAPER)
#define TAPERFUN1(th) (0.0)
#define TAPERFUN2(th) (SLOPE3*((th)-TAPERPOS1))
#define TAPERFUN3(th) (1.0)
#define TAPERFUN4(th) (-SLOPE3*((th)-TAPERPOS4))
#define TAPERFUN5(th) (0.0)
#define COOLTAPER3(th)   (((th)<TAPERPOS1 ? TAPERFUN1(th) : (  ((th)<TAPERPOS2) ? TAPERFUN2(th) : (  ((th)<TAPERPOS3) ? TAPERFUN3(th): (  ((th)<TAPERPOS4) ? TAPERFUN4(th) : TAPERFUN5(th)   )))) )
/* cooling function, if any */

FTYPE rin;

int coolfunc_thindisk(FTYPE h_over_r, FTYPE *pr, struct of_geom *geom, struct of_state *q,FTYPE (*dUcomp)[NPR])
{
        FTYPE X[NDIM],V[NDIM],r,th,R,Wcirc,cs_circ,rho,u,P,w,wcirc,dUcool;
	FTYPE taper0;
	int ii,jj, kk;
	FTYPE pressure;

	ii=geom->i;
	jj=geom->j;
	kk=geom->k;

        /* cooling function for maintaining fixed H/R */
        rho = pr[RHO] ;
        u = pr[UU] ;
	P=pressure_rho0_u(rho,u);
        w = rho + u + P ;

        coord(ii,jj,kk,CENT,X) ;
        bl_coord(X,V) ;
	r=V[1];
	th=V[2];
        R = r*sin(th) ;


	rin = (1. + h_over_r)*Risco;

        /* crude approximation */
        Wcirc = pow(R,-1.5) ;
        cs_circ = THETACOOL/sqrt(R) ;
	//        wcirc = rho*(1. + cs_circ*cs_circ/(gam - 1.)) ;
        wcirc = rho*(1. + cs_circ*cs_circ) ;

        if(t > 0.){
	  dUcool = -(Wcirc/TAUCOOL)*( (w - wcirc)*(q->ucon[TT])*(q->ucov[TT])) ;
	  // shape function to avoid problems near pole
	  //taper0=COOLTAPER(0);
	  //dUcool*=1./(1.-1./taper0)+1./(1.-taper0)*COOLTAPER(th);
	  //dUcool*=COOLTAPER2(th);
	  dUcool*=COOLTAPER3(th);
	  dUcool*=taper_func(R,Rhor); // don't cool inside horizon
	}
        else{
	  dUcool = 0. ;
	}

	dUcomp[RADSOURCE][UU]=dUcool;

        return(0) ;
}





int coolfunc_neutrino(FTYPE *pr, struct of_geom *geom, struct of_state *q,FTYPE (*dUcomp)[NPR])
{
        FTYPE X[NDIM],V[NDIM],r,th,R,w,rho,u,P,dUcool ;
	FTYPE rhocode,ucode,Pcode,Rcode,Tcode;
	int ii,jj,kk;
	FTYPE swsq,cv,ca,cvp,cap,mue,nf,ca2,cv2,cvp2,cap2;
	FTYPE Tk,Tmev,lambda,xi,lambdap5,lambda2,lambda3,lambda4,lambda6,lambda8,lambda9;
	FTYPE xi2,xi3,Tk2,lgTk,lgden;
	FTYPE a0,a1,a2,b1,b2,b3,c1;
	FTYPE phys2code;
	FTYPE glambda,fpair,qpair,Qpair,Qpaircode;
	FTYPE pgam2,pgam,pgamp5,pgam3o2,pgam7o2,pgam6;
	FTYPE ft,fl,xp,yp,minfun,minf,fxy2,fxy,Qv,Qplasma,Qplasmacode;
	FTYPE xnucode,Qecap,etae,Qecapcode,Qbrem,Qbremcode;
	FTYPE Ev,Tv,Tvp,Kt,HOR,ntau,EXT;
	FTYPE Qntot,cofactor;
	FTYPE Tfeok,height,zhere,ztogo;

	ii=geom->i;
	jj=geom->j;
	kk=geom->k;

        rhocode = pr[RHO];
        ucode = pr[UU];    
	Pcode=pressure_rho0_u(rhocode,ucode);
	Tcode=Pcode/rhocode;

        coord(ii,jj,kk,CENT,X) ;
        bl_coord(X,V) ;
	r=V[1];
	th=V[2];
        Rcode = r*sin(th) ;


	// for all processes Itoh 1989/1996
	swsq=0.2319;
	cv=0.5+2.*swsq;
	ca=0.5;
	cvp=1.0-cv;
	cap=1.0-ca;
	mue=2.0; // like most normal matter, including n+p+e
	nf=2.0; // flavors conditioned


	// go ahead and scale primitives used below
	R=Rcode*Lunit;
	rho=rhocode*rhounit;
	u=ucode*rhounit;
	P=Pcode*Pressureunit;
	Tk=Tcode*Tempunit;
	Tmev=kb*Tk/ergPmev;
	// from Itoh 1983 or Itoh 2002
	Tfeok=5.9302E9*(pow(1.+1.018*pow(rho/1E6/mue,2./3.),0.5)-1.0)/Tk; // degen e, then Tfeok >>1

	// should solve for electron chemical potential from Kuhri & Mineshige 2002 in deg, ultrarel. e- case.

	// now should all be dimensional cgs units except for final code result output for each cooling process


	ca2=ca*ca;
	cv2=cv*cv;
	cvp2=cvp*cvp;
	cap2=cap*cap;


	lambda=Tk/(5.9302E9);
	xi=pow(rho*rhounit/mue/1E9,1.0/3.0)/lambda;

	lambdap5=sqrt(lambda);
	lambda2=lambda*lambda;
	lambda3=lambda2*lambda;
	lambda4=lambda2*lambda2;
	lambda6=lambda2*lambda4;
	lambda8=lambda4*lambda4;
	lambda9=lambda8*lambda;

	xi2=xi*xi;
	xi3=xi2*xi;

	Tk2=Tk*Tk;
	lgTk=log10(Tk);
	lgden=log10(2.0*rho/mue);

	phys2code=1.0/(edotunit/pow(Lunit,3.0));

	if(Tfeok<1.0){ // then not e. degen and pair process important ala Korhi

	// for pair capture on neutrons
	a0=6.002E19;
	a1=2.084E20;
	a2=1.872E21;
	b1=(Tk<1E10) ? 9.383E-1 : 1.2383;
	b2=(Tk<1E10) ? -4.141E-1 : -0.8141;
	b3=(Tk<1E10) ? 5.829E-2 : 0.0;
	c1=(Tk<1E10) ? 5.5924 : 4.9924;


	
	glambda=1.0-13.04*lambda2+133.5*lambda4+1534.0*lambda6+918.6*lambda8;
	fpair=(a0+a1*xi+a2*xi2)*exp(-c1*xi)/(xi3+b1/lambda+b2/lambda2+b3/lambda3);
	qpair=pow(10.7480*lambda2+0.3967*lambdap5+1.0050,(-1.00))*pow(1.0+(rho*rhounit/mue)*pow(7.692E7*lambda3+9.715E6*lambdap5,(-1.0)),(-0.3));
	Qpair=0.5*((cv2+ca2)+nf*(cvp2+cap2))*(1.0+((cv2-ca2)+nf*(cvp2+cap2))*qpair/((cv2+ca2)+nf*(cvp2+cap2)))*glambda*exp(-2.0/lambda)*fpair;
	Qpaircode=Qpair*phys2code;
	//	qpairtot=-SUM(gdet*dV*Qpaircode*ud0*uu0)*edotunit;
	}
	else Qpaircode=0;

	if(Tfeok>1.0){ // then e. degen and plasma process important and formula applicable ala Korhi,Itoh

	// for plasma process, valid only in high electron degeneracy limit
	pgam2=1.1095E11*rho/mue/(Tk2*pow(1+pow(1.019E-6*rho/mue,(2./3.)),(1./2.)));
	pgam=sqrt(pgam2);
	pgamp5=sqrt(pgam);
	pgam3o2=pow(pgamp5,3.0);
	pgam7o2=pow(pgamp5,7.0);
	pgam6=pgam2*pgam2*pgam2;

	ft=2.4+0.6*pgamp5+0.51*pgam+1.25*pgam3o2;
	fl=(8.6*pgam2+1.35*pgam7o2)/(225.-17.*pgam+pgam2);
	xp=1./6.*(+17.5+lgden-3.*lgTk);
	yp=1./6.*(-24.5+lgden+3.*lgTk);
	minfun=yp-1.6+1.25*xp;
	minf=(minfun>0.0) ? 0.0 : minfun;
	fxy2=1.05+(0.39-1.25*xp-0.35*sin(4.5*xp)-0.3*exp(-pow(4.5*xp+0.9,2.0)))*exp(-pow(minf/(0.57-0.25*xp),2.0));
	fxy=(fabs(xp)>0.7 || yp<0) ? 1.0 : fxy2;
	Qv=3.00E21*lambda9*pgam6*exp(-pgam)*(ft+fl)*fxy;
	Qplasma=(cv2+nf*cvp2)*Qv;
	Qplasmacode=Qplasma*phys2code;
	//qplasmatot=-SUM(gdet*dV*Qplasmacode*ud0*uu0)*edotunit;
	}
	else Qplasmacode=0;

	////////////////////////
	// electron capture

	// fraction of free nucleons
	//	xnucode=26.0*pow(Tk*kb/ergPmev,9.0/8.0)/pow(rho*rhounit/1.E10,3.0/4.0)*exp(-7.074/(Tk*kb/ergPmev));
	xnucode=30.97*pow(Tk/1E10,9.0/8.0)*pow(rho*rhounit/1.E10,-3.0/4.0)*exp(-6.096/(Tk/1E10));
	if(xnucode>1.0) xnucode=1.0;

	if(Tfeok<1.0){ // then not e. degen and e cap process important and formula applicable ala Korhi,Itoh,PWF99
	  Qecap=9.2E33*(rho/1E10)*pow(Tk/1E11,6.0); // PWF99 or Kohri
	}
	else{ // then otherwise use Korhi's eq46
	  etae=Tfeok; // this is WRONG, made up, but best guess
	  Qecap=1.1E31*pow(etae*Tk/1E11,9.0);
	}
	Qecapcode=Qecap*xnucode*phys2code; // final code answer

	/////////////////////////
	// n-n scattering
	if(Tfeok>1.0){ // degen
	  Qbrem=3.4E33*pow(Tk/1E11,8.0)*pow(rho/1E13,1.0/3.0);
	}
	else Qbrem=1.5E33*pow(Tk/1E11,5.5)*pow(rho/1E13,2.0); // non degen
	Qbremcode=Qbrem*phys2code; // final code answer

	///////////////////////////
	// consider optical depth

	// my estimate from Fryer & Meszaros 2003
        /*
        Ev=44.0*pow(rho/1E12/mue,1.0/3.0)*ergPmev;
        Tv=Ev/(21.0*kb);
        Tvp=Tv*kb/ergPmev;
        Kt=1.5E-17*(rho)*pow(Tvp*0.25,2.0);
        HOR=0.26;
        ntau=Kt*HOR*R;
        */

	// from Di Matteo, Perna, Narayan 2002
	// simplified, should use full form, as like also in Popham & Narayan 1995
	// should use HOR calculator rather than fixed value, but per radius!!  and limited to disk, not plunging region
	// consider flavor counting?  not sure what to do there
	HOR=0.26;
	zhere=fabs(R*cos(th));
	ztogo=HOR*fabs(R);
	height=fabs(ztogo-zhere);
	ntau=2.5E-7*pow(Tk/1E11,5.0)*height;
	ntau+=(4.5E-7*xnucode+2.7E-7)*pow(Tk/1E11,2.0)*(rho/1E10)*height;
	// extinction of neutrinos (should add pressure due to neutrinos, but small anyways)
	EXT=exp(-ntau);

	// should perhaps use emission from surface of neutrinosphere ala eq15 of Di Matteo, Perna, Narayan 2002.  nah...just let EXT handle that.

	// include cooling due to photodisintegration, which should then be included in the pressure/internal energy
	// dXdr is radial derivative of Xnuc, must be determined numerically, say using the same method to interpolate, then just use dq/dr.  Might avoid numerical problems in dXdr
	//qphoto=1E29*pow(rho/1E10)*(q->ucon[RR])/(q->ucon[TT])*dXdr;


	///////////////////////////
	// comoving cooling factor

	cofactor=-(q->ucon[TT])*(q->ucov[TT]) ;

        if(t > 0.){
	  // provide energy density per second cooled
	  // multiplying by -uu0*ud0 gives the conserved energy at infinity
	  dUcomp[NEUTRINOANN][UU]=Qpaircode*cofactor*EXT;
	  dUcomp[NEUTRINOPLASMA][UU]=Qplasmacode*cofactor*EXT;
	  dUcomp[NEUTRINOECAP][UU]=Qecapcode*cofactor*EXT;
	  dUcomp[NEUTRINOBREM][UU]=Qbremcode*cofactor*EXT;
	}
	else{
	  dUcomp[NEUTRINOANN][UU]=0;
	  dUcomp[NEUTRINOPLASMA][UU]=0;
	  dUcomp[NEUTRINOECAP][UU]=0;
	  dUcomp[NEUTRINOBREM][UU]=0;
	}

	// rho>~ 10^(11) g/cm^3, optically thick to neutrinos

	// consider limiting dt based upon cooling time scale, or at
	// least to limit changes in internal energy to less than few
	// percent per time step (MW99 265, 2nd to last paragraph)
	// Woosley&Baron 1992: neutrino energy of 15MeV.
	// assumed isotropic emission in comoving frame, so no momentum change
	// mass-energy (rho) is lost?
	

	// consider mass and momentum losses?!
	// see M99 page266, just before section 4.

	// add pair annihilation? MW265 last paragraph
	
	// why is MW99's and Korhi's approximations so bad?
	// MW99 off by factors of 2-3X for Xnuc=1.
	
	// what is the effect of Xnuc in MW99?

	// What about Woosley&Baron1992
	// affect on internal energy
	// 1) added energy by v capture on nucleons (eq4) in optically thin limit outside the neutrinosphere (0 inside)
	// 2) added by electron scatter of neutrinos (eq7)
	// 3) added by v v* annihilation


        return(0) ;
}


#define TRELAX	0.05

void misc_source(FTYPE *pr, struct of_geom *geom, 
		struct of_state *q, FTYPE *dU, FTYPE Dt) 
{
	FTYPE X[NDIM],V[NDIM],r,th ;
	FTYPE rhoscal,uuscal ;
	FTYPE drhodt,dudt,W,isqr,ir ;
	FTYPE ucon_norm[NDIM],ucov_norm[NDIM],alpha ;
	int ii,jj,kk,j ;
	FTYPE dpdrho0,dpdu;


	ii=geom->i;
	jj=geom->j;
	kk=geom->k;

        coord(ii,jj,kk,CENT,X) ;
        bl_coord(X,V) ;
	r=V[1];
	th=V[2];

	ir = 1./r ;
	isqr = 1./sqrt(r) ;
	rhoscal = ir*isqr ;
	uuscal = rhoscal*ir ;

	/* set up normal observer four-velocities, if needed */
	if(pr[RHO] < RHOMIN*rhoscal || pr[UU] < UUMIN*uuscal) {
		alpha = 1./sqrt(-geom->gcon[0][0]) ;
		ucon_norm[TT] = 1./alpha ;
		SLOOPA(j) ucon_norm[j] = geom->gcon[0][j]*alpha ;
		ucov_norm[TT] = -alpha ;
		SLOOPA(j) ucov_norm[j] = 0. ;
	}

	if(pr[RHO] < RHOMIN*rhoscal) {
		drhodt = -0.05*(pr[RHO] - RHOMIN*rhoscal)/dt ;

		dpdrho0 = dpdrho0_rho0_u(pr[RHO],pr[UU]);

		dU[RHO] += (ucon_norm[TT])                                               *drhodt ;
		dU[UU]  += (((1.0+dpdrho0)*ucon_norm[TT]*ucov_norm[TT]) + dpdrho0)       *drhodt ;
		dU[U1]  += ((1.0+dpdrho0)*ucon_norm[TT]*ucov_norm[1] )                   *drhodt ;
		dU[U2]  += ((1.0+dpdrho0)*ucon_norm[TT]*ucov_norm[2] )                   *drhodt ;
		dU[U3]  += ((1.0+dpdrho0)*ucon_norm[TT]*ucov_norm[3] )                   *drhodt ;
	}

	if(pr[UU] < UUMIN*uuscal) {
		dudt = -0.05*(pr[UU] - UUMIN*uuscal)/dt ;

		dpdu = dpdu_rho0_u(pr[RHO],pr[UU]);

		dU[UU] += ((1.0+dpdu)*ucon_norm[TT]*ucov_norm[TT] +    dpdu)  * dudt ;
		dU[U1] += ((1.0+dpdu)*ucon_norm[TT]*ucov_norm[1])             * dudt ;
		dU[U2] += ((1.0+dpdu)*ucon_norm[TT]*ucov_norm[2])             * dudt ;
		dU[U3] += ((1.0+dpdu)*ucon_norm[TT]*ucov_norm[3])             * dudt ;
	}
}
