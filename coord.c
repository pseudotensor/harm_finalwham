
#include "decs.h"

/** 
 *
 * this file contains all the coordinate dependent
 * parts of the code, except the initial and boundary
 * conditions 
 *
 **/

// static variables with global scope to this file
// could make any or all of these true global if want to change them in, say, init.c

// for defcoord==COMPLEX1TH
static FTYPE der0=9;
//static FTYPE Ri=8.0; // too unresolved for 256x128
static FTYPE Ri=20.0;

// for defcoord==COMPLEX2TH
static FTYPE x2trans;
static FTYPE m2,d2,c2,m3,b3,myx2,flip1,flip2,thetatores;

// for defcoord==JET1COORDS
static FTYPE h0,hf,rh0,myrout,dmyhslope1dr,myhslope1,myhslope2,myhslope,dmyhslope2dx1,dmyhslopedx1,x1in,x1out;
static FTYPE npow;

// for defcoord=JET2COORDS
static FTYPE r1jet,njet,rpjet,dmyhslopedr;

// for defcoord=JET3COORDS
static FTYPE r0jet,rsjet,Qjet;

// for defcoord=PULSARCOORDS
static FTYPE hinner,houter;


void set_coord_parms(void)
{
  // assumes R0, Rin, Rout, and hslope are so general that are set in init.c
  if (defcoord == LOGRSINTH) {
  }
  else if (defcoord == COMPLEX0TH) {
  }
  else if(defcoord == UNIRSINTH){
  }
  else if (defcoord == EQMIRROR) {
  }
  else if(defcoord == COMPLEX1TH) {
  }
  else if(defcoord == COMPLEX2TH) {
    x2trans=0.1; // user settable, must be same as below in dxdxp
    thetatores=2.5*h_over_r;

    // fixed coefficients
    m2=(3.*(-2.*thetatores + M_PI))/(2.*x2trans) + (4.*thetatores)/(-1. + 2.*x2trans);
    d2=(2.*thetatores - M_PI + 2.*M_PI*x2trans)/(-2.*pow(x2trans,3.) + 4.*pow(x2trans,4.));
    c2=(6.*thetatores - 3.*M_PI + 6.*M_PI*x2trans)/(2.*pow(x2trans,2.) - 4.*pow(x2trans, 3.));
    m3=(2.*thetatores)/(1. - 2.*x2trans);
    b3=M_PI/2. + thetatores/(-1. + 2.*x2trans);
  }
  else if (defcoord == LOGRUNITH) { // uniform theta and log in radius
  }
  else if (defcoord == JET1COORDS) {
    // optimal is npow=10 R0=-3
    npow=1.0;
    //R0=0.0;

    // must be same as in dxdxp()
    h0=hslope;
    hf=2.0-0.22;
    rh0=40.0;
    myrout=Rout;
    dmyhslope1dr = (hf-h0)/(myrout-rh0);
    dmyhslope2dx1=(hf-h0)/(x1out-x1in);
    x1in=log(Rin-R0);
    x1out=log(Rout-R0);

  }
  else if (defcoord == JET2COORDS) {
    npow=1.0;

    // must be same as in dxdxp()
    if(1){
      r1jet=16.0;
      njet=0.3;
      rpjet=0.9;
    }
    else{
      r1jet=9.0;
      njet=0.3;
      rpjet=.9;
    }
  }
  else if (defcoord == JET3COORDS) {
    npow=1.0;

    // must be same as in dxdxp()
    if(0){ // first attempt
      r1jet=2.8;
      njet=0.3;
      r0jet=7.0;
      rsjet=21.0;
      Qjet=1.7;
    }
    else if(0){ // chosen to resolve disk then resolve jet
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.8;
    }
    else if(1){
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.3; // chosen to help keep jet resolved even within disk region
    }
  }
  else if (defcoord == PULSARCOORDS) {

    if(0){// pulsar in force free
      // pulsar_grid.nb for theta part and for the radial part:
      // see pulsar_gridnew.nb
      // for Rout=10^6 and R0=0.786*Rin Rin=4.84, npow=10 gives same dr/r as npow=1 R0=0.9*Rin at r=Rin
      npow=1.0;
      
      // must be same as in dxdxp()
      hinner=hslope; // hslope specifies inner hslope
      houter=hslope*0.05; // reduce by some arbitrary factor (currently 1/20)
      r0jet=5.0; // spread in radius over which hslope changes
      rsjet=18.0; // location of current sheet beginning for NS pulsar
    }
    else if(1){ // NS-pulsar in GRMHD
      npow=10.0;
      hinner=1.9*hslope; // hslope specifies inner hslope
      //houter=hslope*0.001; // reduce by some arbitrary factor (currently 1/20)
      houter=hslope*1.5; // increase houter up to 2.0
      r0jet=5.0; // spread in radius over which hslope changes
      rsjet=15.0; // location of current sheet beginning for NS pulsar
    }

  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
  }
  else if (defcoord == BILOGCYLCOORDS) {
    npow=10.0; // exponential rate
  }
  else if (defcoord == RAMESHCOORDS) {
    //    myhslope=pow( (*r-rsjet)/r0jet , njet);
    npow=1.0;
    //npow=3.0;

    // for Ramesh problem with paraboloidal field (nu=1)
    r0jet=10.0; // divisor
    njet=jetalpha; // power \theta_j \propto r^{-njet}
    //njet=1.0;
    rsjet=0.0; // subtractor
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of set_coord_parms\n");
    myexit(1);
  }

}


void write_coord_parms(void)
{
  FILE *out;
  int dimen;

  if(myid==0){
    if((out=fopen("coordparms.dat","wt"))==NULL){
      dualfprintf(fail_file,"Couldn't write coordparms.dat file\n");
      myexit(1);
    }
    else{

      // same for all coords (notice no carraige return)
      fprintf(out,"%21.15g %21.15g %21.15g %21.15g ",R0,Rin,Rout,hslope);

      if (defcoord == LOGRSINTH) {
      }
      else if (defcoord == COMPLEX0TH) {
      }
      else if(defcoord == UNIRSINTH){
      }
      else if (defcoord == EQMIRROR) {
      }
      else if(defcoord == COMPLEX1TH) {
      }
      else if(defcoord == COMPLEX2TH) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",x2trans,thetatores,m2,d2,c2,m3,b3,h_over_r);
      }
      else if (defcoord == LOGRUNITH) { // uniform theta and log in radius
      }
      else if (defcoord == JET1COORDS) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",npow,h0,hf,rh0,myrout,dmyhslope1dr,dmyhslope2dx1,x1in,x1out);
      }
      else if (defcoord == JET2COORDS) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",npow,r1jet,njet,rpjet);
      }
      else if (defcoord == JET3COORDS) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",npow,r1jet,njet,r0jet,rsjet,Qjet);
      }
      else if (defcoord == PULSARCOORDS) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g\n",npow,hinner,houter,r0jet,rsjet);
      }
      else if (defcoord == UNIFORMCOORDS) {
	//uniform grid for Cartesian coordinates
	DIMENLOOP(dimen) fprintf(out,"%21.15g ",Rin_array[dimen]);
	DIMENLOOP(dimen) fprintf(out,"%21.15g ",Rout_array[dimen]);
	fprintf(out,"\n");
      }
      else if (defcoord == BILOGCYLCOORDS) {
	fprintf(out,"%21.15g\n",npow);
      }
      else if (defcoord == RAMESHCOORDS) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",npow,r0jet,njet,rsjet);
      }
      else{
	dualfprintf(fail_file,"Shouldn't reach end of write_coord_parms\n");
	myexit(1);
      }

      fclose(out);
    }
  }
}



void read_coord_parms(void)
{
  FILE *in;
  FTYPE ftemp;
  int dimen;
  
  if(myid==0){
    in=fopen("coordparms.dat","rt");
    if(in==NULL){
      dualfprintf(fail_file,"Couldn't read coordparms.dat file.  I'll assume coded coordinates and let restart header overwrite any global restart parameters\n");
      set_coord_parms();
    }
    else{
      // don't want to overwrite since restart file sets this
      //      fscanf(in,HEADER4IN,&R0,&Rin,&Rout,&hslope);
      fscanf(in,HEADER4IN,&ftemp,&ftemp,&ftemp,&ftemp);

      if (defcoord == LOGRSINTH) {
      }
      else if (defcoord == COMPLEX0TH) {
      }
      else if(defcoord == UNIRSINTH){
      }
      else if (defcoord == EQMIRROR) {
      }
      else if(defcoord == COMPLEX1TH) {
      }
      else if(defcoord == COMPLEX2TH) {
	fscanf(in,HEADER8IN,&x2trans,&thetatores,&m2,&d2,&c2,&m3,&b3,&h_over_r);
      }
      else if (defcoord == LOGRUNITH) { // uniform theta and log in radius
      }
      else if (defcoord == JET1COORDS) {
	fscanf(in,HEADER9IN,&npow,&h0,&hf,&rh0,&myrout,&dmyhslope1dr,&dmyhslope2dx1,&x1in,&x1out);
      }
      else if (defcoord == JET2COORDS) {
	fscanf(in,HEADER4IN,&npow,&r1jet,&njet,&rpjet);
      }
      else if (defcoord == JET3COORDS) {
	fscanf(in,HEADER6IN,&npow,&r1jet,&njet,&r0jet,&rsjet,&Qjet);
      }
      else if (defcoord == PULSARCOORDS) {
	fscanf(in,HEADER5IN,&npow,&hinner,&houter,&r0jet,&rsjet);
      }
      else if (defcoord == UNIFORMCOORDS) {
	//uniform grid for Cartesian coordinates
	DIMENLOOP(dimen) fscanf(in,HEADERONEIN,&Rin_array[dimen]);
	DIMENLOOP(dimen) fscanf(in,HEADERONEIN,&Rout_array[dimen]);
      }
      else if (defcoord == BILOGCYLCOORDS) {
	fscanf(in,HEADERONEIN,&npow);
      }
      else if (defcoord == RAMESHCOORDS) {
	fscanf(in,HEADER4IN,&npow,&r0jet,&njet,&rsjet);
      }
      else{
	dualfprintf(fail_file,"Shouldn't reach end of read_coord_parms\n");
	myexit(1);
      }

      fclose(in);
    }
  }

#if(USEMPI)
  // broadcast
  MPI_Bcast(&R0, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Rin, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Rout, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hslope, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);

  if (defcoord == LOGRSINTH) {
  }
  else if (defcoord == COMPLEX0TH) {
  }
  else if(defcoord == UNIRSINTH){
  }
  else if (defcoord == EQMIRROR) {
  }
  else if(defcoord == COMPLEX1TH) {
  }
  else if(defcoord == COMPLEX2TH) {
    MPI_Bcast(&x2trans, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&thetatores, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m2, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d2, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c2, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m3, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b3, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    //  MPI_Bcast(&h_over_r, 1, MPI_FTYPE, 0, MPI_COMM_WORLD); // set by pre_init_specific_init() in init.c
  }
  else if (defcoord == LOGRUNITH) { // uniform theta and log in radius
  }
  else if (defcoord == JET1COORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h0, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&hf, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rh0, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&myrout, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dmyhslope1dr, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dmyhslope2dx1, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x1in, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x1out, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == JET2COORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rpjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == JET3COORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Qjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == PULSARCOORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&hinner, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&houter, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    DIMENLOOP(dimen) MPI_Bcast(&Rin_array[dimen], 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    DIMENLOOP(dimen) MPI_Bcast(&Rout_array[dimen], 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == BILOGCYLCOORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == RAMESHCOORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of read_coord_parms\n");
    myexit(1);
  }
  
#endif

}





/* Returns boyer-lindquist coordinte of point */
void bl_coord(FTYPE *X, FTYPE *V)
{
  extern FTYPE mysin(FTYPE th);

  // in spherical polar coords: t=V[0] r=V[1] th=V[2] phi=V[3]

  if (defcoord == LOGRSINTH) {
    //    V[1] = Rin*exp(X[1]) ;
    V[1] = R0+exp(X[1]) ;
    //V[1] = Rin * exp(X[1]);
    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - hslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      //      V[2] = 0.5*M_PI + M_PI * fabs(X[2]-0.5) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }

    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];

  } else if (defcoord == COMPLEX0TH) {
    V[1] = R0+Rin * exp(X[1] * log(Rout / Rin));
    V[2] =
      ((-49. * hslope + 60. * M_PI) * X[2]) / 12. +
      ((247. * hslope - 240. * M_PI) * pow(X[2],2)) / 12. +
      ((-83. * hslope + 80. * M_PI) * pow(X[2],3)) / 2. -
      (5. * (-25. * hslope + 24. * M_PI) * pow(X[2], 4)) / 3. +
      (2. * (-25. * hslope + 24. * M_PI) * pow(X[2], 5)) / 3.;
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if(defcoord==UNIRSINTH){
    V[1]=X[1];
    V[2] = M_PI* X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == EQMIRROR) {
    // MIRROR at equator, equator is outer theta edge
    V[1] = R0+exp(X[1]) ;
    V[2] = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if(defcoord == COMPLEX1TH) {
    V[1] = R0+exp(X[1]) ;

    V[2] = (der0*X[2]*(-32.*pow(-1. + X[2],3.)*pow(X[2],2.)*(-1. + 2.*X[2]) - 
		       Ri*(-1. + X[2])*pow(-1. + 2.*X[2],3.)*
		       (-1. + 7.*(-1. + X[2])*X[2])) + 
	    M_PI*Ri*pow(X[2],3.)*(70. + 
				  3.*X[2]*(-105. + 2.*X[2]*(91. + 10.*X[2]*(-7. + 2.*X[2])))))/Ri;

    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if(defcoord == COMPLEX2TH) {

    V[1] = R0+exp(X[1]) ;

    // now assign values
    if(X[2]<0.5){ myx2=X[2]; flip1=0.0; flip2=1.0;}
    else{ myx2=1.0-X[2]; flip1=M_PI; flip2=-1.0;}

    if(myx2<=x2trans){
      V[2] = flip1+flip2*(d2*pow(myx2,3.0)+c2*pow(myx2,2.0)+m2*myx2);
    }
    else{
      V[2] = flip1+flip2*(m3*myx2+b3);
    }
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];

  }
  else if (defcoord == LOGRUNITH) { // uniform theta and log in radius
    V[1] = R0+exp(X[1]) ;
    V[2] = M_PI * X[2] ;
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == JET1COORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope1 = h0+dmyhslope1dr*((V[1])-rh0);
    myhslope2 = h0+(hf-h0)*(X[1]-x1in)/(x1out-x1in);

    myhslope=(myhslope1+myhslope2)*0.5;
    V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == JET2COORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope=2.0-pow(V[1]/r1jet,njet*(-1.0+exp(1.0)/exp(V[1]+rpjet)));

    V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == JET3COORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));

    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == PULSARCOORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope=(0.5+1.0/M_PI*atan((V[1]-rsjet)/r0jet))*(houter-hinner)+hinner;

    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    V[1] = Rin_array[1] + X[1] * ( Rout_array[1] - Rin_array[1] );
    V[2] = Rin_array[2] + X[2] * ( Rout_array[2] - Rin_array[2] );
    V[3] = Rin_array[3] + X[3] * ( Rout_array[3] - Rin_array[3] );
  }
  else if (defcoord == BILOGCYLCOORDS) {
    // R : cylindrical radius, assumes X[1]=0..1
    // exponential grid
    V[1] = ((Rout-Rin)*exp(npow*X[1])+Rin*exp(npow)-Rout)/(exp(npow)-1.0);
    
    // z : cylindrical height, assumes X[2]=-1..1
    // bi-exponential grid
    // here the grid goes from Zin to Zout in a bi-log way, and X[2]=0 is Z=0
    if(X[2]>0.0) V[2] = ((Zout-0)*exp(npow*fabs(X[2])) + 0*exp(npow)-Zout)/(exp(npow)-1.0);
    else V[2] = ((Zin-0)*exp(npow*fabs(X[2])) + 0*exp(npow)-Zin)/(exp(npow)-1.0);
    
  } 
  else if (defcoord == RAMESHCOORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope=pow( (V[1]-rsjet)/r0jet , njet);

    V[2] = 0.5*M_PI*(1.0 + atan(myhslope*(X[2]-0.5))/atan(myhslope*0.5));

    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];


  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of bl_coord\n");
    myexit(1);
  }

  // don't allow to be smaller to avoid singularity
  // noted this caused problems with jon_interp in calculating jacobian
  if(POSDEFMETRIC){
    if(V[2]<0) V[2] = -V[2];
    if(V[2]>M_PI) V[2]=2.0*M_PI-V[2];

    if(V[3]<0) V[3] = -V[3];
    if(V[3]>2.0*M_PI) V[3]=4.0*M_PI-V[3];

  }

	#if( COORDSINGFIXCYL )   //SUPERSASMARK fix the singularity for the cylinrical coordinates
		if (fabs(V[1]) < SINGSMALL){
			if(V[1]>=0) V[1]=SINGSMALL;
			if(V[1]<0) V[1]=-SINGSMALL;
		}
	#endif

  if(COORDSINGFIX){
			if (fabs(V[2]) < SINGSMALL){
				if(V[2]>=0) V[2]=SINGSMALL;
				if(V[2]<0) V[2]=-SINGSMALL;
			}
			if (fabs(M_PI-V[2]) < SINGSMALL){
				if(V[2]>=M_PI) V[2]=M_PI+SINGSMALL;
				if(V[2]<M_PI) V[2]=M_PI-SINGSMALL;
			}
  }



}



// Jacobian for dx uniform per dx nonuniform (dx/dr / dx/dr')
// i.e. Just take d(bl-coord)/d(ksp uniform coord)
// e.g. dr/dx1 d\theta/dx2

// take note of the ordering of indicies
// dxdxp[j][k]=dxdxp[mu][nu]=(dx^\mu_{BL}/dx^\nu_{KSP uni})

// should make this numerical like connection, then to conserve CPU, would need all over grid
void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM])
{
  void dxdxp_numerical(FTYPE *X, FTYPE (*dxdxp)[NDIM]);
  void dxdxp_analytic(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);

  if(defcoord<=ANALYTICSWITCH){ // then have analytic dxdxp
    dxdxp_analytic(X,V,dxdxp);
  }
  else{
    dxdxp_numerical(X,dxdxp);
  }

}

// should make this numerical like connection, then to conserve CPU, would need all over grid
void dxdxp_analytic(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM])
{
  int j,k;
  extern FTYPE mycos(FTYPE th);
  extern FTYPE mysin(FTYPE th);


  // default identity transformation
  DLOOP(j,k) dxdxp[j][k]=0.0;
  DLOOPA(j) dxdxp[j][j]=1.0;

  if (defcoord == LOGRSINTH) {
    dxdxp[1][1] = V[1]-R0;
    if(X[2]<0.5){
      dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * mycos(2. * M_PI * X[2]);
    }
    else{
      dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * mycos(2. * M_PI * (1.0-X[2]) );
    }
    dxdxp[3][3] = 2.0*M_PI;
  } else if (defcoord == COMPLEX0TH) {
    dxdxp[1][1] = (V[1]-R0) * log(Rout / Rin);

    dxdxp[2][2] = (-49. * hslope + 60. * M_PI) / 12. +
      ((247. * hslope - 240. * M_PI) * X[2]) / 6. +
      (3. * (-83. * hslope + 80. * M_PI) * pow(X[2], 2)) / 2. -
      (20. * (-25. * hslope + 24. * M_PI) * pow(X[2], 3)) / 3. +
      (10. * (-25. * hslope + 24. * M_PI) * pow(X[2], 4)) / 3.;

    dxdxp[3][3] = 2.0*M_PI;

  } else if(defcoord == UNIRSINTH){
    dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[3][3] = 2.0*M_PI;
  }
  else if (defcoord == EQMIRROR) {
    dxdxp[1][1] = V[1]-R0;
    dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[3][3] = 2.0*M_PI;
  } 
  else if (defcoord == COMPLEX1TH) {
    dxdxp[1][1] = V[1]-R0;
    dxdxp[2][2] = (210.*M_PI*Ri*pow(1. - 2.*X[2],2.)*pow(-1. + X[2],2.)*
		   pow(X[2],2.) + der0*
		   (-32.*pow(-1. + X[2],2.)*pow(X[2],2.)*
		    (3. + 14.*(-1. + X[2])*X[2]) - 
		    Ri*pow(1. - 2.*X[2],2.)*
		    (-1. + 2.*(-1. + X[2])*X[2]*(2. + 49.*(-1. + X[2])*X[2]))))/Ri;
    dxdxp[3][3] = 2.0*M_PI;
  } 
  else if(defcoord == COMPLEX2TH) {
    dxdxp[1][1] = V[1]-R0;

    // now assign values
    if(X[2]<0.5){ myx2=X[2]; flip1=0.0; flip2=1.0;}
    else{ myx2=1.0-X[2]; flip1=M_PI; flip2=-1.0;}

    if(myx2<=x2trans){
      dxdxp[2][2] = (3.0*d2*pow(myx2,2.0)+2.0*c2*pow(myx2,1.0)+m2);
    }
    else{
      dxdxp[2][2] = (m3);
    }


    dxdxp[3][3] = 2.0*M_PI;


  }
  else if (defcoord == LOGRUNITH) {
    dxdxp[1][1] = V[1]-R0;
    dxdxp[2][2] = M_PI;
    dxdxp[3][3] = 2.0*M_PI;
  }
  else if (defcoord == JET1COORDS) {

    //dxdxp[1][1] = npow*(V[1]-R0)*pow(log(V[1]-R0),(npow-1.0)/npow);
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);

    myhslope1 = h0+dmyhslope1dr*(V[1]-rh0);
    myhslope2 = h0+dmyhslope2dx1*(X[1]-x1in);

    dmyhslopedx1=0.5*(dmyhslope1dr*dxdxp[1][1]+dmyhslope2dx1);
    myhslope=0.5*(myhslope1+myhslope2);

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    // d\theta/dx1 not 0
    // d\theta/dx1 = (d\theta/dr)*(dr/dx1)  (is this generally true or just when V[1](x1)?
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);

    dxdxp[3][3] = 2.0*M_PI;
  }
  else if (defcoord == JET2COORDS) {
    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);


    myhslope=2.0-pow(V[1]/r1jet,njet*(-1.0+exp(1.0)/exp(V[1]+rpjet)));

    dmyhslopedr=-((pow(exp(1.0),-V[1] - rpjet)*myhslope*njet*(-exp(1.0) + pow(exp(1.0),V[1] + rpjet) + exp(1.0)*V[1]*log(V[1]/r1jet)))/V[1]);
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);
    dxdxp[3][3] = 2.0*M_PI;
  }
  else if(defcoord == JET3COORDS){
    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);


    myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));

    dmyhslopedr=-((Qjet*(-((njet*(0.5 + atan(V[1]/r0jet - rsjet/r0jet)/M_PI))/V[1]) - (njet*r0jet*log(V[1]/r1jet))/(M_PI*(pow(r0jet,2) + pow(V[1] - rsjet,2)))))/pow(V[1]/r1jet,njet*(0.5 + atan(V[1]/r0jet - rsjet/r0jet)/M_PI)));

    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    if(X[2]<0.5){
      dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * mycos(2. * M_PI * X[2]);
      dxdxp[2][1] = -0.5*dmyhslopedx1* mysin(2. * M_PI * X[2]);
    }
    else{
      dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * mycos(2. * M_PI * (1.0-X[2]));
      dxdxp[2][1] = -0.5*dmyhslopedx1* (-mysin(2. * M_PI * (1.0-X[2])));
    }


    dxdxp[3][3] = 2.0*M_PI;

    

  }
  else if(defcoord == PULSARCOORDS){
    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);

    myhslope=(0.5+1.0/M_PI*atan((V[1]-rsjet)/r0jet))*(houter-hinner)+hinner;
    dmyhslopedr=(houter-hinner)*r0jet/(M_PI*(r0jet*r0jet+(V[1]-rsjet)*(V[1]-rsjet)));
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);



    dxdxp[3][3] = 2.0*M_PI;



    

  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    dxdxp[1][1] = Rout_array[1] - Rin_array[1];
    dxdxp[2][2] = Rout_array[2] - Rin_array[2];
    dxdxp[3][3] = Rout_array[3] - Rin_array[3];
    //    dualfprintf(fail_file,"COORD.c: dxdxp[1][1]=%21.15g\n",dxdxp[1][1]);
  }
  else if (defcoord == BILOGCYLCOORDS) {
    myexit(6666);
  }
  else if(defcoord==RAMESHCOORDS){


    //    V[1] = R0+exp(pow(X[1],npow)) ;
    //    myhslope=pow( (*r-rsjet)/r0jet , njet);
    //    V[2] = 0.5*M_PI*(1.0 + atan(myhslope*(X[2]-0.5))/atan(myhslope*0.5));


    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);
    // drdx2 = 0


    myhslope=pow( (V[1]-rsjet)/r0jet , njet);
    dmyhslopedr=(njet/r0jet)*pow( (V[1]-rsjet)/r0jet , njet-1.0);

    if(!finite(dmyhslopedr)){
      dualfprintf(fail_file,"Problem with dmyhslopedr=%g\n",dmyhslopedr);
      dualfprintf(fail_file,"njet=%g r=%g rsjet=%g r0jet=%g\n",njet,V[1],rsjet,r0jet);
      myexit(1);
    }

    // dhslope/dx1
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];


    // d\theta/dx2
    // (2*Pi*h(r))/(ArcTan(h(r)/2.)*(4 + Power(1 - 2*x2,2)*Power(h(r),2)))
    
    dxdxp[2][2] = (2.0*M_PI*myhslope)/(atan(myhslope*0.5)*(4.0 + pow(1.0 - 2.0*X[2],2.0)*pow(myhslope,2.0)));


    // d\theta/dr
    //(Pi*((-4*ArcTan((-0.5 + x2)*h(r)))/(4 + Power(h(r),2)) + 
    //       (4*(-1 + 2*x2)*ArcTan(h(r)/2.))/
    //        (4 + Power(1 - 2*x2,2)*Power(h(r),2)))*
    //     Derivative(1)(h)(r))/(4.*Power(ArcTan(h(r)/2.),2))

    
    // d\theta/dx1  = d\theta/dr dr/dx1
    dxdxp[2][1] =     (M_PI*dmyhslopedx1*(
					  (-4.0*atan((-0.5 + X[2])*myhslope))/(4. + pow(myhslope,2.)) + 
					  (4.*(-1. + 2.*X[2])*atan(myhslope/2.))/(4. + pow(1. - 2.*X[2],2.)*pow(myhslope,2.))
					  )
		       )/(4.*pow(atan(myhslope/2.),2.));

    dxdxp[3][3] = 2.0*M_PI;


  }

  else{
    dualfprintf(fail_file,"Shouldn't reach end of dxdxp\n");
    myexit(1);
  }




}



#define GAMMIEDERIVATIVE 0
#define NUMREC 1  // at present this is too slow since dxdxp is called many times.  Could setup permenant memory space, but kinda stupid

// Note that can't use NUMREC for connection if using numerical dxdxp.  Any other form of connection and then any dxdxp can be used.

// For example, if both connection and dxdxp are computed using GAMMIEDERIVATIVE, seems to work fine as a decent numerical approximation

// Also, one can use NUMREC for dxdxp if using analytic connection or GAMMIEDERIVATIVE connection.

//#define DXDERTYPE NUMREC 
#define DXDERTYPE GAMMIEDERIVATIVE

// see conn_func() for notes
#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
#define DXDELTA 1E-5
#elif(REALTYPE==LONGDOUBLETYPE)
#define DXDELTA 1E-6
#endif

void dxdxp_numerical(FTYPE *X, FTYPE (*dxdxp)[NDIM])
{
  int j,k,l;
  FTYPE Xh[NDIM], Xl[NDIM];
  FTYPE Vh[NDIM],Vl[NDIM];
  FTYPE blcoordsimple(FTYPE*X, int i, int j);
  extern FTYPE dfridr(FTYPE (*func)(FTYPE*,int,int), FTYPE *X,int ii, int jj, int kk);
  void donothing(volatile FTYPE *temp);
  volatile FTYPE temp;
  volatile FTYPE dxmachine[NDIM];


  if(DXDERTYPE==GAMMIEDERIVATIVE){

    for(k=0;k<NDIM;k++){
      for(j=0;j<NDIM;j++){

	if((j==TT)||(k==TT)){
	  if(j!=k) dxdxp[j][k]=0.0;
	  else dxdxp[j][k]=1.0;
	}
	else if((j==PH)||(k==PH)){
	  if(j!=k) dxdxp[j][k]=0.0;
	  else dxdxp[j][k]=1.0;
	}
	else{
	  for(l=0;l<NDIM;l++){
	    Xl[l]=Xh[l]=X[l]; // location of derivative
	    temp = X[l]-DXDELTA;
	    donothing(&temp);
	    X[l] = temp+DXDELTA;
	    temp = X[l]+DXDELTA;
	    donothing(&temp);
	    dxmachine[l] = temp-X[l];
	  }
	  
	  Xh[k]+=dxmachine[k]; // shift up
	  Xl[k]-=dxmachine[k]; // shift down

	  // below 2 lines redundant because gets both coordinates, but ok
	  bl_coord(Xh, Vh);
	  bl_coord(Xl, Vl);
	  dxdxp[j][k] = (Vh[j] - Vl[j]) / (Xh[k] - Xl[k]);
	  // GODMARK: unless N is a power of 2, X causes V to not be machine representable
	  // So even for a uniform grid dxdxp can vary near machine level
	  //  dualfprintf(fail_file,"(Vh[%d] - Vl[%d])=%21.15g (Xh[%d] - Xl[%d])=%21.15g\n",j,j,(Vh[j] - Vl[j]),k,k,(Xh[k] - Xl[k]));
	}
      }
    }

  }
  else if(DXDERTYPE==NUMREC){

    for(k=0;k<NDIM;k++) for(j=0;j<NDIM;j++){
      dxdxp[j][k]=dfridr(blcoordsimple,X,0,j,k);
    }

  }
}

#undef GAMMIEDERIVATIVE
#undef NUMREC
#undef DXDERTYPE
#undef DXDELTA

void donothing(volatile FTYPE *temp)
{
  *temp=*temp;

}

FTYPE blcoordsimple(FTYPE*X, int i, int j) // i not used
{
  FTYPE V[NDIM];
  
  if((j==TT)||(j==PH)) return(X[j]); // dummy linear relationship
  else{
    bl_coord(X, V);
    return(V[j]);
  }

}



// /////////////////////////////////////////////////////////////////
// 
// Below set X uniform grid -- usually doesn't change.
// Can usually force startx[1]=startx[2]=0. and dx[1]=1./N1 dx[2]=1./N2
// 
// /////////////////////////////////////////////////////////////////


/* some grid location, dxs */
// could find this by root finding.  Needed if no obvious bounds
// alternatively, could always define grid so x1=0..1 and x2=0..1 (likely more reasonable!)
void set_points()
{

  if (defcoord == LOGRSINTH) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } else if (defcoord == COMPLEX0TH) {
    startx[1] = 0.;
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = 1. / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } else if (defcoord == UNIRSINTH) {
    startx[1] = Rin;
    startx[2] = 0.324;
    startx[3] = 0.;
    dx[1] = (Rout-Rin) / totalsize[1];
    dx[2] = (1.0-2*.324) / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == EQMIRROR) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 0.5 / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == COMPLEX1TH) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == COMPLEX2TH) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == LOGRUNITH) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } 
  else if (defcoord == JET1COORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } 
  else if (defcoord == JET2COORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } 
  else if (defcoord == JET3COORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } 
  else if (defcoord == PULSARCOORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];

    dx[3] = 1.0 / totalsize[3];
    //    dx[3] = 2.0*M_PI;
  } 
  else if (defcoord == BILOGCYLCOORDS) {
    startx[1] = 0.;
    startx[2] = -1.0;
    startx[3] = 0.0;
    dx[1] = 1.0 / totalsize[1];
    dx[2] = 2.0 / totalsize[2];
    dx[3] = 1.0 / totalsize[3];
  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    startx[1] = 0.;
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = 1. / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1. / totalsize[3];
  }
  else if (defcoord == RAMESHCOORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1. / totalsize[3];
  } 
  else{
    dualfprintf(fail_file,"Shouldn't reach end of set_points\n");
    myexit(1);
  }


}

#define MAXIHOR 10
#define FRACN1 (0.1)
#define ADJUSTFRACT (0.25)

int setihor(void)
{
  // set to smaller of either totalsize[1]*0.1 or MAXIHOR
  if(totalsize[1]*FRACN1>MAXIHOR) return((int)MAXIHOR);
  else return((int)((FTYPE)totalsize[1]*(FTYPE)FRACN1));
}


// there's probably a way to do this in general
// probably can root find to get this

// set Rin so horizon exactly on FACE1 at i=ihor
FTYPE setRin(int ihor)
{
 
  FTYPE ftemp;
  FTYPE ihoradjust;


  ihoradjust=((FTYPE)ihor)+ADJUSTFRACT; // can't have grid edge exactly on horizon due to ucon_calc()

  //  fprintf(stderr,"ihoradjust = %21.15g\n",ihoradjust);

  if(defcoord == LOGRSINTH){
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == COMPLEX0TH){ // even though form appears different in X1, same Rin results
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == UNIRSINTH){ // uniform
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return((Rhor-ftemp*Rout)/(1.0-ftemp));
  }
  else if(defcoord == EQMIRROR){ // as defcoord=0
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == COMPLEX1TH){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == COMPLEX2TH){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == LOGRUNITH){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == JET1COORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord == JET2COORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord == JET3COORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord == PULSARCOORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    //no horizon in Cartesian coords -> do not set anything
    return( 1.0 );
  }
  else if(defcoord==RAMESHCOORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of setRin\n");
    myexit(1);
    return(-1);
  }

}



void coord(int i, int j, int k, int loc, FTYPE *X)
{

  if (loc == CENT) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE1) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE2) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE3) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN1) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN2) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN3) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == CORNT) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }

  X[0]=0.0;

  return;
}

// identical to coord except declaration using FTYPEs
void coordf(FTYPE i, FTYPE j, FTYPE k, int loc, FTYPE *X)
{

  if (loc == CENT) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE1) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE2) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE3) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN1) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN2) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN3) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == CORNT) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }

  X[0]=0.0;

  return;
}

void icoord(FTYPE *X,int loc, int *i, int *j, int *k)
{
  if(loc == CENT){
    *i = (int)((X[1]-startx[1])/dx[1] - 0.5) ;
    *j = (int)((X[2]-startx[2])/dx[2] - 0.5) ;
    *k = (int)((X[3]-startx[3])/dx[3] - 0.5) ;
  }

  if(startpos[1]+ *i>=totalsize[1]+N1BND) *i=totalsize[1]-1+N1BND;
  if(startpos[2]+ *j>=totalsize[2]+N2BND) *j=totalsize[2]-1+N2BND;
  if(startpos[3]+ *k>=totalsize[3]+N3BND) *k=totalsize[3]-1+N3BND;
  
  if(startpos[1]+ *i<-N1BND) *i=-N1BND;
  if(startpos[2]+ *j<-N2BND) *j=-N2BND;
  if(startpos[3]+ *k<-N3BND) *k=-N3BND;

}
