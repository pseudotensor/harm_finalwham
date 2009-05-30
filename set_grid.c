
#include "decs.h"

// GODMARK: may want to make grid functions explicitly 2D for axisymmetric space-times when in axisymmetry with space-time axis aligned with grid.

/* set up all grid functions */
void set_grid()
{
  int i, j, k, l, m;
  int ii, jj, kk;
  FTYPE X[NDIM];
  struct of_geom geom;
  int loc;
  FTYPE V[NDIM];
  extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  extern void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
  extern void eomfunc_func(int getprim, int whichcoord, FTYPE *X, FTYPE *eomfunc);

  /* set up boundaries, steps in coordinate grid */
  set_points();
  dV = dx[1] * dx[2] * dx[3]; // computational volume
  dVF = dV ; // full 3d volume (used for diagnostics only)



  DLOOPA(j) X[j] = 0.;


  ///////////////////
  //
  // Grid functions that only exist at many locations and are assigned
  // values on all points INCLUDING another value at the outer edges
  // so have edge grid data there -- makes setting up certain things
  // easier
  //
  // Notice that coord() and bl_coord() work without this.  So those
  // functions that only require those functions can do without this
  // extra grid stuff.
  //
  //////////////////

#if(MCOORD!=CARTMINKMETRIC)
  FULLLOOPP1
#else
    // doesn't depend on position, only store/use 1 value
  i=j=k=0;
#endif
  {
    

    // over grid locations needing these quantities
    for (loc = NPG - 1; loc >= 0; loc--) {

      coord(i, j, k, loc, X);

      /////////////////
      //
      // (1,MCOORD) here actually means PRIMCOORDS since the "1" means convert MCOORD to PRIMCOORDS.

      //      dxdxp_func(X,dxdxp[i][j][k][loc]); // future numerical version
      gcov_func(1,MCOORD,X, gcov[i][j][k][loc],gcovpert[i][j][k][loc]);
      gdet[i][j][k][loc] = gdet_func(gcov[i][j][k][loc]);
      gcon_func(1,MCOORD,X,gcov[i][j][k][loc],gcon[i][j][k][loc]);
      eomfunc_func(1,MCOORD,X,&eomfunc[i][j][k][loc]);

      // check if near static limit since can't divide by the below in ucon_calc
      // GODMARK
      if (fabs(gcon[i][j][k][loc][TT][TT]) < SMALL) {
	bl_coord(X,V);
	dualfprintf(fail_file, "grid location too near g_{tt}==0: %d %d %d : r=%21.15g th=%21.15g phi=%21.15g : Rin=%21.15g %21.15g\n", i,j,k,V[1],V[2],V[3],Rin,gcon[i][j][k][loc][TT][TT]);
	myexit(1);
      }
      if (fabs(gcon[i][j][k][loc][RR][RR]) < SMALL) {
	bl_coord(X,V);
	dualfprintf(fail_file, "grid location too near g^{rr}==0:  %d %d %d : r=%21.15g th=%21.15g phi=%21.15g :  Rin=%21.15g %21.15g\n", i,j,k,V[1],V[2],V[3],Rin,gcon[i][j][k][loc][RR][RR]);
	myexit(1);
      }
      if (fabs(gcon[i][j][k][loc][TH][TH]) < SMALL) {
	bl_coord(X,V);
	dualfprintf(fail_file,"grid location too near g^{\\theta\\theta}==0:  %d %d %d : r=%21.15g th=%21.15g phi=%21.15g :  Rin=%21.15g %21.15g\n", i,j,k,V[1],V[2],V[3],Rin,gcon[i][j][k][loc][TH][TH]);
	myexit(1);
      }
      if (fabs(gcon[i][j][k][loc][PH][PH]) < SMALL) {
	bl_coord(X,V);
	dualfprintf(fail_file,"grid location too near g^{\\phi\\phi}==0:  %d %d %d : r=%21.15g th=%21.15g phi=%21.15g :  Rin=%21.15g %21.15g\n", i,j,k,V[1],V[2],V[3],Rin,gcon[i][j][k][loc][PH][PH]);
	myexit(1);
      }
      // what about g_{tt}==0? Do I ever divide by g_{tt}?
      // Yes, for ucon[TT] for 4 velocity, which is done if WHICHVEL==VEL4 or init.c
      // what about g^{rr}==0? Do I ever divide by g^{rr}?
      // Doesn't appear so
      // g^{pp} inf taken care of in metric.c by avoiding theta=0,Pi

    }
  }


  ///////////////////
  //
  // Grid functions that only exist at one location for all grid points
  //
  //////////////////
#if(MCOORD!=CARTMINKMETRIC)
  FULLLOOP // connection only needed at center, and only has memory on full normal grid (not like gcon/gcov that have extra upper edge)
#else
    // doesn't depend on position, only store/use 1 value
  i=j=k=0;
#endif
  {

    // over grid locations needing these quantities
    for (loc = NPG - 1; loc >= 0; loc--) {

      coord(i, j, k, loc, X);

      if (loc == CENT) {
	get_geometry(i, j, k, loc, &geom);
	conn_func(MCOORD,X, &geom, conn[i][j][k],conn2[i][j][k]);
      }
    }
  }

  ///////////////////
  //
  // Grid functions that only exist at one location AND only on active grid
  //
  //////////////////
  
  if(VOLUMEDIFF){
#if(MCOORD!=CARTMINKMETRIC)
  ZLOOP
#else
    // doesn't depend on position, only store/use 1 value
  i=j=k=0;
#endif
  {
      // only at 1 location, centered, using surrounding edge values
      if((defcoord==LOGRUNITH)&&(MCOORD==KSCOORDS)){
	mks_unitheta_idxvol_func(i,j,k,idxvol[i][j][k]);
      }
      else{
	idxvol[i][j][k][TT]=1.0; // really 1/dt, but changes in time      
	idxvol[i][j][k][RR]=1.0/dx[1];
	idxvol[i][j][k][TH]=1.0/dx[2];
	idxvol[i][j][k][PH]=1.0/dx[3];
      }
    }
  }

  /* done! */
}
