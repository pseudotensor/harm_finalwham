
#include "decs.h"

// this is the presently used function

int bound_pflag_user(int boundstage, int prim[][N2M][N3M][NUMPFLAGS])
{
  int i, j, k, pl;
  int failreturn;
  int ri, rj, rk; // reference i,j,k

  // includes all conditions from global.h except NSSURFACE and FIXED

  // OUTFLOW: direct copy here since value is discrete and has discrete meaning
  // that is, properties are just copied.

  /* if fixed BC: do nothing */
  // GODMARK: no, since doing multi-steps, one should assign some fixed primitive to prim, as done for real primitives
  // general assume flag doesn't care about FIXED, since represents a failure.  One's own FIXED conditions shouldn't fail.


  /////////////////////
  //
  // PERIODIC
  //
  /////////////////////

  // periodic x1
  if ( (mycpupos[1] == 0)&&(mycpupos[1] == ncpux1 - 1) ) {
    if( (BCtype[X1DN]==PERIODIC)&&(BCtype[X1UP]==PERIODIC) ){
      // just copy from one side to another
	
      LOOPF2 LOOPF3{

	// copy from upper side to lower boundary zones
	ri=N1;
	rj=j;
	rk=k;
	LOOPBOUND1IN FLOOP(pl) prim[i][j][k][pl] = prim[ri+i][rj][rk][pl];

	// copy from lower side to upper boundary zones
	ri=0;
	rj=j;
	rk=k;
	LOOPBOUND1OUT FLOOP(pl) prim[i][j][k][pl] = prim[ri+(i-N1)][rj][rk][pl];
      }
    }
  }

  // periodic x2
  if ( (mycpupos[2] == 0)&&(mycpupos[2] == ncpux2 - 1) ) {
    if( (BCtype[X2DN]==PERIODIC)&&(BCtype[X2UP]==PERIODIC) ){
      // just copy from one side to another
	
      LOOPF1 LOOPF3{

	// copy from upper side to lower boundary zones
	ri=i;
	rj=N2;
	rk=k;
	LOOPBOUND2IN FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj+j][rk][pl];

	// copy from lower side to upper boundary zones
	ri=i;
	rj=0;
	rk=k;
	LOOPBOUND2OUT FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj+(j-N2)][rk][pl];
      }
    }
  }

  // periodic x3
  if ( (mycpupos[3] == 0)&&(mycpupos[3] == ncpux3 - 1) ) {
    if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){
      // just copy from one side to another
	
      LOOPF1 LOOPF2{

	// copy from upper side to lower boundary zones
	ri=i;
	rj=j;
	rk=N3;
	LOOPBOUND3IN FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk+k][pl];

	// copy from lower side to upper boundary zones
	ri=i;
	rj=j;
	rk=0;
	LOOPBOUND3OUT FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk+(k-N3)][pl];
      }
    }
  }

  /////////////////////
  //
  // OUTFLOW/FIXEDOUTFLOW (for events, just copy)
  //
  /////////////////////


  // outflow inner x1
  if (mycpupos[1] == 0) {
    if( (BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)||(BCtype[X1DN]==OUTFLOWNOINFLOW) ){
      /* inner r boundary condition: u, just copy */
      LOOPF2 LOOPF3{
	ri=0;
	rj=j;
	rk=k;
	LOOPBOUND1IN FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
      }
    }
  }

  // outflow outer x1
  if (mycpupos[1] == ncpux1 - 1) {
    if( (BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)||(BCtype[X1UP]==OUTFLOWNOINFLOW) ){
      /* outer r BC: outflow */
      LOOPF2 LOOPF3{
	ri=N1-1;
	rj=j;
	rk=k;
	LOOPBOUND1OUT FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
      }
    }
  }

  // outflow inner x2
  if (mycpupos[2] == 0) {
    if( (BCtype[X2DN]==OUTFLOW)||(BCtype[X2DN]==FIXEDOUTFLOW)||(BCtype[X2DN]==OUTFLOWNOINFLOW) ){
      /* inner r boundary condition: u, just copy */
      LOOPF1 LOOPF3{
	ri=i;
	rj=0;
	rk=k;
	LOOPBOUND2IN FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
      }
    }
  }

  // outflow outer x2
  if (mycpupos[2] == ncpux2 - 1) {
    if( (BCtype[X2UP]==OUTFLOW)||(BCtype[X2UP]==FIXEDOUTFLOW)||(BCtype[X2UP]==OUTFLOWNOINFLOW) ){
      /* outer r BC: outflow */
      LOOPF1 LOOPF3{
	ri=i;
	rj=N2-1;
	rk=k;
	LOOPBOUND2OUT FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
      }
    }
  }

  // outflow inner x3
  if (mycpupos[3] == 0) {
    if( (BCtype[X3DN]==OUTFLOW)||(BCtype[X3DN]==FIXEDOUTFLOW)||(BCtype[X3DN]==OUTFLOWNOINFLOW) ){
      /* inner r boundary condition: u, just copy */
      LOOPF1 LOOPF2{
	ri=i;
	rj=j;
	rk=0;
	LOOPBOUND3IN FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
      }
    }
  }

  // outflow outer x3
  if (mycpupos[3] == ncpux3 - 1) {
    if( (BCtype[X3UP]==OUTFLOW)||(BCtype[X3UP]==FIXEDOUTFLOW)||(BCtype[X3UP]==OUTFLOWNOINFLOW) ){
      /* outer r BC: outflow */
      LOOPF1 LOOPF2{
	ri=i;
	rj=j;
	rk=N3-1;
	LOOPBOUND3OUT FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
      }
    }
  }


  /////////////////////
  //
  // POLARAXIS/SYMM/ASYMM (for events (not values) these are the same)
  //
  /////////////////////


  // symmetry on inner x1
  if (mycpupos[1] == 0) {
    if( (BCtype[X1DN]==POLARAXIS)||(BCtype[X1DN]==SYMM)||(BCtype[X1DN]==ASYMM)) {
      LOOPF2 LOOPF3{
	ri=0;
	rj=j;
	rk=k;
	// symmetric copy
	LOOPBOUND1IN FLOOP(pl)  prim[i][j][k][pl] = prim[ri+(ri-i-1)][rj][rk][pl];
      }
    }
  }

  // symmetry on outer x1
  if (mycpupos[1] == ncpux1 - 1) {
    if( (BCtype[X1UP]==POLARAXIS)||(BCtype[X1UP]==SYMM)||(BCtype[X1UP]==ASYMM)) {
      LOOPF2 LOOPF3{
	ri=N1-1;
	rj=j;
	rk=k;
	LOOPBOUND1OUT FLOOP(pl)  prim[i][j][k][pl] = prim[ri+(ri-i+1)][rj][rk][pl];
      }
    }
  }


  // symmetry on inner x2
  if (mycpupos[2] == 0) {
    if( (BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM)) {
      LOOPF1 LOOPF3{
	ri=i;
	rj=0;
	rk=k;
	// symmetric copy
	LOOPBOUND2IN FLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj+(rj-j-1)][rk][pl];
      }
    }
  }

  // symmetry on outer x2
  if (mycpupos[2] == ncpux2 - 1) {
    if( (BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM)) {
      LOOPF1 LOOPF3{
	ri=i;
	rj=N2-1;
	rk=k;
	LOOPBOUND2OUT FLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj+(rj-j+1)][rk][pl];
      }
    }
  }

  // symmetry on inner x3
  if (mycpupos[3] == 0) {
    if( (BCtype[X3DN]==POLARAXIS)||(BCtype[X3DN]==SYMM)||(BCtype[X3DN]==ASYMM)) {
      LOOPF1 LOOPF2{
	ri=i;
	rj=j;
	rk=0;
	// symmetric copy
	LOOPBOUND3IN FLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj][rk+(rk-k-1)][pl];
      }
    }
  }

  // symmetry on outer x3
  if (mycpupos[3] == ncpux3 - 1) {
    if( (BCtype[X3UP]==POLARAXIS)||(BCtype[X3UP]==SYMM)||(BCtype[X3UP]==ASYMM)) {
      LOOPF1 LOOPF2{
	ri=i;
	rj=j;
	rk=N3-1;
	LOOPBOUND3OUT FLOOP(pl)  prim[i][j][k][pl] = prim[ri][rj][rk+(rk-k+1)][pl];
      }
    }
  }



  return (0);
}
