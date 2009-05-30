
#include "decs.h"
#include "reconstructeno.h"  

// monoindicator[0]: -1,0,1 for rough, ambiguous, and monotonic
// monoindicator[1]: whether set cell's left value (or central value, for a2c/c2a reconstruction)
// monoindicator[2]: whether set cell's right value
// yout[0][i] is the left interface value for c2e (centered value for a2c & c2a) for grid cell i
// yout[1][i] is the right interface value for c2e for grid cell i, does not make sense for a2c & c2a
void compute_monotonicity_line(
																 int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
																 int *minorder, int *maxorder, int *shift,   
																 FTYPE *shockindicator, FTYPE (*df)[NBIGM], 
																 FTYPE (*monoindicator)[NBIGM], FTYPE *yin,  FTYPE (*yout)[NBIGM] )
{ 

	int i; 

#if( DOMONOINTERP == JMONOINTERP ) 
  extern void compute_jmonotonicity_line(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *shockindicator, FTYPE (*df)[NBIGM],  FTYPE (*monoindicator)[NBIGM] , FTYPE *yin, FTYPE (*yout)[NBIGM]);
#elif( DOMONOINTERP == SMONOINTERP ) 
	#include "smonointerp.h"
#elif( DOMONOINTERP == NOMONOINTERP )
#else
	#error Unknown monointerp type
#endif

#if( DOMONOINTERP == JMONOINTERP ) 
	compute_jmonotonicity_line(
																 recontype, whichreduce, preforder, pl, bs, ps, pe, be, 
																 minorder, maxorder, shift,   
																 shockindicator, df, 
																 monoindicator, yin,  yout );
#elif( DOMONOINTERP == SMONOINTERP ) 
	compute_smonotonicity_line(
																 recontype, whichreduce, preforder, pl, bs, ps, pe, be, 
																 minorder, maxorder, shift,   
																 shockindicator, df, 
																 monoindicator, yin,  yout );
#elif( DOMONOINTERP == NOMONOINTERP )
//no monointerp used; so set the monotonicity indicators to zero
	if(recontype!=CVT_C2E){
    for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
      monoindicator[0][i]=0;
      monoindicator[1][i]=0;
    }
  }
  else{
    for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
      monoindicator[0][i]=0;
      monoindicator[1][i]=0;
      monoindicator[2][i]=0;
    }
  }
#else
	#error Unknown monointerp type
#endif

}
