
#include "decs.h"


#define CORNGDETVERSION 1
// whether to specify gdet at end when setting EMF or to have internal to variables before averaging.
// point is that gdet at end is probably better, esp. at coordinate singularities.


int flux_ct(int stage, FTYPE pb[][N2M][N3M][NPR],FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR])
{
  int i, j, k, pl, l;
  // Below stuff is for Athena 1 and Athena 2
  FTYPE ucon[NDIM];
  struct of_state state;
  int dir;
  struct of_geom geomf1,geomf2,geomf3,geomc,geomco1,geomco2,geomco3;
  FTYPE cmax1,cmin1,cmax2,cmin2,ctop1,ctop2;
  int ignorecourant;
  // Gammie stuff
  FTYPE v1_av,w1,v2_av,w2,w_12 ;
  FTYPE rho_av,sqrtrho,b1_av,b2_av,va1,va2 ;
  FTYPE emfmm[NDIM],emfpm[NDIM],emfmp[NDIM],emfpp[NDIM],alpha[NDIM] ;
  FTYPE B1pp,B1pm,B1mp,B1mm;
  FTYPE B2pp,B2pm,B2mp,B2mm ;
  FTYPE B3pp,B3pm,B3mp,B3mm ;
  FTYPE U1pp,U1pm,U1mp,U1mm;
  FTYPE U2pp,U2pm,U2mp,U2mm ;
  FTYPE U3pp,U3pm,U3mp,U3mm ;
  //  FTYPE cms_func(FTYPE *prim_var) ;
  FTYPE B1d,B1u,B1l,B1r;
  FTYPE B2d,B2u,B2l,B2r;
  FTYPE B3d,B3u,B3l,B3r;
  FTYPE pbavg[NPR];
  FTYPE diffusiveterm[NDIM];
  FTYPE coefemf[NDIM];



  ///////////////////////
  //
  //  COMPUTE PRE-FUNCTIONS used by Athena method
  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)){
    // compute v^i

    // loop must go over EMFZLOOP's range minus 1 for i and j and k
    //    PREEMFZLOOP{
    FULLLOOP{ // GODMARK: could try to be more constrained
      // use ucon_calc() to get v^i 

      // EMF below is based upon averaging of zone-centered quantities, so use CENT here (i.e. not CORN)
      get_geometry(i, j, k, CENT, &geomc);
      MYFUN(ucon_calc(pb[i][j][k], &geomc, ucon),"fluxct.c:flux_ct()", "ucon_calc() dir=0", 1);


      // geom.g is \detg for EMF flux (geom.e is EOM factor for flux equation)
#if(CORNGDETVERSION)
      for(l=U1;l<=U3;l++) vconemf[i][j][k][l]=(ucon[l-U1+1]/ucon[TT]); // put in at end
#else
      for(l=U1;l<=U3;l++) vconemf[i][j][k][l]=(ucon[l-U1+1]/ucon[TT])*(geomc.e[l]);
#endif
    }
  }

#if(CORNGDETVERSION)

  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)||(FLUXB==FLUXCTTOTH)||(FLUXB==FLUXCD)){

    /////////////
    //
    // strip off geometry factor, which is added when computing the EMF
    //
    ////////////
    FULLLOOP{ // GODMARK: could try to be more constrained

#if(N1>1)
      get_geometry(i,j,k,FACE1,&geomf1);
      F1[i][j][k][B1]/=(geomf1.e[B1]);
      F1[i][j][k][B2]/=(geomf1.e[B2]);
      F1[i][j][k][B3]/=(geomf1.e[B3]);
#endif

#if(N2>1)
      get_geometry(i,j,k,FACE2,&geomf2);
      F2[i][j][k][B1]/=(geomf2.e[B1]);
      F2[i][j][k][B2]/=(geomf2.e[B2]);
      F2[i][j][k][B3]/=(geomf2.e[B3]);
#endif
    
#if(N3>1)
      get_geometry(i,j,k,FACE3,&geomf3);
      F3[i][j][k][B1]/=(geomf3.e[B1]);
      F3[i][j][k][B2]/=(geomf3.e[B2]);
      F3[i][j][k][B3]/=(geomf3.e[B3]);
#endif
    }
  }
  // don't need to put back geometry factor since in the end the EMF defines the Flux completely (except for FLUXB==FLUXCTHLL)
#endif



  // GODMARK: strange one has to do this.  Related to FLUXCT not reducing correctly for plane-parallel grid-aligned flows?
  if((N2>1)&&(N3>1)) coefemf[1]=0.25;
  else coefemf[1]=0.5; // or 0 if neither as controlled by below
  if((N1>1)&&(N3>1)) coefemf[2]=0.25;
  else coefemf[2]=0.5; // or 0 if neither as controlled by below
  if((N1>1)&&(N2>1)) coefemf[3]=0.25;
  else coefemf[3]=0.5; // or 0 if neither as controlled by below

  /////////////////////////
  //
  // COMPUTE EMF
  //
  //////////////////////////  

  if(FLUXB==FLUXCTHLL){
    // do nothing
  }
  else if((FLUXB==FLUXCTTOTH)||(FLUXB==ATHENA1)||(FLUXB==ATHENA2)){
    /* calculate EMFs */
    /* Toth approach: just average */

    // flux part is just average of same emf term at 4 different edge locations of (B^2 v^1 - B^1 v^2)
    //    EMFZLOOP{
    //LOOPINFP1{ // constrain or control better? GODMARK
    // B^i = \dF^{it}
    // E_i = - [ijk] v^j B^k  , such that (\detg B^i),t = - (\detg(B^i v^j - B^j v^i)),j = - (\detg [ijk] E_k),j = ([ijk] emf[k]),j
      
    // -> E_1 = v^3 B^2 - v^2 B^3
    // -> E_2 = v^1 B^3 - v^3 B^1
    // -> E_3 = v^2 B^1 - v^1 B^2

    // emf[i] = - \detg E_i

    // And notice that Fj[Bi] = \dF^{ij} = B^i v^j - B^j v^i , where j=dir

    // so:
    // emf_1 = B^3 v^2 - B^2 v^3 = F2[B3] or -F3[B2]
    // emf_2 = B^1 v^3 - B^3 v^1 = F3[B1] or -F1[B3]
    // emf_3 = B^2 v^1 - B^1 v^2 = F1[B2] or -F2[B1]

    // Notice only 6 independent ways.  The diagonal terms vanish (e.g. Fi[Bi]=0).

    // emf_i centered on corner of plane with normal direction i


    LOOPINFP1dir1full{
#if((N2>1)||(N3>1))
      emf[1][i][j][k] =
	coefemf[1] * (
#if(N2>1)
		      + F2[i][j][k][B3] + F2[i][j][km1][B3]
#endif
#if(N3>1)
		      - F3[i][j][k][B2] - F3[i][jm1][k][B2]
#endif
		      );
#else
      emf[1][i][j][k]=0.0; // not really 0, but differences in emf will be 0, and that's all that matters
#endif
    }



    LOOPINFP1dir2full{
#if((N1>1)||(N3>1))
      emf[2][i][j][k] =
	coefemf[2] * (
#if(N3>1)
		      + F3[i][j][k][B1] + F3[im1][j][k][B1]
#endif
#if(N1>1)
		      - F1[i][j][k][B3] - F1[i][j][km1][B3]
#endif
		      );
#else
      emf[2][i][j][k]=0.0; // not really 0, but differences in emf will be 0, and that's all that matters
#endif
    }



    LOOPINFP1dir3full{
#if((N1>1)||(N2>1))
      emf[3][i][j][k] =
	coefemf[3] * (
#if(N1>1)
		      + F1[i][j][k][B2] + F1[i][jm1][k][B2]
#endif
#if(N2>1)
		      - F2[i][j][k][B1] - F2[im1][j][k][B1]
#endif
		      );
#else
      emf[3][i][j][k]=0.0; // not really 0, but differences in emf will be 0, and that's all that matters
#endif
    }



#if(CORNGDETVERSION)// then tack on geometry


    FULLLOOP{
      
      get_geometry(i,j,k,CORN1,&geomf1);
      get_geometry(i,j,k,CORN2,&geomf2);
      get_geometry(i,j,k,CORN3,&geomf3);
      
      // obviously geom.e[B2] has to be equal to geom.e[B3] for this method
      emf[1][i][j][k] *=(geomf1.e[B2]);
      // obviously geom.e[B1] has to be equal to geom.e[B2] for this method
      emf[2][i][j][k] *=(geomf2.e[B1]);
      // obviously geom.e[B1] has to be equal to geom.e[B2] for this method
      emf[3][i][j][k] *=(geomf3.e[B1]);

    }


#endif

    //    }// end EMFZLOOP
  }// end if FLUXCT TOTH
  else if(FLUXB==FLUXCD){
    // Centered Difference (CD)
    // here emf is actually electric field
    // in Toth 2000, the F1=f^x F2=f^y
    // see Toth 2000 equation 28/29
    // 0.125 here takes care of larger differencing used in advance() in advance.c
    // all emf's end up at CENT
    //    EMFZLOOP{
    //    LOOPOUTFM1{ // constrain or control better? GODMARK



    // GODMARK: Why did Charles change sign in front of F's (see old code) ?  He changed it as if real sign such that emf[i] = \detg E_i but then later flux was wrong sign since he set flux=emf

    // GODMARK: should just compute F with diffusive term at CENT directly instead of averaging flux

    
#if((N2>1)||(N3>1))
    LOOPOUTFM1dir1full{
      emf[1][i][j][k] =
	0.125 * (
		 + F2[i][j][k][B3] + F2[i][jp1][k][B3]
		 - F3[i][j][k][B2] - F3[i][j][kp1][B2]
		 );
    }
#endif
#if((N1>1)||(N3>1))
    LOOPOUTFM1dir2full{
      emf[2][i][j][k] =
	0.125 * (
		 + F3[i][j][k][B1] + F3[i][j][kp1][B1]
		 - F1[i][j][k][B3] - F1[ip1][j][k][B3]
		 );
    }
#endif
#if((N1>1)||(N2>1))
    LOOPOUTFM1dir3full{
      emf[3][i][j][k] =
	0.125 * (
		 + F1[i][j][k][B2] + F1[ip1][j][k][B2]
		 - F2[i][j][k][B1] - F2[i][jp1][k][B1]
		 );
    }
#endif

    FULLLOOP{
#if(CORNGDETVERSION)// then tack on geometry
      
      // EMF's all at CENT
      get_geometry(i,j,k,CENT,&geomc);
      
      // obviously geom.e[B2] has to be equal to geom.e[B3] for this method
      emf[1][i][j][k] *=(geomc.e[B2]);
      // obviously geom.e[B1] has to be equal to geom.e[B3] for this method
      emf[2][i][j][k] *=(geomc.e[B1]);
      // obviously geom.e[B1] has to be equal to geom.e[B2] for this method
      emf[3][i][j][k] *=(geomc.e[B1]);
#endif
    }

    //    }// end EMFZLOOP
  }// end if FLUXCD






  ///////////////////////////////////////////////////////
  //
  // ADD DIFFUSIVE CORRECTIONS
  //
  ///////////////////////////////////////////////////////




  // add diffusive term
  // must come after FLUXCT
  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)){

    // Stone & Gardiner point out that Toth FLUXCT and FLUXCD not consistent with underlying integration algorithm for plane-parallel, grid-aligned flows.
    // fix is to change 0.25 to 0.5 and use a diffusive term as in ATHENA1

    /* Stone & Gardiner eq. 39 */
    // Charles Gammie (8/17/05) says this is simple, but repairs HARM defect that 2D does not reduce to 1D when waves are along coordinate lines


    LOOPINFP1{ // constrain or control better? GODMARK
      //    EMFZLOOP{
      // {emf}_i=-\epsilon_{ijk} v^i B^k


      // average of the below results in averaged emf located at corner (CORN)
      // same sign as F's (B^2 v^1 - B^1 v^2) with gdet built into vconemf

      // could remove gdet from vconemf and F1/F2 and put gdet at CORN for final result!  Avoids axis problems?
      // applies to above Toth version as well!
      // GODMARK

#if((N2>1)||(N3>1))
      // emf_1
      emfmp[1] = 
	pb[i  ][jm1][k  ][B3]*vconemf[i  ][jm1][k  ][U2] -
	pb[i  ][jm1][k  ][B2]*vconemf[i  ][jm1][k  ][U3] ;
      emfmm[1] = 
	pb[i  ][jm1][km1][B3]*vconemf[i  ][jm1][km1][U2] -
	pb[i  ][jm1][km1][B2]*vconemf[i  ][jm1][km1][U3] ;
      emfpm[1] = 
	pb[i  ][j  ][km1][B3]*vconemf[i  ][j  ][km1][U2] -
	pb[i  ][j  ][km1][B2]*vconemf[i  ][j  ][km1][U3] ;
      emfpp[1] = 
	pb[i  ][j  ][k  ][B3]*vconemf[i  ][j  ][k  ][U2] -
	pb[i  ][j  ][k  ][B2]*vconemf[i  ][j  ][k  ][U3] ;
#endif

#if((N1>1)||(N3>1))
      // emf_2
      emfmp[2] = 
	pb[im1][j  ][k  ][B1]*vconemf[im1][j  ][k  ][U3] -
	pb[im1][j  ][k  ][B3]*vconemf[im1][j  ][k  ][U1] ;
      emfmm[2] = 
	pb[im1][j  ][km1][B1]*vconemf[im1][j  ][km1][U3] -
	pb[im1][j  ][km1][B3]*vconemf[im1][j  ][km1][U1] ;
      emfpm[2] = 
	pb[i  ][j  ][km1][B1]*vconemf[i  ][j  ][km1][U3] -
	pb[i  ][j  ][km1][B3]*vconemf[i  ][j  ][km1][U1] ;
      emfpp[2] = 
	pb[i  ][j  ][k  ][B1]*vconemf[i  ][j  ][k  ][U3] -
	pb[i  ][j  ][k  ][B3]*vconemf[i  ][j  ][k  ][U1] ;
#endif

#if((N1>1)||(N2>1))
      // emf_3 (same sign as FLUXCT method as emf[3])
      emfmp[3] = 
	pb[im1][j  ][k  ][B2]*vconemf[im1][j  ][k  ][U1] -
	pb[im1][j  ][k  ][B1]*vconemf[im1][j  ][k  ][U2] ;
      emfmm[3] = 
	pb[im1][jm1][k  ][B2]*vconemf[im1][jm1][k  ][U1] -
	pb[im1][jm1][k  ][B1]*vconemf[im1][jm1][k  ][U2] ;
      emfpm[3] = 
	pb[i  ][jm1][k  ][B2]*vconemf[i  ][jm1][k  ][U1] -
	pb[i  ][jm1][k  ][B1]*vconemf[i  ][jm1][k  ][U2] ;
      emfpp[3] = 
	pb[i  ][j  ][k  ][B2]*vconemf[i  ][j  ][k  ][U1] -
	pb[i  ][j  ][k  ][B1]*vconemf[i  ][j  ][k  ][U2] ;
#endif

      for(l=1;l<=3;l++){
	diffusiveterm[l]= 0.25*(emfmp[l] + emfmm[l] + emfpm[l] + emfpp[l]);
      }

#if(CORNGDETVERSION)// then tack on geometry
      
      get_geometry(i,j,k,CORN1,&geomf1);
      get_geometry(i,j,k,CORN2,&geomf2);
      get_geometry(i,j,k,CORN3,&geomf3);
      
      // obviously geom.e[B2] has to be equal to geom.e[B3] for this method
      diffusiveterm[1] *=(geomf1.e[B2]);
      // obviously geom.e[B1] has to be equal to geom.e[B2] for this method
      diffusiveterm[2] *=(geomf2.e[B1]);
      // obviously geom.e[B1] has to be equal to geom.e[B2] for this method
      diffusiveterm[3] *=(geomf3.e[B1]);
#endif

      // now add diffusive term to emf
      // notice original emf multiplied by 2 to account for diffusive term being subtracted, so result is consistent
      for(l=1;l<=3;l++){
	emf[l][i][j][k] = 2.0*emf[l][i][j][k] - diffusiveterm[l];
      }


    }// end EMFZLOOP
  }// end ATHENA1
      







  // add another diffusive flux correction
  // must come after ATHENA1
  if(FLUXB==ATHENA2){

    if(LIMADJUST>0){
      dualfprintf(fail_file,"Cannot use Athena2 with limadjust since dq's not defined\n");
      myexit(11);
    }
    // GODMARK
    // should just use unrescaled p2interp stored, but need one for each direction.
    // should also use wave speeds from edges?

    /* Stone & Gardiner eq. 48 */
    // Charles Gammie (8/17/05) says does much better on flux loop advection test than ordinary HARM
    LOOPINFP1{ // constrain or control better? GODMARK
      //    EMFZLOOP{

      // dq1 and dq2 and dq3 are well-defined from fluxcalc() (both directions)
     

      // simple average of 1-D linear extrapolation here, where could use p2interp data saved from step_ch.c?  still available by this function call?


#if((N2>1)||(N3>1))
      // for emf_1

      // B2 located at FACE2 @ (k-1)
      B2d = 0.5*(
		 pb[i][jm1][km1][B2] + 0.5*dq2[i][jm1][km1][B2] +
		 pb[i][j  ][km1][B2] - 0.5*dq2[i][j  ][km1][B2]
		 ) ;

      // B2 located at FACE2 @ (k)
      B2u = 0.5*(
		 pb[i][jm1][k][B2] + 0.5*dq2[i][jm1][k][B2] +
		 pb[i][j  ][k][B2] - 0.5*dq2[i][j  ][k][B2]
		 ) ;

      // B3 located at FACE3 @ (j-1)
      B3l = 0.5*(
		 pb[i][jm1][km1][B3] + 0.5*dq3[i][jm1][km1][B3] +
		 pb[i][jm1][k  ][B3] - 0.5*dq3[i][jm1][k  ][B3]
		 ) ;
      // B3 located at FACE3 @ j
      B3r = 0.5*(
		 pb[i][j][km1][B3] + 0.5*dq3[i][j][km1][B3] +
		 pb[i][j][k  ][B3] - 0.5*dq3[i][j][k  ][B3]
		 ) ;

      // B2 for all centers around CORN1
      // B2[j - mp][k - mp]
      B2mm = pb[i][jm1][km1][B2] ;
      B2mp = pb[i][jm1][k  ][B2] ;
      B2pm = pb[i][j  ][km1][B2] ;
      B2pp = pb[i][j  ][k  ][B2] ;

      // B3 for all centers around CORN1
      // B3[j - mp][k - mp]
      B3mm = pb[i][jm1][km1][B3] ;
      B3mp = pb[i][jm1][k  ][B3] ;
      B3pm = pb[i][j  ][km1][B3] ;
      B3pp = pb[i][j  ][k  ][B3] ;

      // compute characteristic velocity -- only for Athena2 method


      // average pb to CORN1 for average phase speed there
      PLOOP(pl) pbavg[pl]=0.25*(pb[i][j][k][pl]+pb[i][jm1][k][pl]+pb[i][j][km1][pl]+pb[i][jm1][km1][pl]);
      get_geometry(i, j, k, CORN1, &geomco1); // used here and below emf's
      MYFUN(get_state(pbavg, &geomco1, &state),"step_ch.c:flux_ct()", "get_state()", 1);
      dir=2; MYFUN(vchar(pbavg, &state, dir, &geomco1, &cmax1, &cmin1,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=1", 1);
      dir=3; MYFUN(vchar(pbavg, &state, dir, &geomco1, &cmax2, &cmin2,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=2", 2);
      ctop1 = max(fabs(cmax1), fabs(cmin1));
      ctop2 = max(fabs(cmax2), fabs(cmin2));
      //      alpha=0.5*(ctop1+ctop2); // use average?
      //      alpha=max(ctop1,ctop2); // use maximum?
      // seems alpha can be arbitrary since 0 is ATHENA1

      //      alpha = dx1/dt ;	/* crude approx */


      // GODMARK: seems to have left/right and up/down asymmetry due to subtraction
      // if fabs were around each sutracted term, then would be ok (e.g. fabs(B1d-B1u))

      // notice that ctop1 and ctop2 have different "units", so cannot use with B2/B3 arbitrarily, must be consistent.
      diffusiveterm[1] =  0.125*(
				 +ctop1*(
					 + B2d - B2mm - B2u + B2mp
					 + B2d - B2pm - B2u + B2pp
					 )
				 +ctop2*(
					 + B3r - B3pm - B3l + B3mm
					 + B3r - B3pp - B3l + B3mp
					 )
				 ) ;
#endif

#if((N1>1)||(N3>1))
      // for emf_2

      // B3 located at FACE3 @ (i-1)
      B3d = 0.5*(
		 pb[im1][j][km1][B3] + 0.5*dq3[im1][j][km1][B3] +
		 pb[im1][j][k  ][B3] - 0.5*dq3[im1][j][k  ][B3]
		 ) ;

      // B3 located at FACE3 @ (i)
      B3u = 0.5*(
		 pb[i][j][km1][B3] + 0.5*dq3[i][j][km1][B3] +
		 pb[i][j][k  ][B3] - 0.5*dq3[i][j][k  ][B3]
		 ) ;

      // B3 located at FACE1 @ (k-1)
      B1l = 0.5*(
		 pb[im1][j][km1][B1] + 0.5*dq1[im1][j][km1][B1] +
		 pb[i  ][j][km1][B1] - 0.5*dq1[i  ][j][km1][B1]
		 ) ;
      // B1 located at FACE1 @ k
      B1r = 0.5*(
		 pb[im1][j][k][B1] + 0.5*dq1[im1][j][k][B1] +
		 pb[i  ][j][k][B1] - 0.5*dq1[i  ][j][k][B1]
		 ) ;

      // B3 for all centers around CORN1
      // B3[k - mp][i - mp]
      B3mm = pb[im1][j][km1][B3] ;
      B3mp = pb[i  ][j][km1][B3] ;
      B3pm = pb[im1][j][k  ][B3] ;
      B3pp = pb[i  ][j][k  ][B3] ;

      // B1 for all centers around CORN1
      // B1[k - mp][i - mp]
      B1mm = pb[im1][j][km1][B1] ;
      B1mp = pb[i  ][j][km1][B1] ;
      B1pm = pb[im1][j][k  ][B1] ;
      B1pp = pb[i  ][j][k  ][B1] ;

      // compute characteristic velocity -- only for Athena2 method


      // average pb to CORN1 for average phase speed there
      PLOOP(pl) pbavg[pl]=0.25*(pb[i][j][k][pl]+pb[im1][j][k][pl]+pb[i][j][km1][pl]+pb[im1][j][km1][pl]);
      get_geometry(i, j, k, CORN2, &geomco2); // used here and below emf's
      MYFUN(get_state(pbavg, &geomco2, &state),"step_ch.c:flux_ct()", "get_state()", 1);
      dir=3; MYFUN(vchar(pbavg, &state, dir, &geomco2, &cmax1, &cmin1,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=1", 1);
      dir=1; MYFUN(vchar(pbavg, &state, dir, &geomco2, &cmax2, &cmin2,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=2", 2);
      ctop1 = max(fabs(cmax1), fabs(cmin1));
      ctop2 = max(fabs(cmax2), fabs(cmin2));
      //      alpha=0.5*(ctop1+ctop2); // use average?
      //     alpha=max(ctop1,ctop2); // use maximum?
      // seems alpha can be arbitrary since 0 is ATHENA1

      //      alpha = dx1/dt ;	/* crude approx */


      // GODMARK: seems to have left/right and up/down asymmetry due to subtraction
      // if fabs were around each sutracted term, then would be ok (e.g. fabs(B3d-B3u))

      // notice that ctop1 and ctop2 have different "units", so cannot use with B2/B3 arbitrarily, must be consistent.
      diffusiveterm[2] =  0.125*(
				 +ctop1*(
					 + B3d - B3mm - B3u + B3mp
					 + B3d - B3pm - B3u + B3pp
					 )
				 +ctop2*(
					 + B1r - B1pm - B1l + B1mm
					 + B1r - B1pp - B1l + B1mp
					 )
				 ) ;
#endif

#if((N1>1)||(N2>1))
      // for emf_3

      // B1 located at FACE1 @ (j-1)
      B1d = 0.5*(
		 pb[im1][jm1][k][B1] + 0.5*dq1[im1][jm1][k][B1] +
		 pb[i  ][jm1][k][B1] - 0.5*dq1[i  ][jm1][k][B1]
		 ) ;

      // B1 located at FACE1 @ j
      B1u = 0.5*(
		 pb[im1][j][k][B1] + 0.5*dq1[im1][j][k][B1] +
		 pb[i  ][j][k][B1] - 0.5*dq1[i  ][j][k][B1]
		 ) ;

      // B2 located at FACE2 @ (i-1)
      B2l = 0.5*(
		 pb[im1][jm1][k][B2] + 0.5*dq2[im1][jm1][k][B2] +
		 pb[im1][j  ][k][B2] - 0.5*dq2[im1][j  ][k][B2]
		 ) ;
      // B2 located at FACE2 @ i
      B2r = 0.5*(
		 pb[i][jm1][k][B2] + 0.5*dq2[i][jm1][k][B2] +
		 pb[i][j  ][k][B2] - 0.5*dq2[i][j  ][k][B2]
		 ) ;

      // B1 for all centers around CORN3
      // B1[i - mp][j - mp]
      B1mm = pb[im1][jm1][k][B1] ;
      B1mp = pb[im1][j  ][k][B1] ;
      B1pm = pb[i  ][jm1][k][B1] ;
      B1pp = pb[i  ][j  ][k][B1] ;

      // B2 for all centers around CORN3
      // B2[i - mp][j - mp]
      B2mm = pb[im1][jm1][k][B2] ;
      B2mp = pb[im1][j  ][k][B2] ;
      B2pm = pb[i  ][jm1][k][B2] ;
      B2pp = pb[i  ][j  ][k][B2] ;

      // compute characteristic velocity -- only for Athena2 method


      // average pb to CORN3 for average phase speed there
      PLOOP(pl) pbavg[pl]=0.25*(pb[i][j][k][pl]+pb[i][jm1][k][pl]+pb[im1][j][k][pl]+pb[im1][jm1][k][pl]);
      get_geometry(i, j, k, CORN3, &geomco3); // used here and below emf's
      MYFUN(get_state(pbavg, &geomco3, &state),"step_ch.c:flux_ct()", "get_state()", 1);
      dir=1; MYFUN(vchar(pbavg, &state, dir, &geomco3, &cmax1, &cmin1,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=1", 1);
      dir=2; MYFUN(vchar(pbavg, &state, dir, &geomco3, &cmax2, &cmin2,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=2", 2);
      ctop1 = max(fabs(cmax1), fabs(cmin1));
      ctop2 = max(fabs(cmax2), fabs(cmin2));
      //      alpha=0.5*(ctop1+ctop2); // use average?
      //      alpha=max(ctop1,ctop2); // use maximum?
      // seems alpha can be arbitrary since 0 is ATHENA1

      //      alpha = dx1/dt ;	/* crude approx */


      // GODMARK: seems to have left/right and up/down asymmetry due to subtraction
      // if fabs were around each sutracted term, then would be ok (e.g. fabs(B1d-B1u))

      // notice that ctop1 and ctop2 have different "units", so cannot use with B2/B3 arbitrarily, must be consistent.
      diffusiveterm[3] =  0.125*(
				 +ctop1*(
					 + B1d - B1mm - B1u + B1mp
					 + B1d - B1pm - B1u + B1pp
					 )
				 +ctop2*(
					 + B2r - B2pm - B2l + B2mm
					 + B2r - B2pp - B2l + B2mp
					 )
				 ) ;
#endif




      /////////////////////////
      //
      // add geometry
      //
      ////////////////////////

      // must add geometry (no choice since above diffusive correction never had geometry)
      get_geometry(i,j,k,CORN1,&geomf1);
      get_geometry(i,j,k,CORN2,&geomf2);
      get_geometry(i,j,k,CORN3,&geomf3);
      
      // obviously geom.e[B2] has to be equal to geom.e[B3] for this method
      diffusiveterm[1] *=(geomf1.e[B2]);
      // obviously geom.e[B1] has to be equal to geom.e[B2] for this method
      diffusiveterm[2] *=(geomf2.e[B1]);
      // obviously geom.e[B1] has to be equal to geom.e[B2] for this method
      diffusiveterm[3] *=(geomf3.e[B1]);

      //////////////
      //
      // now add diffusive term to emf
      //
      //////////////
      for(l=1;l<=3;l++){
	emf[l][i][j][k]+= - diffusiveterm[l];
      }


    }// end EMF loop
  }// end if athena2

















  ////////////////////////////////////
  ///////////////////////////////////
  //
  // compute flux from EMF (signs must be consistent so that similar signed quantity in end)
  // e.g. flux terms originally involving B^2 v^1 - B^1 v^2 must come out with same sign at end
  //
  // This calculation takes care of case when EMF is not needed, so above need not worry about resetting values of emf to 0 or something
  //
  // we set flux to 0 when that dimension doesn't matter, but note that flux really is not 0.  Rather the differences in the flux are 0 across that direction.  However, this difference never enters since this flux is only used as flux differences.
  // These flux differences are *assumed* to vanish if that dimension has N=1.  However, a problem could be setup such that those differences would NOT be 0.
  // For example, a wedge out of an axisymmetric slice that has boundaries that are not symmetric around the equator.  This would have flux differences that would lead to changes in the quantity.  This would be a quasi-2D problem modelled in 1D.
  //
  //////////////////////////////////
  //////////////////////////////////

  //  dualfprintf(fail_file,"got here: FLUXB=%d\n",FLUXB);

  if(FLUXB==FLUXCTHLL){
    // no change

  }
  else if((FLUXB==FLUXCTTOTH)||(FLUXB==ATHENA2)||(FLUXB==ATHENA1)){
    /* rewrite EMFs as fluxes, after Toth */

    /////////////////////////////////////
    //
    // F1
    //
    ////////////////////////////////////

#if(N1>1)
    LOOPOUTFM1dir1full{ // constrain or control better? GODMARK
      //    F1CTZLOOP {

      // below line always true
      F1[i][j][k][B1] = 0.;
      F1[i][j][k][B2] = 0.5 * (emf[3][i][j][k] + emf[3][i][jp1][k]); // put emf3 back to FACE1
      F1[i][j][k][B3] = - 0.5 * (emf[2][i][j][k] + emf[2][i][j][kp1]); // put emf2 back to FACE1

    }
#endif

    /////////////////////////////////////
    //
    // F2
    //
    ////////////////////////////////////

#if(N2>1)
    LOOPOUTFM1dir2full{ // constrain or control better? GODMARK
      //    F2CTZLOOP {
      F2[i][j][k][B1] = - 0.5 * (emf[3][i][j][k] + emf[3][ip1][j][k]);
      // below line always true
      F2[i][j][k][B2] = 0.;
      F2[i][j][k][B3] = 0.5 * (emf[1][i][j][k] + emf[1][i][j][kp1]);

      
    }
#endif

    /////////////////////////////////////
    //
    // F3
    //
    ////////////////////////////////////

#if(N3>1)
    LOOPOUTFM1dir3full{ // constrain or control better? GODMARK
      //    F3CTZLOOP {

      F3[i][j][k][B1] = 0.5 * (emf[2][i][j][k] + emf[2][ip1][j][k]);
      F3[i][j][k][B2] = - 0.5 * (emf[1][i][j][k] + emf[1][i][jp1][k]);
      // below line always true
      F3[i][j][k][B3] = 0.;

    }
#endif

  } // end if FLUXCT
  else if(FLUXB==FLUXCD){

    // F's are emf's, where dF/dx is (F(i+1)-F(i-1))/dx
    // fluxes remain at center!

    /////////////////////////////////////
    //
    // F1
    //
    ////////////////////////////////////
#if(N1>1)
    //    F1CTZLOOP {
    FULLLOOP{ // SHOULD try to be more constrained GODMARK

      // always below line
      F1[i][j][k][B1] = 0.0;
      F1[i][j][k][B2] = emf[3][i][j][k];
      F1[i][j][k][B3] = -emf[2][i][j][k];

    }
#endif

    /////////////////////////////////////
    //
    // F2
    //
    ////////////////////////////////////
#if(N2>1)
    //    F2CTZLOOP {
    FULLLOOP{ // SHOULD try to be more constrained // GODMARK

      F2[i][j][k][B1] = -emf[3][i][j][k];
      // always below line
      F2[i][j][k][B2] = 0.;
      F2[i][j][k][B3] = emf[1][i][j][k];
    }
#endif

    /////////////////////////////////////
    //
    // F3
    //
    ////////////////////////////////////
#if(N3>1)
    //    F3CTZLOOP {
    FULLLOOP{ // SHOULD try to be more constrained // GODMARK

      F3[i][j][k][B1] = emf[2][i][j][k];
      F3[i][j][k][B2] = -emf[1][i][j][k];
      // always below line
      F3[i][j][k][B3] = 0.;


    }
#endif


  } // end if FLUXCD



  return(0);
}
