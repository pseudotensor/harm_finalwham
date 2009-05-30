
#include "decs.h"

void set_arrays()
{
  int i, j, k, pl, l, m;
  int ii;
  int floor,pf, tscale,dtstage;
  FTYPE valueinit;
  int dir,interpi,enodebugi;
  int dimen;


  p = (FTYPE (*)[N2M][N3M][NPR]) (&(a_p[N1BND][N2BND][N3BND][0]));
  panalytic = (FTYPE (*)[N2M][N3M][NPR]) (&(a_panalytic[N1BND][N2BND][N3BND][0]));
  omegafanalytic = (FTYPE (*)[N2M][N3M][NPR]) (&(a_omegafanalytic[N1BND][N2BND][N3BND][0]));


  emf = (FTYPE (*)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3]) (&(a_emf[-1][N1BND][N2BND][N3BND]));// inner shift still same
  vconemf = (FTYPE (*)[N2M][N3M][NDIM-1]) (&(a_vconemf[N1BND][N2BND][N3BND][-U1]));
#if(STOREWAVESPEEDS)
  wspeed = (FTYPE (*)[2][N1M][N2M][N3M]) (&(a_wspeed[-1][0][N1BND][N2BND][N3BND])); // shifted so wspeed[1,2,3] accesses the memory
#endif


  unew = (FTYPE (*)[N2M][N3M][NPR]) (&(a_unew[N1BND][N2BND][N3BND][0]));
  ulast = (FTYPE (*)[N2M][N3M][NPR]) (&(a_ulast[N1BND][N2BND][N3BND][0]));

#if(DOENOFLUXMEMORY)
  uinitial = (FTYPE (*)[N2M][N3M][NPR]) (&(a_uinitial[N1BND][N2BND][N3BND][0]));
  upoint = (FTYPE (*)[N2M][N3M][NPR]) (&(a_upoint[N1BND][N2BND][N3BND][0]));
  dUgeomarray = (FTYPE (*)[N2M][N3M][NPR]) (&(a_dUgeomarray[N1BND][N2BND][N3BND][0]));
#endif

  pleft = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_pleft[N1BND][N2BND][N3BND][0]));
  pright = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_pright[N1BND][N2BND][N3BND][0]));

  prc = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_prc[N1BND][N2BND][N3BND][0]));

#if(DODQMEMORY)
#if(N1>1)
  dq1 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_dq1[N1BND][N2BND][N3BND][0]));
#endif
#if(N2>1)
  dq2 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_dq2[N1BND][N2BND][N3BND][0]));
#endif
#if(N3>1)
  dq3 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_dq3[N1BND][N2BND][N3BND][0]));
#endif
#else
  dq1 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (0);
  dq2 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (0);
  dq3 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (0);
#endif

  ptemparray = (FTYPE (*)[N2M][N3M][NPR]) (&(a_ptemparray[N1BND][N2BND][N3BND][0]));
  
#if(N1>1)
  F1 = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F1[N1BND][N2BND][N3BND][0]));
#endif
#if(N2>1)
  F2 = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F2[N1BND][N2BND][N3BND][0]));
#endif
#if(N3>1)
  F3 = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F3[N1BND][N2BND][N3BND][0]));
#endif

#if( (FLUXDIMENSPLIT==UNSPLIT)&&(N1NOT1*N2NOT1*N3NOT1) )
  Fa = (FTYPE (*)[N2M][N3M][NPR]) (&(a_Fa[N1BND][N2BND][N3BND][0]));
  Fb = (FTYPE (*)[N2M][N3M][NPR]) (&(a_Fb[N1BND][N2BND][N3BND][0]));
  trifprintf("Allocating Fa/Fb for UNSPLIT flux method\n");
#else
  Fa = (FTYPE (*)[N2M][N3M][NPR]) (0);
  Fb = (FTYPE (*)[N2M][N3M][NPR]) (0);
#endif

#if( LIMIT_FLUXC2A_PRIM_CHANGE )
  fluxvectemp = (FTYPE (*)[N2M][N3M][NPR]) (&(a_fluxvectemp[N1BND][N2BND][N3BND][0]));	/* termporary storage for flux */  //atch
#endif


#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
  weno_lower_order_fraction = (FTYPE (*)[N2M][N3M][NPR]) (&(a_weno_lower_order_fraction[N1BND][N2BND][N3BND][0]));  //atch
  
  //init the array with zeroes
  FULLLOOP PLOOP(pl) weno_lower_order_fraction[i][j][k][pl] = 0.0;
#endif
  
#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )
  weno_prim_lower_order_fraction = (FTYPE (*)[N1M][N2M][N3M]) (&(a_weno_prim_lower_order_fraction[0][N1BND][N2BND][N3BND]));  //atch
  
  //init the array with zeroes
  FULLLOOP  DIMENLOOP(dimen) weno_prim_lower_order_fraction[dimen][i][j][k] = 0.0;
#endif



#if((A2CDIMENSPLIT==UNSPLIT)&&( (N1NOT1*N2NOT1*N3NOT1)||(N1NOT1&&N2NOT1)||(N1NOT1&&N3NOT1)||(N2NOT1&&N3NOT1) ) )
  a2cin  = (FTYPE (*)[N2M][N3M][NPR]) (&(a_Fa[N1BND][N2BND][N3BND][0]));
  a2cout = (FTYPE (*)[N2M][N3M][NPR]) (&(a_Fb[N1BND][N2BND][N3BND][0]));
  trifprintf("Allocating a2cin/a2cout for UNSPLIT a2c method\n");
#else
  a2cin  = (FTYPE (*)[N2M][N3M][NPR]) (0);
  a2cout = (FTYPE (*)[N2M][N3M][NPR]) (0);
#endif

#if(N1>1)
  //  F1CT = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F1CT[N1BND][N2BND][N3BND][0]));
#endif
#if(N2>1)
  //  F2CT = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F2CT[N1BND][N2BND][N3BND][0]));
#endif
#if(N3>1)
  //  F3CT = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F3CT[N1BND][N2BND][N3BND][0]));
#endif

  pk = (FTYPE (*)[N1M][N2M][N3M][NPR]) (&(a_pk[0][N1BND][N2BND][N3BND][0]));

  pflag = (int (*)[N2M][N3M][NUMPFLAGS]) (&(a_pflag[N1BND][N2BND][N3BND][0]));

#if(DOENODEBUG)
  enodebugarray = (CTYPE (*)[N2M][N3M][3-1][NUMINTERPTYPES][NPR-4][NUMENODEBUGS]) (&(a_enodebugarray[N1BND][N2BND][N3BND][0][0][0][0])); //SASMARKK (make dir count from 1 to 2 by changing 0 to -1)  //atch debug
#endif

#if(DODEBUG)
    failfloorcount = (CTYPE (*)[N2M][N3M][NUMTSCALES][NUMFAILFLOORFLAGS]) (&(a_failfloorcount[N1BND][N2BND][N3BND][0][0]));  //SASMARK: added another [0] on the right -- was a memory leak? atch
#endif


#if(DODISS)
    dissfunpos = (FTYPE (*)[N2M][N3M][NUMDISSFUNPOS]) (&(a_dissfunpos[N1BND][N2BND][N3BND][0]));
#endif


#if(CALCFARADAYANDCURRENTS)
  // this faraday needed for current calculation
  cfaraday =  (FTYPE (*)[N2M][N3M][NUMCURRENTSLOTS][3]) (&(a_cfaraday[N1BND][N2BND][N3BND][0][0]));
  fcon =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_fcon[N1BND][N2BND][N3BND][0]));
  jcon = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_jcon[N1BND][N2BND][N3BND][0]));
#endif

  // assume time average stuff gets zeroed in avg routine
#if(DOAVG)
  normalvarstavg =  (FTYPE (*)[N2M][N3M][NUMNORMDUMP]) (&(a_normalvarstavg[N1BND][N2BND][N3BND][0]));
  anormalvarstavg =  (FTYPE (*)[N2M][N3M][NUMNORMDUMP]) (&(a_anormalvarstavg[N1BND][N2BND][N3BND][0]));

  fcontavg =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_fcontavg[N1BND][N2BND][N3BND][0]));
  fcovtavg =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_fcovtavg[N1BND][N2BND][N3BND][0]));

  afcontavg =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_afcontavg[N1BND][N2BND][N3BND][0]));
  afcovtavg =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_afcovtavg[N1BND][N2BND][N3BND][0]));

  massfluxtavg =  (FTYPE (*)[N2M][N3M][NDIM]) (&(a_massfluxtavg[N1BND][N2BND][N3BND][0]));
  amassfluxtavg =  (FTYPE (*)[N2M][N3M][NDIM]) (&(a_amassfluxtavg[N1BND][N2BND][N3BND][0]));

  othertavg =  (FTYPE (*)[N2M][N3M][NUMOTHER]) (&(a_othertavg[N1BND][N2BND][N3BND][0]));
  aothertavg =  (FTYPE (*)[N2M][N3M][NUMOTHER]) (&(a_aothertavg[N1BND][N2BND][N3BND][0]));

#if(CALCFARADAYANDCURRENTS)
  jcontavg = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_jcontavg[N1BND][N2BND][N3BND][0]));
  jcovtavg = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_jcovtavg[N1BND][N2BND][N3BND][0]));

  ajcontavg = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_ajcontavg[N1BND][N2BND][N3BND][0]));
  ajcovtavg = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_ajcovtavg[N1BND][N2BND][N3BND][0]));
#endif

  tudtavg = (FTYPE (*)[N2M][N3M][NUMSTRESSTERMS]) (&(a_tudtavg[N1BND][N2BND][N3BND][0]));
  atudtavg = (FTYPE (*)[N2M][N3M][NUMSTRESSTERMS]) (&(a_atudtavg[N1BND][N2BND][N3BND][0]));
#endif  

#if(DOENODEBUG)
  FULLLOOP DIMENLOOP(dir) INTERPLOOP(interpi) PLOOP(pl) ENODEBUGLOOP(enodebugi){
    if(dir<=2 && pl<=U2){
      enodebugarray[i][j][k][dir-1][interpi][pl][enodebugi]=0;
    }
  }
#endif


  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  valueinit=sqrt(-1.0);
  FULLLOOP {
    if(DODEBUG) TSCALELOOP(tscale) FLOORLOOP(floor) failfloorcount[i][j][k][tscale][floor]=valueinit;
    PFLAGLOOP(pf) pflag[i][j][k][pf] = -100;

    PLOOPINTERP(pl){
      prc[i][j][k][pl] = valueinit;
#if(DODQMEMORY)
#if(N1>1)
      dq1[i][j][k][pl] = valueinit;
#endif
#if(N2>1)
      dq2[i][j][k][pl] = valueinit;
#endif
#if(N3>1)
      dq3[i][j][k][pl] = valueinit;
#endif
#endif
    }

    PLOOP(pl) {
#if(STOREWAVESPEEDS)
      for(l=1;l<=COMPDIM;l++) for(m=0;m<2;m++) wspeed[l][m][i][j][k] = valueinit;
#endif

      p[i][j][k][pl] = valueinit;
      panalytic[i][j][k][pl] = valueinit;
      omegafanalytic[i][j][k][pl] = valueinit;
      pleft[i][j][k][pl] = valueinit;
      pright[i][j][k][pl] = valueinit;
      unew[i][j][k][pl] = valueinit;
      DTSTAGELOOP(dtstage) pk[dtstage][i][j][k][pl] = valueinit;


#if(N1>1)
      F1[i][j][k][pl] = valueinit;
#endif
#if(N2>1)
      F2[i][j][k][pl] = valueinit;
#endif
#if(N3>1)
      F3[i][j][k][pl] = valueinit;
#endif

#if(N1>1)
      //      F1CT[i][j][k][pl] = valueinit;
#endif
#if(N2>1)
      //      F2CT[i][j][k][pl] = valueinit;
#endif
#if(N3>1)
      //      F3CT[i][j][k][pl] = valueinit;
#endif
    }
#if(CALCFARADAYANDCURRENTS)
    for(pl=0;pl<NUMCURRENTSLOTS;pl++) for(l=0;l<3;l++){
      cfaraday[i][j][k][pl][l]=valueinit;
    }
    for(pl=0;pl<NUMFARADAY;pl++){
      fcon[i][j][k][pl]=valueinit;
    }
    for(pl=0;pl<NDIM;pl++){
      jcon[i][j][k][pl]=valueinit;
    }
#endif

  }

#if(MCOORD!=CARTMINKMETRIC)
  /* grid functions */
  // GODMARK: for axisymmetric space-times, may want to keep metric functions 2D to save memory

  // these have 1 extra value on outer edge.  Shift for real pointer no different
  gcon = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
      (&(a_gcon[N1BND][N2BND][N3BND][0][0][0]));
  gcov = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
      (&(a_gcov[N1BND][N2BND][N3BND][0][0][0]));
  gcovpert = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
      (&(a_gcovpert[N1BND][N2BND][N3BND][0][0]));
  gdet = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
      (&(a_gdet[N1BND][N2BND][N3BND][0]));
  eomfunc = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
      (&(a_eomfunc[N1BND][N2BND][N3BND][0]));
  /*
  dxdxp = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
      (&(a_dxdxp[N1BND][N2BND][N3BND][0][0][0]));
  */

  // rest are always located at CENT
  conn = (FTYPE (*)[N2M][N3M][NDIM][NDIM][NDIM])
      (&(a_conn[N1BND][N2BND][N3BND][0][0][0]));
  conn2 = (FTYPE (*)[N2M][N3M][NDIM])
      (&(a_conn2[N1BND][N2BND][N3BND][0]));

#if(VOLUMEDIFF)
  idxvol = (FTYPE (*)[N2M][N3M][NDIM])(&(a_idxvol[N1BND][N2BND][N3BND][0]));
#endif

#else// MCOORD==CARTMINKMETRIC
  // pointers on left remain same, but if access "0" element of that pointer, should get really allocated memory.

  gcon = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
      (&(a_gcon[0][0][0][0][0][0]));
  gcov = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
      (&(a_gcov[0][0][0][0][0][0]));
  gcovpert = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
      (&(a_gcovpert[0][0][0][0][0]));
  gdet = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
      (&(a_gdet[0][0][0][0]));
  eomfunc = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
      (&(a_eomfunc[0][0][0][0]));
  /*
  dxdxp = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
      (&(a_dxdxp[0][0][0][0][0][0]));
  */

  // rest are always located at CENT
  conn = (FTYPE (*)[N2M][N3M][NDIM][NDIM][NDIM])
      (&(a_conn[0][0][0][0][0][0]));
  conn2 = (FTYPE (*)[N2M][N3M][NDIM])
      (&(a_conn2[0][0][0][0]));

#if(VOLUMEDIFF)
  idxvol = (FTYPE (*)[N2M][N3M][NDIM])(&(a_idxvol[0][0][0][0]));
#endif

#endif



#if(DOLUMVSR)
    // yes, for each cpu
    lumvsr=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(lumvsr==NULL){
      dualfprintf(fail_file,"Couldn't open lumvsr memory\n");
      myexit(1);
    }
  
    lumvsr_tot=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(lumvsr_tot==NULL){
      dualfprintf(fail_file,"Couldn't open lumvsr_tot memory\n");
      myexit(1);
    }
#endif
#if(DODISSVSR)
    // yes, for each cpu
    dissvsr=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(dissvsr==NULL){
      dualfprintf(fail_file,"Couldn't open dissvsr memory\n");
      myexit(1);
    }
  
    dissvsr_tot=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(dissvsr_tot==NULL){
      dualfprintf(fail_file,"Couldn't open dissvsr_tot memory\n");
      myexit(1);
    }
#endif
    //for(ii=0;ii<ncpux1*N1;ii++) dissvsr[ii]=0;
    //for(ii=0;ii<ncpux1*N1;ii++) dissvsr_tot[ii]=0;
}
