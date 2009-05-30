
/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
			    void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int), 
			    FTYPE (*res_func) (FTYPE []) )
{
  FTYPE f, f_old, df, df_old, dx[NEWT_DIM], dx_old[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, errx_old, errx_oldest, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE randtmp, tmp;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;
  FTYPE resid_norm, resid_check, grad_check;

  FTYPE res_func_val, res_func_old, res_func_new;
  FTYPE dn[NEWT_DIM], del_f[NEWT_DIM];

  static void my_lnsrch(int, FTYPE [], FTYPE, FTYPE [], FTYPE [], FTYPE [], FTYPE *, 
			  FTYPE, FTYPE, int *, FTYPE (*res_func) (FTYPE []));

  static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  ) ;

  int   keep_iterating, i_increase, retval2,retval = 0;
  const int ltrace  = 0;
  const int ltrace2 = 1;


  retval = 0;


  errx = 1. ; 
  errx_old = 2.;
  df = df_old = f = f_old = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;


  vsq_old = vsq = W = W_old = 0.;


  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 
     nstroke++;
     lntries++;

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */
      

#if(!OPTIMIZED)
    /*  Check for bad untrapped divergences : */
    if( (finite(f)==0) ||  (finite(df)==0) ) {
      if( debugfail >= 2 ) { 
	dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
	dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
      }
      return(1);
    }
#endif


#if(!OPTIMIZED)
    /* Randomly rescale Newton step to break out of iteration cycles: */
    if( ((n_iter+1) % CYCLE_BREAK_PERIOD) == 0 ) {
      randtmp = ( (1.*rand())/(1.*RAND_MAX) );
      for( id = 0; id < n ; id++) dx[id] *= randtmp;
      //	for( id = 0; id < n ; id++) dx[id] *= ( (1.*rand())/(1.*RAND_MAX) );
    }
#endif

    /* Save old values before calculating the new: */
    errx_oldest = errx_old;
    errx_old = errx;
    lerrx=errx;
    errx = 0.;
    f_old = f;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    if( do_line_search == 1 ) { 

      /* Compare the residual to its initial value */ 
      if( n_iter == 0 ) { 
	resid_norm = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_norm += fabs(resid[id]);
	}
	resid_norm /= 1.0*n ;
	if( resid_norm == 0.0 ) resid_norm = 1.0;
      }

      for( id = 0; id < n ; id++) {
	tmp = 0.;
	for( jd = 0; jd < n ; jd++) {
	  tmp += jac[jd][id] * resid[jd];
	}
	del_f[id] = tmp;
      }
      for( id = 0; id < n ; id++) {
	dn[id] = dx[id];
      }

      my_lnsrch(n, x_old-1, f_old, del_f-1, dn-1, x-1, &f, TOL_LINE_STEP, SCALEMAX, &retval, res_func);

      /* dx is needed for errx calculation below: */
      for( id = 0; id < n ; id++) {
	dx[id] = x[id] - x_old[id];
      }

#if(!OPTIMIZED)
      if( ltrace ) { 
	res_func_val = res_func(x);
	res_func_old = res_func(x_old);
	dualfprintf(fail_file,"gnr(): f_old, f, res_func_old, res_func_val = %21.15g  %21.15g  %21.15g  %21.15g  \n",
		f_old, f, res_func_old, res_func_val );
	dualfprintf(fail_file,"gnr(): x_old = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",x_old[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): x     = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",x[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): dn    = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",dn[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): del_f = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",del_f[id]);
	}
	dualfprintf(fail_file,"\n ");
      }
#endif

      /* Check to see if line search problem is because the residual vector is already small enough */
      if( retval == 1 ) {
	resid_check = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_check += fabs(resid[id]);
	}
	resid_check /= 1.0*n;
	
	if( resid_check <= resid_norm * NEWT_FUNC_TOL ) {
	  retval = 0;
	}
	if( ltrace && retval ) { 
	  dualfprintf(fail_file,"general_newton_raphson():  retval, resid_check = %4i  %21.15g \n",retval, resid_check);
	}	  
      }
      /* If initial Newton step is bad, then try again without line searching: */
      if( (retval == 2) && (USE_LINE_SEARCH == do_line_search) ) { 
#if(!OPTIMIZED)	  
	if( ltrace ) { 
	  dualfprintf(fail_file,"gnr(): bad first step: retval, f_old, f  = %4i  %21.15g  %21.15g  \n",retval,f_old,f);
	  dualfprintf(fail_file,"gnr: doing recursive call, retval, errx = %4i  %21.15g \n", retval, errx );
	}
#endif
	retval = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
	for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
	return( retval );
      }

      /* Check to see if it is trapped in a local minimum, i.e. gradient is too small */ 
      if( retval == 1 ) { 
	grad_check = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_check = (x[id] == 0.) ? 1.0 : fabs(x[id]) ;
	  grad_check  +=  del_f[id] * resid_check ;
	}
	resid_check = (f == 0.) ? 1.0 : fabs(f) ;
	grad_check /= resid_check;
	
	/* Then we've most likely found a solution: */
	if( grad_check > GRADMIN ) { 
	  retval = -1;
	}
	else if( ltrace ) { 
	  dualfprintf(fail_file,"general_newton_raphson():  retval, grad_check = %4i  %21.15g \n",retval, grad_check);
	}
      }
    }
    else {
      /* don't use line search : */
      for( id = 0; id < n ; id++) {
	x[id] += dx[id]  ;
      }

    } /* End of "to do line search or not to do line search..." */


    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific: 

#if( NEWCONVERGE == 1 )
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);

    /* For the old criterion, look at errors in each indep. variable(s) (except for 5D) : */
#else
    for( id = 0; id < n ; id++) {
      errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]);
    }
    errx /= 1.*n;
#endif


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    // METHOD specific:

    validate_x( x, x_old ) ;


    /****************************************/
    /* Check to see if we're in a infinite loop with error function: */
    /****************************************/
#if( CHECK_FOR_STALL )
    if( ( (errx_old == errx) || (errx_oldest == errx) ) && (errx <= MIN_NEWT_TOL) )  errx = -errx;
#endif 

    /****************************************/
    /* If there's a problem with line search, then stop iterating: */
    /****************************************/
    if( (retval == 1) || (retval == -1) ) errx = -errx;


#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file," general_newton_raphson(): niter,f_old,f,errx_old,errx = %4i  %21.15g  %21.15g  %21.15g  %21.15g\n",  
	      n_iter,f_old,f,errx_old,errx );
      dualfprintf(fail_file,"gnr(): x_old = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",x_old[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): x     = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",x[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): dx     = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",dx[id]);
      }
      dualfprintf(fail_file,"\n ");
      
    }
#endif

    /****************************************/
    /* Prepare for the next iteration, set the "old" variables: */
    /****************************************/
    for( id = 0; id < n ; id++)  dx_old[id] = dx[id] ;
    f_old  = f;
    df_old = df;


    /****************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations before stopping */
    /****************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

#if(CRAZYDEBUG)
    if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"n_iter=%d errx=%21.15g %21.15g\n",n_iter,errx,MIN_NEWT_TOL);
    }
#endif

  }   // END of while(keep_iterating)

  
    /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) ||  (finite(df)==0) || (finite(x[0])==0) || (finite(x[1])==0)) {
#if(!OPTIMIZED)
    if( debugfail >= 2 ) { 
      dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
      dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
    }
#endif
    return(1);
  }


  if( fabs(errx) > MIN_NEWT_TOL){
    if( (do_line_search != USE_LINE_SEARCH) || (USE_LINE_SEARCH < 0) ) { 
#if(DOHISTOGRAM)
      bin_newt_data( errx, n_iter, 0, 0 );
#endif

#if(!OPTIMIZED)
      if(ltrace2) {
	dualfprintf(fail_file," totalcount = %d   0   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      }
      if(ltrace) {
	dualfprintf(fail_file,"general_newton_raphson():  did not find solution \n");
	if( retval == -1 ) {
	  dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x = ");
	  for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x[id]);
	  dualfprintf(fail_file,"\n");
	  dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x_old = ");
	  for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x_old[id]);
	  dualfprintf(fail_file,"\n");
	}
	
      }
      //      dualfprintf(fail_file,"gnr retval2 = %4i \n", 1); 
#endif
      return(1);
    } 
    else {
      /* If bad return and we tried line searching, try it without before giving up: */
      //      dualfprintf(fail_file,"gnr: doing recursive call, do_line_search, retval, errx = %4i  %4i  %21.15g \n", do_line_search, retval, errx );
      //      
      retval2 = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
      for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
      //      dualfprintf(fail_file,"gnr retval3 = %4i \n", retval2); 
      return( retval2 );
    }
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
#if(DOHISTOGRAM)
    bin_newt_data( errx, n_iter, 1, 0 );
#endif
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   1   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      
    }
    if(ltrace) {
      dualfprintf(fail_file,"general_newton_raphson(): found minimal solution \n");
      
    }
    //    dualfprintf(fail_file,"gnr retval4 = %4i \n", 0); 
#endif
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
#if(DOHISTOGRAM)
    bin_newt_data( errx, n_iter, 2, 0 );
#endif
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   2   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra, errx); 
      
    }
    //    dualfprintf(fail_file,"gnr retval5 = %4i \n", 0); 
#endif
    return(0);
  }

#if(!OPTIMIZED)
  dualfprintf(fail_file,"gnr retval6 = %4i \n", 0);
#endif
  return(0);

}
