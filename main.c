// #include "decs.h"
#include "defs.h"


int main(int argc, char *argv[])
{
  int size;
  extern int step_ch(void);
  extern int init(int argc, char *argv[]);
  //  extern void testffdeinversion(void);

  /* perform initializations */
  if (init(argc,argv) >= 1) {
    dualfprintf(fail_file, "main:init: failure\n");
    myexit(1);
  }
	  

  /* Do initial diagnostics */
  if (DODIAGS) {
    trifprintf("proc: %04d : Start initial diagnostics\n", myid);
    // no error_check since if init passed, diag(0) should pass
    diag(0);
    trifprintf("proc: %04d : End initial diagnostics\n", myid);
  }

  //  testffdeinversion();
  //  myexit(0);


  trifprintf("proc: %04d : Start computation\n", myid);

  while (t < tf) {


    /* step variables forward in time */
    nstroke = 0;
    
    step_ch();

    // must check before MPI operation (since asymmetries would
    // desynchronize cpus)
    MYFUN(error_check(2),"main.c","error_check",1);
    //^^ otherwise ok
    
    // eventually all cpus come here, either in failure mode or not,
    // and cleanly tell others if failed/exit/dump/etc.

    postdt(); // here one can alter variables and try to restart, or implement any post step operations

    /* restart dump */
    // if(nstep == 130) restart_write(1) ;
    //    if(nstep==5) break;

    nstep++;
    // restartsteps[whichrestart]=realnstep;


    /* perform diagnostics */  
    // no error check since assume if step_ch passed, diag(1) will pass

    if (DODIAGS && !DODIAGEVERYSUBSTEP){
      pdump = p;
      diag(1);
#if(PRODUCTION==0)
      trifprintf( "D");
#endif
    }

#if(PRODUCTION)
    if((!(nstep%100))||(nstep<100))
#else
      if(1)
#endif
	{
	  fprintf(log_file,"step: %21.15g %21.15g %21.15g %8ld %8ld %21.15g\n", t, dt,
		  cour, nstep, realnstep, 1.0*nstroke/(2.0*N1*N2)); fflush(log_file);
	  mpiisum0(&nstroke,0);
	  myfprintf(logfull_file,"step: %21.15g %21.15g %21.15g %8ld %8ld %21.15g\n", t, dt,
		    cour, nstep, realnstep, 1.0*nstroke/(2.0*totalzones));
	  myfprintf(stderr,"step: %21.15g %21.15g %21.15g %8ld %8ld %21.15g\n", t, dt,
		    cour, nstep, realnstep, 1.0*nstroke/(2.0*totalzones));
	}
  }
  trifprintf("proc: %04d : End computation\n", myid);

  /* do final diagnostics */
  if (DODIAGS){
    pdump = p;
    diag(2);
  }


  trifprintf("ns,ts: %ld %lld\n", nstep, (long long)(nstep) *(long long)totalsize[1] * (long long)totalsize[2]);

  myexit(0);
  return (0);
}
