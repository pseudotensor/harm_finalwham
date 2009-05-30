
#define PRIMECOORDS -1 // whatever the prime coordinate/metric is, used in transforms.c
#define CARTMINKMETRIC 0 // cartesian that is
#define BLCOORDS 1 // just set Rin oustide horizon for normal star
#define KSCOORDS 2 // as with BLCOORDS
#define HTMETRIC 3 // this defines exterior only, need specific Rin[\theta]
#define CYLMINKMETRIC 4 // cylindrical minkowski
#define HTMETRICACCURATE 5 // consistent expansion form of HT metric
#define SPCMINKMETRIC 6
#define UNIGRAVITY 7 // for RT test problem

// 0 : Boyer-Lindquist (based on r theta)
// 1 : Kerr-Schild (based on bl coords)
// contains metric definitions
// Kerr-Schild: future regularized, ep=-1, k=1

