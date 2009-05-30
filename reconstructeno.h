#ifndef RECONSTRUCTENO_H
#define RECONSTRUCTENO_H

#define DO_ENO_STENCIL_REDUCTION 1  //atch

#define DO_NORMALIZE_SMOOTHNESS_INDICATORS 1 //atch normbeta
#define DO_RESCALE_WEIGHTS 0

#define USE_NORMU 1 // use WENO epsilon that is proportional to the value of the function

// c2e stuff
#define DO_LIMIT_C2E_CORRECTION 0
#define WENO_C2E_MAX_FRAC_DIFFERENCE (0.5)
#define DO_CENT2EDGE 1 // whether to do c2e; if set to 0, edge values are equal to the center value
#define WENO_C2E_REDUCE_NEAR_CUSPS 0 // let mono handle this


// a2c/c2a stuff
#define DO_AVG2CEN 1 // whether to do a2c
#define DO_CEN2AVG 1 // whether to do c2a
#define WENO_AC_REDUCE_NEAR_CUSPS 0 // let mono handle this
#define DO_LIMIT_AC_CORRECTION 1
#define WENO_AC_MAX_FRAC_DIFFERENCE (0.5)  //maximum fractional change allowed during the a2c/c2a conversion; should be >= 0 and < 1
#define WENO_DO_FLIP_CONS_SIGN 1
#define WENO_DIR_FLIP_CONS_SIGN_DN ( (MCOORD == KSCOORDS&&WENO_DO_FLIP_CONS_SIGN)?(2):( (MCOORD == CYLMINKMETRIC&&WENO_DO_FLIP_CONS_SIGN)?(1):(0) ) )
#define WENO_DIR_FLIP_CONS_SIGN_UP ( (MCOORD == KSCOORDS&&WENO_DO_FLIP_CONS_SIGN)?(2):( (MCOORD == CYLMINKMETRIC&&WENO_DO_FLIP_CONS_SIGN)?(0):(0) ) )


#define SYMMETRIZE_SIMPLE_WENO 1
#define SYMMETRIZE_WENO 1
#define SSYMCHECKS  0 //atch symmetry checks flag

#define DO_OR_STENCIL_REDUCTION (0)  //if equals 1, uses the "OR" version of stencil reduction, if 0 uses the "AND" (formerly l&r version), if (-1) -- the combined (AND+OR)/2 version
#define DO_OR_STENCIL_REDUCTION_REDUCE_AC_EXTRA_POINT (0)  //helps to keep the hard test from RAM paper much more stable by reducing the a2c more, only appies when DO_OR_STENCIL_REDUCTION == 1
#define WENO_REDUCE_COMBINED_FACTOR (0.1)  //only applies for DO_OR_STENCIL_REDUCTION = -1

////////
//
//  possible settings for do_weight_or_recon
//
#define WEIGHT_CALC 1
#define RECON_CALC 2
#define ALL_CALC (WEIGHT_CALC | RECON_CALC)  //compute everything

//0 -- don't do weighs minimization
// -- do 1st version of weights minimization
#define NOSPLITA2C 0
#define MINIMIZE_ALL_WEIGHTS 1
#define ENERGY_CONTROLS_ALL_WEIGHTS 2

#define DO_SPLITA2C ENERGY_CONTROLS_ALL_WEIGHTS
//#define DO_SPLITA2C NOSPLITA2C






////////
//
//  MONOINTERP 
//
//#define DOMONOINTERP NOMONOINTERP
//#define DOMONOINTERP JMONOINTERP
#define DOMONOINTERP SMONOINTERP
#define NOMONOINTERP 0
#define JMONOINTERP 1
#define SMONOINTERP 2

//whether to smoothly rescale the SWENO weights in accordance with smoothness indicators to allow for smooth change in the behaviour of SWENO stencil reduction
#define DO_SMOOTH_ADJUSTMENT_UNOPTIMIZED_WEIGHTS 1
#define SMONO_DO_SMOOTH_TRANSITION 1
#define SMONO_EPSILON SQRT_WENO_EPSILON


//SMONO settings (if SMONO is enabled; these should not be here but it is nice to have them all in one file)
#define DO_SMONO_C2A 0
#define DO_SMONO_A2C 0
#define DO_SMONO_C2E 1

#define DO_MONO_1ST_DERIVATIVE 1
#define DO_3RD_DER 1
#define DO_SMONO_C2E_CUSP_INDICATOR 1
#define DO_SMONO_A2C_CUSP_INDICATOR 1
#define DO_SMONO_C2A_CUSP_INDICATOR 1

////////
//
//  POST-/PRESHOCK REDUCTION
//
// reduce extra point downstream of shock
#define DO_REDUCE_POST_SHOCK 1
// reduce extra point upstream of shock
#define DO_REDUCE_PRE_SHOCK 1
// Fractional pressure jump for which the POST/PRESHOCK reduction fully takes place; it is linearly interpolated to smaller pressure jumps
#define REDUCE_DP_FRACTION (0.5)



////////
//
//  EXPERIMENTAL/DIDN'T WORK
//
#define DO_LIMIT_AC_CORRECTION_NEAR_INFLECTIONS 0
#define DO_MONOTONICITY_INDICATOR 0 // using new mono in interpline.c
#define DESENSITISE_STENCIL_REDUCTION 0  //whether to add epsilon in the calculation of the weights used for stencil reduction
#define WENO_HIGH_ORDER_CENTRAL_WEIGHT 0
#define DO_LIMIT_C2E_WEIGHTS 0  //forces weights to be exactly equal if they are sufficiently close; forces those that are small but nonzero to be exactly zero
#define DO_LIMIT_AC_WEIGHTS 0  //forces weights to be exactly equal if they are sufficiently close; forces those that are small but nonzero to be exactly zero
////////


#define FORCE_AC_CA_WEIGHT_TO_BE_OPTIMAL 0
#define FORCE_C2E_WEIGHT_TO_BE_OPTIMAL 0


// maximum order of scheme indicators can handle
#define MAXORDERS MAXSPACEORDER


#define MIN_CVT_ORDER 1
#define MAX_CVT_ORDER 7
#define CVT_A2C 0
#define CVT_C2A 1
#define CVT_C2L 2
#define CVT_C2R 3
#define NUM_CVT_TYPES 4

#define CVT_C2E CVT_C2L //use the same number as CVT_C2L because does not add a new reconstruction type

#define DO_ASSERTS 0

#if(!DO_ASSERTS)
#define assert {}
#else
#define assert assert_func
#endif

FTYPE weno_point_updates[MAX_CVT_ORDER];    //atch symmetrize
FTYPE weno_stencil_updates[MAX_CVT_ORDER];  //atch symmetrize

#define MAX_STENCIL_LENGTH 5  //added by atch (to use as the max. size of an array)
#define MAX_NO_OF_WEIGHTS (3*MAX_STENCIL_LENGTH)  //we maximally store three sets of weights optimized for different points, that's why need to mult. by 3

struct weno_weights_s {
	double weights[MAX_NO_OF_WEIGHTS];		//array that holds the weights
	int len;                            //total number of weights; it is larger than order because both non-optimal and optimal weights are stored in sequence
	int order;													//reconstruction order for which the weights are computed
#if( WENO_C2E_REDUCE_NEAR_CUSPS || WENO_AC_REDUCE_NEAR_CUSPS  )
	int do_reduce;
#endif

#if( WENO_HIGH_ORDER_CENTRAL_WEIGHT )
	FTYPE high_order_central_weight;
#endif

	double w_ratio_min, w_ratio_max, w_ratio_left, w_ratio_right, w_ratio_combined, lower_order_fraction, lower_order_fraction_fordpp;
};

typedef struct weno_weights_s weno_weights_t;

#define LOWER_ORDER_PREFACTOR (1.) //set to zero to switch off stencil reduction
                                         //fewer the number of times when the stencil size is reduced
#define WENO_MIN_NONZERO_WEIGHT (1.) //minimum weight that is considered non-zero

//The interval of weights where a linear combination of reconstructions with different orders are used
//Use different limiting weught ratios for OR and AND versions of stencil reduction.  Because the OR version is more
//sensitive to these ratios, they should be larger for it
#if( DO_OR_STENCIL_REDUCTION )
	#define MAX_TRANSITION_RATIO   (15.)
	#define MIN_TRANSITION_RATIO    (9.)
	#define MAX_TRANSITION_RATIO_AC (15.) //(1.6)  using such small weight thresholds does not work for 2d Noh problem -- lots of negative internal energies in smooth preshock region
	#define MIN_TRANSITION_RATIO_AC (9.)  //(1.3)
#else
	#define MAX_TRANSITION_RATIO (1.6)
	#define MIN_TRANSITION_RATIO (1.3)
	#define MAX_TRANSITION_RATIO_AC (15.)  //using such small weight thresholds does not work for 2d Noh problem -- lots of negative internal energies in smooth preshock region
	#define MIN_TRANSITION_RATIO_AC (9.)
#endif

#define WENO_DELTA_I (1)

#define MONOYIN 0
#define LEFTYOUT 1
#define RIGHTYOUT 2
#define CENTYOUT 1


//Implemented types of cell-centers to cel avg. reconstruction:
#define USEENO 0
#define USEWENO 1

#define WHICHENO USEWENO  //which type of cell-centers to cel avg. reconstruction to use
//#define WENO_EPSILON (1.e-37)
#define WENO_EPSILON (1.e-26)
#define SQRT_WENO_EPSILON (1.e-13)  //used for determining the scale at which the weights change significantly

//#define REDUCEEPSILON (NUMEPSILON * 100.)
//#define REDUCEEPSILON (1E-10)
#define REDUCEEPSILON WENO_EPSILON
#define DESENISTISE_STENCIL_REDUCTION_EPSILON (1E-10)

#define MATRIX_ROW_COL( matrix_name, index_row, index_col, num_columns ) matrix_name[(index_row) * (num_columns) + (index_col)]

//FTYPE a2c( FTYPE uc[][N2M][NPR], FTYPE ua[][N2M][NPR] );
//FTYPE c2a( FTYPE ut[][N2M][NPR], FTYPE uc[][N2M][NPR] );

extern FTYPE get_sum_of_elements( int number_of_elements, FTYPE *array_to_sum_up );

extern int eno_line_a2c( int whichquantity, int do_weight_or_recon, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, FTYPE *shockindicator,
												 FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE *P, FTYPE *yin,  FTYPE *yout);
extern int eno_line_c2a( int whichquantity, int do_weight_or_recon, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, FTYPE *shockindicator, 
												 FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE *P, FTYPE *yin,  FTYPE *yout);
extern int eno_line_c2e( int whichquantity, int do_weight_or_recon, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, FTYPE *shockindicator, 
												 FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE *P, FTYPE *yin,  FTYPE *yout_left, FTYPE *yout_right);
extern void c2e_simple_eno(int full_order, int is_interpolate_to_right, FTYPE *yin, FTYPE *pout);
extern void a2c_simple_eno(int order, FTYPE *yin, FTYPE *pout);
extern void c2a_simple_eno(int order, FTYPE *yin, FTYPE *pout);
extern void c2e_simple_weno(int order, int ii, int bs, int be, FTYPE *yin, FTYPE *pleft, FTYPE *pright);

extern FTYPE limit_ac_correction( int order, int pl, int bs, int bf, FTYPE max_frac_difference, FTYPE *yin, FTYPE *yout );



extern int assert_func( int is_bad_val, char *s, ... );





#endif
