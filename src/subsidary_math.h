/**
    LDPC
    subsidary_math.h
    Purpose: Ñontains subsidary mathematical operations 
	for using by other methods.

    @author Leonov V.R.
    @version 24.04.18
*/


#ifndef SUBSIDARY_MATH_H_INCLUDED
#define SUBSIDARY_MATH_H_INCLUDED


/**
    Fills int array with random permutation.

    @param x Int array.
    @param n Array length.
*/
void fill_with_permutation(int *x, int n);

/**
    Calculate log((exp(x) - 1) / (exp(x) + 1)).

    @param x Some float value.
    @return log((exp(x) - 1) / (exp(x) + 1)).
*/
float log_exp(float x);

/**
    Calculate log((exp(x) + 1) / (exp(x) - 1)).

    @param value Some float value.
    @return log((exp(x) + 1) / (exp(x) - 1)).
*/
float log_tahn(float value);

/**
    Determines the sign of a value.

    @param value Some float value.
    @return 1, if value is positive,
    -1, if value is negative,
    0, if value is zero.
*/
int sign(float value);

/**
    Generates a random number from the normal distribution.

    @return Random number from the normal distribution.
*/
float randn();

#endif // SUBSIDARY_MATH_H_INCLUDED
