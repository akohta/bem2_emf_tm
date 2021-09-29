#if !defined CHANKEL1_01_H
#define CHANKEL1_01_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define EPS_CH1 1.0e-10
#define SMX_CH1	50

// complex argument, ascending series expansions and asymptotic expansions (hankel's expansions)
double complex chankel1_0 (double complex z); // return H_0^(1)(z) ( = J_0(z) + I * Y_0(z) )
double complex chankel1_1 (double complex z); // return H_1^(1)(z) ( = J_1(z) + I * Y_1(z) )
void chankel1_01(double complex z,double complex *H0,double complex *H1); // H_0^(1)(z),

// modified hankel function, argument z must be |z|<5
double complex chankel1_Q_1(double complex z); // return H_1^(1)(z) + 2*i/(M_PI*z)
double complex chankel1_O_1(double complex z); // return -I*M_PI/z*( H_1^(1)(z) + 2*i/(M_PI*z) )= O_1^(1)(z)

// real argument, polynomial approximation (fast)
double besj0(double x);
double besy0(double x);
double besj1(double x);
double besy1(double x);

#endif
