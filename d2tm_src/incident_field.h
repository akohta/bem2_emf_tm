/*
 * incident_field.h
 *
 *  Created on: Sep 17, 2018
 *      Author: ohta
 */

#ifndef INCIDENT_FIELD_H_
#define INCIDENT_FIELD_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "osu_mksa.h"


typedef struct incident_field_data{
  char fname[128];   // datafile name
  double complex E0; // amplitude
  double lambda0;    // wavelength in vacuum
  double ne;         // refractive index
  double angle;      // incident angle

  double k0;         // wave number in vacuum (=angular frequency in this system)
  double kxh,kyh;    // wave number unit vector
}INFD;

// -- incident_field.c (TE mode plane wave) --
void read_infd(char *fname,INFD *wd); // read incident field data
void print_infd(INFD *wd);            // print incident field data
void print_infd_MKSA(INFD *wd);       // print incident field data in MKSA system of units

void infd_grad(double complex *Hz,double complex *gradHz,double *r,INFD *wd);
// outputs
// Hz=H_z, gradHz[0]=dH_z/dx, gradHz[1]=dH_z/dy
// inputs 
// r[0]=x, r[1]=y, pointer of INFD 
double complex infd_Hz(double *r,INFD *wd); 
// return H_z
void infd_Er(double complex *Er,double *r,INFD *wd); 
// outputs
// Er[0]=E_x, H[1]=E_y
void infd_HE(double complex *Hz,double complex *Er,double *r,INFD *wd);
// outputs
// Hz=H_z, Er[0]=E_x, Er[1]=E_y


#endif /* INCIDENT_FIELD_H_ */
