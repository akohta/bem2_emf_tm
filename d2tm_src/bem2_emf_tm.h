/*
 * d2_md_h2.h
 *
 *  Created on: Sep 17, 2018
 *      Author: ohta
 */

#ifndef BEM2_EMF_TM_H_
#define BEM2_EMF_TM_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <mkl.h>

#include "physical_constant.h"
#include "osu_mksa.h"
#include "my_utils.h"
#include "gauleg.h"
#include "incident_field.h"
#include "pb_elem.h"


// -- d2tm_setup.c --
void read_data(int argc,char **argv,DOMD *md);    // read datafile, for bv_solver
void print_data(DOMD *md);                        // print data
void print_data_MKSA(DOMD *md);                   // print data in MKSA system of units
void initialize_domd(DOMD *md);                   // initialize the data, for bv_solver 
void mfree_domd(DOMD *md);                        // free allocated memory
int domain_id(double *rt,DOMD *md);               // return domain id at point rt, rt=(x,y)
void dat_read (char *dname,DOMD *md);             // read datafile output by dat_write()
void dat_write(char *dname,DOMD *md);             // write datafile
void output_node_particles(char *fname,DOMD *md); // outputs the nodes as point cloud data ( .particles file ) 


// -- d2tm_solve_bieq.c --
void solve_bieq(DOMD *md);                        // solve boundary integral equations, for bv_solver


// -- d2tm_field.c --
int Hz_s(double complex *Hz,double *rt,int type,DOMD *md); // scattered or internal field
int Hz_t(double complex *Hz,double *rt,int type,DOMD *md); // total field 
int Hz_i(double complex *Hz,double *rt,int type,DOMD *md); // incident field
int Er_s(double complex *Er,double *rt,int type,DOMD *md); // scattered or internal field
int Er_t(double complex *Er,double *rt,int type,DOMD *md); // total field
int Er_i(double complex *Er,double *rt,int type,DOMD *md); // incident field
int HE_s(double complex *Hz,double complex *Er,double *rt,int type,DOMD *md); // scattered or internal field
int HE_t(double complex *Hz,double complex *Er,double *rt,int type,DOMD *md); // total field
int HE_i(double complex *Hz,double complex *Er,double *rt,int type,DOMD *md); // incident field
int HE_s_dbieq(double complex *Hz,double complex *Er,double *rt,DOMD *md); // for far-field
int HE_t_dbieq(double complex *Hz,double complex *Er,double *rt,DOMD *md); // 
int HE_i_dbieq(double complex *Hz,double complex *Er,double *rt,DOMD *md); //
// optputs
// Hz=H_z, Er[0]=E_x,Er[1]=E_y, return domain id. 
// inputs 
// rt[0]=x, rt[1]=y, type : select integration method, md : pointer of DOMD object.
// type=0 : 3-point Gauss-Legendre, type=1 : 7 point Guass-Kronrod (Gauss-Kronrod extension of 3-point GL)
// type=2 : GLH (defined in d2te_const.h) point Gauss-Legendre, type=3 : DE method (for test)
// type>3 : integration with numerical validation 

// boundary value 
double complex Hz_bv(int did,double eta_t,int t,DOMD *md);
void HE_bv(double complex *Hz,double complex *Er,int did,double eta_t,int t,DOMD *md);
// inputs
// did : domain id, eta_t : parameter of point on each element ( -1 ~ 1 ), t : element id, pointer of DOMD object.

// tangential directinal derivative of boundary vale at the node
double complex dHzdt_bv_node_ndmtd(int did,int t,int tn,DOMD *md); // using numerical differentiation
double complex dHzdt_bv_node_dbieq(int did,int t,int tn,DOMD *md); // using derivative boundary integral equation
// return tangential directinal derivative of Hz at the node.
// inputs
// did : domain id, t : element id, tn : node id (0~2), md : pointer of DOMD object.


// -- d2tm_force.c --
void force_FN(double *Fr,double *Nz,double *rc,int type,DOMD *md);
// outputs
// Fr[0]=F_x, Fr[1]=F_y (radiation force), Nx=N_z (radiation torque)
// inputs
// rc[0]=x,rc[1]=y : coordinate of rotation center,
// type : select integration medthod, md : pointer of DOMD object.
// type=0 : 3-point Gauss-Legendre, type>0 : 7-point Gauss-Kronrod.

#endif
