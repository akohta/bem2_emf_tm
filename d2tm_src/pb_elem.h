/*
 * pb_elem.h
 *
 *  Created on: Sep 12, 2018
 *      Author: ohta
 */

#ifndef PB_ELEM_H_
#define PB_ELEM_H_

#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "incident_field.h"
#include "pb_const.h"
#include "de_int.h"
#include "chankel1_01.h"

typedef struct sub_domain_data{
  int Ne;   // number of elements
  int *eid; // element id

  double complex **u,**dudn; // main field boundary value
  double complex **v,**dvdn; // sub field boundary value
  double complex **w,**dwdn; // sub field boundary value
}SUBD;

typedef struct boundary_data{
  int Nn;     // number of nodes
  int Ne;     // number of elements
  // element constant
  double et[7];       // nodes 
  double wg[7],wk[7]; // wg : weight of gauss node, wk : weight of kronrod node
  double M[3][7];     // shape function value of each nodes
  double *xh,*wh;     // gauss-legender nodes 
  // node data
  double *xn,*yn; // xn[node id]=x-coordinate of node, yn[node id]=y-coordinate of node
  // element data
  int **ed;   // ed[element id][position id(0,1,2)]=node id
  int *md;    // main domain id
  int *sd;    // sub domain id
  int *egd;   // group id (elementary geometrical entity)
  double **x,**y;   // gauss node coordinate, r[element id][position id(0,1,2)]
  double **dx,**dy; // dx/deta and dy/deta value of the gauss nodes
  double complex **ui,**duidn; // incidnt field data
  // sub domain data
  SUBD *sb;
}BOUD;

typedef struct domain_data{
  char med_fn[128],msh_fn[128]; // medium_datafile_name , mesh data file name
  double ra,tx,ty; // rotation angle, translation vector (tx,ty)

  int MN; // number of mediums
  double complex *n; // complex refractive index of medium

  INFD wd; // incident field data
  BOUD bd; // boundary data
}DOMD;

typedef struct tmp_data{
  double eta_t,pad0;
  double complex k;
  double rs[2][3],c1,c2,cr[2];
  double rt[2],vt[2],drt[2];
  int type,node,pad[2];
}TDATA;

// -- pb_elem.c --
double Mn(int n,double eta); // calc shape function Mn(eta), n=0,1,2
double Nn(int n,double eta); // calc alpha=1 shape function Nn(eta), n=0,1,2
void  rn_eta(double  *r,double rn[][3],double eta);  // calc boundary coordinate from the node (eta=-1,1,0)
void drn_eta(double *dr,double rn[][3],double eta);  // dx/deta and dy/deta value of the coordinate eta
void  rs_eta(double  *r,double rs[][3],double eta);  // calc boundary coordinate from the gauss node (eta=-a,a,0)
void drs_eta(double *dr,double rs[][3],double eta);  // dx/deta and dy/deta value of the coordinate eta (rs is gauss node)

// -- pb_coef.c --
void coefficient_GL(double complex *CC,double *rt,int s,double complex k,DOMD *md); // 3 point GL
void coefficient_GK(double complex *CC,double *rt,int s,double complex k,DOMD *md); // 7 point GK
void coefficient_HP(double complex *CC,double *rt,int s,double complex k,DOMD *md); // GLH(defined in d2te_const.h) point GL
void coefficient_DE(double complex *CC,double *rt,int s,double complex k,DOMD *md); // Double-exponential method (test function)
int  coefficient_NV(double complex *CC,double *rt,int s,double complex k,DOMD *md); // with numerical validation using Gauss-Kronrod extension
void coefficient_bd_eta(double complex *CC,double eta_t,int s,double complex k,DOMD *md);

// -- pb_dcoef.c --
void d_coef_grad_GL(double complex *CC,double complex *Cx,double complex *Cy,double *rt,int s,double complex k,DOMD *md);
void d_coef_GL(double complex *Cv,double *vt,double *rt,int s,double complex k,DOMD *md);
void d_coef_GK(double complex *Cv,double *vt,double *rt,int s,double complex k,DOMD *md);
void d_coef_HP(double complex *Cv,double *vt,double *rt,int s,double complex k,DOMD *md);
void d_coef_DE(double complex *Cv,double *vt,double *rt,int s,double complex k,DOMD *md);
int  d_coef_NV(double complex *Cv,double *vt,double *rt,int s,double complex k,DOMD *md);
void d_coef_bd_node_t(double complex *Cv,int s,int sn,double complex k,DOMD *md);


#endif 
