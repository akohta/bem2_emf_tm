/*
 * pb_elm.c
 *
 *  Created on: Sep 12, 2018
 *      Author: ohta
 */
#include "pb_elem.h"

typedef double (*FNAME)(double);

double N0(double eta);
double N1(double eta);
double N2(double eta);
FNAME FN[3]={N0,N1,N2}; // function pointer

double M0(double eta);
double M1(double eta);
double M2(double eta);
FNAME FM[3]={M0,M1,M2}; // function pointer


double Mn(int n,double eta)
{
  return (*FM[n])(eta);
}

double Nn(int n,double eta)
{
  return (*FN[n])(eta);
}

void rn_eta(double *r,double rn[][3],double eta)
{
  r[0]=(0.5*(rn[0][0]+rn[0][1])-rn[0][2])*eta*eta+0.5*(rn[0][1]-rn[0][0])*eta+rn[0][2];
  r[1]=(0.5*(rn[1][0]+rn[1][1])-rn[1][2])*eta*eta+0.5*(rn[1][1]-rn[1][0])*eta+rn[1][2];
}

void drn_eta(double *dr,double rn[][3],double eta)
{
  dr[0]=(rn[0][0]+rn[0][1]-2.0*rn[0][2])*eta+0.5*(rn[0][1]-rn[0][0]);
  dr[1]=(rn[1][0]+rn[1][1]-2.0*rn[1][2])*eta+0.5*(rn[1][1]-rn[1][0]);
}

void rs_eta(double *r,double rs[][3],double eta)
{
  double epa;
  epa=eta*i_A_GL;
  r[0]=(0.5*(rs[0][0]+rs[0][1])-rs[0][2])*epa*epa+0.5*(rs[0][1]-rs[0][0])*epa+rs[0][2];
  r[1]=(0.5*(rs[1][0]+rs[1][1])-rs[1][2])*epa*epa+0.5*(rs[1][1]-rs[1][0])*epa+rs[1][2];
}

void drs_eta(double *dr,double rs[][3],double eta)
{
  double epa;
  epa=eta*i_A_GL;
  dr[0]=i_A_GL*((rs[0][0]+rs[0][1]-2.0*rs[0][2])*epa+0.5*(rs[0][1]-rs[0][0]));
  dr[1]=i_A_GL*((rs[1][0]+rs[1][1]-2.0*rs[1][2])*epa+0.5*(rs[1][1]-rs[1][0]));
}

double N0(double eta)
{
  return -0.5*(1.0-eta)*eta;
}

double N1(double eta)
{
  return  0.5*(1.0+eta)*eta;
}

double N2(double eta)
{
  return (1.0-eta)*(1.0+eta);
}

double M0(double eta)
{
  return -0.5*(1.0-eta/A_GL)*eta/A_GL;
}

double M1(double eta)
{
  return 0.5*(1.0+eta/A_GL)*eta/A_GL;
}

double M2(double eta)
{
  return (1.0-eta/A_GL)*(1.0+eta/A_GL);
}

