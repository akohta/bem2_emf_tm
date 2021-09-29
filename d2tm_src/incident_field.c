/*
 * incident_field.c
 *
 *  Created on: Sep 17, 2018
 *      Author: ohta
 */

#include "incident_field.h"


void read_infd(char *fname,INFD *wd)
{
  FILE *fp;
  if((fp=fopen(fname,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",fname);    exit(1);  }
  strcpy(wd->fname,fname);
  char buf[256]="";   double tmpd,tmpd2;
  fgets(buf,256,fp);  fgets(buf,256,fp);

  fscanf(fp,"%lf",&tmpd);
  fscanf(fp,"%lf",&tmpd2);   wd->E0      =tmpd+I*tmpd2;
  fscanf(fp,"%lf",&tmpd);    wd->lambda0 =tmpd;
  fscanf(fp,"%lf",&tmpd);    wd->ne      =tmpd;
  fscanf(fp,"%lf\n",&tmpd);  wd->angle   =tmpd;
  fclose(fp);

  // init parameter
  wd->k0=2.0*M_PI/wd->lambda0;
  wd->kxh=cos(wd->angle);
  wd->kyh=sin(wd->angle);
}

void print_infd(INFD *wd)
{
  printf("-- incidnt field data --\n");
  printf("filename                 : %s\n",wd->fname);
  printf("planewave amplitude E0   :%7.6g+%7.6gI\n",creal(wd->E0),cimag(wd->E0));
  printf("wavelength in vacuum     : %15.14g\n",wd->lambda0);
  printf("refractive index         : %15.14g\n",wd->ne);
  printf("incident angle      [rad]: %15.14g\n",wd->angle);
}

void print_infd_MKSA(INFD *wd)
{
  printf("-- incidnt field data --\n");
  printf("filename                   : %s\n",wd->fname);
  printf("planewave amplitude E0[N/C]:%7.6g+%7.6gI\n",creal(OSUtoMKSA_ElectricField(wd->E0)),cimag(OSUtoMKSA_ElectricField(wd->E0)));
  printf("wavelength in vacuum    [m]: %15.14g\n",OSUtoMKSA_length(wd->lambda0));
  printf("refractive index           : %15.14g\n",wd->ne);
  printf("incident angle        [rad]: %15.14g\n",wd->angle);
}

void infd_grad(double complex *Hz,double complex *gradHz,double *r,INFD *wd)
{
  double theta=wd->k0*wd->ne*(wd->kxh*r[0]+wd->kyh*r[1]);
  *Hz=wd->ne*wd->E0*(cos(theta)+I*sin(theta));
  gradHz[0]=I*wd->k0*wd->ne*wd->kxh*(*Hz);
  gradHz[1]=I*wd->k0*wd->ne*wd->kyh*(*Hz);
}

double complex infd_Hz(double *r,INFD *wd)
{
  double theta=wd->k0*wd->ne*(wd->kxh*r[0]+wd->kyh*r[1]);
  return wd->ne*wd->E0*(cos(theta)+I*sin(theta));
}

void infd_Er(double complex *Er,double *r,INFD *wd)
{
  double theta=wd->k0*wd->ne*(wd->kxh*r[0]+wd->kyh*r[1]);
  double _Complex ec=wd->E0*(cos(theta)+I*sin(theta));
  Er[0]=-wd->kyh*ec;
  Er[1]= wd->kxh*ec;
}

void infd_HE(double complex *Hz,double complex *Er,double *r,INFD *wd)
{
  double theta=wd->k0*wd->ne*(wd->kxh*r[0]+wd->kyh*r[1]);
  double _Complex ec=wd->E0*(cos(theta)+I*sin(theta));
  *Hz=wd->ne*ec;
  Er[0]=-wd->kyh*ec;
  Er[1]= wd->kxh*ec;
}
