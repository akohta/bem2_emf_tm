/*
 * pb_coef.c
 *
 *  Created on: Sep 13, 2018
 *      Author: ohta
 */
#include "pb_elem.h"

void coefficient_GL(double complex *CC,double *rt,int s,double complex k,DOMD *pc)
{
  void convert_CC(double complex *CC);
  
  double complex H0[3],H1[3],kR;
  double X,Y,R,dx,dy,i_R[3],tD[3],aJ[3];
  int i,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<3;i++){
    X=pc->bd.x[s][i]-rt[0];
    Y=pc->bd.y[s][i]-rt[1];
    R=sqrt(X*X+Y*Y);
    i_R[i]=1.0/R;
    dx=pc->bd.dx[s][i];
    dy=pc->bd.dy[s][i];
    tD[i]=X*i_R[i]*dy-Y*i_R[i]*dx;
    aJ[i]=sqrt(dx*dx+dy*dy);
    kR=k*R;
    chankel1_01(kR,&(H0[i]),&(H1[i]));
  }

  CC[6]=0.0;
  for(i=0;i<3;i++){
    CC[  i]=I/4.0*pc->bd.wg[i]*H0[i]*aJ[i];
    CC[3+i]=-I*k/4.0*pc->bd.wg[i]*H1[i]*tD[i];
    CC[6]+=pc->bd.wg[i]*i_R[i]*tD[i];
  }
  CC[6]*=1.0/(2.0*M_PI);

  if(sig==1) convert_CC(CC);
}

void coefficient_GK(double complex *CGK,double *rt,int s,double complex k,DOMD *pc)
{
  void convert_CC(double complex *CC);
  
  double complex H0[7],H1[7],kR;
  double dr[2],X,Y,R,i_R[7],tD[7],aJ[7];
  int i,j,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<7;i++){
    X=pc->bd.x[s][i]-rt[0]; 
    Y=pc->bd.y[s][i]-rt[1]; 
    R=sqrt(X*X+Y*Y);
    i_R[i]=1.0/R;
    dr[0]=pc->bd.dx[s][i];
    dr[1]=pc->bd.dy[s][i];
    tD[i]=X*i_R[i]*dr[1]-Y*i_R[i]*dr[0]; 
    aJ[i]=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
    kR=k*R;
    chankel1_01(kR,&(H0[i]),&(H1[i]));
  }

  // gauss-kronrod
  for(i=0;i<7;i++)   CGK[i]=0.0;
  for(j=0;j<3;j++){
    for(i=0;i<7;i++){
      CGK[  j]+=pc->bd.wk[i]*H0[i]*pc->bd.M[j][i]*aJ[i];
      CGK[3+j]+=pc->bd.wk[i]*H1[i]*tD[i]*pc->bd.M[j][i];
      if(j==0) CGK[6]+=pc->bd.wk[i]*i_R[i]*tD[i];
    }
  }
  for(j=0;j<3;j++){
    CGK[  j]*=I/4.0;
    CGK[j+3]*=-I*k/4.0;
  }
  CGK[6]*=1.0/(2.0*M_PI);

  if(sig==1) convert_CC(CGK);
}

void coefficient_HP(double complex *CC,double *rt,int s,double complex k,DOMD *pc)
{
  void convert_CC(double complex *CC);
  
  double complex H0,H1,kR;
  double r[2],dr[2],rs[2][3],X,Y,R,i_R,tD,aJ;
  int i,j,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<3;i++){
    rs[0][i]=pc->bd.x[s][i];
    rs[1][i]=pc->bd.y[s][i];
  }

  for(i=0;i<7;i++)   CC[i]=0.0;
  for(i=0;i<GLH;i++){
    rs_eta(r,rs,pc->bd.xh[i]);
    drs_eta(dr,rs,pc->bd.xh[i]);
    X=r[0]-rt[0]; 
    Y=r[1]-rt[1];
    R=sqrt(X*X+Y*Y);
    i_R=1.0/R;
    tD=X*i_R*dr[1]-Y*i_R*dr[0]; 
    aJ=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
    kR=k*R;
    chankel1_01(kR,&H0,&H1);
    
    for(j=0;j<3;j++){
      CC[  j]+=pc->bd.wh[i]*H0*Mn(j,pc->bd.xh[i])*aJ;
      CC[3+j]+=pc->bd.wh[i]*H1*tD*Mn(j,pc->bd.xh[i]);
      if(j==0) CC[6]+=pc->bd.wh[i]*i_R*tD;
    }
  }

  for(j=0;j<3;j++){
    CC[  j]*=I/4.0;
    CC[j+3]*=-I*k/4.0;
  }
  CC[6]*=1.0/(2.0*M_PI);

  if(sig==1) convert_CC(CC); 
}

int coefficient_NV(double complex *CGK,double *rt,int s,double complex k,DOMD *pc)
{
  void convert_CC(double complex *CC);
  
  double complex H0[7],H1[7],CGL[7],kR;
  double dr[2],X,Y,R,i_R[7],tD[7],aJ[7],CD;
  int i,j,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<7;i++){
    X=pc->bd.x[s][i]-rt[0]; if(2.0*fabs(X)/fabs(pc->bd.x[s][i]+rt[0])<MEPS) X=0.0;
    Y=pc->bd.y[s][i]-rt[1]; if(2.0*fabs(Y)/fabs(pc->bd.y[s][i]+rt[1])<MEPS) Y=0.0;
    R=sqrt(X*X+Y*Y);
    i_R[i]=1.0/R;
    dr[0]=pc->bd.dx[s][i];
    dr[1]=pc->bd.dy[s][i];
    tD[i]=X*i_R[i]*dr[1]-Y*i_R[i]*dr[0]; if(2.0*fabs(tD[i])/fabs(X*i_R[i]*dr[1]+Y*i_R[i]*dr[0])<MEPS) tD[i]=0.0;
    aJ[i]=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
    kR=k*R;
    chankel1_01(kR,&(H0[i]),&(H1[i]));
  }

  // gauss-kronrod
  for(i=0;i<7;i++)   CGK[i]=0.0;
  for(j=0;j<3;j++){
    for(i=0;i<7;i++){
      CGK[  j]+=pc->bd.wk[i]*H0[i]*pc->bd.M[j][i]*aJ[i];
      CGK[3+j]+=pc->bd.wk[i]*H1[i]*tD[i]*pc->bd.M[j][i];
      if(j==0) CGK[6]+=pc->bd.wk[i]*i_R[i]*tD[i];
    }
  }
  for(j=0;j<3;j++){
    CGK[  j]*=I/4.0;
    CGK[j+3]*=-I*k/4.0;
  }
  CGK[6]*=1.0/(2.0*M_PI);

  // gauss-legendre
  CGL[6]=0.0;
  for(i=0;i<3;i++){
    CGL[  i]= I/4.0*pc->bd.wg[i]*H0[i]*aJ[i];
    CGL[3+i]=-I*k/4.0*pc->bd.wg[i]*H1[i]*tD[i];
    CGL[6]+=pc->bd.wg[i]*i_R[i]*tD[i];
  }
  CGL[6]*=1.0/(2.0*M_PI);

  // check
  j=0;
  for(i=0;i<7;i++){
    CD=2.0*cabs(CGK[i]-CGL[i])/cabs(CGK[i]+CGL[i]);
    if(CD>IEPS) j++;
  }
  if(j!=0){
    if(sig==0) coefficient_HP(CGK,rt, s,k,pc);
    else       coefficient_HP(CGK,rt,-s,k,pc);
    return 1;
  }
  else{
    if(sig==1) convert_CC(CGK);
    return 0;
  }
}

void coefficient_DE(double complex *CC,double *rt,int s,double complex k,DOMD *pc)
{
  double complex sG(double eta,void *tmp);
  double complex sH(double eta,void *tmp);
  double sF(double eta,void *tmp);
  void convert_CC(double complex *CC);
  
  TDATA td;
  double complex G0,G1,G2,H0,H1,H2,F;
  double err0;
  int i,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  td.k=k;
  td.rt[0]=rt[0];  td.rt[1]=rt[1];
  for(i=0;i<3;i++){
    td.rs[0][i]=pc->bd.x[s][i];
    td.rs[1][i]=pc->bd.y[s][i];
  }

  td.type=0;
  G0=deintz(sG,-1.0,1.0,&td,DEPS,&err0);
  if(err0<0.0){ printf("DE integration error. coefficient_DE(), G0()! Exit...\n");  exit(1);  }
  td.type=1;
  G1=deintz(sG,-1.0,1.0,&td,DEPS,&err0);
  if(err0<0.0){ printf("DE integration error. coefficient_DE(), G1()! Exit...\n");  exit(1);  }
  td.type=2;
  G2=deintz(sG,-1.0,1.0,&td,DEPS,&err0);
  if(err0<0.0){ printf("DE integration error. coefficient_DE(), G2()! Exit...\n");  exit(1);  }

  if( ((2.0*fabs(td.rs[0][0]-rt[0])/fabs(td.rs[0][0]+rt[0])<MEPS) && (2.0*fabs(td.rs[0][1]-rt[0])/fabs(td.rs[0][1]+rt[0])<MEPS)
      && (2.0*fabs(td.rs[0][2]-rt[0])/fabs(td.rs[0][2]+rt[0])<MEPS)) ||
      ((2.0*fabs(td.rs[1][0]-rt[1])/fabs(td.rs[1][0]+rt[1])<MEPS) && (2.0*fabs(td.rs[1][1]-rt[1])/fabs(td.rs[1][1]+rt[1])<MEPS)
          && (2.0*fabs(td.rs[1][2]-rt[1])/fabs(td.rs[1][2]+rt[1])<MEPS))){
    H0=0.0;    H1=0.0;    H2=0.0;
    F=0.0;
  }
  else {
    td.type=0;
    H0=deintz(sH,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. coefficient_DE(), H0()! Exit...\n");  exit(1);  }
    td.type=1;
    H1=deintz(sH,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. coefficient_DE(), H1()! Exit...\n");  exit(1);  }
    td.type=2;
    H2=deintz(sH,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. coefficient_DE(), H2()! Exit...\n");  exit(1);  }

    F=deintd(sF,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. coefficient_DE(), F()! Exit...\n");  exit(1);  }
  }

  CC[0]=I/4.0*G0;
  CC[1]=I/4.0*G1;
  CC[2]=I/4.0*G2;
  CC[3]=-I*td.k/4.0*H0;
  CC[4]=-I*td.k/4.0*H1;
  CC[5]=-I*td.k/4.0*H2;
  CC[6]=1.0/(2.0*M_PI)*F;
  if(sig==1) convert_CC(CC);
}

void coefficient_bd_eta(double complex *CC,double eta_t,int s,double complex k,DOMD *pc)
{
  double complex sG_bd(double xi,void *tmp);
  double complex sH_bd(double xi,void *tmp);
  double sF_bd(double xi,void *tmp);
  void convert_CC(double complex *CC);
  
  TDATA td;
  double err0,err1,F,Cs,ts,tr;
  double complex G0,G1,G2,H0,H1,H2;
  int i,sig;

  if(s>0) sig=0;
  else {
    s=-s;
    eta_t=-eta_t;
    sig=1;
  }

  td.k=k;
  td.eta_t=eta_t;
  for(i=0;i<3;i++){
    td.rs[0][i]=pc->bd.x[s][i];
    td.rs[1][i]=pc->bd.y[s][i];
  }

  td.type=0;
  G0=-deintz(sG_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
        + deintz(sG_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);
  if(err0<0.0||err1<0.0){ printf("DE integration error. coefficient_bd_eta(), G0_bd()! Exit...\n");  exit(1);  }
  td.type=1;
  G1=-deintz(sG_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
        + deintz(sG_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);
  if(err0<0.0||err1<0.0){ printf("DE integration error. coefficient_bd_eta(), G1_bd()! Exit...\n");  exit(1);  }
  td.type=2;
  G2=-deintz(sG_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
        + deintz(sG_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);
  if(err0<0.0||err1<0.0){ printf("DE integration error. coefficient_bd_eta(), G2_bd()! Exit...\n");  exit(1);  }


  ts=(td.rs[0][1]-td.rs[0][2])*td.rs[1][0]+(td.rs[0][2]-td.rs[0][0])*td.rs[1][1]+(td.rs[0][0]-td.rs[0][1])*td.rs[1][2];
  tr=pow(td.rs[0][1]-td.rs[0][0],2)+pow(td.rs[1][1]-td.rs[1][0],2);

  if( fabs(ts)/tr < MEPS){ // Cs=0.0;
    Cs=0.0;    H0=0.0;    H1=0.0;    H2=0.0;    F=0.0;
  }
  else {
    Cs=0.5*i_A_GL*i_A_GL*i_A_GL*ts;
    td.type=0;
    H0= deintz(sH_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
          + deintz(sH_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);
    if(err0<0.0||err1<0.0){ printf("DE integration error. coefficient_bd_eta(), H0_bd()! Exit...\n");    exit(1);  }
    td.type=1;
    H1= deintz(sH_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
          + deintz(sH_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);
    if(err0<0.0||err1<0.0){ printf("DE integration error. coefficient_bd_eta(), H1_bd()! Exit...\n");    exit(1);  }
    td.type=2;
    H2= deintz(sH_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
          + deintz(sH_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);
    if(err0<0.0||err1<0.0){ printf("DE integration error. coefficient_bd_eta(), H2_bd()! Exit...\n");    exit(1);  }

    F=-deintd(sF_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
          +deintd(sF_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);
    if(err0<0.0||err1<0.0){ printf("DE integration error. coefficient_bd_eta(), F_bd()! Exit...\n");    exit(1);  }
  }

  CC[0]= I/4.0*G0;
  CC[1]= I/4.0*G1;
  CC[2]= I/4.0*G2;
  CC[3]=-I*td.k/4.0*Cs*H0;
  CC[4]=-I*td.k/4.0*Cs*H1;
  CC[5]=-I*td.k/4.0*Cs*H2;
  CC[6]= 1.0/(2.0*M_PI)*Cs*F;
  if(sig==1) convert_CC(CC);
}


// ---------------------------------------------------------------------------
void convert_CC(double complex *CC)
{
  double complex tmpg,tmph;
  tmpg=CC[0];
  CC[0]=CC[1];
  CC[1]=tmpg;
  tmph=CC[3];
  CC[3]=-CC[4];
  CC[4]=-tmph;
  CC[5]*=-1.0;
  CC[6]*=-1.0;
}

double complex sG(double eta,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex H0,kR;
  double r[2],dr[2],X,Y,aJ;

  rs_eta(r,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  aJ=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  X=r[0]-td->rt[0];
  Y=r[1]-td->rt[1];
  kR=td->k*sqrt(X*X+Y*Y);
  H0=chankel1_0(kR);

  return H0*Mn(td->type,eta)*aJ;
}

double complex sH(double eta,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex H1,kR;
  double r[2],dr[2],X,Y,R,i_R,tD;

  rs_eta(r,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  X=r[0]-td->rt[0];  if(2.0*fabs(X)/fabs(r[0]+td->rt[0])<MEPS) X=0.0;
  Y=r[1]-td->rt[1];  if(2.0*fabs(Y)/fabs(r[1]+td->rt[1])<MEPS) Y=0.0;
  R=sqrt(X*X+Y*Y);
  i_R=1.0/R;
  kR=td->k*R;
  H1=chankel1_1(kR);
  tD=X*i_R*dr[1]-Y*i_R*dr[0];  if(2.0*fabs(tD)/fabs(X*i_R*dr[1]+Y*i_R*dr[0])<MEPS) tD=0.0;

  return H1*tD*Mn(td->type,eta);
}

double sF(double eta,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double r[2],dr[2],X,Y,i_R,tD;

  rs_eta(r,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  X=r[0]-td->rt[0];  if(2.0*fabs(X)/fabs(r[0]+td->rt[0])<MEPS) X=0.0;
  Y=r[1]-td->rt[1];  if(2.0*fabs(Y)/fabs(r[1]+td->rt[1])<MEPS) Y=0.0;
  i_R=1.0/sqrt(X*X+Y*Y);
  tD=X*i_R*dr[1]-Y*i_R*dr[0];  if(2.0*fabs(tD)/fabs(X*i_R*dr[1]+Y*i_R*dr[0])<MEPS) tD=0.0;
  return i_R*tD;
}

double complex sG_bd(double xi,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex H0,kR;
  double dr[2],aJ;

  drs_eta(dr,td->rs,xi+td->eta_t);
  aJ=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);

  drs_eta(dr,td->rs,0.5*xi+td->eta_t);
  kR=td->k*fabs(xi)*sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  H0=chankel1_0(kR);

  return H0*Mn(td->type,xi+td->eta_t)*aJ;
}

double complex sH_bd(double xi,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex H1,kR;
  double dr[2],aJR;

  drs_eta(dr,td->rs,0.5*xi+td->eta_t);
  aJR=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  kR=td->k*fabs(xi)*aJR;
  H1=chankel1_1(kR);

  return H1*Mn(td->type,xi+td->eta_t)*xi/aJR;
}

double sF_bd(double xi,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double dr[2];

  drs_eta(dr,td->rs,0.5*xi+td->eta_t);
  return 1.0/(dr[0]*dr[0]+dr[1]*dr[1]);
}
