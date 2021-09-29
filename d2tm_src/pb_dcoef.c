/*
 * pb_dcoef.c
 *
 *  Created on: Sep 13, 2018
 *      Author: ohta
 */
#include "pb_elem.h"


void d_coef_grad_GL(double complex *CC,double complex *Cx,double complex *Cy,double *rt,int s,double complex k,DOMD *md)
{
  void convert_CC(double complex *CC); // coef.c
  void convert_Cv(double complex *Cv);
  
  double complex H0[3],H1[3],kR;
  double X,Y,R,dx,dy,i_R[3],tD[3],vDx[3],vDy[3],aJ[3],nvx[3],nvy[3];
  int i,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<3;i++){
    X=md->bd.x[s][i]-rt[0];
    Y=md->bd.y[s][i]-rt[1];
    R=sqrt(X*X+Y*Y);
    i_R[i]=1.0/R;
    dx=md->bd.dx[s][i];
    dy=md->bd.dy[s][i];
    tD[i]=X*i_R[i]*dy-Y*i_R[i]*dx;
    vDx[i]=X*i_R[i];
    vDy[i]=Y*i_R[i];
    aJ[i]=sqrt(dx*dx+dy*dy);
    nvx[i]=dy;
    nvy[i]=-dx;
    kR=k*R;
    chankel1_01(kR,&(H0[i]),&(H1[i]));
  }

  for(i=0;i<3;i++){
    CC[  i]= I/4.0*md->bd.wg[i]*H0[i]*aJ[i];
    Cx[  i]= I*k/4.0*md->bd.wg[i]*H1[i]*vDx[i]*aJ[i];
    Cy[  i]= I*k/4.0*md->bd.wg[i]*H1[i]*vDy[i]*aJ[i];
    CC[3+i]=-I*k/4.0*md->bd.wg[i]*H1[i]*tD[i];
    Cx[3+i]=I*k/4.0*md->bd.wg[i]*(k*H0[i]*tD[i]*vDx[i]-2.0*i_R[i]*H1[i]*(tD[i]*vDx[i]-0.5*nvx[i]));
    Cy[3+i]=I*k/4.0*md->bd.wg[i]*(k*H0[i]*tD[i]*vDy[i]-2.0*i_R[i]*H1[i]*(tD[i]*vDy[i]-0.5*nvy[i]));
  }
  CC[6]=0.0;
  Cx[6]=0.0;
  Cy[6]=0.0;
  for(i=0;i<3;i++){
    CC[6]+=md->bd.wg[i]*i_R[i]*tD[i];
      Cx[6]+=md->bd.wg[i]*i_R[i]*i_R[i]*(tD[i]*vDx[i]-0.5*nvx[i]);
      Cy[6]+=md->bd.wg[i]*i_R[i]*i_R[i]*(tD[i]*vDy[i]-0.5*nvy[i]);
  }
  CC[6]*=1.0/(2.0*M_PI);
  Cx[6]*=1.0/M_PI;
  Cy[6]*=1.0/M_PI;

  if(sig==1){
    convert_CC(CC);
    convert_Cv(Cx);
    convert_Cv(Cy);
  }
}

void d_coef_GL(double complex *Cv,double *vt,double *rt,int s,double complex k,DOMD *md)
{
  void convert_Cv(double complex *Cv);
  
  double complex H0[3],H1[3],kR;
  double X,Y,R,dx,dy,i_R[3],tD[3],vD[3],aJ[3],nv[3];
  double wt[3]={WA_GL,WA_GL,W0_GL};
  int i,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<3;i++){
    X=md->bd.x[s][i]-rt[0];
    Y=md->bd.y[s][i]-rt[1];
    R=sqrt(X*X+Y*Y);
    i_R[i]=1.0/R;
    dx=md->bd.dx[s][i];
    dy=md->bd.dy[s][i];
    tD[i]=X*i_R[i]*dy-Y*i_R[i]*dx;
    vD[i]=X*i_R[i]*vt[0]+Y*i_R[i]*vt[1];
    aJ[i]=sqrt(dx*dx+dy*dy);
    nv[i]=dy*vt[0]-dx*vt[1];
    kR=k*R;
    chankel1_01(kR,&(H0[i]),&(H1[i]));
  }

  for(i=0;i<3;i++){
    Cv[  i]= I*k/4.0*wt[i]*H1[i]*vD[i]*aJ[i];
    Cv[3+i]= I*k/4.0*wt[i]*(k*H0[i]*tD[i]*vD[i]-2.0*i_R[i]*H1[i]*(tD[i]*vD[i]-0.5*nv[i]));
  }
  Cv[6]=0.0;
  for(i=0;i<3;i++){
    Cv[6]+=wt[i]*i_R[i]*tD[i];
  }
  Cv[6]*=1.0/(2.0*M_PI); // F_{t,s}
  if(sig==1) convert_Cv(Cv);
}

void d_coef_GK(double complex *CvGK,double *vt,double *rt,int s,double complex k,DOMD *md)
{
  void convert_Cv(double complex *Cv);
  
  double complex H0[7],H1[7],kR;
  double X,Y,R,dx,dy,i_R[7],tD[7],vD[7],aJ[7],nv[7];
  int i,j,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<7;i++){
    X=md->bd.x[s][i]-rt[0];
    Y=md->bd.y[s][i]-rt[1];
    R=sqrt(X*X+Y*Y);
    i_R[i]=1.0/R;
    dx=md->bd.dx[s][i];
    dy=md->bd.dy[s][i];
    tD[i]=X*i_R[i]*dy-Y*i_R[i]*dx;
    vD[i]=X*i_R[i]*vt[0]+Y*i_R[i]*vt[1];
    aJ[i]=sqrt(dx*dx+dy*dy);
    nv[i]=dy*vt[0]-dx*vt[1];
    kR=k*R;
    chankel1_01(kR,&(H0[i]),&(H1[i]));
  }

  // gauss-kronrod
  for(i=0;i<7;i++) CvGK[i]=0.0;
  for(j=0;j<3;j++){
    for(i=0;i<7;i++){
      CvGK[  j]+=md->bd.wk[i]*H1[i]*vD[i]*md->bd.M[j][i]*aJ[i];
      CvGK[3+j]+=md->bd.wk[i]*(k*H0[i]*tD[i]*vD[i]-2.0*i_R[i]*H1[i]*(tD[i]*vD[i]-0.5*nv[i]))*md->bd.M[j][i];
      if(j==0){
        CvGK[6]+=md->bd.wk[i]*i_R[i]*tD[i]; //F
      }
    }
  }
  for(j=0;j<3;j++){
    CvGK[  j]*=I*k/4.0;
    CvGK[3+j]*=I*k/4.0;
  }
  CvGK[6]*=0.5/M_PI;

  if(sig==1) convert_Cv(CvGK);
}

void d_coef_HP(double complex *Cv,double *vt,double *rt,int s,double complex k,DOMD *md)
{
  void convert_Cv(double complex *Cv);
  
  double complex H0,H1,kR;
  double X,Y,R,i_R,dr[2],r[2],rs[2][3],tD,vD,aJ,nv;
  int i,j,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<3;i++){
    rs[0][i]=md->bd.x[s][i];
    rs[1][i]=md->bd.y[s][i];
  }

  for(i=0;i<7;i++) Cv[i]=0.0;
  for(i=0;i<GLH;i++){
    rs_eta(r,rs,md->bd.xh[i]);
    drs_eta(dr,rs,md->bd.xh[i]);
    X=r[0]-rt[0];
    Y=r[1]-rt[1];
    R=sqrt(X*X+Y*Y);
    i_R=1.0/R;
    tD=X*i_R*dr[1]-Y*i_R*dr[0];
    vD=X*i_R*vt[0]+Y*i_R*vt[1];
    aJ=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
    nv=dr[1]*vt[0]-dr[0]*vt[1];
    kR=k*R;
    chankel1_01(kR,&H0,&H1);
    
    for(j=0;j<3;j++){ 
      Cv[  j]+=md->bd.wh[i]*H1*vD*Mn(j,md->bd.xh[i])*aJ;
      Cv[3+j]+=md->bd.wh[i]*(k*H0*tD*vD-2.0*i_R*H1*(tD*vD-0.5*nv))*Mn(j,md->bd.xh[i]);
      if(j==0) Cv[6]+=md->bd.wh[i]*i_R*tD; //F
    }
  }

  for(j=0;j<3;j++){
    Cv[  j]*=I*k/4.0;
    Cv[3+j]*=I*k/4.0;
  }
  Cv[6]*=0.5/M_PI;

  if(sig==1) convert_Cv(Cv);
}

int d_coef_NV(double complex *CvGK,double *vt,double *rt,int s,double complex k,DOMD *md)
{
  void convert_Cv(double complex *Cv);
  
  double complex H0[7],H1[7],CvGL[7],kR;
  double X,Y,R,dx,dy,i_R[7],tD[7],vD[7],aJ[7],nv[7],CD;
  int i,j,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  for(i=0;i<7;i++){
    X=md->bd.x[s][i]-rt[0];
    Y=md->bd.y[s][i]-rt[1];
    R=sqrt(X*X+Y*Y);
    i_R[i]=1.0/R;
    dx=md->bd.dx[s][i];
    dy=md->bd.dy[s][i];
    tD[i]=X*i_R[i]*dy-Y*i_R[i]*dx;
    vD[i]=X*i_R[i]*vt[0]+Y*i_R[i]*vt[1];
    aJ[i]=sqrt(dx*dx+dy*dy);
    nv[i]=dy*vt[0]-dx*vt[1];
    kR=k*R;
    chankel1_01(kR,&(H0[i]),&(H1[i]));
  }

  // gauss-kronrod
  for(i=0;i<7;i++) CvGK[i]=0.0;
  for(j=0;j<3;j++){
    for(i=0;i<7;i++){
      CvGK[  j]+=md->bd.wk[i]*H1[i]*vD[i]*md->bd.M[j][i]*aJ[i];
      CvGK[3+j]+=md->bd.wk[i]*(k*H0[i]*tD[i]*vD[i]-2.0*i_R[i]*H1[i]*(tD[i]*vD[i]-0.5*nv[i]))*md->bd.M[j][i];
      if(j==0){
        CvGK[6]+=md->bd.wk[i]*i_R[i]*tD[i]; //F
      }
    }
  }
  for(j=0;j<3;j++){
    CvGK[  j]*=I*k/4.0;
    CvGK[3+j]*=I*k/4.0;
  }
  CvGK[6]*=0.5/M_PI;

  // gauss-legendre
  for(i=0;i<3;i++){
    CvGL[  i]= I*k/4.0*md->bd.wg[i]*H1[i]*vD[i]*aJ[i];
    CvGL[3+i]= I*k/4.0*md->bd.wg[i]*(k*H0[i]*tD[i]*vD[i]-2.0*i_R[i]*H1[i]*(tD[i]*vD[i]-0.5*nv[i]));
  }
  CvGL[6]=0.0; // F
  for(i=0;i<3;i++){
    CvGL[6]+=md->bd.wg[i]*i_R[i]*tD[i];
  }
  CvGL[6]*=0.5/M_PI;

  // check
  j=0;
  for(i=0;i<7;i++){
    CD=2.0*cabs(CvGK[i]-CvGL[i])/cabs(CvGK[i]+CvGL[i]);
    if(CD>IEPS) j++;
  }
  if(j!=0){
    if(sig==0) d_coef_HP(CvGK,vt,rt, s,k,md);
    else       d_coef_HP(CvGK,vt,rt,-s,k,md);
    return 1;
  }
  else{
    if(sig==1) convert_Cv(CvGK);
    return 0;
  }
}

void d_coef_DE(double complex *Cv,double *vt,double *rt,int s,double complex k,DOMD *md)
{
  int dc_check_vector1(TDATA *td,double *vt);
  int dc_check_vector2(TDATA *td,double *vt);
  int dc_check_vector3(TDATA *td);
  double complex dc_dsG(double eta,void *tmp);
  double complex dc_dsH(double eta,void *tmp);
  double dc_sF(double eta,void *tmp);
  void convert_Cv(double complex *Cv);
  
  TDATA td;
  double complex dG0,dG1,dG2,dH0,dH1,dH2,F;
  double err0;
  int i,sig;

  if(s>0) sig=0;
  else{
    s=-s;
    sig=1;
  }

  td.k=k;
  td.rt[0]=rt[0];
  td.rt[1]=rt[1];
  td.vt[0]=vt[0];
  td.vt[1]=vt[1];

  for(i=0;i<3;i++){
    td.rs[0][i]=md->bd.x[s][i];
    td.rs[1][i]=md->bd.y[s][i];
  }

  if(dc_check_vector1(&td,vt)==0){
    td.type=0;
    dG0=deintz(dc_dsG,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. d_coef_DE(), dc_dsG()! Exit...\n");  exit(1);  }
    td.type=1;
    dG1=deintz(dc_dsG,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. d_coef_DE(), dc_dsG()! Exit...\n");  exit(1);  }
    td.type=2;
    dG2=deintz(dc_dsG,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. d_coef_DE(), dc_dsG()! Exit...\n");  exit(1);  }
  }
  else {
    dG0=0.0;
    dG1=0.0;
    dG2=0.0;
  }

  if(dc_check_vector2(&td,vt)==0){
    td.type=0;
    dH0=deintz(dc_dsH,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. d_coef_DE(), dc_dsH()! Exit...\n");  exit(1);  }
    td.type=1;
    dH1=deintz(dc_dsH,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. d_coef_DE(), dc_dsH()! Exit...\n");  exit(1);  }
    td.type=2;
    dH2=deintz(dc_dsH,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. d_coef_DE(), dc_dsH()! Exit...\n");  exit(1);  }

  }
  else {
    dH0=0.0;
    dH1=0.0;
    dH2=0.0;
  }

  if(dc_check_vector3(&td)==0){
    F=deintd(dc_sF,-1.0,1.0,&td,DEPS,&err0);
    if(err0<0.0){ printf("DE integration error. d_coef_DE(), dc_F()! Exit...\n");  exit(1);  }
  }
  else F=0.0;

  Cv[0]=I*td.k/4.0*dG0;
  Cv[1]=I*td.k/4.0*dG1;
  Cv[2]=I*td.k/4.0*dG2;
  Cv[3]=-I*td.k*td.k/4.0*dH0;
  Cv[4]=-I*td.k*td.k/4.0*dH1;
  Cv[5]=-I*td.k*td.k/4.0*dH2;
  Cv[6]=0.5*F/M_PI;
  if(sig==1) convert_Cv(Cv);
}

void d_coef_bd_node_t(double complex *Cv,int s,int sn,double complex k,DOMD *md)
{
  double iGF1(double eta,void *tmp);
  double iGF1mf(double eta,void *tmp);
  double complex ipvGF(void *tmp);
  double complex iGF2m(double xi,void *tmp);
  double iHF1(double eta,void *tmp);
  double complex iHF2m(double xi,void *tmp);
  double complex iHF3m(double xi,void *tmp);
  double complex ipvGF(void *td);
  double isF(double eta,void *tmp);
  void convert_Cv(double complex *Cv);
  void convert_tv(double complex *Cv);
  
  TDATA td;
  double complex tc1,tc2;
  double t1;
  double err,err0;
  int i,t,sig,csn;

  double teps=DEPS;

  // check sn;
  if(sn<0 || sn>2){
    printf("d_coef_bd_node_t(), node id error. node id =%d. Exit...\n",sn);
    exit(1);
  }

  if(s>0){
     sig=0;
    csn=sn;
  }
  else {
    sig=1;
    s=-s;
    if(sn==0) csn=1;
    else if(sn==1) csn=0;
    else csn=sn;
  }

  td.node=csn;
  td.k=k;
  td.eta_t=md->bd.et[csn];
  for(i=0;i<3;i++){
    td.rs[0][i]=md->bd.x[s][i];
    td.rs[1][i]=md->bd.y[s][i];
  }
  td.drt[0]=md->bd.dx[s][csn];
  td.drt[1]=md->bd.dy[s][csn];
  td.c1=sqrt(td.drt[0] * td.drt[0] + td.drt[1] * td.drt[1]); // |J_s(\eta_t)|
  td.c2=0.5*i_A_GL*i_A_GL*i_A_GL*( (td.rs[0][1]-td.rs[0][2])*td.rs[1][0]
                                  +(td.rs[0][2]-td.rs[0][0])*td.rs[1][1]
                                  +(td.rs[0][0]-td.rs[0][1])*td.rs[1][2]); // C_s
  td.cr[0]=0.5*i_A_GL*i_A_GL*(td.rs[0][0]+td.rs[0][1]-2.0*td.rs[0][2]); // C_x
  td.cr[1]=0.5*i_A_GL*i_A_GL*(td.rs[1][0]+td.rs[1][1]-2.0*td.rs[1][2]); // C_y
  // dG/dt
  for(t=0;t<3;t++){
    td.type=t;
    if(csn==t){ // principal value
      t1=deintd(iGF1mf,-1.0,1.0,&td,teps,&err);
      if(err<0.0 && fabs(err)>teps){
        printf("DE integration error. d_coef_bd_node_t() iGF1mf\n");
        printf("error code %15.14e. Exit...\n",err);
        exit(1);
      }
       tc1=-deintz(iGF2m,0.0,-1.0-td.eta_t,&td,teps,&err);
      tc2= deintz(iGF2m,0.0, 1.0-td.eta_t,&td,teps,&err0);
      if(err<0.0 || err0<0.0){
        printf("DE integration error. d_coef_bd_node_t() iGF2m\n");
        printf("error code %15.14e, %15.14e. Exit...\n",err,err0);
        exit(1);
      }
      Cv[t]=1.0/(2.0*M_PI)*(t1+ipvGF(&td))-td.k*td.k/(4.0*M_PI)*(tc1+tc2);
    }
    else { // normal integral
      t1=deintd(iGF1,-1.0,1.0,&td,teps,&err);
      if(err<0.0){
        printf("DE integration error. d_coef_bd_node_t() iGF1\n");
         printf("error code %15.14e. Exit...\n",err);
        exit(1);
      }
       tc1=-deintz(iGF2m,0.0,-1.0-td.eta_t,&td,teps,&err);
      tc2= deintz(iGF2m,0.0, 1.0-td.eta_t,&td,teps,&err0);
      if(err<0.0 || err0<0.0){
        printf("DE integration error. d_coef_bd_node_t() iGF2m\n");
        printf("error code %15.14e, %15.14e. Exit...\n",err,err0);
        exit(1);
      }
      Cv[t]=1.0/(2.0*M_PI)*t1-td.k*td.k/(4.0*M_PI)*(tc1+tc2);
     }
  }
  // dH/dt
  for(t=0;t<3;t++){
    td.type=t;
    
     t1=td.c2/(M_PI*td.c1)*deintd(iHF1,-1.0,1.0,&td,teps,&err);
    if(err<0.0){
      printf("DE integration error. d_coef_bd_node_t() iHF1\n");
      printf("error code %15.14e. Exit...\n",err);
      exit(1);
    }
    Cv[t+3]=t1;
     tc1=-deintz(iHF2m,0.0,-1.0-td.eta_t,&td,teps,&err);
    tc2= deintz(iHF2m,0.0, 1.0-td.eta_t,&td,teps,&err0);
    if(err<0.0 || err0<0.0){
      printf("DE integration error. d_coef_bd_node_t() iHF2m\n");
      printf("error code %15.14e, %15.14e. Exit...\n",err,err0);
      exit(1);
    }
    Cv[t+3]+=-td.k*td.k*td.c2/(2.0*M_PI*td.c1)*(tc1+tc2);
     tc1=-deintz(iHF3m,0.0,-1.0-td.eta_t,&td,teps,&err);
    tc2= deintz(iHF3m,0.0, 1.0-td.eta_t,&td,teps,&err0);
    if(err<0.0 || err0<0.0){
      printf("DE integration error. d_coef_bd_node_t() iHF3m\n");
      printf("error code %15.14e, %15.14e. Exit...\n",err,err0);
      exit(1);
    }
    Cv[t+3]+=I*td.k*td.k*td.c2/(4.0*td.c1)*(tc1+tc2);
  }
  // F
  Cv[6]=td.c2/(2.0*M_PI)*deintd(isF,-1.0,1.0,&td,teps,&err);
  if(err<0.0){ printf("DE integration error d_coef_bd_node_t(), isF. Exit...\n");  exit(1);  }
  
  if(sig==1){
    convert_Cv(Cv);
    convert_tv(Cv); // convert tangential vector direction
  }
}

//------------------------------------------------------------------
int dc_check_vector1(TDATA *td,double *vt)
{
  double A[2],B[2],nr0,nr1;

  A[0]=td->rs[0][0]-td->rt[0];
  A[1]=td->rs[1][0]-td->rt[1];
  nr0=1.0/sqrt(A[0]*A[0]+A[1]*A[1]);
  A[0]*=nr0;
  A[1]*=nr0;
  B[0]=td->rs[0][1]-td->rt[0];
  B[1]=td->rs[1][1]-td->rt[1];
  nr1=1.0/sqrt(B[0]*B[0]+B[1]*B[1]);
  B[0]*=nr1;
  B[1]*=nr1;
  if(fabs(A[0]*B[1]-A[1]*B[0])>MEPS) return 0;
  B[0]=td->rs[0][2]-td->rt[0];
  B[1]=td->rs[1][2]-td->rt[1];
  nr1=1.0/sqrt(B[0]*B[0]+B[1]*B[1]);
  B[0]*=nr1;
  B[1]*=nr1;
  if(fabs(A[0]*B[1]-A[1]*B[0])>MEPS) return 0;
  if(fabs(A[0]*vt[0]+A[1]*vt[1])>MEPS) return 0;

  return 1;
}

int dc_check_vector2(TDATA *td,double *vt)
{
  double A[2],B[2],nr0,nr1;
  A[0]=td->rs[0][0]-td->rt[0];
  A[1]=td->rs[1][0]-td->rt[1];
  nr0=1.0/sqrt(A[0]*A[0]+A[1]*A[1]);
  A[0]*=nr0;
  A[1]*=nr0;
  B[0]=td->rs[0][1]-td->rt[0];
  B[1]=td->rs[1][1]-td->rt[1];
  nr1=1.0/sqrt(B[0]*B[0]+B[1]*B[1]);
  B[0]*=nr1;
  B[1]*=nr1;
  if(fabs(A[0]*B[1]-A[1]*B[0])>MEPS) return 0;
  B[0]=td->rs[0][2]-td->rt[0];
  B[1]=td->rs[1][2]-td->rt[1];
  nr1=1.0/sqrt(B[0]*B[0]+B[1]*B[1]);
  B[0]*=nr1;
  B[1]*=nr1;
  if(fabs(A[0]*B[1]-A[1]*B[0])>MEPS) return 0;
  if(fabs(A[0]*vt[1]-A[1]*vt[0])>MEPS) return 0;

  return 1;
}

int dc_check_vector3(TDATA *td)
{
  double A[2],B[2],nr0,nr1;
  A[0]=td->rs[0][0]-td->rt[0];
  A[1]=td->rs[1][0]-td->rt[1];
  nr0=1.0/sqrt(A[0]*A[0]+A[1]*A[1]);
  A[0]*=nr0;
  A[1]*=nr0;
  B[0]=td->rs[0][1]-td->rt[0];
  B[1]=td->rs[1][1]-td->rt[1];
  nr1=1.0/sqrt(B[0]*B[0]+B[1]*B[1]);
  B[0]*=nr1;
  B[1]*=nr1;
  if(fabs(A[0]*B[1]-A[1]*B[0])>MEPS) return 0;
  B[0]=td->rs[0][2]-td->rt[0];
  B[1]=td->rs[1][2]-td->rt[1];
  nr1=1.0/sqrt(B[0]*B[0]+B[1]*B[1]);
  B[0]*=nr1;
  B[1]*=nr1;
  if(fabs(A[0]*B[1]-A[1]*B[0])>MEPS) return 0;

  return 1;
}

void convert_Cv(double complex *Cv)
{
  double complex tdG,tdH;
  tdG=Cv[0];
  Cv[0]=Cv[1];
  Cv[1]=tdG;
  tdH=Cv[3];
  Cv[3]=-Cv[4];
  Cv[4]=-tdH;
  Cv[5]=-Cv[5];
  Cv[6]=-Cv[6];
}

double complex dc_dsG(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex H1,kR;
  double r[2],dr[2],X,Y,R,i_R,vD,aJ;

  rs_eta(r,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  aJ=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  X=r[0]-td->rt[0];
  Y=r[1]-td->rt[1];
  R=sqrt(X*X+Y*Y);
  i_R=1.0/R;
  vD=X*i_R*td->vt[0]+Y*i_R*td->vt[1];
  kR=td->k*R;
  H1=chankel1_1(kR);
  return H1*vD*Mn(td->type,eta)*aJ;
}

double complex dc_dsH(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex H0,H1,kR;
  double r[2],dr[2],X,Y,R,i_R,tD,vD,nv;

  rs_eta(r,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  X=r[0]-td->rt[0];
  Y=r[1]-td->rt[1];
  R=sqrt(X*X+Y*Y);
  i_R=1.0/R;
  tD=X*i_R*dr[1]-Y*i_R*dr[0];
  vD=X*i_R*td->vt[0]+Y*i_R*td->vt[1];
  nv=dr[1]*td->vt[0]-dr[0]*td->vt[1];
  kR=td->k*R;
  chankel1_01(kR,&H0,&H1);

  return -(H0*tD*vD-2.0/kR*H1*(tD*vD-0.5*nv))*Mn(td->type,eta);
}

double dc_dsF(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double r[2],dr[2],X,Y,R,i_R,tD,vD,nv;

  rs_eta(r,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  X=r[0]-td->rt[0];
  Y=r[1]-td->rt[1];
  R=sqrt(X*X+Y*Y);
  i_R=1.0/R;
  tD=X*i_R*dr[1]-Y*i_R*dr[0];
  vD=X*i_R*td->vt[0]+Y*i_R*td->vt[1];
  nv=dr[1]*td->vt[0]-dr[0]*td->vt[1];
  return i_R*i_R*(tD*vD-0.5*nv);
}


double dc_sF(double eta,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double r[2],dr[2],X,Y,i_R,tD;

  rs_eta(r,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  X=r[0]-td->rt[0];
  Y=r[1]-td->rt[1];
  i_R=1.0/sqrt(X*X+Y*Y);
  tD=X*i_R*dr[1]-Y*i_R*dr[0];
  return i_R*tD;
}

double iGF1(double eta,void *tmp)
{
  double Mn_de(int type,int node,double eta);
  
  TDATA *td=(TDATA *)tmp;
  double dr[2],drR[2],aJ,aJR2;

  drs_eta(dr ,td->rs,eta);
  drs_eta(drR,td->rs,0.5*(eta+td->eta_t));
  aJ =sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  aJR2=drR[0]*drR[0]+drR[1]*drR[1];

  return Mn_de(td->type,td->node,eta)*(drR[0]*td->drt[0]+drR[1]*td->drt[1])/aJR2*(aJ/td->c1);
}

double iGF1mf(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double dr[2],drR[2],aJ,aJR2,rj;

  drs_eta(dr ,td->rs,eta);
  drs_eta(drR,td->rs,0.5*(eta+td->eta_t));
  aJ =sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  aJR2=drR[0]*drR[0]+drR[1]*drR[1];
  rj=aJ/td->c1;

  return (eta-td->eta_t)*Mn(td->type,eta)*(drR[0]*(td->drt[0]*rj-drR[0])+drR[1]*(td->drt[1]*rj-drR[1]))/aJR2;
}

double complex iGF2m(double xi,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex O1,kR;
  double dr[2],drR[2],aJ,aJR,eta;

  eta=xi+td->eta_t;
  drs_eta(dr ,td->rs,eta);
  drs_eta(drR,td->rs,0.5*(eta+td->eta_t));
  aJ =sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  aJR=sqrt(drR[0]*drR[0]+drR[1]*drR[1]);
  kR=td->k*fabs(xi)*aJR;
  O1=chankel1_O_1(kR);

  return xi*Mn(td->type,eta)*O1*(drR[0]*td->drt[0]+drR[1]*td->drt[1])*(aJ/td->c1);
}

double iHF1(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double drR[2],aJR2;

  drs_eta(drR,td->rs,0.5*(eta+td->eta_t));
  aJR2=drR[0]*drR[0]+drR[1]*drR[1];

  return (td->cr[0]*drR[0]+td->cr[1]*drR[1])/aJR2*Mn(td->type,eta);
}

double complex iHF2m(double xi,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex O1,kR;
  double drR[2],aJR2,eta;

  eta=xi+td->eta_t;
  drs_eta(drR,td->rs,0.5*(eta+td->eta_t));
  aJR2=drR[0]*drR[0]+drR[1]*drR[1];
  kR=td->k*fabs(xi)*sqrt(aJR2);
  O1=chankel1_O_1(kR);

  return xi*xi*O1*(td->cr[0]*drR[0]+td->cr[1]*drR[1])/aJR2*Mn(td->type,eta);
}

double complex iHF3m(double xi,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex H0,kR;
  double drR[2],aJR2,eta;

  eta=xi+td->eta_t;
  drs_eta(drR,td->rs,0.5*(eta+td->eta_t));
  aJR2=drR[0]*drR[0]+drR[1]*drR[1];
  kR=td->k*fabs(xi)*sqrt(aJR2);
  H0=chankel1_0(kR);

  return xi*H0*(drR[0]*td->drt[0]+drR[1]*td->drt[1])/aJR2*Mn(td->type,eta);
}

double complex ipvGF(void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double etpa,lnet;

  etpa=td->eta_t*i_A_GL;
  lnet=log((1.0-td->eta_t)/(1.0+td->eta_t));

  if(td->type==0){
    return -i_A_GL*(1.0-etpa)*(0.5*td->eta_t*lnet+1.0);
  }
  else if(td->type==1){
    return i_A_GL*(1.0+etpa)*(0.5*td->eta_t*lnet+1.0);
  }
  else if(td->type==2){
    return (1.0+etpa)*(1.0-etpa)*lnet-2.0*i_A_GL*etpa;
  }
  else {
    printf("pb_dcoef.c, type number error in ipvGF(). type=%d. Exit...\n",td->type);
    exit(1);
  }
  return 0;
}

double isF(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double drR[2],aJR2;

  drs_eta(drR,td->rs,0.5*(eta+td->eta_t));
  aJR2=drR[0]*drR[0]+drR[1]*drR[1];
  return 1.0/aJR2;
}

double Mn_de(int type,int node,double eta)
{
  if(type==0){
    if(node==1)             return 0.5*eta*i_A_GL*i_A_GL;
    else if(node==2)return -0.5*(1.0-eta*i_A_GL)*i_A_GL;
  }
  else if(type==1){
    if(node==0)             return 0.5*eta*i_A_GL*i_A_GL;
    else if(node==2)return 0.5*(1.0+eta*i_A_GL)*i_A_GL;
  }
  else if(type==2){
    if(node==0)             return (1.0-eta*i_A_GL)*i_A_GL;
    else if(node==1)return -(1.0+eta*i_A_GL)*i_A_GL;
  }
  printf("pb_dcoef.c, Mn_de() parameter error. type=%d, node=%d. Exit...\n",type,node);
  exit(1);
}

void convert_tv(double complex *Cv)
{
  int i;
  for(i=0;i<6;i++) Cv[i]*=-1.0;
}
