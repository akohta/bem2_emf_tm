/*
 * d2tm_solve_bieq.c
 *
 *  Created on: Sep 18, 2018
 *      Author: ohta
 */
#include "bem2_emf_tm.h"

typedef struct coef_matrix{
  int M; // total domain (material) number
	int N; // total matrix size
	int tt,pad;
	lapack_complex_double *A,*B; // total coefficient matrix
	lapack_int *piv;

	char **tgfn; // G matrix name of each domain
	char **thfn; // H matrix name of each domain
}CMD;

void solve_bieq(DOMD *md)
{
  void init_cmd(CMD *cm);
  void create_matrix(CMD *cm,DOMD *md);
  void solve_mkl_lapack(CMD *cm);
  void calc_bv_u_dudn(CMD *cm,DOMD *md);
  void calc_bv_vw(DOMD *md);
  void malloc_cmatrix_cmd(int N,CMD *cm);
  void create_amatrix_subdomain(int did,CMD *cm,DOMD *md);
  void create_bmatrix_subdomain(int did,CMD *cm,double complex **u);
  void LU_fact_amatrix(CMD *cm);
  void solve_LU(CMD *cm);
  void input_bmatrix(CMD *cm,DOMD *md,double complex **dudn);
  void mfree_cmd(CMD *cm);
  
  time_t start,end;
  CMD cm;
  int d;

  printf("solve electromagnetic field boundary value \n");
  time(&start);

  printf("  coefficient matrix          "); fflush(stdout);
  cm.M=md->MN;
  cm.N=md->bd.Ne*3*2;
  init_cmd(&cm);
  create_matrix(&cm,md);
  printf(" finished\n");

  printf("  solve boundary value Hz     "); fflush(stdout);
  solve_mkl_lapack(&cm);
  calc_bv_u_dudn(&cm,md);
  printf(" finished\n");

  printf("  boundary value Ex Ey        "); fflush(stdout);
  calc_bv_vw(md);
  printf(" finished\n");

  for(d=0;d<=md->MN;d++){
    if(md->bd.sb[d].Ne==0) continue;
    printf("  boundary value of domain %-3d",d); fflush(stdout);
    malloc_cmatrix_cmd(md->bd.sb[d].Ne*3,&cm);
    create_amatrix_subdomain(d,&cm,md);
    LU_fact_amatrix(&cm);
    // v
    create_bmatrix_subdomain(d,&cm,md->bd.sb[d].v);
    solve_LU(&cm);
    input_bmatrix(&cm,md,md->bd.sb[d].dvdn);
    // w
    create_bmatrix_subdomain(d,&cm,md->bd.sb[d].w);
    solve_LU(&cm);
    input_bmatrix(&cm,md,md->bd.sb[d].dwdn);
    printf(" finished\n");
  }

  mfree_cmd(&cm);
  time(&end);
  printf("elapsed time : %g (sec)\n",difftime(end,start));
}

///////////////////////////////////////////////////////////////////////
void init_cmd(CMD *cm)
{
  int i,MN,N;

  MN=cm->M;

  cm->tgfn=(char **)m_alloc2(MN+1,sizeof(char *),"init_cmd(),cm->tgfn");
  cm->thfn=(char **)m_alloc2(MN+1,sizeof(char *),"init_cmd(),cm->thfn");
  for(i=0;i<=MN;i++){
    cm->tgfn[i]=(char *)m_alloc2(16,sizeof(char ),"init_cmd(),tgfn[i]");
    cm->thfn[i]=(char *)m_alloc2(16,sizeof(char ),"init_cmd(),thfn[i]");
    sprintf(cm->tgfn[i],"tmpG_%05d.dat",i);
    sprintf(cm->thfn[i],"tmpH_%05d.dat",i);
  }

  // memory allocation
  N=cm->N;
  cm->A=(lapack_complex_double *)m_alloc2(N*N,sizeof(lapack_complex_double),"init_cmd(),cm->A");
  cm->B=(lapack_complex_double *)m_alloc2(  N,sizeof(lapack_complex_double),"init_cmd(),cm->B");
  cm->piv=(lapack_int *)m_alloc2(N,sizeof(lapack_int),"init_cmd(),piv");
}

void mfree_cmd(CMD *cm)
{
  int i;

  // delete temporary file
  for(i=0;i<=cm->M;i++){
    remove(cm->tgfn[i]);
    remove(cm->thfn[i]);
  }

  // free memory
  for(i=0;i<=cm->M;i++){
    free(cm->tgfn[i]);    free(cm->thfn[i]);
  }
  free(cm->tgfn);  free(cm->thfn);
  cm->M=0;

  free(cm->A);
  free(cm->B);
  free(cm->piv);
  cm->N=0;
}

void create_matrix(CMD *cm,DOMD *md)
{
  void create_matrix_domain(int did,CMD *cm,DOMD *md);
  
  int i;

  cm->tt=0;
  for(i=0;i<=md->MN;i++) create_matrix_domain(i,cm,md);
}

void create_matrix_domain(int did,CMD *cm,DOMD *md)
{
  void swap_CC(double complex *CC);
  
  FILE *fg,*fh;
  double complex *tG,*tH,CC[7],k,ctp,ceps2;
  double rt[2],F;
  int i,j,jc,t,td,atd,s,sd,asd,Ne,N;

  Ne=md->bd.sb[did].Ne;
  N=3*Ne;

  tG=(double complex *)m_alloc2(N,sizeof(double complex),"create_matrix_domain(),tG");
  tH=(double complex *)m_alloc2(N,sizeof(double complex),"create_matrix_domain(),tH");
  if((fg=fopen(cm->tgfn[did],"wb"))==NULL){    printf("Can not open the %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"wb"))==NULL){    printf("Can not open the %s file.\n",cm->thfn[did]);    exit(1);  }

  k=md->wd.k0*md->n[did];
  for(t=1;t<=Ne;t++){
    td=md->bd.sb[did].eid[t];
    atd=abs(td);
    for(j=0;j<3;j++){
      if(td>0) jc=j;
      else {
        if     (j==0) jc=1;
        else if(j==1) jc=0;
        else          jc=2;
      }
      rt[0]=md->bd.x[atd][jc];   rt[1]=md->bd.y[atd][jc];

      F=0.0;
      #pragma omp parallel for schedule(dynamic) reduction(+:F) private(sd,asd,CC,i,ctp,ceps2)
      for(s=1;s<=Ne;s++){
        sd=md->bd.sb[did].eid[s];
        asd=abs(sd);
        if(t!=s) coefficient_NV(CC,rt,sd,k,md);
        else coefficient_bd_eta(CC,md->bd.et[j],sd,k,md);

        for(i=0;i<3;i++){
          tG[(s-1)*3+i]=CC[  i];
          tH[(s-1)*3+i]=CC[3+i];
        }

        if(did==md->bd.md[asd]){
          for(i=0;i<3;i++){
            cm->A[cm->tt*cm->N+cm->N/2*0+(asd-1)*3+i].real=creal( CC[3+i]); // H
            cm->A[cm->tt*cm->N+cm->N/2*0+(asd-1)*3+i].imag=cimag( CC[3+i]);
            cm->A[cm->tt*cm->N+cm->N/2*1+(asd-1)*3+i].real=creal(-CC[  i]); // G
            cm->A[cm->tt*cm->N+cm->N/2*1+(asd-1)*3+i].imag=cimag(-CC[  i]);
          }
        }
        else if(did==md->bd.sd[asd]){
          swap_CC(CC);
          ceps2=cpow(md->n[did]/md->n[md->bd.md[asd]],2);
          for(i=0;i<3;i++){
            cm->A[cm->tt*cm->N+cm->N/2*0+(asd-1)*3+i].real=creal( CC[3+i]);
            cm->A[cm->tt*cm->N+cm->N/2*0+(asd-1)*3+i].imag=cimag( CC[3+i]);
            cm->A[cm->tt*cm->N+cm->N/2*1+(asd-1)*3+i].real=creal( CC[  i]*ceps2);
            cm->A[cm->tt*cm->N+cm->N/2*1+(asd-1)*3+i].imag=cimag( CC[  i]*ceps2);
            if(md->bd.md[asd]==0){
              ctp=-CC[3+i]*md->bd.ui[asd][i]-ceps2*CC[i]*md->bd.duidn[asd][i];
              #pragma omp atomic
              cm->B[cm->tt].real+=creal(ctp);
              #pragma omp atomic
              cm->B[cm->tt].imag+=cimag(ctp);
            }
          }
        }
        else {
          printf("error in create_matrix_domain()! Exit...\n");
          printf("did=%d,md=%d,sd=%d\n",did,md->bd.md[asd],md->bd.sd[asd]);
          exit(1);
        }
        F+=creal(CC[6]);
      } // end parallel
      if(did==0){
        tH[(t-1)*3+j]+=1.0+F;
        cm->A[cm->tt*cm->N+(atd-1)*3+jc].real+=1.0+F;
      }
      else{
        tH[(t-1)*3+j]+=F;
        cm->A[cm->tt*cm->N+(atd-1)*3+jc].real+=F;
        if(did==md->bd.sd[atd] && 0==md->bd.md[atd]){
          cm->B[cm->tt].real+=creal(-F*md->bd.ui[atd][jc]);
          cm->B[cm->tt].imag+=cimag(-F*md->bd.ui[atd][jc]);
        }
      }
      // fwrite
      fwrite(tG,sizeof(double complex),N,fg);
      fwrite(tH,sizeof(double complex),N,fh);
      cm->tt++;
    }
  }

  fclose(fg);  fclose(fh);
  free(tG);  free(tH);
}

void swap_CC(double complex *CC)
{
  double complex tmp1,tmp2;
  tmp1=CC[0];
  CC[0]=CC[1];
  CC[1]=tmp1;
  tmp2=CC[3];
  CC[3]=CC[4];
  CC[4]=tmp2;
}

void solve_mkl_lapack(CMD *cm)
{
  lapack_int m,n,lda,info,nrhs,ldb;
  
  m=(lapack_int)cm->N;
  n=m;
  lda=m;

  info=LAPACKE_zgetrf(LAPACK_ROW_MAJOR,m,n,cm->A,lda,cm->piv); // LU factorization
  if(info!=0){
    printf("error in solve_mkl_lapack(), LAPACKE_zgetrf(). info=%lld\n",info);
    exit(1);
  }
  nrhs=1;
  ldb=1;
  info=LAPACKE_zgetrs(LAPACK_ROW_MAJOR,'N',n,nrhs,cm->A,lda,cm->piv,cm->B,ldb);
  if(info!=0){
    printf("error in solve_mkl_lapack(), LAPACKE_zgetrs(). info=%lld\n",info);
    exit(1);
  }
}

void calc_bv_u_dudn(CMD *cm,DOMD *md)
{
  double complex ceps2;
  int d,s,eid,aeid;

  for(d=0;d<=md->MN;d++){
    for(s=1;s<=md->bd.sb[d].Ne;s++){
      eid=md->bd.sb[d].eid[s];
      aeid=abs(eid);
      if(md->bd.md[aeid]==d){
        md->bd.sb[d].u[s][0]=cm->B[(aeid-1)*3+0].real+cm->B[(aeid-1)*3+0].imag*I;
        md->bd.sb[d].u[s][1]=cm->B[(aeid-1)*3+1].real+cm->B[(aeid-1)*3+1].imag*I;
        md->bd.sb[d].u[s][2]=cm->B[(aeid-1)*3+2].real+cm->B[(aeid-1)*3+2].imag*I;
        md->bd.sb[d].dudn[s][0]=cm->B[cm->N/2+(aeid-1)*3+0].real+cm->B[cm->N/2+(aeid-1)*3+0].imag*I;
        md->bd.sb[d].dudn[s][1]=cm->B[cm->N/2+(aeid-1)*3+1].real+cm->B[cm->N/2+(aeid-1)*3+1].imag*I;
        md->bd.sb[d].dudn[s][2]=cm->B[cm->N/2+(aeid-1)*3+2].real+cm->B[cm->N/2+(aeid-1)*3+2].imag*I;
      }
      else {
        ceps2=cpow(md->n[d]/md->n[md->bd.md[aeid]],2);
        if(md->bd.md[aeid]!=0){
          md->bd.sb[d].u[s][0]=cm->B[(aeid-1)*3+1].real+cm->B[(aeid-1)*3+1].imag*I;
          md->bd.sb[d].u[s][1]=cm->B[(aeid-1)*3+0].real+cm->B[(aeid-1)*3+0].imag*I;
          md->bd.sb[d].u[s][2]=cm->B[(aeid-1)*3+2].real+cm->B[(aeid-1)*3+2].imag*I;
          md->bd.sb[d].dudn[s][0]=-ceps2*(cm->B[cm->N/2+(aeid-1)*3+1].real+cm->B[cm->N/2+(aeid-1)*3+1].imag*I);
          md->bd.sb[d].dudn[s][1]=-ceps2*(cm->B[cm->N/2+(aeid-1)*3+0].real+cm->B[cm->N/2+(aeid-1)*3+0].imag*I);
          md->bd.sb[d].dudn[s][2]=-ceps2*(cm->B[cm->N/2+(aeid-1)*3+2].real+cm->B[cm->N/2+(aeid-1)*3+2].imag*I);
        }
        else {
          md->bd.sb[d].u[s][0]=cm->B[(aeid-1)*3+1].real+cm->B[(aeid-1)*3+1].imag*I+md->bd.ui[aeid][1];
          md->bd.sb[d].u[s][1]=cm->B[(aeid-1)*3+0].real+cm->B[(aeid-1)*3+0].imag*I+md->bd.ui[aeid][0];
          md->bd.sb[d].u[s][2]=cm->B[(aeid-1)*3+2].real+cm->B[(aeid-1)*3+2].imag*I+md->bd.ui[aeid][2];
          md->bd.sb[d].dudn[s][0]=-ceps2*(cm->B[cm->N/2+(aeid-1)*3+1].real+cm->B[cm->N/2+(aeid-1)*3+1].imag*I+md->bd.duidn[aeid][1]);
          md->bd.sb[d].dudn[s][1]=-ceps2*(cm->B[cm->N/2+(aeid-1)*3+0].real+cm->B[cm->N/2+(aeid-1)*3+0].imag*I+md->bd.duidn[aeid][0]);
          md->bd.sb[d].dudn[s][2]=-ceps2*(cm->B[cm->N/2+(aeid-1)*3+2].real+cm->B[cm->N/2+(aeid-1)*3+2].imag*I+md->bd.duidn[aeid][2]);
        }
      }
    }
  }
}

void calc_bv_vw(DOMD *md)
{
  double complex i_k1,dudt;
  double dx,dy,i_aJ;
  int d,s,i,ic,eid,aeid,sig;

  for(d=0;d<=md->MN;d++){
    i_k1=1.0/(I*md->wd.k0*md->n[d]*md->n[d]);
    #pragma omp parallel for schedule(dynamic) private(eid,aeid,i,ic,sig,dx,dy,i_aJ,dudt)
    for(s=1;s<=md->bd.sb[d].Ne;s++){
      eid=md->bd.sb[d].eid[s];
      aeid=abs(eid);
      for(i=0;i<3;i++){
        if(eid>0) {
          ic=i;
          sig=1;
        }
        else {
          if     (i==0) ic=1;
          else if(i==1) ic=0;
          else          ic=2;
          sig=-1;
        }
        dx=(double)sig*md->bd.dx[aeid][ic];
        dy=(double)sig*md->bd.dy[aeid][ic];

        i_aJ=1.0/sqrt(dx*dx+dy*dy);
        if      (SFDM==0) dudt=dHzdt_bv_node_dbieq(d,s,i,md);
        else if (SFDM==1) dudt=dHzdt_bv_node_ndmtd(d,s,i,md);
        else {
          printf("SFDM parameter error. SFDM=%d. check pb_const.h Exit...\n",SFDM);
          exit(1);
        }
        md->bd.sb[d].v[s][i]=-i_k1*i_aJ*(-md->bd.sb[d].dudn[s][i]*dx+dudt*dy);
        md->bd.sb[d].w[s][i]= i_k1*i_aJ*( md->bd.sb[d].dudn[s][i]*dy+dudt*dx);
      }
    }
  }
}

void malloc_cmatrix_cmd(int N,CMD *cm)
{
  free(cm->A);
  free(cm->B);
  free(cm->piv);

  cm->N=N;
  cm->A=(lapack_complex_double *)m_alloc2(N*N,sizeof(lapack_complex_double),"malloc_cmatrix_cmd(),cm->A");
  cm->B=(lapack_complex_double *)m_alloc2(  N,sizeof(lapack_complex_double),"malloc_cmatrix_cmd(),cm->B");
  cm->piv=(lapack_int *)m_alloc2(N,sizeof(lapack_int),"malloc_cmatrix_cmd(),cm->piv");
}

void create_amatrix_subdomain(int did,CMD *cm,DOMD *md)
{
  FILE *fg;
  double complex *tG;
  int j,s;

  if((fg=fopen(cm->tgfn[did],"rb"))==NULL){    printf("Can not open the %s file.\n",cm->tgfn[did]);    exit(1);  }
  tG=(double complex *)m_alloc2(cm->N,sizeof(double complex),"create_amatrix_subdomain(),tG");
  
  for(j=0;j<cm->N;j++){
    fread(tG,cm->N,sizeof(double complex),fg);
    for(s=0;s<cm->N;s++){
      cm->A[j*cm->N+s].real=creal(tG[s]);
      cm->A[j*cm->N+s].imag=cimag(tG[s]);  
    }
  }
  
  fclose(fg);
  free(tG);
}

void create_bmatrix_subdomain(int did,CMD *cm,double complex **u)
{
  FILE *fh;
  double complex *tH,tc;
  int j,s,sn;

  if((fh=fopen(cm->thfn[did],"rb"))==NULL){      printf("Can not open the %s file. Exit...\n",cm->thfn[did]);   exit(1);        }
  tH=(double complex *)m_alloc2(cm->N,sizeof(double complex),"create_bmatrix_subdomain(),tH");

  for(j=0;j<cm->N;j++){
    cm->B[j].real=0.0;
    cm->B[j].imag=0.0;
    fread(tH,cm->N,sizeof(double complex),fh);
    for(s=1;s<=cm->N/3;s++){
      for(sn=0;sn<3;sn++){
        tc=tH[3*(s-1)+sn]*u[s][sn];
        cm->B[j].real+=creal(tc);
        cm->B[j].imag+=cimag(tc);
      }
    }
  }
  
  fclose(fh);
  free(tH);
}

void LU_fact_amatrix(CMD *cm)
{
  lapack_int m,n,lda,info;
  
  m=(lapack_int)cm->N;
  n=m;
  lda=m;

  info=LAPACKE_zgetrf(LAPACK_ROW_MAJOR,m,n,cm->A,lda,cm->piv); // LU factorization
  if(info!=0){
    printf("error in LU_fact_amatrix(), LAPACKE_zgetrf(). info=%lld\n",info);
    exit(1);
  }
}

void solve_LU(CMD *cm)
{
  lapack_int m,n,lda,info,nrhs,ldb;
  
  m=(lapack_int)cm->N;
  n=m;
  lda=m;
  nrhs=1;
  ldb=1;
  info=LAPACKE_zgetrs(LAPACK_ROW_MAJOR,'N',n,nrhs,cm->A,lda,cm->piv,cm->B,ldb);
    if(info!=0){
    printf("error in solve_LU(), LAPACKE_zgetrs(). info=%lld\n",info);
    exit(1);
  }
}

void input_bmatrix(CMD *cm,DOMD *md,double complex **dudn)
{
  int j,s,sn;

  for(j=0;j<cm->N;j++){
    for(s=1;s<=cm->N/3;s++){
      for(sn=0;sn<3;sn++){
        dudn[s][sn]=cm->B[3*(s-1)+sn].real+cm->B[3*(s-1)+sn].imag*I;
      }
    }
  }
}
