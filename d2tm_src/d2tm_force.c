/*
 * d2tm_force.c
 *
 *  Created on: Sep 18, 2018
 *      Author: ohta
 */
#include "bem2_emf_tm.h"


void force_FN(double *Fr,double *Nz,double *rc,int type,DOMD *md)
{
  void force_FN_GL(double *Fr,double *Nz,double *rc,DOMD *md);
  void force_FN_GK(double *Fr,double *Nz,double *rc,DOMD *md);
  
  if(type==0) force_FN_GL(Fr,Nz,rc,md);
  else        force_FN_GK(Fr,Nz,rc,md);
}

//////////////////////////////////////////////////////////////////
void force_FN_GL(double *Fr,double *Nz,double *rc,DOMD *md)
{
	double complex Hz,Er[2];
	double r[2],dx,dy,Txx,Txy,Tyx,Tyy,ahz2,aex2,aey2,dFx,dFy,epsd,tE,tH;
	int s,i,eid;

	for(i=0;i<2;i++) Fr[i]=0.0;
	*Nz=0.0;

	epsd=md->n[0]*md->n[0];
	for(s=1;s<=md->bd.sb[0].Ne;s++){
		eid=md->bd.sb[0].eid[s];
		for(i=0;i<3;i++){
			r[0]=md->bd.x [eid][i];       r[1]=md->bd.y [eid][i];
			dx  =md->bd.dx[eid][i];       dy  =md->bd.dy[eid][i];
			infd_HE(&Hz,Er,r,&(md->wd));
			Hz   +=md->bd.sb[0].u[s][i];
			Er[0]+=md->bd.sb[0].v[s][i];
			Er[1]+=md->bd.sb[0].w[s][i];
			ahz2=creal(Hz*conj(Hz));
			aex2=creal(Er[0]*conj(Er[0]));
			aey2=creal(Er[1]*conj(Er[1]));
			tE=0.25*epsd*(aex2-aey2);
			tH=0.25*ahz2;
			Txx=-tH+tE;
			Txy=0.50*epsd*creal(Er[0]*conj(Er[1]));
			Tyx=Txy;
			Tyy=-tH-tE;
			dFx=Txx*dy-Txy*dx;
			dFy=Tyx*dy-Tyy*dx;
			Fr[0]+=md->bd.wg[i]*dFx;
			Fr[1]+=md->bd.wg[i]*dFy;
			*Nz+=md->bd.wg[i]*((r[0]-rc[0])*dFy-(r[1]-rc[1])*dFx);
		}
	}

	Fr[0]*=-1.0;
	Fr[1]*=-1.0;
	*Nz*=-1.0;
}

void force_FN_GK(double *Fr,double *Nz,double *rc,DOMD *md)
{
	double _Complex Hz,Er[2],Hsz,Esr[2];
	double r[2],dx,dy,Txx,Txy,Tyx,Tyy,ahz2,aex2,aey2,dFx,dFy,epsd,tE,tH;
	int s,i,eid;

	for(i=0;i<2;i++) Fr[i]=0.0;
	*Nz=0.0;

	epsd=md->n[0]*md->n[0];
	for(s=1;s<=md->bd.sb[0].Ne;s++){
		eid=md->bd.sb[0].eid[s];
		for(i=0;i<7;i++){
			r[0]=md->bd.x [eid][i];       r[1]=md->bd.y [eid][i];
			dx  =md->bd.dx[eid][i];       dy  =md->bd.dy[eid][i];
			infd_HE(&Hz,Er,r,&(md->wd));
			if(i<3){
				Hz   +=md->bd.sb[0].u[s][i];
				Er[0]+=md->bd.sb[0].v[s][i];
				Er[1]+=md->bd.sb[0].w[s][i];
			}
			else {
				HE_bv(&Hsz,Esr,0,md->bd.et[i],s,md);
				Hz+=Hsz;
				Er[0]+=Esr[0];
				Er[1]+=Esr[1];
			}
			ahz2=creal(Hz*conj(Hz));
			aex2=creal(Er[0]*conj(Er[0]));
			aey2=creal(Er[1]*conj(Er[1]));
			tE=0.25*epsd*(aex2-aey2);
			tH=0.25*ahz2;
			Txx=-tH+tE;
			Txy=0.50*epsd*creal(Er[0]*conj(Er[1]));
			Tyx=Txy;
			Tyy=-tH-tE;
			dFx=Txx*dy-Txy*dx;
			dFy=Tyx*dy-Tyy*dx;
			Fr[0]+=md->bd.wk[i]*dFx;
			Fr[1]+=md->bd.wk[i]*dFy;
			*Nz+=md->bd.wk[i]*((r[0]-rc[0])*dFy-(r[1]-rc[1])*dFx);
		}
	}

	Fr[0]*=-1.0;
	Fr[1]*=-1.0;
	*Nz*=-1.0;
}
