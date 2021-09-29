#include "bem2_emf_tm.h"

int main(int argc,char *argv[])
{
  DOMD md;
  double complex Hz,Er[2];
  double Fr[2],Nz,r[2],rc[2];
  int type;

  dat_read(argv[1],&md);      // read datafile output by 'd2te_bv_solver'
  print_data(&md);            // print data
  //print_data_MKSA(&md);     // for print MKSA system of units

  r[0]=0.3;                   // set x-coordinate
  r[1]=0.0;                   // set y-coordinate
  type=0;                     // selsect 3 point Gauss-Legendre 
  HE_t(&Hz,Er,r,type,&md);    // calclation of total field ( add incident field to scattered field )
  printf("Total electromagnetic field at r=( % g,% g )\n",r[0],r[1]);
  printf("type=0 setting ( 3 point Gauss-Legendre )\n");
  printf("Ex = % 15.14e %+15.14e I \n",creal(Er[0]),cimag(Er[0]));
  printf("Ey = % 15.14e %+15.14e I \n",creal(Er[1]),cimag(Er[1]));
  printf("Hz = % 15.14e %+15.14e I \n",creal(Hz),cimag(Hz)); 
  
  type=1;                     // selsect 7 point Gauss-Kronrod 
  HE_t(&Hz,Er,r,type,&md);    // calclation of total field 
  printf("type=1 setting ( 7 point Gauss-Kronrod )\n");
  printf("Ex = % 15.14e %+15.14e I \n",creal(Er[0]),cimag(Er[0]));
  printf("Ey = % 15.14e %+15.14e I \n",creal(Er[1]),cimag(Er[1]));
  printf("Hz = % 15.14e %+15.14e I \n",creal(Hz),cimag(Hz)); 


  rc[0]=0.0;                  // set x-coordinate of rotation center 
  rc[1]=0.0;                  // set y-coordinate of rotation center 
  printf("\nRadiation force and torque\n");
  printf("type=0 ( 3 point Gauss-Legendre ) \n"); 
  type=0;
  force_FN(Fr,&Nz,rc,type,&md);
  printf("(Fx,Fy) = (% 15.14e,% 15.14e)\n",Fr[0],Fr[1]);
  printf("   Nz   =  % 15.14e, center of rotation (% g,% g)\n",Nz,rc[0],rc[1]); 
  printf("type=1 ( 7 point Gauss-Kronrod ) \n");
  type=1;                   
  force_FN(Fr,&Nz,rc,type,&md);
  printf("(Fx,Fy) = (% 15.14e,% 15.14e)\n",Fr[0],Fr[1]);
  printf("   Nz   =  % 15.14e, center of rotation (% g,% g)\n",Nz,rc[0],rc[1]); 

  mfree_domd(&md);
  return 0;
}
