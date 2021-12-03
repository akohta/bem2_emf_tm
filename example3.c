#include "bem2_emf_tm.h"
#include <sys/stat.h>
#include <errno.h>  
#include <png.h>

typedef struct image_data{
  char dir_name[64];      // directory name to output image
  int scale;              // number for enlarge the output image
  
  int m;                  // sampling number 
  double rang;            // range of sampling

  int ts;                 // time step per cycle
  
  int type;               // type setting for surface integral
  
  double complex *ve,*vh; // electromagnetic field data
  double mhz,mer[2];      // maximum amplitude of each field component
}IMD;

void directory_name(char *src,char *nn);
void make_directory(char *dir_name);

void eh_field_z(IMD *id,DOMD *md);
void output_field(char *pl,IMD *id,DOMD *md);

// color table
png_byte ct1[9][3]={{0x00,0x00,0x90},{0x00,0x0f,0xff},{0x00,0x90,0xff},{0x0f,0xff,0xee},
                    {0xff,0xff,0xff},{0xff,0xee,0x00},{0xff,0x70,0x00},{0xee,0x00,0x00},{0x7f,0x00,0x00}};
/*                    
png_byte ct1[9][3]={{0x00,0x00,0x90},{0x00,0x0f,0xff},{0x00,0x90,0xff},{0x0f,0xff,0xee},
                    {0x90,0xff,0x70},{0xff,0xee,0x00},{0xff,0x70,0x00},{0xee,0x00,0x00},{0x7f,0x00,0x00}};  
*/

int main(int argc,char *argv[])
{
  DOMD md;
  IMD id;

  dat_read(argv[1],&md); // read data file
  print_data(&md);       // print data
  
  directory_name(argv[1],id.dir_name); // remove file-extension from argv[1] and add "_images"
  id.scale=1;                          // number for enlarge the output image
  id.m=200;                            // sampling number 
  id.rang=1.5*md.wd.lambda0;           // range of sampling
  id.ts=40;                            // time step per cycle
  id.type=1;                           // type=1 : 7 point GK
  
  make_directory(id.dir_name);

  id.ve=(double complex *)m_alloc2(id.m*id.m*2,sizeof(double complex),"example3.c, ve");
  id.vh=(double complex *)m_alloc2(id.m*id.m*1,sizeof(double complex),"example3.c, vh");

  // z=0 plane
  eh_field_z(&id,&md);
  output_field("xy",&id,&md);
  printf("Image output is finished.\n");

  free(id.ve);
  free(id.vh);
  
  mfree_domd(&md);
  return 0;
}

void directory_name(char *src,char *nn)
{
  int s1,s2;
  char *sd,fo[64]={},buf[54]={};
  
  s1=strlen(src);
  if(s1>54){
    printf("example3.c, directory_name(), directory name is too long. exit...\n");
    exit(1);
  }
  sprintf(fo,"%s",src);
  sd=strrchr(fo,'.');
  if(sd!=NULL){
    s2=strlen(sd);
    strncpy(buf,src,s1-s2);
    sprintf(fo,"%s_images",buf);
  }
  sprintf(nn,"%s",fo);
}

void make_directory(char *dir_name)
{
  int ret;
  
  ret=mkdir(dir_name,S_IRWXU|S_IRWXG);
  if(ret!=0 && errno!=EEXIST){
    printf("failed to make directory. Exit..");
    exit(1);
  }
}

void eh_field_z(IMD *id,DOMD *md)
{
  double complex hz,er[2];
  double x[2],dr;
  int i,j;
  
  dr=id->rang*2.0/(double)(id->m-1);
  
  id->mhz=0.0;
  id->mer[0]=0.0;
  id->mer[1]=0.0;
   
  // z=0 plane  
  #pragma omp parallel for schedule(dynamic) private(x,j,hz,er) 
  for(i=0;i<id->m;i++){
    x[1]=id->rang-(double)i*dr;
    for(j=0;j<id->m;j++){
      x[0]=-id->rang+(double)j*dr;
      HE_t(&hz,er,x,id->type,md); // total field
      
      #pragma omp critical
      {
        if(cabs(hz)>id->mhz) id->mhz=cabs(hz);
        if(cabs(er[0])>id->mer[0]) id->mer[0]=cabs(er[0]);
        if(cabs(er[1])>id->mer[1]) id->mer[1]=cabs(er[1]);
      }
      
      id->vh[i*id->m+j]=hz;
      id->ve[i*id->m*2+j*2+0]=er[0];
      id->ve[i*id->m*2+j*2+1]=er[1];
    }
  }
}

void output_field(char *pl,IMD *id,DOMD *md)
{
  void output_png(int nt,double complex cet,char *pl,IMD *id);
  void output_color_bar(IMD *id);
  
  FILE *fp;
  char fn[128];
  double dt;
  int n;
  
  dt=md->wd.lambda0/(double)id->ts;
  
  #pragma omp parallel for schedule(dynamic) 
  for(n=0;n<id->ts;n++){
    output_png(n,cexp(-I*md->wd.k0*dt*(double)n),pl,id);
  }

  // print info
  sprintf(fn,"%s/%s_info.txt",id->dir_name,pl);
  fp=fopen(fn,"wt");
  if(fp==NULL){
    printf("Failed to open the %s file. Exit...\n",fn);
    exit(1);
  }
  fprintf(fp,"the range of color bar\n");
  fprintf(fp,"Ex is %8e to %8e\n",-id->mer[0],id->mer[0]);
  fprintf(fp,"Ey is %8e to %8e\n",-id->mer[1],id->mer[1]);
  fprintf(fp,"Hz is %8e to %8e\n",-id->mhz,id->mhz);
  fclose(fp);
  
  // output color bar image
  output_color_bar(id);
}

void output_png(int nt,double complex cet,char *pl,IMD *id)
{
  int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b); // -1 <= x <= 1
  
  FILE *fp[3];
  char fname[256],*sf[3]={"Hz","Ex","Ey"};
  int j,i,sj,si,d,m,scale;
  png_uint_32 width,height;
  png_structp png[3];
  png_infop info[3];
  png_bytepp pd[3];
  png_byte r,g,b;

  m=id->m;
  scale=id->scale;

  width =m*(scale+1);
  height=m*(scale+1);

  for(d=0;d<3;d++){
    png[d] =png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info[d]=png_create_info_struct(png[d]);
    sprintf(fname,"%s/%s_%s_%03d.png",id->dir_name,pl,sf[d],nt);
    fp[d]=fopen(fname,"wb");
    if(fp[d]==NULL){
      printf("Failed to open the %s file. Exit...\n",fname);
      exit(1);
    }
    
    png_init_io(png[d],fp[d]);
    png_set_IHDR(png[d],info[d],width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
    pd[d]=(png_bytepp)png_malloc(png[d],sizeof(png_bytep)*height);
    png_set_rows(png[d],info[d],pd[d]);
    
    for(j=0;j<height;j++) pd[d][j]=(png_bytep)png_malloc(png[d],sizeof(png_byte)*width*3);
  }

  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      color_rgb(creal(cet*id->vh[i*m+j])/id->mhz,&r,&g,&b);
      for(si=0;si<=scale;si++){
        for(sj=0;sj<=scale;sj++){
          pd[0][i*(scale+1)+si][(j*(scale+1)+sj)*3+0]=r;
          pd[0][i*(scale+1)+si][(j*(scale+1)+sj)*3+1]=g;
          pd[0][i*(scale+1)+si][(j*(scale+1)+sj)*3+2]=b;
        }
      }

      color_rgb(creal(cet*id->ve[i*m*2+j*2+0])/id->mer[0],&r,&g,&b);
      for(si=0;si<=scale;si++){
        for(sj=0;sj<=scale;sj++){
          pd[1][i*(scale+1)+si][(j*(scale+1)+sj)*3+0]=r;
          pd[1][i*(scale+1)+si][(j*(scale+1)+sj)*3+1]=g;
          pd[1][i*(scale+1)+si][(j*(scale+1)+sj)*3+2]=b;
        }
      }
      
      color_rgb(creal(cet*id->ve[i*m*2+j*2+1])/id->mer[1],&r,&g,&b);
      for(si=0;si<=scale;si++){
        for(sj=0;sj<=scale;sj++){
          pd[2][i*(scale+1)+si][(j*(scale+1)+sj)*3+0]=r;
          pd[2][i*(scale+1)+si][(j*(scale+1)+sj)*3+1]=g;
          pd[2][i*(scale+1)+si][(j*(scale+1)+sj)*3+2]=b;
        }
      }      
    }
  }
  
  for(d=0;d<3;d++){
    png_write_png(png[d],info[d],PNG_TRANSFORM_IDENTITY,NULL);
    
    for(j=0;j<height;j++) png_free(png[d],pd[d][j]);
    png_free(png[d],pd[d]);
    
    fclose(fp[d]);
  }
}

void output_color_bar(IMD *id)
{
  int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b); // -1 <= x <= 1
  
  FILE *fp;
  char fname[128];
  int j,i;
  
  png_uint_32 width,height;
  png_structp png;
  png_infop info;
  png_bytepp pdata;
  png_byte r,g,b;

  sprintf(fname,"%s/color_bar.png",id->dir_name);

  height=id->m*(id->scale+1);
  width=height/16;
  
  png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info= png_create_info_struct(png);
  
  fp=fopen(fname,"wb");
  if(fp==NULL){
    printf("Failed to open the %s file. Exit...\n",fname);
    exit(1);
  }
  
  png_init_io(png, fp);
  png_set_IHDR(png,info,width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
  pdata=(png_bytepp)png_malloc(png, sizeof(png_bytep)*height);
  png_set_rows(png,info,pdata);

  for(j=0;j<height;j++){
    pdata[j]=(png_bytep)png_malloc(png,sizeof(png_byte)*width*3);
  }
  
  for(i=0;i<height;i++){
    color_rgb(1.0-(2.0/(double)height)*(double)i,&r,&g,&b);
    for(j=0;j<width;j++){
      pdata[i][j*3+0]=r;
      pdata[i][j*3+1]=g;
      pdata[i][j*3+2]=b;
    }
  }
  
  png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);
  
  for(j=0;j<height;j++){
    png_free(png,pdata[j]);
  }
  png_free(png,pdata);
  fclose(fp);
}

int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b) // -1 <= x <= 1
{
  double i_nc,dr,dg,db;
  unsigned int i,n,nc,nd;

  if(x<-1.0 || x>1.0){
    *r=0x00;    *g=0x00;    *b=0x00;
    return -1;
  }
  
  n=(unsigned int)floor(pow(2,23)*(x+1.0));
  nc=(unsigned int)pow(2,21);
  i_nc=1.0/(double)nc;
  
  if(n<nc*1)      i=1;
  else if(n<nc*2) i=2;
  else if(n<nc*3) i=3;
  else if(n<nc*4) i=4;
  else if(n<nc*5) i=5;
  else if(n<nc*6) i=6;
  else if(n<nc*7) i=7;
  else if(n<nc*8) i=8;
  else {
    *r=ct1[8][0];    *g=ct1[8][1];    *b=ct1[8][2];
    return 0;
  }
    
  nd=n-nc*(i-1);
  dr=(double)(ct1[i][0]-ct1[i-1][0])*i_nc;
  dg=(double)(ct1[i][1]-ct1[i-1][1])*i_nc;
  db=(double)(ct1[i][2]-ct1[i-1][2])*i_nc;
  *r=(png_byte)floor((double)ct1[i-1][0]+dr*(double)nd);
  *g=(png_byte)floor((double)ct1[i-1][1]+dg*(double)nd);
  *b=(png_byte)floor((double)ct1[i-1][2]+db*(double)nd);
  
  return 0;  
}
