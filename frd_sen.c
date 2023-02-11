/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void frd_sen(double *co,ITG *nk,double *dstn,ITG *inum,ITG *nmethod,
         ITG *kode,char *filab,double *time,ITG *nstate_,
	 ITG *istep,ITG *iinc,ITG *mode,ITG *noddiam,char *description,
	 ITG *mi,ITG *ngraph,ITG *ne,double *cs,char *set,ITG *nset,
	 ITG *istartset,ITG *iendset,ITG *ialset,
	 char *jobnamec,char *output,double *v,ITG *iobject,
	 char *objectset,ITG *ntrans,ITG *inotr,double *trab,
	 ITG *idesvar,char *orname,ITG *icoordinate,ITG *inorm,
         ITG *irand,ITG *ishape,ITG *ifeasd){
	 
     /* stores the results in frd format

     iselect selects which nodes are to be stored:
          iselect=-1 means only those nodes for which inum negative
                     ist, i.e. network nodes
          iselect=+1 means only those nodes for which inum positive
                     ist, i.e. structural nodes
          iselect=0  means both of the above */
  
  FILE *f1;
  
  char m1[4]=" -1",m2[4]=" -2",m3[4]=" -3",fneig[132]="",text[132];

  static ITG icounter=0,nkcoords;

  ITG null,one,i,noutloc,iset,iselect,two,three,nout,noutplus,noutmin,
      mt=mi[1]+1,ioutall=0;

  ITG ncompscalar=1,ifieldscalar[1]={1},icompscalar[1]={0},
      nfieldscalar[2]={1,0};
  ITG ncomptensoro=6,ifieldtensoro[6]={1,1,1,1,1,1},
      icomptensoro[6]={0,1,2,3,5,4},nfieldtensoro[2]={6,0},
      ncomptensord=2,ifieldtensord[4]={1,1},icomptensord[2]={0,1},
      nfieldtensord[2]={2,0};
  ITG ncompvector=3,ifieldvector[3]={1,1,1},icompvector[3]={0,1,2},
      nfieldvector1[2]={3,0};

  double pi,oner;

  strcpy2(fneig,jobnamec,132);
  strcat(fneig,".frd");

  if((f1=fopen(fneig,"ab"))==NULL){
    printf("*ERROR in frd: cannot open frd file for writing...");
    exit(0);
  }

  /* check whether all results have to be stored (also those
     corresponding to inactive nodes or elements) */
  
  if(strcmp1(&output[3],"a")==0) ioutall=1;

  pi=4.*atan(1.);
  null=0;
  one=1;two=2;three=3;
  oner=1.;

  /* determining nout, noutplus and noutmin 
              nout: number of structural and network nodes
              noutplus: number of structural nodes
              noutmin: number of network nodes */

  if(*nmethod!=0){
      nout=0;
      noutplus=0;
      noutmin=0;
      if(ioutall==0){
	  for(i=0;i<*nk;i++){
	      if(inum[i]==0) continue;
	      nout++;
	      if(inum[i]>0) noutplus++;
	      if(inum[i]<0) noutmin++;   
	  }
      }else{
	  for(i=0;i<*nk;i++){
	      nout++;
	      if(inum[i]>0) noutplus++;
	      if(inum[i]<0) noutmin++;   
	  }
      }
  }else{
      nout=*nk;
  }  

  nkcoords=*nk;
  iselect=1;

  if(*inorm==1){

      /* storing the normals to the structure */
      
      frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc); 
      
      fprintf(f1," -4  NORM        4    1\n");
      fprintf(f1," -5  NORMX       1    2	 1    0\n");
      fprintf(f1," -5  NORMY       1    2	 2    0\n");
      fprintf(f1," -5  NORMZ       1    2	 3    0\n");
      fprintf(f1," -5  ALL         1    2	 0    0    1ALL\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompvector,ifieldvector,icompvector,
                nfieldvector1,&iselect,m2,f1,output,m3);

  }else if(*ishape==1){

      /* storing the shape updates to the structure */
      
      frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc); 
      
      fprintf(f1," -4  SEVE        4    1\n");
      fprintf(f1," -5  SEVEX       1    2	 1    0\n");
      fprintf(f1," -5  SEVEY       1    2	 2    0\n");
      fprintf(f1," -5  SEVEZ       1    2	 3    0\n");
      fprintf(f1," -5  ALL         1    2	 0    0    1ALL\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompvector,ifieldvector,icompvector,
                nfieldvector1,&iselect,m2,f1,output,m3);

  }else if(*irand==1){

      /* storing the normals to the structure */
      
      iselect=0;
    
      frdset(&filab[87],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	        &noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  RANDSCA     1    1\n");
      fprintf(f1," -5  RAND        1    1    0    0\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
                nfieldscalar,&iselect,m2,f1,output,m3);

  }else if(*irand==2){

      /* storing the random vectors in the structure */

      frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc); 
      
      fprintf(f1," -4  RANDVEC     4    1\n");
      fprintf(f1," -5  RANDX       1    2	 1    0\n");
      fprintf(f1," -5  RANDY       1    2	 2    0\n");
      fprintf(f1," -5  RANDZ       1    2	 3    0\n");
      fprintf(f1," -5  ALL         1    2	 0    0    1ALL\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompvector,ifieldvector,icompvector,
                nfieldvector1,&iselect,m2,f1,output,m3);
      
  }else if(*irand==3){

      /* storing the normals to the structure */
      
      iselect=0;
    
      frdset(&filab[87],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	        &noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  MEANSCA     1    1\n");
      fprintf(f1," -5  MEAN        1    1    0    0\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
                nfieldscalar,&iselect,m2,f1,output,m3);

  }else if(*irand==4){

      /* storing the random vectors in the structure */

      frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc); 
      
      fprintf(f1," -4  MEANVEC     4    1\n");
      fprintf(f1," -5  MEANX       1    2	 1    0\n");
      fprintf(f1," -5  MEANY       1    2	 2    0\n");
      fprintf(f1," -5  MEANZ       1    2	 3    0\n");
      fprintf(f1," -5  ALL         1    2	 0    0    1ALL\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompvector,ifieldvector,icompvector,
                nfieldvector1,&iselect,m2,f1,output,m3);
      
  }else if(*irand==5){

      /* storing the normals to the structure */
      
      iselect=0;
    
      frdset(&filab[87],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	        &noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  LOCREL      1    1\n");
      fprintf(f1," -5  LOCREL      1    1    0    0\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
                nfieldscalar,&iselect,m2,f1,output,m3);

  }else if(*ifeasd==1){

      /* storing the sensitivities of the gradient projection method 
         w.r.t. to the control field s */
      
      frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc); 
      
      fprintf(f1," -4  DFDS        3    1\n");
      fprintf(f1," -5  DFDSCONS    1    1	 1    0\n");
      fprintf(f1," -5  DFDSOBJE    1    1	 2    0\n");
      fprintf(f1," -5  DFDSFEAS    1    1	 3    0\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompvector,ifieldvector,icompvector,
                nfieldvector1,&iselect,m2,f1,output,m3);   

  }else if(*ifeasd==2){

      /* storing the sensitivities of the gradient projection method 
         w.r.t. to the control field s */
      
      frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc); 
      
      fprintf(f1," -4  DFDS        3    1\n");
      fprintf(f1," -5  DFDSCORR    1    1	 1    0\n");
      fprintf(f1," -5  DFDSPROJ    1    1	 2    0\n");
      fprintf(f1," -5  DFDSFEAS    1    1	 3    0\n");
      
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompvector,ifieldvector,icompvector,
                nfieldvector1,&iselect,m2,f1,output,m3);   

  }else if(*ifeasd==3){

      /* storing the feasible direction snsitivities of the FE nodes */
      
      frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc); 
      
      fprintf(f1," -4  FEASDIR     2    1\n");
      fprintf(f1," -5  DFDN        1    1    1    0\n");
      fprintf(f1," -5  DFDNFIL     1    1    2    0\n");

      frdselect(&v[2**nk**iobject],v,&iset,&nkcoords,inum,m1,istartset,
	    iendset,ialset,ngraph,&ncomptensord,ifieldtensord,icomptensord,
	    nfieldtensord,&iselect,m2,f1,output,m3);
      
  }else if(*icoordinate!=1){

      /* storing the orientation sensitivities in the nodes */
  
      if((strcmp1(&objectset[(*iobject-1)*405],"ALL-DISP")==0)||
	 (strcmp1(&objectset[(*iobject-1)*405],"EIGENFREQUENCY")==0)||
	 (strcmp1(&objectset[(*iobject-1)*405],"GREEN")==0)){
	  
	  frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
	  
	  frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		    &noutloc,description,kode,nmethod,f1,output,istep,iinc); 
	  
	  strcpy1(&text[0]," -4              4    1",23);
	  strcpy1(&text[5],"D",1);
	  strcpy1(&text[6],&orname[80*(*idesvar/3)],5);
	  if(*idesvar-(*idesvar/3)*3==0){
	      strcpy1(&text[11],"Rx",2);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dD1dRx      1    2    1    0\n");
	      fprintf(f1," -5  dD2dRx      1    2    1    0\n");
	      fprintf(f1," -5  dD3dRx      1    2    1    0\n");
	  }else if(*idesvar-(*idesvar/3)*3==1){
	      strcpy1(&text[11],"Ry",2);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dD1dRy      1    2    1    0\n");
	      fprintf(f1," -5  dD2dRy      1    2    1    0\n");
	      fprintf(f1," -5  dD3dRy      1    2    1    0\n");
	  }else{
	      strcpy1(&text[11],"Rz",2);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dD1dRz      1    2    1    0\n");
	      fprintf(f1," -5  dD2dRz      1    2    1    0\n");
	      fprintf(f1," -5  dD3dRz      1    2    1    0\n");
	  }
	  fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");
	  
	  frdvector(v,&iset,ntrans,&filab[4002],&nkcoords,inum,m1,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);

	  if(*noddiam>=0){
	  
	    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		      &noutloc,description,kode,nmethod,f1,output,istep,iinc); 
	    
	    strcpy1(&text[0]," -4              4    1",23);
	    strcpy1(&text[5],"D",1);
	    strcpy1(&text[6],&orname[80*(*idesvar/3)],5);
	    if(*idesvar-(*idesvar/3)*3==0){
	      strcpy1(&text[10],"RxI",3);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dD1dRx      1    2    1    0\n");
	      fprintf(f1," -5  dD2dRx      1    2    1    0\n");
	      fprintf(f1," -5  dD3dRx      1    2    1    0\n");
	    }else if(*idesvar-(*idesvar/3)*3==1){
	      strcpy1(&text[10],"RyI",3);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dD1dRy      1    2    1    0\n");
	      fprintf(f1," -5  dD2dRy      1    2    1    0\n");
	      fprintf(f1," -5  dD3dRy      1    2    1    0\n");
	    }else{
	      strcpy1(&text[10],"RzI",3);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dD1dRz      1    2    1    0\n");
	      fprintf(f1," -5  dD2dRz      1    2    1    0\n");
	      fprintf(f1," -5  dD3dRz      1    2    1    0\n");
	    }
	    fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");
	    
	    frdvector(&v[*nk*mt],&iset,ntrans,&filab[4002],&nkcoords,inum,m1,inotr,
		    trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);
	  }
	  
      }else if(strcmp1(&objectset[(*iobject-1)*405],"MISESSTRESS")==0){
	  
	  frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
		 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
		 ngraph);
	  
	  frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		    &noutloc,description,kode,nmethod,f1,output,istep,iinc); 
	  
	  strcpy1(&text[0]," -4              6    1",23);
	  strcpy1(&text[5],"S",1);
	  strcpy1(&text[6],&orname[80*(*idesvar/3)],5);
	  if(*idesvar-(*idesvar/3)*3==0){
	      strcpy1(&text[11],"Rx",2);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dSXXdRx     1    4    1    1\n");
	      fprintf(f1," -5  dSYYdRx     1    4    2    2\n");
	      fprintf(f1," -5  dSZZdRx     1    4    3    3\n");
	      fprintf(f1," -5  dSXYdRx     1    4    1    2\n");
	      fprintf(f1," -5  dSYZdRx     1    4    2    3\n");
	      fprintf(f1," -5  dSZXdRx     1    4    3    1\n");
	  }else if(*idesvar-(*idesvar/3)*3==1){
	      strcpy1(&text[11],"Ry",2);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dSXXdRy     1    4    1    1\n");
	      fprintf(f1," -5  dSYYdRy     1    4    2    2\n");
	      fprintf(f1," -5  dSZZdRy     1    4    3    3\n");
	      fprintf(f1," -5  dSXYdRy     1    4    1    2\n");
	      fprintf(f1," -5  dSYZdRy     1    4    2    3\n");
	      fprintf(f1," -5  dSZXdRy     1    4    3    1\n");
	  }else{
	      strcpy1(&text[11],"Rz",2);
	      text[23]='\0';
	      fprintf(f1,"%s\n",text);
	      fprintf(f1," -5  dSXXdRz     1    4    1    1\n");
	      fprintf(f1," -5  dSYYdRz     1    4    2    2\n");
	      fprintf(f1," -5  dSZZdRz     1    4    3    3\n");
	      fprintf(f1," -5  dSXYdRz     1    4    1    2\n");
	      fprintf(f1," -5  dSYZdRz     1    4    2    3\n");
	      fprintf(f1," -5  dSZXdRz     1    4    3    1\n");
	  }
	  
	  frdselect(dstn,dstn,&iset,&nkcoords,inum,m1,istartset,iendset,
		    ialset,ngraph,&ncomptensoro,ifieldtensoro,icomptensoro,
		    nfieldtensoro,&iselect,m2,f1,output,m3);	  
      } 
      
  }else{

      /* storing the coordinate sensitivities in the nodes */
  
      frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc); 
      
      if(strcmp1(&objectset[*iobject*405],"STRAINENERGY")==0){
	  fprintf(f1," -4  SENENER     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"MASS")==0){
	  fprintf(f1," -4  SENMASS     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"ALL-DISP")==0){
	  fprintf(f1," -4  SENDISA     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"X-DISP")==0){
	  fprintf(f1," -4  SENDISX     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"Y-DISP")==0){
	  fprintf(f1," -4  SENDISY     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"Z-DISP")==0){
	  fprintf(f1," -4  SENDISZ     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"MISESSTRESS")==0){
	  fprintf(f1," -4  SENSTRE     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"PS1STRESS")==0){
	  fprintf(f1," -4  SENSTRE     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"PS3STRESS")==0){
	  fprintf(f1," -4  SENSTRE     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"EQPLASTICSTRAIN")==0){
	  fprintf(f1," -4  SENSPEQ     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"MODALSTRESS")==0){
	  fprintf(f1," -4  SENSTRE     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"EIGENFREQUENCY")==0){
	  fprintf(f1," -4  SENFREQ     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"MAXMEMBERSIZE")==0){
	  fprintf(f1," -4  SENMAXM     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"MINMEMBERSIZE")==0){
	  fprintf(f1," -4  SENMINM     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"MAXGROWTH")==0){
	  fprintf(f1," -4  SENGROW     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"MAXSHRINKAGE")==0){
	  fprintf(f1," -4  SENSHRN     2    1\n");
      }else if(strcmp1(&objectset[*iobject*405],"PACKAGING")==0){
	  fprintf(f1," -4  SENPACK     2    1\n");
      }
      
      fprintf(f1," -5  DFDN        1    1    1    0\n");
      fprintf(f1," -5  DFDNFIL     1    1    2    0\n");

      frdselect(&v[2**nk**iobject],v,&iset,&nkcoords,inum,m1,istartset,
	    iendset,ialset,ngraph,&ncomptensord,ifieldtensord,icomptensord,
	    nfieldtensord,&iselect,m2,f1,output,m3);

  }
  
  fclose(f1);
  return;
  
}
