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

void frd(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne0,
	 double *v,double *stn,ITG *inum,ITG *nmethod,ITG *kode,
	 char *filab,double *een,double *t1,double *fn,double *time,
	 double *epn,ITG *ielmat,char *matname,double *enern,
	 double *xstaten,ITG *nstate_,ITG *istep,ITG *iinc,
	 ITG *ithermal,double *qfn,ITG *mode,ITG *noddiam,
	 double *trab,ITG *inotr,ITG *ntrans,double *orab,
	 ITG *ielorien,ITG *norien,char *description,ITG *ipneigh,
	 ITG *neigh,ITG *mi,double *stx,double *vr,double *vi,
	 double *stnr,double *stni,double *vmax,double *stnmax,
	 ITG *ngraph,double *veold,double *ener,ITG *ne,double *cs,
	 char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
	 double *eenmax,double *fnr,double *fni,double *emn,
	 double *thicke,char *jobnamec,char *output,double *qfx,
         double *cdn,ITG *mortar,double *cdnr,double *cdni,ITG *nmat,
         ITG *ielprop,double *prop,double *sti){

  /* stores the results in frd format

     iselect selects which nodes are to be stored:
     iselect=-1 means only those nodes for which inum negative
     ist, i.e. network nodes
     iselect=+1 means only those nodes for which inum positive
     ist, i.e. structural nodes
     iselect=0  means both of the above */
  
  FILE *f1;
  
  char c[2]="C",m1[4]=" -1",m2[4]=" -2",m3[4]=" -3",
    p0[6]="    0",p1[6]="    1",p2[6]="    2",p3[6]="    3",p4[6]="    4",
    p5[6]="    5",p6[6]="    6",p7[6]="    7",p8[6]="    8",p9[6]="    9",
    p10[6]="   10",p11[6]="   11",
    p12[6]="   12", fneig[132]="",date[8],clock[10],newdate[21],newclock[9],
    material[59]="                                                          ",
    text[2]=" ";

  static ITG icounter=0,nkcoords,iaxial;

  ITG null,one,i,j,k,indexe,nemax,nlayer,noutloc,iset,iselect,ncomp,nope,
    nodes,ifield[7],nfield[2],icomp[7],ifieldstate[*nstate_],two,three,
    icompstate[*nstate_],ip0=0,ip1=1,ip2=2,ip3=3,ip4=4,ip5=5,ip6=6,ip7=7,
    ip8=8,ip9=9,ip10=10,ip11=11,ip12=12,imat,nelout,ioutall=0,*inumshell=NULL,
    nterms,nout,noutplus,noutmin,mt=mi[1]+1,iflag,numfield,iorienloc;

  ITG ncompscalar=1,ifieldscalar[1]={1},icompscalar[1]={0},
    nfieldscalar[2]={1,0};
  ITG ncompvector=3,ifieldvector[3]={1,1,1},icompvector[3]={0,1,2},
    nfieldvector1[2]={3,0},nfieldvector0[2]={mi[1]+1,0},
    icompvectorlast[3]={3,4,5};
  ITG ncomptensor=6,ifieldtensor[6]={1,1,1,1,1,1},icomptensor[6]={0,1,2,3,5,4},
    nfieldtensor[2]={6,0};
  ITG ncompscalph=2,ifieldscalph[2]={1,2},icompscalph[2]={0,0},
    nfieldscalph[2]={0,0};
  ITG ncompvectph=6,ifieldvectph[6]={1,1,1,2,2,2},icompvectph[6]={1,2,3,1,2,3},
    nfieldvectph[2]={mi[1]+1,mi[1]+1};
  ITG ncomptensph=12,ifieldtensph[12]={1,1,1,1,1,1,2,2,2,2,2,2},
    icomptensph[12]={0,1,2,3,5,4,0,1,2,3,5,4},nfieldtensph[2]={6,6};
      
  int iw;

  float fl;

  double pi,oner,*errn=NULL,*ethn=NULL;

  strcpy2(fneig,jobnamec,132);
  strcat(fneig,".frd");

  if((f1=fopen(fneig,"ab"))==NULL){
    printf(" *ERROR in frd: cannot open frd file for writing...");
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
	if(inum[i]>=0) noutplus++;
	if(inum[i]<=0) noutmin++;
      }
    }
  }else{
    nout=*nk;
  }

  /* first time something is written in the frd-file: store
     computational metadata, the nodal coordinates and the
     topology */

  if((*kode==1)&&((*nmethod!=5)||(*mode!=0))){
    iaxial=0.;

    /* date and time */

    FORTRAN(dattime,(date,clock));

    for(i=0;i<20;i++) newdate[i]=' ';
    newdate[20]='\0';

    strcpy1(newdate,&date[6],2);
    strcpy1(&newdate[2],".",1);
    if(strcmp1(&date[4],"01")==0){
      strcpy1(&newdate[3],"january.",8);
      strcpy1(&newdate[11],&date[0],4);
    }else if(strcmp1(&date[4],"02")==0){
      strcpy1(&newdate[3],"february.",9);
      strcpy1(&newdate[12],&date[0],4);
    }else if(strcmp1(&date[4],"03")==0){
      strcpy1(&newdate[3],"march.",6);
      strcpy1(&newdate[9],&date[0],4);
    }else if(strcmp1(&date[4],"04")==0){
      strcpy1(&newdate[3],"april.",6);
      strcpy1(&newdate[9],&date[0],4);
    }else if(strcmp1(&date[4],"05")==0){
      strcpy1(&newdate[3],"may.",4);
      strcpy1(&newdate[7],&date[0],4);
    }else if(strcmp1(&date[4],"06")==0){
      strcpy1(&newdate[3],"june.",5);
      strcpy1(&newdate[8],&date[0],4);
    }else if(strcmp1(&date[4],"07")==0){
      strcpy1(&newdate[3],"july.",5);
      strcpy1(&newdate[8],&date[0],4);
    }else if(strcmp1(&date[4],"08")==0){
      strcpy1(&newdate[3],"august.",7);
      strcpy1(&newdate[10],&date[0],4);
    }else if(strcmp1(&date[4],"09")==0){
      strcpy1(&newdate[3],"september.",10);
      strcpy1(&newdate[13],&date[0],4);
    }else if(strcmp1(&date[4],"10")==0){
      strcpy1(&newdate[3],"october.",8);
      strcpy1(&newdate[11],&date[0],4);
    }else if(strcmp1(&date[4],"11")==0){
      strcpy1(&newdate[3],"november.",9);
      strcpy1(&newdate[12],&date[0],4);
    }else if(strcmp1(&date[4],"12")==0){
      strcpy1(&newdate[3],"december.",9);
      strcpy1(&newdate[12],&date[0],4);
    }

    strcpy1(newclock,clock,2);
    strcpy1(&newclock[2],":",1);
    strcpy1(&newclock[3],&clock[2],2);
    strcpy1(&newclock[5],":",1);
    strcpy1(&newclock[6],&clock[4],2);
    newclock[8]='\0';

    fprintf(f1,"%5sUUSER                                                              \n",p1);
    fprintf(f1,"%5sUDATE              %20s                            \n",p1,newdate);
    fprintf(f1,"%5sUTIME              %8s                                        \n",p1,newclock);
    fprintf(f1,"%5sUHOST                                                              \n",p1);
    fprintf(f1,"%5sUPGM               CalculiX                                        \n",p1);
    fprintf(f1,"%5sUVERSION           Version DEVELOPMENT                             \n",p1);
    fprintf(f1,"%5sUCOMPILETIME       Tue Mar 14 10:00:37 CET 2023                    \n",p1);
    fprintf(f1,"%5sUDIR                                                               \n",p1);
    fprintf(f1,"%5sUDBN                                                               \n",p1);
    
    for(i=0;i<*nmat;i++){
      strcpy1(material,&matname[80*i],58);
      fprintf(f1,"%5sUMAT%5" ITGFORMAT "%58s\n",p1,i+1,material);
    }

    /* storing the header of the coordinates */

    if(strcmp1(output,"asc")==0){
      fprintf(f1,"%5s%1s                  %12" ITGFORMAT "%38" ITGFORMAT "\n",p2,c,nout,one);
    }else{
      fprintf(f1,"%5s%1s                  %12" ITGFORMAT "%38" ITGFORMAT "\n",p2,c,nout,three);
    }

    /* storing the coordinates themselves */

    if(*nmethod!=0){
      for(i=0;i<*nk;i++){
	if((inum[i]==0)&&(ioutall!=1)) continue;
	if(strcmp1(output,"asc")==0){
	  fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E\n",m1,i+1,(float)co[3*i],
		  (float)co[3*i+1],(float)co[3*i+2]);
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  fwrite(&co[3*i],sizeof(double),1,f1);
	  fwrite(&co[3*i+1],sizeof(double),1,f1);
	  fwrite(&co[3*i+2],sizeof(double),1,f1);
	}
      }
    }else{
      for(i=0;i<*nk;i++){
	if(strcmp1(output,"asc")==0){
	  fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E\n",m1,i+1,(float)co[3*i],
		  (float)co[3*i+1],(float)co[3*i+2]);
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  fwrite(&co[3*i],sizeof(double),1,f1);
	  fwrite(&co[3*i+1],sizeof(double),1,f1);
	  fwrite(&co[3*i+2],sizeof(double),1,f1);
	}
      }
    }

    /* nkcoords is the number of nodes at the time when 
       the nodal coordinates are stored in the frd file */

    nkcoords=*nk;
    if(strcmp1(output,"asc")==0)fprintf(f1,"%3s\n",m3);

    /* determining the number of elements */

      nelout=0;
      for(i=0;i<*ne0;i++){
	if(((ipkon[i]<=-1)&&(ioutall==0))||(ipkon[i]==-1)){
	  continue;

	  /* the following elements are not stored in the .frd file: */

	  /* contact spring element */

	}else if(strcmp1(&lakon[8*i],"ESPRNGC")==0){
	  continue;

	  /* film advection element */

	}else if(strcmp1(&lakon[8*i],"ESPRNGF")==0){
	  continue;

	  /* one-noded spring element */

	}else if((strcmp1(&lakon[8*i],"E")==0)&&
		 (strcmp1(&lakon[8*i+6],"1")==0)){
	  continue;

	  /* coupling element */

	}else if(strcmp1(&lakon[8*i],"DCOUP3D")==0){
	  continue;

	  /* mass element */

	}else if(strcmp1(&lakon[8*i],"MASS")==0){
	  continue;

	  /* user element */

	}else if(strcmp1(&lakon[8*i],"U1")==0){
	  continue;
	}
	nelout++;
      }

    /* storing the topology */

    if(strcmp1(output,"asc")==0){
      fprintf(f1,"%5s%1s                  %12" ITGFORMAT "%38" ITGFORMAT "\n",p3,c,nelout,one);
    }else{
      fprintf(f1,"%5s%1s                  %12" ITGFORMAT "%38" ITGFORMAT "\n",p3,c,nelout,two);
    }
    nemax=*ne0;

    for(i=0;i<*ne0;i++){
      if(ipkon[i]<=-1){

	if(ipkon[i]==-1){
	  continue;
	}else if(ioutall!=0){
	  indexe=-2-ipkon[i];
	}else{
	  continue;
	}

	/* next element types are not stored */

      }else if(strcmp1(&lakon[8*i],"F")==0){
	continue;
      }else if(strcmp1(&lakon[8*i],"ESPRNGC")==0){
	continue;
      }else if(strcmp1(&lakon[8*i],"ESPRNGF")==0){
	continue;
      }else if((strcmp1(&lakon[8*i],"E")==0)&&
               (strcmp1(&lakon[8*i+6],"1")==0)){
	continue;
      }else if(strcmp1(&lakon[8*i],"DCOUP3D")==0){
	continue;
      }else if(strcmp1(&lakon[8*i],"MASS")==0){
	continue;
      }else if(strcmp1(&lakon[8*i],"U1")==0){
	continue;
      }else{
	indexe=ipkon[i];
      }
      imat=ielmat[i*mi[2]];
      if(strcmp1(&lakon[8*i+3],"2")==0){

	/* 20-node brick element */

	if(((strcmp1(&lakon[8*i+6]," ")==0)||
            (strcmp1(&filab[4],"E")==0)||
	    (strcmp1(&lakon[8*i+6],"I")==0))&&
           (strcmp2(&lakon[8*i+6],"LC",2)!=0)){
	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		    m1,i+1,p4,p0,imat,m2);
	    for(j=0;j<10;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    fprintf(f1,"\n%3s",m2);
	    for(j=10;j<12;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    for(j=16;j<19;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    for(j=19;j<20;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    for(j=12;j<16;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    fprintf(f1,"\n");
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip4;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    for(j=0;j<10;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	    for(j=10;j<12;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	    for(j=16;j<19;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	    for(j=19;j<20;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	    for(j=12;j<16;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	  }
	}else if(strcmp2(&lakon[8*i+6],"LC",2)==0){

          /* composite material */

          /* 20-node brick elements */

	  nlayer=0;
	  for(k=0;k<mi[2];k++){
	    if(ielmat[i*mi[2]+k]==0) break;
	    nlayer++;
	  }
	  for(k=0;k<nlayer;k++){
	    nemax++;
	    if(strcmp1(output,"asc")==0){
	      fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		      m1,nemax,p4,p0,imat,m2);
	      for(j=0;j<10;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+28+20*k+j]);
	      fprintf(f1,"\n%3s",m2);
	      for(j=10;j<12;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+28+20*k+j]);
	      for(j=16;j<19;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+28+20*k+j]);
	      for(j=19;j<20;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+28+20*k+j]);
	      for(j=12;j<16;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+28+20*k+j]);
	      fprintf(f1,"\n");
	    }else{
	      iw=(int)nemax;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip4;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	      for(j=0;j<10;j++){iw=(int)kon[indexe+28+20*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	      for(j=10;j<12;j++){iw=(int)kon[indexe+28+20*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	      for(j=16;j<19;j++){iw=(int)kon[indexe+28+20*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	      for(j=19;j<20;j++){iw=(int)kon[indexe+28+20*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	      for(j=12;j<16;j++){iw=(int)kon[indexe+28+20*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	    }
	  }
	}else if(strcmp1(&lakon[8*i+6],"B")==0){

	  /* 3-node beam element */

	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n",m1,i+1,p12,p0,imat);
	    fprintf(f1,"%3s%10" ITGFORMAT "%10" ITGFORMAT "%10" ITGFORMAT "\n",m2,kon[indexe+20],
		    kon[indexe+22],kon[indexe+21]);
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip12;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe+20];fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe+22];fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe+21];fwrite(&iw,sizeof(int),1,f1);
	  }
	}else{

	  /* 8-node 2d element */
    
	  if(strcmp1(&lakon[8*i+6],"A")==0) iaxial=1;
	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
                    m1,i+1,p10,p0,imat,m2);
	    for(j=0;j<8;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+20+j]);
	    fprintf(f1,"\n");
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip10;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    for(j=0;j<8;j++){iw=(int)kon[indexe+20+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	  }
	}
      }else if(strcmp1(&lakon[8*i+3],"8")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){

	  /* 8-node brick element */

	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		    m1,i+1,p1,p0,imat,m2);
	    for(j=0;j<8;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    fprintf(f1,"\n");
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip1;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    for(j=0;j<8;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	  }
	}else if(strcmp1(&lakon[8*i+6],"B")==0){

	  /* 2-node 1d element */

	  if(strcmp1(&lakon[8*i+4],"R")==0){
	    if(strcmp1(output,"asc")==0){
	      fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n",m1,i+1,p11,p0,imat);
	      fprintf(f1,"%3s%10" ITGFORMAT "%10" ITGFORMAT "\n",m2,kon[indexe+8],
		      kon[indexe+9]);
	    }else{
	      iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip11;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)kon[indexe+8];fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)kon[indexe+9];fwrite(&iw,sizeof(int),1,f1);
	    }
	  }else if(strcmp1(&lakon[8*i+4],"I")==0){
	    if(strcmp1(output,"asc")==0){
	      fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n",m1,i+1,p11,p0,imat);
	      fprintf(f1,"%3s%10" ITGFORMAT "%10" ITGFORMAT "\n",m2,kon[indexe+11],
		      kon[indexe+12]);
	    }else{
	      iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip11;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)kon[indexe+11];fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)kon[indexe+12];fwrite(&iw,sizeof(int),1,f1);
	    }
	  }
	}else{

	  /* 4-node 2d element */

	  if(strcmp1(&lakon[8*i+6],"A")==0) iaxial=1;
	  if((strcmp1(&lakon[8*i+4],"R")==0)||
	     (strcmp1(&lakon[8*i+4]," ")==0)){
	    if(strcmp1(output,"asc")==0){
	      fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		      m1,i+1,p9,p0,imat,m2);
	      for(j=0;j<4;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+8+j]);
	      fprintf(f1,"\n");
	    }else{
	      iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip9;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	      for(j=0;j<4;j++){iw=(int)kon[indexe+8+j];
		fwrite(&iw,sizeof(int),1,f1);}
	    }
	  }else if(strcmp1(&lakon[8*i+4],"I")==0){
	    if(strcmp1(output,"asc")==0){
	      fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		      m1,i+1,p9,p0,imat,m2);
	      for(j=0;j<4;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+11+j]);
	      fprintf(f1,"\n");
	    }else{
	      iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip9;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	      for(j=0;j<4;j++){iw=(int)kon[indexe+11+j];
		fwrite(&iw,sizeof(int),1,f1);}
	    }
	  }
	}
      }else if((strcmp1(&lakon[8*i+3],"10")==0)||
               (strcmp1(&lakon[8*i+3],"14")==0)){

	/* 10-node tetrahedral element */

	if(strcmp1(output,"asc")==0){
	  fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		  m1,i+1,p6,p0,imat,m2);
	  for(j=0;j<10;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	  fprintf(f1,"\n");
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip6;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	  for(j=0;j<10;j++){iw=(int)kon[indexe+j];
	    fwrite(&iw,sizeof(int),1,f1);}
	}
      }else if(strcmp1(&lakon[8*i+3],"4")==0){

	/* 4-node tetrahedral element */

	if(strcmp1(output,"asc")==0){
	  fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		  m1,i+1,p3,p0,imat,m2);
	  for(j=0;j<4;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	  fprintf(f1,"\n");
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip3;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	  for(j=0;j<4;j++){iw=(int)kon[indexe+j];
	    fwrite(&iw,sizeof(int),1,f1);}
	}
      }else if(strcmp1(&lakon[8*i+3],"15")==0){
	if(((strcmp1(&lakon[8*i+6]," ")==0)||
	    (strcmp1(&filab[4],"E")==0))&&
           (strcmp2(&lakon[8*i+6],"LC",2)!=0)){

          /* 15-node wedge element */

	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		    m1,i+1,p5,p0,imat,m2);
	    for(j=0;j<9;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    for(j=12;j<13;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    fprintf(f1,"\n%3s",m2);
	    for(j=13;j<15;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    for(j=9;j<12;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    fprintf(f1,"\n");
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip5;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    for(j=0;j<9;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	    for(j=12;j<13;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	    for(j=13;j<15;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	    for(j=9;j<12;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	  }

	}else if(strcmp2(&lakon[8*i+6],"LC",2)==0){

          /* composite material */

          /* 15-node wedge elements */

	  nlayer=0;
	  for(k=0;k<mi[2];k++){
	    if(ielmat[i*mi[2]+k]==0) break;
	    nlayer++;
	  }
	  for(k=0;k<nlayer;k++){
	    nemax++;
	    if(strcmp1(output,"asc")==0){
	      fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		      m1,nemax,p5,p0,imat,m2);
	      for(j=0;j<9;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+21+15*k+j]);
	      for(j=12;j<13;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+21+15*k+j]);
	      fprintf(f1,"\n%3s",m2);
	      for(j=13;j<15;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+21+15*k+j]);
	      for(j=9;j<12;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+21+15*k+j]);
	      fprintf(f1,"\n");
	    }else{
	      iw=(int)nemax;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip5;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	      iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	      for(j=0;j<9;j++){iw=(int)kon[indexe+21+15*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	      for(j=12;j<13;j++){iw=(int)kon[indexe+21+15*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	      for(j=13;j<15;j++){iw=(int)kon[indexe+21+15*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	      for(j=9;j<12;j++){iw=(int)kon[indexe+21+15*k+j];
		fwrite(&iw,sizeof(int),1,f1);}
	    }
	  }
	}else{

	  /* 6-node 2d element */

	  if(strcmp1(&lakon[8*i+6],"A")==0) iaxial=1;
	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		    m1,i+1,p8,p0,imat,m2);
	    for(j=0;j<6;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+15+j]);
	    fprintf(f1,"\n");
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip8;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    for(j=0;j<6;j++){iw=(int)kon[indexe+15+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	  }
	}
      }else if(strcmp1(&lakon[8*i+3],"6")==0){
	if((strcmp1(&lakon[8*i+6]," ")==0)||
	   (strcmp1(&filab[4],"E")==0)){

	  /* 6-node wedge element */

	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		    m1,i+1,p2,p0,imat,m2);
	    for(j=0;j<6;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	    fprintf(f1,"\n");
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip2;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    for(j=0;j<6;j++){iw=(int)kon[indexe+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	  }
	}else{

	  /* 3-node 2d element */

	  if(strcmp1(&lakon[8*i+6],"A")==0) iaxial=1;
	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		    m1,i+1,p7,p0,imat,m2);
	    for(j=0;j<3;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+6+j]);
	    fprintf(f1,"\n");
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip7;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    for(j=0;j<3;j++){iw=(int)kon[indexe+6+j];
	      fwrite(&iw,sizeof(int),1,f1);}
	  }
	}
      }else if((strcmp1(&lakon[8*i],"D")==0)&&(ithermal[1]>1)){
	if(kon[indexe]==0){

	  /* 2-node 1d element (network entry element) */

	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n",m1,i+1,p11,p0,imat);
	    fprintf(f1,"%3s%10" ITGFORMAT "%10" ITGFORMAT "\n",m2,
		    kon[indexe+1],kon[indexe+2]);
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip11;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe+1];fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe+2];fwrite(&iw,sizeof(int),1,f1);
	  }
	}else if(kon[indexe+2]==0){

	  /* 2-node 1d element (network exit element) */

	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n",m1,i+1,p11,p0,imat);
	    fprintf(f1,"%3s%10" ITGFORMAT "%10" ITGFORMAT "\n",m2,kon[indexe],
		    kon[indexe+1]);
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip11;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe];fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe+1];fwrite(&iw,sizeof(int),1,f1);
	  }
	}else{

	  /* 3-node 1d element (genuine network element) */

	  if(strcmp1(output,"asc")==0){
	    fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n",m1,i+1,p12,p0,imat);
	    fprintf(f1,"%3s%10" ITGFORMAT "%10" ITGFORMAT "%10" ITGFORMAT "\n",m2,kon[indexe],
		    kon[indexe+2],kon[indexe+1]);
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip12;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe];fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe+2];fwrite(&iw,sizeof(int),1,f1);
	    iw=(int)kon[indexe+1];fwrite(&iw,sizeof(int),1,f1);
	  }
	}
      }else if((strcmp1(&lakon[8*i],"E")==0)&&
               ((strcmp1(&lakon[8*i+6],"A")==0)||
                (strcmp1(&lakon[8*i+6],"2")==0))){

	/* 2-node 1d element (spring element) */

	if(strcmp1(output,"asc")==0){
	  fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n",m1,i+1,p11,p0,imat);
	  fprintf(f1,"%3s%10" ITGFORMAT "%10" ITGFORMAT "\n",m2,kon[indexe],kon[indexe+1]);
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip11;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)kon[indexe];fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)kon[indexe+1];fwrite(&iw,sizeof(int),1,f1);
	}
	//      }else if(strcmp1(&lakon[8*i],"MASS")==0){

	//  /* MASS: store nothing */

	/* user elements */
	
      }else if(strcmp1(&lakon[8*i],"US45")==0){
	if(strcmp1(output,"asc")==0){
	  fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		  m1,i+1,p9,p0,imat,m2);
	  for(j=0;j<4;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	  fprintf(f1,"\n");
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip9;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	  for(j=0;j<4;j++){iw=(int)kon[indexe+j];
	    fwrite(&iw,sizeof(int),1,f1);}
	}
      
      }else if(strcmp1(&lakon[8*i],"US3")==0){
	if(strcmp1(output,"asc")==0){
	  fprintf(f1,"%3s%10" ITGFORMAT "%5s%5s%5" ITGFORMAT "\n%3s",
		  m1,i+1,p7,p0,imat,m2);
	  for(j=0;j<3;j++)fprintf(f1,"%10" ITGFORMAT "",kon[indexe+j]);
	  fprintf(f1,"\n");
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip7;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)ip0;fwrite(&iw,sizeof(int),1,f1);
	  iw=(int)imat;fwrite(&iw,sizeof(int),1,f1);
	  for(j=0;j<3;j++){iw=(int)kon[indexe+j];
	    fwrite(&iw,sizeof(int),1,f1);}
	}
      }else{

	/* not treated element type: may lead to an inconsistency
	   in the element count and element output, which may 
	   cause a crash while reading a binary output file */

	FORTRAN(writeelem,(&i,lakon));

      }
    }
    if(strcmp1(output,"asc")==0)fprintf(f1,"%3s\n",m3);

    if(*nmethod==0){fclose(f1);return;}
  }

  /*  for cyclic symmetry frequency calculations only results for
      even numbers (= odd modes, numbering starts at 0) are stored */

  if((*nmethod==2)&&(((*mode/2)*2!=*mode)&&(*noddiam>=0))){fclose(f1);return;}

  /* storing the displacements in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(filab,"U ")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(filab,set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      if(mi[1]==3){

	fprintf(f1," -4  DISP        4    1\n");
	fprintf(f1," -5  D1          1    2    1    0\n");
	fprintf(f1," -5  D2          1    2    2    0\n");
	fprintf(f1," -5  D3          1    2    3    0\n");
	fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");
      
	frdvector(v,&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
		  trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);

      }else if((mi[1]>3)&&(mi[1]<7)){

	fprintf(f1," -4  DISP        %1" ITGFORMAT "    1\n",mi[1]+1);
	fprintf(f1," -5  D1          1    2    1    0\n");
	fprintf(f1," -5  D2          1    2    2    0\n");
	fprintf(f1," -5  D3          1    2    3    0\n");
	for(j=4;j<=mi[1];j++){
	  fprintf(f1," -5  D%1" ITGFORMAT "          1    1    0    0\n",j);
	}
	fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

	frdgeneralvector(v,&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
			 trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);
      }else{
	printf(" *WARNING in frd:\n");
	printf("          for output purposes only 4, 5 or 6\n");
	printf("          degrees of freedom are allowed\n");
	printf("          for generalized vectors;\n");
	printf("          actual degrees of freedom = %" ITGFORMAT "\n",mi[1]);
	printf("          output request ist not performed;\n");
      }
    }
  }

  /*     storing the imaginary part of displacements in the nodes
         for the odd modes of cyclic symmetry calculations */

  if(*noddiam>=0){
    if((strcmp1(filab,"U ")==0)&&(*ithermal!=2)){
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  DISPI       4    1\n");
      fprintf(f1," -5  D1          1    2    1    0\n");
      fprintf(f1," -5  D2          1    2    2    0\n");
      fprintf(f1," -5  D3          1    2    3    0\n");
      fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");
      
      frdvector(&v[*nk*mt],&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
		trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);
    }
  }

  /*     storing the imaginary part of displacements in the nodes
         for steady state calculations */

  if((*nmethod==5)&&(*mode==0)){
    if((strcmp1(filab,"U ")==0)&&(*ithermal!=2)){
      iselect=1;

      frdset(filab,set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  DISPI       4    1\n");
      fprintf(f1," -5  D1          1    2    1    0\n");
      fprintf(f1," -5  D2          1    2    2    0\n");
      fprintf(f1," -5  D3          1    2    3    0\n");
      fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");
      
      frdvector(v,&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
		trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);
    }
  }

  /* storing the velocities in the nodes */
  
  if((strcmp1(&filab[1740],"V   ")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[1740],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  VELO        4    1\n");
    fprintf(f1," -5  V1          1    2    1    0\n");
    fprintf(f1," -5  V2          1    2    2    0\n");
    fprintf(f1," -5  V3          1    2    3    0\n");
    fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

    frdvector(veold,&iset,ntrans,&filab[1740],&nkcoords,inum,m1,inotr,
	      trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);
  }

  /* storing the temperatures in the nodes */
  
  if(strcmp1(&filab[87],"NT  ")==0){
    iselect=0;
    
    frdset(&filab[87],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  NDTEMP      1    1\n");
    fprintf(f1," -5  T           1    1    0    0\n");

    if(*ithermal<=1){
      frdselect(t1,t1,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
                nfieldscalar,&iselect,m2,f1,output,m3);
    }else{
      frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
                nfieldvector0,&iselect,m2,f1,output,m3);
    }
  }

  /* storing the electrical potential in the nodes */
  
  if((strcmp1(&filab[3654],"POT ")==0)&&(*ithermal==2)){
    iselect=0;
    
    frdset(&filab[3654],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  ELPOT       1    1\n");
    fprintf(f1," -5  V           1    1    0    0\n");
    
    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
	      nfieldvector0,&iselect,m2,f1,output,m3);
  }

  /* storing the stresses in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[174],"S   ")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[174],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  STRESS      6    1\n");
      fprintf(f1," -5  SXX         1    4    1    1\n");
      fprintf(f1," -5  SYY         1    4    2    2\n");
      fprintf(f1," -5  SZZ         1    4    3    3\n");
      fprintf(f1," -5  SXY         1    4    1    2\n");
      fprintf(f1," -5  SYZ         1    4    2    3\n");
      fprintf(f1," -5  SZX         1    4    3    1\n");
  
      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[4350],"SNEG")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[4350],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  STRNEG      6    1\n");
      fprintf(f1," -5  SXX         1    4    1    1\n");
      fprintf(f1," -5  SYY         1    4    2    2\n");
      fprintf(f1," -5  SZZ         1    4    3    3\n");
      fprintf(f1," -5  SXY         1    4    1    2\n");
      fprintf(f1," -5  SYZ         1    4    2    3\n");
      fprintf(f1," -5  SZX         1    4    3    1\n");

      iflag=-1;
      numfield=6;
      if((norien>0)&&(strcmp1(&filab[4355],"L")==0)){
	iorienloc=1;
      }else{
	iorienloc=0;
      }
      NNEW(inumshell,ITG,*nk);
      FORTRAN(extrapolateshell,(stx,stn,ipkon,inumshell,kon,lakon,&numfield,
				nk,ne,mi,&numfield,orab,ielorien,co,
				&iorienloc,&filab[4354],ielmat,
				thicke,ielprop,prop,&iflag));
      SFREE(inumshell);
      
      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[4437],"SMID")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[4437],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  STRMID      6    1\n");
      fprintf(f1," -5  SXX         1    4    1    1\n");
      fprintf(f1," -5  SYY         1    4    2    2\n");
      fprintf(f1," -5  SZZ         1    4    3    3\n");
      fprintf(f1," -5  SXY         1    4    1    2\n");
      fprintf(f1," -5  SYZ         1    4    2    3\n");
      fprintf(f1," -5  SZX         1    4    3    1\n");

      iflag=0;
      numfield=6;
      if((norien>0)&&(strcmp1(&filab[4442],"L")==0)){
	iorienloc=1;
      }else{
	iorienloc=0;
      }
      NNEW(inumshell,ITG,*nk);
      FORTRAN(extrapolateshell,(stx,stn,ipkon,inumshell,kon,lakon,&numfield,
				nk,ne,mi,&numfield,orab,ielorien,co,
				&iorienloc,&filab[4441],ielmat,
				thicke,ielprop,prop,&iflag));
      SFREE(inumshell);
      
      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[4524],"SPOS")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[4524],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  STRPOS      6    1\n");
      fprintf(f1," -5  SXX         1    4    1    1\n");
      fprintf(f1," -5  SYY         1    4    2    2\n");
      fprintf(f1," -5  SZZ         1    4    3    3\n");
      fprintf(f1," -5  SXY         1    4    1    2\n");
      fprintf(f1," -5  SYZ         1    4    2    3\n");
      fprintf(f1," -5  SZX         1    4    3    1\n");

      iflag=1;
      numfield=6;
      if((norien>0)&&(strcmp1(&filab[4529],"L")==0)){
	iorienloc=1;
      }else{
	iorienloc=0;
      }
      NNEW(inumshell,ITG,*nk);
      FORTRAN(extrapolateshell,(stx,stn,ipkon,inumshell,kon,lakon,&numfield,
				nk,ne,mi,&numfield,orab,ielorien,co,
				&iorienloc,&filab[4528],ielmat,
				thicke,ielprop,prop,&iflag));
      SFREE(inumshell);
      
      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }

  /* storing the imaginary part of the stresses in the nodes
     for the odd modes of cyclic symmetry calculations */
  
  if(*noddiam>=0){
    if((strcmp1(&filab[174],"S   ")==0)&&(*ithermal!=2)){
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
      
      fprintf(f1," -4  STRESSI     6    1\n");
      fprintf(f1," -5  SXX         1    4    1    1\n");
      fprintf(f1," -5  SYY         1    4    2    2\n");
      fprintf(f1," -5  SZZ         1    4    3    3\n");
      fprintf(f1," -5  SXY         1    4    1    2\n");
      fprintf(f1," -5  SYZ         1    4    2    3\n");
      fprintf(f1," -5  SZX         1    4    3    1\n");
      
      frdselect(&stn[6**nk],stn,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
                nfieldtensor,&iselect,m2,f1,output,m3);
      
    }
  }

  /* storing the imaginary part of the stresses in the nodes
     for steady state calculations */
  
  if((*nmethod==5)&&(*mode==0)){
    if((strcmp1(&filab[174],"S   ")==0)&&(*ithermal!=2)){
      iselect=1;

      frdset(&filab[174],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
      
      fprintf(f1," -4  STRESSI     6    1\n");
      fprintf(f1," -5  SXX         1    4    1    1\n");
      fprintf(f1," -5  SYY         1    4    2    2\n");
      fprintf(f1," -5  SZZ         1    4    3    3\n");
      fprintf(f1," -5  SXY         1    4    1    2\n");
      fprintf(f1," -5  SYZ         1    4    2    3\n");
      fprintf(f1," -5  SZX         1    4    3    1\n");
      
      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
                nfieldtensor,&iselect,m2,f1,output,m3);
      
    }
  }

  /* storing the electromagnetic field E in the nodes */
  
  if((strcmp1(&filab[3741],"EMFE")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[3741],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  EMFE        4    1\n");
    fprintf(f1," -5  E1          1    2    1    0\n");
    fprintf(f1," -5  E2          1    2    2    0\n");
    fprintf(f1," -5  E3          1    2    3    0\n");
    fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

    frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompvector,ifieldvector,icompvector,
	      nfieldtensor,&iselect,m2,f1,output,m3);

    if(*nmethod==2){
      iselect=1;
    
      frdset(&filab[3741],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  EMFEI       4    1\n");
      fprintf(f1," -5  E1          1    2    1    0\n");
      fprintf(f1," -5  E2          1    2    2    0\n");
      fprintf(f1," -5  E3          1    2    3    0\n");
      fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

      frdselect(&stn[6**nk],stn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncompvector,ifieldvector,icompvector,
		nfieldtensor,&iselect,m2,f1,output,m3);
    }    
  }

  /* storing the electromagnetic field B in the nodes */
  
  if((strcmp1(&filab[3828],"EMFB")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[3828],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  EMFB        4    1\n");
    fprintf(f1," -5  B1          1    2    1    0\n");
    fprintf(f1," -5  B2          1    2    2    0\n");
    fprintf(f1," -5  B3          1    2    3    0\n");
    fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

    frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompvector,ifieldvector,icompvectorlast,
	      nfieldtensor,&iselect,m2,f1,output,m3);

    if(*nmethod==2){
      iselect=1;
    
      frdset(&filab[3828],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  EMFBI       4    1\n");
      fprintf(f1," -5  B1          1    2    1    0\n");
      fprintf(f1," -5  B2          1    2    2    0\n");
      fprintf(f1," -5  B3          1    2    3    0\n");
      fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

      frdselect(&stn[6**nk],stn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncompvector,ifieldvector,icompvectorlast,
		nfieldtensor,&iselect,m2,f1,output,m3);
    }
    
  }

  /* storing the total strains in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[261],"E   ")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[261],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  TOSTRAIN    6    1\n");
      fprintf(f1," -5  EXX         1    4    1    1\n");
      fprintf(f1," -5  EYY         1    4    2    2\n");
      fprintf(f1," -5  EZZ         1    4    3    3\n");
      fprintf(f1," -5  EXY         1    4    1    2\n");
      fprintf(f1," -5  EYZ         1    4    2    3\n");
      fprintf(f1," -5  EZX         1    4    3    1\n");
  
      frdselect(een,een,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }

  /* storing the imaginary part of the total strains in the nodes
     for the odd modes of cyclic symmetry calculations */
  
  if(*noddiam>=0){
    if((strcmp1(&filab[261],"E   ")==0)&&(*ithermal!=2)){
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
      
      fprintf(f1," -4  TOSTRAII    6    1\n");
      fprintf(f1," -5  EXX         1    4    1    1\n");
      fprintf(f1," -5  EYY         1    4    2    2\n");
      fprintf(f1," -5  EZZ         1    4    3    3\n");
      fprintf(f1," -5  EXY         1    4    1    2\n");
      fprintf(f1," -5  EYZ         1    4    2    3\n");
      fprintf(f1," -5  EZX         1    4    3    1\n");
      
      frdselect(&een[6**nk],een,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
                nfieldtensor,&iselect,m2,f1,output,m3);
      
    }
  }

  /* storing the imaginary part of the total strains in the nodes
     for steady state calculations */
  
  if((*nmethod==5)&&(*mode==0)){
    if((strcmp1(&filab[261],"E   ")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[261],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  TOSTRAII    6    1\n");
      fprintf(f1," -5  EXX         1    4    1    1\n");
      fprintf(f1," -5  EYY         1    4    2    2\n");
      fprintf(f1," -5  EZZ         1    4    3    3\n");
      fprintf(f1," -5  EXY         1    4    1    2\n");
      fprintf(f1," -5  EYZ         1    4    2    3\n");
      fprintf(f1," -5  EZX         1    4    3    1\n");
  
      frdselect(een,een,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }

  /* storing the mechanical strains in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[2697],"ME  ")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[2697],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  MESTRAIN    6    1\n");
      fprintf(f1," -5  MEXX        1    4    1    1\n");
      fprintf(f1," -5  MEYY        1    4    2    2\n");
      fprintf(f1," -5  MEZZ        1    4    3    3\n");
      fprintf(f1," -5  MEXY        1    4    1    2\n");
      fprintf(f1," -5  MEYZ        1    4    2    3\n");
      fprintf(f1," -5  MEZX        1    4    3    1\n");
  
  
      frdselect(emn,emn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }

  /* storing the imaginary part of the mechanical strains in the nodes
     for the odd modes of cyclic symmetry calculations */
  
  if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
    if((strcmp1(&filab[2697],"ME  ")==0)&&(*ithermal!=2)){
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
      
      fprintf(f1," -4  MESTRAII    6    1\n");
      fprintf(f1," -5  MEXX        1    4    1    1\n");
      fprintf(f1," -5  MEYY        1    4    2    2\n");
      fprintf(f1," -5  MEZZ        1    4    3    3\n");
      fprintf(f1," -5  MEXY        1    4    1    2\n");
      fprintf(f1," -5  MEYZ        1    4    2    3\n");
      fprintf(f1," -5  MEZX        1    4    3    1\n");
      
      frdselect(&emn[6**nk],een,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
                nfieldtensor,&iselect,m2,f1,output,m3);
      
    }
  }

  /* storing the thermal strains in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[4698],"THE ")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[4698],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  THSTRAIN    6    1\n");
      fprintf(f1," -5  THXX        1    4    1    1\n");
      fprintf(f1," -5  THYY        1    4    2    2\n");
      fprintf(f1," -5  THZZ        1    4    3    3\n");
      fprintf(f1," -5  THXY        1    4    1    2\n");
      fprintf(f1," -5  THYZ        1    4    2    3\n");
      fprintf(f1," -5  THZX        1    4    3    1\n");
  
      NNEW(ethn,double,6**nk);
      for(i=0;i<6**nk;i++){
	ethn[i]=een[i]-emn[i];
      }
      
      frdselect(ethn,ethn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);

      SFREE(ethn);
      
    }
  }

  /* storing the forces in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[348],"RF  ")==0)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[348],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      if(mi[1]==3){
  
	fprintf(f1," -4  FORC        4    1\n");
	fprintf(f1," -5  F1          1    2    1    0\n");
	fprintf(f1," -5  F2          1    2    2    0\n");
	fprintf(f1," -5  F3          1    2    3    0\n");
	fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");
  
	if((iaxial==1)&&(strcmp1(&filab[352],"I")==0)){for(i=0;i<*nk;i++){fn[1+i*mt]*=180.;fn[2+i*mt]*=180.;fn[3+i*mt]*=180.;}}
	frdvector(fn,&iset,ntrans,&filab[348],&nkcoords,inum,m1,inotr,
		  trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);
	if((iaxial==1)&&(strcmp1(&filab[352],"I")==0)){for(i=0;i<*nk;i++){fn[1+i*mt]/=180.;fn[2+i*mt]/=180.;fn[3+i*mt]/=180.;}}

      }else if((mi[1]>3)&&(mi[1]<7)){

	fprintf(f1," -4  FORC        %1" ITGFORMAT "    1\n",mi[1]+1);
	fprintf(f1," -5  F1          1    2    1    0\n");
	fprintf(f1," -5  F2          1    2    2    0\n");
	fprintf(f1," -5  F3          1    2    3    0\n");
	for(j=4;j<=mi[1];j++){
	  fprintf(f1," -5  F%1" ITGFORMAT "          1    1    0    0\n",j);
	}
	fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

	frdgeneralvector(fn,&iset,ntrans,&filab[348],&nkcoords,inum,m1,inotr,
			 trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);
      }else{
	printf(" *WARNING in frd:\n");
	printf("          for output purposes only 4, 5 or 6\n");
	printf("          degrees of freedom are allowed\n");
	printf("          for generalized vectors;\n");
	printf("          actual degrees of freedom = %" ITGFORMAT "\n",mi[1]);
	printf("          output request ist not performed;\n");
      }
    }
  }

  /*     storing the imaginary part of the forces in the nodes
         for the odd modes of cyclic symmetry calculations */

  if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
    if((strcmp1(&filab[348],"RF  ")==0)&&(*ithermal!=2)){
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  FORCI       4    1\n");
      fprintf(f1," -5  F1          1    2    1    0\n");
      fprintf(f1," -5  F2          1    2    2    0\n");
      fprintf(f1," -5  F3          1    2    3    0\n");
      fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");
      
      frdvector(&fn[*nk*mt],&iset,ntrans,filab,&nkcoords,inum,m1,inotr,
		trab,co,istartset,iendset,ialset,mi,ngraph,f1,output,m3);
    }
  }

  /* storing the equivalent plastic strains in the nodes */
  
  if((strcmp1(&filab[435],"PEEQ")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[435],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  PE          1    1\n");
    fprintf(f1," -5  PE          1    1    0    0\n");

    frdselect(epn,epn,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
	      nfieldscalar,&iselect,m2,f1,output,m3);

  }

  /* storing the energy in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[522],"ENER")==0)&&(*ithermal!=2)){
      iselect=1;
    
      frdset(&filab[522],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
    
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  ENER        1    1\n");
      fprintf(f1," -5  ENER        1    1    0    0\n");

      frdselect(enern,enern,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
                nfieldscalar,&iselect,m2,f1,output,m3);

    }
  }
  
  /* storing the contact displacements and stresses at the slave nodes */

  /* node-to-face penalty */
  
  if((strcmp1(&filab[2175],"CONT")==0)&&(*mortar!=1)&&(*ithermal!=2)&&((*nmethod!=2)&&(*nmethod!=13))){
    
    for(i=*ne-1;i>=0;i--){
      if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	break;
    }
    noutloc=*ne-i-1;
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    
    fprintf(f1," -4  CONTACT     6    1\n");
    fprintf(f1," -5  COPEN       1    4    1    1\n");
    fprintf(f1," -5  CSLIP1      1    4    2    2\n");
    fprintf(f1," -5  CSLIP2      1    4    3    3\n");
    fprintf(f1," -5  CPRESS      1    4    1    2\n");
    fprintf(f1," -5  CSHEAR1     1    4    2    3\n");
    fprintf(f1," -5  CSHEAR2     1    4    3    1\n");
    
    for(i=*ne-1;i>=0;i--){
      if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	break;
      strcpy1(text,&lakon[8*i+7],1);
      nope=atoi(text)+1;
      nodes=kon[ipkon[i]+nope-1];
      if(strcmp1(output,"asc")==0){
	fprintf(f1,"%3s%10" ITGFORMAT "",m1,nodes);
	for(j=0;j<6;j++)fprintf(f1,"%12.5E",(float)stx[6*mi[0]*i+j]);
      }else{
	iw=(int)(nodes);fwrite(&iw,sizeof(int),1,f1);
	if(strcmp1(output,"bin")==0){
	  for(j=0;j<6;j++){
	    fl=(float)stx[6*mi[0]*i+j];
	    fwrite(&fl,sizeof(float),1,f1);
	  }
	}else{
	  for(j=0;j<6;j++){
	    fwrite(&stx[6*mi[0]*i+j],sizeof(double),1,f1);
	  }
	}
      }
      if(strcmp1(output,"asc")==0)fprintf(f1,"\n");
    }
    
    if(strcmp1(output,"asc")==0)fprintf(f1,"%3s\n",m3);
  }

  /* face-to-face penalty */

  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[2175],"CONT")==0)&&(*mortar==1)&&(*ithermal!=2)){
      iselect=1;
  
      frdset(&filab[2175],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
      fprintf(f1," -4  CONTACT     6    1\n");
      fprintf(f1," -5  COPEN       1    4    1    1\n");
      fprintf(f1," -5  CSLIP1      1    4    2    2\n");
      fprintf(f1," -5  CSLIP2      1    4    3    3\n");
      fprintf(f1," -5  CPRESS      1    4    1    2\n");
      fprintf(f1," -5  CSHEAR1     1    4    2    3\n");
      fprintf(f1," -5  CSHEAR2     1    4    3    1\n");
  
      frdselect(cdn,cdn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }

  /* storing imaginary part of the differential contact displacements 
     and the contact stresses for the odd modes of cyclic symmetry
     calculations */

  if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
    if((strcmp1(&filab[2175],"CONT")==0)&&(*mortar==1)&&(*ithermal!=2)){
      iselect=1;
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
      fprintf(f1," -4  CONTACTI    6    1\n");
      fprintf(f1," -5  COPEN       1    4    1    1\n");
      fprintf(f1," -5  CSLIP1      1    4    2    2\n");
      fprintf(f1," -5  CSLIP2      1    4    3    3\n");
      fprintf(f1," -5  CPRESS      1    4    1    2\n");
      fprintf(f1," -5  CSHEAR1     1    4    2    3\n");
      fprintf(f1," -5  CSHEAR2     1    4    3    1\n");
  
      frdselect(&cdn[6**nk],cdn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }
  /* storing the contact energy at the slave nodes */
  
  if((strcmp1(&filab[2262],"CELS")==0)&&(*ithermal!=2)){
    
    for(i=*ne-1;i>=0;i--){
      if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	break;
    }
    noutloc=*ne-i-1;
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);
    
    fprintf(f1," -4  CELS        1    1\n");
    fprintf(f1," -5  CELS        1    1    0    0\n");
    
    for(i=*ne-1;i>=0;i--){
      if((strcmp1(&lakon[8*i+1],"S")!=0)||(strcmp1(&lakon[8*i+6],"C")!=0))
	break;
      nope=atoi(&lakon[8*i+7])+1;
      nodes=kon[ipkon[i]+nope-1];
      if(strcmp1(output,"asc")==0){
	fprintf(f1,"%3s%10" ITGFORMAT "%12.5E\n",m1,nodes,(float)ener[2*i*mi[0]]);
      }else{
	iw=(int)(nodes);fwrite(&iw,sizeof(int),1,f1);
	if(strcmp1(output,"bin")==0){
	  fl=(float)ener[i*mi[0]];
	  fwrite(&fl,sizeof(float),1,f1);
	}else{
	  fwrite(&ener[i*mi[0]],sizeof(double),1,f1);
	}
      }
    }
    
    if(strcmp1(output,"asc")==0)fprintf(f1,"%3s\n",m3);
  }
  
  /* storing the internal state variables in the nodes */
  
  if(strcmp1(&filab[609],"SDV ")==0){
    iselect=1;
    
    frdset(&filab[609],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  SDV       %3" ITGFORMAT "    1\n",*nstate_);
    for(j=1;j<=*nstate_;j++){
      fprintf(f1," -5  SDV%-3" ITGFORMAT "      1    1    0    0\n",j);
    }

    for(i=0;i<*nstate_;i++){
      ifieldstate[i]=1;icompstate[i]=i;
    }
    nfield[0]=*nstate_;

    frdselect(xstaten,xstaten,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,nstate_,ifieldstate,icompstate,
	      nfield,&iselect,m2,f1,output,m3);

  }
  
  /* storing the heat flux in the nodes
     the heat flux has been extrapolated from the integration points
     in subroutine extrapolate.f, taking into account whether the 
     results are requested in the global system or in a local system.
     Therefore, subroutine frdvector cannot be used, since it assumes
     the values are stored in the global system */
  
  if((strcmp1(&filab[696],"HFL ")==0)&&(*ithermal>1)){
    iselect=1;
    
    frdset(&filab[696],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  FLUX        4    1\n");
    fprintf(f1," -5  F1          1    2    1    0\n");
    fprintf(f1," -5  F2          1    2    2    0\n");
    fprintf(f1," -5  F3          1    2    3    0\n");
    fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

    frdselect(qfn,qfn,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompvector,ifieldvector,icompvector,
	      nfieldvector1,&iselect,m2,f1,output,m3);

  }
  
  /* storing the electrical current in the nodes
     (cf. heat flux HFL above)  */
  
  if((strcmp1(&filab[3567],"ECD ")==0)&&(*ithermal==2)){
    iselect=1;
    
    frdset(&filab[3567],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  CURR        4    1\n");
    fprintf(f1," -5  j1          1    2    1    0\n");
    fprintf(f1," -5  j2          1    2    2    0\n");
    fprintf(f1," -5  j3          1    2    3    0\n");
    fprintf(f1," -5  ALL         1    2    0    0    1ALL\n");

    frdselect(qfn,qfn,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompvector,ifieldvector,icompvector,
	      nfieldvector1,&iselect,m2,f1,output,m3);

  }
  
  /* storing the heat generation in the nodes */

  if((strcmp1(&filab[783],"RFL ")==0)&&(*ithermal>1)){
    iselect=1;
    
    frdset(&filab[783],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  RFL         1    1\n");
    fprintf(f1," -5  RFL         1    1    0    0\n");

    frdselect(fn,fn,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
	      nfieldvector0,&iselect,m2,f1,output,m3);

  }
  
  /* storing the Zienkiewicz-Zhu improved stresses in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[1044],"ZZS")==0)&&(*ithermal!=2)){
  
      FORTRAN(zienzhu,(co,nk,kon,ipkon,lakon,ne0,stn,ipneigh,neigh,
		       stx,&mi[0]));
  
      iselect=1;
  
      frdset(&filab[1044],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  ZZSTR       6    1\n");
      fprintf(f1," -5  SXX         1    4    1    1\n");
      fprintf(f1," -5  SYY         1    4    2    2\n");
      fprintf(f1," -5  SZZ         1    4    3    3\n");
      fprintf(f1," -5  SXY         1    4    1    2\n");
      fprintf(f1," -5  SYZ         1    4    2    3\n");
      fprintf(f1," -5  SZX         1    4    3    1\n");
  
      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }

  /* storing the imaginary part of the Zienkiewicz-Zhu 
     improved stresses in the nodes
     for the odd modes of cyclic symmetry calculations */
  
  if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
    if((strcmp1(&filab[1044],"ZZS")==0)&&(*ithermal!=2)){

      FORTRAN(zienzhu,(co,nk,kon,ipkon,lakon,ne0,stn,ipneigh,neigh,
		       &stx[6*mi[0]**ne],&mi[0]));
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
      
      fprintf(f1," -4  ZZSTRI      6    1\n");
      fprintf(f1," -5  SXX         1    4    1    1\n");
      fprintf(f1," -5  SYY         1    4    2    2\n");
      fprintf(f1," -5  SZZ         1    4    3    3\n");
      fprintf(f1," -5  SXY         1    4    1    2\n");
      fprintf(f1," -5  SYZ         1    4    2    3\n");
      fprintf(f1," -5  SZX         1    4    3    1\n");
      
      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
                nfieldtensor,&iselect,m2,f1,output,m3);
      
    }
  }
  
  /* storing the error estimator in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2)){

      NNEW(errn,double,6**nk);
  
      nterms=6;
      FORTRAN(errorestimator,(stx,errn,ipkon,kon,lakon,nk,ne,
			      mi,ielmat,&nterms,inum,co,v,&filab[1048],
			      ielprop,prop));
  
      iselect=1;
  
      frdset(&filab[1044],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  ERROR       1    1\n");
      fprintf(f1," -5  STR(%%)      1    1    0    0\n");

      ncomp=1;
      ifield[0]=1;
      icomp[0]=0;
  
      frdselect(errn,errn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomp,ifield,icomp,
		nfieldtensor,&iselect,m2,f1,output,m3);
  
    }
  }

  /* storing the imaginary part of the error estimator in the nodes
     for the odd modes of cyclic symmetry calculations */
  
  if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
    if((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2)){

      nterms=6;
      FORTRAN(errorestimator,(&stx[6*mi[0]**ne],stn,ipkon,kon,lakon,nk,ne,
			      mi,ielmat,&nterms,inum,co,v,&filab[1048],
			      ielprop,prop));
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  ERRORI      1    1\n");
      fprintf(f1," -5  STR(%%)      1    1    0    0\n");
      
      ncomp=1;
      ifield[0]=1;
      icomp[0]=0;

      frdselect(stn,stn,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncomp,ifield,icomp,
                nfieldtensor,&iselect,m2,f1,output,m3);
      
    }
  }
  
  /* storing the thermal error estimator in the nodes */
  
  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[2784],"HER")==0)&&(*ithermal>1)){
  
      nterms=3;
      FORTRAN(errorestimator,(qfx,qfn,ipkon,kon,lakon,nk,ne,
			      mi,ielmat,&nterms,inum,co,v,&filab[2788],
			      ielprop,prop));
  
      iselect=1;
  
      frdset(&filab[2784],set,&iset,istartset,iendset,ialset,
	     inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	     ngraph);
  
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);
  
      fprintf(f1," -4  HERROR      1    1\n");
      fprintf(f1," -5  TEM(%%)      1    1    0    0\n");
  
      ncomp=1;
      ifield[0]=1;
      icomp[0]=0;
  
      frdselect(qfn,qfn,&iset,&nkcoords,inum,m1,istartset,iendset,
		ialset,ngraph,&ncomp,ifield,icomp,
		nfieldvector1,&iselect,m2,f1,output,m3);
  
    }
  }

  /* storing the imaginary part of the thermal error estimator in the nodes
     for the odd modes of cyclic symmetry calculations */
  
  if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
    if((strcmp1(&filab[2784],"HER")==0)&&(*ithermal>1)){

      nterms=3;
      FORTRAN(errorestimator,(&qfx[3*mi[0]**ne],qfn,ipkon,kon,lakon,nk,ne,
			      mi,ielmat,&nterms,inum,co,v,&filab[2788],
			      ielprop,prop));
      
      frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
		&noutloc,description,kode,nmethod,f1,output,istep,iinc);

      fprintf(f1," -4  HERRORI     1    1\n");
      fprintf(f1," -5  TEM(%%)      1    1    0    0\n");
      
      ncomp=1;
      ifield[0]=1;
      icomp[0]=0;

      frdselect(qfn,qfn,&iset,&nkcoords,inum,m1,istartset,iendset,
                ialset,ngraph,&ncomp,ifield,icomp,
                nfieldtensor,&iselect,m2,f1,output,m3);
      
    }
  }

  /* storing the total temperatures in the network nodes */
  
  if((strcmp1(&filab[1131],"TT  ")==0)&&(*ithermal>1)){

    iselect=-1;
    frdset(&filab[1131],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  TOTEMP      1    1\n");
    fprintf(f1," -5  TT          1    1    0    0\n");

    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icompscalar,
	      nfieldvector0,&iselect,m2,f1,output,m3);

  }

  /* storing the mass flow in the network nodes */
  
  if((strcmp1(&filab[1218],"MF  ")==0)&&(*ithermal>1)){

    iselect=-1;
    frdset(&filab[1218],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  MAFLOW      1    1\n");
    fprintf(f1," -5  MF          1    1    0    0\n");

    icomp[0]=1;
    if((iaxial==1)&&(strcmp1(&filab[1222],"I")==0)){for(i=0;i<*nk;i++)v[1+i*mt]*=180.;}
    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
	      nfieldvector0,&iselect,m2,f1,output,m3);
    if((iaxial==1)&&(strcmp1(&filab[1222],"I")==0)){for(i=0;i<*nk;i++)v[1+i*mt]/=180.;}

  }

  /* storing the total pressure in the network nodes */
  
  if((strcmp1(&filab[1305],"PT  ")==0)&&(*ithermal>1)){

    iselect=-1;
    frdset(&filab[1305],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  TOPRES      1    1\n");
    fprintf(f1," -5  PT          1    1    0    0\n");

    icomp[0]=2;
    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
	      nfieldvector0,&iselect,m2,f1,output,m3);

  }

  /* storing the static pressure in the liquid network nodes */
  
  if((strcmp1(&filab[1827],"PS  ")==0)&&(*ithermal>1)){

    iselect=-1;
    frdset(&filab[1827],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  STPRES      1    1\n");
    fprintf(f1," -5  PS          1    1    0    0\n");

    icomp[0]=2;
    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
	      nfieldvector0,&iselect,m2,f1,output,m3);

  }

  /* storing the liquid depth in the channel nodes */
  
  if((strcmp1(&filab[2349],"DEPT")==0)&&(*ithermal>1)){

    iselect=-1;
    frdset(&filab[2349],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  DEPTH       1    1\n");
    fprintf(f1," -5  DEPTH       1    1    0    0\n");

    icomp[0]=2;
    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
	      nfieldvector0,&iselect,m2,f1,output,m3);

  }

  /* storing the critical depth in the channel nodes */
  
  if((strcmp1(&filab[2436],"HCRI")==0)&&(*ithermal>1)){

    iselect=-1;
    frdset(&filab[2436],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  HCRIT       1    1\n");
    fprintf(f1," -5  HCRIT       1    1    0    0\n");

    icomp[0]=3;
    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
	      nfieldvector0,&iselect,m2,f1,output,m3);

  }

  /* storing the static temperature in the network nodes */
  
  if((strcmp1(&filab[1392],"TS  ")==0)&&(*ithermal>1)){

    iselect=-1;
    frdset(&filab[1392],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  STTEMP      1    1\n");
    fprintf(f1," -5  TS          1    1    0    0\n");

    icomp[0]=3;
    frdselect(v,v,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalar,ifieldscalar,icomp,
	      nfieldvector0,&iselect,m2,f1,output,m3);

  }

  /* mesh refinement */
  
  if(strcmp1(&filab[4089],"RM")==0){
    refinemesh(nk,ne,co,ipkon,kon,v,veold,stn,een,emn,epn,enern,
	       qfn,errn,filab,mi,lakon,jobnamec,istartset,iendset,
	       ialset,set,nset,matname,ithermal,output,nmat);
  }

  /* remove auxiliary field for the error estimator at the nodes */  

  if((*nmethod!=5)||(*mode==-1)){
    if((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2)){
      SFREE(errn);
    }
  }

  /*  the remaining lines only apply to frequency calculations
      with cyclic symmetry, complex frequency and steady state calculations */

  if((*nmethod!=2)&&(*nmethod!=13)&&(*nmethod!=5)&&(*nmethod!=6)&&(*nmethod!=7)){fclose(f1);return;}
  if((*nmethod==5)&&(*mode==-1)){fclose(f1);return;}

  /* storing the displacements in the nodes (magnitude, phase) */
  
  if((strcmp1(&filab[870],"PU  ")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[870],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  PDISP       6    1\n");
    fprintf(f1," -5  MAG1        1   12    1    0\n");
    fprintf(f1," -5  MAG2        1   12    2    0\n");
    fprintf(f1," -5  MAG3        1   12    3    0\n");
    fprintf(f1," -5  PHA1        1   12    4    0\n");
    fprintf(f1," -5  PHA2        1   12    5    0\n");
    fprintf(f1," -5  PHA3        1   12    6    0\n");

    frdselect(vr,vi,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompvectph,ifieldvectph,icompvectph,
	      nfieldvectph,&iselect,m2,f1,output,m3);

  }

  /* storing the temperatures in the nodes (magnitude, phase) */
  
  if((strcmp1(&filab[957],"PNT ")==0)&&(*ithermal>1)){
    iselect=1;
    
    frdset(&filab[957],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  PNDTEMP     2    1\n");
    fprintf(f1," -5  MAG1        1    1    1    0\n");
    fprintf(f1," -5  PHA1        1    1    2    0\n");

    frdselect(vr,vi,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompscalph,ifieldscalph,icompscalph,
	      nfieldscalph,&iselect,m2,f1,output,m3);

  }

  /* storing the stresses in the nodes (magnitude, phase) */
  
  if((strcmp1(&filab[1479],"PHS ")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[1479],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  PSTRESS    12    1\n");
    fprintf(f1," -5  MAGXX       1   14    1    1\n");
    fprintf(f1," -5  MAGYY       1   14    2    2\n");
    fprintf(f1," -5  MAGZZ       1   14    3    3\n");
    fprintf(f1," -5  MAGXY       1   14    1    2\n");
    fprintf(f1," -5  MAGYZ       1   14    2    3\n");
    fprintf(f1," -5  MAGZX       1   14    3    1\n");
    fprintf(f1," -5  PHAXX       1   14    1    1\n");
    fprintf(f1," -5  PHAYY       1   14    2    2\n");
    fprintf(f1," -5  PHAZZ       1   14    3    3\n");
    fprintf(f1," -5  PHAXY       1   14    1    2\n");
    fprintf(f1," -5  PHAYZ       1   14    2    3\n");
    fprintf(f1," -5  PHAZX       1   14    3    1\n");

    frdselect(stnr,stni,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncomptensph,ifieldtensph,icomptensph,
	      nfieldtensph,&iselect,m2,f1,output,m3);

  }

  /* storing the differential contact displacements and
     the contact stresses in the nodes (magnitude, phase)
     only for face-to-face penalty contact */
  
  if((strcmp1(&filab[3915],"PCON")==0)&&(*ithermal!=2)&&(*mortar==1)){
    iselect=1;
    
    frdset(&filab[3915],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  PCONTAC    12    1\n");
    fprintf(f1," -5  MAGO        1   14    1    1\n");
    fprintf(f1," -5  MAGSL1      1   14    2    2\n");
    fprintf(f1," -5  MAGSL2      1   14    3    3\n");
    fprintf(f1," -5  MAGP        1   14    1    2\n");
    fprintf(f1," -5  MAGSH1      1   14    2    3\n");
    fprintf(f1," -5  MAGSH2      1   14    3    1\n");
    fprintf(f1," -5  PHAO        1   14    1    1\n");
    fprintf(f1," -5  PHASL1      1   14    2    2\n");
    fprintf(f1," -5  PHASL2      1   14    3    3\n");
    fprintf(f1," -5  PHAP        1   14    1    2\n");
    fprintf(f1," -5  PHASH1      1   14    2    3\n");
    fprintf(f1," -5  PHASH2      1   14    3    1\n");

    frdselect(cdnr,cdni,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncomptensph,ifieldtensph,icomptensph,
	      nfieldtensph,&iselect,m2,f1,output,m3);

  }

  /* storing the forces in the nodes (magnitude, phase) */
  
  if((strcmp1(&filab[2610],"PRF ")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[2610],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  PFORC       6    1\n");
    fprintf(f1," -5  MAG1        1   12    1    0\n");
    fprintf(f1," -5  MAG2        1   12    2    0\n");
    fprintf(f1," -5  MAG3        1   12    3    0\n");
    fprintf(f1," -5  PHA1        1   12    4    0\n");
    fprintf(f1," -5  PHA2        1   12    5    0\n");
    fprintf(f1," -5  PHA3        1   12    6    0\n");

    frdselect(fnr,fni,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncompvectph,ifieldvectph,icompvectph,
	      nfieldvectph,&iselect,m2,f1,output,m3);

  }

  /* the remaining parts are for frequency calculations with cyclic symmetry only */

  if((*nmethod!=2)&&(*nmethod!=13)){fclose(f1);return;}

  /* storing the maximum displacements of the nodes in the base sector
     (components, magnitude) */
  
  if((strcmp1(&filab[1566],"MAXU")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[1566],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  MDISP       4    1\n");
    fprintf(f1," -5  DX          1    4    1    0\n");
    fprintf(f1," -5  DY          1    4    2    0\n");
    fprintf(f1," -5  DZ          1    4    3    0\n");
    fprintf(f1," -5  ANG         1    4    4    0\n");
    
    ncomp=4;
    ifield[0]=1;icomp[0]=1;
    ifield[1]=1;icomp[1]=2;
    ifield[2]=1;icomp[2]=3;
    ifield[3]=1;icomp[3]=0;
    nfield[0]=4;nfield[1]=4;

    frdselect(vmax,vmax,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncomp,ifield,icomp,
	      nfield,&iselect,m2,f1,output,m3);

  }

  /* storing the worst principal stress at the nodes
     in the basis sector (components, magnitude)

     the worst principal stress is the maximum of the
     absolute value of all principal stresses, times
     its original sign */
  
  if((strcmp1(&filab[1653],"MAXS")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[1653],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  MSTRESS     7    1\n");
    fprintf(f1," -5  SXX         1    4    1    1\n");
    fprintf(f1," -5  SYY         1    4    2    2\n");
    fprintf(f1," -5  SZZ         1    4    3    3\n");
    fprintf(f1," -5  SXY         1    4    1    2\n");
    fprintf(f1," -5  SYZ         1    4    2    3\n");
    fprintf(f1," -5  SZX         1    4    3    1\n");
    fprintf(f1," -5  MAG         1    4    0    0\n");
    
    ncomp=7;
    ifield[0]=1;icomp[0]=1;
    ifield[1]=1;icomp[1]=2;
    ifield[2]=1;icomp[2]=3;
    ifield[3]=1;icomp[3]=4;
    ifield[4]=1;icomp[4]=6;
    ifield[5]=1;icomp[5]=5;
    ifield[6]=1;icomp[6]=0;
    nfield[0]=7;nfield[1]=7;

    frdselect(stnmax,stnmax,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncomp,ifield,icomp,
	      nfield,&iselect,m2,f1,output,m3);

  }

  /* storing the worst principal strain at the nodes
     in the basis sector (components, magnitude)

     the worst principal strain is the maximum of the
     absolute value of all principal strains, times
     its original sign */
  
  if((strcmp1(&filab[2523],"MAXE")==0)&&(*ithermal!=2)){
    iselect=1;
    
    frdset(&filab[2523],set,&iset,istartset,iendset,ialset,
	   inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	   ngraph);
    
    frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	      &noutloc,description,kode,nmethod,f1,output,istep,iinc);

    fprintf(f1," -4  MSTRAIN     7    1\n");
    fprintf(f1," -5  EXX         1    4    1    1\n");
    fprintf(f1," -5  EYY         1    4    2    2\n");
    fprintf(f1," -5  EZZ         1    4    3    3\n");
    fprintf(f1," -5  EXY         1    4    1    2\n");
    fprintf(f1," -5  EYZ         1    4    2    3\n");
    fprintf(f1," -5  EZX         1    4    3    1\n");
    fprintf(f1," -5  MAG         1    4    0    0\n");
    
    ncomp=7;
    ifield[0]=1;icomp[0]=1;
    ifield[1]=1;icomp[1]=2;
    ifield[2]=1;icomp[2]=3;
    ifield[3]=1;icomp[3]=4;
    ifield[4]=1;icomp[4]=6;
    ifield[5]=1;icomp[5]=5;
    ifield[6]=1;icomp[6]=0;
    nfield[0]=7;nfield[1]=7;

    frdselect(eenmax,eenmax,&iset,&nkcoords,inum,m1,istartset,iendset,
	      ialset,ngraph,&ncomp,ifield,icomp,
	      nfield,&iselect,m2,f1,output,m3);

  }
  
  fclose(f1);
  return;
  
}
