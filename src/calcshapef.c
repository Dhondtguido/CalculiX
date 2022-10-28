/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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

void calcshapef(ITG *nvar_,ITG *ipvar,double **varp,ITG *ne,
		char *lakon,double *co,ITG *ipkon,ITG *kon,
                ITG *nelemface,char *sideface,ITG *nface,
                ITG *nvarf_,ITG *ipvarf,double **varfp,
		ITG *iturbulent,double *yy){
    
    /*   determines the shape functions and their derivatives for
	 a fluid mesh and stores the results in ipvar and var */
    
    ITG i,j,k,kk,nope,mint3d,iflag=3,nmethod=1,nvar,indexe,nelem,
        mint2d,nopes,ig,nvarf,id;
    
    ITG ifaceq[48]=
	{4,3,2,1,11,10,9,12,
	 5,6,7,8,13,14,15,16,
	 1,2,6,5,9,18,13,17,
	 2,3,7,6,10,19,14,18,
	 3,4,8,7,11,20,15,19,
	 4,1,5,8,12,17,16,20};

    ITG ifacet[24]=
	{1,3,2,7,6,5,
	 1,2,4,5,9,8,
	 2,3,4,6,10,9,
	 1,4,3,8,10,7};

    ITG ifacew[40]=
	{1,3,2,9,8,7,0,0,
	 4,5,6,10,11,12,0,0,
	 1,2,5,4,7,14,10,13,
	 2,3,6,5,8,15,11,14,
	 4,6,3,1,12,15,9,13};
    
    double *var=NULL,xl[60],xi,et,ze,xsj,xl2[24],xs2[21],shp2[56],xsj2[3],
      xi3d,et3d,ze3d,weight,*varf=NULL,yyl[8],y;
    
    /* gauss points */
    
    double gauss2d1[2]={0.,0.};
    
    double gauss2d2[8]={
	-0.577350269189626,-0.577350269189626,
	0.577350269189626,-0.577350269189626,
	-0.577350269189626,0.577350269189626,
	0.577350269189626,0.577350269189626};
    
    double gauss2d4[2]={0.333333333333333,0.333333333333333};
    
    double gauss3d1[3]={0.,0.,0.};
    
    double gauss3d2[24]={
	-0.577350269189626,-0.577350269189626,-0.577350269189626,
	    0.577350269189626,-0.577350269189626,-0.577350269189626,
	    -0.577350269189626,0.577350269189626,-0.577350269189626,
	    0.577350269189626,0.577350269189626,-0.577350269189626,
	    -0.577350269189626,-0.577350269189626,0.577350269189626,
	    0.577350269189626,-0.577350269189626,0.577350269189626,
	    -0.577350269189626,0.577350269189626,0.577350269189626,
	    0.577350269189626,0.577350269189626,0.577350269189626};
    
    double gauss3d4[3]={0.25,0.25,0.25};
    
    double gauss3d7[6]={
	0.333333333333333,0.333333333333333,-0.577350269189626,
	    0.333333333333333,0.333333333333333,0.577350269189626};
 
    double weight2d1[1]={4.};
 
    double weight2d2[4]={1.,1.,1.,1.};
 
    double weight2d4[1]={0.5};
 
    double weight3d1[1]={8.};
 
    double weight3d2[8]={1.,1.,1.,1.,1.,1.,1.,1.};
 
    double weight3d4[1]={0.166666666666667};
 
    double weight3d7[2]={0.5,0.5};

    double xlocal8r[18]={
      0.000000000000000e+0, 0.000000000000000e+0,-0.100000000000000e+1
     , 0.000000000000000e+0, 0.000000000000000e+0, 0.100000000000000e+1
     , 0.000000000000000e+0,-0.100000000000000e+1, 0.000000000000000e+0
     , 0.100000000000000e+1, 0.000000000000000e+0, 0.000000000000000e+0
     , 0.000000000000000e+0, 0.100000000000000e+1, 0.000000000000000e+0
     ,-0.100000000000000e+1, 0.000000000000000e+0,0.000000000000000e+0};

    double xlocal8[72]={
     -0.577350269189626e+0, 0.577350269189626e+0,-0.100000000000000e+1
     , 0.577350269189626e+0, 0.577350269189626e+0,-0.100000000000000e+1
     ,-0.577350269189626e+0,-0.577350269189626e+0,-0.100000000000000e+1
     , 0.577350269189626e+0,-0.577350269189626e+0,-0.100000000000000e+1
     ,-0.577350269189626e+0,-0.577350269189626e+0, 0.100000000000000e+1
     , 0.577350269189626e+0,-0.577350269189626e+0, 0.100000000000000e+1
     ,-0.577350269189626e+0, 0.577350269189626e+0, 0.100000000000000e+1
     , 0.577350269189626e+0, 0.577350269189626e+0, 0.100000000000000e+1
     ,-0.577350269189626e+0,-0.100000000000000e+1,-0.577350269189626e+0
     , 0.577350269189626e+0,-0.100000000000000e+1,-0.577350269189626e+0
     ,-0.577350269189626e+0,-0.100000000000000e+1, 0.577350269189626e+0
     , 0.577350269189626e+0,-0.100000000000000e+1, 0.577350269189626e+0
     , 0.100000000000000e+1,-0.577350269189626e+0,-0.577350269189626e+0
     , 0.100000000000000e+1, 0.577350269189626e+0,-0.577350269189626e+0
     , 0.100000000000000e+1,-0.577350269189626e+0, 0.577350269189626e+0
     , 0.100000000000000e+1, 0.577350269189626e+0, 0.577350269189626e+0
     , 0.577350269189626e+0, 0.100000000000000e+1,-0.577350269189626e+0
     ,-0.577350269189626e+0, 0.100000000000000e+1,-0.577350269189626e+0
     , 0.577350269189626e+0, 0.100000000000000e+1, 0.577350269189626e+0
     ,-0.577350269189626e+0, 0.100000000000000e+1, 0.577350269189626e+0
     ,-0.100000000000000e+1, 0.577350269189626e+0,-0.577350269189626e+0
     ,-0.100000000000000e+1,-0.577350269189626e+0,-0.577350269189626e+0
     ,-0.100000000000000e+1, 0.577350269189626e+0, 0.577350269189626e+0
     ,-0.100000000000000e+1,-0.577350269189626e+0,0.577350269189626e+0};

    double xlocal4[12]={
      0.333333333333333e+0, 0.333333333333333e+0, 0.000000000000000e+0
     , 0.333333333333333e+0, 0.000000000000000e+0, 0.333333333333333e+0
     , 0.333333333333334e+0, 0.333333333333333e+0, 0.333333333333333e+0
     , 0.000000000000000e+0, 0.333333333333333e+0,0.333333333333333e+0};

    double xlocal6[15]={
      0.333333333333333e+0, 0.333333333333333e+0,-0.100000000000000e+1
     , 0.333333333333333e+0, 0.333333333333333e+0, 0.100000000000000e+1
     , 0.500000000000000e+0, 0.000000000000000e+0, 0.000000000000000e+0
     , 0.500000000000000e+0, 0.500000000000000e+0, 0.000000000000000e+0
     , 0.000000000000000e+0, 0.500000000000000e+0,0.000000000000000e+0};
    
    var=*varp;
    varf=*varfp;
    nvar=0;
    nvarf=0;
    
    /* loop over all elements */
    
    for(i=0;i<*ne;i++){
	
	/* check for fluid elements */
	
	if(strcmp1(&lakon[8*i],"F")!=0) continue;

        /* check whether element is actif */

	if(ipkon[i]<0) continue;

	/* storing the beginning of the field for element i */

	ipvar[i]=nvar;
	ipvarf[i]=nvarf;
	//	printf("calcshapef.c nelem %d index %d\n",i+1,nvarf);
	
	/* determining the number of nodes belonging to the element
	   (nope) and the number of integration points (mint3d) */
	
	if(strcmp1(&lakon[8*i+3],"8 ")==0){
	    nope=8;
	    nopes=4;
	    mint2d=4;
	    mint3d=8;
	}else if(strcmp1(&lakon[8*i+3],"8R")==0){
	    nope=8;
	    nopes=4;
	    mint2d=1;
	    mint3d=1;
	}else if(strcmp1(&lakon[8*i+3],"4")==0){
	    nope=4;
	    nopes=3;
	    mint2d=1;
	    mint3d=1;
	}else if(strcmp1(&lakon[8*i+3],"6 ")==0){
	    nope=6;
	    mint3d=2;
	}else{
	    continue;
	}

	/* copying the coordinates in a local field */
	
	indexe=ipkon[i];
	for(j=0;j<nope;j++){
	    for(k=0;k<3;k++){
		xl[3*j+k]=co[3*(kon[indexe+j]-1)+k];
	    }
	}

	/* copying the distance from a solid surface in a local field */

	if(*iturbulent>0){
	  for(j=0;j<nope;j++){
	    yyl[j]=yy[kon[indexe+j]-1];
	  }
	}

	/* loop over the integration points */

	for(kk=0;kk<mint3d;kk++){
	
	/* check size of var */
	
	    if(nvar+4*nope+2>*nvar_){
		*nvar_=(ITG)(1.1**nvar_+4*nope+2);
		RENEW(var,double,*nvar_);
	    }

	    if(strcmp1(&lakon[8*i+3],"8R")==0){
		xi=gauss3d1[3*kk];		
		et=gauss3d1[3*kk+1];		
		ze=gauss3d1[3*kk+2];
		weight=weight3d1[kk];
	    }else if(strcmp1(&lakon[8*i+3],"8")==0){
		xi=gauss3d2[3*kk];		
		et=gauss3d2[3*kk+1];		
		ze=gauss3d2[3*kk+2];
		weight=weight3d2[kk];
	    }else if(strcmp1(&lakon[8*i+3],"4")==0){
		xi=gauss3d4[3*kk];		
		et=gauss3d4[3*kk+1];		
		ze=gauss3d4[3*kk+2];
		weight=weight3d4[kk];
	    }else if(strcmp1(&lakon[8*i+3],"6")==0){
		xi=gauss3d7[3*kk];		
		et=gauss3d7[3*kk+1];		
		ze=gauss3d7[3*kk+2];
		weight=weight3d7[kk];
	    }

	    /* calculating the shape functions and their
               derivatives */

	    if(nope==8){
		FORTRAN(shape8h,(&xi,&et,&ze,xl,&xsj,&var[nvar],&iflag));
	    }else if(nope==4){
		FORTRAN(shape4tet,(&xi,&et,&ze,xl,&xsj,&var[nvar],&iflag));
	    }else if(nope==6){
		FORTRAN(shape6w,(&xi,&et,&ze,xl,&xsj,&var[nvar],&iflag));
	    }

	    /* check the Jacobian determinant */

	    if(xsj<1.e-20){
		printf(" *ERROR in calcshapef: nonpositive Jacobian\n");
		printf("        determinant in element %d\n\n",i);
		xsj=fabs(xsj);
		nmethod=0;
	    }

	    /* storing the Jacobian determinant*weight */
	    
	    var[nvar+4*nope]=xsj*weight;

	    /* calculating the distance of the integration point
               from the next solid surface */

	    if(*iturbulent>0){
	      y=0.;
	      for(j=0;j<nope;j++){
		y+=var[nvar+4*j+3]*yyl[j];
	      }
	      var[nvar+4*nope+1]=y;
	    }
	    
            /* updating the pointer */

	    nvar+=4*nope+2;

	}
    
	/* free stream or solid surface boundaries */
	
	if(*nface!=0){
	    
	    nelem=i+1;
	    FORTRAN(nident,(nelemface,&nelem,nface,&id));
	    
	    do{
		if(id==0) break;
		if(nelemface[id-1]!=nelem) break;
		ig=sideface[id-1]-'0';
		
		/* treatment of wedge faces */
		
		if(strcmp1(&lakon[8*i+3],"6")==0){
		    mint2d=1;
		    if(ig<=2){
			nopes=3;
		    }else{
			nopes=4;
		    }
		}
		
		/* storing the coordinates of the face nodes */
		
		if(nope==8){
		    for(j=0;j<nopes;j++){
			for(k=0;k<3;k++){
			    xl2[3*j+k]=co[3*(kon[indexe+
						 ifaceq[8*(ig-1)+j]-1]-1)+k];
			}
		    }
		}else if(nope==4){
		    for(j=0;j<nopes;j++){
			for(k=0;k<3;k++){
			    xl2[3*j+k]=co[3*(kon[indexe+
						 ifacet[6*(ig-1)+j]-1]-1)+k];
			}
		    }
		}else{
		    for(j=0;j<nopes;j++){
			for(k=0;k<3;k++){
			    xl2[3*j+k]=co[3*(kon[indexe+
						 ifacew[8*(ig-1)+j]-1]-1)+k];
			}
		    }
		}
		
		//		printf("calcshapef.c %d  %d\n",nvarf,4*nope+nopes+4);
		for(kk=0;kk<mint2d;kk++){
		    
		    /* check size of varf */
		    
		    if(nvarf+4*nope+nopes+4>*nvarf_){
			*nvarf_=(ITG)(1.1**nvarf_+4*nope+nopes+4);
			RENEW(varf,double,*nvarf_);
		    }
		    
		    if((strcmp1(&lakon[8*i+3],"8R")==0)||
		       ((strcmp1(&lakon[8*i+3],"6")==0)&&(nopes==4))){
			xi=gauss2d1[2*kk];
			et=gauss2d1[2*kk+1];
			weight=weight2d1[kk];
		    }else if(strcmp1(&lakon[8*i+3],"8")==0){
			xi=gauss2d2[2*kk];
			et=gauss2d2[2*kk+1];
			weight=weight2d2[kk];
		    }else if((strcmp1(&lakon[8*i+3],"4")==0)||
			     ((strcmp1(&lakon[8*i+3],"6")==0)&&(nopes==3))){
			xi=gauss2d4[2*kk];
			et=gauss2d4[2*kk+1];
			weight=weight2d4[kk];
		    }
		    
		    /* local surface normal */
		    
		    iflag=2;
		    if(nopes==4){
			FORTRAN(shape4q,(&xi,&et,xl2,xsj2,xs2,shp2,&iflag));
		    }else{
			FORTRAN(shape3tri,(&xi,&et,xl2,xsj2,xs2,shp2,&iflag));
		    }
		    
		    /* copying the shape function values and the Jacobian
		       determinant into varf */
		    
		    for(j=0;j<nopes;j++){
			varf[nvarf+j]=shp2[7*j+3];
		    }
		    nvarf+=nopes;
		    for(j=0;j<3;j++){
			varf[nvarf+j]=xsj2[j]*weight;
		    }
		    nvarf+=3;
		    
		    /* local coordinates of the surface integration point
		       within the local element coordinate system */
		    
		    iflag=3;
		    if(strcmp1(&lakon[8*i+3],"8R")==0){
			xi3d=xlocal8r[3*(ig-1)+3*kk];
			et3d=xlocal8r[3*(ig-1)+3*kk+1];
			ze3d=xlocal8r[3*(ig-1)+3*kk+2];
			FORTRAN(shape8h,(&xi3d,&et3d,&ze3d,xl,&xsj,
					 &varf[nvarf],&iflag));
		    }else if(strcmp1(&lakon[8*i+3],"8")==0){
			xi3d=xlocal8[12*(ig-1)+3*kk];
			et3d=xlocal8[12*(ig-1)+3*kk+1];
			ze3d=xlocal8[12*(ig-1)+3*kk+2];
			FORTRAN(shape8h,(&xi3d,&et3d,&ze3d,xl,&xsj,
					 &varf[nvarf],&iflag));
		    }else if(strcmp1(&lakon[8*i+3],"4")==0){
			xi3d=xlocal4[3*(ig-1)+3*kk];
			et3d=xlocal4[3*(ig-1)+3*kk+1];
			ze3d=xlocal4[3*(ig-1)+3*kk+2];
			FORTRAN(shape4tet,(&xi3d,&et3d,&ze3d,xl,&xsj,
					   &varf[nvarf],&iflag));
		    }else if(strcmp1(&lakon[8*i+3],"6")==0){
			xi3d=xlocal6[3*(ig-1)+3*kk];
			et3d=xlocal6[3*(ig-1)+3*kk+1];
			ze3d=xlocal6[3*(ig-1)+3*kk+2];
			FORTRAN(shape6w,(&xi3d,&et3d,&ze3d,xl,&xsj,
					 &varf[nvarf],&iflag));
		    }

		    /* calculating the distance of the integration point
		       from the next solid surface */

		    if(*iturbulent>0){
		      y=0.;
		      for(j=0;j<nope;j++){
			y+=varf[nvarf+4*j+3]*yyl[j];
		      }
		      varf[nvarf+4*nope]=y;
		    }
		    
		    nvarf+=4*nope+1;
		}
		id--;
	    }while(1);
	}
    }
    
    if(nmethod==0){FORTRAN(stop,());}
    
    RENEW(var,double,nvar);
    RENEW(varf,double,nvarf);
    
    *varp=var;
    *varfp=varf;
    
}
