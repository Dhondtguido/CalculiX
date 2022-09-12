/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

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

void getlocalresults(ITG **integerglobp,double **doubleglobp,ITG *nktet,
                     double *cotet,double *h,ITG *netet_,ITG *kontet,
                     ITG *ifatet,double *planfa,ITG *kontetor){
    
  ITG *ifatetl=NULL,*kon=NULL,*ipkon=NULL,*kontyp=NULL,*iparent=NULL,
    *ielemnr=NULL,*integerglob=NULL,nfaces,netet,*nnx=NULL,
    *nny=NULL,*nnz=NULL,i,j,nkon,ne,n1,n2,n3,n4,kflag=2;
    
  double *doubleglob=NULL,*cgtet=NULL,*x=NULL,*y=NULL,*z=NULL,
    *xo=NULL,*yo=NULL,*zo=NULL;
    
  integerglob=*integerglobp;doubleglob=*doubleglobp;

  nfaces=0;
  netet=0;

  /* catalogueing all tetrahedral elements in successive order */

  //  NNEW(kon,ITG,10**netet_);
  NNEW(kon,ITG,4**netet_);
  NNEW(ifatetl,ITG,4**netet_);
  NNEW(ipkon,ITG,*netet_);
  NNEW(kontyp,ITG,*netet_);
  NNEW(iparent,ITG,*netet_);
  NNEW(ielemnr,ITG,*netet_);

  for(i=0;i<*netet_;i++){
    if(kontet[4*i]==0) continue;
    for(j=0;j<4;j++){
      //      kon[10*netet+j]=kontet[4*i+j];
      kon[4*netet+j]=kontet[4*i+j];
      ifatetl[4*netet+j]=ifatet[4*i+j];
      if(ifatet[4*i+j]>nfaces) nfaces=ifatet[4*i+j];
    }

    /* storing the midnodes */

    //	for(j=0;j<6;j++){kon[10*netet+4+j]=kontetor[6*i+j];}

    //    ipkon[netet]=10*netet;
    ipkon[netet]=4*netet;
    kontyp[netet]=3;
    iparent[netet]=netet+1;
    ielemnr[netet]=i+1;
    netet++;
  }
  //  nkon=10*netet;
  nkon=4*netet;
  ne=netet;

  RENEW(kon,ITG,nkon);
  RENEW(ifatetl,ITG,4*netet);
  RENEW(ipkon,ITG,ne);
  RENEW(kontyp,ITG,ne);
  RENEW(iparent,ITG,ne);
  RENEW(ielemnr,ITG,ne);
    
  /* calculating the center of gravity of the tetrahedra */
    
  NNEW(cgtet,double,3*netet);
  for(i=0;i<netet;i++){
    n1=kon[4*i]-1;
    n2=kon[4*i+1]-1;
    n3=kon[4*i+2]-1;
    n4=kon[4*i+3]-1;
    cgtet[3*i]=(cotet[3*n1]+cotet[3*n2]+cotet[3*n3]+cotet[3*n4])/4.;
    cgtet[3*i+1]=(cotet[3*n1+1]+cotet[3*n2+1]+cotet[3*n3+1]+cotet[3*n4+1])/4.;
    cgtet[3*i+2]=(cotet[3*n1+2]+cotet[3*n2+2]+cotet[3*n3+2]+cotet[3*n4+2])/4.;
  }
    
  /* initialization of additional fields */
    
  NNEW(x,double,netet);
  NNEW(y,double,netet);
  NNEW(z,double,netet);
  NNEW(xo,double,netet);
  NNEW(yo,double,netet);
  NNEW(zo,double,netet);
  NNEW(nnx,ITG,netet);
  NNEW(nny,ITG,netet);
  NNEW(nnz,ITG,netet);
  for(i=0;i<netet;i++){
    nnx[i]=i+1;
    nny[i]=i+1;
    nnz[i]=i+1;
    x[i]=cgtet[3*i];
    y[i]=cgtet[3*i+1];
    z[i]=cgtet[3*i+2];
    xo[i]=x[i];
    yo[i]=y[i];
    zo[i]=z[i];
  }
  FORTRAN(dsort,(x,nnx,&netet,&kflag));
  FORTRAN(dsort,(y,nny,&netet,&kflag));
  FORTRAN(dsort,(z,nnz,&netet,&kflag));
  SFREE(cgtet);

  /* storing the global data in a common block */
    
  NNEW(integerglob,ITG,5+3*ne+nkon+8*netet);
    
  integerglob[0]=*nktet;
  integerglob[1]=netet;
  integerglob[2]=ne;
  integerglob[3]=nkon;
  integerglob[4]=nfaces;
  memcpy(&integerglob[5],&nnx[0],sizeof(ITG)*netet);
  memcpy(&integerglob[netet+5],&nny[0],sizeof(ITG)*netet);
  memcpy(&integerglob[2*netet+5],&nnz[0],sizeof(ITG)*netet);
  memcpy(&integerglob[3*netet+5],&ifatetl[0],sizeof(ITG)*4*netet);
  memcpy(&integerglob[7*netet+5],&kontyp[0],sizeof(ITG)*ne);
  memcpy(&integerglob[ne+7*netet+5],&ipkon[0],sizeof(ITG)*ne);
  memcpy(&integerglob[2*ne+7*netet+5],&kon[0],sizeof(ITG)*nkon);
  memcpy(&integerglob[nkon+2*ne+7*netet+5],&iparent[0],sizeof(ITG)*netet);
  memcpy(&integerglob[nkon+2*ne+8*netet+5],&ielemnr[0],sizeof(ITG)*ne);
    
  NNEW(doubleglob,double,4**nktet+4*nfaces+6*netet);
    
  memcpy(&doubleglob[0],&x[0],sizeof(double)*netet);
  memcpy(&doubleglob[netet],&y[0],sizeof(double)*netet);
  memcpy(&doubleglob[2*netet],&z[0],sizeof(double)*netet);
  memcpy(&doubleglob[3*netet],&xo[0],sizeof(double)*netet);
  memcpy(&doubleglob[4*netet],&yo[0],sizeof(double)*netet);
  memcpy(&doubleglob[5*netet],&zo[0],sizeof(double)*netet);
  memcpy(&doubleglob[6*netet],&planfa[0],sizeof(double)*4*nfaces);
  memcpy(&doubleglob[4*nfaces+6*netet],&h[0],sizeof(double)**nktet);
  memcpy(&doubleglob[*nktet+4*nfaces+6*netet],&cotet[0],sizeof(double)*3**nktet);
    
  SFREE(nnx);SFREE(nny);SFREE(nnz);SFREE(kontyp);SFREE(ipkon);
  SFREE(kon);SFREE(iparent);SFREE(ielemnr);SFREE(ifatetl);

  SFREE(x);SFREE(y);SFREE(z);SFREE(xo);SFREE(yo);SFREE(zo);

  *integerglobp=integerglob;*doubleglobp=doubleglob;
 
}
