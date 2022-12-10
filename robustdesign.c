/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#ifdef SPOOLES
#include "spooles.h"
#endif
#ifdef SGI
#include "sgi.h"
#endif
#ifdef TAUCS
#include "tau.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif

void robustdesign(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
		  ITG *ne,
		  ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
		  ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
		  ITG *nmpc,
		  ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
		  ITG *nelemload,char *sideload,double *xload,
		  ITG *nload,ITG *nactdof,
		  ITG *icol,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
		  ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
		  ITG *ilboun,
		  double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
		  double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
		  ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
		  double *t0,double *t1,double *t1old,
		  ITG *ithermal,double *prestr,ITG *iprestr,
		  double *vold,ITG *iperturb,double *sti,ITG *nzs, 
		  ITG *kode,char *filab,double *eme,
		  ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
		  ITG *nplkcon,
		  double **xstatep,ITG *npmat_,char *matname,ITG *isolver,
		  ITG *mi,ITG *ncmat_,ITG *nstate_,double *cs,ITG *mcs,
		  ITG *nkon,double **enerp,double *xbounold,
		  double *xforcold,double *xloadold,
		  char *amname,double *amta,ITG *namta,
		  ITG *nam,ITG *iamforc,ITG *iamload,
		  ITG *iamt1,ITG *iamboun,double *ttime,char *output,
		  char *set,ITG *nset,ITG *istartset,
		  ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
		  char *prset,ITG *nener,double *trab,
		  ITG *inotr,ITG *ntrans,double *fmpc,ITG *ipobody,
		  ITG *ibody,
		  double *xbody,ITG *nbody,double *xbodyold,double *timepar,
		  double *thicke,char *jobnamec,char *tieset,ITG *ntie,
		  ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
		  ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,
		  ITG *nobject,char **objectsetp,ITG *istat,char *orname,
		  ITG *nzsprevstep,ITG *nlabel,double *physcon,char *jobnamef,
		  ITG *iponor2d,ITG *knor2d,ITG *ne2d,ITG *iponoel2d,
		  ITG *inoel2d,
		  ITG *mpcend,ITG *irobustdesign,ITG *irandomtype,
		  double *randomval){
	     
  char description[13]="            ",*lakon=NULL,cflag[1]=" ",
    *lakonfa=NULL,*objectset=NULL,filabnew[5]="    ";
       
  ITG *inum=NULL,k,*irow=NULL,iinc=1,mode=-1,noddiam=-1,ngraph=1,ne0,
    *integerglob=NULL,*kon=NULL,*ipkon=NULL,*ielmat=NULL,i,
    ndesi,iobject,*iponoel=NULL,node,*nodedesi=NULL,*ipoface=NULL,*nodface=NULL,
    *inoel=NULL,icoordinate=0,*istartdesi=NULL,*ialdesi=NULL,
    *istartelem=NULL,*ialelem=NULL,inoelsize,*itmp=NULL,*ielorien=NULL,
    iglob=0,idesvar=0,inorm=0,irand=0,*nodedesiinv=NULL,
    iregion=0,*konfa=NULL,*ipkonfa=NULL,nsurfs,*iponoelfa=NULL,
    *inoelfa=NULL,*iponor=NULL,*iponexp=NULL,ifreemax,*ipretinfo=NULL,
    nfield,iforce,*iponod2dto3d=NULL,*iponk2dto3d=NULL,ishape=0,ndesibou,
    *nodedesibou=NULL,*nodedesiinvbou=NULL,nmethodnew=0,*neigh=NULL,
    *ipneigh=NULL,ifeasd=0;
      
  double *stn=NULL,*tper,*xdesi=NULL,ptime=0.,*doubleglob=NULL,*xstate=NULL,
    *ener=NULL,sigma=0,*extnor=NULL,dtime,time,*xnor=NULL,*cdni=NULL,
    *cdnr=NULL,*cdn=NULL,*qfx=NULL,*emn=NULL,*fni=NULL,*fnr=NULL,
    *eenmax=NULL,*veold=NULL,*stnmax=NULL,*vmax=NULL,*stni=NULL,
    *stnr=NULL,*vi=NULL,*vr=NULL,*qfn=NULL,*xstaten=NULL,*enern=NULL,
    *epn=NULL,*fn=NULL,*een=NULL,*v=NULL;
  
#ifdef SGI
  ITG token;
#endif

  irow=*irowp;ener=*enerp;xstate=*xstatep;ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;ielmat=*ielmatp;ielorien=*ielorienp;objectset=*objectsetp;

  tper=&timepar[1];

  time=*tper;
  dtime=*tper;

  ne0=*ne;

  /* determining the global values to be used as boundary conditions
     for a submodel */

  ITG irefine=0;
  getglobalresults(&jobnamec[396],&integerglob,&doubleglob,nboun,iamboun,xboun,
		   nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
                   ithermal,nk,t1,iamt1,&sigma,&irefine);
  
  /* check which design variables are active */
  
  for(i=0;i<*ntie;i++){
      if(strcmp1(&tieset[i*243+80],"D")==0){
	  if(strcmp1(&tieset[i*243],"COORDINATE")==0){
	      icoordinate=1;
	      break;
	  }else if(strcmp1(&tieset[i*243],"ORIENTATION")==0){
	      printf(" *ERROR in robustdesign: the ORIENTATION sensitivity was requested,\n");
	      FORTRAN(stop,());
	      break;
	  }
      }
  }
  
  /* in robust design anaylsis "edge preservation" option is always active */
  iregion=1;
  
  /* determining the elements belonging to a given node */
  
  NNEW(iponoel,ITG,*nk);
  NNEW(inoel,ITG,2**nkon);
  FORTRAN(elementpernode,(iponoel,inoel,lakon,ipkon,kon,ne));

  /* find the 
     - external faces belonging to a given node
     - nodes belonging to a given external surface */

  NNEW(ipoface,ITG,*nk);
  NNEW(nodface,ITG,5*6**ne);
  NNEW(konfa,ITG,8*6**ne);
  NNEW(ipkonfa,ITG,6**ne);
  NNEW(lakonfa,char,8*6**ne);
  FORTRAN(findextsurface,(nodface,ipoface,ne,ipkon,lakon,kon,
  			  konfa,ipkonfa,nk,lakonfa,&nsurfs,
  			  &ifreemax));
  RENEW(nodface,ITG,5*ifreemax);
  RENEW(konfa,ITG,8*nsurfs);
  RENEW(ipkonfa,ITG,nsurfs);
  RENEW(lakonfa,char,8*nsurfs);

  /* find the external faces belonging to a given node */

  NNEW(iponoelfa,ITG,*nk);
  NNEW(inoelfa,ITG,3*8*6**ne);
  FORTRAN(extfacepernode,(iponoelfa,inoelfa,lakonfa,ipkonfa,konfa,
  			  &nsurfs,&inoelsize));
  RENEW(inoelfa,ITG,3*inoelsize);
  
  /* determining the design variables */
  
  NNEW(nodedesi,ITG,*nk);
  NNEW(itmp,ITG,*nk);
  NNEW(nodedesiinv,ITG,*nk);
  
  if(*ne2d!=0){

    NNEW(iponod2dto3d,ITG,3**nk);
    NNEW(iponk2dto3d,ITG,*nk);
     
    FORTRAN(getdesiinfo2d,(set,istartset,iendset,ialset,nset,
			   mi,nactdof,&ndesi,nodedesi,ntie,tieset,
			   nodedesiinv,lakon,ipkon,kon,iponoelfa,
			   iponod2dto3d,iponor2d,knor2d,iponoel2d,
			   inoel2d,nobject,objectset,iponk2dto3d,ne,
			   jobnamef));
    						 
  
  }else{
  
    FORTRAN(getdesiinfo3d_robust,(set,istartset,iendset,ialset,nset,
				  mi,nactdof,&ndesi,nodedesi,ntie,tieset,
				  itmp,nmpc,nodempc,ipompc,nodedesiinv,
				  iponoel,inoel,lakon,ipkon,
				  kon,&iregion,ipoface,nodface,nk,
				  irandomtype,jobnamef));  
  }
  
  SFREE(itmp);
  RENEW(nodedesi,ITG,ndesi);

  /* determining the designvariables at the boundary of the variable field */
  /* in case of a conditional random field */ 
  
  if(irobustdesign[2]==1){
  
    NNEW(nodedesibou,ITG,*nk);
    NNEW(nodedesiinvbou,ITG,*nk);
  
    FORTRAN(getdesiinfobou,(&ndesibou,nodedesibou,nodedesiinv,
			    lakon,ipkon,kon,ipoface,nodface,
			    nodedesiinvbou,&ndesi,nodedesi,nk));

    RENEW(nodedesibou,ITG,ndesibou);
  
  }
  
  /* storing the elements to which each design variable belongs
     in field ialdesi */

  NNEW(istartdesi,ITG,ndesi+1);
  NNEW(ialdesi,ITG,*nkon);
  FORTRAN(elemperdesi,(&ndesi,nodedesi,iponoel,inoel,
		       istartdesi,ialdesi,lakon,ipkon,kon,
		       nodedesiinv,&icoordinate,&iregion));
  RENEW(ialdesi,ITG,istartdesi[ndesi]-1);
  	
  /* calculating the normal direction for every designvariable */
  
  NNEW(extnor,double,3**nk);
  
  FORTRAN(normalsonsurface_robust,(ipkon,kon,lakon,extnor,co,nk,ipoface,
    			           nodface,nactdof,mi,nodedesiinv,&iregion,
    			           iponoelfa,&ndesi,nodedesi,iponod2dto3d,
    			           ikboun,nboun,ne2d)); 
  
  /* if the sensitivity calculation is used in a optimization script
     this script usually contains a loop consisting of:
     1. a call to CalculiX to define the sensitivities
     2. a small modification of the surface geometry in a direction which
     decrease the objective function (only the design variables)
     3. a modification of the internal mesh in order to preserve
     mesh quality
     The latter point can be done by performing a linear elastic
     calculation in which the small modification in 2. is applied
     a *boundary condition and all other nodes (on the external 
     surface but no design variables) are fixed by *equation's
     in a direction normal to the surface. At corners and edges
     there my be more than one normal. The necessary equations are
     calculated in normalsforeq_se.f and stored in jobname.equ */

  NNEW(iponor,ITG,8*nsurfs);
  for(i=0;i<8*nsurfs;i++) iponor[i]=-1;
  NNEW(xnor,double,24*nsurfs);
  NNEW(iponexp,ITG,2**nk);
  NNEW(ipretinfo,ITG,*nk);
  
  FORTRAN(normalsforequ_se,(nk,co,iponoelfa,inoelfa,konfa,ipkonfa,lakonfa,
			    &nsurfs,iponor,xnor,nodedesiinv,jobnamef,
			    iponexp,nmpc,labmpc,ipompc,nodempc,ipretinfo,
			    kon,ipkon,lakon,iponoel,inoel,iponor2d,knor2d,
			    iponod2dto3d,ipoface,nodface));
    	  
  SFREE(konfa);SFREE(ipkonfa);SFREE(lakonfa);SFREE(iponor);SFREE(xnor);
  SFREE(iponoelfa);SFREE(inoelfa);SFREE(iponexp);SFREE(ipretinfo);
       
  /* createinum is called in order to determine the nodes belonging
     to elements; this information is needed in frd_se */
  
  NNEW(inum,ITG,*nk);
  FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
		      nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
		      ielprop,prop));
    				       
  /* storing the normal information in the frd-file for the optimizer */
      
  ++*kode;

  inorm=1;
  nfield=3;
  iforce=0;
  
  if(strcmp1(&filab[4],"I")==0){
  
    FORTRAN(map3dto1d2d,(extnor,ipkon,inum,kon,lakon,&nfield,nk,
			 ne,cflag,co,vold,&iforce,mi,ielprop,prop));
  }

  /* storing the coordinates and topology (if not already done so) */

  frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,&nmethodnew,
      kode,filabnew,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
      nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
      ntrans,orab,ielorien,norien,description,ipneigh,neigh,
      mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
      cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
      thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
      ielprop,prop,sti);
  
  frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,
    	  istep,
    	  &iinc,&mode,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
    	  istartset,iendset,ialset,jobnamec,output,
    	  extnor,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
    	  &icoordinate,&inorm,&irand,&ishape,&ifeasd); 
  inorm=0;

  /* storing the normal direction for every design variable */

  NNEW(xdesi,double,3*ndesi);
  for(k=0;k<ndesi;k++){
    node=nodedesi[k]-1;
    memcpy(&xdesi[3*k],&extnor[3*node],sizeof(double)*3);
  }
  
  /* calculation of gaussian random fields for robust optimization */
  
  randomfieldmain(kon,ipkon,lakon,ne,nmpc,nactdof,mi,nodedesi,&ndesi,
		  istartdesi,ialdesi,co,physcon,isolver,ntrans,nk,inotr,trab,jobnamec,
		  nboun,cs,mcs,inum,nmethod,kode,filab,nstate_,istep,description,set,
		  nset,iendset,output,istartset,ialset,extnor,irandomtype,randomval,
		  irobustdesign,&ndesibou,nodedesibou,nodedesiinvbou); 
    		       
  SFREE(inum);SFREE(extnor);
  if(irobustdesign[2]==1){SFREE(nodedesibou);SFREE(nodedesiinvbou);}	  
            
  SFREE(iponoel);SFREE(inoel);SFREE(nodedesiinv);
  
  if(*ne2d!=0){SFREE(iponod2dto3d);SFREE(iponk2dto3d);}
     
  // if(*nbody>0) SFREE(ipobody);

  SFREE(istartdesi);SFREE(ialdesi);SFREE(istartelem);SFREE(ialelem);

  if(icoordinate==1){
    SFREE(nodedesi);SFREE(xdesi);SFREE(ipoface);SFREE(nodface);
  }

  *irowp=irow;*enerp=ener;*xstatep=xstate;*ipkonp=ipkon;*lakonp=lakon;
  *konp=kon;*ielmatp=ielmat;*ielorienp=ielorien;*objectsetp=objectset;

  (*ttime)+=(*tper);
 
  return;
}
