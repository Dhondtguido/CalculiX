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
#include <ctype.h>
#include "CalculiX.h"

void readinput(char *jobnamec, char **inpcp, ITG *nline, ITG *nset,
	       ITG *ipoinp, ITG **inpp, ITG **ipoinpcp, ITG *ithermal,
               ITG *nuel,ITG *inp_size){

  /*   reads and stores the input deck in inpcp; determines the
       number of sets  */

  FILE *f1[10];

  char buff[1320]="", fninp[132]="", includefn[132]="", *inpc=NULL,
       textpart[2112]="",*set=NULL;

  ITG i,j,k,n,in=0,nlinemax=100000,irestartread,irestartstep,
    icntrl,nload,nforc,nboun,nk,ne,nmpc,nalset,nmat,ntmat,npmat,
    norien,nam,nprint,mi[3],ntrans,ncs,namtot,ncmat,memmpc,ne1d,
    ne2d,nflow,*meminset=NULL,*rmeminset=NULL, *inp=NULL,ntie,
    nener,nstate,nentries=18,ifreeinp,ikey,lincludefn,nslavs,
    nbody,ncharmax=1000000,*ipoinpc=NULL,ichangefriction=0,nkon,
    ifile,mcs,initialtemperature=0,nprop,mortar,ifacecount,
    nintpoint,infree[4],iheading=0,ichangesurfacebehavior=0,
    nef,mpcend,nheading_,network,nfc,ndc,iprestr; 

  /* initialization */

  /* nentries is the number of different keyword cards for which
     the input deck order is important, cf keystart.f */

  NNEW(inpc,char,ncharmax);
  NNEW(ipoinpc,ITG,nlinemax+1);
  NNEW(inp,ITG,3*nlinemax);
  *nline=0;
  for(i=0;i<2*nentries;i++){ipoinp[i]=0;}
  ifreeinp=1;
  ikey=0;

  /* opening the input file */

  strcpy2(fninp,jobnamec,132);
  strcat(fninp,".inp");
  if((f1[in]=fopen(fninp,"r"))==NULL){
      printf(" *ERROR in readinput: cannot open file %s\n",fninp);
      exit(0);
  }

  /* starting to read the input file */

  do{
      if(fgets(buff,1320,f1[in])==NULL){
	  fclose(f1[in]);
	  if(in!=0){
	      in--;
	      continue;
	  }
	  else{break;}
      }
	  
      /* check for heading lines: should not be changed */

      if(iheading==1){
	  if((buff[0]=='*')&&(buff[1]!='*')){
	      iheading=0;
	  }
      }

      /* storing the significant characters */
      /* get rid of blanks  */
	
      k=0;
      i=-1;
      if(iheading==0){
	  do{
	      i++;
	      if((buff[i]=='\0')||(buff[i]=='\n')||(buff[i]=='\r')||(k==1320)) break;
	      if((buff[i]==' ')||(buff[i]=='\t')) continue;
	      buff[k]=buff[i];
	      k++;
	  }while(1);
      }else{
	  do{
	      i++;
	      if((buff[i]=='\0')||(buff[i]=='\n')||(buff[i]=='\r')||(k==1320)) break;
	      buff[k]=buff[i];
	      k++;
	  }while(1);
      }
	
      /* check for blank lines and comments */

      if(k==0) continue;
      if(strcmp1(&buff[0],"**")==0) continue;

      /* changing to uppercase except filenames */

      if(iheading==0){
	  j=0;
	  ifile=0;
	  do{
	      if(j>=5){
		  if(strcmp1(&buff[j-5],"FILE=")==0) ifile=1;
	      }
	      if(j>=6){
		  if(strcmp1(&buff[j-6],"INPUT=")==0) ifile=1;
	      }
	      if(j>=7){
		  if(strcmp1(&buff[j-7],"OUTPUT=")==0) ifile=1;
	      }
	      if(j>=9){
		  if(strcmp1(&buff[j-9],"FILENAME=")==0) ifile=1;
	      }
	      if(ifile==1){
		  do{
		      if(strcmp1(&buff[j],",")!=0){
			  j++;
		      }else{
			  ifile=0;
			  break;
		      }
		  }while(j<k);
	      }else{
		  buff[j]=toupper(buff[j]);
	      }
	      j++;
	  }while(j<k);
      }

      /* check for a *HEADING card */

      if(strcmp1(&buff[0],"*HEADING")==0){
	  iheading=1;
      }

      /* check for a *KINEMATIC or *DISTRIBUTING card and change
         the asterisk into a C (is a "dependent" card of the 
         *COUPLING card */

      if((strcmp1(&buff[0],"*KINEMATIC")==0)||
         ((strcmp1(&buff[0],"*DISTRIBUTING")==0)&&
          (strcmp1(&buff[0],"*DISTRIBUTINGCOUPLING")!=0)))
         {
	  buff[0]='C';
      }
	  
      /* check for include statements */
	  
      if(strcmp1(&buff[0],"*INCLUDE")==0){
	  lincludefn=k;
	  FORTRAN(includefilename,(buff,includefn,&lincludefn));
          includefn[lincludefn]='\0';
	  in++;
	  if(in>9){
	      printf(" *ERROR in readinput: include statements can \n not be cascaded over more than 9 levels\n");
	  }
	  if((f1[in]=fopen(includefn,"r"))==NULL){
	      printf(" *ERROR in readinput: cannot open file %s\n",includefn);
	      exit(0);
	  }
          continue;
      }

      /* adding a line */
	  
      (*nline)++;
      if(*nline>nlinemax){
	  nlinemax=(ITG)(1.1*nlinemax);
	  RENEW(ipoinpc,ITG,nlinemax+1);
	  RENEW(inp,ITG,3*nlinemax);
      }

      /* checking the total number of characters */

      if(ipoinpc[*nline-1]+k>ncharmax){
	  ncharmax=(ITG)(1.1*ncharmax);
	  RENEW(inpc,char,ncharmax);
      }
	  
      /* copying into inpc */

      for(j=0;j<k;j++){
	  inpc[ipoinpc[*nline-1]+j]=buff[j];
      }
      ipoinpc[*nline]=ipoinpc[*nline-1]+k;

      /* counting sets */
      
      if(strcmp1(&buff[0],"*AMPLITUDE")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"AMPLITUDE",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CHANGEFRICTION")==0){
	ichangefriction=1;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CHANGESURFACEBEHAVIOR")==0){
	ichangesurfacebehavior=1;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CONDUCTIVITY")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CONTACTDAMPING")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CONTACTPAIR")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"CONTACTPAIR",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*COUPLING")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"COUPLING",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CREEP")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CYCLICHARDENING")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DAMPING")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DISTRIBUTION")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"DISTRIBUTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ELASTIC")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DEFORMATIONPLASTICITY")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DENSITY")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DEPVAR")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ELASTIC")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ELECTRICALCONDUCTIVITY")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if((strcmp1(&buff[0],"*ELEMENT")==0)&&
              (strcmp1(&buff[0],"*ELEMENTOUTPUT")!=0)){
        (*nset)++;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"ELEMENT",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ELSET")==0){
        (*nset)++;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"ELSET",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*EXPANSION")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*FLUIDCONSTANTS")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if((strcmp1(&buff[0],"*FRICTION")==0)&&(ichangefriction==0)){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*GAPCONDUCTANCE")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*GAPHEATGENERATION")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*HYPERELASTIC")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*HYPERFOAM")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*INITIALCONDITIONS")==0){
	  FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INITIALCONDITIONS",
			    nline,&ikey));
	  FORTRAN(splitline,(buff,textpart,&n));
	  for(i=0;i<n;i++){
	      if(strcmp1(&textpart[(long long)132*i],"TYPE=TEMPERATURE")==0){
		  initialtemperature=1;
	      }
          }
      }
      else if(strcmp1(&buff[0],"*MAGNETICPERMEABILITY")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*MATERIAL")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if((strcmp1(&buff[0],"*NODE")==0)&&
	      (strcmp1(&buff[0],"*NODEPRINT")!=0)&&
	      (strcmp1(&buff[0],"*NODEOUTPUT")!=0)&&
	      (strcmp1(&buff[0],"*NODEFILE")!=0)){
        (*nset)++;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"NODE",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*NSET")==0){
        (*nset)++;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"NSET",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ORIENTATION")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"ORIENTATION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*PLASTIC")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*RATEDEPENDENT")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*RESTART")==0){
	  irestartread=0;
	  irestartstep=0;
	  strcpy1(&buff[k]," ",1);
	  FORTRAN(splitline,(buff,textpart,&n));
	  for(i=0;i<n;i++){
	      if(strcmp1(&textpart[(long long)132*i],"READ")==0){
		  irestartread=1;
	      }
	      if(strcmp1(&textpart[(long long)132*i],"STEP")==0){
		  irestartstep=atoi(&textpart[(long long)132*i+5]);
	      }
          }
          if(irestartread==1){
            icntrl=0;
            FORTRAN(restartshort,(nset,&nload,&nbody,&nforc,&nboun,&nk,
              &ne,&nmpc,&nalset,&nmat,&ntmat,&npmat,&norien,&nam,
              &nprint,mi,&ntrans,&ncs,&namtot,&ncmat,&memmpc,
              &ne1d,&ne2d,&nflow,set,meminset,rmeminset,jobnamec,
	      &irestartstep,&icntrl,ithermal,&nener,&nstate,&ntie,
	      &nslavs,&nkon,&mcs,&nprop,&mortar,&ifacecount,&nintpoint,
	      infree,&nef,&mpcend,&nheading_,&network,&nfc,&ndc,&iprestr));
            FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"RESTART,READ",
                              nline,&ikey));
	  }
          else{
            FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                              nline,&ikey));
          }

      }
      else if(strcmp1(&buff[0],"*SPECIFICGASCONSTANT")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*SPECIFICHEAT")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*SUBMODEL")==0){
	(*nset)+=2;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*SURFACEINTERACTION")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*SURFACEBEHAVIOR")==0){
	  if(ichangesurfacebehavior==0){
	      FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
	  }else{
	      FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
				nline,&ikey));
	  }
      }
      else if(strcmp1(&buff[0],"*SURFACE")==0){
        (*nset)++;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"SURFACE",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*TIE")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"TIE",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*TRANSFORM")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"TRANSFORM",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*USERELEMENT")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"USERELEMENT",
                          nline,&ikey));
	(*nuel)++;
      }
      else if(strcmp1(&buff[0],"*USERMATERIAL")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*")==0){
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                          nline,&ikey));

        /* checking whether the calculation is mechanical,
           thermal or thermomechanical: needed to know
           which mpc's to apply to 2-D elements */

	if((strcmp1(&buff[0],"*STATIC")==0)||
	   (strcmp1(&buff[0],"*VISCO")==0)||
	   (strcmp1(&buff[0],"*DYNAMIC")==0)){
	    if(ithermal[1]==0){
		if(initialtemperature==1)ithermal[1]=1;
	    }else if(ithermal[1]==2){
		ithermal[1]=3;
	    }
	}else if(strcmp1(&buff[0],"*HEATTRANSFER")==0){
	    if(ithermal[1]<2) ithermal[1]=ithermal[1]+2;
	}else if(strcmp1(&buff[0],"*COUPLEDTEMPERATURE-DISPLACEMENT")==0){
	    ithermal[1]=3;
	}else if(strcmp1(&buff[0],"*UNCOUPLEDTEMPERATURE-DISPLACEMENT")==0){
	    ithermal[1]=3;
	}else if(strcmp1(&buff[0],"*ELECTROMAGNETICS")==0){
	    ithermal[1]=3;
	}
      }
  }while(1);

  inp[3*ipoinp[2*ikey-1]-2]=*nline;
  RENEW(inpc,char,(long long)132**nline);
  
  *inp_size=3*ipoinp[2*ikey-1];
  RENEW(inp,ITG,*inp_size);
  
  *inpcp=inpc;
  *ipoinpcp=ipoinpc;
  *inpp=inp;
  
//  FORTRAN(writeinput,(inpc,ipoinp,inp,nline,&ipoinp[2*ikey-1],ipoinpc));

  return;

}






