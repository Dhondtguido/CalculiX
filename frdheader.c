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

void frdheader(ITG *icounter,double *oner,double *time,double *pi,
	       ITG *noddiam,double *cs,ITG *null,ITG *mode,
	       ITG *noutloc,char *description,ITG *kode,ITG *nmethod,
               FILE *f1, char *output,ITG *istep,ITG *iinc){

  char tmp[132],text[132];

  ITG i,ncomma;
  
  /* icounter counts the number of loadcases in the frd-file 
     kode counts the number of increments in the frd-file */

  strcpy1(&text[0],"    1PSTEP",10);
  for(i=10;i<70;i++)text[i]=' ';text[70]='\0';
  (*icounter)++;
  sprintf(&text[24],"%12" ITGFORMAT "",*icounter);
  sprintf(&text[36],"%12" ITGFORMAT "",*iinc);
  sprintf(&text[48],"%12" ITGFORMAT "",*istep);text[60]=' ';
  fprintf(f1,"%s\n",text);

  /* additional headers for frequency calculations */

  if((*nmethod==2)||(*nmethod==6)||(*nmethod==7)||(*nmethod==13)||
     ((*nmethod==12)&&(*noddiam>-1))){
    strcpy1(&text[0],"    1PGM",8);
    for(i=8;i<70;i++)text[i]=' ';text[70]='\0';
    sprintf(&text[24],"%12.6E",*oner);text[36]=' ';
    fprintf(f1,"%s\n",text);

    strcpy1(&text[0],"    1PGK",8);
    for(i=8;i<70;i++)text[i]=' ';text[70]='\0';
    sprintf(&text[24],"%12.6E",(*time*2.**pi)*(*time*2.**pi));text[36]=' ';
    fprintf(f1,"%s\n",text);

    strcpy1(&text[0],"    1PHID",9);
    for(i=9;i<70;i++)text[i]=' ';text[70]='\0';
    sprintf(&text[24],"%12" ITGFORMAT "",*noddiam);text[36]=' ';
    fprintf(f1,"%s\n",text);

    /* additional headers for cyclic symmetry calculations */

    if((*noddiam>=0)&&(cs!=NULL)){
      strcpy1(&text[0],"    1PAX",8);
      for(i=8;i<24;i++)text[i]=' ';
      sprintf(&text[24],"%12.5E",cs[5]);
      sprintf(&text[36],"%12.5E",cs[6]);
      sprintf(&text[48],"%12.5E",cs[7]);
      sprintf(&text[60],"%12.5E",cs[8]);
      sprintf(&text[72],"%12.5E",cs[9]);
      sprintf(&text[84],"%12.5E",cs[10]);
      fprintf(f1,"%s\n",text);
    }

    strcpy1(&text[0],"    1PSUBC",10);
    for(i=10;i<70;i++)text[i]=' ';text[70]='\0';
    sprintf(&text[24],"%12" ITGFORMAT "",*null);text[36]=' ';
    fprintf(f1,"%s\n",text);

    strcpy1(&text[0],"    1PMODE",10);
    for(i=10;i<70;i++)text[i]=' ';text[70]='\0';
    sprintf(&text[24],"%12" ITGFORMAT "",*mode+1);text[36]=' ';
    fprintf(f1,"%s\n",text);
  }

#ifdef COMPANY
  writeBasisParameter(f1,istep,iinc);
#endif	  
     
  /* 100CL line */
 
  for(i=0;i<75;i++)text[i]=' ';
  if(abs(*nmethod)==1){
    strcpy1(&text[0],"  100CL       .00000E+00                                 0    1",63);
  }else if((*nmethod==2)||(*nmethod==13)){
    strcpy1(&text[0],"  100CL       .00000E+00                                 2    1",63);
  }else if(*nmethod==3){
    strcpy1(&text[0],"  100CL       .00000E+00                                 4    1",63);
  }else if((*nmethod==4)||(*nmethod==5)){
    strcpy1(&text[0],"  100CL       .00000E+00                                 1    1",63);
  }else{
    strcpy1(&text[0],"  100CL       .00000E+00                                 3    1",63);
  }

  sprintf(tmp,"%12" ITGFORMAT "",*noutloc);
  strcpy1(&text[24],tmp,12);
  strcpy1(&text[36],description,12);
  if((*nmethod==2)||(*nmethod==13))strcpy1(&text[63],"MODAL",5);
  if(strcmp1(output,"asc")==0){
    strcpy1(&text[74],"1",1);
  }else if(strcmp1(output,"bin")==0){
    strcpy1(&text[74],"2",1);
  }else{
    strcpy1(&text[74],"3",1);
  }
  sprintf(tmp,"%5" ITGFORMAT "",100+(*kode));
  strcpy1(&text[7],tmp,5);
//  sprintf(tmp,"%12.5E",*time);

// change 21 Oct 2018  if((*time<=0.)||(*nmethod==2)){
  if(*time<=0.){
      sprintf(tmp,"%12.5E",*time);
  }else if((log10(*time)>=0)&&(log10(*time)<10.)){
      ncomma=10-floor(log10(*time)+1.);
      if(ncomma==0){
	  sprintf(tmp,"%12.0f",*time);
      }else if(ncomma==1){
	  sprintf(tmp,"%12.1f",*time);
      }else if(ncomma==2){
	  sprintf(tmp,"%12.2f",*time);
      }else if(ncomma==3){
	  sprintf(tmp,"%12.3f",*time);
      }else if(ncomma==4){
	  sprintf(tmp,"%12.4f",*time);
      }else if(ncomma==5){
	  sprintf(tmp,"%12.5f",*time);
      }else if(ncomma==6){
	  sprintf(tmp,"%12.6f",*time);
      }else if(ncomma==7){
	  sprintf(tmp,"%12.7f",*time);
      }else if(ncomma==8){
	  sprintf(tmp,"%12.8f",*time);
      }else{
	  sprintf(tmp,"%12.9f",*time);
      }
  }else{
      sprintf(tmp,"%12.5E",*time);
  }

  strcpy1(&text[12],tmp,12);
  sprintf(tmp,"%5" ITGFORMAT "",*kode);
  strcpy1(&text[58],tmp,5);
  text[75]='\0';
  fprintf(f1,"%s\n",text);

}
