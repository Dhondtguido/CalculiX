/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                     */

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

#ifdef CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT

#include<stdio.h>
#include<ctype.h>
#include<stdlib.h>
#include<string.h>

#if (defined _WIN32) && (! defined __CYGWIN__)
#include<windows.h>
#define dlopen(name, mode) LoadLibrary (TEXT (name))
#define dlsym(handle, func) GetProcAddress (handle, func)
typedef HINSTANCE* lib_handler;
#else
#include"dlfcn.h"
typedef void* lib_handler;
#endif

#include"CalculiX.h"

static const char* search(const char *b,
			  const char *e,
			  const char  c){
  const char* p=b;
  while(p!=e){
    if(*p==c){
      return p;
    }
    ++p;
  }
  return e;
}

static const char* search_space(const char *b,
				const char *e){
  const char* c=b;
  while(c!=e){
    if(isspace(*c)){
      return c;
    }
    ++c;
  }
  return e;
}


static void report_failure(const char* n,
			   const char* e){
  char b[81] = {'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'};
  memcpy(b,n,search_space(n,n+80)-n);
  printf(" *ERROR: invalid material name '%s' (%s)\n",b,e);
  exit(-1);
}  

static void check(const int b,
		  const char* n,
		  const char* e){
  if(!b){
    report_failure(n,e);
  }
}


static void getExternalBehaviours(CalculixExternalBehaviour*** p,
				  size_t **s)
{
  static CalculixExternalBehaviour* rf = NULL;
  static size_t nrf = 0;
  *p = &rf;
  *s = &nrf;
}

static void treatExternalBehaviour(const char *n,
				   const CalculixInterface i,
				   const char *l,
				   const char *f){
#ifdef __CYGWIN__
  char b[80] = {'c','y','g','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'};
#else
  char b[80] = {'l','i','b','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
		'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'};
#endif
#ifdef CALCULIX_EXTERNAL_BEHAVIOUR_DEBUG
  fprintf(stdout,"treatExternalBehaviour: loading function '%s' from library '%s'\n",f,l);
#endif
  CalculixExternalBehaviour** rf;
  size_t *nrf;
  // full library name
  const size_t s = strlen(l);
  memcpy(b+3,l,s);
#if (! defined _WIN32) && (! defined __CYGWIN__)
  b[3+s] = '.';
  b[4+s] = 's';
  b[5+s] = 'o';
#else /* _WIN32 */
  b[3+s] = '.';
  b[4+s] = 'd';
  b[5+s] = 'l';
  b[6+s] = 'l';
#endif /* _WIN32 */
  getExternalBehaviours(&rf,&nrf);
  lib_handler lib = dlopen(b,RTLD_NOW);
#ifndef _WIN32
  check(lib!=NULL,n,dlerror());
#else
  check(lib!=NULL,n,"unable to load library");
#endif /* _WIN32 */
  void * ptr = dlsym(lib,f);
  check(ptr!=NULL,n,"unable to load function");
  if(calculix_searchExternalBehaviour(n)==NULL){
#ifdef CALCULIX_EXTERNAL_BEHAVIOUR_DEBUG
    fprintf(stdout,"treatExternalBehaviour: registring material (library '%s', function '%s')\n",l,f);
#endif
    RENEW(*rf,CalculixExternalBehaviour,*nrf+1);
    char *nn = malloc(81*sizeof(char));;
    check(nn!=NULL,n,"no memory left");
    memcpy(nn,n,80*sizeof(char));
    nn[80]='\0';
    (*rf)[*nrf]. n  = nn;
    (*rf)[*nrf]. i   = i;
    (*rf)[*nrf]. ptr = ptr;
    *nrf+=1;
#ifdef CALCULIX_EXTERNAL_BEHAVIOUR_DEBUG
  } else {
    fprintf(stdout,"searchExternalBehaviour: material is already registred\n");
#endif
  }
}

void calculix_freeExternalBehaviours()
{
  CalculixExternalBehaviour** rf;
  size_t *nrf;
  size_t i;
  getExternalBehaviours(&rf,&nrf);
  for(i=0;i!=*nrf;++i){
    free((void*) (*rf)[i].n);
  }
  free(*rf);
  *rf  = NULL;
  *nrf = 0;
}

static int compareMaterialName(const char *m1,
			       const char *m2){
  size_t i = 0;
  while(i!=80){
    if((isspace(m1[i]))||(isspace(m2[i]))){
      return m1[i]==m2[i];
    }
    if(m1[i]!=m2[i]){
      return 0;
    }
    ++i;
  }
  return 1;
}

const CalculixExternalBehaviour*
calculix_searchExternalBehaviour(const char* n){
  CalculixExternalBehaviour** rf;
  size_t *nrf;
  getExternalBehaviours(&rf,&nrf);
  size_t i = 0;
  while(i!=*nrf){
#ifdef CALCULIX_EXTERNAL_BEHAVIOUR_DEBUG
    fprintf(stdout,"searchExternalBehaviour: considering material %s\n",(*rf)[i].n);
#endif
    if(compareMaterialName((*rf)[i].n,n)==1){
      return (*rf)+i;
    }
    ++i;
  }
#ifdef CALCULIX_EXTERNAL_BEHAVIOUR_DEBUG
  fprintf(stdout,"searchExternalBehaviour: material not registred yet\n");
#endif
  return NULL;
}

void calculix_registerExternalBehaviour(const char* n) {
  // buffers to avoid memory allocation
  char b1[80] = {'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'};
  char b2[80] = {'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0',
         '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'};
  // library name
  const char * l;
  // function name
  const char * f;
  // interface name
  CalculixInterface i = CALCULIX_STANDARD_INTERFACE;
  // current position
  const char * c = n+1;
  // past-the-end position
  const char * e = search_space(n,n+80);
  const size_t s = e-n;
  const char *p1,*p2,*p3;
  check(n!=e,n,"empty string");
  if(*n!='@'){
    return;
  }
  if((s>=10)&&(strncmp(c,"ABAQUSNL_",9)==0)){
    i = CALCULIX_ABAQUSNL_INTERFACE;
    c += 9;
  } else if((s>=8)&&(strncmp(c,"ABAQUS_",7)==0)){
    i = CALCULIX_ABAQUS_INTERFACE;
    c += 7;
  }
  p1 = c;
  p2 = c = search(c,e,'_');
  if(c==e){
    p3 = search(p1,e,'@');
    l = memcpy(b1,p1,p3-p1);
    check(strlen(l)!=0,n,"empty library name");
    treatExternalBehaviour(n,i,l,"umat_");
  } else {
    // library and function
    p2 = c = search(c,e,'_');
    p3 = search(c,e,'@');
    l = memcpy(b1,p1,p2-p1);
    f = memcpy(b2,p2+1,p3-p2-1);
    check(strlen(l)!=0,n,"empty library name");
    check(strlen(f)!=0,n,"empty function name");
    treatExternalBehaviour(n,i,l,f);
  }
}

#endif /* CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT */
