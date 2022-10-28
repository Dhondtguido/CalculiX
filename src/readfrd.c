/* --------------------------------------------------------------------  */
/*                          CALCULIX                                     */
/*                   - GRAPHICAL INTERFACE -                             */
/*                                                                       */
/*     A 3-dimensional pre- and post-processor for finite elements       */
/*              Copyright (C) 1996 Klaus Wittig                          */
/*                                                                       */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; version 2 of           */
/*     the License.                                                      */
/*                                                                       */
/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */
/*                                                                       */
/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */
/* --------------------------------------------------------------------  */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "readfrd.h"

#define TEST     0

#define INI_FIELD_SIZE 1000000


/* ToDo:
*/


void freeDatasets(Datasets *lcase, int nr)
{
  register int i;

  printf(" free lc[%d] ncomps:%d\n",nr,lcase[nr].ncomps);
  if(lcase[nr].loaded)
  {
    for(i=0; i<lcase[nr].ncomps; i++) SFREE(lcase[nr].dat[i]);
  }
  /* always allocated */
  SFREE(lcase[nr].dat);
  lcase[nr].dat=NULL;

  if(lcase[nr].npheader)
  {
    for(i=0; i<lcase[nr].npheader; i++) SFREE(lcase[nr].pheader[i]);
    SFREE(lcase[nr].pheader);
    lcase[nr].pheader=NULL;
  }
  for(i=0; i<lcase[nr].ncomps; i++)
  {
    SFREE(lcase[nr].compName[i]);
    SFREE(lcase[nr].icname[i]);
  }
  /* always allocated */
  SFREE(lcase[nr].compName);
  SFREE(lcase[nr].icname);
  lcase[nr].compName=NULL;
  lcase[nr].icname=NULL;

  SFREE(lcase[nr].ictype);
  SFREE(lcase[nr].icind1);
  SFREE(lcase[nr].icind2);
  SFREE(lcase[nr].iexist);
  SFREE(lcase[nr].max);
  SFREE(lcase[nr].menu);
  SFREE(lcase[nr].min);
  SFREE(lcase[nr].nmax);
  SFREE(lcase[nr].nmin);
  SFREE(lcase[nr].fileptr);
  lcase[nr].fileptr=NULL;
  lcase[nr].loaded=0;

  /* edat not propper implemented or deleted */
  // for(i=0; i<3; i++) for(e=0; e<anz->e; e++) SFREE(lcase[nr].edat[i][e]);
}


/* read_mode=0: jump '100C' data-blocks and read them later on demand */
/* in any case for each results of a new step the first data-block is readed to see how much nodes are included */

int readfrd( char *datin, Summen *anz, Nodes **nptr, Elements **eptr, Datasets **lptr, int read_mode )
{
  FILE *handle;
  int i=0, j=0;
  int  nodeflag=0, elemflag=0, errFlag=0, firsttime=1;
  int n;  /* used in format_flag */
  long offset=0;
  fpos_t *filepntr=NULL;
  int elem_data=0,nod_data=0, nod_1st_block=0; /* nodes in resultblock, nodes in 1st block (if no "nr of nodes" are given in frd file, 100C-line) */

  int  ncomps, indx, maxcomps=0, nvals, nentities;
  char rec_str[MAX_LINE_LENGTH];
  int  node_field_size, elem_field_size;
  int  e_nmax=1, e_nmin=1;
  int  length, flag, format_flag;
  int  ipuf, nodenr=0;
  static float *value=NULL;
  static double *dvalue=NULL;

  char **dat, **compName;
  int          *menu, *ictype, *icind1, *icind2, *iexist; 

  int anz_p=-1;
  char **pheader=NULL;

  Nodes     *node=NULL;
  Elements  *elem=NULL;
  Datasets  *lcase=NULL;


  if ( (lcase = (Datasets *)malloc( 1 * sizeof(Datasets))) == NULL )
    printf("\n\n ERROR: malloc failed\n\n") ;

  anz->u=anz->n=anz->e=anz->l=-1;
  anz->nmax=0;  anz->nmin=MAX_INTEGER;
  anz->emax=0;  anz->emin=MAX_INTEGER;
  length = 1;
  format_flag=0;

  /* Open the files and check to see that it was opened correctly */
  handle = fopen (datin, "rb");
  if ( handle== NULL )  { printf ("ERROR in readfrd: The input file \"%s\" could not be opened.\n\n", datin); return(-1); }
  else  printf (" file:%s opened\n", datin);


  printf (" reading frd format\n");
  length = frecord( handle, rec_str);
  rec_str[length]='\0';
  flag = stoi(rec_str,4,5);
  if (flag == 1 )
  {
    stos(rec_str,7,12,anz->model);
    printf (" MODEL NAME:  %s", anz->model);
  }
  else
  {
    printf ("\n\nFATAL ERROR: no proper file-format found.\n\n");
    return (-1);
  }

  while(length)
  {

    /* store the beginning of the data-block for later reading */
    if(filepntr==NULL)
    {  if( (filepntr=(fpos_t *)malloc(1*sizeof(fpos_t))) == NULL ) printf(" ERROR: malloc failed\n"); }

    if(fgetpos( handle, (fpos_t *)filepntr)!=0) { printf("error in fgetpos"); return(-1); }

    read_again:;
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;
    else rec_str[length] =(char)0;
    printf ("record:%s\n", rec_str);


    flag = stoi(rec_str,1,5);
    format_flag = stoi(rec_str,74,75);
    //printf ("OPCODE:%d IFORMT:%d\n", flag, format_flag );

    if(flag == 9999) goto read_again;
    if(( (nodeflag==1)&&(flag == 2) ) || ( (elemflag==1)&&(flag == 3) ))
    {
      printf ("found a second mesh. This mesh will be ignored\n");
      if(format_flag < 2) flag=-1;
      else flag*=-1;
    }
    if(flag == 1)
    {
      /* User Header used to store general information */
      if(rec_str[5]=='U')
      {
        anz->u++;
        if(!anz->u)
        { if(( anz->uheader=(char **)malloc( sizeof(char *))) == NULL )
          printf("\n\n ERROR: malloc failed\n\n") ; }
        else if(( anz->uheader=(char **)realloc((char **)anz->uheader, (anz->u+1)*sizeof(char *))) == NULL )
          printf("\n\n ERROR: realloc failed\n\n") ;
        if(( anz->uheader[anz->u]=(char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
          printf("\n\n ERROR: malloc failed\n\n") ;
        strcpy(anz->uheader[anz->u],rec_str);
      }

      /* Project Header used to store additional Dataset information */
      if(rec_str[5]=='P')
      {
        anz_p++;
        if(!anz_p)
        { if(( pheader=(char **)malloc( sizeof(char *))) == NULL )
          printf("\n\n ERROR: malloc failed\n\n") ; }
        else if(( pheader=(char **)realloc((char **)pheader, (anz_p+1)*sizeof(char *))) == NULL )
          printf("\n\n ERROR: realloc failed\n\n") ;
        if(( pheader[anz_p]=(char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
          printf("\n\n ERROR: malloc failed\n\n") ;
        strcpy(pheader[anz_p],rec_str);
      }
    }

    else if(flag == 2)
    {
      /* store the pheaders which are leading this block */
      anz->p=anz_p+1;
      anz->pheader=pheader;
      anz_p=-1;
      pheader=NULL;

      printf ("reading Nodes\n");
      nodeflag=1;

      /* nr of nodes per block can be read from the frd file, this is not documented in the original frd-spec. */
      nod_data=stoi( rec_str, 25, 36 );
      if(nod_data>0) node_field_size=nod_data;
      else node_field_size=INI_FIELD_SIZE;
      do
      {
        if ( (node = (Nodes *)realloc( (Nodes *)node, (node_field_size+1) * sizeof(Nodes))) == NULL )
        {
          printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", node_field_size );
          node_field_size/=2;
        }
        if(node_field_size<0)
        {
          printf("\n\n ERROR: not enough memory in readfrd()\n\n");
          exit(-1);
        }
      }while(!node);
      for(i=0; i<node_field_size; i++) node[i].indx=-1;

      if (format_flag < 2)
      { 
       do
       {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        flag = stoi(rec_str,1,3);
        anz->n++;
        if (flag == -3) break;
        if (!format_flag) node[anz->n].nr = stoi(rec_str,4,8);
        else              node[anz->n].nr = stoi(rec_str,4,13);
        if (node[anz->n].nr>=node_field_size)
	{
          if(node[anz->n].nr<MAX_INTEGER/2) node_field_size=node[anz->n].nr*2+1; else node_field_size=MAX_INTEGER-2;
          nodenr=node[anz->n].nr;
          do
          {
            if ( (node = (Nodes *)realloc( (Nodes *)node, (node_field_size+1) * sizeof(Nodes))) == NULL )
            {
              printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", node_field_size );
              node_field_size=nodenr+(node_field_size-nodenr)/2;
            }
            if(node_field_size<=nodenr)
            {
              printf("\n\n ERROR: not enough memory in readfrd() for node-nr:%d available\n\n", nodenr);
              exit(-1);
            }
          }while(!node);
          for(i=anz->nmax+1; i<node_field_size; i++) node[i].indx=-1;
        }
        /* save only nodes which are not already stored */
        if(node[node[anz->n].nr].indx<0)
        {
          node[node[anz->n].nr].indx=anz->n;
          if (!format_flag)
          {
            node[node[anz->n].nr].nx = stof(&rec_str[8],1,12);
            node[node[anz->n].nr].ny = stof(&rec_str[20],1,12);
            node[node[anz->n].nr].nz = stof(&rec_str[32],1,12);
          }
          else
          {
            node[node[anz->n].nr].nx = stof(&rec_str[13],1,12);
            node[node[anz->n].nr].ny = stof(&rec_str[25],1,12);
            node[node[anz->n].nr].nz = stof(&rec_str[37],1,12);
          }
          if (node[anz->n].nr >  anz->nmax)  anz->nmax=node[anz->n].nr;
          if (node[anz->n].nr <  anz->nmin)  anz->nmin=node[anz->n].nr;
#if TEST
        printf (" n=%d x=%lf y=%lf z=%lf \n",  node[anz->n].nr,
          node[node[anz->n].nr].nx, node[node[anz->n].nr].ny,
          node[node[anz->n].nr].nz); 
#endif
        } 
       } while(flag != -3);
      }

      /* binary format */
      else
      { 
       if ( (value = (float *)realloc((float *)value, (3) * sizeof(float))) == NULL )
         printf("\n\n ERROR: realloc failed, value\n\n") ;
       for(i=0; i<nod_data; i++)
       {
        anz->n++;
        length=fread((int *)&node[anz->n].nr,sizeof(int),1,handle);
	//printf("n:%d\n", node[anz->n].nr);
        if (node[anz->n].nr>=node_field_size)
	{
          if(node[anz->n].nr<MAX_INTEGER/2) node_field_size=node[anz->n].nr*2+1; else node_field_size=MAX_INTEGER-2;
          nodenr= node[anz->n].nr;
          do
          {
            if ( (node = (Nodes *)realloc( (Nodes *)node, (node_field_size+1) * sizeof(Nodes))) == NULL )
            {
              printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", node_field_size );
              node_field_size=nodenr+(node_field_size-nodenr)/2;
            }
            if(node_field_size<=nodenr)
            {
              printf("\n\n ERROR: not enough memory in readfrd() for the node-nr:%d available\n\n", nodenr);
              exit(-1);
            }
          }while(!node);
          for(n=anz->nmax+1; n<node_field_size; n++) node[n].indx=-1;
        }
        /* save only nodes which are not already stored */
        if (format_flag == 2)
	{
          length=fread((float *)value,sizeof(float),3,handle);
	  //printf("n:%f %f %f\n", value[0],value[1],value[2]);
          if(node[node[anz->n].nr].indx<0)
          {
            node[node[anz->n].nr].indx=anz->n;
	    node[node[anz->n].nr].nx = value[0];
            node[node[anz->n].nr].ny = value[1];
            node[node[anz->n].nr].nz = value[2];
          }
          else fseeko(handle, 3*sizeof(float), SEEK_CUR);
	}
        else
	{
          if(node[node[anz->n].nr].indx<0)
          {
            length=fread((double *)&node[node[anz->n].nr].nx,sizeof(double),3,handle);
            node[node[anz->n].nr].indx=anz->n;
          }
          else fseeko(handle, 3*sizeof(double), SEEK_CUR);
	}
        if (node[anz->n].nr >  anz->nmax)  anz->nmax=node[anz->n].nr;
        if (node[anz->n].nr <  anz->nmin)  anz->nmin=node[anz->n].nr;
#if TEST
        printf (" n=%d x=%lf y=%lf z=%lf \n",  node[anz->n].nr,
          node[node[anz->n].nr].nx, node[node[anz->n].nr].ny,
          node[node[anz->n].nr].nz); 
#endif
       }
       anz->n++;
      }
      node_field_size=anz->nmax+1;
      if((node = (Nodes *)realloc( (Nodes *)node, node_field_size * sizeof(Nodes))) == NULL )
        printf("\n\n ERROR: realloc failed\n\n") ;
      else
        printf ("\n %d nodes reallocated \n",anz->nmax);
      //nod_1st_block=anz->n;
    }

    else if(flag == 3)
    {
      printf ("reading Elements\n");
      // anz->emax=-MAX_INTEGER;  anz->emin=MAX_INTEGER;
      elemflag=1;
      e_nmax=-MAX_INTEGER;  e_nmin=MAX_INTEGER;

      /* nr of elems per block can be read from the frd file, this is not documented in the original frd-spec. */
      elem_data=stoi( rec_str, 25, 36 );
      if(elem_data>0) elem_field_size=elem_data;
      else elem_field_size=INI_FIELD_SIZE;
      do
      {
        if((elem = (Elements *)realloc( (Elements *)elem, (elem_field_size+1) * sizeof(Elements))) == NULL )
        {
          printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", elem_field_size );
          elem_field_size/=2;
        }
        if(elem_field_size<0)
        {
          printf("\n\n ERROR: not enough memory in readfrd()\n\n");
          exit(-1);
        }
      }while(!elem);

      /* binary format */
      if (format_flag == 2)
      { 
        if ( (elem = (Elements *)realloc((Elements *)elem, elem_data * sizeof(Elements))) == NULL )
          printf("\n\n ERROR: in readfrd realloc failed\n\n") ;
        else
          printf ("\n %d elements allocated \n", elem_data);
        for (i=0; i<elem_data; i++)
        {
          anz->e++;
          length=fread((int *)&elem[anz->e].nr,sizeof(int),1,handle);
          length=fread((int *)&elem[anz->e].type,sizeof(int),1,handle);
          length=fread((int *)&elem[anz->e].group,sizeof(int),1,handle);
          length=fread((int *)&elem[anz->e].mat,sizeof(int),1,handle);
	  elem[anz->e].attr = 0;
          anz->etype[elem[anz->e].type]++;
          if (elem[anz->e].nr >  anz->emax)  anz->emax=elem[anz->e].nr;
          if (elem[anz->e].nr <  anz->emin)  anz->emin=elem[anz->e].nr;
          if (elem[anz->e].type == 1)      ipuf = 8;   /* HEXA8  */
          else if (elem[anz->e].type == 2) ipuf = 6;   /* PE6   */
          else if (elem[anz->e].type == 3) ipuf = 4;   /* TET4   */
          else if (elem[anz->e].type == 4) ipuf = 20;  /* HEXA20 */
          else if (elem[anz->e].type == 5) ipuf = 15;  /* PE15  */
          else if (elem[anz->e].type == 6) ipuf = 10;  /* TET10  */
          else if (elem[anz->e].type == 7) ipuf = 3;   /* TRI3   */
          else if (elem[anz->e].type == 8) ipuf = 6;   /* TRI6   */
          else if (elem[anz->e].type == 9) ipuf = 4;   /* QUAD4  */
          else if (elem[anz->e].type == 10) ipuf = 8; /* QUAD8  */
          else if (elem[anz->e].type == 11) ipuf = 2;  /* BEAM2   */
          else if (elem[anz->e].type == 12) ipuf = 3;  /* BEAM3   */
	  //printf("el:%d t:%d g:%d m:%d n:%d\n", elem[anz->e].nr, elem[anz->e].type,elem[anz->e].group,elem[anz->e].mat, ipuf);
          length=fread((int *)elem[anz->e].nod,sizeof(int),ipuf,handle);
	  //for(j=0;j<ipuf; j++) printf(" %d",elem[anz->e].nod[j]); printf("\n"); 
        }
        anz->e++;
      }
      else
      {
       do
       {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        flag = stoi(rec_str,1,3);
        anz->e++;

        if (flag == -3) break;
        else if ((flag == -1)||(flag == -2))
        {
          if (anz->e>=elem_field_size)
          {
            if(anz->e<MAX_INTEGER/2) elem_field_size=anz->e*2+1; else elem_field_size=MAX_INTEGER-2;
            do
            {
              if((elem = (Elements *)realloc( (Elements *)elem, (elem_field_size+1) * sizeof(Elements))) == NULL )
              {
                printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", elem_field_size );
                elem_field_size=anz->e+(elem_field_size-anz->e)/2;
              }
              if(elem_field_size<=anz->e)
              {
                printf("\n\n ERROR: not enough memory in readfrd()\n\n");
                exit(-1);
              }
            }while(!elem);
          }
          if (!format_flag)
          {
	    elem[anz->e].nr = stoi(&rec_str[3], 1, 5);
	    elem[anz->e].type = stoi(&rec_str[8], 1, 5);
	    elem[anz->e].group = stoi(&rec_str[13], 1, 5);
	    elem[anz->e].mat = stoi(&rec_str[18], 1, 5);
          }
          else
          {
	    elem[anz->e].nr = stoi(&rec_str[3], 1, 10);
	    elem[anz->e].type = stoi(&rec_str[13], 1, 5);
	    elem[anz->e].group = stoi(&rec_str[18], 1, 5);
	    elem[anz->e].mat = stoi(&rec_str[23], 1, 5);
          }
	  elem[anz->e].attr = 0;
          ipuf=0;
          if (elem[anz->e].nr >  anz->emax)  anz->emax=elem[anz->e].nr;
          if (elem[anz->e].nr <  anz->emin)  anz->emin=elem[anz->e].nr;
          if (elem[anz->e].type == 1)      ipuf = 8;   /* HEXA8  */
          else if (elem[anz->e].type == 2) ipuf = 6;   /* PE6   */
          else if (elem[anz->e].type == 3) ipuf = 4;   /* TET4   */
          else if (elem[anz->e].type == 4) ipuf = 20;  /* HEXA20 */
          else if (elem[anz->e].type == 5) ipuf = 15;  /* PE15  */
          else if (elem[anz->e].type == 6) ipuf = 10;  /* TET10  */
          else if (elem[anz->e].type == 7) ipuf = 3;   /* TRI3   */
          else if (elem[anz->e].type == 8) ipuf = 6;   /* TRI6   */
          else if (elem[anz->e].type == 9) ipuf = 4;   /* QUAD4  */
          else if (elem[anz->e].type == 10) ipuf = 8; /* QUAD8  */
          else if (elem[anz->e].type == 11) ipuf = 2;  /* BEAM2   */
          else if (elem[anz->e].type == 12) ipuf = 3;  /* BEAM3   */
#if TEST
          printf ("\n%d e=%d typ=%d grp=%d mat=%d \n", flag, elem[anz->e].nr,
                    elem[anz->e].type, elem[anz->e].group, elem[anz->e].mat );
#endif
          length = frecord( handle, rec_str );
          if (ipuf==0)
          {
            printf (" element:%d is from unknown type:%d\n", elem[anz->e].nr, elem[anz->e].type);
          }
          else
          {
            anz->etype[elem[anz->e].type]++;
            /* read the node-lines */
            if (!format_flag)
            {
              j=0;
              for (i=0; i<ipuf; i++)
              {
		elem[anz->e].nod[i] = stoi(&rec_str[3+j*5], 1, 5);
                if (j<14) j++;
                else
                {
                  if (i<ipuf-1) length = frecord( handle, rec_str ); j=0;
                }
              }
            }
            else
            {
              j=0;
              for (i=0; i<ipuf; i++)
              {
                elem[anz->e].nod[i] = stoi(&rec_str[3+j*10], 1, 10);
                if (j<9) j++;
                else
                {
                  if (i<ipuf-1) length = frecord( handle, rec_str ); j=0;
                }
              }
            }
          }
        }
        else
        {
          printf ("ERROR: flag:%d is not expected, must be -1 or -2!\n%s", flag, rec_str );
          exit(-1);
        }
       } while(flag != -3);
       elem_field_size=anz->e+1;
       if ( (elem = (Elements *)realloc((Elements *)elem, elem_field_size * sizeof(Elements))) == NULL )
         printf("\n\n ERROR: in readfrd realloc failed\n\n") ;
       else
         printf ("\n %d elements reallocated \n", anz->e);
      }
    }

    else if(flag == 100)
    {
      anz->l++;

      printf ("reading Dataset No:%d\n",anz->l+1);
      if ( (lcase = (Datasets *)realloc((Datasets *)lcase, (anz->l+2) * sizeof(Datasets))) == NULL )
      { printf("\n\n ERROR: malloc failure\n\n" ); exit(1); }

      lcase[anz->l].handle=(FILE *)NULL;

      /* store the pheaders which are leading this block */
      lcase[anz->l].npheader=anz_p+1;
      lcase[anz->l].pheader=pheader;
      lcase[anz->l].fileptr=NULL;
      lcase[anz->l].loaded=1;
      lcase[anz->l].format_flag=format_flag;
      anz_p=-1;
      pheader=NULL;
      offset=0;

      stos( rec_str, 7, 12, lcase[anz->l].dataset_name);
      lcase[anz->l].value=stof( rec_str, 13, 24 );

      /* nr of nodes per block can be read from the frd file, this is not documented in the original frd-spec. */
      nod_data=stoi( rec_str, 25, 36 );
      /* because of a bug in ccx2.0 this nr can be wrong. In this case it is higher than the actual nr of nodes. */
      if(nod_data>anz->n)
      {
        printf(" WARNING: in this result-block are more nodes announced:%d than in the model defined:%d\n Please inform the program-admin of the originator of the frd-file\n\n", nod_data, anz->n);
        //exit(0);
        nod_data=anz->n;
      }
#ifdef DEVEL
      if(!nod_data)
      { 
        nod_data=nod_1st_block;
        printf("nods in block assumed:%d\n",nod_data );
      }
#endif
      stos( rec_str, 37, 56, lcase[anz->l].dataset_text);
      lcase[anz->l].analysis_type=stoi( rec_str, 57, 58 );
      lcase[anz->l].step_number=stoi( rec_str, 59, 63 );
      if(strlen(rec_str)>72 )
      {
        stos(rec_str,64,73,lcase[anz->l].analysis_name);
      }
      else
      {
        strcpy(lcase[anz->l].analysis_name,"");
      }
      ncomps=nentities=indx=0;
      if (!format_flag) n=8;
      else n=13;
      errFlag=0;
      firsttime=1;

      ipuf=-1;
      //if(lcase[anz->l].analysis_type==2) //in the moment ccx writes the wrong number, therefore:
      if(lcase[anz->l].analysis_type>=2)
      {
        for(i=0;i<lcase[anz->l].npheader; i++)
        {
          if(compare(&lcase[anz->l].pheader[i][5],"PHID", 4)==4)
          {
            sscanf(lcase[anz->l].pheader[i],"%*s %d", &ipuf);
            if(ipuf>-1) sprintf(lcase[anz->l].dataset_text,"ND:%d",ipuf);
          }
	}
        if(ipuf!=-1) for(i=0;i<lcase[anz->l].npheader; i++)
        {
          if(compare(&lcase[anz->l].pheader[i][5],"PMODE", 5)==5)
          {
            sscanf(lcase[anz->l].pheader[i],"%*s %d", &ipuf);
            if(ipuf>-1) sprintf(&lcase[anz->l].dataset_text[strlen(lcase[anz->l].dataset_text)]," MODE:%d",ipuf);
          }
        }
      }

      /* initialize */
      lcase[anz->l].ncomps=0;
      do
      {
        /* bin mode, active after last column definition was read (-5 lines). Attention flag:-6 not permitted so far! */
        if (( format_flag>=2)&&(lcase[anz->l].ncomps>0)&&(ncomps==lcase[anz->l].ncomps))
        {
          //printf("format_flag=%d ncomps:%d lcncomps:%d\n",format_flag,ncomps,lcase[anz->l].ncomps);

	  /* if offset is known jump the filepointer before the next block or else continue reading assuming values for all nodes are provided */
	  if(offset)
          {
            lcase[anz->l].loaded=0;
	    if (firsttime)
	    {
              firsttime=0;

              /* store the beginning of the data-block for later reading */
              lcase[anz->l].fileptr=filepntr;
              filepntr=NULL;
              lcase[anz->l].handle=handle;
              strcpy(lcase[anz->l].filename,datin);

              if( fseeko( handle, offset, SEEK_CUR )!=0) printf("error in fseeko\n");
	    }
	  }
          else if(format_flag==2)
	  {
            if ( (value = (float *)realloc((float *)value, (lcase[anz->l].ncomps) * sizeof(float))) == NULL )
              printf("\n\n ERROR: realloc failed, value\n\n") ;
            for(n=0; n<nod_data; n++)
            {
              length=fread((int *)&nodenr,sizeof(int),1,handle);
              length=fread((float *)value,sizeof(float),lcase[anz->l].ncomps,handle);
	      // printf("n:%d N:%d ",n+1, nodenr); 
              for(i=0; i<lcase[anz->l].ncomps; i++)
              {
	        // printf(" %f",value[i]); 
                lcase[anz->l].dat[i][nodenr]= value[i];
              }
              // printf("\n");
            }
	  }
          else if(format_flag==3)
	  {
            if ( (dvalue = (double *)realloc((double *)dvalue, (lcase[anz->l].ncomps) * sizeof(double))) == NULL )
              printf("\n\n ERROR: realloc failed, dvalue\n\n") ;
            for(n=0; n<nod_data; n++)
            {
              length=fread((int *)&nodenr,sizeof(int),1,handle);
              length=fread((double *)dvalue,sizeof(double),lcase[anz->l].ncomps,handle);
	      // printf("n:%d N:%d ",n+1, nodenr); 
              for(i=0; i<lcase[anz->l].ncomps; i++)
              {
	        // printf(" %f",dvalue[i]); 
                lcase[anz->l].dat[i][nodenr]= dvalue[i];
              }
              // printf("\n");
            }
	  }
          break;
        }

        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        flag = stoi(rec_str,1,3);
	//printf("flag:%d\n", flag);
	//printf("rec in block:%s\n", rec_str);

        if(flag == -1)
        {
	  /* if offset is known jump the filepointer before the next block and continue reading until flag=-3  */
	  if(offset)
          {
            lcase[anz->l].loaded=0;
	    if (firsttime)
	    {
              firsttime=0;

              /* store the beginning of the data-block for later reading */
              lcase[anz->l].fileptr=filepntr;
              filepntr=NULL;
              lcase[anz->l].handle=handle;
              strcpy(lcase[anz->l].filename,datin);

              /* reduce the offset by the current record */
              if( fseeko( handle, offset-length, SEEK_CUR )!=0) printf("error in fseeko\n");
	    }
	  }
          else
	  {
            nod_1st_block++;
            if (format_flag) nodenr = stoi(rec_str,4,13); else nodenr = stoi(rec_str,4,8); 
            if (nodenr>anz->nmax)
            {
              if (!errFlag) { errFlag=1; printf("WARNING: found node:%d in Dataset higher than in geometry allocated:%d\n", nodenr, anz->nmax); }
            }
            else if ( lcase[anz->l].irtype == 1 )
    	    {
              if(maxcomps==6)
              {
                i=6;
                if ( format_flag)
                {
                  lcase[anz->l].dat[0][nodenr]= stof(&rec_str[  13  ], 1, 12);
                  lcase[anz->l].dat[1][nodenr]= stof(&rec_str[  25  ], 1, 12);
                  lcase[anz->l].dat[2][nodenr]= stof(&rec_str[  37  ], 1, 12);
                  lcase[anz->l].dat[3][nodenr]= stof(&rec_str[  49  ], 1, 12);
                  lcase[anz->l].dat[4][nodenr]= stof(&rec_str[  61  ], 1, 12);
                  lcase[anz->l].dat[5][nodenr]= stof(&rec_str[  73  ], 1, 12);
                }
                else
                {
                  lcase[anz->l].dat[0][nodenr]= stof(&rec_str[  8   ], 1, 12);
                  lcase[anz->l].dat[1][nodenr]= stof(&rec_str[  20  ], 1, 12);
                  lcase[anz->l].dat[2][nodenr]= stof(&rec_str[  32  ], 1, 12);
                  lcase[anz->l].dat[3][nodenr]= stof(&rec_str[  44  ], 1, 12);
                  lcase[anz->l].dat[4][nodenr]= stof(&rec_str[  56  ], 1, 12);
                  lcase[anz->l].dat[5][nodenr]= stof(&rec_str[  68  ], 1, 12);
                }
              }
              else
              {
	        for(i=0; i<maxcomps; i++) lcase[anz->l].dat[i][nodenr]= stof(&rec_str[n+i*12], 1, 12);
	      }
    	      /* printf("%d", nodenr); for (i=0; i<maxcomps; i++) printf(" %f",lcase[anz->l].dat[i][nodenr] ); printf("\n"); */
            }
            else i=0;
	  }
	}
        else if(flag == -2)
	{
          if (!format_flag) n=8;
          else n=13;
          j=0;

          /* in case the data should be read directly and not on demand */
	  if(!offset)
          {
            do
            {
              lcase[anz->l].dat[i][nodenr]= stof(&rec_str[n+j*12], 1, 12);
              i++;j++;
            }while((j<6)&&(i<lcase[anz->l].ncomps));
	  }
        }
        else if (flag == -4)
        {
          stos( rec_str, 6, 13, lcase[anz->l].name);
          lcase[anz->l].ncomps = stoi(rec_str,14,18);
          lcase[anz->l].irtype = stoi(rec_str,19,23);

          if( lcase[anz->l].irtype > 2 )
          {
            printf(" Found ELEMENT DATA, this is not suported!\n");
            anz->l--;
            goto next;
          }

          if ( (lcase[anz->l].nmax = (int *)malloc( (lcase[anz->l].ncomps) * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].nmin = (int *)malloc( (lcase[anz->l].ncomps) * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].max = (CGXFLOAT *)malloc( (lcase[anz->l].ncomps) * sizeof(CGXFLOAT))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].min = (CGXFLOAT *)malloc( (lcase[anz->l].ncomps) * sizeof(CGXFLOAT))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].compName = (char **)malloc( (lcase[anz->l].ncomps) * sizeof(char *))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].icname = (char **)malloc( (lcase[anz->l].ncomps) * sizeof(char *))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].menu = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].ictype = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].icind1 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].icind2 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].iexist = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (lcase[anz->l].dat = (CGXFLOAT **)malloc( (lcase[anz->l].ncomps) * sizeof(CGXFLOAT *))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
  printf(" gen lc[%d] ncomps:%d\n",anz->l,lcase[anz->l].ncomps);
          for(i=0; i<(lcase[anz->l].ncomps); i++)
	  {
            if ( (lcase[anz->l].compName[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
               printf("\n\n ERROR: malloc failed\n\n" );
            if ( (lcase[anz->l].icname[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
               printf("\n\n ERROR: malloc failed\n\n" );
            lcase[anz->l].max[i]=-MAX_FLOAT;
            lcase[anz->l].min[i]=MAX_FLOAT;
	  }
        }
        else if(flag == -5)
        {
          if(indx<lcase[anz->l].ncomps)
	  {
            stos(rec_str, 6, 13, lcase[anz->l].compName[indx]);
            lcase[anz->l].menu[indx] = stoi(rec_str,14,18);
            lcase[anz->l].ictype[indx] = stoi(rec_str,19,23);
            lcase[anz->l].icind1[indx] = stoi(rec_str,24,28);
            lcase[anz->l].icind2[indx] = stoi(rec_str,29,33);
            lcase[anz->l].iexist[indx] = stoi(rec_str,34,38);

            /* requests for additional components are not supported so far */
            if(lcase[anz->l].iexist[indx]==1)
	    {
              SFREE(lcase[anz->l].compName[indx]); lcase[anz->l].compName[indx]=NULL;
              SFREE(lcase[anz->l].icname[indx]); lcase[anz->l].icname[indx]=NULL;
	    }
            else ncomps++;
            indx++;
	  }
          else
	  {
	    rec_str[14]='\0';
	    printf(" WARNING: unallocated component:%d \"%s\" %d\n", ncomps, &rec_str[5],lcase[anz->l].ncomps);
	    exit(0);
	  }
          nentities++;

          /* this is the last -5 line, try to figure out an offset (length of data-block) for the file-pointer */
          /* and allocate data if no offset is defined (first time) */
          if(nentities==lcase[anz->l].ncomps)
	  {
            lcase[anz->l].ncomps=ncomps;
            if(lcase[anz->l].ncomps<6) maxcomps=lcase[anz->l].ncomps;
            else maxcomps=6;

            nvals=0;
            for (i=0; i<ncomps; i++) if(lcase[anz->l].iexist[i]!=1) nvals++;
	    printf("ncomps:%d nvals:%d\n", ncomps, nvals);

            if(!read_mode)
            {
              if(nod_data)
              {
                if (format_flag==2)
		{
                  offset= nod_data * (4+nvals*4); 
		}
                else if (format_flag==3)
		{
                  offset= nod_data * (4+nvals*8); 
		}
                else
		{
                  /* just to get an approximate offset: */
                  if (!format_flag) n=8;
                  else n=13;
                  if(nvals<=6) offset= nod_data * (n+nvals*12+1); 
                  else
                  {
                    offset=0;
                    for(i=0; i<nvals/6; i++) 
                      offset+= nod_data * (n+6*12+1);
                    if(nvals%6)
                      offset+= nod_data * (n+(nvals%6)*12+1);
      	          }
		  //printf("offset:%d nod_data:%d n:%d nvals:%d\n", offset,nod_data,n,nvals);
		}
              }
            }

            if(!offset)
            {
              /* in case the data should be read directly and not on demand */
              for(i=0; i<(lcase[anz->l].ncomps); i++)
	      {
                if ( (lcase[anz->l].dat[i] = (CGXFLOAT *)malloc( (anz->nmax+1) * sizeof(CGXFLOAT))) == NULL )
                  printf("\n\n ERROR: malloc failure\n\n" );	               
                for(j=0; j<=anz->nmax; j++) lcase[anz->l].dat[i][j]=0.;
	      }
	    }
	  }

        }
        else if(flag == -6)
        {
          length= strsplt( rec_str, ' ', &dat);
          ipuf=atoi(dat[2]);
          if ( (compName = (char **)malloc( ipuf * sizeof(char *))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          for(i=0; i<ipuf; i++)
	  {
            if ( (compName[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
               printf("\n\n ERROR: malloc failed\n\n" );
	  }
          for (i=0; i<ipuf; i++) strcpy(compName[i], lcase[anz->l].compName[ atoi(dat[i+3])-1 ]);
          for (i=0; i<ipuf; i++) strcpy(lcase[anz->l].compName[i],compName[i]);
          if ( (menu = (int *)malloc( ipuf * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (ictype = (int *)malloc( ipuf * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (icind1 = (int *)malloc( ipuf * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (icind2 = (int *)malloc( ipuf * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          if ( (iexist = (int *)malloc( ipuf * sizeof(int))) == NULL )
            printf("\n\n ERROR: malloc failure\n\n" );
          for (i=0; i<ipuf; i++) menu[i] = lcase[anz->l].menu[atoi(dat[i+3])-1];
          for (i=0; i<ipuf; i++) lcase[anz->l].menu[i] =menu[i];
          for (i=0; i<ipuf; i++) ictype[i] = lcase[anz->l].ictype[atoi(dat[i+3])-1];
          for (i=0; i<ipuf; i++) lcase[anz->l].ictype[i] =ictype[i];
          for (i=0; i<ipuf; i++) icind1[i] = lcase[anz->l].icind1[atoi(dat[i+3])-1];
          for (i=0; i<ipuf; i++) lcase[anz->l].icind1[i] =icind1[i];
          for (i=0; i<ipuf; i++) icind2[i] = lcase[anz->l].icind2[atoi(dat[i+3])-1];
          for (i=0; i<ipuf; i++) lcase[anz->l].icind2[i] =icind2[i];
          for (i=0; i<ipuf; i++) iexist[i] = lcase[anz->l].iexist[atoi(dat[i+3])-1];
          for (i=0; i<ipuf; i++) lcase[anz->l].iexist[i] =iexist[i];
	  
          for(i=0; i<length; i++) SFREE(dat[i]);
          SFREE(dat); 
	  
          for(i=0; i<ipuf; i++) SFREE(compName[i]);
          SFREE(compName); 
          SFREE(menu); 
          SFREE(ictype); 
          SFREE(icind1); 
          SFREE(icind2); 
          SFREE(iexist); 
	}
      }while(flag!=-3);

      /* in case the data should be read directly and not on demand */
      if(!offset)
      {
       for(j=0; j<anz->n; j++)
       {
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          if(lcase[anz->l].dat[i][node[j].nr] > lcase[anz->l].max[i])
          {
            lcase[anz->l].max[i]=lcase[anz->l].dat[i][node[j].nr];
            lcase[anz->l].nmax[i]=node[j].nr;
          }
          if(lcase[anz->l].dat[i][node[j].nr] < lcase[anz->l].min[i])
          {
            lcase[anz->l].min[i]=lcase[anz->l].dat[i][node[j].nr];
            lcase[anz->l].nmin[i]=node[j].nr;
          }
        }
       }
      }
    }

    else if(flag == 4)
    {
      anz->l++;
      printf ("reading Dataset No:%d\n",anz->l+1);
      if ( (lcase = (Datasets *)realloc((Datasets *)lcase, (anz->l+2) * sizeof(Datasets))) == NULL )
      { printf("\n\n ERROR: malloc failure\n\n" ); exit(1); }

      lcase[anz->l].value=stof( rec_str, 13, 25 );
      strcpy(lcase[anz->l].name,"DISP    ");
      lcase[anz->l].ncomps = 3;
      lcase[anz->l].irtype = 1;
      lcase[anz->l].loaded = 1;
      lcase[anz->l].format_flag=format_flag;

      if ( (lcase[anz->l].nmax = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].nmin = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].max = (CGXFLOAT *)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].min = (CGXFLOAT *)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].dat = (CGXFLOAT **)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].compName = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icname = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      for(i=0; i<lcase[anz->l].ncomps; i++)
      {
        if ( (lcase[anz->l].compName[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );
        if ( (lcase[anz->l].icname[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );
        lcase[anz->l].max[i]=-MAX_FLOAT;
        lcase[anz->l].min[i]=MAX_FLOAT;
        if ( (lcase[anz->l].dat[i] = (CGXFLOAT *)malloc( (anz->nmax+1) * sizeof(CGXFLOAT))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );	               
        for(j=0; j<=anz->nmax; j++) lcase[anz->l].dat[i][j]=0.;
      }
      if ( (lcase[anz->l].menu = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].ictype = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icind1 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icind2 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].iexist = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );

      for(i=0; i<lcase[anz->l].ncomps; i++)
      {
        lcase[anz->l].menu[i] = 1;
        lcase[anz->l].ictype[i] = 2;
        lcase[anz->l].icind1[i] = i+1;
        lcase[anz->l].icind2[i] = 0;
        lcase[anz->l].iexist[i] = 0;
      }

      strcpy(lcase[anz->l].compName[0], "x       ");
      strcpy(lcase[anz->l].compName[1], "y       ");
      strcpy(lcase[anz->l].compName[2], "z       ");
      errFlag=0;
      do
      {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        flag = stoi(rec_str,1,3);
        if(flag == -1)
        {
          if (!format_flag) nodenr = stoi(rec_str,4,8);
          else              nodenr = stoi(rec_str,4,13);
          if (nodenr>anz->nmax)
          {
            if (!errFlag) { errFlag=1; printf("WARNING: found node:%d in Dataset higher than in geometry allocated:%d\n", nodenr, anz->nmax); }
          }
          else if (!format_flag) 
          {
            for(i=0; i<lcase[anz->l].ncomps; i++)
              lcase[anz->l].dat[i][nodenr]= stof(&rec_str[8+i*12], 1, 12);
          }
          else  
          {
            for(i=0; i<lcase[anz->l].ncomps; i++)
              lcase[anz->l].dat[i][nodenr]= stof(&rec_str[13+i*12], 1, 12);
          }
        }
      }while(flag!=-3);
      for(n=0; n<anz->n; n++)
      {
        nodenr=node[n].nr;
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          if(lcase[anz->l].dat[i][nodenr] > lcase[anz->l].max[i])
          {
            lcase[anz->l].max[i]=lcase[anz->l].dat[i][nodenr];
            lcase[anz->l].nmax[i]=nodenr;
          }
          if(lcase[anz->l].dat[i][nodenr] < lcase[anz->l].min[i])
          {
            lcase[anz->l].min[i]=lcase[anz->l].dat[i][nodenr];
            lcase[anz->l].nmin[i]=nodenr;
          }
        }
      }
    }

    else if(flag == 5)
    {
      anz->l++;
      printf ("reading Dataset No:%d\n",anz->l+1);
      if ( (lcase = (Datasets *)realloc((Datasets *)lcase, (anz->l+2) * sizeof(Datasets))) == NULL )
      { printf("\n\n ERROR: malloc failure\n\n" ); exit(1); }

      lcase[anz->l].value=stof( rec_str, 13, 25 );
      strcpy(lcase[anz->l].name,"STRESS  ");
      lcase[anz->l].ncomps = 6;
      lcase[anz->l].irtype = 1;
      lcase[anz->l].loaded = 1;
      lcase[anz->l].format_flag=format_flag;

      if ( (lcase[anz->l].nmax = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].nmin = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].max = (CGXFLOAT *)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].min = (CGXFLOAT *)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].dat = (CGXFLOAT **)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].compName = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icname = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      for(i=0; i<lcase[anz->l].ncomps; i++)
      {
        if ( (lcase[anz->l].dat[i] = (CGXFLOAT *)malloc( (anz->nmax+1) * sizeof(CGXFLOAT))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );	               
        if ( (lcase[anz->l].compName[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );
        if ( (lcase[anz->l].icname[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );
        lcase[anz->l].max[i]=-MAX_FLOAT;
        lcase[anz->l].min[i]=MAX_FLOAT;
        for(j=0; j<=anz->nmax; j++) lcase[anz->l].dat[i][j]=0.;
      }
      if ( (lcase[anz->l].menu = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].ictype = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icind1 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icind2 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].iexist = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );

      for(i=0; i<lcase[anz->l].ncomps; i++)
      {
        lcase[anz->l].menu[i] = 1;
        lcase[anz->l].ictype[i] = 4;
        lcase[anz->l].iexist[i] = 0;
      }
      lcase[anz->l].icind1[0] = 1;
      lcase[anz->l].icind2[0] = 1;
      lcase[anz->l].icind1[1] = 2;
      lcase[anz->l].icind2[1] = 2;
      lcase[anz->l].icind1[2] = 3;
      lcase[anz->l].icind2[2] = 3;
      lcase[anz->l].icind1[3] = 1;
      lcase[anz->l].icind2[3] = 2;
      lcase[anz->l].icind1[4] = 2;
      lcase[anz->l].icind2[4] = 3;
      lcase[anz->l].icind1[5] = 3;
      lcase[anz->l].icind2[5] = 1;

      strcpy(lcase[anz->l].compName[0], "xx      ");
      strcpy(lcase[anz->l].compName[1], "yy      ");
      strcpy(lcase[anz->l].compName[2], "zz      ");
      strcpy(lcase[anz->l].compName[3], "xy      ");
      strcpy(lcase[anz->l].compName[4], "yz      ");
      strcpy(lcase[anz->l].compName[5], "zx      ");
      errFlag=0;
      do
      {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        flag = stoi(rec_str,1,3);
        if(flag == -1)
        {
          if (!format_flag) nodenr = stoi(rec_str,4,8);
          else              nodenr = stoi(rec_str,4,13);
          if (nodenr>anz->nmax)
          {
            if (!errFlag) { errFlag=1; printf("WARNING: found node:%d in Dataset higher than in geometry allocated:%d\n", nodenr, anz->nmax); }
          }
          else
          {
            /* new line */
            length = frecord( handle, rec_str);
            if (rec_str[length] == (char)EOF) break;
            if (!format_flag) 
            {
              for(i=0; i<lcase[anz->l].ncomps; i++)
                lcase[anz->l].dat[i][nodenr]= stof(&rec_str[8+i*12], 1, 12);
            }
            else  
            {
              for(i=0; i<lcase[anz->l].ncomps; i++)
                lcase[anz->l].dat[i][nodenr]= stof(&rec_str[13+i*12], 1, 12);
            }
          }
        }
      }while(flag!=-3);
      for(n=0; n<anz->n; n++)
      {
        nodenr=node[n].nr;
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          if(lcase[anz->l].dat[i][nodenr] > lcase[anz->l].max[i])
          {
            lcase[anz->l].max[i]=lcase[anz->l].dat[i][nodenr];
            lcase[anz->l].nmax[i]=nodenr;
          }
          if(lcase[anz->l].dat[i][nodenr] < lcase[anz->l].min[i])
          {
            lcase[anz->l].min[i]=lcase[anz->l].dat[i][nodenr];
            lcase[anz->l].nmin[i]=nodenr;
          }
        }
      }
    }

    else if((flag == 7)||(flag == 9))
    {
      anz->l++;
      printf ("reading Dataset No:%d\n",anz->l+1);
      if ( (lcase = (Datasets *)realloc((Datasets *)lcase, (anz->l+2) * sizeof(Datasets))) == NULL )
      { printf("\n\n ERROR: malloc failure\n\n" ); exit(1); }

      lcase[anz->l].value=stof( rec_str, 13, 25 );
      strcpy(lcase[anz->l].name,"TEMP    ");
      lcase[anz->l].ncomps = 1;
      lcase[anz->l].irtype = 1;
      lcase[anz->l].loaded = 1;
      lcase[anz->l].format_flag=format_flag;

      if ( (lcase[anz->l].nmax = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].nmin = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].max = (CGXFLOAT *)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].min = (CGXFLOAT *)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].dat = (CGXFLOAT **)malloc( lcase[anz->l].ncomps * sizeof(CGXFLOAT *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].compName = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icname = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      for(i=0; i<lcase[anz->l].ncomps; i++)
      {
        if ( (lcase[anz->l].dat[i] = (CGXFLOAT *)malloc( (anz->nmax+1) * sizeof(CGXFLOAT))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );	               
        if ( (lcase[anz->l].compName[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );
        if ( (lcase[anz->l].icname[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );
        lcase[anz->l].max[i]=-MAX_FLOAT;
        lcase[anz->l].min[i]=MAX_FLOAT;
        for(j=0; j<=anz->nmax; j++) lcase[anz->l].dat[i][j]=0.;
      }
      if ( (lcase[anz->l].menu = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].ictype = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icind1 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].icind2 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );
      if ( (lcase[anz->l].iexist = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
        printf("\n\n ERROR: malloc failure\n\n" );

      for(i=0; i<lcase[anz->l].ncomps; i++)
      {
        lcase[anz->l].menu[i] = 1;
        lcase[anz->l].ictype[i] = 1;
        lcase[anz->l].icind1[i] = i+1;
        lcase[anz->l].icind2[i] = 0;
        lcase[anz->l].iexist[i] = 0;
      }

      strcpy(lcase[anz->l].compName[0], "Value   ");
      errFlag=0;
      do
      {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        flag = stoi(rec_str,1,3);
        if(flag == -1)
        {
          if (!format_flag) nodenr = stoi(rec_str,4,8);
          else              nodenr = stoi(rec_str,4,13);
          if (nodenr>anz->nmax)
          {
            if (!errFlag) { errFlag=1; printf("WARNING: found node:%d in Dataset higher than in geometry allocated:%d\n", nodenr, anz->nmax); }
          }
          else if (!format_flag) 
          {
            for(i=0; i<lcase[anz->l].ncomps; i++)
              lcase[anz->l].dat[i][nodenr]= stof(&rec_str[8+i*12], 1, 12);
          }
          else  
          {
            for(i=0; i<lcase[anz->l].ncomps; i++)
              lcase[anz->l].dat[i][nodenr]= stof(&rec_str[13+i*12], 1, 12);
          }
        }
      }while(flag!=-3);
      for(n=0; n<anz->n; n++)
      {
        nodenr=node[n].nr;
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          if(lcase[anz->l].dat[i][nodenr] > lcase[anz->l].max[i])
          {
            lcase[anz->l].max[i]=lcase[anz->l].dat[i][nodenr];
            lcase[anz->l].nmax[i]=nodenr;
          }
          if(lcase[anz->l].dat[i][nodenr] < lcase[anz->l].min[i])
          {
            lcase[anz->l].min[i]=lcase[anz->l].dat[i][nodenr];
            lcase[anz->l].nmin[i]=nodenr;
          }
        }
      }
    }

    else if(flag == -1)
    {
  next:;
      printf (" overread Block: %s\n", rec_str);
      do
      {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        /* printf ("\n record:%d %s\n", length, rec_str);   */
        if (length != 0)
        {
          flag = stoi(rec_str,1,5);
        }
      } while(flag != -3);
    }
    else if(flag == -2)
    {
      /* nr of nodes per block can be read from the frd file, this is not documented in the original frd-spec. */
      nod_data=stoi( rec_str, 25, 36 );
      printf ("scip %d nodes\n",nod_data);
      if(format_flag == 2) fseeko(handle, nod_data*(sizeof(int)+3*sizeof(float)), SEEK_CUR);
      else  fseeko(handle, nod_data*(sizeof(int)+3*sizeof(double)), SEEK_CUR);
    }
    else if(flag == -3)
    {
      /* nr of elems per block can be read from the frd file, this is not documented in the original frd-spec. */
      elem_data=stoi( rec_str, 25, 36 );
      printf ("scip %d elements\n",elem_data);
      for (i=0; i<elem_data; i++)
      {
        fseeko(handle, sizeof(int), SEEK_CUR);
        length=fread((int *)&j,sizeof(int),1,handle);
        fseeko(handle, 2*sizeof(int), SEEK_CUR);
        if (j == 1)      ipuf = 8;   /* HEXA8  */
        else if (j == 2) ipuf = 6;   /* PE6   */
        else if (j == 3) ipuf = 4;   /* TET4   */
        else if (j == 4) ipuf = 20;  /* HEXA20 */
        else if (j == 5) ipuf = 15;  /* PE15  */
        else if (j == 6) ipuf = 10;  /* TET10  */
        else if (j == 7) ipuf = 3;   /* TRI3   */
        else if (j == 8) ipuf = 6;   /* TRI6   */
        else if (j == 9) ipuf = 4;   /* QUAD4  */
        else if (j == 10) ipuf = 8; /* QUAD8  */
        else if (j == 11) ipuf = 2;  /* BEAM2   */
        else if (j == 12) ipuf = 3;  /* BEAM3   */
        fseeko(handle, ipuf*sizeof(int), SEEK_CUR);
      }
    }
    else
    {
      printf (" overread Block: %s\n", rec_str);
      do
      {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        /* printf ("\n record:%d %s\n", length, rec_str);   */
        if (length != 0)
        {
          flag = stoi(rec_str,1,5);
        }
      } while(flag != -3);
    }
  }

  for(i=0; i<anz_p; i++) SFREE(pheader[i]); SFREE(pheader);
  SFREE( filepntr);

  anz->u++;
  anz->l++;
  if(anz->n<0) anz->n=0;
  if(anz->e<0) anz->e=0;

  if ( e_nmax > (anz->nmax) )
  {
    printf ("\nWARNING: element requestes a nodename higher than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }
  if ( e_nmin < 1 )
  {
    printf ("\nWARNING: element requestes a nodename lower than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }

  anz->orign    = anz->n;
  anz->orignmax = anz->nmax;
  anz->olc = anz->l;
  
  *nptr = node; *eptr = elem; *lptr = lcase;
  return(1);
}



/* if a block was skipped during first read of frd-file read it now */
/* return -1 if failure */
int readfrdblock(int lc, Summen *anz,   Nodes     *node, Datasets *lcase )
{
  register int i,j, n;
  int  length, flag, format_flag, nodenr=0, ncomps;
  int  nod_data=0, maxcomps=0;
  int  errFlag=0;
  char rec_str[MAX_LINE_LENGTH];
  static float *value=NULL;
  static double *dvalue=NULL;
  FILE *handle;

  // printf("readfrdblock for file:%s\n", lcase[lc].filename);

  /* Open the files and check to see that it was opened correctly */
  handle = lcase[lc].handle;
  if ( handle== NULL )  { printf ("ERROR in readfrdblock: The input file \"%s\" could not be opened.\n\n", lcase[lc].filename); return(-1); }

  if( fsetpos( handle, (fpos_t *)lcase[lc].fileptr)!=0) { printf("error in fsetpos"); return(-1); }
  lcase[lc].loaded=1;
  
  length = frecord( handle, rec_str);
  if (rec_str[length] == (char)EOF) return(-1);

  flag = stoi(rec_str,1,5);
  format_flag = stoi(rec_str,74,75);

  if(lcase[lc].ncomps<6) maxcomps=lcase[lc].ncomps;
  else maxcomps=6;

  if( lcase[lc].irtype > 2 )
  {
    printf(" ERROR: Found ELEMENT DATA, this is not suported!\n");
    return(-1);
  }

  if ( (lcase[lc].dat = (CGXFLOAT **)malloc( (lcase[lc].ncomps) * sizeof(CGXFLOAT *))) == NULL )
    printf("\n\n ERROR: malloc failure\n\n" );
  for(i=0; i<(lcase[lc].ncomps); i++)
  {
    if ( (lcase[lc].dat[i] = (CGXFLOAT *)malloc( (anz->nmax+1) * sizeof(CGXFLOAT))) == NULL )
      printf("\n\n ERROR: malloc failure\n\n" );	               
    for(j=0; j<=anz->nmax; j++) lcase[lc].dat[i][j]=0.;
  }


  if(flag == 100)
  {
    nod_data=stoi( rec_str, 25, 36 );
    if (!format_flag) n=8;
    else n=13;
    do
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
      flag = stoi(rec_str,1,3);

      /* bin mode, active after last column definition was read (-5 lines). Attention flag:-6 not permitted so far! */
      if ( format_flag>=2)
      {
        if (flag == -4)
        {
          ncomps = stoi(rec_str,14,18);
        }
        else continue;

        //printf("format_flag=%d ncomps:%d\n",format_flag,ncomps);

        /* skip the meta-data */
        for(i=0; i<ncomps; i++) length = frecord( handle, rec_str);

        if ( format_flag==2)
        {
          if ( (value = (float *)realloc((float *)value, (lcase[lc].ncomps) * sizeof(float))) == NULL )
            printf("\n\n ERROR: realloc failed, value\n\n") ;
          for(n=0; n<nod_data; n++)
          {
            length=fread((int *)&nodenr,sizeof(int),1,handle);
            length=fread((float *)value,sizeof(float),lcase[lc].ncomps,handle);
  	  // printf("N:%d ",nodenr);
            for(i=0; i<lcase[lc].ncomps; i++)
            {
  	    // printf(" %f",value[i]); 
              lcase[lc].dat[i][nodenr]= value[i];
            }
  	  // printf("\n");
          }
        }
        if ( format_flag==3)
        {
          if ( (dvalue = (double *)realloc((double *)dvalue, (lcase[lc].ncomps) * sizeof(double))) == NULL )
            printf("\n\n ERROR: realloc failed, dvalue\n\n") ;
          for(n=0; n<nod_data; n++)
          {
            length=fread((int *)&nodenr,sizeof(int),1,handle);
            length=fread((double *)dvalue,sizeof(double),lcase[lc].ncomps,handle);
  	  // printf("N:%d ",nodenr);
            for(i=0; i<lcase[lc].ncomps; i++)
            {
  	    // printf(" %f",dvalue[i]); 
              lcase[lc].dat[i][nodenr]= dvalue[i];
            }
  	  // printf("\n");
          }
        }
        break;        
      }

      else if(flag == -1)
      {
        if (format_flag) nodenr = stoi(rec_str,4,13); else nodenr = stoi(rec_str,4,8); 
        if (nodenr>anz->nmax)
        {
          if (!errFlag) { errFlag=1; printf("WARNING: found node:%d in Dataset higher than in geometry allocated:%d\n", nodenr, anz->nmax); }
        }
        else if ( lcase[lc].irtype == 1 )
	{
          if(maxcomps==6)
          {
            i=6;
            if ( format_flag)
            {
              lcase[lc].dat[0][nodenr]= stof(&rec_str[  13  ], 1, 12);
              lcase[lc].dat[1][nodenr]= stof(&rec_str[  25  ], 1, 12);
              lcase[lc].dat[2][nodenr]= stof(&rec_str[  37  ], 1, 12);
              lcase[lc].dat[3][nodenr]= stof(&rec_str[  49  ], 1, 12);
              lcase[lc].dat[4][nodenr]= stof(&rec_str[  61  ], 1, 12);
              lcase[lc].dat[5][nodenr]= stof(&rec_str[  73  ], 1, 12);
            }
            else
            {
              lcase[lc].dat[0][nodenr]= stof(&rec_str[  8   ], 1, 12);
              lcase[lc].dat[1][nodenr]= stof(&rec_str[  20  ], 1, 12);
              lcase[lc].dat[2][nodenr]= stof(&rec_str[  32  ], 1, 12);
              lcase[lc].dat[3][nodenr]= stof(&rec_str[  44  ], 1, 12);
              lcase[lc].dat[4][nodenr]= stof(&rec_str[  56  ], 1, 12);
              lcase[lc].dat[5][nodenr]= stof(&rec_str[  68  ], 1, 12);
            }
          }
          else for(i=0; i<maxcomps; i++) lcase[lc].dat[i][nodenr]= stof(&rec_str[n+i*12], 1, 12);
        }
        else i=0;
      }
      else if(flag == -2)
      {
        if (!format_flag) n=8;
        else n=13;
        j=0;
        do
        {
          lcase[lc].dat[i][nodenr]= stof(&rec_str[n+j*12], 1, 12);
          i++;j++;
        }while((j<6)&&(i<lcase[lc].ncomps));
      }
    }while(flag!=-3);

    for(j=0; j<anz->orign; j++)
    {
      for(i=0; i<lcase[lc].ncomps; i++)
      {
        if(lcase[lc].dat[i][node[j].nr] > lcase[lc].max[i])
        {
          lcase[lc].max[i]=lcase[lc].dat[i][node[j].nr];
          lcase[lc].nmax[i]=node[j].nr;
        }
        if(lcase[lc].dat[i][node[j].nr] < lcase[lc].min[i])
        {
          lcase[lc].min[i]=lcase[lc].dat[i][node[j].nr];
          lcase[lc].nmin[i]=node[j].nr;
        }
      }
    }
  }

  return(0);
}


/* regula falsi to find the matching record fast */
/* not finished */
char *getRecord(FILE *handle, int n, int x0 )
{
    int ii, m, n1,n2, x=0, offset=0;

    /* search the intersection */
    n1=0;                              
    n2=n;                            
    for(ii=0; ii<n; ii++)
    {                     
      m=(n2+n1)/2;

      
      if( fseeko( handle, offset, SEEK_CUR )!=0) printf("error in fseeko\n");
                          
      if(x0>= x ) n1=m;              
      if(x0 < x ) n2=m;              
      if((n2-n1) == 1) break;           
    }                                 
#if TEST
    printf("x:%d x0:%d\n", x,x0); 
#endif
  return(NULL);
}


/* return -1 if failure */
/* return 0 if successfull */
int readOneNode( int lc, Summen *anz, Datasets *lcase, int nodenr, double **vptr, long *byte_offset )
{
  register int i=0,j, n;
  int length, flag, inodenr=0;
  int   maxcomps=0, ncomps;
  long offset=0;
  int  errFlag=0, readFlag=0, bailout=0;
  char rec_str[MAX_LINE_LENGTH];
  FILE *handle;
  static float *value=NULL;
  static double *dvalue=NULL;
  double *dat;
  int  nod_data=0, nvals=0;

  // printf("readOneNode for file:%s *byte_offset:%d\n", lcase[lc].filename, *byte_offset);

  /* Open the files and check to see that it was opened correctly */
  handle = lcase[lc].handle;
  if ( handle== NULL )  { printf ("ERROR in readOneNode: The input file \"%s\" could not be opened.\n\n", lcase[lc].filename); return(-1); }

  if( fsetpos( handle, (fpos_t *)lcase[lc].fileptr)!=0) { printf("error in fsetpos"); return(-1); }

  /* the header-lines must be skipped in case byte_offset==0 and format == 2(bin) */
  if((lcase[lc].format_flag>=2)&&(!*byte_offset))
  {
    do
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) return(-1);
      *byte_offset+=length+1;
      printf ("record:%s\n", rec_str);
      flag = stoi(rec_str,1,5);
      if (flag == -4)
      {
        ncomps = stoi(rec_str,14,18);
      }
      else continue;

      /* skip the meta-data */
      for(i=0; i<ncomps; i++) *byte_offset+= frecord( handle, rec_str)+1;
      break;        
    }while(1);
  }
  else{ if( fseeko( handle, *byte_offset, SEEK_CUR )!=0) printf("error in fseeko\n"); }

  offset=*byte_offset;

  if ( (dat = (double *)malloc( (lcase[lc].ncomps) * sizeof(double ))) == NULL )
    printf("\n\n ERROR: malloc failure\n\n" );
  *vptr=dat;

  /* if bin mode */
  if(lcase[lc].format_flag==2)
  {
    if ( (value = (float *)realloc((float *)value, (lcase[lc].ncomps) * sizeof(float))) == NULL )
      printf("\n\n ERROR: realloc failed, value\n\n") ;
    do
    {
      length=fread((int *)&inodenr,sizeof(int),1,handle)*sizeof(int);
      length+=fread((float *)value,sizeof(float),lcase[lc].ncomps,handle)*sizeof(float);
      //printf("N:%d ",inodenr);
      //for(i=0; i<lcase[lc].ncomps; i++) printf(" %f",value[i]);
      //printf("\n");
      if(inodenr==nodenr) break;
      else offset+=length;
    }while(1);
    for(i=0; i<lcase[lc].ncomps; i++)
    {
	//printf(" %f",value[i]);
      dat[i]= value[i];
    }

    *byte_offset=offset;
    //printf("offset:%d\n", offset);
    return(0);
  }
  if(lcase[lc].format_flag==3)
  {
    if ( (dvalue = (double *)realloc((double *)dvalue, (lcase[lc].ncomps) * sizeof(double))) == NULL )
      printf("\n\n ERROR: realloc failed, dvalue\n\n") ;
    do
    {
      length=fread((int *)&inodenr,sizeof(int),1,handle)*sizeof(int);
      length+=fread((double *)dvalue,sizeof(double),lcase[lc].ncomps,handle)*sizeof(double);
      //printf("N:%d ",inodenr);
      //for(i=0; i<lcase[lc].ncomps; i++) printf(" %f",dvalue[i]);
      //printf("\n");
      if(inodenr==nodenr) break;
      else offset+=length;
    }while(1);
    for(i=0; i<lcase[lc].ncomps; i++)
    {
	//printf(" %f",dvalue[i]);
      dat[i]= dvalue[i];
    }

    *byte_offset=offset;
    //printf("offset:%d\n", offset);
    return(0);
  }

  if(lcase[lc].ncomps<6) maxcomps=lcase[lc].ncomps;
  else maxcomps=6;
  if (!lcase[lc].format_flag) n=8;
  else n=13;

  if( lcase[lc].irtype > 2 )
  {
    printf(" ERROR: Found ELEMENT DATA, this is not suported!\n");
    return(-1);
  }

 repeat_search:;
  do
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF)
    {
      /* offset obviously false */
      offset=0;
      if( fsetpos( handle, (fpos_t *)lcase[lc].fileptr)!=0) { printf("error in fsetpos"); return(-1); }
    }
    flag = stoi(rec_str,1,3);

    if((readFlag)&&(flag == -2))
    {
      //printf("-2 node found:%d %d at pos:%d\n", inodenr,nodenr, nod_data);
      if (!lcase[lc].format_flag) n=8;
      else n=13;
      j=0;
      do
      {
        dat[i]= stof(&rec_str[n+j*12], 1, 12);
        i++;j++;
      }while((j<6)&&(i<lcase[lc].ncomps));
    }
    else if(readFlag)
    {
      /* leave after all records of that node are read */

      nod_data--;

      /* get the offset for the next call: */
      for (i=0; i<lcase[lc].ncomps; i++) if(lcase[lc].iexist[i]!=1) nvals++;

      if (!lcase[lc].format_flag) n=8;
      else n=13;
      if(nvals<=6) offset+= nod_data * (n+nvals*12); 
      else
      {
        for(i=0; i<nvals/6; i++) 
          offset+= nod_data * (n+6*12+1);
        if(nvals%6)
          offset+= nod_data * (n+(nvals%6)*12+1);
      }
      *byte_offset=offset;
      //printf("offset:%d nod_data:%d n:%d nvals:%d dat:%f\n", offset,nod_data,n,nvals, dat[0]);
      return(0);
    }
    else if(flag == -1)
    {
      if (lcase[lc].format_flag) inodenr = stoi(rec_str,4,13); else inodenr = stoi(rec_str,4,8);
      nod_data++; 
      //printf("node:%d %d at pos:%d ir:%d\n", inodenr,nodenr, nod_data, lcase[lc].irtype);
      if (inodenr>anz->nmax)
      {
        if (!errFlag) { errFlag=1; printf("WARNING: found node:%d in Dataset higher than in geometry allocated:%d\n", inodenr, anz->nmax); }
      }
      else if (( inodenr==nodenr)&&(lcase[lc].irtype == 1 ))
      {
        //printf("node found:%d %d at pos:%d\n", inodenr,nodenr, nod_data);
        readFlag=1;
        if(maxcomps==6)
        {
          i=6;
          if ( lcase[lc].format_flag)
          {
            dat[0]= stof(&rec_str[  13  ], 1, 12);
            dat[1]= stof(&rec_str[  25  ], 1, 12);
            dat[2]= stof(&rec_str[  37  ], 1, 12);
            dat[3]= stof(&rec_str[  49  ], 1, 12);
            dat[4]= stof(&rec_str[  61  ], 1, 12);
            dat[5]= stof(&rec_str[  73  ], 1, 12);
          }
          else
          {
            dat[0]= stof(&rec_str[  8   ], 1, 12);
            dat[1]= stof(&rec_str[  20  ], 1, 12);
            dat[2]= stof(&rec_str[  32  ], 1, 12);
            dat[3]= stof(&rec_str[  44  ], 1, 12);
            dat[4]= stof(&rec_str[  56  ], 1, 12);
            dat[5]= stof(&rec_str[  68  ], 1, 12);
          }
        }
        else for(i=0; i<maxcomps; i++) dat[i]= stof(&rec_str[n+i*12], 1, 12);
      }
      else i=0;
    }
    
  }while(flag!=-3);

  if(!bailout)
  {
    bailout=1;
    offset=0;
    //printf("repeat search\n");
    if( fsetpos( handle, (fpos_t *)lcase[lc].fileptr)!=0) { printf("error in fsetpos"); return(-1); }
    goto repeat_search;
  }
  return(-1);
}

