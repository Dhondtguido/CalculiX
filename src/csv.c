/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2023 Guido Dhondt                          */

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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

#define min(a, b) ((a) <= (b) ? (a) : (b))
#define max(a, b) ((a) >= (b) ? (a) : (b))
char *get_set(char *s)
{
  int l = strlen(s);

  while (isspace(s[l - 1]))
    --l;
  while (*s && isspace(*s))
    ++s, --l;

  return strndup(s, l - 1);
}
void csv(double *co, ITG *nk, ITG *kon, ITG *ipkon, char *lakon, ITG *ne0,
         double *v, double *stn, ITG *inum, ITG *nmethod, ITG *kode,
         char *filab, double *een, double *t1, double *fn, double *time,
         double *epn, ITG *ielmat, char *matname, double *enern,
         double *xstaten, ITG *nstate_, ITG *istep, ITG *iinc,
         ITG *ithermal, double *qfn, ITG *mode, ITG *noddiam,
         double *trab, ITG *inotr, ITG *ntrans, double *orab,
         ITG *ielorien, ITG *norien, char *description, ITG *ipneigh,
         ITG *neigh, ITG *mi, double *stx, double *vr, double *vi,
         double *stnr, double *stni, double *vmax, double *stnmax,
         ITG *ngraph, double *veold, double *ener, ITG *ne, double *cs,
         char *set, ITG *nset, ITG *istartset, ITG *iendset, ITG *ialset,
         double *eenmax, double *fnr, double *fni, double *emn,
         double *thicke, char *jobnamec, char *output, double *qfx,
         double *cdn, ITG *mortar, double *cdnr, double *cdni, ITG *nmat,
         ITG *ielprop, double *prop, double *sti)
{

  /* stores the results in frd format

     iselect selects which nodes are to be stored:
     iselect=-1 means only those nodes for which inum negative
     ist, i.e. network nodes
     iselect=+1 means only those nodes for which inum positive
     ist, i.e. structural nodes
     iselect=0  means both of the above */

  char text[2] = " ";

  char nodes_filename[30],
      time_step[12],
      set_name[64],
      csv_filename[64],
      csv_header[240];
  
  FILE *nodes_file;

  static ITG icounter = 0, nkcoords, iaxial;

  ITG i, j, k, noutloc, iset, iselect, ncomp, nope, nodes, ifield[7],
      nfield[2], icomp[7], ifieldstate[*nstate_], icompstate[*nstate_],
      ioutall = 0, *inumshell = NULL, nterms, nout, noutplus, noutmin,
      mt = mi[1] + 1, iflag, numfield, iorienloc;

  ITG ncompscalar = 1, ifieldscalar[1] = {1}, icompscalar[1] = {0},
      nfieldscalar[2] = {1, 0};
  ITG ncompvector = 3, ifieldvector[3] = {1, 1, 1}, icompvector[3] = {0, 1, 2},
      nfieldvector1[2] = {3, 0}, nfieldvector0[2] = {mi[1] + 1, 0},
      icompvectorlast[3] = {3, 4, 5};
  ITG ncomptensor = 6, ifieldtensor[6] = {1, 1, 1, 1, 1, 1}, icomptensor[6] = {0, 1, 2, 3, 5, 4},
      nfieldtensor[2] = {6, 0};
  ITG ncompscalph = 2, ifieldscalph[2] = {1, 2}, icompscalph[2] = {0, 0},
      nfieldscalph[2] = {0, 0};
  ITG ncompvectph = 6, ifieldvectph[6] = {1, 1, 1, 2, 2, 2}, icompvectph[6] = {1, 2, 3, 1, 2, 3},
      nfieldvectph[2] = {mi[1] + 1, mi[1] + 1};
  ITG ncomptensph = 12, ifieldtensph[12] = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2},
      icomptensph[12] = {0, 1, 2, 3, 5, 4, 0, 1, 2, 3, 5, 4}, nfieldtensph[2] = {6, 6};

  double *errn = NULL, *ethn = NULL;

  /* check whether all results have to be stored (also those
     corresponding to inactive nodes or elements) */

  if (strcmp1(&output[3], "a") == 0)
    ioutall = 1;

  strncpy(set_name, filab + 6, 30);
  sprintf(set_name, get_set(set_name));

  /* determining nout, noutplus and noutmin
     nout: number of structural and network nodes
     noutplus: number of structural nodes
     noutmin: number of network nodes */

  if (*nmethod != 0)
  {
    nout = 0;
    noutplus = 0;
    noutmin = 0;
    if (ioutall == 0)
    {
      for (i = 0; i < *nk; i++)
      {
        if (inum[i] == 0)
          continue;
        nout++;
        if (inum[i] > 0)
          noutplus++;
        if (inum[i] < 0)
          noutmin++;
      }
    }
    else
    {
      for (i = 0; i < *nk; i++)
      {
        nout++;
        if (inum[i] >= 0)
          noutplus++;
        if (inum[i] <= 0)
          noutmin++;
      }
    }
  }
  else
  {
    nout = *nk;
  }

  ITG ncomma;

  if (*time <= 0.)
  {
    sprintf(time_step, "%12.5E", *time);
  }
  else if ((log10(*time) >= 0) && (log10(*time) < 10.))
  {
    ncomma = 10 - floor(log10(*time) + 1.);
    if (ncomma == 0)
    {
      sprintf(time_step, "%12.0f", *time);
    }
    else if (ncomma == 1)
    {
      sprintf(time_step, "%12.1f", *time);
    }
    else if (ncomma == 2)
    {
      sprintf(time_step, "%12.2f", *time);
    }
    else if (ncomma == 3)
    {
      sprintf(time_step, "%12.3f", *time);
    }
    else if (ncomma == 4)
    {
      sprintf(time_step, "%12.4f", *time);
    }
    else if (ncomma == 5)
    {
      sprintf(time_step, "%12.5f", *time);
    }
    else if (ncomma == 6)
    {
      sprintf(time_step, "%12.6f", *time);
    }
    else if (ncomma == 7)
    {
      sprintf(time_step, "%12.7f", *time);
    }
    else if (ncomma == 8)
    {
      sprintf(time_step, "%12.8f", *time);
    }
    else
    {
      sprintf(time_step, "%12.9f", *time);
    }
  }
  else
  {
    sprintf(time_step, "%12.5E", *time);
  }

  /* store the nodal coordinates */

  if ((*kode == 1) && ((*nmethod != 5) || (*mode != 0)))
  {
    iaxial = 0.;
    sprintf(nodes_filename, "nodes_%s.csv", jobnamec);
    nodes_file = fopen(nodes_filename, "w");
    if (nodes_file == NULL)
    {
      fprintf(stderr, "Couldn't open %s\n", nodes_filename);
      exit(1);
    }

    /* storing the coordinates themselves */

    if (*nmethod != 0)
    {
      fprintf(nodes_file, "# nodal coordinates (node,x,y,z) for %s\n"
                          "node,x,y,z\n",
              jobnamec);
      for (i = 0; i < *nk; i++)
      {
        if ((inum[i] == 0) && (ioutall != 1))
          continue;
        fprintf(nodes_file, "%" ITGFORMAT ",%12.5E,%12.5E,%12.5E\n", i + 1, (float)co[3 * i],
                (float)co[3 * i + 1], (float)co[3 * i + 2]);
      }
    }
    fclose(nodes_file);
  }
  /* nkcoords is the number of nodes at the time when
   the nodal coordinates are stored in the csv file */

  nkcoords = *nk;

  /*  for cyclic symmetry frequency calculations only results for
      even numbers (= odd modes, numbering starts at 0) are stored */

  if ((*nmethod == 2) && (((*mode / 2) * 2 != *mode) && (*noddiam >= 0)))
    return;

  /* storing the displacements in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(filab, "U ") == 0) && (*ithermal != 2))
    {
      iselect = 1;
      frdset(filab, set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      if (mi[1] == 3)
      {
        sprintf(csv_header, "# displacements as (D1,D2,D3) for set %s and time %s\n"
                            "node,D1,D2,D3\n",
                set_name, time_step);
        csvvector(v, &iset, ntrans, filab, &nkcoords, inum, inotr,
                  trab, co, istartset, iendset, ialset, mi, ngraph,
                  output, time, csv_header, set_name, "DISP", time_step);
      }
      else if ((mi[1] > 3) && (mi[1] < 7))
      {
        sprintf(csv_header, "# displacements as (D1,D2,D3,...F%d) for set %s and time %s\n"
                            "node,D1,D2,D3",
                mi[1], set_name, time_step);

        for (j = 4; j <= mi[1]; j++)
        {
          sprintf(csv_header, "D%" ITGFORMAT ",", j);
        }
        sprintf(csv_header, "\n");

        csvgeneralvector(v, &iset, ntrans, filab, &nkcoords, inum, inotr,
                         trab, co, istartset, iendset, ialset, mi, ngraph,
                         output, csv_header, set_name, "DISP", time_step);
      }
    }
  }

  /*     storing the imaginary part of displacements in the nodes
         for the odd modes of cyclic symmetry calculations */

  if (*noddiam >= 0)
  {
    if ((strcmp1(filab, "U ") == 0) && (*ithermal != 2))
    {
      sprintf(csv_header, "# Imaginary part of displacements (D1,D2,D3) for set %s and time %s\n"
                          "node,D1,D2,D3\n",
              set_name, time_step);
      csvvector(v, &iset, ntrans, filab, &nkcoords, inum, inotr,
                trab, co, istartset, iendset, ialset, mi, ngraph,
                output, time, csv_header, set_name, "DISP", time_step);
    }
  }

  /*     storing the imaginary part of displacements in the nodes
         for steady state calculations */

  if ((*nmethod == 5) && (*mode == 0))
  {
    if ((strcmp1(filab, "U ") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(filab, set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# Imaginary part of displacements (D1,D2,D3) for set %s and time %s\n"
                          "node,D1,D2,D3\n",
              set_name, time_step);
      csvvector(v, &iset, ntrans, filab, &nkcoords, inum, inotr,
                trab, co, istartset, iendset, ialset, mi, ngraph,
                output, time, csv_header, set_name, "DISPI", time_step);
    }
  }

  /* storing the velocities in the nodes */

  if ((strcmp1(&filab[1740], "V   ") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[1740], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# velocities in the nodes as (V1,V2,V3) for set %s and time %s\n"
                        "node,V1,V2,V3\n",
            set_name, time_step);
    csvvector(v, &iset, ntrans, filab, &nkcoords, inum, inotr,
              trab, co, istartset, iendset, ialset, mi, ngraph,
              output, time, csv_header, set_name, "VELO", time_step);
  }

  /* storing the temperatures in the nodes */

  if (strcmp1(&filab[87], "NT  ") == 0)
  {
    iselect = 0;

    frdset(&filab[87], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# temperatures in the nodes as (T) for set %s and time %s\n"
                        "node,T\n",
            set_name, time_step);

    if (*ithermal <= 1)
    {
      csvselect(t1, t1, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncompscalar, ifieldscalar, icompscalar,
                nfieldscalar, &iselect, output, csv_header, set_name,
                "NDTEMP", time_step);
    }
    else
    {
      csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncompscalar, ifieldscalar, icompscalar,
                nfieldscalar, &iselect, output, csv_header, set_name,
                "NDTEMP", time_step);
    }
  }

  /* storing the electrical potential in the nodes */

  if ((strcmp1(&filab[3654], "POT ") == 0) && (*ithermal == 2))
  {
    iselect = 0;

    frdset(&filab[3654], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);
    sprintf(csv_header, "# electrical potential in the nodes as (V) for set %s and time %s\n"
                        "node,V\n",
            set_name, time_step);

    csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icompscalar,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "ELPOT", time_step);
  }

  /* storing the stresses in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[174], "S   ") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[174], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# stresses (SXX,SYY,SZZ,SXY,SYZ,SZX) for set %s and time %s\n"
                          "node,SXX,SYY,SZZ,SXY,SYZ,SZX\n",
              set_name, time_step);

      csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "STRESS", time_step);
    }
  }

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[4350], "SNEG") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[4350], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# stress at the negative surface of a shell as (SXX,SYY,SZZ,SXY,SYZ,SZX) for set %s and time %s\n"
                          "node,SXX,SYY,SZZ,SXY,SYZ,SZX\n",
              set_name, time_step);

      iflag = -1;
      numfield = 6;
      if ((norien > 0) && (strcmp1(&filab[4355], "L") == 0))
      {
        iorienloc = 1;
      }
      else
      {
        iorienloc = 0;
      }
      NNEW(inumshell, ITG, *nk);
      FORTRAN(extrapolateshell, (stx, stn, ipkon, inumshell, kon, lakon, &numfield,
                                 nk, ne, mi, &numfield, orab, ielorien, co,
                                 &iorienloc, &filab[4354], ielmat,
                                 thicke, ielprop, prop, &iflag));
      SFREE(inumshell);

      csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "STRNEG", time_step);
    }
  }

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[4437], "SMID") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[4437], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# stress at the midsurface of a shell as (SXX,SYY,SZZ,SXY,SYZ,SZX) for set %s and time %s\n"
                          "node,SXX,SYY,SZZ,SXY,SYZ,SZX\n",
              set_name, time_step);

      iflag = 0;
      numfield = 6;
      if ((norien > 0) && (strcmp1(&filab[4442], "L") == 0))
      {
        iorienloc = 1;
      }
      else
      {
        iorienloc = 0;
      }
      NNEW(inumshell, ITG, *nk);
      FORTRAN(extrapolateshell, (stx, stn, ipkon, inumshell, kon, lakon, &numfield,
                                 nk, ne, mi, &numfield, orab, ielorien, co,
                                 &iorienloc, &filab[4441], ielmat,
                                 thicke, ielprop, prop, &iflag));
      SFREE(inumshell);

      csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "STRMID", time_step);
    }
  }

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[4524], "SPOS") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[4524], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# stress at the positive surface of a shell as (SXX,SYY,SZZ,SXY,SYZ,SZX) for set %s and time %s\n"
                          "node,SXX,SYY,SZZ,SXY,SYZ,SZX\n",
              set_name, time_step);

      iflag = 1;
      numfield = 6;
      if ((norien > 0) && (strcmp1(&filab[4529], "L") == 0))
      {
        iorienloc = 1;
      }
      else
      {
        iorienloc = 0;
      }
      NNEW(inumshell, ITG, *nk);
      FORTRAN(extrapolateshell, (stx, stn, ipkon, inumshell, kon, lakon, &numfield,
                                 nk, ne, mi, &numfield, orab, ielorien, co,
                                 &iorienloc, &filab[4528], ielmat,
                                 thicke, ielprop, prop, &iflag));
      SFREE(inumshell);
      csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "STRPOS", time_step);
    }
  }

  /* storing the imaginary part of the stresses in the nodes
     for the odd modes of cyclic symmetry calculations */

  if (*noddiam >= 0)
  {
    if ((strcmp1(&filab[174], "S   ") == 0) && (*ithermal != 2))
    {

      sprintf(csv_header, "# imaginary part of the stresses (SXX,SYY,SZZ,SXY,SYZ,SZX) for set %s and time %s\n"
                          "node,SXX,SYY,SZZ,SXY,SYZ,SZX\n",
              set_name, time_step);

      csvselect(&stn[6 * *nk], stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "STRESSI", time_step);
    }
  }

  /* storing the imaginary part of the stresses in the nodes
     for steady state calculations */

  if ((*nmethod == 5) && (*mode == 0))
  {
    if ((strcmp1(&filab[174], "S   ") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[174], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);
      sprintf(csv_header, "# imaginary part of the stresses (SXX,SYY,SZZ,SXY,SYZ,SZX) for set %s and time %s\n"
                          "node,SXX,SYY,SZZ,SXY,SYZ,SZX\n",
              set_name, time_step);

      csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "STRESSI", time_step);
    }
  }

  /* storing the electromagnetic field E in the nodes */

  if ((strcmp1(&filab[3741], "EMFE") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[3741], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# electromagnetic field as (E1, E2, E3) for set %s and time %s\n"
                        "node,E1,E2,E3\n",
            set_name, time_step);

    csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompvector, ifieldvector, icompvector,
              nfieldtensor, &iselect, output, csv_header, set_name,
              "EMFE", time_step);

    if (*nmethod == 2)
    {
      iselect = 1;

      frdset(&filab[3741], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# Imaginary electromagnetic field as (E1, E2, E3) for set %s and time %s\n"
                          "node,E1,E2,E3\n",
              set_name, time_step);

      csvselect(&stn[6 * *nk], stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncompvector, ifieldvector, icompvector,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "EMFEI", time_step);
    }
  }

  /* storing the electromagnetic field B in the nodes */

  if ((strcmp1(&filab[3828], "EMFB") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[3828], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# electromagnetic field B as (B1, B2, B3) for set %s and time %s\n"
                        "node,B1,B2,B3\n",
            set_name, time_step);

    csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompvector, ifieldvector, icompvectorlast,
              nfieldtensor, &iselect, output, csv_header, set_name,
              "EMFB", time_step);

    if (*nmethod == 2)
    {
      iselect = 1;

      frdset(&filab[3828], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# imaginary electromagnetic field B as (B1, B2, B3) for set %s and time %s\n"
                          "node,B1,B2,B3\n",
              set_name, time_step);

      csvselect(&stn[6 * *nk], stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncompvector, ifieldvector, icompvectorlast,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "EMFBI", time_step);
    }
  }

  /* storing the total strains in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[261], "E   ") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[261], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);
      sprintf(csv_header, "# total strains as (EXX, EYY, EZZ, EXY, EYZ, EZX) for set %s and time %s\n"
                          "node,EXX,EYY,EZZ,EXY,EYZ,EZX\n",
              set_name, time_step);

      csvselect(een, een, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "TOSTRAIN", time_step);
    }
  }

  /* storing the imaginary part of the total strains in the nodes
     for the odd modes of cyclic symmetry calculations */

  if (*noddiam >= 0)
  {
    if ((strcmp1(&filab[261], "E   ") == 0) && (*ithermal != 2))
    {

      sprintf(csv_header, "# imaginary part of the total strains as (EXX, EYY, EZZ, EXY, EYZ, EZX) for set %s and time %s\n"
                          "node,EXX,EYY,EZZ,EXY,EYZ,EZX\n",
              set_name, time_step);

      csvselect(&een[6 * *nk], een, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "TOSTRAII", time_step);
    }
  }

  /* storing the imaginary part of the total strains in the nodes
     for steady state calculations */

  if ((*nmethod == 5) && (*mode == 0))
  {
    if ((strcmp1(&filab[261], "E   ") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[261], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# imaginary part of the total strains as (EXX, EYY, EZZ, EXY, EYZ, EZX) for set %s and time %s\n"
                          "node,EXX,EYY,EZZ,EXY,EYZ,EZX\n",
              set_name, time_step);

      csvselect(een, een, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "TOSTRAII", time_step);
    }
  }

  /* storing the mechanical strains in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[2697], "ME  ") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[2697], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# mechanical strains as (MEXX, MEYY, MEZZ, MEXY, MEYZ, MEZX) for set %s and time %s\n"
                          "node,MEXX,MEYY,MEZZ,MEXY,MEYZ,MEZX\n",
              set_name, time_step);

      csvselect(emn, emn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "MESTRAIN", time_step);
    }
  }

  /* storing the imaginary part of the mechanical strains in the nodes
     for the odd modes of cyclic symmetry calculations or for
     steady state calculations */

  if ((*noddiam >= 0) || ((*nmethod == 5) && (*mode == 0)))
  {
    if ((strcmp1(&filab[2697], "ME  ") == 0) && (*ithermal != 2))
    {

      /* check for a set in steady state calculations */

      if ((*nmethod == 5) && (*mode == 0))
      {
        iselect = 1;
        frdset(&filab[2697], set, &iset, istartset, iendset, ialset,
               inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
               ngraph);
      }

      sprintf(csv_header, "# imaginary part of the mechanical strains as (MEXX, MEYY, MEZZ, MEXY, MEYZ, MEZX) for set %s and time %s\n"
                          "node,MEXX,MEYY,MEZZ,MEXY,MEYZ,MEZX\n",
              set_name, time_step);

      if (*noddiam >= 0)
      {
        csvselect(&emn[6 * *nk], emn, &iset, &nkcoords, inum, istartset, iendset,
                  ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                  nfieldtensor, &iselect, output, csv_header, set_name,
                  "MESTRAII", time_step);
      }
      else
      {
        csvselect(emn, emn, &iset, &nkcoords, inum, istartset, iendset,
                  ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                  nfieldtensor, &iselect, output, csv_header, set_name,
                  "MESTRAII", time_step);
      }
    }
  }

  /* storing the thermal strains in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[4698], "THE ") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[4698], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# thermal strains as (THXX, THYY, THZZ, THXY, THYZ, THZX) for set %s and time %s\n"
                          "node,THXX,THYY,THZZ,THXY,THYZ,THZX\n",
              set_name, time_step);

      NNEW(ethn, double, 6 * *nk);
      for (i = 0; i < 6 * *nk; i++)
      {
        ethn[i] = een[i] - emn[i];
      }
      csvselect(ethn, ethn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "THSTRAIN", time_step);
      SFREE(ethn);
    }
  }

  /* storing the forces in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[348], "RF  ") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[348], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      if (mi[1] == 3)
      {
        sprintf(csv_header, "# forces as (F1,F2,F3) for set %s and time %s\n"
                            "node,F1,F2,F3\n",
                set_name, time_step);
        if ((iaxial == 1) && (strcmp1(&filab[352], "I") == 0))
        {
          for (i = 0; i < *nk; i++)
          {
            fn[1 + i * mt] *= 180.;
            fn[2 + i * mt] *= 180.;
            fn[3 + i * mt] *= 180.;
          }
        }
        csvvector(v, &iset, ntrans, &filab[348], &nkcoords, inum, inotr,
                  trab, co, istartset, iendset, ialset, mi, ngraph,
                  output, time, csv_header, set_name, "FORC", time_step);

        if ((iaxial == 1) && (strcmp1(&filab[352], "I") == 0))
        {
          for (i = 0; i < *nk; i++)
          {
            fn[1 + i * mt] /= 180.;
            fn[2 + i * mt] /= 180.;
            fn[3 + i * mt] /= 180.;
          }
        }
      }
      else if ((mi[1] > 3) && (mi[1] < 7))
      {
        sprintf(csv_header, "# forces as (F1,F2,F3,...F%d) for set %s and time %s\n"
                            "node,F1,F2,F3",
                mi[1], set_name, time_step);

        for (j = 4; j <= mi[1]; j++)
        {
          sprintf(csv_header, "F%" ITGFORMAT ",", j);
        }
        sprintf(csv_header, "\n");

        csvgeneralvector(v, &iset, ntrans, &filab[348], &nkcoords, inum, inotr,
                         trab, co, istartset, iendset, ialset, mi, ngraph,
                         output, csv_header, set_name, "FORC", time_step);
      }
      else
      {
        printf(" *WARNING in csv:\n");
        printf("          for output purposes only 4, 5 or 6\n");
        printf("          degrees of freedom are allowed\n");
        printf("          for generalized vectors;\n");
        printf("          actual degrees of freedom = %" ITGFORMAT "\n", mi[1]);
        printf("          output request ist not performed;\n");
      }
    }
  }

  /*     storing the imaginary part of the forces in the nodes
         for the odd modes of cyclic symmetry calculations or for
         steady state calculations */

  if ((*noddiam >= 0) || ((*nmethod == 5) && (*mode == 0)))
  {
    if ((strcmp1(&filab[348], "RF  ") == 0) && (*ithermal != 2))
    {

      /* check for a set in steady state calculations */

      if ((*nmethod == 5) && (*mode == 0))
      {
        iselect = 1;
        frdset(&filab[348], set, &iset, istartset, iendset, ialset,
               inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
               ngraph);
      }
      sprintf(csv_header, "# imaginary part of the forces as (F1, F2, F3) for set %s and time %s\n"
                          "node,F1,F2,F3\n",
              set_name, time_step);

      if (*noddiam >= 0)
      {
        csvvector(&fn[*nk * mt], &iset, ntrans, &filab[348], &nkcoords, inum, inotr,
                  trab, co, istartset, iendset, ialset, mi, ngraph,
                  output, time, csv_header, set_name, "FORCI", time_step);
      }
      else
      {
        csvvector(fn, &iset, ntrans, filab, &nkcoords, inum, inotr,
                  trab, co, istartset, iendset, ialset, mi, ngraph,
                  output, time, csv_header, set_name, "FORCI", time_step);
      }
    }
  }

  /* storing the equivalent plastic strains in the nodes */

  if ((strcmp1(&filab[435], "PEEQ") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[435], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# equivalent plastic strain as (PE) for set %s and time %s\n"
                        "node,PE\n",
            set_name, time_step);

    csvselect(ethn, ethn, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icompscalar,
              nfieldscalar, &iselect, output, csv_header, set_name,
              "PE", time_step);
  }

  /* storing the energy in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[522], "ENER") == 0) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[522], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# energy as (PE) for set %s and time %s\n"
                          "node,ENER\n",
              set_name, time_step);

      csvselect(enern, enern, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncompscalar, ifieldscalar, icompscalar,
                nfieldscalar, &iselect, output, csv_header, set_name,
                "ENER", time_step);
    }
  }

  /* storing the contact displacements and stresses at the slave nodes */

  /* node-to-face penalty */

  if ((strcmp1(&filab[2175], "CONT") == 0) && (*mortar != 1) && (*ithermal != 2) && ((*nmethod != 2) && (*nmethod != 13)))
  {

    for (i = *ne - 1; i >= 0; i--)
    {
      if ((strcmp1(&lakon[8 * i + 1], "S") != 0) || (strcmp1(&lakon[8 * i + 6], "C") != 0))
        break;
    }
    noutloc = *ne - i - 1;

    sprintf(csv_header, "# contact displacements and stresses at the slave as (COPEN, CSLIP1, CSLIP2, CPRESS, CSHEAR1, CSHEAR2) for set %s and time %s\n"
                        "node,COPEN,CSLIP1,CSLIP2,CPRESS,CSHEAR1,CSHEAR2\n",
            set_name, time_step);

    char csv_cont_filename[64];
    FILE *fcontact;

    for (i = *ne - 1; i >= 0; i--)
    {
      if ((strcmp1(&lakon[8 * i + 1], "S") != 0) || (strcmp1(&lakon[8 * i + 6], "C") != 0))
        break;
      sprintf(csv_cont_filename, "FILE_CONTACT_%s_%s.csv", set_name, time_step);
      fcontact = fopen(csv_cont_filename, "w");
      if (fcontact == NULL)
      {
        fprintf(stderr, "Couldn't open %s\n", csv_cont_filename);
        exit(1);
      }
      fprintf(fcontact, csv_header);

      strcpy1(text, &lakon[8 * i + 7], 1);
      nope = atoi(text) + 1;
      nodes = kon[ipkon[i] + nope - 1];

      fprintf(fcontact, "%d" ITGFORMAT ",", nodes);
      for (j = 0; j < 6; j++)
        fprintf(fcontact, ",%12.5E", (float)stx[6 * mi[0] * i + j]);
    }

    fprintf(fcontact, "\n");
    fclose(fcontact);
  }

  /* face-to-face penalty */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[2175], "CONT") == 0) && (*mortar == 1) && (*ithermal != 2))
    {
      iselect = 1;

      frdset(&filab[2175], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# contact displacements and stresses face-to-face penalty at the slave as (COPEN, CSLIP1, CSLIP2, CPRESS, CSHEAR1, CSHEAR2) for set %s and time %s\n"
                          "node,COPEN,CSLIP1,CSLIP2,CPRESS,CSHEAR1,CSHEAR2\n",
              set_name, time_step);

      csvselect(cdn, cdn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "CONTACT", time_step);
    }
  }

  /* storing imaginary part of the differential contact displacements
     and the contact stresses for the odd modes of cyclic symmetry
     calculations (only face-to-face penalty) */

  //  if((*noddiam>=0)||((*nmethod==5)&&(*mode==0))){
  if (*noddiam >= 0)
  {
    if ((strcmp1(&filab[2175], "CONT") == 0) && (*mortar == 1) && (*ithermal != 2))
    {
      iselect = 1;

      sprintf(csv_header, "# imaginary part of the differential contact displacements"
                          "# and the contact stresses for the odd modes of cyclic symmetry"
                          "# calculations (only face-to-face penalty)"
                          "# as (COPEN, CSLIP1, CSLIP2, CPRESS, CSHEAR1, CSHEAR2) for set %s and time %s\n"
                          "node,COPEN,CSLIP1,CSLIP2,CPRESS,CSHEAR1,CSHEAR2\n",
              set_name, time_step);

      csvselect(&cdn[6 * *nk], cdn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "CONTACTI", time_step);
    }
  }
  /* storing the contact energy at the slave nodes */

  if ((strcmp1(&filab[2262], "CELS") == 0) && (*ithermal != 2))
  {

    for (i = *ne - 1; i >= 0; i--)
    {
      if ((strcmp1(&lakon[8 * i + 1], "S") != 0) || (strcmp1(&lakon[8 * i + 6], "C") != 0))
        break;
    }
    noutloc = *ne - i - 1;

    sprintf(csv_header, "# contact Energy contact energy at the slave nodes as (CELS) for set %s and time %s\n"
                        "node,CELS\n",
            set_name, time_step);

    char csv_contactenergy_filename[64];
    FILE *fcontactenergy;

    for (i = *ne - 1; i >= 0; i--)
    {
      if ((strcmp1(&lakon[8 * i + 1], "S") != 0) || (strcmp1(&lakon[8 * i + 6], "C") != 0))
        break;
      nope = atoi(&lakon[8 * i + 7]) + 1;
      nodes = kon[ipkon[i] + nope - 1];

      sprintf(csv_contactenergy_filename, "FILE_CELS_%s_%s.csv", set_name, time_step);
      fcontactenergy = fopen(csv_contactenergy_filename, "w");
      if (fcontactenergy == NULL)
      {
        fprintf(stderr, "Couldn't open %s\n", csv_contactenergy_filename);
        exit(1);
      }

      fprintf(fcontactenergy, csv_header);

      fprintf(fcontactenergy, "%" ITGFORMAT ",%12.5E\n", nodes, (float)ener[2 * i * mi[0]]);
    }

    fprintf(fcontactenergy, "\n");
    fclose(fcontactenergy);
  }

  /* storing the internal state variables in the nodes */

  if (strcmp1(&filab[609], "SDV ") == 0)
  {
    iselect = 1;

    frdset(&filab[609], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# internal state variables for set %s and time %s\n"
                        "node,",
            set_name, time_step);

    for (j = 1; j <= *nstate_; j++)
    {
      sprintf(csv_header, "SDV%-3," ITGFORMAT, j);
    }
    sprintf(csv_header, "\n");

    for (i = 0; i < *nstate_; i++)
    {
      ifieldstate[i] = 1;
      icompstate[i] = i;
    }
    nfield[0] = *nstate_;
    csvselect(xstaten, xstaten, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, nstate_, ifieldstate, icompstate,
              nfield, &iselect, output, csv_header, set_name,
              "SDV", time_step);
  }

  /* storing the heat flux in the nodes
     the heat flux has been extrapolated from the integration points
     in subroutine extrapolate.f, taking into account whether the
     results are requested in the global system or in a local system.
     Therefore, subroutine frdvector cannot be used, since it assumes
     the values are stored in the global system */

  if ((strcmp1(&filab[696], "HFL ") == 0) && (*ithermal > 1))
  {
    iselect = 1;

    frdset(&filab[696], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# heat flux in the nodes (F1, F2, F3) for set %s and time %s\n"
                        "node,F1,F2,F3\n",
            set_name, time_step);

    csvselect(qfn, qfn, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompvector, ifieldvector, icompvector,
              nfieldvector1, &iselect, output, csv_header, set_name,
              "FLUX", time_step);
  }

  /* storing the electrical current in the nodes
     (cf. heat flux HFL above)  */

  if ((strcmp1(&filab[3567], "ECD ") == 0) && (*ithermal == 2))
  {
    iselect = 1;

    frdset(&filab[3567], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# electrical current in the nodes (j1, j2, j3) for set %s and time %s\n"
                        "node,j1,j2,j3\n",
            set_name, time_step);

    csvselect(qfn, qfn, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompvector, ifieldvector, icompvector,
              nfieldvector1, &iselect, output, csv_header, set_name,
              "CURR", time_step);
  }

  /* storing the heat generation in the nodes */

  if ((strcmp1(&filab[783], "RFL ") == 0) && (*ithermal > 1))
  {
    iselect = 1;

    frdset(&filab[783], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# heat generation in the nodes (RFL) for set %s and time %s\n"
                        "node,RFL\n",
            set_name, time_step);

    csvselect(fn, fn, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icompscalar,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "RFL", time_step);
  }

  /* storing the Zienkiewicz-Zhu improved stresses in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[1044], "ZZS") == 0) && (*ithermal != 2))
    {

      FORTRAN(zienzhu, (co, nk, kon, ipkon, lakon, ne0, stn, ipneigh, neigh,
                        stx, &mi[0]));

      iselect = 1;

      frdset(&filab[1044], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# Zienkiewicz-Zhu improved stresses in the nodes (SXX,SYY,SZZ,SXY,SYZ,SZX) for set %s and time %s\n"
                          "node,SXX,SYY,SZZ,SXY,SYZ,SZX\n",
              set_name, time_step);

      csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "ZZSTR", time_step);
    }
  }

  /* storing the imaginary part of the Zienkiewicz-Zhu
     improved stresses in the nodes
     for the odd modes of cyclic symmetry calculations or for
     steady state dynamics calculations */

  if ((*noddiam >= 0) || ((*nmethod == 5) && (*mode == 0)))
  {
    if ((strcmp1(&filab[1044], "ZZS") == 0) && (*ithermal != 2))
    {

      if (*noddiam >= 0)
      {

        /* cyclic symmetry := call from arpackcs.c */

        FORTRAN(zienzhu, (co, nk, kon, ipkon, lakon, ne0, stn, ipneigh, neigh,
                          &stx[6 * mi[0] * *ne], &mi[0]));
      }
      else
      {

        /* steady state := call from steadystate.c */

        FORTRAN(zienzhu, (co, nk, kon, ipkon, lakon, ne0, stn, ipneigh, neigh,
                          stx, &mi[0]));
        iselect = 1;
        frdset(&filab[1044], set, &iset, istartset, iendset, ialset,
               inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
               ngraph);
      }

      sprintf(csv_header, "# imaginary part of the Zienkiewicz-Zhu improved stresses in the nodes (SXX,SYY,SZZ,SXY,SYZ,SZX) for set %s and time %s\n"
                          "node,SXX,SYY,SZZ,SXY,SYZ,SZX\n",
              set_name, time_step);

      csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomptensor, ifieldtensor, icomptensor,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "ZZSTRI", time_step);
    }
  }

  /* storing the error estimator in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[1044], "ERR") == 0) && (*ithermal != 2))
    {

      NNEW(errn, double, 6 * *nk);

      nterms = 6;
      FORTRAN(errorestimator, (stx, errn, ipkon, kon, lakon, nk, ne,
                               mi, ielmat, &nterms, inum, co, v, &filab[1048],
                               ielprop, prop));

      iselect = 1;

      frdset(&filab[1044], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      ncomp = 1;
      ifield[0] = 1;
      icomp[0] = 0;

      sprintf(csv_header, "# error estimator in the nodes (STR(%%)) for set %s and time %s\n"
                          "node,STR(%%)\n",
              set_name, time_step);

      csvselect(errn, errn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomp, ifield, icomp,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "ERROR", time_step);
    }
  }

  /* storing the imaginary part of the error estimator in the nodes
     for the odd modes of cyclic symmetry calculations or for
     steady state dynamics calculations */

  if ((*noddiam >= 0) || ((*nmethod == 5) && (*mode == 0)))
  {
    if ((strcmp1(&filab[1044], "ERR") == 0) && (*ithermal != 2))
    {

      nterms = 6;
      if (*noddiam >= 0)
      {
        FORTRAN(errorestimator, (&stx[6 * mi[0] * *ne], stn, ipkon, kon, lakon, nk, ne,
                                 mi, ielmat, &nterms, inum, co, v, &filab[1048],
                                 ielprop, prop));
      }
      else
      {
        FORTRAN(errorestimator, (stx, stn, ipkon, kon, lakon, nk, ne,
                                 mi, ielmat, &nterms, inum, co, v, &filab[1048],
                                 ielprop, prop));
        iselect = 1;
        frdset(&filab[1044], set, &iset, istartset, iendset, ialset,
               inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
               ngraph);
      }

      sprintf(csv_header, "# Imaginary part of error estimator in the nodes (STR(%%)) for set %s and time %s\n"
                          "node,STR(%%)\n",
              set_name, time_step);

      ncomp = 1;
      ifield[0] = 1;
      icomp[0] = 0;

      csvselect(stn, stn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomp, ifield, icomp,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "ERRORI", time_step);
    }
  }

  /* storing the thermal error estimator in the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[2784], "HER") == 0) && (*ithermal > 1))
    {

      nterms = 3;
      FORTRAN(errorestimator, (qfx, qfn, ipkon, kon, lakon, nk, ne,
                               mi, ielmat, &nterms, inum, co, v, &filab[2788],
                               ielprop, prop));

      iselect = 1;

      frdset(&filab[2784], set, &iset, istartset, iendset, ialset,
             inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
             ngraph);

      sprintf(csv_header, "# thermal error estimator in the nodes (TEM(%%)) for set %s and time %s\n"
                          "node,TEM(%%)\n",
              set_name, time_step);

      ncomp = 1;
      ifield[0] = 1;
      icomp[0] = 0;

      csvselect(qfn, qfn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomp, ifield, icomp,
                nfieldvector1, &iselect, output, csv_header, set_name,
                "HERROR", time_step);
    }
  }

  /* storing the imaginary part of the thermal error estimator in the nodes
     for the odd modes of cyclic symmetry calculations */

  if ((*noddiam >= 0) || ((*nmethod == 5) && (*mode == 0)))
  {
    if ((strcmp1(&filab[2784], "HER") == 0) && (*ithermal > 1))
    {

      nterms = 3;
      FORTRAN(errorestimator, (&qfx[3 * mi[0] * *ne], qfn, ipkon, kon, lakon, nk, ne,
                               mi, ielmat, &nterms, inum, co, v, &filab[2788],
                               ielprop, prop));

      sprintf(csv_header, "# imaginary part of thermal error estimator in the nodes (TEM(%%)) for set %s and time %s\n"
                          "node,TEM(%%)\n",
              set_name, time_step);

      ncomp = 1;
      ifield[0] = 1;
      icomp[0] = 0;

      csvselect(qfn, qfn, &iset, &nkcoords, inum, istartset, iendset,
                ialset, ngraph, &ncomp, ifield, icomp,
                nfieldtensor, &iselect, output, csv_header, set_name,
                "HERRORI", time_step);
    }
  }

  /* storing the total temperatures in the network nodes */

  if ((strcmp1(&filab[1131], "TT  ") == 0) && (*ithermal > 1))
  {

    iselect = -1;
    frdset(&filab[1131], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# total temperatures in the network nodes (TT) for set %s and time %s\n"
                        "node,TT)\n",
            set_name, time_step);

    csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icompscalar,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "TOTEMP", time_step);
  }

  /* storing the mass flow in the network nodes */

  if ((strcmp1(&filab[1218], "MF  ") == 0) && (*ithermal > 1))
  {

    iselect = -1;
    frdset(&filab[1218], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# mass flow in the network nodes as (MF) for set %s and time %s\n"
                        "node,MF)\n",
            set_name, time_step);

    icomp[0] = 1;
    if ((iaxial == 1) && (strcmp1(&filab[1222], "I") == 0))
    {
      for (i = 0; i < *nk; i++)
        v[1 + i * mt] *= 180.;
    }
    csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icomp,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "MAFLOW", time_step);
    if ((iaxial == 1) && (strcmp1(&filab[1222], "I") == 0))
    {
      for (i = 0; i < *nk; i++)
        v[1 + i * mt] /= 180.;
    }
  }

  /* storing the total pressure in the network nodes */

  if ((strcmp1(&filab[1305], "PT  ") == 0) && (*ithermal > 1))
  {

    iselect = -1;
    frdset(&filab[1305], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# total pressure in the network nodes as (PT) for set %s and time %s\n"
                        "node,PT)\n",
            set_name, time_step);
    icomp[0] = 2;
    csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icomp,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "TOPRES", time_step);
  }

  /* storing the static pressure in the liquid network nodes */

  if ((strcmp1(&filab[1827], "PS  ") == 0) && (*ithermal > 1))
  {

    iselect = -1;
    frdset(&filab[1827], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# static pressure in the liquid network nodes as (PS) for set %s and time %s\n"
                        "node,PS)\n",
            set_name, time_step);

    icomp[0] = 2;

    csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icomp,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "STPRES", time_step);
  }

  /* storing the liquid depth in the channel nodes */

  if ((strcmp1(&filab[2349], "DEPT") == 0) && (*ithermal > 1))
  {

    iselect = -1;
    frdset(&filab[2349], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# liquid depth in the channel nodes nodes as (PS) for set %s and time %s\n"
                        "node,DEPTH)\n",
            set_name, time_step);

    icomp[0] = 2;

    csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icomp,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "DEPTH", time_step);
  }

  /* storing the critical depth in the channel nodes */

  if ((strcmp1(&filab[2436], "HCRI") == 0) && (*ithermal > 1))
  {

    iselect = -1;
    frdset(&filab[2436], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# critical depth in the channel nodes as (PS) for set %s and time %s\n"
                        "node,HCRIT)\n",
            set_name, time_step);

    icomp[0] = 3;

    csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icomp,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "HCRIT", time_step);
  }

  /* storing the static temperature in the network nodes */

  if ((strcmp1(&filab[1392], "TS  ") == 0) && (*ithermal > 1))
  {

    iselect = -1;
    frdset(&filab[1392], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# static temperature in the network nodes as (TS) for set %s and time %s\n"
                        "node,TS)\n",
            set_name, time_step);

    icomp[0] = 3;

    csvselect(v, v, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalar, ifieldscalar, icomp,
              nfieldvector0, &iselect, output, csv_header, set_name,
              "STTEMP", time_step);
  }

  /* mesh refinement */

  if (strcmp1(&filab[4089], "RM") == 0)
  {
    refinemesh(nk, ne, co, ipkon, kon, v, veold, stn, een, emn, epn, enern,
               qfn, errn, filab, mi, lakon, jobnamec, istartset, iendset,
               ialset, set, nset, matname, ithermal, output, nmat);
  }

  /* remove auxiliary field for the error estimator at the nodes */

  if ((*nmethod != 5) || (*mode == -1))
  {
    if ((strcmp1(&filab[1044], "ERR") == 0) && (*ithermal != 2))
    {
      SFREE(errn);
    }
  }

  /* storing the displacements in the nodes (magnitude, phase) */

  if ((strcmp1(&filab[870], "PU  ") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[870], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# displacements in the nodes (magnitude, phase) as (MAG1,MAG2,MAG3,PHA1,PHA2,PHA3) for set %s and time %s\n"
                        "node,MAG1,MAG2,MAG3,PHA1,PHA2,PHA3)\n",
            set_name, time_step);

    csvselect(vr, vi, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompvectph, ifieldvectph, icompvectph,
              nfieldvectph, &iselect, output, csv_header, set_name,
              "PDISP", time_step);
  }

  /* storing the temperatures in the nodes (magnitude, phase) */

  if ((strcmp1(&filab[957], "PNT ") == 0) && (*ithermal > 1))
  {
    iselect = 1;

    frdset(&filab[957], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# temperatures in the nodes (magnitude, phase) as (MAG1,MAG2,MAG3,PHA1,PHA2,PHA3) for set %s and time %s\n"
                        "node,MAG1,PHA1)\n",
            set_name, time_step);

    csvselect(vr, vi, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompscalph, ifieldscalph, icompscalph,
              nfieldscalph, &iselect, output, csv_header, set_name,
              "PNDTEMP", time_step);
  }

  /* storing the stresses in the nodes (magnitude, phase) */

  if ((strcmp1(&filab[1479], "PHS ") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[1479], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# stresses in the nodes (magnitude, phase)\n"
                        "# as (MAGXX,MAGYY,MAGZZ,MAGXY,MAGYZ,MAGZX,PHAXX,PHAYY,PHAZZ,PHAXY,PHAYZ,PHAZX)\n"
                        "# for set %s and time %s\n"
                        "node,MAGXX,MAGYY,MAGZZ,MAGXY,MAGYZ,MAGZX,PHAXX,PHAYY,PHAZZ,PHAXY,PHAYZ,PHAZX\n",
            set_name, time_step);

    csvselect(stnr, stni, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncomptensph, ifieldtensph, icomptensph,
              nfieldtensph, &iselect, output, csv_header, set_name,
              "PSTRESS", time_step);
  }

  /* storing the differential contact displacements and
     the contact stresses in the nodes (magnitude, phase)
     only for face-to-face penalty contact */

  if ((strcmp1(&filab[3915], "PCON") == 0) && (*ithermal != 2) && (*mortar == 1))
  {
    iselect = 1;

    frdset(&filab[3915], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# differential contact displacements and the contact stresses in the nodes (magnitude, phase)\n"
                        "# as (MAGO,MAGSL1,MAGSL2,MAGP,MAGSH1,MAGSH2,PHAO,PHASL1,PHASL2,PHAP,PHASH1,PHASH2)\n"
                        "# for set %s and time %s\n"
                        "node,MAGO,MAGSL1,MAGSL2,MAGP,MAGSH1,MAGSH2,PHAO,PHASL1,PHASL2,PHAP,PHASH1,PHASH2\n",
            set_name, time_step);

    csvselect(cdnr, cdni, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncomptensph, ifieldtensph, icomptensph,
              nfieldtensph, &iselect, output, csv_header, set_name,
              "PCONTAC", time_step);
  }

  /* storing the forces in the nodes (magnitude, phase) */

  if ((strcmp1(&filab[2610], "PRF ") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[2610], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# forces in the nodes (magnitude, phase) as (MAG1,MAG2,MAG3,PHA1,PHA2,PHA3)\n"
                        "# for set %s and time %s\n"
                        "node,MAG1,MAG2,MAG3,PHA1,PHA2,PHA3\n",
            set_name, time_step);

    csvselect(fnr, fni, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncompvectph, ifieldvectph, icompvectph,
              nfieldvectph, &iselect, output, csv_header, set_name,
              "PFORC", time_step);
  }

  /* the remaining parts are for frequency calculations with cyclic symmetry only */

  if ((*nmethod != 2) && (*nmethod != 13))
    return;

  /* storing the maximum displacements of the nodes in the base sector
     (components, magnitude) */

  if ((strcmp1(&filab[1566], "MAXU") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[1566], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "#  maximum displacements of the nodes in the base sector"
                        "# (components, magnitude) as (DX,DY,DZ,ANG)\n"
                        "# for set %s and time %s\n"
                        "node,DX,DY,DZ,ANG\n",
            set_name, time_step);

    ncomp = 4;
    ifield[0] = 1;
    icomp[0] = 1;
    ifield[1] = 1;
    icomp[1] = 2;
    ifield[2] = 1;
    icomp[2] = 3;
    ifield[3] = 1;
    icomp[3] = 0;
    nfield[0] = 4;
    nfield[1] = 4;

    csvselect(vmax, vmax, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncomp, ifield, icomp,
              nfield, &iselect, output, csv_header, set_name,
              "MDISP", time_step);
  }

  /* storing the worst principal stress at the nodes
     in the basis sector (components, magnitude)

     the worst principal stress is the maximum of the
     absolute value of all principal stresses, times
     its original sign */

  if ((strcmp1(&filab[1653], "MAXS") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[1653], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);
    sprintf(csv_header, "# worst principal stress at the nodes in the basis sector\n"
                        "# (components, magnitude). The worst principal stress is the\n"
                        "# maximum of the absolute value of all principal stresses, times\n"
                        "# its original sign (SXX,SYY,SZZ,SXY,SYZ,SZX,MAG)\n"
                        "# for set %s and time %s\n"
                        "node,SXX,SYY,SZZ,SXY,SYZ,SZX,MAG\n",
            set_name, time_step);

    ncomp = 7;
    ifield[0] = 1;
    icomp[0] = 1;
    ifield[1] = 1;
    icomp[1] = 2;
    ifield[2] = 1;
    icomp[2] = 3;
    ifield[3] = 1;
    icomp[3] = 4;
    ifield[4] = 1;
    icomp[4] = 6;
    ifield[5] = 1;
    icomp[5] = 5;
    ifield[6] = 1;
    icomp[6] = 0;
    nfield[0] = 7;
    nfield[1] = 7;

    csvselect(stnmax, stnmax, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncomp, ifield, icomp,
              nfield, &iselect, output, csv_header, set_name,
              "MSTRESS", time_step);
  }

  /* storing the worst principal strain at the nodes
     in the basis sector (components, magnitude)

     the worst principal strain is the maximum of the
     absolute value of all principal strains, times
     its original sign */

  if ((strcmp1(&filab[2523], "MAXE") == 0) && (*ithermal != 2))
  {
    iselect = 1;

    frdset(&filab[2523], set, &iset, istartset, iendset, ialset,
           inum, &noutloc, &nout, nset, &noutmin, &noutplus, &iselect,
           ngraph);

    sprintf(csv_header, "# Worst principal strain at the nodes in the basis sector\n"
                        "# (components, magnitude). The worst principal strain is the\n"
                        "# maximum of the absolute value of all principal stresses, times\n"
                        "# its original sign (EXX,EYY,EZZ,EXY,EYZ,EZX,MAG)\n"
                        "# for set %s and time %s\n"
                        "node,EXX,EYY,EZZ,EXY,EYZ,EZX,MAG\n",
            set_name, time_step);

    ncomp = 7;
    ifield[0] = 1;
    icomp[0] = 1;
    ifield[1] = 1;
    icomp[1] = 2;
    ifield[2] = 1;
    icomp[2] = 3;
    ifield[3] = 1;
    icomp[3] = 4;
    ifield[4] = 1;
    icomp[4] = 6;
    ifield[5] = 1;
    icomp[5] = 5;
    ifield[6] = 1;
    icomp[6] = 0;
    nfield[0] = 7;
    nfield[1] = 7;

    csvselect(eenmax, eenmax, &iset, &nkcoords, inum, istartset, iendset,
              ialset, ngraph, &ncomp, ifield, icomp,
              nfield, &iselect, output, csv_header, set_name,
              "MSTRAIN", time_step);
  }

  return;
}
