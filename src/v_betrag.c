
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>


double v_betrag(double *a)
/* ********************************************************* */
/*      laenge von Vektor a                                  */
/* ********************************************************* */
{
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}


