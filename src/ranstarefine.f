      SUBROUTINE RANSTArefine ( N, X, FU, EPS, F0, IER,cotet,
     &     kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C
C Random search to determine a point with a function value that is
C different from the one found at the starting point X
C This is an auxiliary subroutine to support the minimization
C subroutine FMINSI
C
C   FMINSI - Fortran subroutines for unconstrained function minimization
C   Copyright (C) 1986, 1993, 2001  Hugo Pfoertner
C
C   This library is free software; you can redistribute it and/or
C   modify it under the terms of the GNU Lesser General Public
C   License as published by the Free Software Foundation; either
C   version 2.1 of the License, or (at your option) any later version.
C
C   This library is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C   Lesser General Public License for more details.
C
C   You should have received a copy of the GNU Lesser General Public
C   License along with this library; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C   USA
C
C   Contact info: mailto:hugo@pfoertner.org
C   or use the information provided at http://www.pfoertner.org/
C
C Author: Hugo Pfoertner, Oberhaching, Germany
C
C Version History (left in German):
! 22.06.21 adapted fminsi for the mesh refinement procedure in CalculiX
!          Author: Guido Dhondt
C 02.06.01 English translation of comments, LGPL header added
C 27.11.93 ZUFALLSZAHLENGENERATOR RANEWR UND INIRAN EINGEBAUT
C 16.11.92 DIMENSIONIERUNG AUF 128 ERHOEHT
C 15.01.91 ABFRAGE AUF AENDERUNG STATT AUF VERMINDERUNG DES
C          FUNKTIONSWERTES
C 09.04.86 BASISVERSION
C
C The meaning of the parameters is the same as in subroutine FMINSI
      implicit doubleprecision (a-h,o-z)
      INTEGER N, IER,kontet(4,*),ipoeln(*),ieln(2,*),node,iedge,
     &     ipoeled(*),ieled(2,*),iedgmid(*),iedtet(6,*)
      doubleprecision FU, F0, X(N), EPS(N),cotet(3,*)
      EXTERNAL FU
C Output:
C F0 ...  Changed function value
C X ...   Variable vector of point with different function value
C         if such a point has been found, unchanged otherwise
C IER ... 0, if a point with changed function value has been found
C         3 otherwise (search stopped without success after 100*N
C         function evaluations.
C
C Dimension of local arrays has to be compliant with corresponding
C size within FMINSI
      PARAMETER ( NMAX=216 )
      doubleprecision RS(NMAX), XS(NMAX)
C Uniform random number generator and corresponding initialization
      REAL RANUWH
      EXTERNAL INIRAN, RANUWH
C
C
C Preset return code with "No success"
      IER = 3
C
C Initialize uniform random number generator
C
      CALL INIRAN
C
C Set initial step size
      DO 10 K = 1, N
      RS(K) = 10. * EPS(K)
10    CONTINUE
C
C After N2 function evaluations, an increase in variance is tried
      N2 = 2**N
C
C Maximum number of function evaluations
      NT = 100 * N
C
C Function value at starting point
      F0 = FU ( N, X,cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
      LL = 0
C
C Loop over maximum number of trials
      DO 20 L = 1, NT
      DO 30 NN = 1, N
C
C Set co-ordinates to starting point + random increment
C in the range +-RS(i)
      XS(NN) = X(NN) + RS(NN) * 2. * ( dble(RANUWH()) - 0.5 )
30    CONTINUE
C
C Corresponding function value
      H = FU ( N, XS,cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C
C Check for change
      IF ( H .NE. F0 ) THEN
C
C Loop is terminated at first occurrence of changed value
        F0 = H
        DO 40 K = 1, N
40      X(K) = XS(K)
        IER = 0
        GOTO 999
      ENDIF
C
C After N2 trials without success, the variance is increased
C
      LL = LL + 1
      IF ( LL .GE. N2 ) THEN
        DO 50 K = 1, N
50      RS(K) = RS(K) + RS(K)
        LL = 0
      ENDIF
C
C End of loop over maximum number of trials
20    CONTINUE
C
999   CONTINUE
      RETURN
C End of subroutine RANSTA
      END
