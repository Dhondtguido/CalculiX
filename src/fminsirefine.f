      SUBROUTINE FMINSIrefine ( N, X, FU, EPS, FMIN, IER,cotet,
     &     kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C
C Minimization of a function of N variables with a polytope method,
C using only function values.
C
C   FMINSI - Fortran subroutines for unconstrained function minimization
C   Copyright (C) 1992, 2001  Hugo Pfoertner
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
C Version History (German Language kept for authenticity ;-)
C
! 22.06.21 adapted fminsi for the mesh refinement procedure in CalculiX
!          Author: Guido Dhondt
C 21.05.05 Tentative doubleprecision version
C 03.06.01 Translation finished, additional comments added, LGPL
C 21.04.01 Start English translation of comments
C 30.11.92 KOMPLETTE KONTRAKTION NUR BEI VERSAGEN DES EINFACHEN
C          KONTRAKTIONSSCHRITTES
C 27.11.92 KLEINSTER STARTSCHRITT: RMULT * EPS(I)
C 16.11.92 MAXIMALE ITERATIONSZAHL = 400*N + N**3/2,
C          ABFRAGE FUER WAHL DES EXPANDIERTEN PUNKTES GEAENDERT,
C          MITTELPUNKTSBERECHNUNG UEBER AKTUALISIERUNG, REDUKTION DES
C          EXPANSIONSFAKTORS BEI MEHRFACHER ERFOLGLOSER EXPANSION.
C 16.01.91 BESCHRAENKUNG DER STARTSCHRITTEXPANSION AUF 10MAL
C 15.01.91 BEI STARTWERTSUCHE MITTELS RANSTA WIRD NUR NOCH
C          NACH VERAENDERTEM, NICHT MEHR NACH KLEINEREM F GESUCHT.
C 13.06.89 EPSI-BERECHNUNG MODIFIZIERT, VERSUCH EINER WIEDERHERSTELLUNG
C          DER VOM ANWENDER ANGEGEBENEN GENAUIGKEITSSCHRANKEN,
C          ZUSAETZLICHES KONVERGENZKRITERIUM : 3*(N+1) ERFOLGLOSE
C          KONTRAKTIONSVERSUCHE.
C 09.06.89 EXPANSIONSFAKTOR GAMMA = 1 + 3/N
C 28.10.87 ZUSAETZLICHE DIAGNOSEMELDUNGEN BEI GRADIENTENBERECHNUNG
C          FUER STARTSCHRITT
C 10.04.86 RANDOM-STARTWERTSUCHE UND DIAGNOSEMELDUNGEN
C 1982     DIVERSE ERGAENZUNGEN (BUCH GILL/MURRAY),
C          VERBESSERTE KONVERGENZKRITERIEN (EINFUEHRUNG EPSI)
C 1978     BASISVERSION (ORIGINAL NELDER/MEAD ALGORITHMUS)
C
C References:
C NELDER J.A., MEAD R.: A SIMPLEX METHOD FOR FUNCTION MINIMIZATION.
C COMPUTER JOURNAL 7, 308-313, 1965
C GILL G.E., MURRAY W., WRIGHT M.H.: PRACTICAL OPTIMIZATION.
C ACADEMIC PRESS, LONDON, 1981 (S. 94-96)
C BARABINO G.P. ET AL : A STUDY ON THE PERFORMANCE OF SIMPLEX METHODS
C FOR FUNCTION MINIMIZATION. IN PROCEEDINGS OF THE IEEE INTERNATIONAL
C CONFERENCE ON CIRCUITS AND COMPUTERS ICCC80, IEEE, NEW YORK, 1980
C (S. 1150-1153)
C SMITH J.D., HOLLOMAN M.E., OTTO W.E.: SIMMIN - A PROGRAM FOR
C MINIMIZING FUNCTIONS BY SIMPLEX METHODS FOR IBM PC.
C U.S. ARMY MISSILE COMMAND, REDSTONE ARSENAL, ALABAMA, DIRECTED
C ENERGY DIRECTORATE,  TECHNICAL REPORT RD-DE-87-4, JULY 1987
C
C Description of parameters:
C
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      INTEGER N, IER,kontet(4,*),ipoeln(*),ieln(2,*),node,iedge,
     &     ipoeled(*),ieled(2,*),iedgmid(*),iedtet(6,*)
      DOUBLEPRECISION FU, FMIN, X(N), EPS(N),cotet(3,*)
      EXTERNAL FU
C N ... Number of variables.
C       Local dimensions currently require N <= 128
C X(N) ... vector of variables. When calling FMINSI the user has to
C       provide the starting point of the minimization, i.e. an
C       assumed or previously determined estimate of the location
C       of a mimimum.
C       After completion of FMINSI X contains the best approximation
C       found for the loacation of a local minimum of the objective
C       function.
C FU(N,X) ... REAL FUNCTION, ojective function to be miminized,
C       need neither be differentiable nor steady.
C       FU has to declared as "EXTERNAL" in the calling program.
C EPS ... vector of desired absolute accuracy for the single components
C       of the result vector. The minimization is terminated, when the
C       change of all components of the variable vector drops below
C       the respective EPSI(I) in one optimization step. EPSI(I) is
C       derived from the input vector EPS(I) by applying some
C       plausibility checks to avoid using convergence limits not
C       attainable with single precision arithmetic.
C FMIN ... Lowest function value found at termination
C IER ... "Return code", when FMINSI is terminated regularly, i.e. if
C       a minimum has been found satisfying the convergence criterion,
C       IER=0 is set on return.
C       Situations returning IER>0:
C       IER=1: N greater than max allowable dimension NMAX
C       IER=2: No convergence reached within ITMAX (see below) steps.
C              X contains the best approximation to the location of
C              the minimum found when the iteration was terminated.
C       IER=3: Objective function constant in the vicinity of the start
C              location. Sometimes this situation is simply caused by
C              the programmers failure to assign a value to the
C              objective function
C  The input value of IER is used to activate a printout of an
C  execution trace to standard output using write (*,...) ...
C  For a full visibility of the execution flow it is advisable to
C  include printouts of the variable vector and of the computed
C  function value into the source of the objective function.
C  As FMINSI overwrites IER on exit, it is necessary to restore IER
C  to a negative value before calling FMINSI again, when trace output
C  is required for consecutive executions.
C
C Declaration of local variables (undeclared variables use the
C Fortran standard implicit type convention (I..N Integer), else Real
C
C FMINSI uses single precision reals, which I've found to be always
C sufficient within the optimization routine in 25 years of optimization
C experience, also with "industrial grade" problems.
C Double precision is usually only needed for some intermediate
C calculations within the objective functions. With the spreading use of
C 64 bit processors this problem gradually vanishes in the long term
C
      LOGICAL TEST, MODEPS, AVGUPD
C
C
C This limits the problem dimension, using the method for higher N
C requires a corresponding modification of NMAX.
C It is not recommended to use FMINSI for such high problem dimensions,
C mainly due to the inevitable progressive collapse of the search
C polytope during the iteration (currently there is no mechanism to
C restore the polytope to full rank other than restarting by a new call)
C
      PARAMETER ( NMAX=216, NMP1=NMAX+1 )
      DIMENSION EPSI(NMAX), SIM(NMAX,NMP1), F(NMP1)
     &         ,XAV(NMAX), XST(NMAX), XSTST(NMAX)
C
C The function name is also passed into the call sequence of
C subroutine RANSTA, therefore
C**** EXTERNAL FU (declaration already made above)
C
C Factors for shape adjustment of the search polytope
C
      DATA ALPHA, BETA, GAMM0, GAMM1 / 1.0, 0.5, 1.0, 3.0 /,
     &     RMULT / 8.0 /,  ISTMAX / 10 /
C
C Some other real constants
C
      DATA ZERO, ONE, HALF / 0.0, 1.0, 0.5 /
     &    ,SMALL, VSMALL / 1.E-9, 1.E-30 /
C
C Activate printout of execution trace
      TEST = IER .LT. 0
C
C Maximum number of iterations, found to be a good compromise except
C for really pathological problems
C
C      ITMAX = 400 * N  +  N**3 / 2
      ITMAX = 1000 * N  +  2*N**3
C
C Check for admissible range of problem dimension
C
      IER = 0
      IF ( TEST ) WRITE (*,*) ' START FMINSI, VERSION 1992/11/30'
      IF ( N .LT. 1 ) THEN
        WRITE ( *, 1003 ) N
1003    FORMAT (' +++++FMINSI : N =',I3,' ILLEGAL')
        IER = 1
        GOTO 999
      ENDIF
C
      IF ( N .GT. NMAX ) THEN
        WRITE ( *, 1000 ) N
1000    FORMAT(' +++++FMINSI: DIMENSION N =',I4,' TOO LARGE')
        IER = 1
        GOTO 999
      ENDIF
C
C Expansion factor variable and dependent on problem dimension:
C This is the most important enhancement of the original method.
C The idea is due to Barabino et al. (see references above),
C but they did not provide a method to determine the factor.
C The formula used here was found to be the best compromise in
C a parameter study performed in 1989 with quadratic problems for N up
C to 100, using random coefficients plus a modification of the diagonal
C elements to get a prescribed problem condition number.
C
      GAMMA = GAMM0 + GAMM1 / REAL ( N )
      IF ( TEST ) WRITE (*,*) ' EXPANSION FACTOR GAMMA = ', GAMMA
C
C Starting value of the expansion factor. The factor is dynamically
C adjusted. It is reduced if successive expansion steps fail.
C Following successul expansion steps it is restored.
C
      GAMMD = GAMMA
C
C Number of vertices of the simplex
      NP = N + 1
C
C Auxiliary value 1 / N
      EDN = 1.0 / REAL(N)
C
C Initial choice of the method to determine average of simplex points
      AVGUPD = .FALSE.
C
C Check plausibility of user supplied absolute accuracy vector.
C Internally a (modified if necessary) vector EPSI is used.
C
      MODEPS = .FALSE.
      DO 305 I = 1, N
      EPSI(I) = MAX ( ABS (EPS(I)), 1.0E-8 * ABS (X(I)) )
      IF ( EPSI(I) .LT. VSMALL ) EPSI(I) = MAX (SMALL,ABS(X(I)*SMALL))
      IF ( EPS(I) .NE. EPSI(I) ) MODEPS = .TRUE.
305   CONTINUE
      IF ( TEST .AND. MODEPS )
     &      WRITE (*,*) ' CONVERGENCE LIMITS MODIFIED'
C
C Use a heuristic method to determine step sizes to compute
C an approximation of the gradient direction at the starting point
      DO 306 I = 1, N
      XAV(I) = MAX ( 0.00001*ABS(X(I)), EPSI(I) )
306   CONTINUE
      IF ( TEST ) WRITE (*,*)
     & ' STEP SIZE TO ESTIMATE GRADIENT : ', (XAV(K),K=1,N)
C
C Determine the possible expansion of an EPSI sized simplex
C by performing an approximate steepest descent step
C
C Difference approximation of -G at starting point, stored on XSTST
      IF ( TEST ) WRITE (*,*) ' DIFFERENCE APPROXIMATION FOR -G :'
      F0 = FU ( N, X(1),cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
      IF ( TEST ) WRITE (*,*) ' FUNCTION VALUE AT STARTING POINT : ', F0
      DO 310 I = 1, N
      DO 320 J = 1, N
320   XST(J) = X(J)
      XST(I) = XST(I) + XAV(I)
      IF ( TEST ) WRITE (*,*) ' COMPONENT ', I,
     &            ' X = ', (XST(K),K=1,N)
      H = FU ( N, XST(1),cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
      IF ( TEST ) WRITE (*,*) ' FUNCTION VALUE : ', H
      F(I) = XAV(I)
C
C Determine positive or negative direction along co-ordinate axes
C changed 27.11.92, old version:
C***  IF ( F0 .LT. H ) F(I) = -HALF*F(I)
      IF ( F0 .LT. H ) F(I) = - F(I)
      XSTST(I) = ( F0 - H ) / XAV(I)
310   CONTINUE
C
C Normalize G to the maximum norm of EPS
      HG = ZERO
      HE = ZERO
      DO 330 I = 1, N
      HG = MAX ( HG, ABS(XSTST(I)) )
      HE = MAX ( HE, XAV(I) )
330   CONTINUE
C
C If the norm of the gradient is too small, i.e. the function
C is approximately constant in the vicinity of the user supplied
C starting point, then a random search with repeated increase
C of variance is tried around the starting point, with the aim
C of to get away from the flat region.
C
      IF ( HG .LT. VSMALL ) THEN
        IF ( TEST )
     &  WRITE (*,*) ' RANDOM SEARCH FOR BETTER STARTING POINT'
        DO 333 I = 1, N
333     XAV(I)= X(I)
C
C The search is performed by a separate subroutine RANSTA.
C It either returns a point with a different function value
C or IER=3 if the search failed.
C
        CALL RANSTArefine ( N, XAV, FU, EPSI, F1, IER,cotet,
     &       kontet,ipoeln,ieln,node,iedge,
     &       ipoeled,ieled,iedgmid,iedtet )
        IF ( IER .EQ. 0 ) THEN
          IF ( TEST )
     &    WRITE (*,*) ' SEARCH FOR START POINT SUCCESSFUL, F1 = ', F1
C
C Check for increase of F1
C
      IF ( F1 .GT. F0 ) THEN
        IF ( TEST )
     &  WRITE (*,*) ' SEARCH DIRECTION REVERSED'
        DO 336 I = 1, N
        F(I) = X(I)
        X(I) = XAV(I)
        XAV(I) = F(I)
336     CONTINUE
        H = F1
        F1 = F0
        F0 = H
      ENDIF
C
C Difference vector from starting point
C
          DO 337 I = 1, N
          XSTST(I) = XAV(I) - X(I)
C
C Steps along the co-ordinate axes are  limited to a minimum
C size of EPSI(i)
C
          IF ( ABS (XSTST(I)) .LT. EPSI(I) ) THEN
            F(I) = SIGN ( EPSI(I), XSTST(I) )
          ELSE
            F(I) = XSTST(I)
          ENDIF
337       CONTINUE
C
C Skip forward, if starting steps have been successfully determined
C by the random search (Label 345 approx 16 lines below)
C
C A short note on software quality:
C Sorry to all those software gurus for the GOTOs, but the roots
C of FMINSI date back to the mid seventies, using Fortran IV.
C You know "IF(I-J)10,20,10" ....
C At that time I didn't bother about maintainability of software.
C BTW, I have written Ada programs used in aircraft on-board systems,
C where software quality assurance would have sent me to hell for using
C only 5% of the control flow complexity of this program.
C I never tried to determine the McCabe of this code ;-)
C
          GOTO 345
        ENDIF
        WRITE ( *, 1002 ) F0
1002    FORMAT(' +++++FMINSI: Function approximately constant',
     &         ' around starting point',/,14X,'F= ', E13.6 )
        FMIN = F0
        GOTO 999
      ENDIF
C
      HG = HE / HG
C
      DO 340 I = 1, N
340   XSTST(I) = XSTST(I) * HG
C
C Label 345 is the target location, if the first search direction
C was the result of a random search
C
345   CONTINUE
C
C Perform a quasi steepest descent step
      START = ONE / RMULT
      J = 0
C
C Start of expansion loop
C
350   CONTINUE
C
C Increase counter J for the number of expansion steps
C
      J = J + 1
      START = START * RMULT
      DO 360 I = 1, N
360   X(I) = X(I) + START * XSTST(I)
      F1 = FU ( N, X(1),cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C
C Expansion is terminated after ISTMAX steps (skip downwards)
C
      IF ( F1 .GT. F0  .OR.  J .GT. ISTMAX ) GOTO 370
C
C Size of starting step can be increased further
C
      F0 = F1
C
C Skip back to the next expansion step
C
      GOTO 350
C
C The point before the last one was the best
C
370   CONTINUE
      F(NP) = F0
C
C Restore the corresponding location
C
      DO 380 I = 1, N
380   X(I) = X(I) - START * XSTST(I)
      H = RMULT
      START = START / H
C
      START = MAX ( START, RMULT )
C
390   CONTINUE
C
      IF ( TEST )
     &  WRITE (*,*) ' MULT=',START,' EPSI:',(EPSI(I),I=1,N)
C
C Create initial simplex, using the start expansion factor START
C and the co-ordinate step sizes stored in vector F
C
      DO 10 I = 1, N
      SIM(I,NP) = X(I)
10    CONTINUE
C
      DO 20 I = 1, N
      DO 25 J = 1, N
25    SIM(J,I) = X(J)
      SIM(I,I) = SIM(I,I) + START * F(I)
20    CONTINUE
C
C Compute function values at vertices of simplex
C
      DO 30 I = 1, N
30    F(I) = FU ( N, SIM(1,I),cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C
C Initialize iteration count L
C
      L = 0
C
C Initialize no-success counters for expansion and contraction steps
C
      NSUCC = 0
      NOXSUC = 0
C
C Start of main iteration loop
C ============================
C
35    CONTINUE
C
      L = L + 1
C
C Test for termination due to exceeding the maximum iteration count
C
      IF ( L .GT. ITMAX ) THEN
        IER = 2
        WRITE (*,1001) L
1001    FORMAT(' +++++FMINSI: STOPPED AFTER',I7,' ITERATIONS',
     &  /      '              NO CONVERGENCE REACHED')
C
C Copy best approximation to output and jump to exit
C
        DO 34 I = 1, N
34      X(I) = SIM(I,1)
        FMIN = F(1)
        GOTO 999
      ENDIF
C
C Search for minimum and maximum function value
C
36    CONTINUE
      ILO = 1
      IHI = 1
      FLO = F(1)
      FHI = F(1)
C
      DO 50 I = 2, NP
      IF ( F(I) .LT. FLO ) THEN
        FLO = F(I)
        ILO = I
      ELSE IF ( F(I) .GT. FHI ) THEN
        FHI = F(I)
        IHI = I
      ENDIF
50    CONTINUE
      IF ( TEST )
     & WRITE (*,*) ' FLO=',FLO, ' IHI=',IHI, ' FHI=',FHI
C
C Convergence is assumed if function is constant within simplex
C
      IF ( ILO .EQ. IHI ) THEN
        DO 55 I = 1, N
55      X(I) = SIM(I,ILO)
        FMIN = F(ILO)
        IF ( TEST ) WRITE (*,*) ' FUNCTION CONSTANT =', FMIN
        GOTO 999
      ENDIF
C
C Move point with highest function value to index NP
C
      IF ( ILO .EQ. NP ) ILO = IHI
      IF ( IHI .NE. NP ) THEN
        DO 60 I = 1, N
        H = SIM(I,NP)
        SIM(I,NP) = SIM(I,IHI)
        SIM(I,IHI) = H
60      CONTINUE
        H = F(NP)
        F(NP) = FHI
        F(IHI) = H
C
C Update location of center of gravity
C (recomputing has been found to be the most expensive loop for
C higher N in a performance analysis of the program)
C
        IF ( AVGUPD ) THEN
          DO 65 I = 1, N
          XAV(I) = XAV(I) + EDN * ( SIM(I,IHI) - SIM(I,NP) )
65        CONTINUE
        ENDIF
      ENDIF
C
C Move point with lowest function value to index 1
C
      IF ( ILO .NE. 1 ) THEN
        DO 70 I = 1, N
        H = SIM(I,1)
        SIM(I,1) = SIM(I,ILO)
        SIM(I,ILO) = H
70      CONTINUE
        H = F(1)
        F(1) = FLO
        F(ILO) = H
      ENDIF
C
C Auxiliary convergence test
C (proposed by Smith et al., see references above)
C Check for many unsuccessful contraction steps
C
      IF ( NSUCC .GT. 3*NP ) THEN
        IF ( TEST ) WRITE (*,*) ' STOPPED AFTER ', NSUCC,
     &    ' FAILED CONTRACTION STEPS'
        GOTO 712
      ENDIF
C
C Target label for repeating the convergence test, if the
C user supplied EPS has been modified during the initial
C plausibility check.
C
701   CONTINUE
C
C Test for convergence: Assume convergence, if all components
C of the difference vector between the maximum and minimum points
C of the simplex have an absolute value less than EPSI(i)
C
      DO 71 I = 1, N
      IF ( ABS ( SIM(I,NP)-SIM(I,1) ) .GT. EPSI(I) ) GOTO 75
71    CONTINUE
C
C Check if the convergence test was passed against modified limits
C
      IF ( MODEPS ) THEN
C
C Try to restore user supplied limits, unless too demanding
C for single precision arithmetic
C
        DO 711 I = 1, N
        EPSI(I) = MAX ( ABS(EPS(I)), ABS(SIM(I,1))*1.0E-5 )
711     CONTINUE
        MODEPS = .FALSE.
        IF ( TEST ) WRITE (*,*) ' CONVERGENCE LIMITS RESTORED'
        GOTO 701
      ENDIF
C
C Target label after exceeding NSUCC
C
712   CONTINUE
C
C Convergence reached
C
      IF ( TEST )
     &    WRITE (*,*) ' CONVERGENCE ACHIEVED, FMIN=', F(1)
      DO 72 I = 1, N
72    X(I) = SIM(I,1)
      FMIN = F(1)
      GOTO 999
C
75    CONTINUE
C
C Compute the center of gravity of all points, but
C omitting the point with the highest function value.
C Normally this is done by an update formula, but after
C a contraction of the whole simplex it is necessary to
C recompute the whole sum.
C
      IF ( .NOT. AVGUPD ) THEN
C
        DO 80 I = 1, N
        H = ZERO
        DO 85 J = 1, N
85      H = H + SIM(I,J)
        XAV(I) = H * EDN
80      CONTINUE
C
C Reset state to "use update formula"
        AVGUPD = .TRUE.
      ENDIF
C
C Reflexion: Reverse the principal search direction
C
      DO 90 I = 1, N
90    XST(I) = ( ONE + ALPHA ) * XAV(I) - ALPHA * SIM(I,NP)
      FST = FU ( N, XST(1),cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C
C Check if the reflected point is a new minimum (skip forward)
C
      IF ( FST .GE. FLO ) GOTO 110
C
C Reset counter for unsuccessful contractions
C
      NSUCC = 0
      IF ( TEST )
     &  WRITE (*,*) ' REFLEXION SUCCESSFUL, FST =',FST
C
C Try an expansion into the direction of the reflected point
C
      DO 120 I = 1, N
120   XSTST(I) = ( ONE - GAMMD ) * XAV(I) + GAMMD * XST(I)
      FSTST = FU ( N, XSTST(1),cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C
C Changed 16.11.92 : Result of expansion is compared against
C result of reflexion and not against minimum
C Instead of:  IF ( FSTST .GE. FLO ) THEN
C New:
      IF ( FSTST .GT. FST ) THEN
C
C Increase counter for unsuccessful expansions
C
        NOXSUC = NOXSUC + 1
        IF ( NOXSUC .GT. N ) THEN
C
C Reduce expansion factor
C
          GAMMD = HALF * ( GAMM0 + GAMMD )
          IF ( TEST )
     &      WRITE (*,*) ' EXPANSION FACTOR REDUCED TO ', GAMMD
        ENDIF
C Skip forward
        GOTO 140
      ENDIF
C
C A new minimum was found by the expansion step
C
      IF ( TEST )
     &  WRITE (*,*) ' EXPANSION SUCCESSFUL, FSTST =', FSTST
C
C Reset all non-success counters
C
      NSUCC = 0
      NOXSUC = 0
C
C Restore expansion factor to the (dimension dependent)
C default value
C
      GAMMD = GAMMA
C
C Target label after successful contraction
C
130   CONTINUE
C
      DO 150 I = 1, N
150   SIM(I,NP) = XSTST(I)
      F(NP) = FSTST
C
C Jump back to the next iteration step
C ====================================
C
      GOTO 35
C
C Reflexion did not find a new minimum
C
110   CONTINUE
C
C Check if reflexion result is at least better than one of the
C remaining points
C
      DO 170 I = 2, N
      IF ( FST .LT. F(I) ) GOTO 140
170   CONTINUE
C
C Reflected point does not lead to an exchange of maximum
C
      GOTO 171
C
C Target label, if the reflexion result can be used to replace
C the worst point
C
140   CONTINUE
C
C Replace worst point by XST
C
      DO 220 I = 1, N
220   SIM(I,NP) = XST(I)
      F(NP) = FST
C Reset counter for unsuccessful contractions
      NSUCC = 0
C
C Skip back to next iteration step
C ================================
C
      GOTO 35
C
C Target label, if reflected point did not lead to an exchange
C of the worst point
C
171   CONTINUE
C
C Check if reflected point is at least better than the previous
C maximum
C
      IF ( FST .LT. FHI ) THEN
        DO 173 I = 1, N
173     SIM(I,NP) = XST(I)
        F(NP) = FST
        FHI = FST
      ENDIF
C
C Reflected point could not replace maximum.
C Try to contract point with highest F towards center of gravity.
C
      IF ( TEST )
     & WRITE (*,*) ' CONTRACTION TOWARDS CENTER'
      DO 180 I = 1, N
180   XSTST(I) = BETA * SIM(I,NP) + ( ONE - BETA ) * XAV(I)
      FSTST = FU ( N, XSTST(1),cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C Changed 30.11.92, Instead of :
C***  IF ( FSTST .LT. FHI .AND. FSTST .GT. FLO ) GOTO 130
      IF ( FSTST .LT. FLO ) THEN
C
C Contracted point is new minimum
C
        IF ( TEST )
     &  WRITE (*,*) ' NEW MINIMUM FOUND BY CONTRACTION:', FSTST
        NSUCC = 0
      ELSE
C
C Increase counter for unsucessful contractions
C
        NSUCC = NSUCC + 1
      ENDIF
C
C Check if contracted point is better than maximum (go up)
C
      IF ( FSTST .LT. FHI ) GOTO 130
C
C After an unsuccessful contraction all points of
C the simplex are contracted towards the minimum.
C Set flag indicating the need for a full recomputation
C of the center of gravity in the next step
C
      AVGUPD = .FALSE.
C
      IF ( TEST )
     &   WRITE (*,*) ' SIMPLEX CONTRACTED'
C
      IMIN = 1
C
      DO 190 I = 2, NP
      DO 200 J = 1, N
200   SIM(J,I) = HALF * ( SIM(J,I) + SIM(J,IMIN) )
      F(I) = FU ( N, SIM(1,I),cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet )
C
C Check if contracted point is a new minimum,
C the remaining points will be moved towards the new target
C
      IF ( F(I) .LT. FLO ) THEN
        NSUCC = 0
        IF ( TEST )
     &  WRITE (*,*) ' MINIMUM EXCHANGED, I=',I, ' F=',F(I)
        FLO = F(I)
        IMIN = I
      ENDIF
C
190   CONTINUE
C
C Skip back to next iteration step
C ================================
C
      GOTO 35
C
C End of main iteration loop
C
999   RETURN
C
C End of subroutine FMINSI
C
      END
