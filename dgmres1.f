!
!    SLATEC: public domain
!
!    Change on 20180508 (eight of May 2018): matvec and msolve were
!                       removed from all argument lists, i.e. these
!                       surroutine names are fixed now (G. Dhondt)
!                       Needed for the multithreading parallellization
!
*DECK DGMRES1
      SUBROUTINE DGMRES1 (N, B, X, NELT, IA, JA, A, ISYM,
     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,
     +   IGWK, LIGW, RWORK, IWORK)
C***BEGIN PROLOGUE  DGMRES
C***PURPOSE  Preconditioned GMRES iterative sparse Ax=b solver.
C            This routine uses the generalized minimum residual
C            (GMRES) method with preconditioning to solve
C            non-symmetric linear systems of the form: Ax = b.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SGMRES-S, DGMRES-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C
C *Usage:
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C      INTEGER   ITER, IERR, IUNIT, LRGW, IGWK(LIGW), LIGW
C      INTEGER   IWORK(USER DEFINED)
C      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N)
C      DOUBLE PRECISION RGWK(LRGW), RWORK(USER DEFINED)
C      EXTERNAL  MATVEC, MSOLVE
C
C      CALL DGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX,
C     $     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK)
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for the solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description", below,
C         for more details.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MATVEC :EXT      External.
C         Name of a routine which performs the matrix vector multiply
C         Y = A*X given A and X.  The name of the MATVEC routine must
C         be declared external in the calling program.  The calling
C         sequence to MATVEC is:
C             CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)
C         where N is the number of unknowns, Y is the product A*X
C         upon return, X is an input vector, and NELT is the number of
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.
C         ISYM is a flag which, if non-zero, denotes that A is
C         symmetric and only the lower or upper triangle is stored.
C MSOLVE :EXT      External.
C         Name of the routine which solves a linear system Mz = r for
C         z given r with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays.  The name of the MSOLVE routine must
C         be declared external in the calling program.  The calling
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is the right-hand side
C         vector and Z is the solution upon return.  NELT, IA, JA, A and
C         ISYM are defined as above.  RWORK is a double precision array
C         that can be used to pass necessary preconditioning information
C         and/or workspace to MSOLVE.  IWORK is an integer work array
C         for the same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore  much less  efficient.  See  ISDGMR (the
C                 stop test routine) for more information.
C         ITOL=1  Means   the  iteration stops   when the first test
C                 described below on  the residual RL  is satisfied,
C                 and there  is either right  or  no preconditioning
C                 being used.
C         ITOL=2  Implies     that   the  user    is   using    left
C                 preconditioning, and the second stopping criterion
C                 below is used.
C         ITOL=3  Means the  iteration stops   when  the  third test
C                 described below on Minv*Residual is satisfied, and
C                 there is either left  or no  preconditioning being
C                 used.
C         ITOL=11 is    often  useful  for   checking  and comparing
C                 different routines.  For this case, the  user must
C                 supply  the  "exact" solution or  a  very accurate
C                 approximation (one with  an  error much less  than
C                 TOL) through a common block,
C                     COMMON /DSLBLK/ SOLN( )
C                 If ITOL=11, iteration stops when the 2-norm of the
C                 difference between the iterative approximation and
C                 the user-supplied solution  divided by the  2-norm
C                 of the  user-supplied solution  is  less than TOL.
C                 Note that this requires  the  user to  set up  the
C                 "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling
C                 routine.  The routine with this declaration should
C                 be loaded before the stop test so that the correct
C                 length is used by  the loader.  This procedure  is
C                 not standard Fortran and may not work correctly on
C                 your   system (although  it  has  worked  on every
C                 system the authors have tried).  If ITOL is not 11
C                 then this common block is indeed standard Fortran.
C TOL    :INOUT    Double Precision.
C         Convergence criterion, as described below.  If TOL is set
C         to zero on input, then a default value of 500*(the smallest
C         positive magnitude, machine epsilon) is used.
C ITMAX  :DUMMY    Integer.
C         Maximum number of iterations in most SLAP routines.  In
C         this routine this does not make sense.  The maximum number
C         of iterations here is given by ITMAX = MAXL*(NRMAX+1).
C         See IGWK for definitions of MAXL and NRMAX.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or
C         ITMAX if convergence criterion could not be achieved in
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows..
C
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               for right or no preconditioning, and
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               for left preconditioning.
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               since right or no preconditioning
C                               being used.
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               since left preconditioning is being
C                               used.
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
C                               i=1,n
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for
C                           RGWK or IGWK.
C               IERR = 2 => Routine DGMRES failed to reduce the norm
C                           of the current residual on its last call,
C                           and so the iteration has stalled.  In
C                           this case, X equals the last computed
C                           approximation.  The user must either
C                           increase MAXL, or choose a different
C                           initial guess.
C               IERR =-1 => Insufficient length for RGWK array.
C                           IGWK(6) contains the required minimum
C                           length of the RGWK array.
C               IERR =-2 => Illegal value of ITOL, or ITOL and JPRE
C                           values are inconsistent.
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
C         left-hand-side of the relevant stopping test defined
C         below associated with the residual for the current
C         approximation X(L).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration,
C         if this is desired for monitoring convergence.  If unit
C         number is 0, no writing will occur.
C SB     :IN       Double Precision SB(N).
C         Array of length N containing scale factors for the right
C         hand side vector B.  If JSCAL.eq.0 (see below), SB need
C         not be supplied.
C SX     :IN       Double Precision SX(N).
C         Array of length N containing scale factors for the solution
C         vector X.  If JSCAL.eq.0 (see below), SX need not be
C         supplied.  SB and SX can be the same array in the calling
C         program if desired.
C RGWK   :INOUT    Double Precision RGWK(LRGW).
C         Double Precision array used for workspace by DGMRES.
C         On return, RGWK(1) = RHOL.  See IERR for definition of RHOL.
C LRGW   :IN       Integer.
C         Length of the double precision workspace, RGWK.
C         LRGW >= 1 + N*(MAXL+6) + MAXL*(MAXL+3).
C         See below for definition of MAXL.
C         For the default values, RGWK has size at least 131 + 16*N.
C IGWK   :INOUT    Integer IGWK(LIGW).
C         The following IGWK parameters should be set by the user
C         before calling this routine.
C         IGWK(1) = MAXL.  Maximum dimension of Krylov subspace in
C            which X - X0 is to be found (where, X0 is the initial
C            guess).  The default value of MAXL is 10.
C         IGWK(2) = KMP.  Maximum number of previous Krylov basis
C            vectors to which each new basis vector is made orthogonal.
C            The default value of KMP is MAXL.
C         IGWK(3) = JSCAL.  Flag indicating whether the scaling
C            arrays SB and SX are to be used.
C            JSCAL = 0 => SB and SX are not used and the algorithm
C               will perform as if all SB(I) = 1 and SX(I) = 1.
C            JSCAL = 1 =>  Only SX is used, and the algorithm
C               performs as if all SB(I) = 1.
C            JSCAL = 2 =>  Only SB is used, and the algorithm
C               performs as if all SX(I) = 1.
C            JSCAL = 3 =>  Both SB and SX are used.
C         IGWK(4) = JPRE.  Flag indicating whether preconditioning
C            is being used.
C            JPRE = 0  =>  There is no preconditioning.
C            JPRE > 0  =>  There is preconditioning on the right
C               only, and the solver will call routine MSOLVE.
C            JPRE < 0  =>  There is preconditioning on the left
C               only, and the solver will call routine MSOLVE.
C         IGWK(5) = NRMAX.  Maximum number of restarts of the
C            Krylov iteration.  The default value of NRMAX = 10.
C            if IWORK(5) = -1,  then no restarts are performed (in
C            this case, NRMAX is set to zero internally).
C         The following IWORK parameters are diagnostic information
C         made available to the user after this routine completes.
C         IGWK(6) = MLWK.  Required minimum length of RGWK array.
C         IGWK(7) = NMS.  The total number of calls to MSOLVE.
C LIGW   :IN       Integer.
C         Length of the integer workspace, IGWK.  LIGW >= 20.
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).
C         Double Precision array that can be used for workspace in
C         MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINED).
C         Integer array that can be used for workspace in MSOLVE.
C
C *Description:
C       DGMRES solves a linear system A*X = B rewritten in the form:
C
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
C
C       with right preconditioning, or
C
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
C
C       with left preconditioning, where A is an N-by-N double precision
C       matrix, X and B are N-vectors,  SB and SX  are diagonal scaling
C       matrices,   and M is  a preconditioning    matrix.   It uses
C       preconditioned  Krylov   subpace  methods  based     on  the
C       generalized minimum residual  method (GMRES).   This routine
C       optionally performs  either  the  full     orthogonalization
C       version of the  GMRES  algorithm or an incomplete variant of
C       it.  Both versions use restarting of the linear iteration by
C       default, although the user can disable this feature.
C
C       The GMRES  algorithm generates a sequence  of approximations
C       X(L) to the  true solution of the above  linear system.  The
C       convergence criteria for stopping the  iteration is based on
C       the size  of the  scaled norm of  the residual  R(L)  =  B -
C       A*X(L).  The actual stopping test is either:
C
C               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),
C
C       for right preconditioning, or
C
C               norm(SB*(M-inverse)*(B-A*X(L))) .le.
C                       TOL*norm(SB*(M-inverse)*B),
C
C       for left preconditioning, where norm() denotes the Euclidean
C       norm, and TOL is  a positive scalar less  than one  input by
C       the user.  If TOL equals zero  when DGMRES is called, then a
C       default  value  of 500*(the   smallest  positive  magnitude,
C       machine epsilon) is used.  If the  scaling arrays SB  and SX
C       are used, then  ideally they  should be chosen  so  that the
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components
C       approximately equal  to  one in  magnitude.  If one wants to
C       use the same scaling in X  and B, then  SB and SX can be the
C       same array in the calling program.
C
C       The following is a list of the other routines and their
C       functions used by DGMRES:
C       DPIGMR  Contains the main iteration loop for GMRES.
C       DORTH   Orthogonalizes a new vector against older basis vectors.
C       DHEQR   Computes a QR decomposition of a Hessenberg matrix.
C       DHELS   Solves a Hessenberg least-squares system, using QR
C               factors.
C       DRLCAL  Computes the scaled residual RL.
C       DXLCAL  Computes the solution XL.
C       ISDGMR  User-replaceable stopping routine.
C
C       This routine does  not care  what matrix data   structure is
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE
C       routines, with  the arguments as  described above.  The user
C       could write any type of structure and the appropriate MATVEC
C       and MSOLVE routines.  It is assumed  that A is stored in the
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
C       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP
C       routines DSDCG and DSICCG are examples of this procedure.
C
C       Two  examples  of  matrix  data structures  are the: 1) SLAP
C       Triad  format and 2) SLAP Column format.
C
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       =================== S L A P Column format ==================
C
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, Reduced Storage
C                  Matrix Methods in Stiff ODE Systems, Lawrence Liver-
C                  more National Laboratory Report UCRL-95088, Rev. 1,
C                  Livermore, California, June 1987.
C               2. Mark K. Seager, A SLAP for the Masses, in
C                  G. F. Carey, Ed., Parallel Supercomputing: Methods,
C                  Algorithms and Applications, Wiley, 1989, pp.135-155.
C***ROUTINES CALLED  D1MACH, DCOPY, DNRM2, DPIGMR
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   891004  Added new reference.
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910506  Corrected errors in C***ROUTINES CALLED list.  (FNF)
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   920929  Corrected format of references.  (FNF)
C   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)
C   921026  Added check for valid value of ITOL.  (FNF)
C***END PROLOGUE  DGMRES
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      DOUBLE PRECISION ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LIGW, LRGW, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RGWK(LRGW), RWORK(*), SB(N),
     +                 SX(N), X(N)
      INTEGER IA(NELT), IGWK(LIGW), IWORK(*), JA(NELT)
C     .. Local Scalars ..
      DOUBLE PRECISION BNRM, RHOL, SUM
      INTEGER I, IFLAG, JPRE, JSCAL, KMP, LDL, LGMR, LHES, LQ, LR, LV,
     +        LW, LXL, LZ, LZM1, MAXL, MAXLP1, NMS, NMSL, NRMAX, NRSTS
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2
      EXTERNAL D1MACH, DNRM2
C     .. External Subroutines ..
      EXTERNAL DCOPY, DPIGMR1
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C***FIRST EXECUTABLE STATEMENT  DGMRES
      IERR = 0
C   ------------------------------------------------------------------
C         Load method parameters with user values or defaults.
C   ------------------------------------------------------------------
      MAXL = IGWK(1)
      IF (MAXL .EQ. 0) MAXL = 10
      IF (MAXL .GT. N) MAXL = N
      KMP = IGWK(2)
      IF (KMP .EQ. 0) KMP = MAXL
      IF (KMP .GT. MAXL) KMP = MAXL
      JSCAL = IGWK(3)
      JPRE = IGWK(4)
C         Check for valid value of ITOL.
      IF( (ITOL.LT.0) .OR. ((ITOL.GT.3).AND.(ITOL.NE.11)) ) GOTO 650
C         Check for consistent values of ITOL and JPRE.
      IF( ITOL.EQ.1 .AND. JPRE.LT.0 ) GOTO 650
      IF( ITOL.EQ.2 .AND. JPRE.GE.0 ) GOTO 650
      NRMAX = IGWK(5)
      IF( NRMAX.EQ.0 ) NRMAX = 10
C         If NRMAX .eq. -1, then set NRMAX = 0 to turn off restarting.
      IF( NRMAX.EQ.-1 ) NRMAX = 0
C         If input value of TOL is zero, set it to its default value.
      IF( TOL.EQ.0.0D0 ) TOL = 500*D1MACH(3)
C
C         Initialize counters.
      ITER = 0
      NMS = 0
      NRSTS = 0
C   ------------------------------------------------------------------
C         Form work array segment pointers.
C   ------------------------------------------------------------------
      MAXLP1 = MAXL + 1
      LV = 1
      LR = LV + N*MAXLP1
      LHES = LR + N + 1
      LQ = LHES + MAXL*MAXLP1
      LDL = LQ + 2*MAXL
      LW = LDL + N
      LXL = LW + N
      LZ = LXL + N
C
C         Load IGWK(6) with required minimum length of the RGWK array.
      IGWK(6) = LZ + N - 1
      IF( LZ+N-1.GT.LRGW ) GOTO 640
C   ------------------------------------------------------------------
C         Calculate scaled-preconditioned norm of RHS vector b.
C   ------------------------------------------------------------------
      IF (JPRE .LT. 0) THEN
         CALL MSOLVE(N, B, RGWK(LR), NELT, IA, JA, A, ISYM,
     $        RWORK, IWORK)
         NMS = NMS + 1
      ELSE
         CALL DCOPY(N, B, 1, RGWK(LR), 1)
      ENDIF
      IF( JSCAL.EQ.2 .OR. JSCAL.EQ.3 ) THEN
         SUM = 0
         DO 10 I = 1,N
            SUM = SUM + (RGWK(LR-1+I)*SB(I))**2
 10      CONTINUE
         BNRM = SQRT(SUM)
      ELSE
         BNRM = DNRM2(N,RGWK(LR),1)
      ENDIF
C   ------------------------------------------------------------------
C         Calculate initial residual.
C   ------------------------------------------------------------------
      CALL MATVEC(N, X, RGWK(LR), NELT, IA, JA, A, ISYM)
      DO 50 I = 1,N
         RGWK(LR-1+I) = B(I) - RGWK(LR-1+I)
 50   CONTINUE
C   ------------------------------------------------------------------
C         If performing restarting, then load the residual into the
C         correct location in the RGWK array.
C   ------------------------------------------------------------------
 100  CONTINUE
      IF( NRSTS.GT.NRMAX ) GOTO 610
      IF( NRSTS.GT.0 ) THEN
C         Copy the current residual to a different location in the RGWK
C         array.
         CALL DCOPY(N, RGWK(LDL), 1, RGWK(LR), 1)
      ENDIF
C   ------------------------------------------------------------------
C         Use the DPIGMR algorithm to solve the linear system A*Z = R.
C   ------------------------------------------------------------------
      CALL DPIGMR1(N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP,
     $       NRSTS, JPRE, NMSL, RGWK(LZ), RGWK(LV),
     $       RGWK(LHES), RGWK(LQ), LGMR, RWORK, IWORK, RGWK(LW),
     $       RGWK(LDL), RHOL, NRMAX, B, BNRM, X, RGWK(LXL), ITOL,
     $       TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
      ITER = ITER + LGMR
      NMS = NMS + NMSL
C
C         Increment X by the current approximate solution Z of A*Z = R.
C
      LZM1 = LZ - 1
      DO 110 I = 1,N
         X(I) = X(I) + RGWK(LZM1+I)
 110  CONTINUE
      IF( IFLAG.EQ.0 ) GOTO 600
      IF( IFLAG.EQ.1 ) THEN
         NRSTS = NRSTS + 1
         GOTO 100
      ENDIF
      IF( IFLAG.EQ.2 ) GOTO 620
C   ------------------------------------------------------------------
C         All returns are made through this section.
C   ------------------------------------------------------------------
C         The iteration has converged.
C
 600  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 0
      RETURN
C
C         Max number((NRMAX+1)*MAXL) of linear iterations performed.
 610  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 1
      RETURN
C
C         GMRES failed to reduce last residual in MAXL iterations.
C         The iteration has stalled.
 620  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 2
      RETURN
C         Error return.  Insufficient length for RGWK array.
 640  CONTINUE
      ERR = TOL
      IERR = -1
      RETURN
C         Error return.  Inconsistent ITOL and JPRE values.
 650  CONTINUE
      ERR = TOL
      IERR = -2
      RETURN
C------------- LAST LINE OF DGMRES FOLLOWS ----------------------------
      END
*DECK DPIGMR
      SUBROUTINE DPIGMR1(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, NRSTS,
     +   JPRE,  NMSL, Z, V, HES, Q, LGMR, RPAR, IPAR, WK,
     +   DL, RHOL, NRMAX, B, BNRM, X, XL, ITOL, TOL, NELT, IA, JA, A,
     +   ISYM, IUNIT, IFLAG, ERR)
C***BEGIN PROLOGUE  DPIGMR
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SPIGMR-S, DPIGMR-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C         This routine solves the linear system A * Z = R0 using a
C         scaled preconditioned version of the generalized minimum
C         residual method.  An initial guess of Z = 0 is assumed.
C
C *Usage:
C      INTEGER N, JSCAL, MAXL, MAXLP1, KMP, NRSTS, JPRE, NMSL, LGMR
C      INTEGER IPAR(USER DEFINED), NRMAX, ITOL, NELT, IA(NELT), JA(NELT)
C      INTEGER ISYM, IUNIT, IFLAG
C      DOUBLE PRECISION R0(N), SR(N), SZ(N), Z(N), V(N,MAXLP1),
C     $                 HES(MAXLP1,MAXL), Q(2*MAXL), RPAR(USER DEFINED),
C     $                 WK(N), DL(N), RHOL, B(N), BNRM, X(N), XL(N),
C     $                 TOL, A(NELT), ERR
C      EXTERNAL MATVEC, MSOLVE
C
C      CALL DPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP,
C     $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR,
C     $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL,
C     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
C
C *Arguments:
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C R0     :IN       Double Precision R0(N)
C         R0 = the right hand side of the system A*Z = R0.
C         R0 is also used as workspace when computing
C         the final approximation.
C         (R0 is the same as V(*,MAXL+1) in the call to DPIGMR.)
C SR     :IN       Double Precision SR(N)
C         SR is a vector of length N containing the non-zero
C         elements of the diagonal scaling matrix for R0.
C SZ     :IN       Double Precision SZ(N)
C         SZ is a vector of length N containing the non-zero
C         elements of the diagonal scaling matrix for Z.
C JSCAL  :IN       Integer
C         A flag indicating whether arrays SR and SZ are used.
C         JSCAL=0 means SR and SZ are not used and the
C                 algorithm will perform as if all
C                 SR(i) = 1 and SZ(i) = 1.
C         JSCAL=1 means only SZ is used, and the algorithm
C                 performs as if all SR(i) = 1.
C         JSCAL=2 means only SR is used, and the algorithm
C                 performs as if all SZ(i) = 1.
C         JSCAL=3 means both SR and SZ are used.
C MAXL   :IN       Integer
C         The maximum allowable order of the matrix H.
C MAXLP1 :IN       Integer
C         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to.  (KMP .le. MAXL)
C NRSTS  :IN       Integer
C         Counter for the number of restarts on the current
C         call to DGMRES.  If NRSTS .gt. 0, then the residual
C         R0 is already scaled, and so scaling of it is
C         not necessary.
C JPRE   :IN       Integer
C         Preconditioner type flag.
C MATVEC :EXT      External.
C         Name of a routine which performs the matrix vector multiply
C         Y = A*X given A and X.  The name of the MATVEC routine must
C         be declared external in the calling program.  The calling
C         sequence to MATVEC is:
C             CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)
C         where N is the number of unknowns, Y is the product A*X
C         upon return, X is an input vector, and NELT is the number of
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.
C         ISYM is a flag which, if non-zero, denotes that A is
C         symmetric and only the lower or upper triangle is stored.
C MSOLVE :EXT      External.
C         Name of the routine which solves a linear system Mz = r for
C         z given r with the preconditioning matrix M (M is supplied via
C         RPAR and IPAR arrays.  The name of the MSOLVE routine must
C         be declared external in the calling program.  The calling
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
C         Where N is the number of unknowns, R is the right-hand side
C         vector and Z is the solution upon return.  NELT, IA, JA, A and
C         ISYM are defined as below.  RPAR is a double precision array
C         that can be used to pass necessary preconditioning information
C         and/or workspace to MSOLVE.  IPAR is an integer work array
C         for the same purpose as RPAR.
C NMSL   :OUT      Integer
C         The number of calls to MSOLVE.
C Z      :OUT      Double Precision Z(N)
C         The final computed approximation to the solution
C         of the system A*Z = R0.
C V      :OUT      Double Precision V(N,MAXLP1)
C         The N by (LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C HES    :OUT      Double Precision HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,I)
C         and V(*,K).
C Q      :OUT      Double Precision Q(2*MAXL)
C         A double precision array of length 2*MAXL containing the
C         components of the Givens rotations used in the QR
C         decomposition of HES.  It is loaded in DHEQR and used in
C         DHELS.
C LGMR   :OUT      Integer
C         The number of iterations performed and
C         the current order of the upper Hessenberg
C         matrix HES.
C RPAR   :IN       Double Precision RPAR(USER DEFINED)
C         Double Precision workspace passed directly to the MSOLVE
C         routine.
C IPAR   :IN       Integer IPAR(USER DEFINED)
C         Integer workspace passed directly to the MSOLVE routine.
C WK     :IN       Double Precision WK(N)
C         A double precision work array of length N used by routines
C         MATVEC and MSOLVE.
C DL     :INOUT    Double Precision DL(N)
C         On input, a double precision work array of length N used for
C         calculation of the residual norm RHO when the method is
C         incomplete (KMP.lt.MAXL), and/or when using restarting.
C         On output, the scaled residual vector RL.  It is only loaded
C         when performing restarts of the Krylov iteration.
C RHOL   :OUT      Double Precision
C         A double precision scalar containing the norm of the final
C         residual.
C NRMAX  :IN       Integer
C         The maximum number of restarts of the Krylov iteration.
C         NRMAX .gt. 0 means restarting is active, while
C         NRMAX = 0 means restarting is not being used.
C B      :IN       Double Precision B(N)
C         The right hand side of the linear system A*X = b.
C BNRM   :IN       Double Precision
C         The scaled norm of b.
C X      :IN       Double Precision X(N)
C         The current approximate solution as of the last
C         restart.
C XL     :IN       Double Precision XL(N)
C         An array of length N used to hold the approximate
C         solution X(L) when ITOL=11.
C ITOL   :IN       Integer
C         A flag to indicate the type of convergence criterion
C         used.  See the driver for its description.
C TOL    :IN       Double Precision
C         The tolerance on residuals R0-A*Z in scaled norm.
C NELT   :IN       Integer
C         The length of arrays IA, JA and A.
C IA     :IN       Integer IA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C JA     :IN       Integer JA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C A      :IN       Double Precision A(NELT)
C         A double precision array of length NELT containing matrix
C         data. It is passed directly to the MATVEC and MSOLVE routines.
C ISYM   :IN       Integer
C         A flag to indicate symmetric matrix storage.
C         If ISYM=0, all non-zero entries of the matrix are
C         stored.  If ISYM=1, the matrix is symmetric and
C         only the upper or lower triangular part is stored.
C IUNIT  :IN       Integer
C         The i/o unit number for writing intermediate residual
C         norm values.
C IFLAG  :OUT      Integer
C         An integer error flag..
C         0 means convergence in LGMR iterations, LGMR.le.MAXL.
C         1 means the convergence test did not pass in MAXL
C           iterations, but the residual norm is .lt. norm(R0),
C           and so Z is computed.
C         2 means the convergence test did not pass in MAXL
C           iterations, residual .ge. norm(R0), and Z = 0.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as
C         defined by ITOL.
C
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DAXPY, DCOPY, DHELS, DHEQR, DNRM2, DORTH, DRLCAL,
C                    DSCAL, ISDGMR
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DPIGMR
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      DOUBLE PRECISION BNRM, ERR, RHOL, TOL
      INTEGER IFLAG, ISYM, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR, MAXL,
     +        MAXLP1, N, NELT, NMSL, NRMAX, NRSTS
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(*), DL(*), HES(MAXLP1,*), Q(*), R0(*),
     +                 RPAR(*), SR(*), SZ(*), V(N,*), WK(*), X(*),
     +                 XL(*), Z(*)
      INTEGER IA(NELT), IPAR(*), JA(NELT)
C     .. Local Scalars ..
      DOUBLE PRECISION C, DLNRM, PROD, R0NRM, RHO, S, SNORMW, TEM
      INTEGER I, I2, INFO, IP1, ITER, ITMAX, J, K, LL, LLP1
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER ISDGMR1
      EXTERNAL DNRM2, ISDGMR1
C     .. External Subroutines ..
      EXTERNAL DAXPY, DCOPY, DHELS, DHEQR, DORTH, DRLCAL, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C***FIRST EXECUTABLE STATEMENT  DPIGMR
C
C         Zero out the Z array.
C
      DO 5 I = 1,N
         Z(I) = 0
 5    CONTINUE
C
      IFLAG = 0
      LGMR = 0
      NMSL = 0
C         Load ITMAX, the maximum number of iterations.
      ITMAX =(NRMAX+1)*MAXL
C   -------------------------------------------------------------------
C         The initial residual is the vector R0.
C         Apply left precon. if JPRE < 0 and this is not a restart.
C         Apply scaling to R0 if JSCAL = 2 or 3.
C   -------------------------------------------------------------------
      IF ((JPRE .LT. 0) .AND.(NRSTS .EQ. 0)) THEN
         CALL DCOPY(N, R0, 1, WK, 1)
         CALL MSOLVE(N, WK, R0, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      IF (((JSCAL.EQ.2) .OR.(JSCAL.EQ.3)) .AND.(NRSTS.EQ.0)) THEN
         DO 10 I = 1,N
            V(I,1) = R0(I)*SR(I)
 10      CONTINUE
      ELSE
         DO 20 I = 1,N
            V(I,1) = R0(I)
 20      CONTINUE
      ENDIF
      R0NRM = DNRM2(N, V, 1)
      ITER = NRSTS*MAXL
C
C         Call stopping routine ISDGMR.
C
      IF (ISDGMR1(N, B, X, XL, NELT, IA, JA, A, ISYM, 
     $    NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, V(1,1), Z, WK,
     $    RPAR, IPAR, R0NRM, BNRM, SR, SZ, JSCAL,
     $    KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
     $    HES, JPRE) .NE. 0) RETURN
      TEM = 1.0D0/R0NRM
      CALL DSCAL(N, TEM, V(1,1), 1)
C
C         Zero out the HES array.
C
      DO 50 J = 1,MAXL
         DO 40 I = 1,MAXLP1
            HES(I,J) = 0
 40      CONTINUE
 50   CONTINUE
C   -------------------------------------------------------------------
C         Main loop to compute the vectors V(*,2) to V(*,MAXL).
C         The running product PROD is needed for the convergence test.
C   -------------------------------------------------------------------
      PROD = 1
      DO 90 LL = 1,MAXL
         LGMR = LL
C   -------------------------------------------------------------------
C        Unscale  the  current V(LL)  and store  in WK.  Call routine
C        MSOLVE    to   compute(M-inverse)*WK,   where    M   is  the
C        preconditioner matrix.  Save the answer in Z.   Call routine
C        MATVEC to compute  VNEW  = A*Z,  where  A is  the the system
C        matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call
C        routine DORTH  to  orthogonalize the    new vector VNEW   =
C        V(*,LL+1).  Call routine DHEQR to update the factors of HES.
C   -------------------------------------------------------------------
        IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
           DO 60 I = 1,N
              WK(I) = V(I,LL)/SZ(I)
 60        CONTINUE
        ELSE
           CALL DCOPY(N, V(1,LL), 1, WK, 1)
        ENDIF
        IF (JPRE .GT. 0) THEN
           CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
           NMSL = NMSL + 1
           CALL MATVEC(N, Z, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ELSE
           CALL MATVEC(N, WK, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ENDIF
        IF (JPRE .LT. 0) THEN
           CALL DCOPY(N, V(1,LL+1), 1, WK, 1)
           CALL MSOLVE(N,WK,V(1,LL+1),NELT,IA,JA,A,ISYM,RPAR,IPAR)
           NMSL = NMSL + 1
        ENDIF
        IF ((JSCAL .EQ. 2) .OR.(JSCAL .EQ. 3)) THEN
           DO 65 I = 1,N
              V(I,LL+1) = V(I,LL+1)*SR(I)
 65        CONTINUE
        ENDIF
        CALL DORTH(V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)
        HES(LL+1,LL) = SNORMW
        CALL DHEQR(HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
C   -------------------------------------------------------------------
C         Update RHO, the estimate of the norm of the residual R0-A*ZL.
C         If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
C         necessarily orthogonal for LL > KMP.  The vector DL must then
C         be computed, and its norm used in the calculation of RHO.
C   -------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*R0NRM)
        IF ((LL.GT.KMP) .AND.(KMP.LT.MAXL)) THEN
           IF (LL .EQ. KMP+1) THEN
              CALL DCOPY(N, V(1,1), 1, DL, 1)
              DO 75 I = 1,KMP
                 IP1 = I + 1
                 I2 = I*2
                 S = Q(I2)
                 C = Q(I2-1)
                 DO 70 K = 1,N
                    DL(K) = S*DL(K) + C*V(K,IP1)
 70              CONTINUE
 75           CONTINUE
           ENDIF
           S = Q(2*LL)
           C = Q(2*LL-1)/SNORMW
           LLP1 = LL + 1
           DO 80 K = 1,N
              DL(K) = S*DL(K) + C*V(K,LLP1)
 80        CONTINUE
           DLNRM = DNRM2(N, DL, 1)
           RHO = RHO*DLNRM
        ENDIF
        RHOL = RHO
C   -------------------------------------------------------------------
C         Test for convergence.  If passed, compute approximation ZL.
C         If failed and LL < MAXL, then continue iterating.
C   -------------------------------------------------------------------
        ITER = NRSTS*MAXL + LGMR
        IF (ISDGMR1(N, B, X, XL, NELT, IA, JA, A, ISYM, 
     $      NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, DL, Z, WK,
     $      RPAR, IPAR, RHOL, BNRM, SR, SZ, JSCAL,
     $      KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
     $      HES, JPRE) .NE. 0) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
C   -------------------------------------------------------------------
C         Rescale so that the norm of V(1,LL+1) is one.
C   -------------------------------------------------------------------
        TEM = 1.0D0/SNORMW
        CALL DSCAL(N, TEM, V(1,LL+1), 1)
 90   CONTINUE
 100  CONTINUE
      IF (RHO .LT. R0NRM) GO TO 150
 120  CONTINUE
      IFLAG = 2
C
C         Load approximate solution with zero.
C
      DO 130 I = 1,N
         Z(I) = 0
 130  CONTINUE
      RETURN
 150  IFLAG = 1
C
C         Tolerance not met, but residual norm reduced.
C
      IF (NRMAX .GT. 0) THEN
C
C        If performing restarting (NRMAX > 0)  calculate the residual
C        vector RL and  store it in the DL  array.  If the incomplete
C        version is being used (KMP < MAXL) then DL has  already been
C        calculated up to a scaling factor.   Use DRLCAL to calculate
C        the scaled residual vector.
C
         CALL DRLCAL(N, KMP, MAXL, MAXL, V, Q, DL, SNORMW, PROD,
     $        R0NRM)
      ENDIF
C   -------------------------------------------------------------------
C         Compute the approximation ZL to the solution.  Since the
C         vector Z was used as workspace, and the initial guess
C         of the linear iteration is zero, Z must be reset to zero.
C   -------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
         R0(K) = 0
 210  CONTINUE
      R0(1) = R0NRM
      CALL DHELS(HES, MAXLP1, LL, Q, R0)
      DO 220 K = 1,N
         Z(K) = 0
 220  CONTINUE
      DO 230 I = 1,LL
         CALL DAXPY(N, R0(I), V(1,I), 1, Z, 1)
 230  CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 240 I = 1,N
            Z(I) = Z(I)/SZ(I)
 240     CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL DCOPY(N, Z, 1, WK, 1)
         CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      RETURN
C------------- LAST LINE OF DPIGMR FOLLOWS ----------------------------
      END
*DECK ISDGMR
      INTEGER FUNCTION ISDGMR1 (N, B, X, XL, NELT, IA, JA, A, ISYM,
     +   NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
     +   RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL, KMP, LGMR, MAXL,
     +   MAXLP1, V, Q, SNORMW, PROD, R0NRM, HES, JPRE)
C***BEGIN PROLOGUE  ISDGMR
C***SUBSIDIARY
C***PURPOSE  Generalized Minimum Residual Stop Test.
C            This routine calculates the stop test for the Generalized
C            Minimum RESidual (GMRES) iteration scheme.  It returns a
C            non-zero if the error estimate (the type of which is
C            determined by ITOL) is less than the user specified
C            tolerance TOL.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (ISSGMR-S, ISDGMR-D)
C***KEYWORDS  GMRES, LINEAR SYSTEM, SLAP, SPARSE, STOP TEST
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C
C *Usage:
C      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NMSL, ITOL
C      INTEGER ITMAX, ITER, IUNIT, IWORK(USER DEFINED), JSCAL
C      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE
C      DOUBLE PRECISION B(N), X(N), XL(MAXL), A(NELT), TOL, ERR,
C     $                 R(N), Z(N), DZ(N), RWORK(USER DEFINED),
C     $                 RNRM, BNRM, SB(N), SX(N), V(N,MAXLP1),
C     $                 Q(2*MAXL), SNORMW, PROD, R0NRM,
C     $                 HES(MAXLP1,MAXL)
C      EXTERNAL MSOLVE
C
C      IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
C     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
C     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL,
C     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
C     $     HES, JPRE) .NE. 0) THEN ITERATION DONE
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand-side vector.
C X      :IN       Double Precision X(N).
C         Approximate solution vector as of the last restart.
C XL     :OUT      Double Precision XL(N)
C         An array of length N used to hold the approximate
C         solution as of the current iteration.  Only computed by
C         this routine when ITOL=11.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description", in the DGMRES,
C         DSLUGM and DSDGMR routines for more details.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system Mz = r for  z
C         given r with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays.  The name of the MSOLVE routine must
C         be declared external in the calling program.  The calling
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is the right-hand side
C         vector and Z is the solution upon return.  NELT, IA, JA, A and
C         ISYM are defined as above.  RWORK is a double precision array
C         that can be used to pass necessary preconditioning information
C         and/or workspace to MSOLVE.  IWORK is an integer work array
C         for the same purpose as RWORK.
C NMSL   :INOUT    Integer.
C         A counter for the number of calls to MSOLVE.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore much less efficient.
C         ITOL=1  Means   the  iteration stops   when the first test
C                 described below on  the residual RL  is satisfied,
C                 and there  is either right  or  no preconditioning
C                 being used.
C         ITOL=2  Implies     that   the  user    is   using    left
C                 preconditioning, and the second stopping criterion
C                 below is used.
C         ITOL=3  Means the  iteration stops   when  the  third test
C                 described below on Minv*Residual is satisfied, and
C                 there is either left  or no  preconditioning begin
C                 used.
C         ITOL=11 is    often  useful  for   checking  and comparing
C                 different routines.  For this case, the  user must
C                 supply  the  "exact" solution or  a  very accurate
C                 approximation (one with  an  error much less  than
C                 TOL) through a common block,
C                     COMMON /DSLBLK/ SOLN( )
C                 If ITOL=11, iteration stops when the 2-norm of the
C                 difference between the iterative approximation and
C                 the user-supplied solution  divided by the  2-norm
C                 of the  user-supplied solution  is  less than TOL.
C                 Note that this requires  the  user to  set up  the
C                 "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling
C                 routine.  The routine with this declaration should
C                 be loaded before the stop test so that the correct
C                 length is used by  the loader.  This procedure  is
C                 not standard Fortran and may not work correctly on
C                 your   system (although  it  has  worked  on every
C                 system the authors have tried).  If ITOL is not 11
C                 then this common block is indeed standard Fortran.
C TOL    :IN       Double Precision.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :IN       Integer.
C         The iteration for which to check for convergence.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows..
C
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               for right or no preconditioning, and
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               for left preconditioning.
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               since right or no preconditioning
C                               being used.
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               since left preconditioning is being
C                               used.
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
C                               i=1,n
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration,
C         if this is desired for monitoring convergence.  If unit
C         number is 0, no writing will occur.
C R      :INOUT    Double Precision R(N).
C         Work array used in calling routine.  It contains
C         information necessary to compute the residual RL = B-A*XL.
C Z      :WORK     Double Precision Z(N).
C         Workspace used to hold the pseudo-residual M z = r.
C DZ     :WORK     Double Precision DZ(N).
C         Workspace used to hold temporary vector(s).
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).
C         Double Precision array that can be used by MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINED).
C         Integer array that can be used by MSOLVE.
C RNRM   :IN       Double Precision.
C         Norm of the current residual.  Type of norm depends on ITOL.
C BNRM   :IN       Double Precision.
C         Norm of the right hand side.  Type of norm depends on ITOL.
C SB     :IN       Double Precision SB(N).
C         Scaling vector for B.
C SX     :IN       Double Precision SX(N).
C         Scaling vector for X.
C JSCAL  :IN       Integer.
C         Flag indicating if scaling arrays SB and SX are being
C         used in the calling routine DPIGMR.
C         JSCAL=0 means SB and SX are not used and the
C                 algorithm will perform as if all
C                 SB(i) = 1 and SX(i) = 1.
C         JSCAL=1 means only SX is used, and the algorithm
C                 performs as if all SB(i) = 1.
C         JSCAL=2 means only SB is used, and the algorithm
C                 performs as if all SX(i) = 1.
C         JSCAL=3 means both SB and SX are used.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to.  (KMP .le. MAXL)
C LGMR   :IN       Integer
C         The number of GMRES iterations performed on the current call
C         to DPIGMR (i.e., # iterations since the last restart) and
C         the current order of the upper Hessenberg
C         matrix HES.
C MAXL   :IN       Integer
C         The maximum allowable order of the matrix H.
C MAXLP1 :IN       Integer
C         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
C V      :IN       Double Precision V(N,MAXLP1)
C         The N by (LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C Q      :IN       Double Precision Q(2*MAXL)
C         A double precision array of length 2*MAXL containing the
C         components of the Givens rotations used in the QR
C         decomposition of HES.
C SNORMW :IN       Double Precision
C         A scalar containing the scaled norm of VNEW before it
C         is renormalized in DPIGMR.
C PROD   :IN       Double Precision
C         The product s1*s2*...*sl = the product of the sines of the
C         Givens rotations used in the QR factorization of the
C         Hessenberg matrix HES.
C R0NRM  :IN       Double Precision
C         The scaled norm of initial residual R0.
C HES    :IN       Double Precision HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,I)
C         and V(*,K).
C JPRE   :IN       Integer
C         Preconditioner type flag.
C         (See description of IGWK(4) in DGMRES.)
C
C *Description
C       When using the GMRES solver,  the preferred value  for ITOL
C       is 0.  This is due to the fact that when ITOL=0 the norm of
C       the residual required in the stopping test is  obtained for
C       free, since this value is already  calculated  in the GMRES
C       algorithm.   The  variable  RNRM contains the   appropriate
C       norm, which is equal to norm(SB*(RL - A*XL))  when right or
C       no   preconditioning is  being  performed,   and equal   to
C       norm(SB*Minv*(RL - A*XL))  when using left preconditioning.
C       Here, norm() is the Euclidean norm.  Nonzero values of ITOL
C       require  additional work  to  calculate the  actual  scaled
C       residual  or its scaled/preconditioned  form,  and/or   the
C       approximate solution XL.  Hence, these values of  ITOL will
C       not be as efficient as ITOL=0.
C
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C     This routine does not verify that ITOL has a valid value.
C     The calling routine should make such a test before calling
C     ISDGMR, as is done in DGMRES.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  D1MACH, DCOPY, DNRM2, DRLCAL, DSCAL, DXLCAL
C***COMMON BLOCKS    DSLBLK
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Corrected conversion errors, etc.  (FNF)
C   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   921026  Corrected D to E in output format.  (FNF)
C   921113  Corrected C***CATEGORY line.  (FNF)
C***END PROLOGUE  ISDGMR
C     .. Scalar Arguments ..
      DOUBLE PRECISION BNRM, ERR, PROD, R0NRM, RNRM, SNORMW, TOL
      INTEGER ISYM, ITER, ITMAX, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR,
     +        MAXL, MAXLP1, N, NELT, NMSL
C     .. Array Arguments ..
      DOUBLE PRECISION A(*), B(*), DZ(*), HES(MAXLP1, MAXL), Q(*), R(*),
     +                 RWORK(*), SB(*), SX(*), V(N,*), X(*), XL(*), Z(*)
      INTEGER IA(*), IWORK(*), JA(*)
C     .. Arrays in Common ..
      DOUBLE PRECISION SOLN(1)
C     .. Local Scalars ..
      DOUBLE PRECISION DXNRM, FUZZ, RAT, RATMAX, SOLNRM, TEM
      INTEGER I, IELMAX
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2
      EXTERNAL D1MACH, DNRM2
C     .. External Subroutines ..
      EXTERNAL DCOPY, DRLCAL, DSCAL, DXLCAL1
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON /DSLBLK/ SOLN
C     .. Save statement ..
      SAVE SOLNRM
C***FIRST EXECUTABLE STATEMENT  ISDGMR
      ISDGMR1 = 0
      IF ( ITOL.EQ.0 ) THEN
C
C       Use input from DPIGMR to determine if stop conditions are met.
C
         ERR = RNRM/BNRM
      ENDIF
      IF ( (ITOL.GT.0) .AND. (ITOL.LE.3) ) THEN
C
C       Use DRLCAL to calculate the scaled residual vector.
C       Store answer in R.
C
         IF ( LGMR.NE.0 ) CALL DRLCAL(N, KMP, LGMR, MAXL, V, Q, R,
     $                                SNORMW, PROD, R0NRM)
         IF ( ITOL.LE.2 ) THEN
C         err = ||Residual||/||RightHandSide||(2-Norms).
            ERR = DNRM2(N, R, 1)/BNRM
C
C         Unscale R by R0NRM*PROD when KMP < MAXL.
C
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0D0/(R0NRM*PROD)
               CALL DSCAL(N, TEM, R, 1)
            ENDIF
         ELSEIF ( ITOL.EQ.3 ) THEN
C         err = Max |(Minv*Residual)(i)/x(i)|
C         When JPRE .lt. 0, R already contains Minv*Residual.
            IF ( JPRE.GT.0 ) THEN
               CALL MSOLVE(N, R, DZ, NELT, IA, JA, A, ISYM, RWORK,
     $              IWORK)
               NMSL = NMSL + 1
            ENDIF
C
C         Unscale R by R0NRM*PROD when KMP < MAXL.
C
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0D0/(R0NRM*PROD)
               CALL DSCAL(N, TEM, R, 1)
            ENDIF
C
            FUZZ = D1MACH(1)
            IELMAX = 1
            RATMAX = ABS(DZ(1))/MAX(ABS(X(1)),FUZZ)
            DO 25 I = 2, N
               RAT = ABS(DZ(I))/MAX(ABS(X(I)),FUZZ)
               IF( RAT.GT.RATMAX ) THEN
                  IELMAX = I
                  RATMAX = RAT
               ENDIF
 25         CONTINUE
            ERR = RATMAX
            IF( RATMAX.LE.TOL ) ISDGMR1 = 1
            IF( IUNIT.GT.0 ) WRITE(IUNIT,1020) ITER, IELMAX, RATMAX
            RETURN
         ENDIF
      ENDIF
      IF ( ITOL.EQ.11 ) THEN
C
C       Use DXLCAL to calculate the approximate solution XL.
C
         IF ( (LGMR.NE.0) .AND. (ITER.GT.0) ) THEN
            CALL DXLCAL1(N, LGMR, X, XL, XL, HES, MAXLP1, Q, V, R0NRM,
     $           DZ, SX, JSCAL, JPRE, NMSL, RWORK, IWORK,
     $           NELT, IA, JA, A, ISYM)
         ELSEIF ( ITER.EQ.0 ) THEN
C         Copy X to XL to check if initial guess is good enough.
            CALL DCOPY(N, X, 1, XL, 1)
         ELSE
C         Return since this is the first call to DPIGMR on a restart.
            RETURN
         ENDIF
C
         IF ((JSCAL .EQ. 0) .OR.(JSCAL .EQ. 2)) THEN
C         err = ||x-TrueSolution||/||TrueSolution||(2-Norms).
            IF ( ITER.EQ.0 ) SOLNRM = DNRM2(N, SOLN, 1)
            DO 30 I = 1, N
               DZ(I) = XL(I) - SOLN(I)
 30         CONTINUE
            ERR = DNRM2(N, DZ, 1)/SOLNRM
         ELSE
            IF (ITER .EQ. 0) THEN
               SOLNRM = 0
               DO 40 I = 1,N
                  SOLNRM = SOLNRM + (SX(I)*SOLN(I))**2
 40            CONTINUE
               SOLNRM = SQRT(SOLNRM)
            ENDIF
            DXNRM = 0
            DO 50 I = 1,N
               DXNRM = DXNRM + (SX(I)*(XL(I)-SOLN(I)))**2
 50         CONTINUE
            DXNRM = SQRT(DXNRM)
C         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms).
            ERR = DXNRM/SOLNRM
         ENDIF
      ENDIF
C
      IF( IUNIT.NE.0 ) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL, MAXL, KMP
         ENDIF
         WRITE(IUNIT,1010) ITER, RNRM/BNRM, ERR
      ENDIF
      IF ( ERR.LE.TOL ) ISDGMR1 = 1
C
      RETURN
 1000 FORMAT(' Generalized Minimum Residual(',I3,I3,') for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Natural Err Est','   Error Estimate')
 1010 FORMAT(1X,I4,1X,D16.7,1X,D16.7)
 1020 FORMAT(1X,' ITER = ',I5, ' IELMAX = ',I5,
     $     ' |R(IELMAX)/X(IELMAX)| = ',D12.5)
C------------- LAST LINE OF ISDGMR FOLLOWS ----------------------------
      END
*DECK DXLCAL
      SUBROUTINE DXLCAL1 (N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
     +   WK, SZ, JSCAL, JPRE, NMSL, RPAR, IPAR, NELT, IA, JA, A,
     +   ISYM)
C***BEGIN PROLOGUE  DXLCAL
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SXLCAL-S, DXLCAL-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C        This  routine computes the solution  XL,  the current DGMRES
C        iterate, given the  V(I)'s and  the  QR factorization of the
C        Hessenberg  matrix HES.   This routine  is  only called when
C        ITOL=11.
C
C *Usage:
C      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, NMSL, IPAR(USER DEFINED)
C      INTEGER NELT, IA(NELT), JA(NELT), ISYM
C      DOUBLE PRECISION X(N), XL(N), ZL(N), HES(MAXLP1,MAXL), Q(2*MAXL),
C     $                 V(N,MAXLP1), R0NRM, WK(N), SZ(N),
C     $                 RPAR(USER DEFINED), A(NELT)
C      EXTERNAL MSOLVE
C
C      CALL DXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
C     $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR,
C     $     NELT, IA, JA, A, ISYM)
C
C *Arguments:
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C LGMR   :IN       Integer
C         The number of iterations performed and
C         the current order of the upper Hessenberg
C         matrix HES.
C X      :IN       Double Precision X(N)
C         The current approximate solution as of the last restart.
C XL     :OUT      Double Precision XL(N)
C         An array of length N used to hold the approximate
C         solution X(L).
C         Warning: XL and ZL are the same array in the calling routine.
C ZL     :IN       Double Precision ZL(N)
C         An array of length N used to hold the approximate
C         solution Z(L).
C HES    :IN       Double Precision HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,i) and V(*,k).
C MAXLP1 :IN       Integer
C         MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
C         MAXL is the maximum allowable order of the matrix HES.
C Q      :IN       Double Precision Q(2*MAXL)
C         A double precision array of length 2*MAXL containing the
C         components of the Givens rotations used in the QR
C         decomposition of HES.  It is loaded in DHEQR.
C V      :IN       Double Precision V(N,MAXLP1)
C         The N by(LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C R0NRM  :IN       Double Precision
C         The scaled norm of the initial residual for the
C         current call to DPIGMR.
C WK     :IN       Double Precision WK(N)
C         A double precision work array of length N.
C SZ     :IN       Double Precision SZ(N)
C         A vector of length N containing the non-zero
C         elements of the diagonal scaling matrix for Z.
C JSCAL  :IN       Integer
C         A flag indicating whether arrays SR and SZ are used.
C         JSCAL=0 means SR and SZ are not used and the
C                 algorithm will perform as if all
C                 SR(i) = 1 and SZ(i) = 1.
C         JSCAL=1 means only SZ is used, and the algorithm
C                 performs as if all SR(i) = 1.
C         JSCAL=2 means only SR is used, and the algorithm
C                 performs as if all SZ(i) = 1.
C         JSCAL=3 means both SR and SZ are used.
C JPRE   :IN       Integer
C         The preconditioner type flag.
C MSOLVE :EXT      External.
C         Name of the routine which solves a linear system Mz = r for
C         z given r with the preconditioning matrix M (M is supplied via
C         RPAR and IPAR arrays.  The name of the MSOLVE routine must
C         be declared external in the calling program.  The calling
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
C         Where N is the number of unknowns, R is the right-hand side
C         vector and Z is the solution upon return.  NELT, IA, JA, A and
C         ISYM are defined as below.  RPAR is a double precision array
C         that can be used to pass necessary preconditioning information
C         and/or workspace to MSOLVE.  IPAR is an integer work array
C         for the same purpose as RPAR.
C NMSL   :IN       Integer
C         The number of calls to MSOLVE.
C RPAR   :IN       Double Precision RPAR(USER DEFINED)
C         Double Precision workspace passed directly to the MSOLVE
C         routine.
C IPAR   :IN       Integer IPAR(USER DEFINED)
C         Integer workspace passed directly to the MSOLVE routine.
C NELT   :IN       Integer
C         The length of arrays IA, JA and A.
C IA     :IN       Integer IA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C JA     :IN       Integer JA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C A      :IN       Double Precision A(NELT)
C         A double precision array of length NELT containing matrix
C         data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C ISYM   :IN       Integer
C         A flag to indicate symmetric matrix storage.
C         If ISYM=0, all non-zero entries of the matrix are
C         stored.  If ISYM=1, the matrix is symmetric and
C         only the upper or lower triangular part is stored.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DAXPY, DCOPY, DHELS
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DXLCAL
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      DOUBLE PRECISION R0NRM
      INTEGER ISYM, JPRE, JSCAL, LGMR, MAXLP1, N, NELT, NMSL
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), HES(MAXLP1,*), Q(*), RPAR(*), SZ(*),
     +                 V(N,*), WK(N), X(N), XL(N), ZL(N)
      INTEGER IA(NELT), IPAR(*), JA(NELT)
C     .. Local Scalars ..
      INTEGER I, K, LL, LLP1
C     .. External Subroutines ..
      EXTERNAL DAXPY, DCOPY, DHELS
C***FIRST EXECUTABLE STATEMENT  DXLCAL
      LL = LGMR
      LLP1 = LL + 1
      DO 10 K = 1,LLP1
         WK(K) = 0
 10   CONTINUE
      WK(1) = R0NRM
      CALL DHELS(HES, MAXLP1, LL, Q, WK)
      DO 20 K = 1,N
         ZL(K) = 0
 20   CONTINUE
      DO 30 I = 1,LL
         CALL DAXPY(N, WK(I), V(1,I), 1, ZL, 1)
 30   CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 40 K = 1,N
            ZL(K) = ZL(K)/SZ(K)
 40      CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL DCOPY(N, ZL, 1, WK, 1)
         CALL MSOLVE(N, WK, ZL, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
C         calculate XL from X and ZL.
      DO 50 K = 1,N
         XL(K) = X(K) + ZL(K)
 50   CONTINUE
      RETURN
C------------- LAST LINE OF DXLCAL FOLLOWS ----------------------------
      END
