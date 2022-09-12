!
!    SLATEC: public domain
!
*DECK DGMRES
      SUBROUTINE DGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
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
C     .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
C     .. Local Scalars ..
      DOUBLE PRECISION BNRM, RHOL, SUM
      INTEGER I, IFLAG, JPRE, JSCAL, KMP, LDL, LGMR, LHES, LQ, LR, LV,
     +        LW, LXL, LZ, LZM1, MAXL, MAXLP1, NMS, NMSL, NRMAX, NRSTS
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2
      EXTERNAL D1MACH, DNRM2
C     .. External Subroutines ..
      EXTERNAL DCOPY, DPIGMR
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
      CALL DPIGMR(N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP,
     $       NRSTS, JPRE, MATVEC, MSOLVE, NMSL, RGWK(LZ), RGWK(LV),
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
*DECK DNRM2
      DOUBLE PRECISION FUNCTION DNRM2 (N, DX, INCX)
C***BEGIN PROLOGUE  DNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3B
C***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX with storage
C     increment INCX.
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four phase method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm.
C
C     Phase 1 scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNRM2
      INTEGER NEXT
      DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO,
     +                 ONE
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
C
      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C***FIRST EXECUTABLE STATEMENT  DNRM2
      IF (N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C
C                                                 BEGIN MAIN LOOP
C
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
C
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = ABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (ABS(DX(I)) .GT. CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (ABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = ABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI / N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J = I,NN,INCX
      IF (ABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = SQRT(SUM)
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END
*DECK DPIGMR
      SUBROUTINE DPIGMR (N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, NRSTS,
     +   JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR, RPAR, IPAR, WK,
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
C     .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
C     .. Local Scalars ..
      DOUBLE PRECISION C, DLNRM, PROD, R0NRM, RHO, S, SNORMW, TEM
      INTEGER I, I2, INFO, IP1, ITER, ITMAX, J, K, LL, LLP1
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER ISDGMR
      EXTERNAL DNRM2, ISDGMR
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
      IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
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
        IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
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
*DECK DRLCAL
      SUBROUTINE DRLCAL (N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD,
     +   R0NRM)
C***BEGIN PROLOGUE  DRLCAL
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SRLCAL-S, DRLCAL-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C         This routine calculates the scaled residual RL from the
C         V(I)'s.
C *Usage:
C      INTEGER N, KMP, LL, MAXL
C      DOUBLE PRECISION V(N,LL), Q(2*MAXL), RL(N), SNORMW, PROD, R0NORM
C
C      CALL DRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD, R0NRM)
C
C *Arguments:
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C KMP    :IN       Integer
C         The number of previous V vectors the new vector VNEW
C         must be made orthogonal to. (KMP .le. MAXL)
C LL     :IN       Integer
C         The current dimension of the Krylov subspace.
C MAXL   :IN       Integer
C         The maximum dimension of the Krylov subspace.
C V      :IN       Double Precision V(N,LL)
C         The N x LL array containing the orthogonal vectors
C         V(*,1) to V(*,LL).
C Q      :IN       Double Precision Q(2*MAXL)
C         A double precision array of length 2*MAXL containing the
C         components of the Givens rotations used in the QR
C         decomposition of HES.  It is loaded in DHEQR and used in
C         DHELS.
C RL     :OUT      Double Precision RL(N)
C         The residual vector RL.  This is either SB*(B-A*XL) if
C         not preconditioning or preconditioning on the right,
C         or SB*(M-inverse)*(B-A*XL) if preconditioning on the
C         left.
C SNORMW :IN       Double Precision
C         Scale factor.
C PROD   :IN       Double Precision
C         The product s1*s2*...*sl = the product of the sines of the
C         Givens rotations used in the QR factorization of
C         the Hessenberg matrix HES.
C R0NRM  :IN       Double Precision
C         The scaled norm of initial residual R0.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DCOPY, DSCAL
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DRLCAL
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      DOUBLE PRECISION PROD, R0NRM, SNORMW
      INTEGER KMP, LL, MAXL, N
C     .. Array Arguments ..
      DOUBLE PRECISION Q(*), RL(N), V(N,*)
C     .. Local Scalars ..
      DOUBLE PRECISION C, S, TEM
      INTEGER I, I2, IP1, K, LLM1, LLP1
C     .. External Subroutines ..
      EXTERNAL DCOPY, DSCAL
C***FIRST EXECUTABLE STATEMENT  DRLCAL
      IF (KMP .EQ. MAXL) THEN
C
C         calculate RL.  Start by copying V(*,1) into RL.
C
         CALL DCOPY(N, V(1,1), 1, RL, 1)
         LLM1 = LL - 1
         DO 20 I = 1,LLM1
            IP1 = I + 1
            I2 = I*2
            S = Q(I2)
            C = Q(I2-1)
            DO 10 K = 1,N
               RL(K) = S*RL(K) + C*V(K,IP1)
 10         CONTINUE
 20      CONTINUE
         S = Q(2*LL)
         C = Q(2*LL-1)/SNORMW
         LLP1 = LL + 1
         DO 30 K = 1,N
            RL(K) = S*RL(K) + C*V(K,LLP1)
 30      CONTINUE
      ENDIF
C
C         When KMP < MAXL, RL vector already partially calculated.
C         Scale RL by R0NRM*PROD to obtain the residual RL.
C
      TEM = R0NRM*PROD
      CALL DSCAL(N, TEM, RL, 1)
      RETURN
C------------- LAST LINE OF DRLCAL FOLLOWS ----------------------------
      END
*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
      use omp_lib
C***BEGIN PROLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DCOPY
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
      use omp_lib
C***BEGIN PROLOGUE  DCOPY
C***PURPOSE  Copy a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy double precision DX to double precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCOPY
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END
*DECK DHELS
      SUBROUTINE DHELS (A, LDA, N, Q, B)
C***BEGIN PROLOGUE  DHELS
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SHELS-S, DHELS-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C        This routine is extracted from the LINPACK routine SGESL with
C        changes due to the fact that A is an upper Hessenberg matrix.
C
C        DHELS solves the least squares problem:
C
C                   MIN(B-A*X,B-A*X)
C
C        using the factors computed by DHEQR.
C
C *Usage:
C      INTEGER LDA, N
C      DOUBLE PRECISION A(LDA,N), Q(2*N), B(N+1)
C
C      CALL DHELS(A, LDA, N, Q, B)
C
C *Arguments:
C A       :IN       Double Precision A(LDA,N)
C          The output from DHEQR which contains the upper
C          triangular factor R in the QR decomposition of A.
C LDA     :IN       Integer
C          The leading dimension of the array A.
C N       :IN       Integer
C          A is originally an (N+1) by N matrix.
C Q       :IN       Double Precision Q(2*N)
C          The coefficients of the N Givens rotations
C          used in the QR factorization of A.
C B       :INOUT    Double Precision B(N+1)
C          On input, B is the right hand side vector.
C          On output, B is the solution vector X.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DAXPY
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DHELS
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      INTEGER LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
C     .. Local Scalars ..
      DOUBLE PRECISION C, S, T, T1, T2
      INTEGER IQ, K, KB, KP1
C     .. External Subroutines ..
      EXTERNAL DAXPY
C***FIRST EXECUTABLE STATEMENT  DHELS
C
C         Minimize(B-A*X,B-A*X).  First form Q*B.
C
      DO 20 K = 1, N
         KP1 = K + 1
         IQ = 2*(K-1) + 1
         C = Q(IQ)
         S = Q(IQ+1)
         T1 = B(K)
         T2 = B(KP1)
         B(K) = C*T1 - S*T2
         B(KP1) = S*T1 + C*T2
 20   CONTINUE
C
C         Now solve  R*X = Q*B.
C
      DO 40 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL DAXPY(K-1, T, A(1,K), 1, B(1), 1)
 40   CONTINUE
      RETURN
C------------- LAST LINE OF DHELS FOLLOWS ----------------------------
      END
*DECK DHEQR
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
C***BEGIN PROLOGUE  DHEQR
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SHEQR-S, DHEQR-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C        This   routine  performs  a QR   decomposition  of an  upper
C        Hessenberg matrix A using Givens  rotations.  There  are two
C        options  available: 1)  Performing  a fresh decomposition 2)
C        updating the QR factors by adding a row and  a column to the
C        matrix A.
C
C *Usage:
C      INTEGER LDA, N, INFO, IJOB
C      DOUBLE PRECISION A(LDA,N), Q(2*N)
C
C      CALL DHEQR(A, LDA, N, Q, INFO, IJOB)
C
C *Arguments:
C A      :INOUT    Double Precision A(LDA,N)
C         On input, the matrix to be decomposed.
C         On output, the upper triangular matrix R.
C         The factorization can be written Q*A = R, where
C         Q is a product of Givens rotations and R is upper
C         triangular.
C LDA    :IN       Integer
C         The leading dimension of the array A.
C N      :IN       Integer
C         A is an (N+1) by N Hessenberg matrix.
C Q      :OUT      Double Precision Q(2*N)
C         The factors c and s of each Givens rotation used
C         in decomposing A.
C INFO   :OUT      Integer
C         = 0  normal value.
C         = K  if  A(K,K) .eq. 0.0 .  This is not an error
C           condition for this subroutine, but it does
C           indicate that DHELS will divide by zero
C           if called.
C IJOB   :IN       Integer
C         = 1     means that a fresh decomposition of the
C                 matrix A is desired.
C         .ge. 2  means that the current decomposition of A
C                 will be updated by the addition of a row
C                 and a column.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DHEQR
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      INTEGER IJOB, INFO, LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), Q(*)
C     .. Local Scalars ..
      DOUBLE PRECISION C, S, T, T1, T2
      INTEGER I, IQ, J, K, KM1, KP1, NM1
C     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
C***FIRST EXECUTABLE STATEMENT  DHEQR
      IF (IJOB .GT. 1) GO TO 70
C   -------------------------------------------------------------------
C         A new factorization is desired.
C   -------------------------------------------------------------------
C         QR decomposition without pivoting.
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
C           Compute K-th column of R.
C           First, multiply the K-th column of A by the previous
C           K-1 Givens rotations.
C
         IF (KM1 .LT. 1) GO TO 20
         DO 10 J = 1, KM1
            I = 2*(J-1) + 1
            T1 = A(J,K)
            T2 = A(J+1,K)
            C = Q(I)
            S = Q(I+1)
            A(J,K) = C*T1 - S*T2
            A(J+1,K) = S*T1 + C*T2
 10      CONTINUE
C
C         Compute Givens components C and S.
C
 20      CONTINUE
         IQ = 2*KM1 + 1
         T1 = A(K,K)
         T2 = A(KP1,K)
         IF( T2.EQ.0.0D0 ) THEN
            C = 1
            S = 0
         ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
            T = T1/T2
            S = -1.0D0/SQRT(1.0D0+T*T)
            C = -S*T
         ELSE
            T = T2/T1
            C = 1.0D0/SQRT(1.0D0+T*T)
            S = -C*T
         ENDIF
         Q(IQ) = C
         Q(IQ+1) = S
         A(K,K) = C*T1 - S*T2
         IF( A(K,K).EQ.0.0D0 ) INFO = K
 60   CONTINUE
      RETURN
C   -------------------------------------------------------------------
C         The old factorization of a will be updated.  A row and a
C         column has been added to the matrix A.  N by N-1 is now
C         the old size of the matrix.
C   -------------------------------------------------------------------
 70   CONTINUE
      NM1 = N - 1
C   -------------------------------------------------------------------
C         Multiply the new column by the N previous Givens rotations.
C   -------------------------------------------------------------------
      DO 100 K = 1,NM1
         I = 2*(K-1) + 1
         T1 = A(K,N)
         T2 = A(K+1,N)
         C = Q(I)
         S = Q(I+1)
         A(K,N) = C*T1 - S*T2
         A(K+1,N) = S*T1 + C*T2
 100  CONTINUE
C   -------------------------------------------------------------------
C         Complete update of decomposition by forming last Givens
C         rotation, and multiplying it times the column
C         vector(A(N,N),A(NP1,N)).
C   -------------------------------------------------------------------
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF ( T2.EQ.0.0D0 ) THEN
         C = 1
         S = 0
      ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
         T = T1/T2
         S = -1.0D0/SQRT(1.0D0+T*T)
         C = -S*T
      ELSE
         T = T2/T1
         C = 1.0D0/SQRT(1.0D0+T*T)
         S = -C*T
      ENDIF
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C------------- LAST LINE OF DHEQR FOLLOWS ----------------------------
      END
*DECK ISDGMR
      INTEGER FUNCTION ISDGMR (N, B, X, XL, NELT, IA, JA, A, ISYM,
     +   MSOLVE, NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
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
C     .. Subroutine Arguments ..
      EXTERNAL MSOLVE
C     .. Arrays in Common ..
      DOUBLE PRECISION SOLN(1)
C     .. Local Scalars ..
      DOUBLE PRECISION DXNRM, FUZZ, RAT, RATMAX, SOLNRM, TEM
      INTEGER I, IELMAX
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2
      EXTERNAL D1MACH, DNRM2
C     .. External Subroutines ..
      EXTERNAL DCOPY, DRLCAL, DSCAL, DXLCAL
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON /DSLBLK/ SOLN
C     .. Save statement ..
      SAVE SOLNRM
C***FIRST EXECUTABLE STATEMENT  ISDGMR
      ISDGMR = 0
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
            IF( RATMAX.LE.TOL ) ISDGMR = 1
            IF( IUNIT.GT.0 ) WRITE(IUNIT,1020) ITER, IELMAX, RATMAX
            RETURN
         ENDIF
      ENDIF
      IF ( ITOL.EQ.11 ) THEN
C
C       Use DXLCAL to calculate the approximate solution XL.
C
         IF ( (LGMR.NE.0) .AND. (ITER.GT.0) ) THEN
            CALL DXLCAL(N, LGMR, X, XL, XL, HES, MAXLP1, Q, V, R0NRM,
     $           DZ, SX, JSCAL, JPRE, MSOLVE, NMSL, RWORK, IWORK,
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
      IF ( ERR.LE.TOL ) ISDGMR = 1
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
*DECK DORTH
      SUBROUTINE DORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
C***BEGIN PROLOGUE  DORTH
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SORTH-S, DORTH-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C        This routine  orthogonalizes  the  vector  VNEW  against the
C        previous KMP  vectors in the   V array.  It uses  a modified
C        Gram-Schmidt   orthogonalization procedure with  conditional
C        reorthogonalization.
C
C *Usage:
C      INTEGER N, LL, LDHES, KMP
C      DOUBLE PRECISION VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW
C
C      CALL DORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
C
C *Arguments:
C VNEW   :INOUT    Double Precision VNEW(N)
C         On input, the vector of length N containing a scaled
C         product of the Jacobian and the vector V(*,LL).
C         On output, the new vector orthogonal to V(*,i0) to V(*,LL),
C         where i0 = max(1, LL-KMP+1).
C V      :IN       Double Precision V(N,LL)
C         The N x LL array containing the previous LL
C         orthogonal vectors V(*,1) to V(*,LL).
C HES    :INOUT    Double Precision HES(LDHES,LL)
C         On input, an LL x LL upper Hessenberg matrix containing,
C         in HES(I,K), K.lt.LL, the scaled inner products of
C         A*V(*,K) and V(*,i).
C         On return, column LL of HES is filled in with
C         the scaled inner products of A*V(*,LL) and V(*,i).
C N      :IN       Integer
C         The order of the matrix A, and the length of VNEW.
C LL     :IN       Integer
C         The current order of the matrix HES.
C LDHES  :IN       Integer
C         The leading dimension of the HES array.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to (KMP .le. MAXL).
C SNORMW :OUT      DOUBLE PRECISION
C         Scalar containing the l-2 norm of VNEW.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DAXPY, DDOT, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DORTH
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      DOUBLE PRECISION SNORMW
      INTEGER KMP, LDHES, LL, N
C     .. Array Arguments ..
      DOUBLE PRECISION HES(LDHES,*), V(N,*), VNEW(*)
C     .. Local Scalars ..
      DOUBLE PRECISION ARG, SUMDSQ, TEM, VNRM
      INTEGER I, I0
C     .. External Functions ..
      DOUBLE PRECISION DDOT, DNRM2
      EXTERNAL DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL DAXPY
C     .. Intrinsic Functions ..
      INTRINSIC MAX, SQRT
C***FIRST EXECUTABLE STATEMENT  DORTH
C
C         Get norm of unaltered VNEW for later use.
C
      VNRM = DNRM2(N, VNEW, 1)
C   -------------------------------------------------------------------
C         Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL).
C         Scaled inner products give new column of HES.
C         Projections of earlier vectors are subtracted from VNEW.
C   -------------------------------------------------------------------
      I0 = MAX(1,LL-KMP+1)
      DO 10 I = I0,LL
         HES(I,LL) = DDOT(N, V(1,I), 1, VNEW, 1)
         TEM = -HES(I,LL)
         CALL DAXPY(N, TEM, V(1,I), 1, VNEW, 1)
 10   CONTINUE
C   -------------------------------------------------------------------
C         Compute SNORMW = norm of VNEW.  If VNEW is small compared
C         to its input value (in norm), then reorthogonalize VNEW to
C         V(*,1) through V(*,LL).  Correct if relative correction
C         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using
C         the dot products involved.
C   -------------------------------------------------------------------
      SNORMW = DNRM2(N, VNEW, 1)
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0
      DO 30 I = I0,LL
         TEM = -DDOT(N, V(1,I), 1, VNEW, 1)
         IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
         HES(I,LL) = HES(I,LL) - TEM
         CALL DAXPY(N, TEM, V(1,I), 1, VNEW, 1)
         SUMDSQ = SUMDSQ + TEM**2
 30   CONTINUE
      IF (SUMDSQ .EQ. 0.0D0) RETURN
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
C
      RETURN
C------------- LAST LINE OF DORTH FOLLOWS ----------------------------
      END
*DECK DXLCAL
      SUBROUTINE DXLCAL (N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
     +   WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR, NELT, IA, JA, A,
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
C     .. Subroutine Arguments ..
      EXTERNAL MSOLVE
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
*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
      use omp_lib
C***BEGIN PROLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END
