!
!    SLATEC: public domain
!
*DECK DDEABM
      SUBROUTINE DDEABM (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,
     +   RWORK, LRW, IWORK, LIW, RPAR, IPAR)
C***BEGIN PROLOGUE  DDEABM
C***PURPOSE  Solve an initial value problem in ordinary differential
C            equations using an Adams-Bashforth method.
C***LIBRARY   SLATEC (DEPAC)
C***CATEGORY  I1A1B
C***TYPE      DOUBLE PRECISION (DEABM-S, DDEABM-D)
C***KEYWORDS  ADAMS-BASHFORTH METHOD, DEPAC, INITIAL VALUE PROBLEMS,
C             ODE, ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
C***AUTHOR  Shampine, L. F., (SNLA)
C           Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   This is the Adams code in the package of differential equation
C   solvers DEPAC, consisting of the codes DDERKF, DDEABM, and DDEBDF.
C   Design of the package was by L. F. Shampine and H. A. Watts.
C   It is documented in
C        SAND79-2374 , DEPAC - Design of a User Oriented Package of ODE
C                              Solvers.
C   DDEABM is a driver for a modification of the code ODE written by
C             L. F. Shampine and M. K. Gordon
C             Sandia Laboratories
C             Albuquerque, New Mexico 87185
C
C **********************************************************************
C * ABSTRACT *
C ************
C
C   Subroutine DDEABM uses the Adams-Bashforth-Moulton
C   Predictor-Corrector formulas of orders one through twelve to
C   integrate a system of NEQ first order ordinary differential
C   equations of the form
C                         DU/DX = DF(X,U)
C   when the vector Y(*) of initial values for U(*) at X=T is given.
C   The subroutine integrates from T to TOUT. It is easy to continue the
C   integration to get results at additional TOUT.  This is the interval
C   mode of operation.  It is also easy for the routine to return with
C   the solution at each intermediate step on the way to TOUT.  This is
C   the intermediate-output mode of operation.
C
C   DDEABM uses subprograms DDES, DSTEPS, DINTP, DHSTRT, DHVNRM,
C   D1MACH, and the error handling routine XERMSG.  The only machine
C   dependent parameters to be assigned appear in D1MACH.
C
C **********************************************************************
C * Description of The Arguments To DDEABM (An Overview) *
C **********************************************************************
C
C   The Parameters are
C
C      DF -- This is the name of a subroutine which you provide to
C             define the differential equations.
C
C      NEQ -- This is the number of (first order) differential
C             equations to be integrated.
C
C      T -- This is a DOUBLE PRECISION value of the independent
C           variable.
C
C      Y(*) -- This DOUBLE PRECISION array contains the solution
C              components at T.
C
C      TOUT -- This is a DOUBLE PRECISION point at which a solution is
C              desired.
C
C      INFO(*) -- The basic task of the code is to integrate the
C             differential equations from T to TOUT and return an
C             answer at TOUT.  INFO(*) is an INTEGER array which is used
C             to communicate exactly how you want this task to be
C             carried out.
C
C      RTOL, ATOL -- These DOUBLE PRECISION quantities represent
C                    relative and absolute error tolerances which you
C                    provide to indicate how accurately you wish the
C                    solution to be computed.  You may choose them to be
C                    both scalars or else both vectors.
C
C      IDID -- This scalar quantity is an indicator reporting what
C             the code did.  You must monitor this INTEGER variable to
C             decide what action to take next.
C
C      RWORK(*), LRW -- RWORK(*) is a DOUBLE PRECISION work array of
C             length LRW which provides the code with needed storage
C             space.
C
C      IWORK(*), LIW -- IWORK(*) is an INTEGER work array of length LIW
C             which provides the code with needed storage space and an
C             across call flag.
C
C      RPAR, IPAR -- These are DOUBLE PRECISION and INTEGER parameter
C             arrays which you can use for communication between your
C             calling program and the DF subroutine.
C
C  Quantities which are used as input items are
C             NEQ, T, Y(*), TOUT, INFO(*),
C             RTOL, ATOL, RWORK(1), LRW and LIW.
C
C  Quantities which may be altered by the code are
C             T, Y(*), INFO(1), RTOL, ATOL,
C             IDID, RWORK(*) and IWORK(*).
C
C **********************************************************************
C * INPUT -- What To Do On The First Call To DDEABM *
C **********************************************************************
C
C   The first call of the code is defined to be the start of each new
C   problem.  Read through the descriptions of all the following items,
C   provide sufficient storage space for designated arrays, set
C   appropriate variables for the initialization of the problem, and
C   give information about how you want the problem to be solved.
C
C
C      DF -- Provide a subroutine of the form
C                               DF(X,U,UPRIME,RPAR,IPAR)
C             to define the system of first order differential equations
C             which is to be solved.  For the given values of X and the
C             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
C             evaluate the NEQ components of the system of differential
C             equations  DU/DX=DF(X,U)  and store the derivatives in the
C             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
C             equations I=1,...,NEQ.
C
C             Subroutine DF must NOT alter X or U(*).  You must declare
C             the name df in an external statement in your program that
C             calls DDEABM.  You must dimension U and UPRIME in DF.
C
C             RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
C             arrays which you can use for communication between your
C             calling program and subroutine DF. They are not used or
C             altered by DDEABM.  If you do not need RPAR or IPAR,
C             ignore these parameters by treating them as dummy
C             arguments. If you do choose to use them, dimension them in
C             your calling program and in DF as arrays of appropriate
C             length.
C
C      NEQ -- Set it to the number of differential equations.
C             (NEQ .GE. 1)
C
C      T -- Set it to the initial point of the integration.
C             You must use a program variable for T because the code
C             changes its value.
C
C      Y(*) -- Set this vector to the initial values of the NEQ solution
C             components at the initial point.  You must dimension Y at
C             least NEQ in your calling program.
C
C      TOUT -- Set it to the first point at which a solution
C             is desired.  You can take TOUT = T, in which case the code
C             will evaluate the derivative of the solution at T and
C             return. Integration either forward in T  (TOUT .GT. T)  or
C             backward in T  (TOUT .LT. T)  is permitted.
C
C             The code advances the solution from T to TOUT using
C             step sizes which are automatically selected so as to
C             achieve the desired accuracy.  If you wish, the code will
C             return with the solution and its derivative following
C             each intermediate step (intermediate-output mode) so that
C             you can monitor them, but you still must provide TOUT in
C             accord with the basic aim of the code.
C
C             The first step taken by the code is a critical one
C             because it must reflect how fast the solution changes near
C             the initial point.  The code automatically selects an
C             initial step size which is practically always suitable for
C             the problem. By using the fact that the code will not step
C             past TOUT in the first step, you could, if necessary,
C             restrict the length of the initial step size.
C
C             For some problems it may not be permissible to integrate
C             past a point TSTOP because a discontinuity occurs there
C             or the solution or its derivative is not defined beyond
C             TSTOP.  When you have declared a TSTOP point (see INFO(4)
C             and RWORK(1)), you have told the code not to integrate
C             past TSTOP.  In this case any TOUT beyond TSTOP is invalid
C             input.
C
C      INFO(*) -- Use the INFO array to give the code more details about
C             how you want your problem solved.  This array should be
C             dimensioned of length 15 to accommodate other members of
C             DEPAC or possible future extensions, though DDEABM uses
C             only the first four entries.  You must respond to all of
C             the following items which are arranged as questions.  The
C             simplest use of the code corresponds to answering all
C             questions as YES ,i.e. setting ALL entries of INFO to 0.
C
C        INFO(1) -- This parameter enables the code to initialize
C               itself.  You must set it to indicate the start of every
C               new problem.
C
C            **** Is this the first call for this problem ...
C                  YES -- set INFO(1) = 0
C                   NO -- not applicable here.
C                         See below for continuation calls.  ****
C
C        INFO(2) -- How much accuracy you want of your solution
C               is specified by the error tolerances RTOL and ATOL.
C               The simplest use is to take them both to be scalars.
C               To obtain more flexibility, they can both be vectors.
C               The code must be told your choice.
C
C            **** Are both error tolerances RTOL, ATOL scalars ...
C                  YES -- set INFO(2) = 0
C                         and input scalars for both RTOL and ATOL
C                   NO -- set INFO(2) = 1
C                         and input arrays for both RTOL and ATOL ****
C
C        INFO(3) -- The code integrates from T in the direction
C               of TOUT by steps.  If you wish, it will return the
C               computed solution and derivative at the next
C               intermediate step (the intermediate-output mode) or
C               TOUT, whichever comes first.  This is a good way to
C               proceed if you want to see the behavior of the solution.
C               If you must have solutions at a great many specific
C               TOUT points, this code will compute them efficiently.
C
C            **** Do you want the solution only at
C                 TOUT (and not at the next intermediate step) ...
C                  YES -- set INFO(3) = 0
C                   NO -- set INFO(3) = 1 ****
C
C        INFO(4) -- To handle solutions at a great many specific
C               values TOUT efficiently, this code may integrate past
C               TOUT and interpolate to obtain the result at TOUT.
C               Sometimes it is not possible to integrate beyond some
C               point TSTOP because the equation changes there or it is
C               not defined past TSTOP.  Then you must tell the code
C               not to go past.
C
C            **** Can the integration be carried out without any
C                 Restrictions on the independent variable T ...
C                  YES -- set INFO(4)=0
C                   NO -- set INFO(4)=1
C                         and define the stopping point TSTOP by
C                         setting RWORK(1)=TSTOP ****
C
C      RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
C             error tolerances to tell the code how accurately you want
C             the solution to be computed.  They must be defined as
C             program variables because the code may change them.  You
C             have two choices --
C                  Both RTOL and ATOL are scalars. (INFO(2)=0)
C                  Both RTOL and ATOL are vectors. (INFO(2)=1)
C             In either case all components must be non-negative.
C
C             The tolerances are used by the code in a local error test
C             at each step which requires roughly that
C                     ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL
C             for each vector component.
C             (More specifically, a Euclidean norm is used to measure
C             the size of vectors, and the error test uses the magnitude
C             of the solution at the beginning of the step.)
C
C             The true (global) error is the difference between the true
C             solution of the initial value problem and the computed
C             approximation.  Practically all present day codes,
C             including this one, control the local error at each step
C             and do not even attempt to control the global error
C             directly.  Roughly speaking, they produce a solution Y(T)
C             which satisfies the differential equations with a
C             residual R(T),    DY(T)/DT = DF(T,Y(T)) + R(T)   ,
C             and, almost always, R(T) is bounded by the error
C             tolerances.  Usually, but not always, the true accuracy of
C             the computed Y is comparable to the error tolerances. This
C             code will usually, but not always, deliver a more accurate
C             solution if you reduce the tolerances and integrate again.
C             By comparing two such solutions you can get a fairly
C             reliable idea of the true error in the solution at the
C             bigger tolerances.
C
C             Setting ATOL=0.D0 results in a pure relative error test on
C             that component. Setting RTOL=0. results in a pure absolute
C             error test on that component.  A mixed test with non-zero
C             RTOL and ATOL corresponds roughly to a relative error
C             test when the solution component is much bigger than ATOL
C             and to an absolute error test when the solution component
C             is smaller than the threshold ATOL.
C
C             Proper selection of the absolute error control parameters
C             ATOL  requires you to have some idea of the scale of the
C             solution components.  To acquire this information may mean
C             that you will have to solve the problem more than once. In
C             the absence of scale information, you should ask for some
C             relative accuracy in all the components (by setting  RTOL
C             values non-zero) and perhaps impose extremely small
C             absolute error tolerances to protect against the danger of
C             a solution component becoming zero.
C
C             The code will not attempt to compute a solution at an
C             accuracy unreasonable for the machine being used.  It will
C             advise you if you ask for too much accuracy and inform
C             you as to the maximum accuracy it believes possible.
C
C      RWORK(*) -- Dimension this DOUBLE PRECISION work array of length
C             LRW in your calling program.
C
C      RWORK(1) -- If you have set INFO(4)=0, you can ignore this
C             optional input parameter.  Otherwise you must define a
C             stopping point TSTOP by setting   RWORK(1) = TSTOP.
C             (for some problems it may not be permissible to integrate
C             past a point TSTOP because a discontinuity occurs there
C             or the solution or its derivative is not defined beyond
C             TSTOP.)
C
C      LRW -- Set it to the declared length of the RWORK array.
C             You must have  LRW .GE. 130+21*NEQ
C
C      IWORK(*) -- Dimension this INTEGER work array of length LIW in
C             your calling program.
C
C      LIW -- Set it to the declared length of the IWORK array.
C             You must have  LIW .GE. 51
C
C      RPAR, IPAR -- These are parameter arrays, of DOUBLE PRECISION and
C             INTEGER type, respectively.  You can use them for
C             communication between your program that calls DDEABM and
C             the  DF subroutine.  They are not used or altered by
C             DDEABM.  If you do not need RPAR or IPAR, ignore these
C             parameters by treating them as dummy arguments.  If you do
C             choose to use them, dimension them in your calling program
C             and in DF as arrays of appropriate length.
C
C **********************************************************************
C * OUTPUT -- After Any Return From DDEABM *
C **********************************************************************
C
C   The principal aim of the code is to return a computed solution at
C   TOUT, although it is also possible to obtain intermediate results
C   along the way.  To find out whether the code achieved its goal
C   or if the integration process was interrupted before the task was
C   completed, you must check the IDID parameter.
C
C
C      T -- The solution was successfully advanced to the
C             output value of T.
C
C      Y(*) -- Contains the computed solution approximation at T.
C             You may also be interested in the approximate derivative
C             of the solution at T.  It is contained in
C             RWORK(21),...,RWORK(20+NEQ).
C
C      IDID -- Reports what the code did
C
C                         *** Task Completed ***
C                   Reported by positive values of IDID
C
C             IDID = 1 -- A step was successfully taken in the
C                       intermediate-output mode.  The code has not
C                       yet reached TOUT.
C
C             IDID = 2 -- The integration to TOUT was successfully
C                       completed (T=TOUT) by stepping exactly to TOUT.
C
C             IDID = 3 -- The integration to TOUT was successfully
C                       completed (T=TOUT) by stepping past TOUT.
C                       Y(*) is obtained by interpolation.
C
C                         *** Task Interrupted ***
C                   Reported by negative values of IDID
C
C             IDID = -1 -- A large amount of work has been expended.
C                       (500 steps attempted)
C
C             IDID = -2 -- The error tolerances are too stringent.
C
C             IDID = -3 -- The local error test cannot be satisfied
C                       because you specified a zero component in ATOL
C                       and the corresponding computed solution
C                       component is zero.  Thus, a pure relative error
C                       test is impossible for this component.
C
C             IDID = -4 -- The problem appears to be stiff.
C
C             IDID = -5,-6,-7,..,-32  -- Not applicable for this code
C                       but used by other members of DEPAC or possible
C                       future extensions.
C
C                         *** Task Terminated ***
C                   Reported by the value of IDID=-33
C
C             IDID = -33 -- The code has encountered trouble from which
C                       it cannot recover.  A message is printed
C                       explaining the trouble and control is returned
C                       to the calling program. For example, this occurs
C                       when invalid input is detected.
C
C      RTOL, ATOL -- These quantities remain unchanged except when
C             IDID = -2.  In this case, the error tolerances have been
C             increased by the code to values which are estimated to be
C             appropriate for continuing the integration.  However, the
C             reported solution at T was obtained using the input values
C             of RTOL and ATOL.
C
C      RWORK, IWORK -- Contain information which is usually of no
C             interest to the user but necessary for subsequent calls.
C             However, you may find use for
C
C             RWORK(11)--which contains the step size H to be
C                        attempted on the next step.
C
C             RWORK(12)--if the tolerances have been increased by the
C                        code (IDID = -2) , they were multiplied by the
C                        value in RWORK(12).
C
C             RWORK(13)--Which contains the current value of the
C                        independent variable, i.e. the farthest point
C                        integration has reached. This will be different
C                        from T only when interpolation has been
C                        performed (IDID=3).
C
C             RWORK(20+I)--Which contains the approximate derivative
C                        of the solution component Y(I).  In DDEABM, it
C                        is obtained by calling subroutine DF to
C                        evaluate the differential equation using T and
C                        Y(*) when IDID=1 or 2, and by interpolation
C                        when IDID=3.
C
C **********************************************************************
C * INPUT -- What To Do To Continue The Integration *
C *             (calls after the first)             *
C **********************************************************************
C
C        This code is organized so that subsequent calls to continue the
C        integration involve little (if any) additional effort on your
C        part. You must monitor the IDID parameter in order to determine
C        what to do next.
C
C        Recalling that the principal task of the code is to integrate
C        from T to TOUT (the interval mode), usually all you will need
C        to do is specify a new TOUT upon reaching the current TOUT.
C
C        Do not alter any quantity not specifically permitted below,
C        in particular do not alter NEQ, T, Y(*), RWORK(*), IWORK(*) or
C        the differential equation in subroutine DF. Any such alteration
C        constitutes a new problem and must be treated as such, i.e.
C        you must start afresh.
C
C        You cannot change from vector to scalar error control or vice
C        versa (INFO(2)) but you can change the size of the entries of
C        RTOL, ATOL.  Increasing a tolerance makes the equation easier
C        to integrate.  Decreasing a tolerance will make the equation
C        harder to integrate and should generally be avoided.
C
C        You can switch from the intermediate-output mode to the
C        interval mode (INFO(3)) or vice versa at any time.
C
C        If it has been necessary to prevent the integration from going
C        past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
C        code will not integrate to any TOUT beyond the currently
C        specified TSTOP.  Once TSTOP has been reached you must change
C        the value of TSTOP or set INFO(4)=0.  You may change INFO(4)
C        or TSTOP at any time but you must supply the value of TSTOP in
C        RWORK(1) whenever you set INFO(4)=1.
C
C        The parameter INFO(1) is used by the code to indicate the
C        beginning of a new problem and to indicate whether integration
C        is to be continued.  You must input the value  INFO(1) = 0
C        when starting a new problem.  You must input the value
C        INFO(1) = 1  if you wish to continue after an interrupted task.
C        Do not set  INFO(1) = 0  on a continuation call unless you
C        want the code to restart at the current T.
C
C                         *** Following A Completed Task ***
C         If
C             IDID = 1, call the code again to continue the integration
C                     another step in the direction of TOUT.
C
C             IDID = 2 or 3, define a new TOUT and call the code again.
C                     TOUT must be different from T. You cannot change
C                     the direction of integration without restarting.
C
C                         *** Following An Interrupted Task ***
C                     To show the code that you realize the task was
C                     interrupted and that you want to continue, you
C                     must take appropriate action and reset INFO(1) = 1
C         If
C             IDID = -1, the code has attempted 500 steps.
C                     If you want to continue, set INFO(1) = 1 and
C                     call the code again. An additional 500 steps
C                     will be allowed.
C
C             IDID = -2, the error tolerances RTOL, ATOL have been
C                     increased to values the code estimates appropriate
C                     for continuing.  You may want to change them
C                     yourself.  If you are sure you want to continue
C                     with relaxed error tolerances, set INFO(1)=1 and
C                     call the code again.
C
C             IDID = -3, a solution component is zero and you set the
C                     corresponding component of ATOL to zero.  If you
C                     are sure you want to continue, you must first
C                     alter the error criterion to use positive values
C                     for those components of ATOL corresponding to zero
C                     solution components, then set INFO(1)=1 and call
C                     the code again.
C
C             IDID = -4, the problem appears to be stiff.  It is very
C                     inefficient to solve such problems with DDEABM.
C                     The code DDEBDF in DEPAC handles this task
C                     efficiently.  If you are absolutely sure you want
C                     to continue with DDEABM, set INFO(1)=1 and call
C                     the code again.
C
C             IDID = -5,-6,-7,..,-32  --- cannot occur with this code
C                     but used by other members of DEPAC or possible
C                     future extensions.
C
C                         *** Following A Terminated Task ***
C         If
C             IDID = -33, you cannot continue the solution of this
C                     problem.  An attempt to do so will result in your
C                     run being terminated.
C
C **********************************************************************
C *Long Description:
C
C **********************************************************************
C *             DEPAC Package Overview           *
C **********************************************************************
C
C ....   You have a choice of three differential equation solvers from
C ....   DEPAC. The following brief descriptions are meant to aid you in
C ....   choosing the most appropriate code for your problem.
C
C ....   DDERKF is a fifth order Runge-Kutta code. It is the simplest of
C ....   the three choices, both algorithmically and in the use of the
C ....   code. DDERKF is primarily designed to solve non-stiff and
C ....   mildly stiff differential equations when derivative evaluations
C ....   are not expensive. It should generally not be used to get high
C ....   accuracy results nor answers at a great many specific points.
C ....   Because DDERKF has very low overhead costs, it will usually
C ....   result in the least expensive integration when solving
C ....   problems requiring a modest amount of accuracy and having
C ....   equations that are not costly to evaluate. DDERKF attempts to
C ....   discover when it is not suitable for the task posed.
C
C ....   DDEABM is a variable order (one through twelve) Adams code.
C ....   Its complexity lies somewhere between that of DDERKF and
C ....   DDEBDF.  DDEABM is primarily designed to solve non-stiff and
C ....   mildly stiff differential equations when derivative evaluations
C ....   are expensive, high accuracy results are needed or answers at
C ....   many specific points are required. DDEABM attempts to discover
C ....   when it is not suitable for the task posed.
C
C ....   DDEBDF is a variable order (one through five) backward
C ....   differentiation formula code. it is the most complicated of
C ....   the three choices. DDEBDF is primarily designed to solve stiff
C ....   differential equations at crude to moderate tolerances.
C ....   If the problem is very stiff at all, DDERKF and DDEABM will be
C ....   quite inefficient compared to DDEBDF. However, DDEBDF will be
C ....   inefficient compared to DDERKF and DDEABM on non-stiff problems
C ....   because it uses much more storage, has a much larger overhead,
C ....   and the low order formulas will not give high accuracies
C ....   efficiently.
C
C ....   The concept of stiffness cannot be described in a few words.
C ....   If you do not know the problem to be stiff, try either DDERKF
C ....   or DDEABM. Both of these codes will inform you of stiffness
C ....   when the cost of solving such problems becomes important.
C
C *********************************************************************
C
C***REFERENCES  L. F. Shampine and H. A. Watts, DEPAC - design of a user
C                 oriented package of ODE solvers, Report SAND79-2374,
C                 Sandia Laboratories, 1979.
C***ROUTINES CALLED  DDES, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891024  Changed references from DVNORM to DHVNRM.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDEABM
C
      INTEGER IALPHA, IBETA, IDELSN, IDID, IFOURU, IG, IHOLD,
     1      INFO, IP, IPAR, IPHI, IPSI, ISIG, ITOLD, ITSTAR, ITWOU,
     2      IV, IW, IWORK, IWT, IYP, IYPOUT, IYY, LIW, LRW, NEQ
      DOUBLE PRECISION ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
      LOGICAL START,PHASE1,NORND,STIFF,INTOUT
C
      DIMENSION Y(*),INFO(15),RTOL(*),ATOL(*),RWORK(*),IWORK(*),
     1          RPAR(*),IPAR(*)
C
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3
C
      EXTERNAL DF
C
C     CHECK FOR AN APPARENT INFINITE LOOP
C
C***FIRST EXECUTABLE STATEMENT  DDEABM
      IF ( INFO(1) .EQ. 0 ) IWORK(LIW) = 0
      IF (IWORK(LIW) .GE. 5) THEN
         IF (T .EQ. RWORK(21 + NEQ)) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDEABM',
     *         'AN APPARENT INFINITE LOOP HAS BEEN DETECTED.$$' //
     *         'YOU HAVE MADE REPEATED CALLS AT T = ' // XERN3 //
     *         ' AND THE INTEGRATION HAS NOT ADVANCED.  CHECK THE ' //
     *         'WAY YOU HAVE SET PARAMETERS FOR THE CALL TO THE ' //
     *         'CODE, PARTICULARLY INFO(1).', 13, 2)
            RETURN
         ENDIF
      ENDIF
C
C     CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION
C
      IDID=0
      IF (LRW .LT. 130+21*NEQ) THEN
         WRITE (XERN1, '(I8)') LRW
         CALL XERMSG ('SLATEC', 'DDEABM', 'THE LENGTH OF THE RWORK ' //
     *      'ARRAY MUST BE AT LEAST 130 + 21*NEQ.$$' //
     *      'YOU HAVE CALLED THE CODE WITH LRW = ' // XERN1, 1, 1)
         IDID=-33
      ENDIF
C
      IF (LIW .LT. 51) THEN
         WRITE (XERN1, '(I8)') LIW
         CALL XERMSG ('SLATEC', 'DDEABM', 'THE LENGTH OF THE IWORK ' //
     *      'ARRAY MUST BE AT LEAST 51.$$YOU HAVE CALLED THE CODE ' //
     *      'WITH LIW = ' // XERN1, 2, 1)
         IDID=-33
      ENDIF
C
C     COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK ARRAY
C
      IYPOUT = 21
      ITSTAR = NEQ + 21
      IYP = 1 + ITSTAR
      IYY = NEQ + IYP
      IWT = NEQ + IYY
      IP = NEQ + IWT
      IPHI = NEQ + IP
      IALPHA = (NEQ*16) + IPHI
      IBETA = 12 + IALPHA
      IPSI = 12 + IBETA
      IV = 12 + IPSI
      IW = 12 + IV
      ISIG = 12 + IW
      IG = 13 + ISIG
      IGI = 13 + IG
      IXOLD = 11 + IGI
      IHOLD = 1 + IXOLD
      ITOLD = 1 + IHOLD
      IDELSN = 1 + ITOLD
      ITWOU = 1 + IDELSN
      IFOURU = 1 + ITWOU
C
      RWORK(ITSTAR) = T
      IF (INFO(1) .EQ. 0) GO TO 50
      START = IWORK(21) .NE. (-1)
      PHASE1 = IWORK(22) .NE. (-1)
      NORND = IWORK(23) .NE. (-1)
      STIFF = IWORK(24) .NE. (-1)
      INTOUT = IWORK(25) .NE. (-1)
C
 50   CALL DDES(DF,NEQ,T,Y,TOUT,INFO,RTOL,ATOL,IDID,RWORK(IYPOUT),
     1         RWORK(IYP),RWORK(IYY),RWORK(IWT),RWORK(IP),RWORK(IPHI),
     2         RWORK(IALPHA),RWORK(IBETA),RWORK(IPSI),RWORK(IV),
     3         RWORK(IW),RWORK(ISIG),RWORK(IG),RWORK(IGI),RWORK(11),
     4         RWORK(12),RWORK(13),RWORK(IXOLD),RWORK(IHOLD),
     5         RWORK(ITOLD),RWORK(IDELSN),RWORK(1),RWORK(ITWOU),
     5         RWORK(IFOURU),START,PHASE1,NORND,STIFF,INTOUT,IWORK(26),
     6         IWORK(27),IWORK(28),IWORK(29),IWORK(30),IWORK(31),
     7         IWORK(32),IWORK(33),IWORK(34),IWORK(35),IWORK(45),
     8         RPAR,IPAR)
C
      IWORK(21) = -1
      IF (START) IWORK(21) = 1
      IWORK(22) = -1
      IF (PHASE1) IWORK(22) = 1
      IWORK(23) = -1
      IF (NORND) IWORK(23) = 1
      IWORK(24) = -1
      IF (STIFF) IWORK(24) = 1
      IWORK(25) = -1
      IF (INTOUT) IWORK(25) = 1
C
      IF (IDID .NE. (-2)) IWORK(LIW) = IWORK(LIW) + 1
      IF (T .NE. RWORK(ITSTAR)) IWORK(LIW) = 0
C
      RETURN
      END
*DECK DDES
      SUBROUTINE DDES (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,
     +   YPOUT, YP, YY, WT, P, PHI, ALPHA, BETA, PSI, V, W, SIG, G, GI,
     +   H, EPS, X, XOLD, HOLD, TOLD, DELSGN, TSTOP, TWOU, FOURU, START,
     +   PHASE1, NORND, STIFF, INTOUT, NS, KORD, KOLD, INIT, KSTEPS,
     +   KLE4, IQUIT, KPREV, IVC, IV, KGI, RPAR, IPAR)
C***BEGIN PROLOGUE  DDES
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEABM
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (DES-S, DDES-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DDEABM merely allocates storage for DDES to relieve the user of the
C   inconvenience of a long call list.  Consequently  DDES  is used as
C   described in the comments for  DDEABM .
C
C***SEE ALSO  DDEABM
C***ROUTINES CALLED  D1MACH, DINTP, DSTEPS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls, cvt GOTOs to
C           IF-THEN-ELSE.  (RWC)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DDES
C
      INTEGER IDID, INFO, INIT, IPAR, IQUIT, IV, IVC, K, KGI, KLE4,
     1      KOLD, KORD, KPREV, KSTEPS, L, LTOL, MAXNUM, NATOLP, NEQ,
     2      NRTOLP, NS
      DOUBLE PRECISION A, ABSDEL, ALPHA, ATOL, BETA, D1MACH,
     1      DEL, DELSGN, DT, EPS, FOURU, G, GI, H,
     2      HA, HOLD, P, PHI, PSI, RPAR, RTOL, SIG, T, TOLD, TOUT,
     3      TSTOP, TWOU, U, V, W, WT, X, XOLD, Y, YP, YPOUT, YY
      LOGICAL STIFF,CRASH,START,PHASE1,NORND,INTOUT
C
      DIMENSION Y(*),YY(*),WT(*),PHI(NEQ,16),P(*),YP(*),
     1  YPOUT(*),PSI(12),ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),
     2  GI(11),IV(10),INFO(15),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
C
      EXTERNAL DF
C
C.......................................................................
C
C  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
C  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER
C  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
C  WORK.
C
      SAVE MAXNUM
      DATA MAXNUM/500/
C
C.......................................................................
C
C***FIRST EXECUTABLE STATEMENT  DDES
      IF (INFO(1) .EQ. 0) THEN
C
C ON THE FIRST CALL , PERFORM INITIALIZATION --
C        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
C        FUNCTION ROUTINE  D1MACH. THE USER MUST MAKE SURE THAT THE
C        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
C
         U=D1MACH(4)
C                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
         TWOU=2.D0*U
         FOURU=4.D0*U
C                       -- SET TERMINATION FLAG
         IQUIT=0
C                       -- SET INITIALIZATION INDICATOR
         INIT=0
C                       -- SET COUNTER FOR ATTEMPTED STEPS
         KSTEPS=0
C                       -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
         INTOUT= .FALSE.
C                       -- SET INDICATOR FOR STIFFNESS DETECTION
         STIFF= .FALSE.
C                       -- SET STEP COUNTER FOR STIFFNESS DETECTION
         KLE4=0
C                       -- SET INDICATORS FOR STEPS CODE
         START= .TRUE.
         PHASE1= .TRUE.
         NORND= .TRUE.
C                       -- RESET INFO(1) FOR SUBSEQUENT CALLS
         INFO(1)=1
      ENDIF
C
C.......................................................................
C
C      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
C
      IF (INFO(1) .NE. 0  .AND.  INFO(1) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(1)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(1) MUST BE ' //
     *      'SET TO 0 FOR THE START OF A NEW PROBLEM, AND MUST BE ' //
     *      'SET TO 1 FOLLOWING AN INTERRUPTED TASK.  YOU ARE ' //
     *      'ATTEMPTING TO CONTINUE THE INTEGRATION ILLEGALLY BY ' //
     *      'CALLING THE CODE WITH INFO(1) = ' // XERN1, 3, 1)
         IDID=-33
      ENDIF
C
      IF (INFO(2) .NE. 0  .AND.  INFO(2) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(2)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(2) MUST BE ' //
     *      '0 OR 1 INDICATING SCALAR AND VECTOR ERROR TOLERANCES, ' //
     *      'RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH INFO(2) = ' //
     *      XERN1, 4, 1)
         IDID=-33
      ENDIF
C
      IF (INFO(3) .NE. 0  .AND.  INFO(3) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(3)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(3) MUST BE ' //
     *      '0 OR 1 INDICATING THE INTERVAL OR INTERMEDIATE-OUTPUT ' //
     *      'MODE OF INTEGRATION, RESPECTIVELY.  YOU HAVE CALLED ' //
     *      'THE CODE WITH  INFO(3) = ' // XERN1, 5, 1)
         IDID=-33
      ENDIF
C
      IF (INFO(4) .NE. 0  .AND.  INFO(4) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(4)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(4) MUST BE ' //
     *      '0 OR 1 INDICATING WHETHER OR NOT THE INTEGRATION ' //
     *      'INTERVAL IS TO BE RESTRICTED BY A POINT TSTOP.  YOU ' //
     *      'HAVE CALLED THE CODE WITH INFO(4) = ' // XERN1, 14, 1)
         IDID=-33
      ENDIF
C
      IF (NEQ .LT. 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM,  THE NUMBER OF ' //
     *      'EQUATIONS NEQ MUST BE A POSITIVE INTEGER.  YOU HAVE ' //
     *      'CALLED THE CODE WITH  NEQ = ' // XERN1, 6, 1)
         IDID=-33
      ENDIF
C
      NRTOLP = 0
      NATOLP = 0
      DO 90 K=1,NEQ
         IF (NRTOLP .EQ. 0 .AND. RTOL(K) .LT. 0.D0) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') RTOL(K)
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, THE RELATIVE ' //
     *         'ERROR TOLERANCES RTOL MUST BE NON-NEGATIVE.  YOU ' //
     *         'HAVE CALLED THE CODE WITH  RTOL(' // XERN1 // ') = ' //
     *         XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' //
     *         'NO FURTHER CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
            IDID = -33
            NRTOLP = 1
         ENDIF
C
         IF (NATOLP .EQ. 0 .AND. ATOL(K) .LT. 0.D0) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') ATOL(K)
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, THE ABSOLUTE ' //
     *         'ERROR TOLERANCES ATOL MUST BE NON-NEGATIVE.  YOU ' //
     *         'HAVE CALLED THE CODE WITH  ATOL(' // XERN1 // ') = ' //
     *         XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' //
     *         'NO FURTHER CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
            IDID = -33
            NATOLP = 1
         ENDIF
C
         IF (INFO(2) .EQ. 0) GO TO 100
         IF (NATOLP.GT.0 .AND. NRTOLP.GT.0) GO TO 100
   90 CONTINUE
C
  100 IF (INFO(4) .EQ. 1) THEN
         IF (SIGN(1.D0,TOUT-T) .NE. SIGN(1.D0,TSTOP-T)
     1      .OR. ABS(TOUT-T) .GT. ABS(TSTOP-T)) THEN
            WRITE (XERN3, '(1PE15.6)') TOUT
            WRITE (XERN4, '(1PE15.6)') TSTOP
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //
     *         'CALLED THE CODE WITH  TOUT = ' // XERN3 // ' BUT ' //
     *         'YOU HAVE ALSO TOLD THE CODE (INFO(4) = 1) NOT TO ' //
     *         'INTEGRATE PAST THE POINT TSTOP = ' // XERN4 //
     *         ' THESE INSTRUCTIONS CONFLICT.', 14, 1)
            IDID=-33
         ENDIF
      ENDIF
C
C     CHECK SOME CONTINUATION POSSIBILITIES
C
      IF (INIT .NE. 0) THEN
         IF (T .EQ. TOUT) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //
     *         'CALLED THE CODE WITH  T = TOUT = ' // XERN3 //
     *         '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 9, 1)
            IDID=-33
         ENDIF
C
         IF (T .NE. TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //
     *         'CHANGED THE VALUE OF T FROM ' // XERN3 // ' TO ' //
     *         XERN4 //'  THIS IS NOT ALLOWED ON CONTINUATION CALLS.',
     *         10, 1)
            IDID=-33
         ENDIF
C
         IF (INIT .NE. 1) THEN
            IF (DELSGN*(TOUT-T) .LT. 0.D0) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, BY ' //
     *            'CALLING THE CODE WITH TOUT = ' // XERN3 //
     *            ' YOU ARE ATTEMPTING TO CHANGE THE DIRECTION OF ' //
     *            'INTEGRATION.$$THIS IS NOT ALLOWED WITHOUT ' //
     *            'RESTARTING.', 11, 1)
               IDID=-33
            ENDIF
         ENDIF
      ENDIF
C
C     INVALID INPUT DETECTED
C
      IF (IDID .EQ. (-33)) THEN
         IF (IQUIT .NE. (-33)) THEN
            IQUIT = -33
            INFO(1) = -1
         ELSE
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INVALID ' //
     *         'INPUT WAS DETECTED ON SUCCESSIVE ENTRIES.  IT IS ' //
     *         'IMPOSSIBLE TO PROCEED BECAUSE YOU HAVE NOT ' //
     *         'CORRECTED THE PROBLEM, SO EXECUTION IS BEING ' //
     *         'TERMINATED.', 12, 2)
         ENDIF
         RETURN
      ENDIF
C
C.......................................................................
C
C     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS
C     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE,
C     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE
C     FOURU WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE
C
      DO 180 K=1,NEQ
        IF (RTOL(K)+ATOL(K) .GT. 0.D0) GO TO 170
        RTOL(K)=FOURU
        IDID=-2
  170   IF (INFO(2) .EQ. 0) GO TO 190
  180   CONTINUE
C
  190 IF (IDID .NE. (-2)) GO TO 200
C                       RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
C                                                SMALL POSITIVE VALUE
      INFO(1)=-1
      RETURN
C
C     BRANCH ON STATUS OF INITIALIZATION INDICATOR
C            INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE
C                   AND DIRECTION NOT YET SET
C            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET
C            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
C
  200 IF (INIT .EQ. 0) GO TO 210
      IF (INIT .EQ. 1) GO TO 220
      GO TO 240
C
C.......................................................................
C
C     MORE INITIALIZATION --
C                         -- EVALUATE INITIAL DERIVATIVES
C
  210 INIT=1
      A=T
      CALL DF(A,Y,YP,RPAR,IPAR)
      IF (T .NE. TOUT) GO TO 220
      IDID=2
      DO 215 L = 1,NEQ
  215    YPOUT(L) = YP(L)
      TOLD=T
      RETURN
C
C                         -- SET INDEPENDENT AND DEPENDENT VARIABLES
C                                              X AND YY(*) FOR STEPS
C                         -- SET SIGN OF INTEGRATION DIRECTION
C                         -- INITIALIZE THE STEP SIZE
C
  220 INIT = 2
      X = T
      DO 230 L = 1,NEQ
  230   YY(L) = Y(L)
      DELSGN = SIGN(1.0D0,TOUT-T)
      H = SIGN(MAX(FOURU*ABS(X),ABS(TOUT-X)),TOUT-X)
C
C.......................................................................
C
C   ON EACH CALL SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL
C   OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT
C
  240 DEL = TOUT - T
      ABSDEL = ABS(DEL)
C
C.......................................................................
C
C   IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN
C
  250 IF(ABS(X-T) .LT. ABSDEL) GO TO 260
      CALL DINTP(X,YY,TOUT,Y,YPOUT,NEQ,KOLD,PHI,IVC,IV,KGI,GI,
     1                                        ALPHA,G,W,XOLD,P)
      IDID = 3
      IF (X .NE. TOUT) GO TO 255
      IDID = 2
      INTOUT = .FALSE.
  255 T = TOUT
      TOLD = T
      RETURN
C
C   IF CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE,
C   EXTRAPOLATE AND RETURN
C
  260 IF (INFO(4) .NE. 1) GO TO 280
      IF (ABS(TSTOP-X) .GE. FOURU*ABS(X)) GO TO 280
      DT = TOUT - X
      DO 270 L = 1,NEQ
  270   Y(L) = YY(L) + DT*YP(L)
      CALL DF(TOUT,Y,YPOUT,RPAR,IPAR)
      IDID = 3
      T = TOUT
      TOLD = T
      RETURN
C
  280 IF (INFO(3) .EQ. 0  .OR.  .NOT.INTOUT) GO TO 300
C
C   INTERMEDIATE-OUTPUT MODE
C
      IDID = 1
      DO 290 L = 1,NEQ
        Y(L)=YY(L)
  290   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INTOUT = .FALSE.
      RETURN
C
C.......................................................................
C
C     MONITOR NUMBER OF STEPS ATTEMPTED
C
  300 IF (KSTEPS .LE. MAXNUM) GO TO 330
C
C                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED
      IDID=-1
      KSTEPS=0
      IF (.NOT. STIFF) GO TO 310
C
C                       PROBLEM APPEARS TO BE STIFF
      IDID=-4
      STIFF= .FALSE.
      KLE4=0
C
  310 DO 320 L = 1,NEQ
        Y(L) = YY(L)
  320   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
C
C.......................................................................
C
C   LIMIT STEP SIZE, SET WEIGHT VECTOR AND TAKE A STEP
C
  330 HA = ABS(H)
      IF (INFO(4) .NE. 1) GO TO 340
      HA = MIN(HA,ABS(TSTOP-X))
  340 H = SIGN(HA,H)
      EPS = 1.0D0
      LTOL = 1
      DO 350 L = 1,NEQ
        IF (INFO(2) .EQ. 1) LTOL = L
        WT(L) = RTOL(LTOL)*ABS(YY(L)) + ATOL(LTOL)
        IF (WT(L) .LE. 0.0D0) GO TO 360
  350   CONTINUE
      GO TO 380
C
C                       RELATIVE ERROR CRITERION INAPPROPRIATE
  360 IDID = -3
      DO 370 L = 1,NEQ
        Y(L) = YY(L)
  370   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
C
  380 CALL DSTEPS(DF,NEQ,YY,X,H,EPS,WT,START,HOLD,KORD,KOLD,CRASH,PHI,P,
     1           YP,PSI,ALPHA,BETA,SIG,V,W,G,PHASE1,NS,NORND,KSTEPS,
     2           TWOU,FOURU,XOLD,KPREV,IVC,IV,KGI,GI,RPAR,IPAR)
C
C.......................................................................
C
      IF(.NOT.CRASH) GO TO 420
C
C                       TOLERANCES TOO SMALL
      IDID = -2
      RTOL(1) = EPS*RTOL(1)
      ATOL(1) = EPS*ATOL(1)
      IF (INFO(2) .EQ. 0) GO TO 400
      DO 390 L = 2,NEQ
        RTOL(L) = EPS*RTOL(L)
  390   ATOL(L) = EPS*ATOL(L)
  400 DO 410 L = 1,NEQ
        Y(L) = YY(L)
  410   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
C
C   (STIFFNESS TEST) COUNT NUMBER OF CONSECUTIVE STEPS TAKEN WITH THE
C   ORDER OF THE METHOD BEING LESS OR EQUAL TO FOUR
C
  420 KLE4 = KLE4 + 1
      IF(KOLD .GT. 4) KLE4 = 0
      IF(KLE4 .GE. 50) STIFF = .TRUE.
      INTOUT = .TRUE.
      GO TO 250
      END
*DECK DINTP
      SUBROUTINE DINTP (X, Y, XOUT, YOUT, YPOUT, NEQN, KOLD, PHI, IVC,
     +   IV, KGI, GI, ALPHA, OG, OW, OX, OY)
C***BEGIN PROLOGUE  DINTP
C***PURPOSE  Approximate the solution at XOUT by evaluating the
C            polynomial computed in DSTEPS at XOUT.  Must be used in
C            conjunction with DSTEPS.
C***LIBRARY   SLATEC (DEPAC)
C***CATEGORY  I1A1B
C***TYPE      DOUBLE PRECISION (SINTRP-S, DINTP-D)
C***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
C             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR,
C             SMOOTH INTERPOLANT
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   The methods in subroutine  DSTEPS  approximate the solution near  X
C   by a polynomial.  Subroutine  DINTP  approximates the solution at
C   XOUT  by evaluating the polynomial there.  Information defining this
C   polynomial is passed from  DSTEPS  so  DINTP  cannot be used alone.
C
C   Subroutine DSTEPS is completely explained and documented in the text
C   "Computer Solution of Ordinary Differential Equations, the Initial
C   Value Problem"  by L. F. Shampine and M. K. Gordon.
C
C   Input to DINTP --
C
C   The user provides storage in the calling program for the arrays in
C   the call list
C      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),OY(NEQN)
C                AND ALPHA(12),OG(13),OW(12),GI(11),IV(10)
C   and defines
C      XOUT -- point at which solution is desired.
C   The remaining parameters are defined in  DSTEPS  and passed to
C   DINTP  from that subroutine
C
C   Output from  DINTP --
C
C      YOUT(*) -- solution at  XOUT
C      YPOUT(*) -- derivative of solution at  XOUT
C   The remaining parameters are returned unaltered from their input
C   values.  Integration with  DSTEPS  may be continued.
C
C***REFERENCES  H. A. Watts, A smoother interpolant for DE/STEP, INTRP
C                 II, Report SAND84-0293, Sandia Laboratories, 1984.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   840201  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DINTP
C
      INTEGER I, IQ, IV, IVC, IW, J, JQ, KGI, KOLD, KP1, KP2,
     1        L, M, NEQN
      DOUBLE PRECISION ALP, ALPHA, C, G, GDI, GDIF, GI, GAMMA, H, HI,
     1       HMU, OG, OW, OX, OY, PHI, RMU, SIGMA, TEMP1, TEMP2, TEMP3,
     2       W, X, XI, XIM1, XIQ, XOUT, Y, YOUT, YPOUT
C
      DIMENSION Y(*),YOUT(*),YPOUT(*),PHI(NEQN,16),OY(*)
      DIMENSION G(13),C(13),W(13),OG(13),OW(12),ALPHA(12),GI(11),IV(10)
C
C***FIRST EXECUTABLE STATEMENT  DINTP
      KP1 = KOLD + 1
      KP2 = KOLD + 2
C
      HI = XOUT - OX
      H = X - OX
      XI = HI/H
      XIM1 = XI - 1.D0
C
C   INITIALIZE W(*) FOR COMPUTING G(*)
C
      XIQ = XI
      DO 10 IQ = 1,KP1
        XIQ = XI*XIQ
        TEMP1 = IQ*(IQ+1)
 10     W(IQ) = XIQ/TEMP1
C
C   COMPUTE THE DOUBLE INTEGRAL TERM GDI
C
      IF (KOLD .LE. KGI) GO TO 50
      IF (IVC .GT. 0) GO TO 20
      GDI = 1.0D0/TEMP1
      M = 2
      GO TO 30
 20   IW = IV(IVC)
      GDI = OW(IW)
      M = KOLD - IW + 3
 30   IF (M .GT. KOLD) GO TO 60
      DO 40 I = M,KOLD
 40     GDI = OW(KP2-I) - ALPHA(I)*GDI
      GO TO 60
 50   GDI = GI(KOLD)
C
C   COMPUTE G(*) AND C(*)
C
 60   G(1) = XI
      G(2) = 0.5D0*XI*XI
      C(1) = 1.0D0
      C(2) = XI
      IF (KOLD .LT. 2) GO TO 90
      DO 80 I = 2,KOLD
        ALP = ALPHA(I)
        GAMMA = 1.0D0 + XIM1*ALP
        L = KP2 - I
        DO 70 JQ = 1,L
 70       W(JQ) = GAMMA*W(JQ) - ALP*W(JQ+1)
        G(I+1) = W(1)
 80     C(I+1) = GAMMA*C(I)
C
C   DEFINE INTERPOLATION PARAMETERS
C
 90   SIGMA = (W(2) - XIM1*W(1))/GDI
      RMU = XIM1*C(KP1)/GDI
      HMU = RMU/H
C
C   INTERPOLATE FOR THE SOLUTION -- YOUT
C   AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT
C
      DO 100 L = 1,NEQN
        YOUT(L) = 0.0D0
 100    YPOUT(L) = 0.0D0
      DO 120 J = 1,KOLD
        I = KP2 - J
        GDIF = OG(I) - OG(I-1)
        TEMP2 = (G(I) - G(I-1)) - SIGMA*GDIF
        TEMP3 = (C(I) - C(I-1)) + RMU*GDIF
        DO 110 L = 1,NEQN
          YOUT(L) = YOUT(L) + TEMP2*PHI(L,I)
 110      YPOUT(L) = YPOUT(L) + TEMP3*PHI(L,I)
 120    CONTINUE
      DO 130 L = 1,NEQN
        YOUT(L) = ((1.0D0 - SIGMA)*OY(L) + SIGMA*Y(L)) +
     1             H*(YOUT(L) + (G(1) - SIGMA*OG(1))*PHI(L,1))
 130    YPOUT(L) = HMU*(OY(L) - Y(L)) +
     1                (YPOUT(L) + (C(1) + RMU*OG(1))*PHI(L,1))
C
      RETURN
      END
*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
*DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
C***BEGIN PROLOGUE  XERPRN
C***SUBSIDIARY
C***PURPOSE  Print error messages processed by XERMSG.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERPRN-A)
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C This routine sends one or more lines to each of the (up to five)
C logical units to which error messages are to be sent.  This routine
C is called several times by XERMSG, sometimes with a single line to
C print and sometimes with a (potentially very long) message that may
C wrap around into multiple lines.
C
C PREFIX  Input argument of type CHARACTER.  This argument contains
C         characters to be put at the beginning of each line before
C         the body of the message.  No more than 16 characters of
C         PREFIX will be used.
C
C NPREF   Input argument of type INTEGER.  This argument is the number
C         of characters to use from PREFIX.  If it is negative, the
C         intrinsic function LEN is used to determine its length.  If
C         it is zero, PREFIX is not used.  If it exceeds 16 or if
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
C         used.  If NPREF is positive and the length of PREFIX is less
C         than NPREF, a copy of PREFIX extended with blanks to length
C         NPREF will be used.
C
C MESSG   Input argument of type CHARACTER.  This is the text of a
C         message to be printed.  If it is a long message, it will be
C         broken into pieces for printing on multiple lines.  Each line
C         will start with the appropriate prefix and be followed by a
C         piece of the message.  NWRAP is the number of characters per
C         piece; that is, after each NWRAP characters, we break and
C         start a new line.  In addition the characters '$$' embedded
C         in MESSG are a sentinel for a new line.  The counting of
C         characters up to NWRAP starts over for each new line.  The
C         value of NWRAP typically used by XERMSG is 72 since many
C         older error messages in the SLATEC Library are laid out to
C         rely on wrap-around every 72 characters.
C
C NWRAP   Input argument of type INTEGER.  This gives the maximum size
C         piece into which to break MESSG for printing on multiple
C         lines.  An embedded '$$' ends a line, and the count restarts
C         at the following character.  If a line break does not occur
C         on a blank (it would split a word) that word is moved to the
C         next line.  Values of NWRAP less than 16 will be treated as
C         16.  Values of NWRAP greater than 132 will be treated as 132.
C         The actual line length will be NPREF + NWRAP after NPREF has
C         been adjusted to fall between 0 and 16 and NWRAP has been
C         adjusted to fall between 16 and 132.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   880621  DATE WRITTEN
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
C           SLASH CHARACTER IN FORMAT STATEMENTS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
C           LINES TO BE PRINTED.
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Added code to break messages between words.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
*DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
C***BEGIN PROLOGUE  XERSVE
C***SUBSIDIARY
C***PURPOSE  Record that an error has occurred.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (XERSVE-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
C
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
C
C *Arguments:
C
C        LIBRAR :IN    is the library that the message is from.
C        SUBROU :IN    is the subroutine that the message is from.
C        MESSG  :IN    is the message to be saved.
C        KFLAG  :IN    indicates the action to be performed.
C                      when KFLAG > 0, the message in MESSG is saved.
C                      when KFLAG=0 the tables will be dumped and
C                      cleared.
C                      when KFLAG < 0, the tables will be dumped and
C                      not cleared.
C        NERR   :IN    is the error number.
C        LEVEL  :IN    is the error severity.
C        ICOUNT :OUT   the number of times this message has been seen,
C                      or zero if the table has overflowed and does not
C                      contain this message specifically.  When KFLAG=0,
C                      ICOUNT will not be altered.
C
C *Description:
C
C   Record that this error occurred and possibly dump and clear the
C   tables.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   800319  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900413  Routine modified to remove reference to KFLAG.  (WRB)
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
C           sequence, use IF-THEN-ELSE, make number of saved entries
C           easily changeable, changed routine name from XERSAV to
C           XERSVE.  (RWC)
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (I)
C***BEGIN PROLOGUE  D1MACH
C***PURPOSE  Return floating point machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   D1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument, and can be referenced as follows:
C
C        D = D1MACH(I)
C
C   where I=1,...,5.  The (output) value of D above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   Assume double precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(14) = T, the number of base-B digits.
C   I1MACH(15) = EMIN, the smallest exponent E.
C   I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890213  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   900911  Added SUN 386i constants.  (WRB)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added CONVEX -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C   010817  Elevated IEEE to highest importance; see next set of
C           comments below.  (DWL)
C***END PROLOGUE  D1MACH
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
C Initial data here correspond to the IEEE standard.  The values for
C DMACH(1), DMACH(3) and DMACH(4) are slight upper bounds.  The value
C for DMACH(2) is a slight lower bound.  The value for DMACH(5) is
C a 20-digit approximation.  If one of the sets of initial data below
C is preferred, do the necessary commenting and uncommenting. (DWL)
      DOUBLE PRECISION DMACH(5)
      DATA DMACH / 2.23D-308, 1.79D+308, 1.111D-16, 2.222D-16,
     1             0.30102999566398119521D0 /
      SAVE DMACH
C
cc      EQUIVALENCE (DMACH(1),SMALL(1))
cc      EQUIVALENCE (DMACH(2),LARGE(1))
cc      EQUIVALENCE (DMACH(3),RIGHT(1))
cc      EQUIVALENCE (DMACH(4),DIVER(1))
cc      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA SMALL(1) / ZC00800000 /
C     DATA SMALL(2) / Z000000000 /
C     DATA LARGE(1) / ZDFFFFFFFF /
C     DATA LARGE(2) / ZFFFFFFFFF /
C     DATA RIGHT(1) / ZCC5800000 /
C     DATA RIGHT(2) / Z000000000 /
C     DATA DIVER(1) / ZCC6800000 /
C     DATA DIVER(2) / Z000000000 /
C     DATA LOG10(1) / ZD00E730E7 /
C     DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O0000000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O0007777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O7770000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O7777777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA SMALL(1) / Z"3001800000000000" /
C     DATA SMALL(2) / Z"3001000000000000" /
C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
C     DATA LARGE(2) / Z"4FFE000000000000" /
C     DATA RIGHT(1) / Z"3FD2800000000000" /
C     DATA RIGHT(2) / Z"3FD2000000000000" /
C     DATA DIVER(1) / Z"3FD3800000000000" /
C     DATA DIVER(2) / Z"3FD3000000000000" /
C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA SMALL(1) / 00564000000000000000B /
C     DATA SMALL(2) / 00000000000000000000B /
C     DATA LARGE(1) / 37757777777777777777B /
C     DATA LARGE(2) / 37157777777777777777B /
C     DATA RIGHT(1) / 15624000000000000000B /
C     DATA RIGHT(2) / 00000000000000000000B /
C     DATA DIVER(1) / 15634000000000000000B /
C     DATA DIVER(2) / 00000000000000000000B /
C     DATA LOG10(1) / 17164642023241175717B /
C     DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn OR -pd8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CC0000000000000' /
C     DATA DMACH(4) / Z'3CD0000000000000' /
C     DATA DMACH(5) / Z'3FF34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3F900000000000000000000000000000' /
C     DATA DMACH(4) / Z'3F910000000000000000000000000000' /
C     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
C
C     MACHINE CONSTANTS FOR THE CRAY
C
C     DATA SMALL(1) / 201354000000000000000B /
C     DATA SMALL(2) / 000000000000000000000B /
C     DATA LARGE(1) / 577767777777777777777B /
C     DATA LARGE(2) / 000007777777777777774B /
C     DATA RIGHT(1) / 376434000000000000000B /
C     DATA RIGHT(2) / 000000000000000000000B /
C     DATA DIVER(1) / 376444000000000000000B /
C     DATA DIVER(2) / 000000000000000000000B /
C     DATA LOG10(1) / 377774642023241175717B /
C     DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC DMACH(5)
C
C     DATA SMALL /    20K, 3*0 /
C     DATA LARGE / 77777K, 3*177777K /
C     DATA RIGHT / 31420K, 3*0 /
C     DATA DIVER / 32020K, 3*0 /
C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA DMACH(1) / '0000000000000010'X /
C     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
C     DATA DMACH(3) / '0000000000003CC0'X /
C     DATA DMACH(4) / '0000000000003CD0'X /
C     DATA DMACH(5) / '79FF509F44133FF3'X /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FORMAT
C
C     DATA DMACH(1) / '0010000000000000'X /
C     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X /
C     DATA DMACH(3) / '3CA0000000000000'X /
C     DATA DMACH(4) / '3CB0000000000000'X /
C     DATA DMACH(5) / '3FD34413509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
C     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
C     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
C     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /        128,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C     DATA DIVER(1), DIVER(2) /       9472,           0 /
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /         16,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C     DATA DIVER(1), DIVER(2) /      15568,           0 /
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) /  40000B,       0 /
C     DATA SMALL(3), SMALL(4) /       0,       1 /
C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
C     DATA RIGHT(3), RIGHT(4) /       0,    225B /
C     DATA DIVER(1), DIVER(2) /  40000B,       0 /
C     DATA DIVER(3), DIVER(4) /       0,    227B /
C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
C
C     DATA SMALL(1) / 2.23D-308  /
C     DATA LARGE(1) / 1.79D+308  /
C     DATA RIGHT(1) / 1.11D-16   /
C     DATA DIVER(1) / 2.22D-16   /
C     DATA LOG10(1) / 0.301029995663981195D0 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    8388608,           0 /
C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
C     DATA DIVER(1), DIVER(2) /  620756992,           0 /
C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
C
C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    128,      0 /
C     DATA SMALL(3), SMALL(4) /      0,      0 /
C     DATA LARGE(1), LARGE(2) /  32767,     -1 /
C     DATA LARGE(3), LARGE(4) /     -1,     -1 /
C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
C     DATA RIGHT(3), RIGHT(4) /      0,      0 /
C     DATA DIVER(1), DIVER(2) /   9472,      0 /
C     DATA DIVER(3), DIVER(4) /      0,      0 /
C     DATA LOG10(1), LOG10(2) /  16282,   8346 /
C     DATA LOG10(3), LOG10(4) / -31493, -12296 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA SMALL(3), SMALL(4) / O000000, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA LARGE(3), LARGE(4) / O177777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
C     DATA DIVER(1), DIVER(2) / O022400, O000000 /
C     DATA DIVER(3), DIVER(4) / O000000, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020232 /
C     DATA LOG10(3), LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
C     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
C     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
C     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
C
C     MACHINE CONSTANTS FOR THE SUN 386i
C
C     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
C     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
C     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'D1MACH',
     +   'I OUT OF BOUNDS', 1, 2)
C
      D1MACH = DMACH(I)
      RETURN
C
      END
*DECK XGETUA
      SUBROUTINE XGETUA (IUNITA, N)
C***BEGIN PROLOGUE  XGETUA
C***PURPOSE  Return unit number(s) to which error messages are being
C            sent.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XGETUA-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
*DECK DSTEPS
      SUBROUTINE DSTEPS (DF, NEQN, Y, X, H, EPS, WT, START, HOLD, K,
     +   KOLD, CRASH, PHI, P, YP, PSI, ALPHA, BETA, SIG, V, W, G,
     +   PHASE1, NS, NORND, KSTEPS, TWOU, FOURU, XOLD, KPREV, IVC, IV,
     +   KGI, GI, RPAR, IPAR)
C***BEGIN PROLOGUE  DSTEPS
C***PURPOSE  Integrate a system of first order ordinary differential
C            equations one step.
C***LIBRARY   SLATEC (DEPAC)
C***CATEGORY  I1A1B
C***TYPE      DOUBLE PRECISION (STEPS-S, DSTEPS-D)
C***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
C             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
C***AUTHOR  Shampine, L. F., (SNLA)
C           Gordon, M. K., (SNLA)
C             MODIFIED BY H.A. WATTS
C***DESCRIPTION
C
C   Written by L. F. Shampine and M. K. Gordon
C
C   Abstract
C
C   Subroutine  DSTEPS  is normally used indirectly through subroutine
C   DDEABM .  Because  DDEABM  suffices for most problems and is much
C   easier to use, using it should be considered before using  DSTEPS
C   alone.
C
C   Subroutine DSTEPS integrates a system of  NEQN  first order ordinary
C   differential equations one step, normally from X to X+H, using a
C   modified divided difference form of the Adams Pece formulas.  Local
C   extrapolation is used to improve absolute stability and accuracy.
C   The code adjusts its order and step size to control the local error
C   per unit step in a generalized sense.  Special devices are included
C   to control roundoff error and to detect when the user is requesting
C   too much accuracy.
C
C   This code is completely explained and documented in the text,
C   Computer Solution of Ordinary Differential Equations, The Initial
C   Value Problem  by L. F. Shampine and M. K. Gordon.
C   Further details on use of this code are available in "Solving
C   Ordinary Differential Equations with ODE, STEP, and INTRP",
C   by L. F. Shampine and M. K. Gordon, SLA-73-1060.
C
C
C   The parameters represent --
C      DF -- subroutine to evaluate derivatives
C      NEQN -- number of equations to be integrated
C      Y(*) -- solution vector at X
C      X -- independent variable
C      H -- appropriate step size for next step.  Normally determined by
C           code
C      EPS -- local error tolerance
C      WT(*) -- vector of weights for error criterion
C      START -- logical variable set .TRUE. for first step,  .FALSE.
C           otherwise
C      HOLD -- step size used for last successful step
C      K -- appropriate order for next step (determined by code)
C      KOLD -- order used for last successful step
C      CRASH -- logical variable set .TRUE. when no step can be taken,
C           .FALSE. otherwise.
C      YP(*) -- derivative of solution vector at  X  after successful
C           step
C      KSTEPS -- counter on attempted steps
C      TWOU -- 2.*U where U is machine unit roundoff quantity
C      FOURU -- 4.*U where U is machine unit roundoff quantity
C      RPAR,IPAR -- parameter arrays which you may choose to use
C            for communication between your program and subroutine F.
C            They are not altered or used by DSTEPS.
C   The variables X,XOLD,KOLD,KGI and IVC and the arrays Y,PHI,ALPHA,G,
C   W,P,IV and GI are required for the interpolation subroutine SINTRP.
C   The remaining variables and arrays are included in the call list
C   only to eliminate local retention of variables between calls.
C
C   Input to DSTEPS
C
C      First call --
C
C   The user must provide storage in his calling program for all arrays
C   in the call list, namely
C
C     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12),
C    1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),
C    2  RPAR(*),IPAR(*)
C
C    **Note**
C
C   The user must also declare  START ,  CRASH ,  PHASE1  and  NORND
C   logical variables and  DF  an EXTERNAL subroutine, supply the
C   subroutine  DF(X,Y,YP)  to evaluate
C      DY(I)/DX = YP(I) = DF(X,Y(1),Y(2),...,Y(NEQN))
C   and initialize only the following parameters.
C      NEQN -- number of equations to be integrated
C      Y(*) -- vector of initial values of dependent variables
C      X -- initial value of the independent variable
C      H -- nominal step size indicating direction of integration
C           and maximum size of step.  Must be variable
C      EPS -- local error tolerance per step.  Must be variable
C      WT(*) -- vector of non-zero weights for error criterion
C      START -- .TRUE.
C      YP(*) -- vector of initial derivative values
C      KSTEPS -- set KSTEPS to zero
C      TWOU -- 2.*U where U is machine unit roundoff quantity
C      FOURU -- 4.*U where U is machine unit roundoff quantity
C   Define U to be the machine unit roundoff quantity by calling
C   the function routine  D1MACH,  U = D1MACH(4), or by
C   computing U so that U is the smallest positive number such
C   that 1.0+U .GT. 1.0.
C
C   DSTEPS  requires that the L2 norm of the vector with components
C   LOCAL ERROR(L)/WT(L)  be less than  EPS  for a successful step.  The
C   array  WT  allows the user to specify an error test appropriate
C   for his problem.  For example,
C      WT(L) = 1.0  specifies absolute error,
C            = ABS(Y(L))  error relative to the most recent value of the
C                 L-th component of the solution,
C            = ABS(YP(L))  error relative to the most recent value of
C                 the L-th component of the derivative,
C            = MAX(WT(L),ABS(Y(L)))  error relative to the largest
C                 magnitude of L-th component obtained so far,
C            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  specifies a mixed
C                 relative-absolute test where  RELERR  is relative
C                 error,  ABSERR  is absolute error and  EPS =
C                 MAX(RELERR,ABSERR) .
C
C      Subsequent calls --
C
C   Subroutine  DSTEPS  is designed so that all information needed to
C   continue the integration, including the step size  H  and the order
C   K , is returned with each step.  With the exception of the step
C   size, the error tolerance, and the weights, none of the parameters
C   should be altered.  The array  WT  must be updated after each step
C   to maintain relative error tests like those above.  Normally the
C   integration is continued just beyond the desired endpoint and the
C   solution interpolated there with subroutine  SINTRP .  If it is
C   impossible to integrate beyond the endpoint, the step size may be
C   reduced to hit the endpoint since the code will not take a step
C   larger than the  H  input.  Changing the direction of integration,
C   i.e., the sign of  H , requires the user set  START = .TRUE. before
C   calling  DSTEPS  again.  This is the only situation in which  START
C   should be altered.
C
C   Output from DSTEPS
C
C      Successful Step --
C
C   The subroutine returns after each successful step with  START  and
C   CRASH  set .FALSE. .  X  represents the independent variable
C   advanced one step of length  HOLD  from its value on input and  Y
C   the solution vector at the new value of  X .  All other parameters
C   represent information corresponding to the new  X  needed to
C   continue the integration.
C
C      Unsuccessful Step --
C
C   When the error tolerance is too small for the machine precision,
C   the subroutine returns without taking a step and  CRASH = .TRUE. .
C   An appropriate step size and error tolerance for continuing are
C   estimated and all other information is restored as upon input
C   before returning.  To continue with the larger tolerance, the user
C   just calls the code again.  A restart is neither required nor
C   desirable.
C
C***REFERENCES  L. F. Shampine and M. K. Gordon, Solving ordinary
C                 differential equations with ODE, STEP, and INTRP,
C                 Report SLA-73-1060, Sandia Laboratories, 1973.
C***ROUTINES CALLED  D1MACH, DHSTRT
C***REVISION HISTORY  (YYMMDD)
C   740101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSTEPS
C
      INTEGER I, IFAIL, IM1, IP1, IPAR, IQ, J, K, KM1, KM2, KNEW,
     1      KOLD, KP1, KP2, KSTEPS, L, LIMIT1, LIMIT2, NEQN, NS, NSM2,
     2      NSP1, NSP2
      DOUBLE PRECISION ABSH, ALPHA, BETA, BIG, D1MACH,
     1      EPS, ERK, ERKM1, ERKM2, ERKP1, ERR,
     2      FOURU, G, GI, GSTR, H, HNEW, HOLD, P, P5EPS, PHI, PSI, R,
     3      REALI, REALNS, RHO, ROUND, RPAR, SIG, TAU, TEMP1,
     4      TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, TWO, TWOU, U, V, W, WT,
     5      X, XOLD, Y, YP
      LOGICAL START,CRASH,PHASE1,NORND
      DIMENSION Y(*),WT(*),PHI(NEQN,16),P(*),YP(*),PSI(12),
     1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),
     2  RPAR(*),IPAR(*)
      DIMENSION TWO(13),GSTR(13)
      EXTERNAL DF
      SAVE TWO, GSTR
C
      DATA TWO(1),TWO(2),TWO(3),TWO(4),TWO(5),TWO(6),TWO(7),TWO(8),
     1     TWO(9),TWO(10),TWO(11),TWO(12),TWO(13)
     2     /2.0D0,4.0D0,8.0D0,16.0D0,32.0D0,64.0D0,128.0D0,256.0D0,
     3      512.0D0,1024.0D0,2048.0D0,4096.0D0,8192.0D0/
      DATA GSTR(1),GSTR(2),GSTR(3),GSTR(4),GSTR(5),GSTR(6),GSTR(7),
     1     GSTR(8),GSTR(9),GSTR(10),GSTR(11),GSTR(12),GSTR(13)
     2     /0.5D0,0.0833D0,0.0417D0,0.0264D0,0.0188D0,0.0143D0,0.0114D0,
     3      0.00936D0,0.00789D0,0.00679D0,0.00592D0,0.00524D0,0.00468D0/
C
C       ***     BEGIN BLOCK 0     ***
C   CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
C   PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
C   STARTING STEP SIZE.
C                   ***
C
C   IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
C
C***FIRST EXECUTABLE STATEMENT  DSTEPS
      CRASH = .TRUE.
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 5
      H = SIGN(FOURU*ABS(X),H)
      RETURN
 5    P5EPS = 0.5D0*EPS
C
C   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE
C
      ROUND = 0.0D0
      DO 10 L = 1,NEQN
 10     ROUND = ROUND + (Y(L)/WT(L))**2
      ROUND = TWOU*SQRT(ROUND)
      IF(P5EPS .GE. ROUND) GO TO 15
      EPS = 2.0D0*ROUND*(1.0D0 + FOURU)
      RETURN
 15   CRASH = .FALSE.
      G(1) = 1.0D0
      G(2) = 0.5D0
      SIG(1) = 1.0D0
      IF(.NOT.START) GO TO 99
C
C   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP
C
C     CALL DF(X,Y,YP,RPAR,IPAR)
C     SUM = 0.0
      DO 20 L = 1,NEQN
        PHI(L,1) = YP(L)
   20   PHI(L,2) = 0.0D0
C20     SUM = SUM + (YP(L)/WT(L))**2
C     SUM = SQRT(SUM)
C     ABSH = ABS(H)
C     IF(EPS .LT. 16.0*SUM*H*H) ABSH = 0.25*SQRT(EPS/SUM)
C     H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
C
      U = D1MACH(4)
      BIG = SQRT(D1MACH(2))
      CALL DHSTRT(DF,NEQN,X,X+H,Y,YP,WT,1,U,BIG,
     1             PHI(1,3),PHI(1,4),PHI(1,5),PHI(1,6),RPAR,IPAR,H)
C
      HOLD = 0.0D0
      K = 1
      KOLD = 0
      KPREV = 0
      START = .FALSE.
      PHASE1 = .TRUE.
      NORND = .TRUE.
      IF(P5EPS .GT. 100.0D0*ROUND) GO TO 99
      NORND = .FALSE.
      DO 25 L = 1,NEQN
 25     PHI(L,15) = 0.0D0
 99   IFAIL = 0
C       ***     END BLOCK 0     ***
C
C       ***     BEGIN BLOCK 1     ***
C   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
C   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED.
C                   ***
C
 100  KP1 = K+1
      KP2 = K+2
      KM1 = K-1
      KM2 = K-2
C
C   NS IS THE NUMBER OF DSTEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
C   ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE
C
      IF(H .NE. HOLD) NS = 0
      IF (NS.LE.KOLD) NS = NS+1
      NSP1 = NS+1
      IF (K .LT. NS) GO TO 199
C
C   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
C   ARE CHANGED
C
      BETA(NS) = 1.0D0
      REALNS = NS
      ALPHA(NS) = 1.0D0/REALNS
      TEMP1 = H*REALNS
      SIG(NSP1) = 1.0D0
      IF(K .LT. NSP1) GO TO 110
      DO 105 I = NSP1,K
        IM1 = I-1
        TEMP2 = PSI(IM1)
        PSI(IM1) = TEMP1
        BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
        TEMP1 = TEMP2 + H
        ALPHA(I) = H/TEMP1
        REALI = I
 105    SIG(I+1) = REALI*ALPHA(I)*SIG(I)
 110  PSI(K) = TEMP1
C
C   COMPUTE COEFFICIENTS G(*)
C
C   INITIALIZE V(*) AND SET W(*).
C
      IF(NS .GT. 1) GO TO 120
      DO 115 IQ = 1,K
        TEMP3 = IQ*(IQ+1)
        V(IQ) = 1.0D0/TEMP3
 115    W(IQ) = V(IQ)
      IVC = 0
      KGI = 0
      IF (K .EQ. 1) GO TO 140
      KGI = 1
      GI(1) = W(2)
      GO TO 140
C
C   IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*)
C
 120  IF(K .LE. KPREV) GO TO 130
      IF (IVC .EQ. 0) GO TO 122
      JV = KP1 - IV(IVC)
      IVC = IVC - 1
      GO TO 123
 122  JV = 1
      TEMP4 = K*KP1
      V(K) = 1.0D0/TEMP4
      W(K) = V(K)
      IF (K .NE. 2) GO TO 123
      KGI = 1
      GI(1) = W(2)
 123  NSM2 = NS-2
      IF(NSM2 .LT. JV) GO TO 130
      DO 125 J = JV,NSM2
        I = K-J
        V(I) = V(I) - ALPHA(J+1)*V(I+1)
 125    W(I) = V(I)
      IF (I .NE. 2) GO TO 130
      KGI = NS - 1
      GI(KGI) = W(2)
C
C   UPDATE V(*) AND SET W(*)
C
 130  LIMIT1 = KP1 - NS
      TEMP5 = ALPHA(NS)
      DO 135 IQ = 1,LIMIT1
        V(IQ) = V(IQ) - TEMP5*V(IQ+1)
 135    W(IQ) = V(IQ)
      G(NSP1) = W(1)
      IF (LIMIT1 .EQ. 1) GO TO 137
      KGI = NS
      GI(KGI) = W(2)
 137  W(LIMIT1+1) = V(LIMIT1+1)
      IF (K .GE. KOLD) GO TO 140
      IVC = IVC + 1
      IV(IVC) = LIMIT1 + 2
C
C   COMPUTE THE G(*) IN THE WORK VECTOR W(*)
C
 140  NSP2 = NS + 2
      KPREV = K
      IF(KP1 .LT. NSP2) GO TO 199
      DO 150 I = NSP2,KP1
        LIMIT2 = KP2 - I
        TEMP6 = ALPHA(I-1)
        DO 145 IQ = 1,LIMIT2
 145      W(IQ) = W(IQ) - TEMP6*W(IQ+1)
 150    G(I) = W(1)
 199    CONTINUE
C       ***     END BLOCK 1     ***
C
C       ***     BEGIN BLOCK 2     ***
C   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED
C   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K,
C   K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED.
C                   ***
C
C   INCREMENT COUNTER ON ATTEMPTED DSTEPS
C
      KSTEPS = KSTEPS + 1
C
C   CHANGE PHI TO PHI STAR
C
      IF(K .LT. NSP1) GO TO 215
      DO 210 I = NSP1,K
        TEMP1 = BETA(I)
        DO 205 L = 1,NEQN
 205      PHI(L,I) = TEMP1*PHI(L,I)
 210    CONTINUE
C
C   PREDICT SOLUTION AND DIFFERENCES
C
 215  DO 220 L = 1,NEQN
        PHI(L,KP2) = PHI(L,KP1)
        PHI(L,KP1) = 0.0D0
 220    P(L) = 0.0D0
      DO 230 J = 1,K
        I = KP1 - J
        IP1 = I+1
        TEMP2 = G(I)
        DO 225 L = 1,NEQN
          P(L) = P(L) + TEMP2*PHI(L,I)
 225      PHI(L,I) = PHI(L,I) + PHI(L,IP1)
 230    CONTINUE
      IF(NORND) GO TO 240
      DO 235 L = 1,NEQN
        TAU = H*P(L) - PHI(L,15)
        P(L) = Y(L) + TAU
 235    PHI(L,16) = (P(L) - Y(L)) - TAU
      GO TO 250
 240  DO 245 L = 1,NEQN
 245    P(L) = Y(L) + H*P(L)
 250  XOLD = X
      X = X + H
      ABSH = ABS(H)
      CALL DF(X,P,YP,RPAR,IPAR)
C
C   ESTIMATE ERRORS AT ORDERS K,K-1,K-2
C
      ERKM2 = 0.0D0
      ERKM1 = 0.0D0
      ERK = 0.0D0
      DO 265 L = 1,NEQN
        TEMP3 = 1.0D0/WT(L)
        TEMP4 = YP(L) - PHI(L,1)
        IF(KM2)265,260,255
 255    ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2
 260    ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2
 265    ERK = ERK + (TEMP4*TEMP3)**2
      IF(KM2)280,275,270
 270  ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2)
 275  ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1)
 280  TEMP5 = ABSH*SQRT(ERK)
      ERR = TEMP5*(G(K)-G(KP1))
      ERK = TEMP5*SIG(KP1)*GSTR(K)
      KNEW = K
C
C   TEST IF ORDER SHOULD BE LOWERED
C
      IF(KM2)299,290,285
 285  IF(MAX(ERKM1,ERKM2) .LE. ERK) KNEW = KM1
      GO TO 299
 290  IF(ERKM1 .LE. 0.5D0*ERK) KNEW = KM1
C
C   TEST IF STEP SUCCESSFUL
C
 299  IF(ERR .LE. EPS) GO TO 400
C       ***     END BLOCK 2     ***
C
C       ***     BEGIN BLOCK 3     ***
C   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) .
C   IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE
C   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
C   TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE
C   PRECISION.
C                   ***
C
C   RESTORE X, PHI(*,*) AND PSI(*)
C
      PHASE1 = .FALSE.
      X = XOLD
      DO 310 I = 1,K
        TEMP1 = 1.0D0/BETA(I)
        IP1 = I+1
        DO 305 L = 1,NEQN
 305      PHI(L,I) = TEMP1*(PHI(L,I) - PHI(L,IP1))
 310    CONTINUE
      IF(K .LT. 2) GO TO 320
      DO 315 I = 2,K
 315    PSI(I-1) = PSI(I) - H
C
C   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP
C   SIZE
C
 320  IFAIL = IFAIL + 1
      TEMP2 = 0.5D0
      IF(IFAIL - 3) 335,330,325
 325  IF(P5EPS .LT. 0.25D0*ERK) TEMP2 = SQRT(P5EPS/ERK)
 330  KNEW = 1
 335  H = TEMP2*H
      K = KNEW
      NS = 0
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 340
      CRASH = .TRUE.
      H = SIGN(FOURU*ABS(X),H)
      EPS = EPS + EPS
      RETURN
 340  GO TO 100
C       ***     END BLOCK 3     ***
C
C       ***     BEGIN BLOCK 4     ***
C   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE
C   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE
C   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP.
C                   ***
 400  KOLD = K
      HOLD = H
C
C   CORRECT AND EVALUATE
C
      TEMP1 = H*G(KP1)
      IF(NORND) GO TO 410
      DO 405 L = 1,NEQN
        TEMP3 = Y(L)
        RHO = TEMP1*(YP(L) - PHI(L,1)) - PHI(L,16)
        Y(L) = P(L) + RHO
        PHI(L,15) = (Y(L) - P(L)) - RHO
 405    P(L) = TEMP3
      GO TO 420
 410  DO 415 L = 1,NEQN
        TEMP3 = Y(L)
        Y(L) = P(L) + TEMP1*(YP(L) - PHI(L,1))
 415    P(L) = TEMP3
 420  CALL DF(X,Y,YP,RPAR,IPAR)
C
C   UPDATE DIFFERENCES FOR NEXT STEP
C
      DO 425 L = 1,NEQN
        PHI(L,KP1) = YP(L) - PHI(L,1)
 425    PHI(L,KP2) = PHI(L,KP1) - PHI(L,KP2)
      DO 435 I = 1,K
        DO 430 L = 1,NEQN
 430      PHI(L,I) = PHI(L,I) + PHI(L,KP1)
 435    CONTINUE
C
C   ESTIMATE ERROR AT ORDER K+1 UNLESS:
C     IN FIRST PHASE WHEN ALWAYS RAISE ORDER,
C     ALREADY DECIDED TO LOWER ORDER,
C     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE
C
      ERKP1 = 0.0D0
      IF(KNEW .EQ. KM1  .OR.  K .EQ. 12) PHASE1 = .FALSE.
      IF(PHASE1) GO TO 450
      IF(KNEW .EQ. KM1) GO TO 455
      IF(KP1 .GT. NS) GO TO 460
      DO 440 L = 1,NEQN
 440    ERKP1 = ERKP1 + (PHI(L,KP2)/WT(L))**2
      ERKP1 = ABSH*GSTR(KP1)*SQRT(ERKP1)
C
C   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER
C   FOR NEXT STEP
C
      IF(K .GT. 1) GO TO 445
      IF(ERKP1 .GE. 0.5D0*ERK) GO TO 460
      GO TO 450
 445  IF(ERKM1 .LE. MIN(ERK,ERKP1)) GO TO 455
      IF(ERKP1 .GE. ERK  .OR.  K .EQ. 12) GO TO 460
C
C   HERE ERKP1 .LT. ERK .LT. MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE
C   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
C
C   RAISE ORDER
C
 450  K = KP1
      ERK = ERKP1
      GO TO 460
C
C   LOWER ORDER
C
 455  K = KM1
      ERK = ERKM1
C
C   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
C
 460  HNEW = H + H
      IF(PHASE1) GO TO 465
      IF(P5EPS .GE. ERK*TWO(K+1)) GO TO 465
      HNEW = H
      IF(P5EPS .GE. ERK) GO TO 465
      TEMP2 = K+1
      R = (P5EPS/ERK)**(1.0D0/TEMP2)
      HNEW = ABSH*MAX(0.5D0,MIN(0.9D0,R))
      HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H)
 465  H = HNEW
      RETURN
C       ***     END BLOCK 4     ***
      END
*DECK FDUMP
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***PURPOSE  Symbolic dump (should be locally written).
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (FDUMP-A)
C***KEYWORDS  ERROR, XERMSG
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
*DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
C***BEGIN PROLOGUE  I1MACH
C***PURPOSE  Return integer machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      INTEGER (I1MACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   I1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument and can be referenced as follows:
C
C        K = I1MACH(I)
C
C   where I=1,...,16.  The (output) value of K above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   I/O unit numbers:
C     I1MACH( 1) = the standard input unit.
C     I1MACH( 2) = the standard output unit.
C     I1MACH( 3) = the standard punch unit.
C     I1MACH( 4) = the standard error message unit.
C
C   Words:
C     I1MACH( 5) = the number of bits per integer storage unit.
C     I1MACH( 6) = the number of characters per integer storage unit.
C
C   Integers:
C     assume integers are represented in the S-digit, base-A form
C
C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C     I1MACH( 7) = A, the base.
C     I1MACH( 8) = S, the number of base-A digits.
C     I1MACH( 9) = A**S - 1, the largest magnitude.
C
C   Floating-Point Numbers:
C     Assume floating-point numbers are represented in the T-digit,
C     base-B form
C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C                where 0 .LE. X(I) .LT. B for I=1,...,T,
C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C     I1MACH(10) = B, the base.
C
C   Single-Precision:
C     I1MACH(11) = T, the number of base-B digits.
C     I1MACH(12) = EMIN, the smallest exponent E.
C     I1MACH(13) = EMAX, the largest exponent E.
C
C   Double-Precision:
C     I1MACH(14) = T, the number of base-B digits.
C     I1MACH(15) = EMIN, the smallest exponent E.
C     I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   891012  Added VAX G-floating constants.  (WRB)
C   891012  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
C           (RWC)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added Convex -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
C           options.  (DWL, RWC and WRB).
C   010817  Elevated IEEE to highest importance; see next set of
C           comments below.  (DWL)
C***END PROLOGUE  I1MACH
C
C Initial data here correspond to the IEEE standard.  If one of the
C sets of initial data below is preferred, do the necessary commenting
C and uncommenting. (DWL)
      INTEGER IMACH(16),OUTPUT
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          6 /
      DATA IMACH( 4) /          6 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -126 /
      DATA IMACH(13) /        127 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1022 /
      DATA IMACH(16) /       1023 /
      SAVE IMACH
cc      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        129 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /          7 /
C     DATA IMACH( 2) /          2 /
C     DATA IMACH( 3) /          2 /
C     DATA IMACH( 4) /          2 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -256 /
C     DATA IMACH(13) /        255 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /       -256 /
C     DATA IMACH(16) /        255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /        -50 /
C     DATA IMACH(16) /         76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /     -32754 /
C     DATA IMACH(16) /      32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -4095 /
C     DATA IMACH(13) /       4094 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -4095 /
C     DATA IMACH(16) /       4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /    6LOUTPUT/
C     DATA IMACH( 5) /         60 /
C     DATA IMACH( 6) /         10 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /       -929 /
C     DATA IMACH(13) /       1070 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /       -929 /
C     DATA IMACH(16) /       1069 /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / Z'7FFFFFFF' /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16383 /
C     DATA IMACH(16) /      16383 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -pd8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 46 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         46 /
C     DATA IMACH( 9) / 1777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 64 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 777777777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /         11 /
C     DATA IMACH( 2) /         12 /
C     DATA IMACH( 3) /          8 /
C     DATA IMACH( 4) /         10 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         24 /
C     DATA IMACH( 6) /          3 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         23 /
C     DATA IMACH( 9) /    8388607 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         38 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /         43 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         63 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         39 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         55 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          7 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1015 /
C     DATA IMACH(16) /       1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) /  Z7FFFFFFF /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         54 /
C     DATA IMACH(15) /       -101 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         62 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1021 /
C     DATA IMACH(13) /       1024 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16381 /
C     DATA IMACH(16) /      16384 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          1 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /      -1024 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /          1 /
C     DATA IMACH( 2) /          1 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
      STOP
      END
*DECK DHSTRT
      SUBROUTINE DHSTRT (DF, NEQ, A, B, Y, YPRIME, ETOL, MORDER, SMALL,
     +   BIG, SPY, PV, YP, SF, RPAR, IPAR, H)
C***BEGIN PROLOGUE  DHSTRT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (HSTART-S, DHSTRT-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DHSTRT computes a starting step size to be used in solving initial
C   value problems in ordinary differential equations.
C
C **********************************************************************
C  ABSTRACT
C
C     Subroutine DHSTRT computes a starting step size to be used by an
C     initial value method in solving ordinary differential equations.
C     It is based on an estimate of the local Lipschitz constant for the
C     differential equation   (lower bound on a norm of the Jacobian) ,
C     a bound on the differential equation  (first derivative) , and
C     a bound on the partial derivative of the equation with respect to
C     the independent variable.
C     (all approximated near the initial point A)
C
C     Subroutine DHSTRT uses a function subprogram DHVNRM for computing
C     a vector norm. The maximum norm is presently utilized though it
C     can easily be replaced by any other vector norm. It is presumed
C     that any replacement norm routine would be carefully coded to
C     prevent unnecessary underflows or overflows from occurring, and
C     also, would not alter the vector or number of components.
C
C **********************************************************************
C  On input you must provide the following
C
C      DF -- This is a subroutine of the form
C                               DF(X,U,UPRIME,RPAR,IPAR)
C             which defines the system of first order differential
C             equations to be solved. For the given values of X and the
C             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
C             evaluate the NEQ components of the system of differential
C             equations  DU/DX=DF(X,U)  and store the derivatives in the
C             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
C             equations I=1,...,NEQ.
C
C             Subroutine DF must not alter X or U(*). You must declare
C             the name DF in an external statement in your program that
C             calls DHSTRT. You must dimension U and UPRIME in DF.
C
C             RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
C             arrays which you can use for communication between your
C             program and subroutine DF. They are not used or altered by
C             DHSTRT. If you do not need RPAR or IPAR, ignore these
C             parameters by treating them as dummy arguments. If you do
C             choose to use them, dimension them in your program and in
C             DF as arrays of appropriate length.
C
C      NEQ -- This is the number of (first order) differential equations
C             to be integrated.
C
C      A -- This is the initial point of integration.
C
C      B -- This is a value of the independent variable used to define
C             the direction of integration. A reasonable choice is to
C             set  B  to the first point at which a solution is desired.
C             You can also use  B, if necessary, to restrict the length
C             of the first integration step because the algorithm will
C             not compute a starting step length which is bigger than
C             ABS(B-A), unless  B  has been chosen too close to  A.
C             (it is presumed that DHSTRT has been called with  B
C             different from  A  on the machine being used. Also see the
C             discussion about the parameter  SMALL.)
C
C      Y(*) -- This is the vector of initial values of the NEQ solution
C             components at the initial point  A.
C
C      YPRIME(*) -- This is the vector of derivatives of the NEQ
C             solution components at the initial point  A.
C             (defined by the differential equations in subroutine DF)
C
C      ETOL -- This is the vector of error tolerances corresponding to
C             the NEQ solution components. It is assumed that all
C             elements are positive. Following the first integration
C             step, the tolerances are expected to be used by the
C             integrator in an error test which roughly requires that
C                        ABS(LOCAL ERROR)  .LE.  ETOL
C             for each vector component.
C
C      MORDER -- This is the order of the formula which will be used by
C             the initial value method for taking the first integration
C             step.
C
C      SMALL -- This is a small positive machine dependent constant
C             which is used for protecting against computations with
C             numbers which are too small relative to the precision of
C             floating point arithmetic.  SMALL  should be set to
C             (approximately) the smallest positive DOUBLE PRECISION
C             number such that  (1.+SMALL) .GT. 1.  on the machine being
C             used. The quantity  SMALL**(3/8)  is used in computing
C             increments of variables for approximating derivatives by
C             differences.  Also the algorithm will not compute a
C             starting step length which is smaller than
C             100*SMALL*ABS(A).
C
C      BIG -- This is a large positive machine dependent constant which
C             is used for preventing machine overflows. A reasonable
C             choice is to set big to (approximately) the square root of
C             the largest DOUBLE PRECISION number which can be held in
C             the machine.
C
C      SPY(*),PV(*),YP(*),SF(*) -- These are DOUBLE PRECISION work
C             arrays of length NEQ which provide the routine with needed
C             storage space.
C
C      RPAR,IPAR -- These are parameter arrays, of DOUBLE PRECISION and
C             INTEGER type, respectively, which can be used for
C             communication between your program and the DF subroutine.
C             They are not used or altered by DHSTRT.
C
C **********************************************************************
C  On Output  (after the return from DHSTRT),
C
C      H -- is an appropriate starting step size to be attempted by the
C             differential equation method.
C
C           All parameters in the call list remain unchanged except for
C           the working arrays SPY(*),PV(*),YP(*), and SF(*).
C
C **********************************************************************
C
C***SEE ALSO  DDEABM, DDEBDF, DDERKF
C***ROUTINES CALLED  DHVNRM
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891024  Changed references from DVNORM to DHVNRM.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DHSTRT
C
      INTEGER IPAR, J, K, LK, MORDER, NEQ
      DOUBLE PRECISION A, ABSDX, B, BIG, DA, DELF, DELY,
     1      DFDUB, DFDXB, DHVNRM,
     2      DX, DY, ETOL, FBND, H, PV, RELPER, RPAR, SF, SMALL, SPY,
     3      SRYDPB, TOLEXP, TOLMIN, TOLP, TOLSUM, Y, YDPB, YP, YPRIME
      DIMENSION Y(*),YPRIME(*),ETOL(*),SPY(*),PV(*),YP(*),
     1          SF(*),RPAR(*),IPAR(*)
      EXTERNAL DF
C
C     ..................................................................
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 160
C***FIRST EXECUTABLE STATEMENT  DHSTRT
         DX = B - A
         ABSDX = ABS(DX)
         RELPER = SMALL**0.375D0
C
C        ...............................................................
C
C             COMPUTE AN APPROXIMATE BOUND (DFDXB) ON THE PARTIAL
C             DERIVATIVE OF THE EQUATION WITH RESPECT TO THE
C             INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW.
C             ALSO COMPUTE A BOUND (FBND) ON THE FIRST DERIVATIVE
C             LOCALLY.
C
         DA = SIGN(MAX(MIN(RELPER*ABS(A),ABSDX),
     1                    100.0D0*SMALL*ABS(A)),DX)
         IF (DA .EQ. 0.0D0) DA = RELPER*DX
         CALL DF(A+DA,Y,SF,RPAR,IPAR)
         DO 10 J = 1, NEQ
            YP(J) = SF(J) - YPRIME(J)
   10    CONTINUE
         DELF = DHVNRM(YP,NEQ)
         DFDXB = BIG
         IF (DELF .LT. BIG*ABS(DA)) DFDXB = DELF/ABS(DA)
         FBND = DHVNRM(SF,NEQ)
C
C        ...............................................................
C
C             COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ
C             CONSTANT FOR THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS
C             ALSO REPRESENTS AN ESTIMATE OF THE NORM OF THE JACOBIAN
C             LOCALLY.  THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO
C             ESTIMATE THE LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES.
C             THE FIRST PERTURBATION VECTOR IS BASED ON THE INITIAL
C             DERIVATIVES AND DIRECTION OF INTEGRATION. THE SECOND
C             PERTURBATION VECTOR IS FORMED USING ANOTHER EVALUATION OF
C             THE DIFFERENTIAL EQUATION.  THE THIRD PERTURBATION VECTOR
C             IS FORMED USING PERTURBATIONS BASED ONLY ON THE INITIAL
C             VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS CHANGED TO
C             NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN
C             INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT
C             COMPONENTS OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE
C             CONSISTENT WITH THE SLOPES OF LOCAL SOLUTION CURVES.
C             ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST
C             DERIVATIVE.
C
C                               PERTURBATION VECTOR SIZE IS HELD
C                               CONSTANT FOR ALL ITERATIONS. COMPUTE
C                               THIS CHANGE FROM THE
C                                       SIZE OF THE VECTOR OF INITIAL
C                                       VALUES.
         DELY = RELPER*DHVNRM(Y,NEQ)
         IF (DELY .EQ. 0.0D0) DELY = RELPER
         DELY = SIGN(DELY,DX)
         DELF = DHVNRM(YPRIME,NEQ)
         FBND = MAX(FBND,DELF)
         IF (DELF .EQ. 0.0D0) GO TO 30
C           USE INITIAL DERIVATIVES FOR FIRST PERTURBATION
            DO 20 J = 1, NEQ
               SPY(J) = YPRIME(J)
               YP(J) = YPRIME(J)
   20       CONTINUE
         GO TO 50
   30    CONTINUE
C           CANNOT HAVE A NULL PERTURBATION VECTOR
            DO 40 J = 1, NEQ
               SPY(J) = 0.0D0
               YP(J) = 1.0D0
   40       CONTINUE
            DELF = DHVNRM(YP,NEQ)
   50    CONTINUE
C
         DFDUB = 0.0D0
         LK = MIN(NEQ+1,3)
         DO 140 K = 1, LK
C           DEFINE PERTURBED VECTOR OF INITIAL VALUES
            DO 60 J = 1, NEQ
               PV(J) = Y(J) + DELY*(YP(J)/DELF)
   60       CONTINUE
            IF (K .EQ. 2) GO TO 80
C              EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED
C              VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES
               CALL DF(A,PV,YP,RPAR,IPAR)
               DO 70 J = 1, NEQ
                  PV(J) = YP(J) - YPRIME(J)
   70          CONTINUE
            GO TO 100
   80       CONTINUE
C              USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE
C                                    IN COMPUTING ONE ESTIMATE
               CALL DF(A+DA,PV,YP,RPAR,IPAR)
               DO 90 J = 1, NEQ
                  PV(J) = YP(J) - SF(J)
   90          CONTINUE
  100       CONTINUE
C           CHOOSE LARGEST BOUNDS ON THE FIRST DERIVATIVE
C                          AND A LOCAL LIPSCHITZ CONSTANT
            FBND = MAX(FBND,DHVNRM(YP,NEQ))
            DELF = DHVNRM(PV,NEQ)
C        ...EXIT
            IF (DELF .GE. BIG*ABS(DELY)) GO TO 150
            DFDUB = MAX(DFDUB,DELF/ABS(DELY))
C     ......EXIT
            IF (K .EQ. LK) GO TO 160
C           CHOOSE NEXT PERTURBATION VECTOR
            IF (DELF .EQ. 0.0D0) DELF = 1.0D0
            DO 130 J = 1, NEQ
               IF (K .EQ. 2) GO TO 110
                  DY = ABS(PV(J))
                  IF (DY .EQ. 0.0D0) DY = DELF
               GO TO 120
  110          CONTINUE
                  DY = Y(J)
                  IF (DY .EQ. 0.0D0) DY = DELY/RELPER
  120          CONTINUE
               IF (SPY(J) .EQ. 0.0D0) SPY(J) = YP(J)
               IF (SPY(J) .NE. 0.0D0) DY = SIGN(DY,SPY(J))
               YP(J) = DY
  130       CONTINUE
            DELF = DHVNRM(YP,NEQ)
  140    CONTINUE
  150    CONTINUE
C
C        PROTECT AGAINST AN OVERFLOW
         DFDUB = BIG
  160 CONTINUE
C
C     ..................................................................
C
C          COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE
C
      YDPB = DFDXB + DFDUB*FBND
C
C     ..................................................................
C
C          DEFINE THE TOLERANCE PARAMETER UPON WHICH THE STARTING STEP
C          SIZE IS TO BE BASED.  A VALUE IN THE MIDDLE OF THE ERROR
C          TOLERANCE RANGE IS SELECTED.
C
      TOLMIN = BIG
      TOLSUM = 0.0D0
      DO 170 K = 1, NEQ
         TOLEXP = LOG10(ETOL(K))
         TOLMIN = MIN(TOLMIN,TOLEXP)
         TOLSUM = TOLSUM + TOLEXP
  170 CONTINUE
      TOLP = 10.0D0**(0.5D0*(TOLSUM/NEQ + TOLMIN)/(MORDER+1))
C
C     ..................................................................
C
C          COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND
C          SECOND DERIVATIVE INFORMATION
C
C                            RESTRICT THE STEP LENGTH TO BE NOT BIGGER
C                            THAN ABS(B-A).   (UNLESS  B  IS TOO CLOSE
C                            TO  A)
      H = ABSDX
C
      IF (YDPB .NE. 0.0D0 .OR. FBND .NE. 0.0D0) GO TO 180
C
C        BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND
C                     DERIVATIVE TERM (YDPB) ARE ZERO
         IF (TOLP .LT. 1.0D0) H = ABSDX*TOLP
      GO TO 200
  180 CONTINUE
C
      IF (YDPB .NE. 0.0D0) GO TO 190
C
C        ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO
         IF (TOLP .LT. FBND*ABSDX) H = TOLP/FBND
      GO TO 200
  190 CONTINUE
C
C        SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO
         SRYDPB = SQRT(0.5D0*YDPB)
         IF (TOLP .LT. SRYDPB*ABSDX) H = TOLP/SRYDPB
  200 CONTINUE
C
C     FURTHER RESTRICT THE STEP LENGTH TO BE NOT
C                               BIGGER THAN  1/DFDUB
      IF (H*DFDUB .GT. 1.0D0) H = 1.0D0/DFDUB
C
C     FINALLY, RESTRICT THE STEP LENGTH TO BE NOT
C     SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF
C     A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO,
C     THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE
C                                     STEP LENGTH.
      H = MAX(H,100.0D0*SMALL*ABS(A))
      IF (H .EQ. 0.0D0) H = SMALL*ABS(B)
C
C     NOW SET DIRECTION OF INTEGRATION
      H = SIGN(H,DX)
C
      RETURN
      END
*DECK DHVNRM
      DOUBLE PRECISION FUNCTION DHVNRM (V, NCOMP)
C***BEGIN PROLOGUE  DHVNRM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (HVNRM-S, DHVNRM-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     Compute the maximum norm of the vector V(*) of length NCOMP and
C     return the result as DHVNRM
C
C***SEE ALSO  DDEABM, DDEBDF, DDERKF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891024  Changed references from DVNORM to DHVNRM.  (WRB)
C   891024  Changed routine name from DVNORM to DHVNRM.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DHVNRM
C
      INTEGER K, NCOMP
      DOUBLE PRECISION V
      DIMENSION V(*)
C***FIRST EXECUTABLE STATEMENT  DHVNRM
      DHVNRM = 0.0D0
      DO 10 K = 1, NCOMP
         DHVNRM = MAX(DHVNRM,ABS(V(K)))
   10 CONTINUE
      RETURN
      END
*DECK J4SAVE
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
C***BEGIN PROLOGUE  J4SAVE
C***SUBSIDIARY
C***PURPOSE  Save or recall global variables needed by error
C            handling routines.
C***LIBRARY   SLATEC (XERROR)
C***TYPE      INTEGER (J4SAVE-I)
C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                = 6 Refers to the 2nd unit for error messages
C                = 7 Refers to the 3rd unit for error messages
C                = 8 Refers to the 4th unit for error messages
C                = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C***SEE ALSO  XERMSG
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900205  Minor modifications to prologue.  (WRB)
C   900402  Added TYPE section.  (WRB)
C   910411  Added KEYWORDS section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
*DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
C***BEGIN PROLOGUE  XERCNT
C***SUBSIDIARY
C***PURPOSE  Allow user control over handling of errors.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCNT-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCNT.
C        If the user has provided his own version of XERCNT, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        LIBRAR - the library that the routine is in.
C        SUBROU - the subroutine that XERMSG is being called from
C        MESSG  - the first 20 characters of the error message.
C        NERR   - same as in the call to XERMSG.
C        LEVEL  - same as in the call to XERMSG.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
C           names, changed routine name from XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
*DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
C***BEGIN PROLOGUE  XERHLT
C***SUBSIDIARY
C***PURPOSE  Abort program execution and print error message.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERHLT-A)
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERHLT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG is as in XERMSG.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to delete length of character
C           and changed routine name from XERABT to XERHLT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END
