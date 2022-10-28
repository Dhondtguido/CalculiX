!
!    SLATEC: public domain
!
*DECK DDEBDF
      SUBROUTINE DDEBDF (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,
     +   RWORK, LRW, IWORK, LIW, RPAR, IPAR, DJAC)
C***BEGIN PROLOGUE  DDEBDF
C***PURPOSE  Solve an initial value problem in ordinary differential
C            equations using backward differentiation formulas.  It is
C            intended primarily for stiff problems.
C***LIBRARY   SLATEC (DEPAC)
C***CATEGORY  I1A2
C***TYPE      DOUBLE PRECISION (DEBDF-S, DDEBDF-D)
C***KEYWORDS  BACKWARD DIFFERENTIATION FORMULAS, DEPAC,
C             INITIAL VALUE PROBLEMS, ODE,
C             ORDINARY DIFFERENTIAL EQUATIONS, STIFF
C***AUTHOR  Shampine, L. F., (SNLA)
C           Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   This is the backward differentiation code in the package of
C   differential equation solvers DEPAC, consisting of the codes
C   DDERKF, DDEABM, and DDEBDF.  Design of the package was by
C   L. F. Shampine and H. A. Watts.  It is documented in
C        SAND-79-2374 , DEPAC - Design of a User Oriented Package of ODE
C                              Solvers.
C   DDEBDF is a driver for a modification of the code LSODE written by
C             A. C. Hindmarsh
C             Lawrence Livermore Laboratory
C             Livermore, California 94550
C
C **********************************************************************
C **             DEPAC PACKAGE OVERVIEW           **
C **********************************************************************
C
C        You have a choice of three differential equation solvers from
C        DEPAC.  The following brief descriptions are meant to aid you
C        in choosing the most appropriate code for your problem.
C
C        DDERKF is a fifth order Runge-Kutta code. It is the simplest of
C        the three choices, both algorithmically and in the use of the
C        code. DDERKF is primarily designed to solve non-stiff and mild-
C        ly stiff differential equations when derivative evaluations are
C        not expensive.  It should generally not be used to get high
C        accuracy results nor answers at a great many specific points.
C        Because DDERKF has very low overhead costs, it will usually
C        result in the least expensive integration when solving
C        problems requiring a modest amount of accuracy and having
C        equations that are not costly to evaluate.  DDERKF attempts to
C        discover when it is not suitable for the task posed.
C
C        DDEABM is a variable order (one through twelve) Adams code. Its
C        complexity lies somewhere between that of DDERKF and DDEBDF.
C        DDEABM is primarily designed to solve non-stiff and mildly
C        stiff differential equations when derivative evaluations are
C        expensive, high accuracy results are needed or answers at
C        many specific points are required.  DDEABM attempts to discover
C        when it is not suitable for the task posed.
C
C        DDEBDF is a variable order (one through five) backward
C        differentiation formula code.  It is the most complicated of
C        the three choices.  DDEBDF is primarily designed to solve stiff
C        differential equations at crude to moderate tolerances.
C        If the problem is very stiff at all, DDERKF and DDEABM will be
C        quite inefficient compared to DDEBDF.  However, DDEBDF will be
C        inefficient compared to DDERKF and DDEABM on non-stiff problems
C        because it uses much more storage, has a much larger overhead,
C        and the low order formulas will not give high accuracies
C        efficiently.
C
C        The concept of stiffness cannot be described in a few words.
C        If you do not know the problem to be stiff, try either DDERKF
C        or DDEABM.  Both of these codes will inform you of stiffness
C        when the cost of solving such problems becomes important.
C
C **********************************************************************
C ** ABSTRACT **
C **********************************************************************
C
C   Subroutine DDEBDF uses the backward differentiation formulas of
C   orders one through five to integrate a system of NEQ first order
C   ordinary differential equations of the form
C                         DU/DX = DF(X,U)
C   when the vector Y(*) of initial values for U(*) at X=T is given.
C   The subroutine integrates from T to TOUT. It is easy to continue the
C   integration to get results at additional TOUT. This is the interval
C   mode of operation. It is also easy for the routine to return with
C   the solution at each intermediate step on the way to TOUT. This is
C   the intermediate-output mode of operation.
C
C **********************************************************************
C * Description of The Arguments To DDEBDF (An Overview) *
C **********************************************************************
C
C   The Parameters are:
C
C      DF -- This is the name of a subroutine which you provide to
C            define the differential equations.
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
C      RTOL, ATOL -- These DOUBLE PRECISION quantities
C             represent relative and absolute error tolerances which you
C             provide to indicate how accurately you wish the solution
C             to be computed.  You may choose them to be both scalars
C             or else both vectors.
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
C             calling program and the DF subroutine (and the DJAC
C             subroutine).
C
C      DJAC -- This is the name of a subroutine which you may choose to
C             provide for defining the Jacobian matrix of partial
C             derivatives DF/DU.
C
C  Quantities which are used as input items are
C             NEQ, T, Y(*), TOUT, INFO(*),
C             RTOL, ATOL, RWORK(1), LRW,
C             IWORK(1), IWORK(2), and LIW.
C
C  Quantities which may be altered by the code are
C             T, Y(*), INFO(1), RTOL, ATOL,
C             IDID, RWORK(*) and IWORK(*).
C
C **********************************************************************
C * INPUT -- What To Do On The First Call To DDEBDF *
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
C             which is to be solved. For the given values of X and the
C             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
C             evaluate the NEQ components of the system of differential
C             equations  DU/DX=DF(X,U)  and store the derivatives in the
C             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
C             equations I=1,...,NEQ.
C
C             Subroutine DF must not alter X or U(*). You must declare
C             the name DF in an external statement in your program that
C             calls DDEBDF. You must dimension U and UPRIME in DF.
C
C             RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
C             arrays which you can use for communication between your
C             calling program and subroutine DF. They are not used or
C             altered by DDEBDF.  If you do not need RPAR or IPAR,
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
C      TOUT -- Set it to the first point at which a solution is desired.
C             You can take TOUT = T, in which case the code
C             will evaluate the derivative of the solution at T and
C             return.  Integration either forward in T  (TOUT .GT. T)
C             or backward in T  (TOUT .LT. T)  is permitted.
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
C             the problem.  By using the fact that the code will not
C             step past TOUT in the first step, you could, if necessary,
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
C             DEPAC or possible future extensions, though DDEBDF uses
C             only the first six entries.  You must respond to all of
C             the following items which are arranged as questions.  The
C             simplest use of the code corresponds to answering all
C             questions as YES ,i.e. setting all entries of INFO to 0.
C
C        INFO(1) -- This parameter enables the code to initialize
C               itself.  You must set it to indicate the start of every
C               new problem.
C
C            **** Is this the first call for this problem ...
C                  YES -- Set INFO(1) = 0
C                   NO -- Not applicable here.
C                         See below for continuation calls.  ****
C
C        INFO(2) -- How much accuracy you want of your solution
C               is specified by the error tolerances RTOL and ATOL.
C               The simplest use is to take them both to be scalars.
C               To obtain more flexibility, they can both be vectors.
C               The code must be told your choice.
C
C            **** Are both error tolerances RTOL, ATOL scalars ...
C                  YES -- Set INFO(2) = 0
C                         and input scalars for both RTOL and ATOL
C                   NO -- Set INFO(2) = 1
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
C                 TOUT (and NOT at the next intermediate step) ...
C                  YES -- Set INFO(3) = 0
C                   NO -- Set INFO(3) = 1 ****
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
C                 restrictions on the independent variable T ...
C                  YES -- Set INFO(4)=0
C                   NO -- Set INFO(4)=1
C                         and define the stopping point TSTOP by
C                         setting RWORK(1)=TSTOP ****
C
C        INFO(5) -- To solve stiff problems it is necessary to use the
C               Jacobian matrix of partial derivatives of the system
C               of differential equations.  If you do not provide a
C               subroutine to evaluate it analytically (see the
C               description of the item DJAC in the call list), it will
C               be approximated by numerical differencing in this code.
C               Although it is less trouble for you to have the code
C               compute partial derivatives by numerical differencing,
C               the solution will be more reliable if you provide the
C               derivatives via DJAC.  Sometimes numerical differencing
C               is cheaper than evaluating derivatives in DJAC and
C               sometimes it is not - this depends on your problem.
C
C               If your problem is linear, i.e. has the form
C               DU/DX = DF(X,U) = J(X)*U + G(X)   for some matrix J(X)
C               and vector G(X), the Jacobian matrix  DF/DU = J(X).
C               Since you must provide a subroutine to evaluate DF(X,U)
C               analytically, it is little extra trouble to provide
C               subroutine DJAC for evaluating J(X) analytically.
C               Furthermore, in such cases, numerical differencing is
C               much more expensive than analytic evaluation.
C
C            **** Do you want the code to evaluate the partial
C                 derivatives automatically by numerical differences ...
C                  YES -- Set INFO(5)=0
C                   NO -- Set INFO(5)=1
C                         and provide subroutine DJAC for evaluating the
C                         Jacobian matrix ****
C
C        INFO(6) -- DDEBDF will perform much better if the Jacobian
C               matrix is banded and the code is told this.  In this
C               case, the storage needed will be greatly reduced,
C               numerical differencing will be performed more cheaply,
C               and a number of important algorithms will execute much
C               faster.  The differential equation is said to have
C               half-bandwidths ML (lower) and MU (upper) if equation I
C               involves only unknowns Y(J) with
C                              I-ML .LE. J .LE. I+MU
C               for all I=1,2,...,NEQ.  Thus, ML and MU are the widths
C               of the lower and upper parts of the band, respectively,
C               with the main diagonal being excluded.  If you do not
C               indicate that the equation has a banded Jacobian,
C               the code works with a full matrix of NEQ**2 elements
C               (stored in the conventional way).  Computations with
C               banded matrices cost less time and storage than with
C               full matrices if  2*ML+MU .LT. NEQ.  If you tell the
C               code that the Jacobian matrix has a banded structure and
C               you want to provide subroutine DJAC to compute the
C               partial derivatives, then you must be careful to store
C               the elements of the Jacobian matrix in the special form
C               indicated in the description of DJAC.
C
C            **** Do you want to solve the problem using a full
C                 (dense) Jacobian matrix (and not a special banded
C                 structure) ...
C                  YES -- Set INFO(6)=0
C                   NO -- Set INFO(6)=1
C                         and provide the lower (ML) and upper (MU)
C                         bandwidths by setting
C                         IWORK(1)=ML
C                         IWORK(2)=MU ****
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
C             (More specifically, a root-mean-square norm is used to
C             measure the size of vectors, and the error test uses the
C             magnitude of the solution at the beginning of the step.)
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
C             Setting ATOL=0. results in a pure relative error test on
C             that component.  Setting RTOL=0. results in a pure abso-
C             lute error test on that component.  A mixed test with non-
C             zero RTOL and ATOL corresponds roughly to a relative error
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
C             (For some problems it may not be permissible to integrate
C             past a point TSTOP because a discontinuity occurs there
C             or the solution or its derivative is not defined beyond
C             TSTOP.)
C
C      LRW -- Set it to the declared length of the RWORK array.
C             You must have
C                  LRW .GE. 250+10*NEQ+NEQ**2
C             for the full (dense) Jacobian case (when INFO(6)=0),  or
C                  LRW .GE. 250+10*NEQ+(2*ML+MU+1)*NEQ
C             for the banded Jacobian case (when INFO(6)=1).
C
C      IWORK(*) -- Dimension this INTEGER work array of length LIW in
C             your calling program.
C
C      IWORK(1), IWORK(2) -- If you have set INFO(6)=0, you can ignore
C             these optional input parameters. Otherwise you must define
C             the half-bandwidths ML (lower) and MU (upper) of the
C             Jacobian matrix by setting    IWORK(1) = ML   and
C             IWORK(2) = MU.  (The code will work with a full matrix
C             of NEQ**2 elements unless it is told that the problem has
C             a banded Jacobian, in which case the code will work with
C             a matrix containing at most  (2*ML+MU+1)*NEQ  elements.)
C
C      LIW -- Set it to the declared length of the IWORK array.
C             You must have LIW .GE. 56+NEQ.
C
C      RPAR, IPAR -- These are parameter arrays, of DOUBLE PRECISION and
C             INTEGER type, respectively. You can use them for
C             communication between your program that calls DDEBDF and
C             the  DF subroutine (and the DJAC subroutine). They are not
C             used or altered by DDEBDF. If you do not need RPAR or
C             IPAR, ignore these parameters by treating them as dummy
C             arguments. If you do choose to use them, dimension them in
C             your calling program and in DF (and in DJAC) as arrays of
C             appropriate length.
C
C      DJAC -- If you have set INFO(5)=0, you can ignore this parameter
C             by treating it as a dummy argument. (For some compilers
C             you may have to write a dummy subroutine named  DJAC  in
C             order to avoid problems associated with missing external
C             routine names.)  Otherwise, you must provide a subroutine
C             of the form
C                          DJAC(X,U,PD,NROWPD,RPAR,IPAR)
C             to define the Jacobian matrix of partial derivatives DF/DU
C             of the system of differential equations   DU/DX = DF(X,U).
C             For the given values of X and the vector
C             U(*)=(U(1),U(2),...,U(NEQ)), the subroutine must evaluate
C             the non-zero partial derivatives  DF(I)/DU(J)  for each
C             differential equation I=1,...,NEQ and each solution
C             component J=1,...,NEQ , and store these values in the
C             matrix PD.  The elements of PD are set to zero before each
C             call to DJAC so only non-zero elements need to be defined.
C
C             Subroutine DJAC must not alter X, U(*), or NROWPD. You
C             must declare the name DJAC in an external statement in
C             your program that calls DDEBDF. NROWPD is the row
C             dimension of the PD matrix and is assigned by the code.
C             Therefore you must dimension PD in DJAC according to
C                              DIMENSION PD(NROWPD,1)
C             You must also dimension U in DJAC.
C
C             The way you must store the elements into the PD matrix
C             depends on the structure of the Jacobian which you
C             indicated by INFO(6).
C             *** INFO(6)=0 -- Full (Dense) Jacobian ***
C                 When you evaluate the (non-zero) partial derivative
C                 of equation I with respect to variable J, you must
C                 store it in PD according to
C                                PD(I,J) = * DF(I)/DU(J) *
C             *** INFO(6)=1 -- Banded Jacobian with ML Lower and MU
C                 Upper Diagonal Bands (refer to INFO(6) description of
C                 ML and MU) ***
C                 When you evaluate the (non-zero) partial derivative
C                 of equation I with respect to variable J, you must
C                 store it in PD according to
C                                IROW = I - J + ML + MU + 1
C                                PD(IROW,J) = * DF(I)/DU(J) *
C
C             RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
C             arrays which you can use for communication between your
C             calling program and your Jacobian subroutine DJAC. They
C             are not altered by DDEBDF. If you do not need RPAR or
C             IPAR, ignore these parameters by treating them as dummy
C             arguments.  If you do choose to use them, dimension them
C             in your calling program and in DJAC as arrays of
C             appropriate length.
C
C **********************************************************************
C * OUTPUT -- After any return from DDEBDF *
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
C             IDID = -4,-5  -- Not applicable for this code but used
C                       by other members of DEPAC.
C
C             IDID = -6 -- DDEBDF had repeated convergence test failures
C                       on the last attempted step.
C
C             IDID = -7 -- DDEBDF had repeated error test failures on
C                       the last attempted step.
C
C             IDID = -8,..,-32  -- Not applicable for this code but
C                       used by other members of DEPAC or possible
C                       future extensions.
C
C                         *** Task Terminated ***
C                   Reported by the value of IDID=-33
C
C             IDID = -33 -- The code has encountered trouble from which
C                       it cannot recover.  A message is printed
C                       explaining the trouble and control is returned
C                       to the calling program.  For example, this
C                       occurs when invalid input is detected.
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
C             RWORK(12)--If the tolerances have been increased by the
C                        code (IDID = -2) , they were multiplied by the
C                        value in RWORK(12).
C
C             RWORK(13)--which contains the current value of the
C                        independent variable, i.e. the farthest point
C                        integration has reached.  This will be
C                        different from T only when interpolation has
C                        been performed (IDID=3).
C
C             RWORK(20+I)--which contains the approximate derivative
C                        of the solution component Y(I).  In DDEBDF, it
C                        is never obtained by calling subroutine DF to
C                        evaluate the differential equation using T and
C                        Y(*), except at the initial point of
C                        integration.
C
C **********************************************************************
C ** INPUT -- What To Do To Continue The Integration **
C **             (calls after the first)             **
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
C        Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2)
C        unless you are going to restart the code.
C
C        The parameter INFO(1) is used by the code to indicate the
C        beginning of a new problem and to indicate whether integration
C        is to be continued.  You must input the value  INFO(1) = 0
C        when starting a new problem.  You must input the value
C        INFO(1) = 1  if you wish to continue after an interrupted task.
C        Do not set  INFO(1) = 0  on a continuation call unless you
C        want the code to restart at the current T.
C
C                         *** Following a Completed Task ***
C         If
C             IDID = 1, call the code again to continue the integration
C                     another step in the direction of TOUT.
C
C             IDID = 2 or 3, define a new TOUT and call the code again.
C                     TOUT must be different from T.  You cannot change
C                     the direction of integration without restarting.
C
C                         *** Following an Interrupted Task ***
C                     To show the code that you realize the task was
C                     interrupted and that you want to continue, you
C                     must take appropriate action and reset INFO(1) = 1
C         If
C             IDID = -1, the code has attempted 500 steps.
C                     If you want to continue, set INFO(1) = 1 and
C                     call the code again.  An additional 500 steps
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
C             IDID = -4,-5  --- cannot occur with this code but used
C                     by other members of DEPAC.
C
C             IDID = -6, repeated convergence test failures occurred
C                     on the last attempted step in DDEBDF.  An inaccu-
C                     rate Jacobian may be the problem.  If you are
C                     absolutely certain you want to continue, restart
C                     the integration at the current T by setting
C                     INFO(1)=0 and call the code again.
C
C             IDID = -7, repeated error test failures occurred on the
C                     last attempted step in DDEBDF.  A singularity in
C                     the solution may be present.  You should re-
C                     examine the problem being solved.  If you are
C                     absolutely certain you want to continue, restart
C                     the integration at the current T by setting
C                     INFO(1)=0 and call the code again.
C
C             IDID = -8,..,-32  --- cannot occur with this code but
C                     used by other members of DDEPAC or possible future
C                     extensions.
C
C                         *** Following a Terminated Task ***
C         If
C             IDID = -33, you cannot continue the solution of this
C                     problem.  An attempt to do so will result in your
C                     run being terminated.
C
C **********************************************************************
C
C         ***** Warning *****
C
C     If DDEBDF is to be used in an overlay situation, you must save and
C     restore certain items used internally by DDEBDF  (values in the
C     common block DDEBD1).  This can be accomplished as follows.
C
C     To save the necessary values upon return from DDEBDF, simply call
C        DSVCO(RWORK(22+NEQ),IWORK(21+NEQ)).
C
C     To restore the necessary values before the next call to DDEBDF,
C     simply call    DRSCO(RWORK(22+NEQ),IWORK(21+NEQ)).
C
C***REFERENCES  L. F. Shampine and H. A. Watts, DEPAC - design of a user
C                 oriented package of ODE solvers, Report SAND79-2374,
C                 Sandia Laboratories, 1979.
C***ROUTINES CALLED  DLSOD, XERMSG
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891024  Changed references from DVNORM to DHVNRM.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900510  Convert XERRWV calls to XERMSG calls, make Prologue comments
C           consistent with DEBDF.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDEBDF
      INTEGER IACOR, IBAND, IBEGIN, ICOMI, ICOMR, IDELSN, IDID, IER,
     1      IEWT, IINOUT, IINTEG, IJAC, ILRW, INFO, INIT,
     2      IOWNS, IPAR, IQUIT, ISAVF, ITOL, ITSTAR, ITSTOP, IWM,
     3      IWORK, IYH, IYPOUT, JSTART, KFLAG, KSTEPS, L, LIW, LRW,
     4      MAXORD, METH, MITER, ML, MU, N, NEQ, NFE, NJE, NQ, NQU,
     5      NST
      DOUBLE PRECISION ATOL, EL0, H, HMIN, HMXI, HU, ROWNS, RPAR,
     1      RTOL, RWORK, T, TN, TOLD, TOUT, UROUND, Y
      LOGICAL INTOUT
      CHARACTER*8 XERN1, XERN2
      CHARACTER*16 XERN3
C
      DIMENSION Y(*),INFO(15),RTOL(*),ATOL(*),RWORK(*),IWORK(*),
     1          RPAR(*),IPAR(*)
C
      COMMON /DDEBD1/ TOLD,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,
     1                IQUIT,INIT,IYH,IEWT,IACOR,ISAVF,IWM,KSTEPS,IBEGIN,
     2                ITOL,IINTEG,ITSTOP,IJAC,IBAND,IOWNS(6),IER,JSTART,
     3                KFLAG,L,METH,MITER,MAXORD,N,NQ,NST,NFE,NJE,NQU
C
      EXTERNAL DF, DJAC
C
C        CHECK FOR AN APPARENT INFINITE LOOP
C
C***FIRST EXECUTABLE STATEMENT  DDEBDF
      IF (INFO(1) .EQ. 0) IWORK(LIW) = 0
C
      IF (IWORK(LIW).GE. 5) THEN
         IF (T .EQ. RWORK(21+NEQ)) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDEBDF',
     *         'AN APPARENT INFINITE LOOP HAS BEEN DETECTED.$$' //
     *         'YOU HAVE MADE REPEATED CALLS AT T = ' // XERN3 //
     *         ' AND THE INTEGRATION HAS NOT ADVANCED.  CHECK THE ' //
     *         'WAY YOU HAVE SET PARAMETERS FOR THE CALL TO THE ' //
     *         'CODE, PARTICULARLY INFO(1).', 13, 2)
            RETURN
         ENDIF
      ENDIF
C
      IDID = 0
C
C        CHECK VALIDITY OF INFO PARAMETERS
C
      IF (INFO(1) .NE. 0 .AND. INFO(1) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(1)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(1) MUST BE SET TO 0 ' //
     *      'FOR THE  START OF A NEW PROBLEM, AND MUST BE SET TO 1 ' //
     *      'FOLLOWING AN INTERRUPTED TASK.  YOU ARE ATTEMPTING TO ' //
     *      'CONTINUE THE INTEGRATION ILLEGALLY BY CALLING THE ' //
     *      'CODE WITH  INFO(1) = ' // XERN1, 3, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(2) .NE. 0 .AND. INFO(2) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(2)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(2) MUST BE 0 OR 1 ' //
     *      'INDICATING SCALAR AND VECTOR ERROR TOLERANCES, ' //
     *      'RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH INFO(2) = ' //
     *      XERN1, 4, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(3) .NE. 0 .AND. INFO(3) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(3)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(3) MUST BE 0 OR 1 ' //
     *      'INDICATING THE INTERVAL OR INTERMEDIATE-OUTPUT MODE OF ' //
     *      'INTEGRATION, RESPECTIVELY.  YOU HAVE CALLED THE CODE ' //
     *      'WITH  INFO(3) = ' // XERN1, 5, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(4) .NE. 0 .AND. INFO(4) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(4)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(4) MUST BE 0 OR 1 ' //
     *      'INDICATING WHETHER OR NOT THE INTEGRATION INTERVAL IS ' //
     *      'TO BE RESTRICTED BY A POINT TSTOP.  YOU HAVE CALLED ' //
     *      'THE CODE WITH INFO(4) = ' // XERN1, 14, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(5) .NE. 0 .AND. INFO(5) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(5)
         CALL XERMSG ('SLATEC',  'DDEBDF', 'INFO(5) MUST BE 0 OR 1 ' //
     *      'INDICATING WHETHER THE CODE IS TOLD TO FORM THE ' //
     *      'JACOBIAN MATRIX BY NUMERICAL DIFFERENCING OR YOU ' //
     *      'PROVIDE A SUBROUTINE TO EVALUATE IT ANALYTICALLY.  ' //
     *      'YOU HAVE CALLED THE CODE WITH INFO(5) = ' // XERN1, 15, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(6) .NE. 0 .AND. INFO(6) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(6)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(6) MUST BE 0 OR 1 ' //
     *      'INDICATING WHETHER THE CODE IS TOLD TO TREAT THE ' //
     *      'JACOBIAN AS A FULL (DENSE) MATRIX OR AS HAVING A ' //
     *      'SPECIAL BANDED STRUCTURE.  YOU HAVE CALLED THE CODE ' //
     *      'WITH INFO(6) = ' // XERN1, 16, 1)
         IDID = -33
      ENDIF
C
      ILRW = NEQ
      IF (INFO(6) .NE. 0) THEN
C
C        CHECK BANDWIDTH PARAMETERS
C
         ML = IWORK(1)
         MU = IWORK(2)
         ILRW = 2*ML + MU + 1
C
         IF (ML.LT.0 .OR. ML.GE.NEQ .OR. MU.LT.0 .OR. MU.GE.NEQ) THEN
            WRITE (XERN1, '(I8)') ML
            WRITE (XERN2, '(I8)') MU
            CALL XERMSG ('SLATEC', 'DDEBDF', 'YOU HAVE SET INFO(6) ' //
     *         '= 1, TELLING THE CODE THAT THE JACOBIAN MATRIX HAS ' //
     *         'A SPECIAL BANDED STRUCTURE.  HOWEVER, THE LOWER ' //
     *         '(UPPER) BANDWIDTHS  ML (MU) VIOLATE THE CONSTRAINTS ' //
     *         'ML,MU .GE. 0 AND  ML,MU .LT. NEQ.  YOU HAVE CALLED ' //
     *         'THE CODE WITH ML = ' // XERN1 // ' AND MU = ' // XERN2,
     *         17, 1)
            IDID = -33
         ENDIF
      ENDIF
C
C        CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION
C
      IF (LRW .LT. 250 + (10 + ILRW)*NEQ) THEN
         WRITE (XERN1, '(I8)') LRW
         IF (INFO(6) .EQ. 0) THEN
            CALL XERMSG ('SLATEC', 'DDEBDF', 'LENGTH OF ARRAY RWORK ' //
     *         'MUST BE AT LEAST 250 + 10*NEQ + NEQ*NEQ.$$' //
     *         'YOU HAVE CALLED THE CODE WITH  LRW = ' // XERN1, 1, 1)
         ELSE
            CALL XERMSG ('SLATEC', 'DDEBDF', 'LENGTH OF ARRAY RWORK ' //
     *         'MUST BE AT LEAST 250 + 10*NEQ + (2*ML+MU+1)*NEQ.$$' //
     *         'YOU HAVE CALLED THE CODE WITH  LRW = ' // XERN1, 18, 1)
         ENDIF
         IDID = -33
      ENDIF
C
      IF (LIW .LT. 56 + NEQ) THEN
         WRITE (XERN1, '(I8)') LIW
         CALL XERMSG ('SLATEC', 'DDEBDF', 'LENGTH OF ARRAY IWORK ' //
     *      'BE AT LEAST  56 + NEQ.  YOU HAVE CALLED THE CODE WITH ' //
     *      'LIW = ' // XERN1, 2, 1)
         IDID = -33
      ENDIF
C
C        COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK
C        ARRAY AND RESTORE COMMON BLOCK DATA
C
      ICOMI = 21 + NEQ
      IINOUT = ICOMI + 33
C
      IYPOUT = 21
      ITSTAR = 21 + NEQ
      ICOMR = 22 + NEQ
C
      IF (INFO(1) .NE. 0) INTOUT = IWORK(IINOUT) .NE. (-1)
C     CALL DRSCO(RWORK(ICOMR),IWORK(ICOMI))
C
      IYH = ICOMR + 218
      IEWT = IYH + 6*NEQ
      ISAVF = IEWT + NEQ
      IACOR = ISAVF + NEQ
      IWM = IACOR + NEQ
      IDELSN = IWM + 2 + ILRW*NEQ
C
      IBEGIN = INFO(1)
      ITOL = INFO(2)
      IINTEG = INFO(3)
      ITSTOP = INFO(4)
      IJAC = INFO(5)
      IBAND = INFO(6)
      RWORK(ITSTAR) = T
C
      CALL DLSOD(DF,NEQ,T,Y,TOUT,RTOL,ATOL,IDID,RWORK(IYPOUT),
     1           RWORK(IYH),RWORK(IYH),RWORK(IEWT),RWORK(ISAVF),
     2           RWORK(IACOR),RWORK(IWM),IWORK(1),DJAC,INTOUT,
     3           RWORK(1),RWORK(12),RWORK(IDELSN),RPAR,IPAR)
C
      IWORK(IINOUT) = -1
      IF (INTOUT) IWORK(IINOUT) = 1
C
      IF (IDID .NE. (-2)) IWORK(LIW) = IWORK(LIW) + 1
      IF (T .NE. RWORK(ITSTAR)) IWORK(LIW) = 0
C     CALL DSVCO(RWORK(ICOMR),IWORK(ICOMI))
      RWORK(11) = H
      RWORK(13) = TN
      INFO(1) = IBEGIN
C
      RETURN
      END
*DECK DLSOD
      SUBROUTINE DLSOD (DF, NEQ, T, Y, TOUT, RTOL, ATOL, IDID, YPOUT,
     +   YH, YH1, EWT, SAVF, ACOR, WM, IWM, DJAC, INTOUT, TSTOP, TOLFAC,
     +   DELSGN, RPAR, IPAR)
C***BEGIN PROLOGUE  DLSOD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (LSOD-S, DLSOD-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   DDEBDF  merely allocates storage for  DLSOD  to relieve the user of
C   the inconvenience of a long call list.  Consequently  DLSOD  is used
C   as described in the comments for  DDEBDF .
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  D1MACH, DHSTRT, DINTYD, DSTOD, DVNRMS, XERMSG
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C***END PROLOGUE  DLSOD
C
      INTEGER IBAND, IBEGIN, IDID, IER, IINTEG, IJAC, INIT, INTFLG,
     1      IOWNS, IPAR, IQUIT, ITOL, ITSTOP, IWM, JSTART, K, KFLAG,
     2      KSTEPS, L, LACOR, LDUM, LEWT, LSAVF, LTOL, LWM, LYH, MAXNUM,
     3      MAXORD, METH, MITER, N, NATOLP, NEQ, NFE, NJE, NQ, NQU,
     4      NRTOLP, NST
      DOUBLE PRECISION ABSDEL, ACOR, ATOL, BIG, D1MACH, DEL,
     1      DELSGN, DT, DVNRMS, EL0, EWT,
     2      H, HA, HMIN, HMXI, HU, ROWNS, RPAR, RTOL, SAVF, T, TOL,
     3      TOLD, TOLFAC, TOUT, TSTOP, U, WM, X, Y, YH, YH1, YPOUT
      LOGICAL INTOUT
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
C
      DIMENSION Y(*),YPOUT(*),YH(NEQ,6),YH1(*),EWT(*),SAVF(*),
     1          ACOR(*),WM(*),IWM(*),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
C
C
      COMMON /DDEBD1/ TOLD,ROWNS(210),EL0,H,HMIN,HMXI,HU,X,U,IQUIT,INIT,
     1                LYH,LEWT,LACOR,LSAVF,LWM,KSTEPS,IBEGIN,ITOL,
     2                IINTEG,ITSTOP,IJAC,IBAND,IOWNS(6),IER,JSTART,
     3                KFLAG,LDUM,METH,MITER,MAXORD,N,NQ,NST,NFE,NJE,NQU
C
      EXTERNAL DF, DJAC
C
C     ..................................................................
C
C       THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
C       NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE
C       COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE
C       EXCESSIVE WORK.
      SAVE MAXNUM
C
      DATA MAXNUM /500/
C
C     ..................................................................
C
C***FIRST EXECUTABLE STATEMENT  DLSOD
      IF (IBEGIN .EQ. 0) THEN
C
C        ON THE FIRST CALL , PERFORM INITIALIZATION --
C        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
C        FUNCTION ROUTINE D1MACH. THE USER MUST MAKE SURE THAT THE
C        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
C
         U = D1MACH(4)
C                          -- SET ASSOCIATED MACHINE DEPENDENT PARAMETER
         WM(1) = SQRT(U)
C                          -- SET TERMINATION FLAG
         IQUIT = 0
C                          -- SET INITIALIZATION INDICATOR
         INIT = 0
C                          -- SET COUNTER FOR ATTEMPTED STEPS
         KSTEPS = 0
C                          -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
         INTOUT = .FALSE.
C                          -- SET START INDICATOR FOR DSTOD CODE
         JSTART = 0
C                          -- SET BDF METHOD INDICATOR
         METH = 2
C                          -- SET MAXIMUM ORDER FOR BDF METHOD
         MAXORD = 5
C                          -- SET ITERATION MATRIX INDICATOR
C
         IF (IJAC .EQ. 0 .AND. IBAND .EQ. 0) MITER = 2
         IF (IJAC .EQ. 1 .AND. IBAND .EQ. 0) MITER = 1
         IF (IJAC .EQ. 0 .AND. IBAND .EQ. 1) MITER = 5
         IF (IJAC .EQ. 1 .AND. IBAND .EQ. 1) MITER = 4
C
C                          -- SET OTHER NECESSARY ITEMS IN COMMON BLOCK
         N = NEQ
         NST = 0
         NJE = 0
         HMXI = 0.0D0
         NQ = 1
         H = 1.0D0
C                          -- RESET IBEGIN FOR SUBSEQUENT CALLS
         IBEGIN = 1
      ENDIF
C
C     ..................................................................
C
C      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
C
      IF (NEQ .LT. 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'DLSOD',
     *      'IN DDEBDF, THE NUMBER OF EQUATIONS MUST BE A ' //
     *      'POSITIVE INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = ' //
     *      XERN1, 6, 1)
         IDID=-33
      ENDIF
C
      NRTOLP = 0
      NATOLP = 0
      DO 60 K = 1, NEQ
         IF (NRTOLP .LE. 0) THEN
            IF (RTOL(K) .LT. 0.) THEN
               WRITE (XERN1, '(I8)') K
               WRITE (XERN3, '(1PE15.6)') RTOL(K)
               CALL XERMSG ('SLATEC', 'DLSOD',
     *            'IN DDEBDF, THE RELATIVE ERROR TOLERANCES MUST ' //
     *            'BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE WITH ' //
     *            'RTOL(' // XERN1 // ') = ' // XERN3 // '$$IN THE ' //
     *            'CASE OF VECTOR ERROR TOLERANCES, NO FURTHER ' //
     *            'CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
               IDID = -33
               IF (NATOLP .GT. 0) GO TO 70
               NRTOLP = 1
            ELSEIF (NATOLP .GT. 0) THEN
               GO TO 50
            ENDIF
         ENDIF
C
         IF (ATOL(K) .LT. 0.) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') ATOL(K)
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, THE ABSOLUTE ERROR ' //
     *         'TOLERANCES MUST BE NON-NEGATIVE.$$YOU HAVE CALLED ' //
     *         'THE CODE WITH ATOL(' // XERN1 // ') = ' // XERN3 //
     *         '$$IN THE CASE OF VECTOR ERROR TOLERANCES, NO FURTHER '
     *         // 'CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
            IDID=-33
            IF (NRTOLP .GT. 0) GO TO 70
            NATOLP=1
         ENDIF
   50    IF (ITOL .EQ. 0) GO TO 70
   60 CONTINUE
C
   70 IF (ITSTOP .EQ. 1) THEN
         IF (SIGN(1.0D0,TOUT-T) .NE. SIGN(1.0D0,TSTOP-T) .OR.
     1      ABS(TOUT-T) .GT. ABS(TSTOP-T)) THEN
            WRITE (XERN3, '(1PE15.6)') TOUT
            WRITE (XERN4, '(1PE15.6)') TSTOP
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, YOU HAVE CALLED THE ' //
     *         'CODE WITH TOUT = ' // XERN3 // '$$BUT YOU HAVE ' //
     *         'ALSO TOLD THE CODE NOT TO INTEGRATE PAST THE POINT ' //
     *         'TSTOP = ' // XERN4 // ' BY SETTING INFO(4) = 1.$$' //
     *         'THESE INSTRUCTIONS CONFLICT.', 14, 1)
            IDID=-33
         ENDIF
      ENDIF
C
C        CHECK SOME CONTINUATION POSSIBILITIES
C
      IF (INIT .NE. 0) THEN
         IF (T .EQ. TOUT) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, YOU HAVE CALLED THE CODE WITH T = TOUT = ' //
     *         XERN3 // '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.',
     *         9, 1)
            IDID=-33
         ENDIF
C
         IF (T .NE. TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, YOU HAVE CHANGED THE VALUE OF T FROM ' //
     *         XERN3 // ' TO ' // XERN4 //
     *         '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
            IDID=-33
         ENDIF
C
         IF (INIT .NE. 1) THEN
            IF (DELSGN*(TOUT-T) .LT. 0.0D0) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'DLSOD',
     *            'IN DDEBDF, BY CALLING THE CODE WITH TOUT = ' //
     *            XERN3 // ' YOU ARE ATTEMPTING TO CHANGE THE ' //
     *            'DIRECTION OF INTEGRATION.$$THIS IS NOT ALLOWED ' //
     *            'WITHOUT RESTARTING.', 11, 1)
               IDID=-33
            ENDIF
         ENDIF
      ENDIF
C
      IF (IDID .EQ. (-33)) THEN
         IF (IQUIT .NE. (-33)) THEN
C                       INVALID INPUT DETECTED
            IQUIT=-33
            IBEGIN=-1
         ELSE
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, INVALID INPUT WAS DETECTED ON ' //
     *         'SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED ' //
     *         'BECAUSE YOU HAVE NOT CORRECTED THE PROBLEM, ' //
     *         'SO EXECUTION IS BEING TERMINATED.', 12, 2)
         ENDIF
         RETURN
      ENDIF
C
C        ...............................................................
C
C             RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED
C             AS ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS
C             CASE, THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE
C             SMALLEST VALUE 100*U WHICH IS LIKELY TO BE REASONABLE FOR
C             THIS METHOD AND MACHINE
C
      DO 180 K = 1, NEQ
         IF (RTOL(K) + ATOL(K) .GT. 0.0D0) GO TO 170
            RTOL(K) = 100.0D0*U
            IDID = -2
  170    CONTINUE
C     ...EXIT
         IF (ITOL .EQ. 0) GO TO 190
  180 CONTINUE
  190 CONTINUE
C
      IF (IDID .NE. (-2)) GO TO 200
C        RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
C                                 SMALL POSITIVE VALUE
         IBEGIN = -1
      GO TO 460
  200 CONTINUE
C        BEGIN BLOCK PERMITTING ...EXITS TO 450
C           BEGIN BLOCK PERMITTING ...EXITS TO 430
C              BEGIN BLOCK PERMITTING ...EXITS TO 260
C                 BEGIN BLOCK PERMITTING ...EXITS TO 230
C
C                    BRANCH ON STATUS OF INITIALIZATION INDICATOR
C                           INIT=0 MEANS INITIAL DERIVATIVES AND
C                           NOMINAL STEP SIZE
C                                  AND DIRECTION NOT YET SET
C                           INIT=1 MEANS NOMINAL STEP SIZE AND
C                           DIRECTION NOT YET SET INIT=2 MEANS NO
C                           FURTHER INITIALIZATION REQUIRED
C
                     IF (INIT .EQ. 0) GO TO 210
C                 ......EXIT
                        IF (INIT .EQ. 1) GO TO 230
C              .........EXIT
                        GO TO 260
  210                CONTINUE
C
C                    ................................................
C
C                         MORE INITIALIZATION --
C                                             -- EVALUATE INITIAL
C                                             DERIVATIVES
C
                     INIT = 1
                     CALL DF(T,Y,YH(1,2),RPAR,IPAR)
                     NFE = 1
C                 ...EXIT
                     IF (T .NE. TOUT) GO TO 230
                     IDID = 2
                     DO 220 L = 1, NEQ
                        YPOUT(L) = YH(L,2)
  220                CONTINUE
                     TOLD = T
C        ............EXIT
                     GO TO 450
  230             CONTINUE
C
C                 -- COMPUTE INITIAL STEP SIZE
C                 -- SAVE SIGN OF INTEGRATION DIRECTION
C                 -- SET INDEPENDENT AND DEPENDENT VARIABLES
C                                      X AND YH(*) FOR DSTOD
C
                  LTOL = 1
                  DO 240 L = 1, NEQ
                     IF (ITOL .EQ. 1) LTOL = L
                     TOL = RTOL(LTOL)*ABS(Y(L)) + ATOL(LTOL)
                     IF (TOL .EQ. 0.0D0) GO TO 390
                     EWT(L) = TOL
  240             CONTINUE
C
                  BIG = SQRT(D1MACH(2))
                  CALL DHSTRT(DF,NEQ,T,TOUT,Y,YH(1,2),EWT,1,U,BIG,
     1                        YH(1,3),YH(1,4),YH(1,5),YH(1,6),RPAR,
     2                        IPAR,H)
C
                  DELSGN = SIGN(1.0D0,TOUT-T)
                  X = T
                  DO 250 L = 1, NEQ
                     YH(L,1) = Y(L)
                     YH(L,2) = H*YH(L,2)
  250             CONTINUE
                  INIT = 2
  260          CONTINUE
C
C              ......................................................
C
C                 ON EACH CALL SET INFORMATION WHICH DETERMINES THE
C                 ALLOWED INTERVAL OF INTEGRATION BEFORE RETURNING
C                 WITH AN ANSWER AT TOUT
C
               DEL = TOUT - T
               ABSDEL = ABS(DEL)
C
C              ......................................................
C
C                 IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND
C                 RETURN
C
  270          CONTINUE
C                 BEGIN BLOCK PERMITTING ...EXITS TO 400
C                    BEGIN BLOCK PERMITTING ...EXITS TO 380
                        IF (ABS(X-T) .LT. ABSDEL) GO TO 290
                           CALL DINTYD(TOUT,0,YH,NEQ,Y,INTFLG)
                           CALL DINTYD(TOUT,1,YH,NEQ,YPOUT,INTFLG)
                           IDID = 3
                           IF (X .NE. TOUT) GO TO 280
                              IDID = 2
                              INTOUT = .FALSE.
  280                      CONTINUE
                           T = TOUT
                           TOLD = T
C        ..................EXIT
                           GO TO 450
  290                   CONTINUE
C
C                       IF CANNOT GO PAST TSTOP AND SUFFICIENTLY
C                       CLOSE, EXTRAPOLATE AND RETURN
C
                        IF (ITSTOP .NE. 1) GO TO 310
                        IF (ABS(TSTOP-X) .GE. 100.0D0*U*ABS(X))
     1                     GO TO 310
                           DT = TOUT - X
                           DO 300 L = 1, NEQ
                              Y(L) = YH(L,1) + (DT/H)*YH(L,2)
  300                      CONTINUE
                           CALL DF(TOUT,Y,YPOUT,RPAR,IPAR)
                           NFE = NFE + 1
                           IDID = 3
                           T = TOUT
                           TOLD = T
C        ..................EXIT
                           GO TO 450
  310                   CONTINUE
C
                        IF (IINTEG .EQ. 0 .OR. .NOT.INTOUT) GO TO 320
C
C                          INTERMEDIATE-OUTPUT MODE
C
                           IDID = 1
                        GO TO 370
  320                   CONTINUE
C
C                       .............................................
C
C                            MONITOR NUMBER OF STEPS ATTEMPTED
C
                        IF (KSTEPS .LE. MAXNUM) GO TO 330
C
C                          A SIGNIFICANT AMOUNT OF WORK HAS BEEN
C                          EXPENDED
                           IDID = -1
                           KSTEPS = 0
                           IBEGIN = -1
                        GO TO 370
  330                   CONTINUE
C
C                          ..........................................
C
C                             LIMIT STEP SIZE AND SET WEIGHT VECTOR
C
                           HMIN = 100.0D0*U*ABS(X)
                           HA = MAX(ABS(H),HMIN)
                           IF (ITSTOP .EQ. 1)
     1                        HA = MIN(HA,ABS(TSTOP-X))
                           H = SIGN(HA,H)
                           LTOL = 1
                           DO 340 L = 1, NEQ
                              IF (ITOL .EQ. 1) LTOL = L
                              EWT(L) = RTOL(LTOL)*ABS(YH(L,1))
     1                                 + ATOL(LTOL)
C                    .........EXIT
                              IF (EWT(L) .LE. 0.0D0) GO TO 380
  340                      CONTINUE
                           TOLFAC = U*DVNRMS(NEQ,YH,EWT)
C                 .........EXIT
                           IF (TOLFAC .LE. 1.0D0) GO TO 400
C
C                          TOLERANCES TOO SMALL
                           IDID = -2
                           TOLFAC = 2.0D0*TOLFAC
                           RTOL(1) = TOLFAC*RTOL(1)
                           ATOL(1) = TOLFAC*ATOL(1)
                           IF (ITOL .EQ. 0) GO TO 360
                              DO 350 L = 2, NEQ
                                 RTOL(L) = TOLFAC*RTOL(L)
                                 ATOL(L) = TOLFAC*ATOL(L)
  350                         CONTINUE
  360                      CONTINUE
                           IBEGIN = -1
  370                   CONTINUE
C           ............EXIT
                        GO TO 430
  380                CONTINUE
C
C                    RELATIVE ERROR CRITERION INAPPROPRIATE
  390                CONTINUE
                     IDID = -3
                     IBEGIN = -1
C           .........EXIT
                     GO TO 430
  400             CONTINUE
C
C                 ...................................................
C
C                      TAKE A STEP
C
                  CALL DSTOD(NEQ,Y,YH,NEQ,YH1,EWT,SAVF,ACOR,WM,IWM,
     1                       DF,DJAC,RPAR,IPAR)
C
                  JSTART = -2
                  INTOUT = .TRUE.
               IF (KFLAG .EQ. 0) GO TO 270
C
C              ......................................................
C
               IF (KFLAG .EQ. -1) GO TO 410
C
C                 REPEATED CORRECTOR CONVERGENCE FAILURES
                  IDID = -6
                  IBEGIN = -1
               GO TO 420
  410          CONTINUE
C
C                 REPEATED ERROR TEST FAILURES
                  IDID = -7
                  IBEGIN = -1
  420          CONTINUE
  430       CONTINUE
C
C           .........................................................
C
C                                  STORE VALUES BEFORE RETURNING TO
C                                  DDEBDF
            DO 440 L = 1, NEQ
               Y(L) = YH(L,1)
               YPOUT(L) = YH(L,2)/H
  440       CONTINUE
            T = X
            TOLD = T
            INTOUT = .FALSE.
  450    CONTINUE
  460 CONTINUE
      RETURN
      END
*DECK DSTOD
      SUBROUTINE DSTOD (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM,
     +   DF, DJAC, RPAR, IPAR)
C***BEGIN PROLOGUE  DSTOD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (STOD-S, DSTOD-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DSTOD integrates a system of first order odes over one step in the
C   integrator package DDEBDF.
C ----------------------------------------------------------------------
C DSTOD  performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note.. DSTOD  is independent of the value of the iteration method
C indicator MITER, when this is .NE. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTOD  is done with the following variables..
C
C Y      = An array of length .GE. N used as the Y argument in
C          all calls to DF and DJAC.
C NEQ    = Integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to DF and DJAC.
C YH     = An NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
C          J-th derivative of Y(I), scaled by H**J/FACTORIAL(J)
C          (J = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = A constant integer .GE. N, the first dimension of YH.
C YH1    = A one-dimensional array occupying the same space as YH.
C EWT    = An array of N elements with which the estimated local
C          errors in YH are compared.
C SAVF   = An array of working storage, of length N.
C ACOR   = A work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(I) contains
C          the estimated one-step local error in Y(I).
C WM,IWM = DOUBLE PRECISION and INTEGER work arrays associated with
C          matrix operations in chord iteration (MITER .NE. 0).
C DPJAC   = Name of routine to evaluate and preprocess Jacobian matrix
C          if a chord method is being used.
C DSLVS   = Name of routine to solve linear system in chord iteration.
C H      = The step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = The minimum absolute value of the step size H to be used.
C HMXI   = Inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = The independent variable. TN is updated on each step taken.
C JSTART = An integer used for input only, with the following
C          values and meanings..
C               0  Perform the first step.
C           .GT.0  Take a new step continuing from the last.
C              -1  Take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  Take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings..
C               0  The step was successful.
C              -1  The requested error could not be achieved.
C              -2  Corrector convergence could not be achieved.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = The maximum order of integration method to be allowed.
C METH/MITER = The method flags.  See description in driver.
C N      = The number of first-order differential equations.
C ----------------------------------------------------------------------
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  DCFOD, DPJAC, DSLVS, DVNRMS
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920422  Changed DIMENSION statement.  (WRB)
C***END PROLOGUE  DSTOD
C
      INTEGER I, I1, IALTH, IER, IOD, IOWND, IPAR, IPUP, IREDO, IRET,
     1      IWM, J, JB, JSTART, KFLAG, KSTEPS, L, LMAX, M, MAXORD,
     2      MEO, METH, MITER, N, NCF, NEQ, NEWQ, NFE, NJE, NQ, NQNYH,
     3      NQU, NST, NSTEPJ, NYH
      DOUBLE PRECISION ACOR, CONIT, CRATE, DCON, DDN,
     1      DEL, DELP, DSM, DUP, DVNRMS, EL, EL0, ELCO,
     2      EWT, EXDN, EXSM, EXUP, H, HMIN, HMXI, HOLD, HU, R, RC,
     3      RH, RHDN, RHSM, RHUP, RMAX, ROWND, RPAR, SAVF, TESCO,
     4      TN, TOLD, UROUND, WM, Y, YH, YH1
      EXTERNAL DF, DJAC
C
      DIMENSION Y(*),YH(NYH,*),YH1(*),EWT(*),SAVF(*),ACOR(*),WM(*),
     1          IWM(*),RPAR(*),IPAR(*)
      COMMON /DDEBD1/ ROWND,CONIT,CRATE,EL(13),ELCO(13,12),HOLD,RC,RMAX,
     1                TESCO(3,12),EL0,H,HMIN,HMXI,HU,TN,UROUND,IOWND(7),
     2                KSTEPS,IOD(6),IALTH,IPUP,LMAX,MEO,NQNYH,NSTEPJ,
     3                IER,JSTART,KFLAG,L,METH,MITER,MAXORD,N,NQ,NST,NFE,
     4                NJE,NQU
C
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 690
C        BEGIN BLOCK PERMITTING ...EXITS TO 60
C***FIRST EXECUTABLE STATEMENT  DSTOD
            KFLAG = 0
            TOLD = TN
            NCF = 0
            IF (JSTART .GT. 0) GO TO 160
            IF (JSTART .EQ. -1) GO TO 10
               IF (JSTART .EQ. -2) GO TO 90
C              ---------------------------------------------------------
C               ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER
C               VARIABLES ARE INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY
C               WHICH H CAN BE INCREASED IN A SINGLE STEP.  IT IS
C               INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL INITIAL H,
C               BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE OCCURS
C               (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT
C               2 FOR THE NEXT INCREASE.
C              ---------------------------------------------------------
               LMAX = MAXORD + 1
               NQ = 1
               L = 2
               IALTH = 2
               RMAX = 10000.0D0
               RC = 0.0D0
               EL0 = 1.0D0
               CRATE = 0.7D0
               DELP = 0.0D0
               HOLD = H
               MEO = METH
               NSTEPJ = 0
               IRET = 3
            GO TO 50
   10       CONTINUE
C              BEGIN BLOCK PERMITTING ...EXITS TO 30
C                 ------------------------------------------------------
C                  THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN
C                  JSTART = -1.  IPUP IS SET TO MITER TO FORCE A MATRIX
C                  UPDATE.  IF AN ORDER INCREASE IS ABOUT TO BE
C                  CONSIDERED (IALTH = 1), IALTH IS RESET TO 2 TO
C                  POSTPONE CONSIDERATION ONE MORE STEP.  IF THE CALLER
C                  HAS CHANGED METH, DCFOD  IS CALLED TO RESET THE
C                  COEFFICIENTS OF THE METHOD.  IF THE CALLER HAS
C                  CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
C                  ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN
C                  ACCORDINGLY.  IF H IS TO BE CHANGED, YH MUST BE
C                  RESCALED.  IF H OR METH IS BEING CHANGED, IALTH IS
C                  RESET TO L = NQ + 1 TO PREVENT FURTHER CHANGES IN H
C                  FOR THAT MANY STEPS.
C                 ------------------------------------------------------
                  IPUP = MITER
                  LMAX = MAXORD + 1
                  IF (IALTH .EQ. 1) IALTH = 2
                  IF (METH .EQ. MEO) GO TO 20
                     CALL DCFOD(METH,ELCO,TESCO)
                     MEO = METH
C              ......EXIT
                     IF (NQ .GT. MAXORD) GO TO 30
                     IALTH = L
                     IRET = 1
C        ............EXIT
                     GO TO 60
   20             CONTINUE
                  IF (NQ .LE. MAXORD) GO TO 90
   30          CONTINUE
               NQ = MAXORD
               L = LMAX
               DO 40 I = 1, L
                  EL(I) = ELCO(I,NQ)
   40          CONTINUE
               NQNYH = NQ*NYH
               RC = RC*EL(1)/EL0
               EL0 = EL(1)
               CONIT = 0.5D0/(NQ+2)
               DDN = DVNRMS(N,SAVF,EWT)/TESCO(1,L)
               EXDN = 1.0D0/L
               RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
               RH = MIN(RHDN,1.0D0)
               IREDO = 3
               IF (H .EQ. HOLD) GO TO 660
               RH = MIN(RH,ABS(H/HOLD))
               H = HOLD
               GO TO 100
   50       CONTINUE
C           ------------------------------------------------------------
C            DCFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS
C            FOR THE CURRENT METH.  THEN THE EL VECTOR AND RELATED
C            CONSTANTS ARE RESET WHENEVER THE ORDER NQ IS CHANGED, OR AT
C            THE START OF THE PROBLEM.
C           ------------------------------------------------------------
            CALL DCFOD(METH,ELCO,TESCO)
   60    CONTINUE
   70    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 680
               DO 80 I = 1, L
                  EL(I) = ELCO(I,NQ)
   80          CONTINUE
               NQNYH = NQ*NYH
               RC = RC*EL(1)/EL0
               EL0 = EL(1)
               CONIT = 0.5D0/(NQ+2)
               GO TO (90,660,160), IRET
C              ---------------------------------------------------------
C               IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
C               RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH
C               IS SET TO L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT
C               MANY STEPS, UNLESS FORCED BY A CONVERGENCE OR ERROR TEST
C               FAILURE.
C              ---------------------------------------------------------
   90          CONTINUE
               IF (H .EQ. HOLD) GO TO 160
               RH = H/HOLD
               H = HOLD
               IREDO = 3
  100          CONTINUE
  110          CONTINUE
                  RH = MIN(RH,RMAX)
                  RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
                  R = 1.0D0
                  DO 130 J = 2, L
                     R = R*RH
                     DO 120 I = 1, N
                        YH(I,J) = YH(I,J)*R
  120                CONTINUE
  130             CONTINUE
                  H = H*RH
                  RC = RC*RH
                  IALTH = L
                  IF (IREDO .NE. 0) GO TO 150
                     RMAX = 10.0D0
                     R = 1.0D0/TESCO(2,NQU)
                     DO 140 I = 1, N
                        ACOR(I) = ACOR(I)*R
  140                CONTINUE
C     ...............EXIT
                     GO TO 690
  150             CONTINUE
C                 ------------------------------------------------------
C                  THIS SECTION COMPUTES THE PREDICTED VALUES BY
C                  EFFECTIVELY MULTIPLYING THE YH ARRAY BY THE PASCAL
C                  TRIANGLE MATRIX.  RC IS THE RATIO OF NEW TO OLD
C                  VALUES OF THE COEFFICIENT  H*EL(1).  WHEN RC DIFFERS
C                  FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER
C                  TO FORCE DPJAC TO BE CALLED, IF A JACOBIAN IS
C                  INVOLVED.  IN ANY CASE, DPJAC IS CALLED AT LEAST
C                  EVERY 20-TH STEP.
C                 ------------------------------------------------------
  160             CONTINUE
  170             CONTINUE
C                    BEGIN BLOCK PERMITTING ...EXITS TO 610
C                       BEGIN BLOCK PERMITTING ...EXITS TO 490
                           IF (ABS(RC-1.0D0) .GT. 0.3D0) IPUP = MITER
                           IF (NST .GE. NSTEPJ + 20) IPUP = MITER
                           TN = TN + H
                           I1 = NQNYH + 1
                           DO 190 JB = 1, NQ
                              I1 = I1 - NYH
                              DO 180 I = I1, NQNYH
                                 YH1(I) = YH1(I) + YH1(I+NYH)
  180                         CONTINUE
  190                      CONTINUE
                           KSTEPS = KSTEPS + 1
C                          ---------------------------------------------
C                           UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A
C                           CONVERGENCE TEST IS MADE ON THE R.M.S. NORM
C                           OF EACH CORRECTION, WEIGHTED BY THE ERROR
C                           WEIGHT VECTOR EWT.  THE SUM OF THE
C                           CORRECTIONS IS ACCUMULATED IN THE VECTOR
C                           ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE
C                           CORRECTOR LOOP.
C                          ---------------------------------------------
  200                      CONTINUE
                              M = 0
                              DO 210 I = 1, N
                                 Y(I) = YH(I,1)
  210                         CONTINUE
                              CALL DF(TN,Y,SAVF,RPAR,IPAR)
                              NFE = NFE + 1
                              IF (IPUP .LE. 0) GO TO 220
C                                ---------------------------------------
C                                 IF INDICATED, THE MATRIX P = I -
C                                 H*EL(1)*J IS REEVALUATED AND
C                                 PREPROCESSED BEFORE STARTING THE
C                                 CORRECTOR ITERATION.  IPUP IS SET TO 0
C                                 AS AN INDICATOR THAT THIS HAS BEEN
C                                 DONE.
C                                ---------------------------------------
                                 IPUP = 0
                                 RC = 1.0D0
                                 NSTEPJ = NST
                                 CRATE = 0.7D0
                                 CALL DPJAC(NEQ,Y,YH,NYH,EWT,ACOR,SAVF,
     1                                      WM,IWM,DF,DJAC,RPAR,IPAR)
C                          ......EXIT
                                 IF (IER .NE. 0) GO TO 440
  220                         CONTINUE
                              DO 230 I = 1, N
                                 ACOR(I) = 0.0D0
  230                         CONTINUE
  240                         CONTINUE
                                 IF (MITER .NE. 0) GO TO 270
C                                   ------------------------------------
C                                    IN THE CASE OF FUNCTIONAL
C                                    ITERATION, UPDATE Y DIRECTLY FROM
C                                    THE RESULT OF THE LAST FUNCTION
C                                    EVALUATION.
C                                   ------------------------------------
                                    DO 250 I = 1, N
                                       SAVF(I) = H*SAVF(I) - YH(I,2)
                                       Y(I) = SAVF(I) - ACOR(I)
  250                               CONTINUE
                                    DEL = DVNRMS(N,Y,EWT)
                                    DO 260 I = 1, N
                                       Y(I) = YH(I,1) + EL(1)*SAVF(I)
                                       ACOR(I) = SAVF(I)
  260                               CONTINUE
                                 GO TO 300
  270                            CONTINUE
C                                   ------------------------------------
C                                    IN THE CASE OF THE CHORD METHOD,
C                                    COMPUTE THE CORRECTOR ERROR, AND
C                                    SOLVE THE LINEAR SYSTEM WITH THAT
C                                    AS RIGHT-HAND SIDE AND P AS
C                                    COEFFICIENT MATRIX.
C                                   ------------------------------------
                                    DO 280 I = 1, N
                                       Y(I) = H*SAVF(I)
     1                                        - (YH(I,2) + ACOR(I))
  280                               CONTINUE
                                    CALL DSLVS(WM,IWM,Y,SAVF)
C                             ......EXIT
                                    IF (IER .NE. 0) GO TO 430
                                    DEL = DVNRMS(N,Y,EWT)
                                    DO 290 I = 1, N
                                       ACOR(I) = ACOR(I) + Y(I)
                                       Y(I) = YH(I,1) + EL(1)*ACOR(I)
  290                               CONTINUE
  300                            CONTINUE
C                                ---------------------------------------
C                                 TEST FOR CONVERGENCE.  IF M.GT.0, AN
C                                 ESTIMATE OF THE CONVERGENCE RATE
C                                 CONSTANT IS STORED IN CRATE, AND THIS
C                                 IS USED IN THE TEST.
C                                ---------------------------------------
                                 IF (M .NE. 0)
     1                              CRATE = MAX(0.2D0*CRATE,DEL/DELP)
                                 DCON = DEL*MIN(1.0D0,1.5D0*CRATE)
     1                                  /(TESCO(2,NQ)*CONIT)
                                 IF (DCON .GT. 1.0D0) GO TO 420
C                                   ------------------------------------
C                                    THE CORRECTOR HAS CONVERGED.  IPUP
C                                    IS SET TO -1 IF MITER .NE. 0, TO
C                                    SIGNAL THAT THE JACOBIAN INVOLVED
C                                    MAY NEED UPDATING LATER.  THE LOCAL
C                                    ERROR TEST IS MADE AND CONTROL
C                                    PASSES TO STATEMENT 500 IF IT
C                                    FAILS.
C                                   ------------------------------------
                                    IF (MITER .NE. 0) IPUP = -1
                                    IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
                                    IF (M .GT. 0)
     1                                 DSM = DVNRMS(N,ACOR,EWT)
     2                                       /TESCO(2,NQ)
                                    IF (DSM .GT. 1.0D0) GO TO 380
C                                      BEGIN BLOCK
C                                      PERMITTING ...EXITS TO 360
C                                         ------------------------------
C                                          AFTER A SUCCESSFUL STEP,
C                                          UPDATE THE YH ARRAY.
C                                          CONSIDER CHANGING H IF IALTH
C                                          = 1.  OTHERWISE DECREASE
C                                          IALTH BY 1.  IF IALTH IS THEN
C                                          1 AND NQ .LT. MAXORD, THEN
C                                          ACOR IS SAVED FOR USE IN A
C                                          POSSIBLE ORDER INCREASE ON
C                                          THE NEXT STEP.  IF A CHANGE
C                                          IN H IS CONSIDERED, AN
C                                          INCREASE OR DECREASE IN ORDER
C                                          BY ONE IS CONSIDERED ALSO.  A
C                                          CHANGE IN H IS MADE ONLY IF
C                                          IT IS BY A FACTOR OF AT LEAST
C                                          1.1.  IF NOT, IALTH IS SET TO
C                                          3 TO PREVENT TESTING FOR THAT
C                                          MANY STEPS.
C                                         ------------------------------
                                          KFLAG = 0
                                          IREDO = 0
                                          NST = NST + 1
                                          HU = H
                                          NQU = NQ
                                          DO 320 J = 1, L
                                             DO 310 I = 1, N
                                                YH(I,J) = YH(I,J)
     1                                                    + EL(J)
     2                                                      *ACOR(I)
  310                                        CONTINUE
  320                                     CONTINUE
                                          IALTH = IALTH - 1
                                          IF (IALTH .NE. 0) GO TO 340
C                                            ---------------------------
C                                             REGARDLESS OF THE SUCCESS
C                                             OR FAILURE OF THE STEP,
C                                             FACTORS RHDN, RHSM, AND
C                                             RHUP ARE COMPUTED, BY
C                                             WHICH H COULD BE
C                                             MULTIPLIED AT ORDER NQ -
C                                             1, ORDER NQ, OR ORDER NQ +
C                                             1, RESPECTIVELY.  IN THE
C                                             CASE OF FAILURE, RHUP =
C                                             0.0 TO AVOID AN ORDER
C                                             INCREASE.  THE LARGEST OF
C                                             THESE IS DETERMINED AND
C                                             THE NEW ORDER CHOSEN
C                                             ACCORDINGLY.  IF THE ORDER
C                                             IS TO BE INCREASED, WE
C                                             COMPUTE ONE ADDITIONAL
C                                             SCALED DERIVATIVE.
C                                            ---------------------------
                                             RHUP = 0.0D0
C                       .....................EXIT
                                             IF (L .EQ. LMAX) GO TO 490
                                             DO 330 I = 1, N
                                                SAVF(I) = ACOR(I)
     1                                                    - YH(I,LMAX)
  330                                        CONTINUE
                                             DUP = DVNRMS(N,SAVF,EWT)
     1                                             /TESCO(3,NQ)
                                             EXUP = 1.0D0/(L+1)
                                             RHUP = 1.0D0
     1                                              /(1.4D0*DUP**EXUP
     2                                                + 0.0000014D0)
C                       .....................EXIT
                                             GO TO 490
  340                                     CONTINUE
C                                      ...EXIT
                                          IF (IALTH .GT. 1) GO TO 360
C                                      ...EXIT
                                          IF (L .EQ. LMAX) GO TO 360
                                          DO 350 I = 1, N
                                             YH(I,LMAX) = ACOR(I)
  350                                     CONTINUE
  360                                  CONTINUE
                                       R = 1.0D0/TESCO(2,NQU)
                                       DO 370 I = 1, N
                                          ACOR(I) = ACOR(I)*R
  370                                  CONTINUE
C     .................................EXIT
                                       GO TO 690
  380                               CONTINUE
C                                   ------------------------------------
C                                    THE ERROR TEST FAILED.  KFLAG KEEPS
C                                    TRACK OF MULTIPLE FAILURES.
C                                    RESTORE TN AND THE YH ARRAY TO
C                                    THEIR PREVIOUS VALUES, AND PREPARE
C                                    TO TRY THE STEP AGAIN.  COMPUTE THE
C                                    OPTIMUM STEP SIZE FOR THIS OR ONE
C                                    LOWER ORDER.  AFTER 2 OR MORE
C                                    FAILURES, H IS FORCED TO DECREASE
C                                    BY A FACTOR OF 0.2 OR LESS.
C                                   ------------------------------------
                                    KFLAG = KFLAG - 1
                                    TN = TOLD
                                    I1 = NQNYH + 1
                                    DO 400 JB = 1, NQ
                                       I1 = I1 - NYH
                                       DO 390 I = I1, NQNYH
                                          YH1(I) = YH1(I) - YH1(I+NYH)
  390                                  CONTINUE
  400                               CONTINUE
                                    RMAX = 2.0D0
                                    IF (ABS(H) .GT. HMIN*1.00001D0)
     1                                 GO TO 410
C                                      ---------------------------------
C                                       ALL RETURNS ARE MADE THROUGH
C                                       THIS SECTION.  H IS SAVED IN
C                                       HOLD TO ALLOW THE CALLER TO
C                                       CHANGE H ON THE NEXT STEP.
C                                      ---------------------------------
                                       KFLAG = -1
C     .................................EXIT
                                       GO TO 690
  410                               CONTINUE
C                    ...............EXIT
                                    IF (KFLAG .LE. -3) GO TO 610
                                    IREDO = 2
                                    RHUP = 0.0D0
C                       ............EXIT
                                    GO TO 490
  420                            CONTINUE
                                 M = M + 1
C                             ...EXIT
                                 IF (M .EQ. 3) GO TO 430
C                             ...EXIT
                                 IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP)
     1                              GO TO 430
                                 DELP = DEL
                                 CALL DF(TN,Y,SAVF,RPAR,IPAR)
                                 NFE = NFE + 1
                              GO TO 240
  430                         CONTINUE
C                             ------------------------------------------
C                              THE CORRECTOR ITERATION FAILED TO
C                              CONVERGE IN 3 TRIES.  IF MITER .NE. 0 AND
C                              THE JACOBIAN IS OUT OF DATE, DPJAC IS
C                              CALLED FOR THE NEXT TRY.  OTHERWISE THE
C                              YH ARRAY IS RETRACTED TO ITS VALUES
C                              BEFORE PREDICTION, AND H IS REDUCED, IF
C                              POSSIBLE.  IF H CANNOT BE REDUCED OR 10
C                              FAILURES HAVE OCCURRED, EXIT WITH KFLAG =
C                              -2.
C                             ------------------------------------------
C                          ...EXIT
                              IF (IPUP .EQ. 0) GO TO 440
                              IPUP = MITER
                           GO TO 200
  440                      CONTINUE
                           TN = TOLD
                           NCF = NCF + 1
                           RMAX = 2.0D0
                           I1 = NQNYH + 1
                           DO 460 JB = 1, NQ
                              I1 = I1 - NYH
                              DO 450 I = I1, NQNYH
                                 YH1(I) = YH1(I) - YH1(I+NYH)
  450                         CONTINUE
  460                      CONTINUE
                           IF (ABS(H) .GT. HMIN*1.00001D0) GO TO 470
                              KFLAG = -2
C     ........................EXIT
                              GO TO 690
  470                      CONTINUE
                           IF (NCF .NE. 10) GO TO 480
                              KFLAG = -2
C     ........................EXIT
                              GO TO 690
  480                      CONTINUE
                           RH = 0.25D0
                           IPUP = MITER
                           IREDO = 1
C                 .........EXIT
                           GO TO 650
  490                   CONTINUE
                        EXSM = 1.0D0/L
                        RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
                        RHDN = 0.0D0
                        IF (NQ .EQ. 1) GO TO 500
                           DDN = DVNRMS(N,YH(1,L),EWT)/TESCO(1,NQ)
                           EXDN = 1.0D0/NQ
                           RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
  500                   CONTINUE
                        IF (RHSM .GE. RHUP) GO TO 550
                           IF (RHUP .LE. RHDN) GO TO 540
                              NEWQ = L
                              RH = RHUP
                              IF (RH .GE. 1.1D0) GO TO 520
                                 IALTH = 3
                                 R = 1.0D0/TESCO(2,NQU)
                                 DO 510 I = 1, N
                                    ACOR(I) = ACOR(I)*R
  510                            CONTINUE
C     ...........................EXIT
                                 GO TO 690
  520                         CONTINUE
                              R = EL(L)/L
                              DO 530 I = 1, N
                                 YH(I,NEWQ+1) = ACOR(I)*R
  530                         CONTINUE
                              NQ = NEWQ
                              L = NQ + 1
                              IRET = 2
C           ..................EXIT
                              GO TO 680
  540                      CONTINUE
                        GO TO 580
  550                   CONTINUE
                        IF (RHSM .LT. RHDN) GO TO 580
                           NEWQ = NQ
                           RH = RHSM
                           IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0)
     1                        GO TO 560
                              IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C                             ------------------------------------------
C                              IF THERE IS A CHANGE OF ORDER, RESET NQ,
C                              L, AND THE COEFFICIENTS.  IN ANY CASE H
C                              IS RESET ACCORDING TO RH AND THE YH ARRAY
C                              IS RESCALED.  THEN EXIT FROM 680 IF THE
C                              STEP WAS OK, OR REDO THE STEP OTHERWISE.
C                             ------------------------------------------
C                 ............EXIT
                              IF (NEWQ .EQ. NQ) GO TO 650
                              NQ = NEWQ
                              L = NQ + 1
                              IRET = 2
C           ..................EXIT
                              GO TO 680
  560                      CONTINUE
                           IALTH = 3
                           R = 1.0D0/TESCO(2,NQU)
                           DO 570 I = 1, N
                              ACOR(I) = ACOR(I)*R
  570                      CONTINUE
C     .....................EXIT
                           GO TO 690
  580                   CONTINUE
                        NEWQ = NQ - 1
                        RH = RHDN
                        IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
                        IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 590
                           IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C                          ---------------------------------------------
C                           IF THERE IS A CHANGE OF ORDER, RESET NQ, L,
C                           AND THE COEFFICIENTS.  IN ANY CASE H IS
C                           RESET ACCORDING TO RH AND THE YH ARRAY IS
C                           RESCALED.  THEN EXIT FROM 680 IF THE STEP
C                           WAS OK, OR REDO THE STEP OTHERWISE.
C                          ---------------------------------------------
C                 .........EXIT
                           IF (NEWQ .EQ. NQ) GO TO 650
                           NQ = NEWQ
                           L = NQ + 1
                           IRET = 2
C           ...............EXIT
                           GO TO 680
  590                   CONTINUE
                        IALTH = 3
                        R = 1.0D0/TESCO(2,NQU)
                        DO 600 I = 1, N
                           ACOR(I) = ACOR(I)*R
  600                   CONTINUE
C     ..................EXIT
                        GO TO 690
  610                CONTINUE
C                    ---------------------------------------------------
C                     CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES
C                     HAVE OCCURRED.  IF 10 FAILURES HAVE OCCURRED, EXIT
C                     WITH KFLAG = -1.  IT IS ASSUMED THAT THE
C                     DERIVATIVES THAT HAVE ACCUMULATED IN THE YH ARRAY
C                     HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST
C                     DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO
C                     1.  THEN H IS REDUCED BY A FACTOR OF 10, AND THE
C                     STEP IS RETRIED, UNTIL IT SUCCEEDS OR H REACHES
C                     HMIN.
C                    ---------------------------------------------------
                     IF (KFLAG .NE. -10) GO TO 620
C                       ------------------------------------------------
C                        ALL RETURNS ARE MADE THROUGH THIS SECTION.  H
C                        IS SAVED IN HOLD TO ALLOW THE CALLER TO CHANGE
C                        H ON THE NEXT STEP.
C                       ------------------------------------------------
                        KFLAG = -1
C     ..................EXIT
                        GO TO 690
  620                CONTINUE
                     RH = 0.1D0
                     RH = MAX(HMIN/ABS(H),RH)
                     H = H*RH
                     DO 630 I = 1, N
                        Y(I) = YH(I,1)
  630                CONTINUE
                     CALL DF(TN,Y,SAVF,RPAR,IPAR)
                     NFE = NFE + 1
                     DO 640 I = 1, N
                        YH(I,2) = H*SAVF(I)
  640                CONTINUE
                     IPUP = MITER
                     IALTH = 5
C              ......EXIT
                     IF (NQ .NE. 1) GO TO 670
                  GO TO 170
  650             CONTINUE
  660             CONTINUE
                  RH = MAX(RH,HMIN/ABS(H))
               GO TO 110
  670          CONTINUE
               NQ = 1
               L = 2
               IRET = 3
  680       CONTINUE
         GO TO 70
  690 CONTINUE
      HOLD = H
      JSTART = 1
      RETURN
C     ----------------------- END OF SUBROUTINE DSTOD
C     -----------------------
      END
*DECK DCFOD
      SUBROUTINE DCFOD (METH, ELCO, TESCO)
C***BEGIN PROLOGUE  DCFOD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (CFOD-S, DCFOD-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   DCFOD defines coefficients needed in the integrator package DDEBDF
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DCFOD
C
C
      INTEGER I, IB, METH, NQ, NQM1, NQP1
      DOUBLE PRECISION AGAMQ, ELCO, FNQ, FNQM1, PC, PINT, RAGQ,
     1      RQ1FAC, RQFAC, TESCO, TSIGN, XPIN
      DIMENSION ELCO(13,12),TESCO(3,12)
C     ------------------------------------------------------------------
C      DCFOD  IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS
C      NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS
C      GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.
C      THE MAXIMUM ORDER ASSUMED HERE IS 12 IF METH = 1 AND 5 IF METH =
C      2.  (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)
C      DCFOD  IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM,
C      AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.
C
C      THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.
C      THE COEFFICIENTS EL(I), 1 .LE. I .LE. NQ+1, FOR THE METHOD OF
C      ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A
C      GENERATING POLYNOMIAL, I.E.,
C          L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.
C      FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY
C          DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) =
C      0.  FOR THE BDF METHODS, L(X) IS GIVEN BY
C          L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,
C      WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
C
C      THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE
C      LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.
C      AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP
C      SIZE AT ORDER NQ - 1 IF K = 1, AT ORDER NQ IF K = 2, AND AT ORDER
C      NQ + 1 IF K = 3.
C     ------------------------------------------------------------------
      DIMENSION PC(12)
C
C***FIRST EXECUTABLE STATEMENT  DCFOD
      GO TO (10,60), METH
C
   10 CONTINUE
         ELCO(1,1) = 1.0D0
         ELCO(2,1) = 1.0D0
         TESCO(1,1) = 0.0D0
         TESCO(2,1) = 2.0D0
         TESCO(1,2) = 1.0D0
         TESCO(3,12) = 0.0D0
         PC(1) = 1.0D0
         RQFAC = 1.0D0
         DO 50 NQ = 2, 12
C           ------------------------------------------------------------
C            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
C                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ-1).
C            INITIALLY, P(X) = 1.
C           ------------------------------------------------------------
            RQ1FAC = RQFAC
            RQFAC = RQFAC/NQ
            NQM1 = NQ - 1
            FNQM1 = NQM1
            NQP1 = NQ + 1
C           FORM COEFFICIENTS OF P(X)*(X+NQ-1).
C           ----------------------------------
            PC(NQ) = 0.0D0
            DO 20 IB = 1, NQM1
               I = NQP1 - IB
               PC(I) = PC(I-1) + FNQM1*PC(I)
   20       CONTINUE
            PC(1) = FNQM1*PC(1)
C           COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X).
C           -----------------------
            PINT = PC(1)
            XPIN = PC(1)/2.0D0
            TSIGN = 1.0D0
            DO 30 I = 2, NQ
               TSIGN = -TSIGN
               PINT = PINT + TSIGN*PC(I)/I
               XPIN = XPIN + TSIGN*PC(I)/(I+1)
   30       CONTINUE
C           STORE COEFFICIENTS IN ELCO AND TESCO.
C           --------------------------------
            ELCO(1,NQ) = PINT*RQ1FAC
            ELCO(2,NQ) = 1.0D0
            DO 40 I = 2, NQ
               ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
   40       CONTINUE
            AGAMQ = RQFAC*XPIN
            RAGQ = 1.0D0/AGAMQ
            TESCO(2,NQ) = RAGQ
            IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
            TESCO(3,NQM1) = RAGQ
   50    CONTINUE
      GO TO 100
C
   60 CONTINUE
         PC(1) = 1.0D0
         RQ1FAC = 1.0D0
         DO 90 NQ = 1, 5
C           ------------------------------------------------------------
C            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
C                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ).
C            INITIALLY, P(X) = 1.
C           ------------------------------------------------------------
            FNQ = NQ
            NQP1 = NQ + 1
C           FORM COEFFICIENTS OF P(X)*(X+NQ).
C           ------------------------------------
            PC(NQP1) = 0.0D0
            DO 70 IB = 1, NQ
               I = NQ + 2 - IB
               PC(I) = PC(I-1) + FNQ*PC(I)
   70       CONTINUE
            PC(1) = FNQ*PC(1)
C           STORE COEFFICIENTS IN ELCO AND TESCO.
C           --------------------------------
            DO 80 I = 1, NQP1
               ELCO(I,NQ) = PC(I)/PC(2)
   80       CONTINUE
            ELCO(2,NQ) = 1.0D0
            TESCO(1,NQ) = RQ1FAC
            TESCO(2,NQ) = NQP1/ELCO(1,NQ)
            TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
            RQ1FAC = RQ1FAC/FNQ
   90    CONTINUE
  100 CONTINUE
      RETURN
C     ----------------------- END OF SUBROUTINE DCFOD
C     -----------------------
      END
*DECK DVNRMS
      DOUBLE PRECISION FUNCTION DVNRMS (N, V, W)
C***BEGIN PROLOGUE  DVNRMS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (VNWRMS-S, DVNRMS-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   DVNRMS computes a weighted root-mean-square vector norm for the
C   integrator package DDEBDF.
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DVNRMS
      INTEGER I, N
      DOUBLE PRECISION SUM, V, W
      DIMENSION V(*),W(*)
C***FIRST EXECUTABLE STATEMENT  DVNRMS
      SUM = 0.0D0
      DO 10 I = 1, N
         SUM = SUM + (V(I)/W(I))**2
   10 CONTINUE
      DVNRMS = SQRT(SUM/N)
      RETURN
C     ----------------------- END OF FUNCTION DVNRMS
C     ------------------------
      END
*DECK DINTYD
      SUBROUTINE DINTYD (T, K, YH, NYH, DKY, IFLAG)
C***BEGIN PROLOGUE  DINTYD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (INTYD-S, DINTYD-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DINTYD approximates the solution and derivatives at T by polynomial
C   interpolation. Must be used in conjunction with the integrator
C   package DDEBDF.
C ----------------------------------------------------------------------
C DINTYD computes interpolated values of the K-th derivative of the
C dependent variable vector Y, and stores it in DKY.
C This routine is called by DDEBDF with K = 0,1 and T = TOUT, but may
C also be called by the user for any K up to the current order.
C (see detailed instructions in LSODE usage documentation.)
C ----------------------------------------------------------------------
C The computed values in DKY are gotten by interpolation using the
C Nordsieck history array YH.  This array corresponds uniquely to a
C vector-valued polynomial of degree NQCUR or less, and DKY is set
C to the K-th derivative of this polynomial at T.
C The formula for DKY is..
C              Q
C  DKY(I)  =  Sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)
C             J=K
C where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.
C The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are
C communicated by common.  The above sum is done in reverse order.
C IFLAG is returned negative if either K or T is out of bounds.
C ----------------------------------------------------------------------
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DINTYD
C
      INTEGER I, IC, IER, IFLAG, IOWND, IOWNS, J, JB, JB2, JJ, JJ1,
     1      JP1, JSTART, K, KFLAG, L, MAXORD, METH, MITER, N, NFE,
     2      NJE, NQ, NQU, NST, NYH
      DOUBLE PRECISION C, DKY, EL0, H, HMIN, HMXI, HU, R, ROWND,
     1      ROWNS, S, T, TN, TP, UROUND, YH
      DIMENSION YH(NYH,*),DKY(*)
      COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,
     1                IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER,
     2                MAXORD,N,NQ,NST,NFE,NJE,NQU
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 130
C***FIRST EXECUTABLE STATEMENT  DINTYD
         IFLAG = 0
         IF (K .LT. 0 .OR. K .GT. NQ) GO TO 110
            TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
            IF ((T - TP)*(T - TN) .LE. 0.0D0) GO TO 10
               IFLAG = -2
C     .........EXIT
               GO TO 130
   10       CONTINUE
C
            S = (T - TN)/H
            IC = 1
            IF (K .EQ. 0) GO TO 30
               JJ1 = L - K
               DO 20 JJ = JJ1, NQ
                  IC = IC*JJ
   20          CONTINUE
   30       CONTINUE
            C = IC
            DO 40 I = 1, N
               DKY(I) = C*YH(I,L)
   40       CONTINUE
            IF (K .EQ. NQ) GO TO 90
               JB2 = NQ - K
               DO 80 JB = 1, JB2
                  J = NQ - JB
                  JP1 = J + 1
                  IC = 1
                  IF (K .EQ. 0) GO TO 60
                     JJ1 = JP1 - K
                     DO 50 JJ = JJ1, J
                        IC = IC*JJ
   50                CONTINUE
   60             CONTINUE
                  C = IC
                  DO 70 I = 1, N
                     DKY(I) = C*YH(I,JP1) + S*DKY(I)
   70             CONTINUE
   80          CONTINUE
C     .........EXIT
               IF (K .EQ. 0) GO TO 130
   90       CONTINUE
            R = H**(-K)
            DO 100 I = 1, N
               DKY(I) = R*DKY(I)
  100       CONTINUE
         GO TO 120
  110    CONTINUE
C
            IFLAG = -1
  120    CONTINUE
  130 CONTINUE
      RETURN
C     ----------------------- END OF SUBROUTINE DINTYD
C     -----------------------
      END
*DECK DPJAC
      SUBROUTINE DPJAC (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, DF,
     +   DJAC, RPAR, IPAR)
C***BEGIN PROLOGUE  DPJAC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (PJAC-S, DPJAC-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DPJAC sets up the iteration matrix (involving the Jacobian) for the
C   integration package DDEBDF.
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  DGBFA, DGEFA, DVNRMS
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920422  Changed DIMENSION statement.  (WRB)
C***END PROLOGUE  DPJAC
C
      INTEGER I, I1, I2, IER, II, IOWND, IOWNS, IPAR, IWM, J, J1,
     1      JJ, JSTART, KFLAG, L, LENP, MAXORD, MBA, MBAND,
     2      MEB1, MEBAND, METH, MITER, ML, ML3, MU, N, NEQ,
     3      NFE, NJE, NQ, NQU, NST, NYH
      DOUBLE PRECISION CON, DI, DVNRMS, EL0, EWT,
     1      FAC, FTEM, H, HL0, HMIN, HMXI, HU, R, R0, ROWND, ROWNS,
     2      RPAR, SAVF, SRUR, TN, UROUND, WM, Y, YH, YI, YJ, YJJ
      EXTERNAL DF, DJAC
      DIMENSION Y(*),YH(NYH,*),EWT(*),FTEM(*),SAVF(*),WM(*),IWM(*),
     1          RPAR(*),IPAR(*)
      COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,
     1                IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER,
     2                MAXORD,N,NQ,NST,NFE,NJE,NQU
C     ------------------------------------------------------------------
C      DPJAC IS CALLED BY DSTOD  TO COMPUTE AND PROCESS THE MATRIX
C      P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN.
C      HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE DJAC IF
C      MITER = 1 OR 4, OR BY FINITE DIFFERENCING IF MITER = 2, 3, OR 5.
C      IF MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED.
C      J IS STORED IN WM AND REPLACED BY P.  IF MITER .NE. 3, P IS THEN
C      SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION
C      OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE
C      BY DGEFA IF MITER = 1 OR 2, AND BY DGBFA IF MITER = 4 OR 5.
C
C      IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION
C      WITH DPJAC USES THE FOLLOWING..
C      Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY.
C      FTEM = WORK ARRAY OF LENGTH N (ACOR IN DSTOD ).
C      SAVF = ARRAY CONTAINING DF EVALUATED AT PREDICTED Y.
C      WM   = DOUBLE PRECISION WORK SPACE FOR MATRICES.  ON OUTPUT IT
C      CONTAINS THE
C             INVERSE DIAGONAL MATRIX IF MITER = 3 AND THE LU
C             DECOMPOSITION OF P IF MITER IS 1, 2 , 4, OR 5.
C             STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
C             WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
C             WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN
C             INCREMENTS.  WM(2) = H*EL0, SAVED FOR LATER USE IF MITER =
C             3.
C      IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING
C             AT IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS
C             THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER
C             IS 4 OR 5.
C      EL0  = EL(1) (INPUT).
C      IER  = OUTPUT ERROR FLAG,  = 0 IF NO TROUBLE, .NE. 0 IF
C             P MATRIX FOUND TO BE SINGULAR.
C      THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND,
C      MITER, N, NFE, AND NJE.
C-----------------------------------------------------------------------
C     BEGIN BLOCK PERMITTING ...EXITS TO 240
C        BEGIN BLOCK PERMITTING ...EXITS TO 220
C           BEGIN BLOCK PERMITTING ...EXITS TO 130
C              BEGIN BLOCK PERMITTING ...EXITS TO 70
C***FIRST EXECUTABLE STATEMENT  DPJAC
                  NJE = NJE + 1
                  HL0 = H*EL0
                  GO TO (10,40,90,140,170), MITER
C                 IF MITER = 1, CALL DJAC AND MULTIPLY BY SCALAR.
C                 -----------------------
   10             CONTINUE
                  LENP = N*N
                  DO 20 I = 1, LENP
                     WM(I+2) = 0.0D0
   20             CONTINUE
                  CALL DJAC(TN,Y,WM(3),N,RPAR,IPAR)
                  CON = -HL0
                  DO 30 I = 1, LENP
                     WM(I+2) = WM(I+2)*CON
   30             CONTINUE
C              ...EXIT
                  GO TO 70
C                 IF MITER = 2, MAKE N CALLS TO DF TO APPROXIMATE J.
C                 --------------------
   40             CONTINUE
                  FAC = DVNRMS(N,SAVF,EWT)
                  R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
                  IF (R0 .EQ. 0.0D0) R0 = 1.0D0
                  SRUR = WM(1)
                  J1 = 2
                  DO 60 J = 1, N
                     YJ = Y(J)
                     R = MAX(SRUR*ABS(YJ),R0*EWT(J))
                     Y(J) = Y(J) + R
                     FAC = -HL0/R
                     CALL DF(TN,Y,FTEM,RPAR,IPAR)
                     DO 50 I = 1, N
                        WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
   50                CONTINUE
                     Y(J) = YJ
                     J1 = J1 + N
   60             CONTINUE
                  NFE = NFE + N
   70          CONTINUE
C              ADD IDENTITY MATRIX.
C              -------------------------------------------------
               J = 3
               DO 80 I = 1, N
                  WM(J) = WM(J) + 1.0D0
                  J = J + (N + 1)
   80          CONTINUE
C              DO LU DECOMPOSITION ON P.
C              --------------------------------------------
               CALL DGEFA(WM(3),N,N,IWM(21),IER)
C     .........EXIT
               GO TO 240
C              IF MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND
C              P. ---------
   90          CONTINUE
               WM(2) = HL0
               IER = 0
               R = EL0*0.1D0
               DO 100 I = 1, N
                  Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
  100          CONTINUE
               CALL DF(TN,Y,WM(3),RPAR,IPAR)
               NFE = NFE + 1
               DO 120 I = 1, N
                  R0 = H*SAVF(I) - YH(I,2)
                  DI = 0.1D0*R0 - H*(WM(I+2) - SAVF(I))
                  WM(I+2) = 1.0D0
                  IF (ABS(R0) .LT. UROUND*EWT(I)) GO TO 110
C           .........EXIT
                     IF (ABS(DI) .EQ. 0.0D0) GO TO 130
                     WM(I+2) = 0.1D0*R0/DI
  110             CONTINUE
  120          CONTINUE
C     .........EXIT
               GO TO 240
  130       CONTINUE
            IER = -1
C     ......EXIT
            GO TO 240
C           IF MITER = 4, CALL DJAC AND MULTIPLY BY SCALAR.
C           -----------------------
  140       CONTINUE
            ML = IWM(1)
            MU = IWM(2)
            ML3 = 3
            MBAND = ML + MU + 1
            MEBAND = MBAND + ML
            LENP = MEBAND*N
            DO 150 I = 1, LENP
               WM(I+2) = 0.0D0
  150       CONTINUE
            CALL DJAC(TN,Y,WM(ML3),MEBAND,RPAR,IPAR)
            CON = -HL0
            DO 160 I = 1, LENP
               WM(I+2) = WM(I+2)*CON
  160       CONTINUE
C        ...EXIT
            GO TO 220
C           IF MITER = 5, MAKE MBAND CALLS TO DF TO APPROXIMATE J.
C           ----------------
  170       CONTINUE
            ML = IWM(1)
            MU = IWM(2)
            MBAND = ML + MU + 1
            MBA = MIN(MBAND,N)
            MEBAND = MBAND + ML
            MEB1 = MEBAND - 1
            SRUR = WM(1)
            FAC = DVNRMS(N,SAVF,EWT)
            R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
            IF (R0 .EQ. 0.0D0) R0 = 1.0D0
            DO 210 J = 1, MBA
               DO 180 I = J, N, MBAND
                  YI = Y(I)
                  R = MAX(SRUR*ABS(YI),R0*EWT(I))
                  Y(I) = Y(I) + R
  180          CONTINUE
               CALL DF(TN,Y,FTEM,RPAR,IPAR)
               DO 200 JJ = J, N, MBAND
                  Y(JJ) = YH(JJ,1)
                  YJJ = Y(JJ)
                  R = MAX(SRUR*ABS(YJJ),R0*EWT(JJ))
                  FAC = -HL0/R
                  I1 = MAX(JJ-MU,1)
                  I2 = MIN(JJ+ML,N)
                  II = JJ*MEB1 - ML + 2
                  DO 190 I = I1, I2
                     WM(II+I) = (FTEM(I) - SAVF(I))*FAC
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
            NFE = NFE + MBA
  220    CONTINUE
C        ADD IDENTITY MATRIX.
C        -------------------------------------------------
         II = MBAND + 2
         DO 230 I = 1, N
            WM(II) = WM(II) + 1.0D0
            II = II + MEBAND
  230    CONTINUE
C        DO LU DECOMPOSITION OF P.
C        --------------------------------------------
         CALL DGBFA(WM(3),MEBAND,N,ML,MU,IWM(21),IER)
  240 CONTINUE
      RETURN
C     ----------------------- END OF SUBROUTINE DPJAC
C     -----------------------
      END
*DECK DSLVS
      SUBROUTINE DSLVS (WM, IWM, X, TEM)
C***BEGIN PROLOGUE  DSLVS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SLVS-S, DSLVS-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DSLVS solves the linear system in the iteration scheme for the
C   integrator package DDEBDF.
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  DGBSL, DGESL
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920422  Changed DIMENSION statement.  (WRB)
C***END PROLOGUE  DSLVS
C
      INTEGER I, IER, IOWND, IOWNS, IWM, JSTART, KFLAG, L, MAXORD,
     1      MEBAND, METH, MITER, ML, MU, N, NFE, NJE, NQ, NQU, NST
      DOUBLE PRECISION DI, EL0, H, HL0, HMIN, HMXI, HU, PHL0,
     1      R, ROWND, ROWNS, TEM, TN, UROUND, WM, X
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
      COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,
     1                IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER,
     2                MAXORD,N,NQ,NST,NFE,NJE,NQU
C     ------------------------------------------------------------------
C      THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR SYSTEM ARISING
C      FROM A CHORD ITERATION.  IT IS CALLED BY DSTOD  IF MITER .NE. 0.
C      IF MITER IS 1 OR 2, IT CALLS DGESL TO ACCOMPLISH THIS.
C      IF MITER = 3 IT UPDATES THE COEFFICIENT H*EL0 IN THE DIAGONAL
C      MATRIX, AND THEN COMPUTES THE SOLUTION.
C      IF MITER IS 4 OR 5, IT CALLS DGBSL.
C      COMMUNICATION WITH DSLVS USES THE FOLLOWING VARIABLES..
C      WM  = DOUBLE PRECISION WORK SPACE CONTAINING THE INVERSE DIAGONAL
C      MATRIX IF MITER
C            IS 3 AND THE LU DECOMPOSITION OF THE MATRIX OTHERWISE.
C            STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
C            WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
C            WM(1) = SQRT(UROUND) (NOT USED HERE),
C            WM(2) = HL0, THE PREVIOUS VALUE OF H*EL0, USED IF MITER =
C            3.
C      IWM = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING
C            AT IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS
C            THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS
C            4 OR 5.
C      X   = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION
C            VECTOR ON OUTPUT, OF LENGTH N.
C      TEM = VECTOR OF WORK SPACE OF LENGTH N, NOT USED IN THIS VERSION.
C      IER = OUTPUT FLAG (IN COMMON).  IER = 0 IF NO TROUBLE OCCURRED.
C            IER = -1 IF A SINGULAR MATRIX AROSE WITH MITER = 3.
C      THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, MITER, AND N.
C-----------------------------------------------------------------------
C     BEGIN BLOCK PERMITTING ...EXITS TO 80
C        BEGIN BLOCK PERMITTING ...EXITS TO 60
C***FIRST EXECUTABLE STATEMENT  DSLVS
            IER = 0
            GO TO (10,10,20,70,70), MITER
   10       CONTINUE
            CALL DGESL(WM(3),N,N,IWM(21),X,0)
C     ......EXIT
            GO TO 80
C
   20       CONTINUE
            PHL0 = WM(2)
            HL0 = H*EL0
            WM(2) = HL0
            IF (HL0 .EQ. PHL0) GO TO 40
               R = HL0/PHL0
               DO 30 I = 1, N
                  DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
C        .........EXIT
                  IF (ABS(DI) .EQ. 0.0D0) GO TO 60
                  WM(I+2) = 1.0D0/DI
   30          CONTINUE
   40       CONTINUE
            DO 50 I = 1, N
               X(I) = WM(I+2)*X(I)
   50       CONTINUE
C     ......EXIT
            GO TO 80
   60    CONTINUE
         IER = -1
C     ...EXIT
         GO TO 80
C
   70    CONTINUE
         ML = IWM(1)
         MU = IWM(2)
         MEBAND = 2*ML + MU + 1
         CALL DGBSL(WM(3),MEBAND,N,ML,MU,IWM(21),X,0)
   80 CONTINUE
      RETURN
C     ----------------------- END OF SUBROUTINE DSLVS
C     -----------------------
      END
*DECK DGBFA
      SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
C***BEGIN PROLOGUE  DGBFA
C***PURPOSE  Factor a band matrix using Gaussian elimination.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2A2
C***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGBFA factors a double precision band matrix by elimination.
C
C     DGBFA is usually called by DGBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT.  N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT.  N .
C                More efficient if  ML .LE. MU .
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGBSL will divide by zero if
C                     called.  Use  RCOND  in DGBCO for a reliable
C                     indication of singularity.
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX(1, J-MU)
C                      I2 = MIN(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGBFA
      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      DOUBLE PRECISION ABD(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C***FIRST EXECUTABLE STATEMENT  DGBFA
      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
C
C        FIND L = PIVOT INDEX
C
         LM = MIN(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN(MAX(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGBSL
      SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
C***BEGIN PROLOGUE  DGBSL
C***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
C            the factors computed by DGBCO or DGBFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2A2
C***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGBSL solves the double precision band system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGBCO or DGBFA.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the output from DGBCO or DGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGBCO or DGBFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B , where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGBCO has set RCOND .GT. 0.0
C        or DGBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGBSL
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      DOUBLE PRECISION ABD(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
C***FIRST EXECUTABLE STATEMENT  DGBSL
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DGEFA
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
C***BEGIN PROLOGUE  DGEFA
C***PURPOSE  Factor a matrix using Gaussian elimination.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
C***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGESL
      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
C***BEGIN PROLOGUE  DGESL
C***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
C            factors computed by DGECO or DGEFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
