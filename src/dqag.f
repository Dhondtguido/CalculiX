c
c     CMLIB - public domain
c
      subroutine dqag(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     *    limit,lenw,last,iwork,work,phi,lambda1,zk0,Pup,Tup,rurd,xflow,
     *     kup)
c***begin prologue  dqag
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             integrand examinator, globally adaptive,
c             gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result)le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c            f      - double precision
c                     function subprogam defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            epsabs - double precision
c                     absolute accoracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            key    - integer
c                     key for choice of local integration rule
c                     a gauss-kronrod pair is used with
c                       7 - 15 points if key.lt.2,
c                      10 - 21 points if key = 2,
c                      15 - 31 points if key = 3,
c                      20 - 41 points if key = 4,
c                      25 - 51 points if key = 5,
c                      30 - 61 points if key.gt.5.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for result and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c                      error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yield no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulaties.
c                             if the position of a local difficulty can
c                             be determined (i.e.singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set
c                             to zero.
c                             except when lenw is invalid, iwork(1),
c                             work(limit*2+1) and work(limit*3+1) are
c                             set to zero, work(1) is set to a and
c                             work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end with
c                    ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdiviosion process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers to the error
c                    estimates over the subintervals, such that
c                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
c                    form a decreasing sequence with k = last if
c                    last.le.(limit/2+2), and k = limit+1-last otherwise
c
c            work  - double precision
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left end
c                    points of the subintervals in the partition of
c                     (a,b),
c                    work(limit+1), ..., work(limit+last) contain the
c                     right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last) contain
c                     the error estimates.
c
c***references  (none)
c***routines called  dqage,xerror
c***end prologue  dqag
      real*8 a,abserr,b,epsabs,epsrel,f,result,work,d1mach(4),
     *     phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup
      integer ier,iwork,key,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of lenw.
c
      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
c
c***first executable statement  dqag
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqage.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,
     *  ier,work(1),work(l1),work(l2),work(l3),iwork,last,phi,lambda1,
     *     zk0,Pup,Tup,rurd,xflow,kup)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
!      if(ier.ne.0) call xerror(26habnormal return from dqag ,26,ier,lvl)
      return
      end
      subroutine dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last,phi,lambda1,zk0,
     *     Pup,Tup,rurd,xflow,kup)
c***begin prologue  dqage
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             integrand examinator, globally adaptive,
c             gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral   i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            key    - integer
c                     key for choice of local integration rule
c                     a gauss-kronrod pair is used with
c                          7 - 15 points if key.lt.2,
c                         10 - 21 points if key = 2,
c                         15 - 31 points if key = 3,
c                         20 - 41 points if key = 4,
c                         25 - 51 points if key = 5,
c                         30 - 61 points if key.gt.5.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for result and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value
c                             of limit.
c                             however, if this yields no improvement it
c                             is rather advised to analyze the integrand
c                             in order to determine the integration
c                             difficulties. if the position of a local
c                             difficulty can be determined(e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling the integrator on the
c                             subranges. if possible, an appropriate
c                             special-purpose integrator should be used
c                             which is designed for handling the type of
c                             difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1) ,
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the left
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            blist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the right
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            rlist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the
c                      integral approximations on the subintervals
c
c            elist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the moduli of the
c                      absolute error estimates on the subintervals
c
c            iord    - integer
c                      vector of dimension at least limit, the first k
c                      elements of which are pointers to the
c                      error estimates over the subintervals,
c                      such that elist(iord(1)), ...,
c                      elist(iord(k)) form a decreasing sequence,
c                      with k = last if last.le.(limit/2+2), and
c                      k = limit+1-last otherwise
c
c            last    - integer
c                      number of subintervals actually produced in the
c                      subdivision process
c
c***references  (none)
c***routines called  d1mach,dqk15,dqk21,dqk31,
c                    dqk41,dqk51,dqk61,dqpsrt
c***end prologue  dqage
c
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,d1mach(4),elist,
     * epmach,epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  rlist(limit)
c
      external f


      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                      (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach  is the largest relative spacing.
c           uflow  is the smallest positive magnitude.
c
c***first executable statement  dqage
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      if(epsabs.le.0.0d+00.and.
     *  epsrel.lt.max(0.5d+02*epmach,0.5d-28)) ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      keyf = key
      if(key.le.0) keyf = 1
      if(key.ge.7) keyf = 6
      neval = 0
      if(keyf.eq.1) call dqk15(f,a,b,result,abserr,defabs,resabs,
     &     phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.2) call dqk21(f,a,b,result,abserr,defabs,resabs,
     &    phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.3) call dqk31(f,a,b,result,abserr,defabs,resabs,
     &     phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.4) call dqk41(f,a,b,result,abserr,defabs,resabs,
     &     phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.5) call dqk51(f,a,b,result,abserr,defabs,resabs,
     &     phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.6) call dqk61(f,a,b,result,abserr,defabs,resabs,
     &     phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
c
c           test on accuracy.
c
      errbnd = max(epsabs,epsrel*dabs(result))
      if(abserr.le.0.5d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs)
     *  .or.abserr.eq.0.0d+00) go to 60
c
c           initialization
c           --------------
c
c
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
c
c           main do-loop
c           ------------
c
      do 30 last = 2,limit
c
c           bisect the subinterval with the largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf.eq.1) call dqk15(f,a1,b1,area1,error1,resabs,defab1,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.2) call dqk21(f,a1,b1,area1,error1,resabs,defab1,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.3) call dqk31(f,a1,b1,area1,error1,resabs,defab1,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.4) call dqk41(f,a1,b1,area1,error1,resabs,defab1,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.5) call dqk51(f,a1,b1,area1,error1,resabs,defab1,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.6) call dqk61(f,a1,b1,area1,error1,resabs,defab1,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.1) call dqk15(f,a2,b2,area2,error2,resabs,defab2,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.2) call dqk21(f,a2,b2,area2,error2,resabs,defab2,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.3) call dqk31(f,a2,b2,area2,error2,resabs,defab2,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.4) call dqk41(f,a2,b2,area2,error2,resabs,defab2,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.5) call dqk51(f,a2,b2,area2,error2,resabs,defab2,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.6) call dqk61(f,a2,b2,area2,error2,resabs,defab2,
     *      phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 5
        if(dabs(rlist(maxerr)-area12).le.0.1d-04*dabs(area12)
     *  .and.erro12.ge.0.99d+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
    5   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 8
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
c
c           set error flag in the case that the number of subintervals
c           equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(max(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*
     *  epmach)*(dabs(a2)+0.1d+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
    8   if(error2.gt.error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
   10   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with the largest error estimate (to be bisected next).
c
   20   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax,phi,
     *       lambda1,zk0,Pup,Tup,rurd,xflow,kup)
c ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 40
   30 continue
c
c           compute final result.
c           ---------------------
c
   40 result = 0.0d+00
      do 50 k=1,last
        result = result+rlist(k)
   50 continue
      abserr = errsum
   60 if(keyf.ne.1) neval = (10*keyf+1)*(2*neval+1)
      if(keyf.eq.1) neval = 30*neval+15
  999 return
      end
      subroutine dqk15(f,a,b,result,abserr,resabs,resasc,phi,lambda1,
     *     zk0,Pup,Tup,rurd,xflow,kup)
c***begin prologue  dqk15
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  15-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the7-point gauss rule(resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,
     *  resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,
     *     zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1 
      external f
c
      dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)

      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21


c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 15-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point kronrod rule
c
c           wg     - weights of the 7-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.1294849661 6886969327 0611432679 082 d0 /
      data wg  (  2) / 0.2797053914 8927666790 1467771423 780 d0 /
      data wg  (  3) / 0.3818300505 0511894495 0369775488 975 d0 /
      data wg  (  4) / 0.4179591836 7346938775 5102040816 327 d0 /
c
      data xgk (  1) / 0.9914553711 2081263920 6854697526 329 d0 /
      data xgk (  2) / 0.9491079123 4275852452 6189684047 851 d0 /
      data xgk (  3) / 0.8648644233 5976907278 9712788640 926 d0 /
      data xgk (  4) / 0.7415311855 9939443986 3864773280 788 d0 /
      data xgk (  5) / 0.5860872354 6769113029 4144838258 730 d0 /
      data xgk (  6) / 0.4058451513 7739716690 6606412076 961 d0 /
      data xgk (  7) / 0.2077849550 0789846760 0689403773 245 d0 /
      data xgk (  8) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0229353220 1052922496 3732008058 970 d0 /
      data wgk (  2) / 0.0630920926 2997855329 0700663189 204 d0 /
      data wgk (  3) / 0.1047900103 2225018383 9876322541 518 d0 /
      data wgk (  4) / 0.1406532597 1552591874 5189590510 238 d0 /
      data wgk (  5) / 0.1690047266 3926790282 6583426598 550 d0 /
      data wgk (  6) / 0.1903505780 6478540991 3256402421 014 d0 /
      data wgk (  7) / 0.2044329400 7529889241 4161999234 649 d0 /
      data wgk (  8) / 0.2094821410 8472782801 2999174891 714 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 15-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      resabs = dabs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk21(f,a,b,result,abserr,resabs,resasc,phi,lambda1,
     &     zk0,Pup,Tup,rurd,xflow,kup)
c***begin prologue  dqk21
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  21-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 21-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 10-point gauss rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk21
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,
     *  resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,
     *  zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
      
      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 21-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 10-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 10-point gauss rule
c
c           wgk    - weights of the 21-point kronrod rule
c
c           wg     - weights of the 10-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0666713443 0868813759 3568809893 332 d0 /
      data wg  (  2) / 0.1494513491 5058059314 5776339657 697 d0 /
      data wg  (  3) / 0.2190863625 1598204399 5534934228 163 d0 /
      data wg  (  4) / 0.2692667193 0999635509 1226921569 469 d0 /
      data wg  (  5) / 0.2955242247 1475287017 3892994651 338 d0 /
c
      data xgk (  1) / 0.9956571630 2580808073 5527280689 003 d0 /
      data xgk (  2) / 0.9739065285 1717172007 7964012084 452 d0 /
      data xgk (  3) / 0.9301574913 5570822600 1207180059 508 d0 /
      data xgk (  4) / 0.8650633666 8898451073 2096688423 493 d0 /
      data xgk (  5) / 0.7808177265 8641689706 3717578345 042 d0 /
      data xgk (  6) / 0.6794095682 9902440623 4327365114 874 d0 /
      data xgk (  7) / 0.5627571346 6860468333 9000099272 694 d0 /
      data xgk (  8) / 0.4333953941 2924719079 9265943165 784 d0 /
      data xgk (  9) / 0.2943928627 0146019813 1126603103 866 d0 /
      data xgk ( 10) / 0.1488743389 8163121088 4826001129 720 d0 /
      data xgk ( 11) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0116946388 6737187427 8064396062 192 d0 /
      data wgk (  2) / 0.0325581623 0796472747 8818972459 390 d0 /
      data wgk (  3) / 0.0547558965 7435199603 1381300244 580 d0 /
      data wgk (  4) / 0.0750396748 1091995276 7043140916 190 d0 /
      data wgk (  5) / 0.0931254545 8369760553 5065465083 366 d0 /
      data wgk (  6) / 0.1093871588 0229764189 9210590325 805 d0 /
      data wgk (  7) / 0.1234919762 6206585107 7958109831 074 d0 /
      data wgk (  8) / 0.1347092173 1147332592 8054001771 707 d0 /
      data wgk (  9) / 0.1427759385 7706008079 7094273138 717 d0 /
      data wgk ( 10) / 0.1477391049 0133849137 4841515972 068 d0 /
      data wgk ( 11) / 0.1494455540 0291690566 4936468389 821 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 10-point gauss formula
c           resk   - result of the 21-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk21
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 21-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0d+00
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resk = wgk(11)*fc
      resabs = dabs(resk)
      do 10 j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(11)*dabs(fc-reskh)
      do 20 j=1,10
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk31(f,a,b,result,abserr,resabs,resasc,phi,lambda1,
     *     zk0,Pup,Tup,rurd,xflow,kup)
c***begin prologue  dqk31
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  31-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 31-point
c                       gauss-kronrod rule (resk), obtained by optimal
c                       addition of abscissae to the 15-point gauss
c                       rule (resg).
c
c              abserr - double precison
c                       estimate of the modulus of the modulus,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk31
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,
     *  resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,
     *  zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)

      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 31-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 15-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 15-point gauss rule
c
c           wgk    - weights of the 31-point kronrod rule
c
c           wg     - weights of the 15-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0307532419 9611726835 4628393577 204 d0 /
      data wg  (  2) / 0.0703660474 8810812470 9267416450 667 d0 /
      data wg  (  3) / 0.1071592204 6717193501 1869546685 869 d0 /
      data wg  (  4) / 0.1395706779 2615431444 7804794511 028 d0 /
      data wg  (  5) / 0.1662692058 1699393355 3200860481 209 d0 /
      data wg  (  6) / 0.1861610000 1556221102 6800561866 423 d0 /
      data wg  (  7) / 0.1984314853 2711157645 6118326443 839 d0 /
      data wg  (  8) / 0.2025782419 2556127288 0620199967 519 d0 /
c
      data xgk (  1) / 0.9980022986 9339706028 5172840152 271 d0 /
      data xgk (  2) / 0.9879925180 2048542848 9565718586 613 d0 /
      data xgk (  3) / 0.9677390756 7913913425 7347978784 337 d0 /
      data xgk (  4) / 0.9372733924 0070590430 7758947710 209 d0 /
      data xgk (  5) / 0.8972645323 4408190088 2509656454 496 d0 /
      data xgk (  6) / 0.8482065834 1042721620 0648320774 217 d0 /
      data xgk (  7) / 0.7904185014 4246593296 7649294817 947 d0 /
      data xgk (  8) / 0.7244177313 6017004741 6186054613 938 d0 /
      data xgk (  9) / 0.6509967412 9741697053 3735895313 275 d0 /
      data xgk ( 10) / 0.5709721726 0853884753 7226737253 911 d0 /
      data xgk ( 11) / 0.4850818636 4023968069 3655740232 351 d0 /
      data xgk ( 12) / 0.3941513470 7756336989 7207370981 045 d0 /
      data xgk ( 13) / 0.2991800071 5316881216 6780024266 389 d0 /
      data xgk ( 14) / 0.2011940939 9743452230 0628303394 596 d0 /
      data xgk ( 15) / 0.1011420669 1871749902 7074231447 392 d0 /
      data xgk ( 16) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0053774798 7292334898 7792051430 128 d0 /
      data wgk (  2) / 0.0150079473 2931612253 8374763075 807 d0 /
      data wgk (  3) / 0.0254608473 2671532018 6874001019 653 d0 /
      data wgk (  4) / 0.0353463607 9137584622 2037948478 360 d0 /
      data wgk (  5) / 0.0445897513 2476487660 8227299373 280 d0 /
      data wgk (  6) / 0.0534815246 9092808726 5343147239 430 d0 /
      data wgk (  7) / 0.0620095678 0067064028 5139230960 803 d0 /
      data wgk (  8) / 0.0698541213 1872825870 9520077099 147 d0 /
      data wgk (  9) / 0.0768496807 5772037889 4432777482 659 d0 /
      data wgk ( 10) / 0.0830805028 2313302103 8289247286 104 d0 /
      data wgk ( 11) / 0.0885644430 5621177064 7275443693 774 d0 /
      data wgk ( 12) / 0.0931265981 7082532122 5486872747 346 d0 /
      data wgk ( 13) / 0.0966427269 8362367850 5179907627 589 d0 /
      data wgk ( 14) / 0.0991735987 2179195933 2393173484 603 d0 /
      data wgk ( 15) / 0.1007698455 2387559504 4946662617 570 d0 /
      data wgk ( 16) / 0.1013300070 1479154901 7374792767 493 d0 /
c
c
c           list of major variables
c           -----------------------
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 15-point gauss formula
c           resk   - result of the 31-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c***first executable statement  dqk31
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 31-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      resabs = dabs(resk)
      do 10 j=1,7
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,8
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(16)*dabs(fc-reskh)
      do 20 j=1,15
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk41(f,a,b,result,abserr,resabs,resasc,phi,lambda1,
     *     zk0,Pup,Tup,rurd,xflow,kup)
c***begin prologue  dqk41
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  41-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 41-point
c                       gauss-kronrod rule (resk) obtained by optimal
c                       addition of abscissae to the 20-point gauss
c                       rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integal of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk41
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,
     *  resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,
     *  zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)
      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 41-point gauss-kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 20-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 20-point gauss rule
c
c           wgk    - weights of the 41-point gauss-kronrod rule
c
c           wg     - weights of the 20-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0176140071 3915211831 1861962351 853 d0 /
      data wg  (  2) / 0.0406014298 0038694133 1039952274 932 d0 /
      data wg  (  3) / 0.0626720483 3410906356 9506535187 042 d0 /
      data wg  (  4) / 0.0832767415 7670474872 4758143222 046 d0 /
      data wg  (  5) / 0.1019301198 1724043503 6750135480 350 d0 /
      data wg  (  6) / 0.1181945319 6151841731 2377377711 382 d0 /
      data wg  (  7) / 0.1316886384 4917662689 8494499748 163 d0 /
      data wg  (  8) / 0.1420961093 1838205132 9298325067 165 d0 /
      data wg  (  9) / 0.1491729864 7260374678 7828737001 969 d0 /
      data wg  ( 10) / 0.1527533871 3072585069 8084331955 098 d0 /
c
      data xgk (  1) / 0.9988590315 8827766383 8315576545 863 d0 /
      data xgk (  2) / 0.9931285991 8509492478 6122388471 320 d0 /
      data xgk (  3) / 0.9815078774 5025025919 3342994720 217 d0 /
      data xgk (  4) / 0.9639719272 7791379126 7666131197 277 d0 /
      data xgk (  5) / 0.9408226338 3175475351 9982722212 443 d0 /
      data xgk (  6) / 0.9122344282 5132590586 7752441203 298 d0 /
      data xgk (  7) / 0.8782768112 5228197607 7442995113 078 d0 /
      data xgk (  8) / 0.8391169718 2221882339 4529061701 521 d0 /
      data xgk (  9) / 0.7950414288 3755119835 0638833272 788 d0 /
      data xgk ( 10) / 0.7463319064 6015079261 4305070355 642 d0 /
      data xgk ( 11) / 0.6932376563 3475138480 5490711845 932 d0 /
      data xgk ( 12) / 0.6360536807 2651502545 2836696226 286 d0 /
      data xgk ( 13) / 0.5751404468 1971031534 2946036586 425 d0 /
      data xgk ( 14) / 0.5108670019 5082709800 4364050955 251 d0 /
      data xgk ( 15) / 0.4435931752 3872510319 9992213492 640 d0 /
      data xgk ( 16) / 0.3737060887 1541956067 2548177024 927 d0 /
      data xgk ( 17) / 0.3016278681 1491300432 0555356858 592 d0 /
      data xgk ( 18) / 0.2277858511 4164507808 0496195368 575 d0 /
      data xgk ( 19) / 0.1526054652 4092267550 5220241022 678 d0 /
      data xgk ( 20) / 0.0765265211 3349733375 4640409398 838 d0 /
      data xgk ( 21) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0030735837 1852053150 1218293246 031 d0 /
      data wgk (  2) / 0.0086002698 5564294219 8661787950 102 d0 /
      data wgk (  3) / 0.0146261692 5697125298 3787960308 868 d0 /
      data wgk (  4) / 0.0203883734 6126652359 8010231432 755 d0 /
      data wgk (  5) / 0.0258821336 0495115883 4505067096 153 d0 /
      data wgk (  6) / 0.0312873067 7703279895 8543119323 801 d0 /
      data wgk (  7) / 0.0366001697 5820079803 0557240707 211 d0 /
      data wgk (  8) / 0.0416688733 2797368626 3788305936 895 d0 /
      data wgk (  9) / 0.0464348218 6749767472 0231880926 108 d0 /
      data wgk ( 10) / 0.0509445739 2372869193 2707670050 345 d0 /
      data wgk ( 11) / 0.0551951053 4828599474 4832372419 777 d0 /
      data wgk ( 12) / 0.0591114008 8063957237 4967220648 594 d0 /
      data wgk ( 13) / 0.0626532375 5478116802 5870122174 255 d0 /
      data wgk ( 14) / 0.0658345971 3361842211 1563556969 398 d0 /
      data wgk ( 15) / 0.0686486729 2852161934 5623411885 368 d0 /
      data wgk ( 16) / 0.0710544235 5344406830 5790361723 210 d0 /
      data wgk ( 17) / 0.0730306903 3278666749 5189417658 913 d0 /
      data wgk ( 18) / 0.0745828754 0049918898 6581418362 488 d0 /
      data wgk ( 19) / 0.0757044976 8455667465 9542775376 617 d0 /
      data wgk ( 20) / 0.0763778676 7208073670 5502835038 061 d0 /
      data wgk ( 21) / 0.0766007119 1799965644 5049901530 102 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 20-point gauss formula
c           resk   - result of the 41-point kronrod formula
c           reskh  - approximation to mean value of f over (a,b), i.e.
c                    to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk41
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 41-point gauss-kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0d+00
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resk = wgk(21)*fc
      resabs = dabs(resk)
      do 10 j=1,10
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,10
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(21)*dabs(fc-reskh)
      do 20 j=1,20
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk51(f,a,b,result,abserr,resabs,resasc,phi,lambda1,
     *     zk0,Pup,Tup,rurd,xflow,kup)
c***begin prologue  dqk51
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  51-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subroutine defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 51-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 25-point gauss rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk51
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,
     *     resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,
     *     zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(25),fv2(25),xgk(26),wgk(26),wg(13)
      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 51-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 25-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 25-point gauss rule
c
c           wgk    - weights of the 51-point kronrod rule
c
c           wg     - weights of the 25-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0113937985 0102628794 7902964113 235 d0 /
      data wg  (  2) / 0.0263549866 1503213726 1901815295 299 d0 /
      data wg  (  3) / 0.0409391567 0130631265 5623487711 646 d0 /
      data wg  (  4) / 0.0549046959 7583519192 5936891540 473 d0 /
      data wg  (  5) / 0.0680383338 1235691720 7187185656 708 d0 /
      data wg  (  6) / 0.0801407003 3500101801 3234959669 111 d0 /
      data wg  (  7) / 0.0910282619 8296364981 1497220702 892 d0 /
      data wg  (  8) / 0.1005359490 6705064420 2206890392 686 d0 /
      data wg  (  9) / 0.1085196244 7426365311 6093957050 117 d0 /
      data wg  ( 10) / 0.1148582591 4571164833 9325545869 556 d0 /
      data wg  ( 11) / 0.1194557635 3578477222 8178126512 901 d0 /
      data wg  ( 12) / 0.1222424429 9031004168 8959518945 852 d0 /
      data wg  ( 13) / 0.1231760537 2671545120 3902873079 050 d0 /
c
      data xgk (  1) / 0.9992621049 9260983419 3457486540 341 d0 /
      data xgk (  2) / 0.9955569697 9049809790 8784946893 902 d0 /
      data xgk (  3) / 0.9880357945 3407724763 7331014577 406 d0 /
      data xgk (  4) / 0.9766639214 5951751149 8315386479 594 d0 /
      data xgk (  5) / 0.9616149864 2584251241 8130033660 167 d0 /
      data xgk (  6) / 0.9429745712 2897433941 4011169658 471 d0 /
      data xgk (  7) / 0.9207471152 8170156174 6346084546 331 d0 /
      data xgk (  8) / 0.8949919978 7827536885 1042006782 805 d0 /
      data xgk (  9) / 0.8658470652 9327559544 8996969588 340 d0 /
      data xgk ( 10) / 0.8334426287 6083400142 1021108693 570 d0 /
      data xgk ( 11) / 0.7978737979 9850005941 0410904994 307 d0 /
      data xgk ( 12) / 0.7592592630 3735763057 7282865204 361 d0 /
      data xgk ( 13) / 0.7177664068 1308438818 6654079773 298 d0 /
      data xgk ( 14) / 0.6735663684 7346836448 5120633247 622 d0 /
      data xgk ( 15) / 0.6268100990 1031741278 8122681624 518 d0 /
      data xgk ( 16) / 0.5776629302 4122296772 3689841612 654 d0 /
      data xgk ( 17) / 0.5263252843 3471918259 9623778158 010 d0 /
      data xgk ( 18) / 0.4730027314 4571496052 2182115009 192 d0 /
      data xgk ( 19) / 0.4178853821 9303774885 1814394594 572 d0 /
      data xgk ( 20) / 0.3611723058 0938783773 5821730127 641 d0 /
      data xgk ( 21) / 0.3030895389 3110783016 7478909980 339 d0 /
      data xgk ( 22) / 0.2438668837 2098843204 5190362797 452 d0 /
      data xgk ( 23) / 0.1837189394 2104889201 5969888759 528 d0 /
      data xgk ( 24) / 0.1228646926 1071039638 7359818808 037 d0 /
      data xgk ( 25) / 0.0615444830 0568507888 6546392366 797 d0 /
      data xgk ( 26) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0019873838 9233031592 6507851882 843 d0 /
      data wgk (  2) / 0.0055619321 3535671375 8040236901 066 d0 /
      data wgk (  3) / 0.0094739733 8617415160 7207710523 655 d0 /
      data wgk (  4) / 0.0132362291 9557167481 3656405846 976 d0 /
      data wgk (  5) / 0.0168478177 0912829823 1516667536 336 d0 /
      data wgk (  6) / 0.0204353711 4588283545 6568292235 939 d0 /
      data wgk (  7) / 0.0240099456 0695321622 0092489164 881 d0 /
      data wgk (  8) / 0.0274753175 8785173780 2948455517 811 d0 /
      data wgk (  9) / 0.0307923001 6738748889 1109020215 229 d0 /
      data wgk ( 10) / 0.0340021302 7432933783 6748795229 551 d0 /
      data wgk ( 11) / 0.0371162714 8341554356 0330625367 620 d0 /
      data wgk ( 12) / 0.0400838255 0403238207 4839284467 076 d0 /
      data wgk ( 13) / 0.0428728450 2017004947 6895792439 495 d0 /
      data wgk ( 14) / 0.0455029130 4992178890 9870584752 660 d0 /
      data wgk ( 15) / 0.0479825371 3883671390 6392255756 915 d0 /
      data wgk ( 16) / 0.0502776790 8071567196 3325259433 440 d0 /
      data wgk ( 17) / 0.0523628858 0640747586 4366712137 873 d0 /
      data wgk ( 18) / 0.0542511298 8854549014 4543370459 876 d0 /
      data wgk ( 19) / 0.0559508112 2041231730 8240686382 747 d0 /
      data wgk ( 20) / 0.0574371163 6156783285 3582693939 506 d0 /
      data wgk ( 21) / 0.0586896800 2239420796 1974175856 788 d0 /
      data wgk ( 22) / 0.0597203403 2417405997 9099291932 562 d0 /
      data wgk ( 23) / 0.0605394553 7604586294 5360267517 565 d0 /
      data wgk ( 24) / 0.0611285097 1705304830 5859030416 293 d0 /
      data wgk ( 25) / 0.0614711898 7142531666 1544131965 264 d0 /
c       note: wgk (26) was calculated from the values of wgk(1..25)
      data wgk ( 26) / 0.0615808180 6783293507 8759824240 066 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 25-point gauss formula
c           resk   - result of the 51-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk51
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 51-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resg = wg(13)*fc
      resk = wgk(26)*fc
      resabs = dabs(resk)
      do 10 j=1,12
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,13
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(26)*dabs(fc-reskh)
      do 20 j=1,25
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk61(f,a,b,result,abserr,resabs,resasc,phi,lambda1,
     *     zk0,Pup,Tup,rurd,xflow,kup)
c***begin prologue  dqk61
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  61-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of dabs(f) over (a,b)
c***description
c
c        integration rule
c        standard fortran subroutine
c        double precision version
c
c
c        parameters
c         on entry
c           f      - double precision
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to be
c                    declared e x t e r n a l in the calling program.
c
c           a      - double precision
c                    lower limit of integration
c
c           b      - double precision
c                    upper limit of integration
c
c         on return
c           result - double precision
c                    approximation to the integral i
c                    result is computed by applying the 61-point
c                    kronrod rule (resk) obtained by optimal addition of
c                    abscissae to the 30-point gauss rule (resg).
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed dabs(i-result)
c
c           resabs - double precision
c                    approximation to the integral j
c
c           resasc - double precision
c                    approximation to the integral of dabs(f-i/(b-a))
c
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk61
c
      double precision a,dabsc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,
     *  resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,
     *  zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)
      d1mach(1)=1E21
      d1mach(2)=0
      d1mach(3)=0
      d1mach(4)=1E-21
c
c           the abscissae and weights are given for the
c           interval (-1,1). because of symmetry only the positive
c           abscissae and their corresponding weights are given.
c
c           xgk   - abscissae of the 61-point kronrod rule
c                   xgk(2), xgk(4)  ... abscissae of the 30-point
c                   gauss rule
c                   xgk(1), xgk(3)  ... optimally added abscissae
c                   to the 30-point gauss rule
c
c           wgk   - weights of the 61-point kronrod rule
c
c           wg    - weigths of the 30-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0079681924 9616660561 5465883474 674 d0 /
      data wg  (  2) / 0.0184664683 1109095914 2302131912 047 d0 /
      data wg  (  3) / 0.0287847078 8332336934 9719179611 292 d0 /
      data wg  (  4) / 0.0387991925 6962704959 6801936446 348 d0 /
      data wg  (  5) / 0.0484026728 3059405290 2938140422 808 d0 /
      data wg  (  6) / 0.0574931562 1761906648 1721689402 056 d0 /
      data wg  (  7) / 0.0659742298 8218049512 8128515115 962 d0 /
      data wg  (  8) / 0.0737559747 3770520626 8243850022 191 d0 /
      data wg  (  9) / 0.0807558952 2942021535 4694938460 530 d0 /
      data wg  ( 10) / 0.0868997872 0108297980 2387530715 126 d0 /
      data wg  ( 11) / 0.0921225222 3778612871 7632707087 619 d0 /
      data wg  ( 12) / 0.0963687371 7464425963 9468626351 810 d0 /
      data wg  ( 13) / 0.0995934205 8679526706 2780282103 569 d0 /
      data wg  ( 14) / 0.1017623897 4840550459 6428952168 554 d0 /
      data wg  ( 15) / 0.1028526528 9355884034 1285636705 415 d0 /
c
      data xgk (  1) / 0.9994844100 5049063757 1325895705 811 d0 /
      data xgk (  2) / 0.9968934840 7464954027 1630050918 695 d0 /
      data xgk (  3) / 0.9916309968 7040459485 8628366109 486 d0 /
      data xgk (  4) / 0.9836681232 7974720997 0032581605 663 d0 /
      data xgk (  5) / 0.9731163225 0112626837 4693868423 707 d0 /
      data xgk (  6) / 0.9600218649 6830751221 6871025581 798 d0 /
      data xgk (  7) / 0.9443744447 4855997941 5831324037 439 d0 /
      data xgk (  8) / 0.9262000474 2927432587 9324277080 474 d0 /
      data xgk (  9) / 0.9055733076 9990779854 6522558925 958 d0 /
      data xgk ( 10) / 0.8825605357 9205268154 3116462530 226 d0 /
      data xgk ( 11) / 0.8572052335 4606109895 8658510658 944 d0 /
      data xgk ( 12) / 0.8295657623 8276839744 2898119732 502 d0 /
      data xgk ( 13) / 0.7997278358 2183908301 3668942322 683 d0 /
      data xgk ( 14) / 0.7677774321 0482619491 7977340974 503 d0 /
      data xgk ( 15) / 0.7337900624 5322680472 6171131369 528 d0 /
      data xgk ( 16) / 0.6978504947 9331579693 2292388026 640 d0 /
      data xgk ( 17) / 0.6600610641 2662696137 0053668149 271 d0 /
      data xgk ( 18) / 0.6205261829 8924286114 0477556431 189 d0 /
      data xgk ( 19) / 0.5793452358 2636169175 6024932172 540 d0 /
      data xgk ( 20) / 0.5366241481 4201989926 4169793311 073 d0 /
      data xgk ( 21) / 0.4924804678 6177857499 3693061207 709 d0 /
      data xgk ( 22) / 0.4470337695 3808917678 0609900322 854 d0 /
      data xgk ( 23) / 0.4004012548 3039439253 5476211542 661 d0 /
      data xgk ( 24) / 0.3527047255 3087811347 1037207089 374 d0 /
      data xgk ( 25) / 0.3040732022 7362507737 2677107199 257 d0 /
      data xgk ( 26) / 0.2546369261 6788984643 9805129817 805 d0 /
      data xgk ( 27) / 0.2045251166 8230989143 8957671002 025 d0 /
      data xgk ( 28) / 0.1538699136 0858354696 3794672743 256 d0 /
      data xgk ( 29) / 0.1028069379 6673703014 7096751318 001 d0 /
      data xgk ( 30) / 0.0514718425 5531769583 3025213166 723 d0 /
      data xgk ( 31) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0013890136 9867700762 4551591226 760 d0 /
      data wgk (  2) / 0.0038904611 2709988405 1267201844 516 d0 /
      data wgk (  3) / 0.0066307039 1593129217 3319826369 750 d0 /
      data wgk (  4) / 0.0092732796 5951776342 8441146892 024 d0 /
      data wgk (  5) / 0.0118230152 5349634174 2232898853 251 d0 /
      data wgk (  6) / 0.0143697295 0704580481 2451432443 580 d0 /
      data wgk (  7) / 0.0169208891 8905327262 7572289420 322 d0 /
      data wgk (  8) / 0.0194141411 9394238117 3408951050 128 d0 /
      data wgk (  9) / 0.0218280358 2160919229 7167485738 339 d0 /
      data wgk ( 10) / 0.0241911620 7808060136 5686370725 232 d0 /
      data wgk ( 11) / 0.0265099548 8233310161 0601709335 075 d0 /
      data wgk ( 12) / 0.0287540487 6504129284 3978785354 334 d0 /
      data wgk ( 13) / 0.0309072575 6238776247 2884252943 092 d0 /
      data wgk ( 14) / 0.0329814470 5748372603 1814191016 854 d0 /
      data wgk ( 15) / 0.0349793380 2806002413 7499670731 468 d0 /
      data wgk ( 16) / 0.0368823646 5182122922 3911065617 136 d0 /
      data wgk ( 17) / 0.0386789456 2472759295 0348651532 281 d0 /
      data wgk ( 18) / 0.0403745389 5153595911 1995279752 468 d0 /
      data wgk ( 19) / 0.0419698102 1516424614 7147541285 970 d0 /
      data wgk ( 20) / 0.0434525397 0135606931 6831728117 073 d0 /
      data wgk ( 21) / 0.0448148001 3316266319 2355551616 723 d0 /
      data wgk ( 22) / 0.0460592382 7100698811 6271735559 374 d0 /
      data wgk ( 23) / 0.0471855465 6929915394 5261478181 099 d0 /
      data wgk ( 24) / 0.0481858617 5708712914 0779492298 305 d0 /
      data wgk ( 25) / 0.0490554345 5502977888 7528165367 238 d0 /
      data wgk ( 26) / 0.0497956834 2707420635 7811569379 942 d0 /
      data wgk ( 27) / 0.0504059214 0278234684 0893085653 585 d0 /
      data wgk ( 28) / 0.0508817958 9874960649 2297473049 805 d0 /
      data wgk ( 29) / 0.0512215478 4925877217 0656282604 944 d0 /
      data wgk ( 30) / 0.0514261285 3745902593 3862879215 781 d0 /
      data wgk ( 31) / 0.0514947294 2945156755 8340433647 099 d0 /
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           dabsc  - abscissa
c           fval*  - function value
c           resg   - result of the 30-point gauss rule
c           resk   - result of the 61-point kronrod rule
c           reskh  - approximation to the mean value of f
c                    over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(b+a)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 61-point kronrod approximation to the
c           integral, and estimate the absolute error.
c
c***first executable statement  dqk61
      resg = 0.0d+00
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resk = wgk(31)*fc
      resabs = dabs(resk)
      do 10 j=1,15
        jtw = j*2
        dabsc = hlgth*xgk(jtw)
        fval1 = f(centr-dabsc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+dabsc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j=1,15
        jtwm1 = j*2-1
        dabsc = hlgth*xgk(jtwm1)
        fval1 = f(centr-dabsc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+dabsc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
  15    continue
      reskh = resk*0.5d+00
      resasc = wgk(31)*dabs(fc-reskh)
      do 20 j=1,30
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax,
     *     phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
c***begin prologue  dqpsrt
c***refer to  dqage,dqagie,dqagpe,dqawse
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  sequential sorting
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine maintains the descending ordering in the
c            list of the local error estimated resulting from the
c            interval subdivision process. at each call two error
c            estimates are inserted using the sequential search
c            method, top-down for the largest error estimate and
c            bottom-up for the smallest error estimate.
c***description
c
c           ordering routine
c           standard fortran subroutine
c           double precision version
c
c           parameters (meaning at output)
c              limit  - integer
c                       maximum number of error estimates the list
c                       can contain
c
c              last   - integer
c                       number of error estimates currently in the list
c
c              maxerr - integer
c                       maxerr points to the nrmax-th largest error
c                       estimate currently in the list
c
c              ermax  - double precision
c                       nrmax-th largest error estimate
c                       ermax = elist(maxerr)
c
c              elist  - double precision
c                       vector of dimension last containing
c                       the error estimates
c
c              iord   - integer
c                       vector of dimension last, the first k elements
c                       of which contain pointers to the error
c                       estimates, such that
c                       elist(iord(1)),...,  elist(iord(k))
c                       form a decreasing sequence, with
c                       k = last if last.le.(limit/2+2), and
c                       k = limit+1-last otherwise
c
c              nrmax  - integer
c                       maxerr = iord(nrmax)
c
c***end prologue  dqpsrt
c
      double precision elist,ermax,errmax,errmin,phi,lambda1,zk0,
     *     Pup,Tup,rurd,xflow,kup
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,
     *  nrmax
      dimension elist(last),iord(last)
c
c           check whether the list contains more than
c           two error estimates.
c
c***first executable statement  dqpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
c
c           this part of the routine is only executed if, due to a
c           difficult integrand, subdivision increased the error
c           estimate. in the normal case the insert procedure should
c           start after the nrmax-th largest error estimate.
c
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
c ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
c
c           compute the number of elements in the list to be maintained
c           in descending order. this number depends on the number of
c           subdivisions still allowed.
c
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
c
c           insert errmax by traversing the list top-down,
c           starting comparison from the element elist(iord(nrmax+1)).
c
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
c ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
c
c           insert errmin by traversing the list bottom-up.
c
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
c ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
c
c           set maxerr and ermax.
c
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
