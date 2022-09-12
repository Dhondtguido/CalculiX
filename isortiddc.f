!
!     SLATEC: public domain
!
*deck isort
      subroutine isortiddc(ix,dy1,dy2,cy,n,kflag)
!
!     modified to make the same interchanges in two 
!     double (dy1 and dy2) and a char*20 aray (cy)
!
C***BEGIN PROLOGUE  ISORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2A
C***TYPE      INTEGER (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Kahaner, D. K., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   ISORT sorts array IX and optionally makes the same interchanges in
C   array IY.  The array IX may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      IX - integer array of values to be sorted
C      IY - integer array to be (optionally) carried along
C      N  - number of values in integer array IX to be sorted
C      KFLAG - control parameter
C            =  2  means sort IX in increasing order and carry IY along.
C            =  1  means sort IX in increasing order (ignoring IY)
C            = -1  means sort IX in decreasing order (ignoring IY)
C            = -2  means sort IX in decreasing order and carry IY along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761118  DATE WRITTEN
C   810801  Modified by David K. Kahaner.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to IX,IY. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!   100411  changed the dimension of IL and IU from 21 to 31.
!
!     field IL and IU have the dimension 31. This is log2 of the largest
!     array size to be sorted. If arrays larger than 2**31 in length have
!     to be sorted, this dimension has to be modified accordingly
!
C***END PROLOGUE  ISORT
!
      implicit none
C     .. Scalar Arguments ..
      integer kflag, n,iside,istat
C     .. Array Arguments ..
      integer ix(2,*)
      real*8 dy1(2,*),dy2(2,*)
      character*20 cy(*)
C     .. Local Scalars ..
      real r
      integer i, ij, j, k, kk, l, m, nn, t, tt
      real*8 tty11,tty12,ty11,ty12,tty21,tty22,ty21,ty22,ttx2,tx2
      character*20 uuy,uy
C     .. Local Arrays ..
      integer il(31), iu(31)
C     .. External Subroutines ..
!      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      intrinsic abs, int
C***FIRST EXECUTABLE STATEMENT  ISORT
!
      do i=1,n
         read(cy(i)(2:2),'(i1)',iostat=istat) iside 
         if(istat.gt.0) iside=0
         ix(1,i)=10*ix(1,i)+iside
      enddo
!
      nn = n
      if (nn .lt. 1) then
!         CALL XERMSG ('SLATEC', 'ISORT',
!     +      'The number of values to be sorted is not positive.', 1, 1)
         return
      endif
C
      kk = abs(kflag)
      if (kk.ne.1 .and. kk.ne.2) then
!         CALL XERMSG ('SLATEC', 'ISORT',
!     +      'The sort control parameter, K, is not 2, 1, -1, or -2.', 2,
!     +      1)
         return
      endif
C
C     Alter array IX to get decreasing order if needed
C
      if (kflag .le. -1) then
         do 10 i=1,nn
            ix(1,i) = -ix(1,i)
   10    continue
      endif
C
      if (kk .eq. 2) go to 100
C
C     Sort IX only
C
      m = 1
      i = 1
      j = nn
      r = 0.375e0
C
   20 if (i .eq. j) go to 60
      if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
      else
         r = r-0.21875e0
      endif
C
   30 k = i
C
C     Select a central element of the array and save it in location T
C
      ij = i + int((j-i)*r)
      t = ix(1,ij)
C
C     If first element of array is greater than T, interchange with T
C
      if (ix(1,i) .gt. t) then
         ix(1,ij) = ix(1,i)
         ix(1,i) = t
         t = ix(1,ij)
      endif
      l = j
C
C     If last element of array is less than than T, interchange with T
C
      if (ix(1,j) .lt. t) then
         ix(1,ij) = ix(1,j)
         ix(1,j) = t
         t = ix(1,ij)
C
C        If first element of array is greater than T, interchange with T
C
         if (ix(1,i) .gt. t) then
            ix(1,ij) = ix(1,i)
            ix(1,i) = t
            t = ix(1,ij)
         endif
      endif
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 l = l-1
      if (ix(1,l) .gt. t) go to 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 k = k+1
      if (ix(1,k) .lt. t) go to 50
C
C     Interchange these elements
C
      if (k .le. l) then
         tt = ix(1,l)
         ix(1,l) = ix(1,k)
         ix(1,k) = tt
         go to 40
      endif
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      if (l-i .gt. j-k) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 70
C
C     Begin again on another portion of the unsorted array
C
   60 m = m-1
      if (m .eq. 0) go to 190
      i = il(m)
      j = iu(m)
C
   70 if (j-i .ge. 1) go to 30
      if (i .eq. 1) go to 20
      i = i-1
C
   80 i = i+1
      if (i .eq. j) go to 60
      t = ix(1,i+1)
      if (ix(1,i) .le. t) go to 80
      k = i
C
   90 ix(1,k+1) = ix(1,k)
      k = k-1
      if (t .lt. ix(1,k)) go to 90
      ix(1,k+1) = t
      go to 80
C
C     Sort IX and carry IY along
C
  100 m = 1
      i = 1
      j = nn
      r = 0.375e0
C
  110 if (i .eq. j) go to 150
      if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
      else
         r = r-0.21875e0
      endif
C
  120 k = i
C
C     Select a central element of the array and save it in location T
C
      ij = i + int((j-i)*r)
      t = ix(1,ij)
      ty11 = dy1(1,ij)
      ty21 = dy1(2,ij)
      ty12 = dy2(1,ij)
      ty22 = dy2(2,ij)
      tx2 = ix(2,ij)
      uy = cy(ij)
C
C     If first element of array is greater than T, interchange with T
C
      if (ix(1,i) .gt. t) then
         ix(1,ij) = ix(1,i)
         ix(1,i) = t
         t = ix(1,ij)
         dy1(1,ij) = dy1(1,i)
         dy1(2,ij) = dy1(2,i)
         dy2(1,ij) = dy2(1,i)
         dy2(2,ij) = dy2(2,i)
         ix(2,ij) = ix(2,i)
         cy(ij) = cy(i)
         dy1(1,i) = ty11
         dy1(2,i) = ty21
         dy2(1,i) = ty12
         dy2(2,i) = ty22
         ix(2,i) = tx2
         cy(i) = uy
         ty11 = dy1(1,ij)
         ty21 = dy1(2,ij)
         ty12 = dy2(1,ij)
         ty22 = dy2(2,ij)
         tx2 = ix(2,ij)
         uy = cy(ij)
      endif
      l = j
C
C     If last element of array is less than T, interchange with T
C
      if (ix(1,j) .lt. t) then
         ix(1,ij) = ix(1,j)
         ix(1,j) = t
         t = ix(1,ij)
         dy1(1,ij) = dy1(1,j)
         dy1(2,ij) = dy1(2,j)
         dy2(1,ij) = dy2(1,j)
         dy2(2,ij) = dy2(2,j)
         ix(2,ij) = ix(2,j)
         cy(ij) = cy(j)
         dy1(1,j) = ty11
         dy1(2,j) = ty21
         dy2(1,j) = ty12
         dy2(2,j) = ty22
         ix(2,j) = tx2
         cy(j) = uy
         ty11 = dy1(1,ij)
         ty21 = dy1(2,ij)
         ty12 = dy2(1,ij)
         ty22 = dy2(2,ij)
         tx2 = ix(2,ij)
         uy = cy(ij)
C
C        If first element of array is greater than T, interchange with T
C
         if (ix(1,i) .gt. t) then
            ix(1,ij) = ix(1,i)
            ix(1,i) = t
            t = ix(1,ij)
            dy1(1,ij) = dy1(1,i)
            dy1(2,ij) = dy1(2,i)
            dy2(1,ij) = dy2(1,i)
            dy2(2,ij) = dy2(2,i)
            ix(2,ij) = ix(2,i)
            cy(ij) = cy(i)
            dy1(1,i) = ty11
            dy1(2,i) = ty21
            dy2(1,i) = ty12
            dy2(2,i) = ty22
            ix(2,i) = tx2
            cy(i) = uy
            ty11 = dy1(1,ij)
            ty21 = dy1(2,ij)
            ty12 = dy2(1,ij)
            ty22 = dy2(2,ij)
            tx2 = ix(2,ij)
            uy = cy(ij)
         endif
      endif
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 l = l-1
      if (ix(1,l) .gt. t) go to 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 k = k+1
      if (ix(1,k) .lt. t) go to 140
C
C     Interchange these elements
C
      if (k .le. l) then
         tt = ix(1,l)
         ix(1,l) = ix(1,k)
         ix(1,k) = tt
         tty11 = dy1(1,l)
         tty21 = dy1(2,l)
         tty12 = dy2(1,l)
         tty22 = dy2(2,l)
         ttx2 = ix(2,l)
         uuy = cy(l)
         dy1(1,l) = dy1(1,k)
         dy1(2,l) = dy1(2,k)
         dy2(1,l) = dy2(1,k)
         dy2(2,l) = dy2(2,k)
         ix(2,l) = ix(2,k)
         cy(l) = cy(k)
         dy1(1,k) = tty11
         dy1(2,k) = tty21
         dy2(1,k) = tty12
         dy2(2,k) = tty22
         ix(2,k) = ttx2
         cy(k) = uuy
         go to 130
      endif
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      if (l-i .gt. j-k) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 160
C
C     Begin again on another portion of the unsorted array
C
  150 m = m-1
      if (m .eq. 0) go to 190
      i = il(m)
      j = iu(m)
C
  160 if (j-i .ge. 1) go to 120
      if (i .eq. 1) go to 110
      i = i-1
C
  170 i = i+1
      if (i .eq. j) go to 150
      t = ix(1,i+1)
      ty11 = dy1(1,i+1)
      ty21 = dy1(2,i+1)
      ty12 = dy2(1,i+1)
      ty22 = dy2(2,i+1)
      tx2 = ix(2,i+1)
      uy = cy(i+1)
      if (ix(1,i) .le. t) go to 170
      k = i
C
  180 ix(1,k+1) = ix(1,k)
      dy1(1,k+1) = dy1(1,k)
      dy1(2,k+1) = dy1(2,k)
      dy2(1,k+1) = dy2(1,k)
      dy2(2,k+1) = dy2(2,k)
      ix(2,k+1) = ix(2,k)
      cy(k+1) = cy(k)
      k = k-1
      if (t .lt. ix(1,k)) go to 180
      ix(1,k+1) = t
      dy1(1,k+1) = ty11
      dy1(2,k+1) = ty21
      dy2(1,k+1) = ty12
      dy2(2,k+1) = ty22
      ix(2,k+1) = tx2
      cy(k+1) = uy
      go to 170
C
C     Clean up
C
  190 if (kflag .le. -1) then
         do 200 i=1,nn
            ix(1,i) = -ix(1,i)
  200    continue
      endif
!
      do i=1,nn
         read(cy(i)(2:2),'(i1)',iostat=istat) iside 
         if(istat.gt.0) iside=0
         ix(1,i)=(ix(1,i)-iside)/10
      enddo
!
      return
      end
