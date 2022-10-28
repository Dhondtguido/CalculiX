!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
C
C  1.TASK          INTERPOLATION OF A TWO DIMENSIONAL FUNCTION DEFINED POINT BY POINT
C  *********       THE X COORDINATES ARE USER SPECIFIED.
c                  THE INTERPOLATION TYPE CAN BE INDEPENDANTLY CHOSEN IN THE TWO DIRECTIONS
C                  EITHER CONSTANT, LINEAR OR DOUBLE QUADRATIC.
C                  BEYOND THE FIELD OF INTERPOLATION AN EXTRAOLATION IS CARRIED OUT.
C                  FOR ALL FOUR EXTRAPOLATION DIRECTIONS DIFFERENT EXTRAPOLATION METHOD
C                 (C ONSTANT,LINEAR,QUADRATIC) CAN BE CHOSEN, WHICH ORDER MUST NOT BE HIGHER
C                  THAN THE IONTERPOLATION ORDER
C
C  2.UP-AUFRUF     CALL TWODINT(T,LSP,IART,XA,YA,ZA,NA,IEXP,IER)
C  ***********                 T = MATRIX OF THE SAMPLE POINTS FORMATED AS FOLLOW
C                                  T(1,1) = NX + NY * 0.001
C                                         NX = NUMBER OF LINES T
C                                         NY = NUMBER OF COLUMNS T
C                                  T(1,2) ... T(1,NY)
C                                        VECTOR OF THE Y COORDINATES OF THE T MATRIX
C                                  T(2,1) ... T(NX,1)
C                                         VECTOR OF THE X COORDINATES OF THE T MATRIX
C                                  REST  OF T-MATRIX:
C                                          POINT(X,Y) OF THE T MATRIX 
C
C                              LSP  = COLUMN STEPOF  T
C                              IART = TYPE OF INTERPOLATION
C                                     IART = INTX * 10 + INTY
C                                     INTX INTERPOLATION TYPE IN X-DIRECTION
C                                     INTY INTERPOLATION TYPE IN Y-DIRECTION
C                              XA = VECTOR OF THE X COORDINATES OF THE VALUE TO BE INTERPOLATED
C                              YA = VECTOR OF THE Y COORDINATES OF THE VALUE TO BE INTERPOLATED
C                              ZA = VECTOR OF THE INTERPOLATED VALUES
C                              NA = ACTUAL LENGTH OF THE 3 PREVIOUS VECTORS
C                              IEXP = TWO ELEMENT VECTOR CONTRAINING THE TYPE OF EXTRAPOLATION 
C                                      CHOSEN BEYOND THE INTERPOLATION DOMAIN
C                                     IEXP(1): EXTRAPOLATION IN X-DIRECTION
C                                     IEXP(1) = IEXPX1 * 10 + IEXPXN
C                                     IEXPX1: EXTRAPOLATION BENEATH THE FIRST POINT
C                                     IEXPXN: EXTRAPOLATION BEYOND THE LAST POINT
C                                     IEXP(2): EXTRAPOLATION IN Y-DIRECTION
C                                     IEXP(2) = IEXPY1 * 10 + IEXPYN
C                                     SAME METHOD AS FOR IEXP(1):
C                              IER = ERROR CODE
C                                     IER = 0: NORMAL PROCEEDING
C                                     IER = -1: ERROR  INPUTDATA
C
C                           REMARK: CHOICE OF THE INTER- EXTRAPOLATION TYPE  IART AND IEXP -
C                         --------             ASSIGNEMENT OF  INTX,INTY,IEXPX1,
C                                       IEXPXN,IEXPY1,IEXPYN:
C                                        = 0 :   CONSTANT
C                                        = 1 :   LINEAR
C                                        = 2 :   DOUBLE QUADRATIC FROM
C                                                THE  SECOND UNTIL PENULTIMATE
C                                                INTERVAL IN THE INTERPOLATION MATRIX T,OTHERWISE QUADRATIC
C
C 3.RESTRICTIONS   THE SAMPLING POINT VECTORS (X UND Y COORDINATES
C ***************  OF THE MATRICX T MUST BE  STRICTLY MONOTONIC INCREASING SORTED
C                  THE PARAMETER FOR THE TYPE OF EXTRAPOLATION
c                  MUST NOT BE GREATER THAN THE ONE FOR TH EINTERPOLATION TYPE
C                  OTHERWISE THE  VALUE IS AUTOMATICALLY ADAPTATED
C                  IF THE NUMBER OF THE SAMPLING POINTS FOR THE REQUIRED TYPE OF INTERPOLATION IS TOO SMALL,
C                  THE DEGREE OF INTERPOLATION WILL BE ACCORDINGLY ADAPTATED
C
C  4.USED UP'S     ONEDINT (ONE DIMENSIONAL INTERPOLATION ANALOG TO THIS PROGRAMM)
C

      SUBROUTINE TWODINT (T,LSP,IART,XA,YA,ZA,NA,IEXP,IER)
      implicit none
      INTEGER IEXP(2),IYU,IYO,IXU,IXO,IDX,IDY,LL,INPY,IEXPX1,IEXPXN,
     &  IEXPY1,IEXPYN,LX,LY,INPX,IART,LSP,IER,NX,NY,L,NA,one
      REAL*8 T(LSP,1),XA(1),YA(1),ZA(1)
      REAL*8 Z1(4),Z2(4)
C      ENTRY ZWEINT (T,LSP,IART,XA,YA,ZA,NA,IEXP,IER)
      IER = 0
      one=1
      NX = T(1,1)
      NY = (T(1,1)-NX)*1000 + 0.1d0
C
C TESTING INPUT
C--------------
      IF ((NX-2).lt.0) then
         go to 900
      elseif((nx-2).eq.0) then
         go to 30
      else
         go to 10
      endif
   10 DO 20 L = 3,NX
   20 IF ((T(L,1)-T(L-1,1)) .LE. 0) GO TO 900
   30 IF ((NY-2).lt.0) then
         go to 900
      elseif((ny-2).eq.0) then
         go to 60
      else
         go to 40
      endif
   40 DO 50 L = 3,NY
   50 IF ((T(1,L)-T(1,L-1)) .LE. 0) GO TO 900
   60 IF (NA .LE. 0) GO TO 900
C
C DEFINING THE CONTROL VALUES
C---------------------------
      INPX = IART/10
      INPY = IART - INPX*10 + 0.1d0
      IEXPX1 = IEXP(1)/10
      IEXPXN = IEXP(1) - IEXPX1*10
      IEXPY1 = IEXP(2)/10
      IEXPYN = IEXP(2) - IEXPY1*10
      IF (NX-2 .LT. INPX) INPX = NX - 2
      IF (NY-2 .LT. INPY) INPY = NY - 2
      IF (IEXPX1 .GT. INPX) IEXPX1 = INPX
      IF (IEXPXN .GT. INPX) IEXPXN = INPX
      IF (IEXPY1 .GT. INPY) IEXPY1 = INPY
      IF (IEXPYN .GT. INPY) IEXPYN = INPY
C
C SUCCESSIVE PROCESSING THE INTERPOLATION EXIGENCES
C-------------------------------------------------------
      DO 400 L = 1,NA
      LX = 2
C
C SETTING REFERENCE POINTS (LX,LY) 
C---------------------------------
  200 IF (XA(L) .LT. T(LX,1)) GO TO 220
      LX = LX + 1
      IF ((LX-NX).le.0) then
         go to 200
      else
         go to 210
      endif
  210 LX = NX
  220 DO 230 LY = 2,NY
  230 IF (YA(L) .LT. T(1,LY)) GO TO 235
      LY = NY
  235 IYU = LY - INPY
      IYO = LY + INPY - 1
      IF (IYU .GE. 2) GO TO 240
      IYU = 2
      IYO = IYU + INPY
  240 IF (IYO .GT. NY) IYO = NY
      IXU = LX - INPX
      IXO = LX + INPX - 1
      IF (IXU .GE. 2) GO TO 245
      IXU = 2
      IXO = IXU + INPX
  245 IF (IXO .GT. NX) IXO = NX
      IDX = IXO - IXU + 1
      IF (IXU .LT. IXO) GO TO 270
      IF (IYU .LT. IYO) GO TO 250
C
C CONSTANT INTERPOLATION
C------------------------
      IF (LX .GT. 2 .AND. XA(L) .LT. T(NX,1)) LX = LX - 1
      IF (LY .GT. 2 .AND. YA(L) .LT. T(1,NY)) LY = LY - 1
      ZA(L) = T(LX,LY)
      GO TO 400
C
C LINEAR AND  QUADRATIC INTERPOLATION USING ONEDINT (ONEDIMENSIONAL)
C---------------------------------------------------------------------
C
C INTERPOLATION ONLY IN Y-DIRECTION
C
  250 IDY = 0
      DO 260 LL = IYU,IYO
      IDY = IDY + 1
      Z1(IDY) = T(1,LL)
  260 Z2(IDY) = T(LX,LL)
      GO TO 300
C
C INTERPOLATION ONLY IN X-DIRECTION
C
  270 IF (IYU .LT. IYO) GO TO 280
      CALL ONEDINT(T(IXU,1),T(IXU,LY),IDX,XA(L),ZA(L),one,INPX,IEXP(1),
     1           IER)
      IF (IER.eq.0) then
         go to 400
      else
         go to 900
      endif
C
C 1.INTERPOLATION STEP IN X-DIRECTION
C
  280 IDY = 0
      DO 290 LL = IYU,IYO
      IDY = IDY + 1
      Z1(IDY) = T(1,LL)
      CALL ONEDINT (T(IXU,1),T(IXU,LL),IDX,XA(L),Z2(IDY),one,INPX,
     1    IEXP(1),IER)
      IF (IER.eq.0) then
         go to 290
      else
         go to 900
      endif
  290 CONTINUE
C
C 1.OR 2.INTERPOLATION STEP IN Y-DIRECTION
C
  300 CALL ONEDINT (Z1,Z2,IDY,YA(L),ZA(L),one,INPY,IEXP(2),IER)
      IF (IER.eq.0) then
         go to 400
      else
         go to 900
      endif
C
C RETURN BY NORMAL PROCEEDING
C--------------------------------
  400 CONTINUE
      IER = 0
      RETURN
C
C ERROR RETURN
C-------------
  900 IER = -1
      RETURN
      END
