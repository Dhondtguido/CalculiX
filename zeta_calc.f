!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2005 Guido Dhondt
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
!     This subroutine enable to compuite the different zeta exponents for
!     the different partial total head loss restrictors. The values of the
!     'zetas' have been found in the following published works
!
!     I.E. IDEL'CHIK 'HANDBOOK OF HYDRAULIC RESISTANCE'
!     2nd edition 1986,HEMISPHERE PUBLISHING CORP.
!     ISBN 0-899116-284-4
! 
!     D.S. MILLER 'INTERNAL FLOW SYSTEMS'
!     1978,vol.5 B.H.R.A FLUID ENGINEERING 
!     ISBN 0-900983-78-7
!
!     author: Yannick Muller
!
      subroutine zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &     isothermal,kon,ipkon,R,kappa,v,mi,iaxial)
!
      implicit none
!
      logical isothermal
!
      character*8 lakon(*)
! 
      integer ielprop(*),nelem,iexp(2),i,j,ier,iwrit1,iexp3(2),
     &     iwrit2,nelem_ref,ipkon(*),kon(*),nelem0,nelem1,nelem2,node10,
     &     node20,nodem0,node11,node21,nodem1,node12,node22,nodem2,
     &     iexpbr1(2) /11,11/,icase,node0,node1,node2,mi(*),n0,n1,n2,
     &     n5,n6,n8,n10,n11,n12,n14,n15,n22,iaxial,index
!
      real*8 zeta,prop(*),lzd,reynolds,ereo,fa2za1,zetap,zeta0,
     &     lambda,thau,a1,a2,dh,dl,a2za1,ldumm,dhdumm,ks,
     &     form_fact,zeta01,zeta02,alpha,rad,delta,a0,b0,azb,rzdh,
     &     A,C,rei,lam,ai,b1,c1,b2,c2,zeta1,re_val,k,ldre,
     &     zetah,cd,cdu,km,Tt0,Ts0,Tt1,Ts1,Tt2,Ts2,
     &     rho0,rho1,rho2,V0,V1,v2,a0a1,a0a2,zetlin,lam10,lam20,pi,
     &     alpha1,alpha2,R,kappa,ang1s,ang2s,cang1s,cang2s,
     &     v(0:mi(2),*),V1V0,V2V0,z1_60,z1_90,
     &     z2_60,z2_90,afakt,V2V0L,kb,ks2,a2a0,Z90LIM11,Z90LIM51,
     &     lam11,lam12,lam21,lam22,W2W0,W1W0,dh0,dh2,hq,z2d390,
     &     z1p090,z90,z60,pt0,pt2,pt1,M0,M1,M2,W0W1,W0W2,
     &     xflow0,xflow1,xflow2,Qred_0, Qred_1, Qred_2,Qred_crit
!
!     THICK EDGED ORIFICE IN STRAIGHT CONDUIT (L/DH > 0.015)
!     I.E. IDEL' CHIK (SECTION III PAGE 140)
!
!     I.E. IDEL'CHIK 'HANDBOOK OF HYDRAULIC RESISTANCE'
!     2nd edition 1986,HEMISPHERE PUBLISHING CORP.
!     ISBN 0-899116-284-4
!
!        ***** long orifice *****
!
!        DIAGRAMS 4-19 p 175 - Reynolds R:epsilon^-_oRe
!
      real*8 XRE (14), YERE (14)
      data XRE /25.d0,40.d0,60.0d0,100.d0,200.d0,400.d0,1000.d0,2000.d0,
     &     4000.d0,10000.d0,20000.d0,100000.d0,200000.d0,1000000./
      data YERE /0.34d0,0.36d0,0.37d0,0.40d0,0.42d0,0.46d0,0.53d0,
     &     0.59d0,0.64d0,0.74d0,0.81d0,0.94d0,0.95d0,0.98/
!     
!     Diagram 4-19 p 175 - Reynolds | A1/A2 R: zeta_phi
!     
      real*8 zzeta (15,11)
      data ((zzeta(i,j),i=1,15),j=1,11) 
     &     /15.011  d0,25.0d0,40.0d0,60.0d0,100.0d0,200.0d0,400.0d0,
     &      1000.0d0,2000.0d0,4000.0d0,10000.0d0,20000.0d0,100000.0d0,
     &      200000.0d0,1000000.0d0,0.d0,1.94d0,1.38d0,1.14d0,0.89d0,
     &      0.69d0,0.64d0,0.39d0,0.30d0,0.22d0,0.15d0,0.11d0,0.04d0,
     &      0.01d0,0.00d0,0.20d0,1.78d0,1.36d0,1.05d0,0.85d0,0.67d0,
     &      0.57d0,0.36d0,0.26d0,0.20d0,0.13d0,0.09d0,0.03d0,0.01d0,
     &      0.00d0,0.30d0,1.57d0,1.16d0,0.88d0,0.75d0,0.57d0,0.43d0,
     &      0.30d0,0.22d0,0.17d0,0.10d0,0.07d0,0.02d0,0.01d0,0.00d0,
     &      0.40d0,1.35d0,0.99d0,0.79d0,0.57d0,0.40d0,0.28d0,0.19d0,
     &      0.14d0,0.10d0,0.06d0,0.04d0,0.02d0,0.01d0,0.00d0,0.50d0,
     &      1.10d0,0.75d0,0.55d0,0.34d0,0.19d0,0.12d0,0.07d0,0.05d0,
     &      0.03d0,0.02d0,0.01d0,0.01d0,0.01d0,0.00d0,0.60d0,0.85d0,
     &      0.56d0,0.30d0,0.19d0,0.10d0,0.06d0,0.03d0,0.02d0,0.01d0,
     &      0.01d0,0.00d0,0.00d0,0.00d0,0.00d0,0.70d0,0.58d0,0.37d0,
     &      0.23d0,0.11d0,0.06d0,0.03d0,0.02d0,0.01d0,0.00d0,0.00d0,
     &      0.00d0,0.00d0,0.00d0,0.00d0,0.80d0,0.40d0,0.24d0,0.13d0,
     &      0.06d0,0.03d0,0.02d0,0.01d0,0.00d0,0.00d0,0.00d0,0.00d0,
     &      0.00d0,0.00d0,0.00d0,0.90d0,0.20d0,0.13d0,0.08d0,0.03d0,
     &      0.01d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
     &      0.00d0,0.00d0,0.95d0,0.03d0,0.03d0,0.02d0,0.00d0,0.00d0,
     &      0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
     &      0.00d0/
!
!     Diagram 4-12 p 169 - l/Dh R: tau
!     
      real*8 XLZD (10), YTOR (10)
      data XLZD /0.0d0,0.2d0,0.4d0,0.6d0,0.8d0,1.0d0,1.2d0,1.6d0,2.0d0,
     &           2.4d0/
      data YTOR /1.35d0,1.22d0,1.10d0,0.84d0,0.42d0,0.24d0,0.16d0,
     &           0.07d0,0.02d0,0.0d0/
      data IEXP /10, 1/
!     
!     ***** wall orifice *****
!     
!     THICK-WALLED ORIFICE IN LARGE WALL (L/DH > 0.015)
!     I.E. IDL'CHIK (page 174)
!     
!     DIAGRAM 4-18 A - l/Dh R: zeta_o
!     
      real*8 XLQD(12)
      DATA XLQD /
     &     0.d0,0.2d0,0.4d0,0.6d0,0.8d0,1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,
     &     2.0d0,10.0d0/
      real*8 YZETA1(12)
      DATA YZETA1 /
     &     2.85d0,2.72d0,2.6d0,2.34d0,1.95d0,1.76d0,1.67d0,1.62d0,1.6d0,
     &     1.58d0,1.55d0,1.55d0/
!     
!     DIAGRAM 4-19 p175 first line - Re (A1/A2=0) R: zeta_phi
!     
      real*8 XRE2(14)
      DATA XRE2 /
     &     25.d0,40.d0,60.d0,100.d0,200.d0,400.d0,1000.d0,2000.d0,
     &     4000.d0,10000.d0,20000.d0,50000.d0,100000.d0,1000000.d0/
      real*8 YZETA2(14)
      DATA YZETA2 /
     &     1.94d0,1.38d0,1.14d0,.89d0,.69d0,.54d0,.39d0,.3d0,.22d0,
     &     .15d0,.11d0,.04d0,.01d0,0.d0/
!     
!     Diagram 4-18 p174 first case * (=multiplication) epsilon^-_oRe p 175
!     
      real*8 YERE2(14)
      DATA YERE2 /
     &     1.d0,1.05d0,1.09d0,1.15d0,1.23d0,1.37d0,1.56d0,1.71d0,1.88d0,
     &     2.17d0,2.38d0,2.56d0,2.72d0,2.85d0/
!
!     ***** expansion *****
!     
!     SUDDEN EXPANSION OF A STREAM WITH UNIFORM VELOCITY DISTRIBUTION
!     I.E. IDL'CHIK (page 160)
!
!     DIAGRAM 4-1 - Re | A1/A2 R:zeta
!
      real*8 ZZETA3(14,8)
      DATA ZZETA3 /
     &     14.008d0, 10.000d0,15.0d0,20.0d0,30.0d0,40.0d0,50.0d0,
     &     100.0d0,200.0d0,500.0d0,1000.0d0,2000.0d0,3000.0d0,3500.0d0,
     &  .01d0,3.10d0,3.20d0,3.00d0,2.40d0,2.15d0,1.95d0,1.70d0,
     &     1.65d0,1.70d0,2.00d0,1.60d0,1.00d0,1.00d0,
     &  0.1d0,3.10d0,3.20d0,3.00d0,2.40d0,2.15d0,1.95d0,1.70d0,1.65d0,
     &     1.70d0,2.00d0,1.60d0,1.00d0,0.81d0,
     &  0.2d0,3.10d0,3.20d0,2.80d0,2.20d0,1.85d0,1.65d0,1.40d0,1.30d0,
     &     1.30d0,1.60d0,1.25d0,0.70d0,0.64d0,
     &  0.3d0,3.10d0,3.10d0,2.60d0,2.00d0,1.60d0,1.40d0,1.20d0,1.10d0,
     &     1.10d0,1.30d0,0.95d0,0.60d0,0.50d0,
     &  0.4d0,3.10d0,3.00d0,2.40d0,1.80d0,1.50d0,1.30d0,1.10d0,1.00d0,
     &     0.85d0,1.05d0,0.80d0,0.40d0,0.36d0,
     &  0.5d0,3.10d0,2.80d0,2.30d0,1.65d0,1.35d0,1.15d0,0.90d0,0.75d0,
     &     0.65d0,0.90d0,0.65d0,0.30d0,0.25d0,
     &  0.6d0,3.10d0,2.70d0,2.15d0,1.55d0,1.25d0,1.05d0,0.80d0,0.60d0,
     &     0.40d0,0.60d0,0.50d0,0.20d0,0.16d0/
!     
      DATA IEXP3 /0,0/
!     
!     ***** contraction *****
!
!     SUDDEN CONTRACTION WITH & WITHOUT CONICAL BELLMOUTH ENTRY
!     I.E. IDL'CHIK  p 168
! 
!     DIAGRAM 4-10 - Re | A1/A2 R: zeta
!
      real*8 ZZETA41(14,7)
      DATA ZZETA41 /
     & 14.007 d0,10.0d0,20.0d0,30.0d0,40.0d0,50.0d0,100.0d0,200.0d0,
     &        500.0d0,1000.0d0,2000.0d0,4000.0d0,5000.0d0,10000.0d0,
     &0.1d0,5.00d0,3.20d0,2.40d0,2.00d0,1.80d0,1.30d0,1.04d0,0.82d0,
     &        0.64d0,0.50d0,0.80d0,0.75d0,0.50d0,
     &0.2d0,5.00d0,3.10d0,2.30d0,1.84d0,1.62d0,1.20d0,0.95d0,0.70d0,
     &        0.50d0,0.40d0,0.60d0,0.60d0,0.40d0,
     &0.3d0,5.00d0,2.95d0,2.15d0,1.70d0,1.50d0,1.10d0,0.85d0,0.60d0,
     &        0.44d0,0.30d0,0.55d0,0.55d0,0.35d0,
     &0.4d0,5.00d0,2.80d0,2.00d0,1.60d0,1.40d0,1.00d0,0.78d0,0.50d0,
     &        0.35d0,0.25d0,0.45d0,0.50d0,0.30d0,
     &0.5d0,5.00d0,2.70d0,1.80d0,1.46d0,1.30d0,0.90d0,0.65d0,0.42d0,
     &        0.30d0,0.20d0,0.40d0,0.42d0,0.25d0,
     &0.6d0,5.00d0,2.60d0,1.70d0,1.35d0,1.20d0,0.80d0,0.56d0,0.35d0,
     &        0.24d0,0.15d0,0.35d0,0.35d0,0.20d0/
!
!      Diagram 3-7 p128  - alpha | l/Dh R: zeta
!
      real*8 ZZETA42(10,7)
      DATA ZZETA42 /
     & 10.007d0,0.d0,10.0d0,20.0d0,30.0d0,40.0d0,60.0d0,100.0d0,140.0d0,
     &    180.0d0,
     &  0.025d0,0.50d0,0.47d0,0.45d0,0.43d0,0.41d0,0.40d0,0.42d0,0.45d0,
     &    0.50d0,
     &  0.050d0,0.50d0,0.45d0,0.41d0,0.36d0,0.33d0,0.30d0,0.35d0,0.42d0,
     &    0.50d0,
     &  0.075d0,0.50d0,0.42d0,0.35d0,0.30d0,0.26d0,0.23d0,0.30d0,0.40d0,
     &    0.50d0,
     &  0.100d0,0.50d0,0.39d0,0.32d0,0.25d0,0.22d0,0.18d0,0.27d0,0.38d0,
     &    0.50d0,
     &  0.150d0,0.50d0,0.37d0,0.27d0,0.20d0,0.16d0,0.15d0,0.25d0,0.37d0,
     &    0.50d0,
     &  0.600d0,0.50d0,0.27d0,0.18d0,0.13d0,0.11d0,0.12d0,0.23d0,0.36d0,
     &    0.50d0/
!
!     ***** bends *****
!     
!     SHARP ELBOW (R/DH=0) AT 0 < DELTA < 180
!     I.E. IDL'CHIK page 294
!     DIAGRAM 6-5  - a0/b0 R: C1
!     
      real*8  XAQB(12)
      DATA XAQB /
     &     0.25d0,0.50d0,0.75d0,1.00d0,1.50d0,2.00d0,3.00d0,4.00d0,
     &     5.00d0,6.00d0,7.00d0,8.00d0/
!     
      real*8 YC(12)
      DATA YC /
     &     1.10d0,1.07d0,1.04d0,1.00d0,0.95d0,0.90d0,0.83d0,0.78d0,
     &     0.75d0,0.72d0,0.71d0,0.70d0/
!
!     DIAGRAM 6-5 - delta R: A
!
      real*8 XDELTA(10)
      DATA XDELTA /
     &     20.0d0,30.0d0,45.0d0,60.0d0,75.0d0,90.0d0,110.d0,130.d0,
     &     150.d0,180.d0/
!     
      real*8 YA(10)
      DATA YA /
     &     2.50d0,2.22d0,1.87d0,1.50d0,1.28d0,1.20d0,1.20d0,1.20d0,
     &     1.20d0,1.20d0/
!     
!     SHARP BENDS 0.5 < R/DH < 1.5 AND 0 < DELTA < 180
!     I.E. IDL'CHIK page 289-290
!     DIAGRAM 6-1  (- delta from diagram 6-5) R: A1
!     
      real*8 YA1(10)
      DATA YA1 /
     &     0.31d0,0.45d0,0.60d0,0.78d0,0.90d0,1.00d0,1.13d0,1.20d0,
     &     1.28d0,1.40d0/
!     
!     DIAGRAM 6-1 - R0/D0 R: B1
!     
      real*8 XRQDH(8)
      DATA XRQDH /
     &     0.50d0,0.60d0,0.70d0,0.80d0,0.90d0,1.00d0,1.25d0,1.50d0/
!     
      real*8 YB1(8)
      DATA YB1 /
     &     1.18d0,0.77d0,0.51d0,0.37d0,0.28d0,0.21d0,0.19d0,0.17d0/
!     
!     DIAGRAM 6-1 (- a0/b0 from diagram 6-5) R: C1
!     
      real*8 YC1(12)
      DATA YC1 /
     &     1.30d0,1.17d0,1.09d0,1.00d0,0.90d0,0.85d0,0.85d0,0.90d0,
     &     0.95d0,0.98d0,1.00d0,1.00d0/
!     
!     SMOOTH BENDS (R/DH > 1.5) AT 0 < DELTA < 180
!     I.E. IDL'CHIK 
!     
!     DIAGRAM 6-1  - R0/D0 R: B1 (continuation of XRQDH)
!
      real*8 XRZDH(14)
      DATA XRZDH /
     &     1.00d0,2.00d0,4.00d0,6.00d0,8.00d0,10.0d0,15.0d0,20.0d0,
     &     25.0d0,30.0d0,35.0d0,40.0d0,45.0d0,50.0d0/
!     
      real*8 YB2(14)
      DATA YB2 /
     &     0.21d0,0.15d0,0.11d0,0.09d0,0.07d0,0.07d0,0.06d0,0.05d0,
     &     0.05d0,0.04d0,0.04d0,0.03d0,0.03d0,0.03d0/
!
!     (- a0/b0 from Diagram 6-5) R: C2
!
      real*8 YC2(12)
      DATA YC2 /
     &     1.80d0,1.45d0,1.20d0,1.00d0,0.68d0,0.45d0,0.40d0,0.43d0,
     &     0.48d0,0.55d0,0.58d0,0.60d0/
!     
!     D.S. MILLER 'INTERNAL FLOW SYSTEMS'
!     1978,vol.5 B.H.R.A FLUID ENGINEERING SERIES
!     ISBN 0-900983-78-7     
!
!        SMOOTH BENDS B.H.R.A HANDBOOK P.141 
!
      REAL*8 ZZETAO(14,15)
      DATA((ZZETAO(I,J),I=1,14),J=1,8) /
     & 14.015d0,0.5d0,0.6d0,0.8d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,
     &     7.0d0,8.0d0,9.0d0,10.d0,
     & 10.00d0, 0.030d0,0.025d0,0.021d0,0.016d0,0.022d0,0.030d0,0.034d0,
     &     0.036d0,0.040d0,0.042d0,0.043d0,0.044d0,0.044d0,
     & 15.00d0, 0.036d0,0.035d0,0.025d0,0.025d0,0.033d0,0.042d0,0.045d0,
     &     0.050d0,0.055d0,0.055d0,0.058d0,0.060d0,0.063d0,
     & 20.00d0, 0.056d0,0.046d0,0.034d0,0.034d0,0.045d0,0.054d0,0.056d0,
     &     0.062d0,0.066d0,0.067d0,0.072d0,0.075d0,0.080d0,
     & 30.00d0, 0.122d0,0.094d0,0.063d0,0.056d0,0.063d0,0.071d0,0.075d0,
     &     0.082d0,0.087d0,0.089d0,0.097d0,0.101d0,0.110d0,
     & 40.00d0, 0.220d0,0.160d0,0.100d0,0.085d0,0.080d0,0.086d0,0.092d0,
     &     0.100d0,0.106d0,0.122d0,0.121d0,0.126d0,0.136d0,
     & 50.00d0, 0.340d0,0.245d0,0.148d0,0.117d0,0.097d0,0.100d0,0.108d0,
     &     0.116d0,0.123d0,0.133d0,0.144d0,0.150d0,0.159d0,
     & 60.00d0, 0.480d0,0.350d0,0.196d0,0.150d0,0.115d0,0.116d0,0.122d0,
     &     0.131d0,0.140d0,0.153d0,0.164d0,0.171d0,0.181d0/
      DATA((ZZETAO(I,J),I=1,14),J=9,15) /
     & 70.00d0, 0.645d0,0.466d0,0.243d0,0.186d0,0.132d0,0.130d0,0.136d0,
     &     0.148d0,0.160d0,0.172d0,0.185d0,0.191d0,0.200d0,
     & 80.00d0, 0.827d0,0.600d0,0.288d0,0.220d0,0.147d0,0.142d0,0.150d0,
     &     0.166d0,0.180d0,0.191d0,0.203d0,0.209d0,0.218d0,
     & 90.00d0, 1.000d0,0.755d0,0.333d0,0.247d0,0.159d0,0.155d0,0.166d0,
     &     0.185d0,0.197d0,0.209d0,0.220d0,0.227d0,0.236d0,
     & 100.0d0, 1.125d0,0.863d0,0.375d0,0.264d0,0.167d0,0.166d0,0.183d0,
     &     0.202d0,0.214d0,0.225d0,0.238d0,0.245d0,0.255d0,
     & 120.0d0, 1.260d0,0.983d0,0.450d0,0.281d0,0.180d0,0.188d0,0.215d0,
     &     0.234d0,0.247d0,0.260d0,0.273d0,0.282d0,0.291d0,
     & 150.0d0, 1.335d0,1.060d0,0.536d0,0.289d0,0.189d0,0.214d0,0.251d0,
     &     0.272d0,0.297d0,0.312d0,0.325d0,0.336d0,0.346d0,
     & 180.0d0, 1.350d0,1.100d0,0.600d0,0.290d0,0.190d0,0.225d0,0.280d0,
     &     0.305d0,0.347d0,0.364d0,0.378d0,0.390d0,0.400d0/
!       
      REAL*8 KRE(22,4)
      DATA KRE  /
     & 22.004d0,1.d+3,2.d+3,3.d+3,4.d+3,5.d+3,6.d+3,7.d+3,8.d+3,9.d+3,
     &        1.d+4,2.d+4,3.d+4,4.d+4,6.d+4,8.d+4,1.d+5,2.d+5,3.d+5,
     &        5.d+5,7.d+5,1.d+6,
     & 1.0d0,3.88d0,3.06d0,2.77d0,2.60d0,2.49d0,2.40d0,2.33d0,2.27d0,
     &        2.22d0,2.18d0,1.86d0,1.69d0,1.57d0,1.41d0,1.30d0,1.22d0,
     &        5*1.00d0,
     & 1.5d0,3.88d0,3.06d0,2.77d0,2.60d0,2.49d0,2.40d0,2.33d0,2.27d0,
     &        2.22d0,2.18d0,1.90d0,1.76d0,1.67d0,1.54d0,1.46d0,1.40d0,
     &        1.22d0,1.12d0,3*1.00d0,
     & 2.0d0,3.88d0,3.06d0,2.77d0,2.60d0,2.49d0,2.40d0,2.33d0,2.27d0,
     &        2.22d0,2.18d0,1.93d0,1.80d0,1.71d0,1.60d0,1.53d0,1.47d0,
     &        1.32d0,1.23d0,1.13d0,1.06d0,1.00d0/
!
      integer iexp6(2)
      DATA iexp6 /0,0/
!
!     Campbell, Slattery
!     "Flow in the entrance of a tube"
!     Journal of Basic Engineering, 1963
!
!     EXIT LOSS COEFFICIENT FOR LAMINAR FLOWS DEPENDING ON THE
!     ACTUAL VELOCITY DISTRIBUTION AT THE EXIT
!
      real*8 XDRE(12)
      DATA XDRE /
     &        0.000d0,0.001d0,0.0035d0,0.0065d0,0.010d0,0.0150d0,
     &        0.020d0,0.025d0,0.035d0,0.045d0,0.056d0,0.065d0/
!
      real*8 ZETAEX(12)
      DATA ZETAEX /
     &        1.00d0,1.200d0,1.40d0,1.54d0,1.63d0,1.73d0,1.80d0,1.85d0,
     &        1.93d0,1.97d0,2.00d0,2.00d0/
!
!     Branch Joint Genium 
!     Branching Flow Part IV - TEES
!     Fluid Flow Division
!     Section 404.2 page 4 December 1986
!     Genium Publishing (see www.genium.com)
!
!     n.b: the values of this table have been scaled by a factor 64.
!
      real*8 XANG(11),YANG(11)
      data (XANG(i),YANG(i),i=1,11)
     &       /0.0d0,62.d0,
     &        15.d0,62.d0,
     &        30.d0,61.d0,
     &        45.d0,61.d0,
     &        60.d0,58.d0,
     &        75.d0,52.d0,
     &        90.d0,40.d0,
     &        105.d0,36.d0,
     &        120.d0,34.d0,
     &        135.d0,33.d0,
     &        150.d0,32.5d0/
!
!     Branch Joint Idelchik 1
!     Diagrams of resistance coefficients 
!     I.E. IDEL'CHIK 'HANDBOOK OF HYDRAULIC RESISTANCE'
!     2nd edition 1986,HEMISPHERE PUBLISHING CORP.
!     ISBN 0-899116-284-4
!
      real*8 TA2A0(12),TAFAKT(12)
      data (TA2A0(i),TAFAKT(i),i=1,12)
     &        /0.d0  ,1.d0  ,
     &        0.16d0 ,1.d0  ,
     &        0.20d0 ,0.99d0,
     &        0.25d0 ,0.95d0,
     &        0.29d0 ,0.90d0,
     &        0.31d0 ,0.85d0,
     &        0.33d0 ,0.80d0,
     &        0.35d0 ,0.78d0,
     &        0.4d0  ,0.75d0,
     &        0.6d0  ,0.70d0,
     &        0.8d0  ,0.65d0,
     &        1.d0   ,0.60d0/   
!
!     Branch Joint Idelchik 2
!     Diagrams of resistance coefficients p348-351 section VII
!     I.E. IDEL'CHIK 'HANDBOOK OF HYDRAULIC RESISTANCE'
!     2nd edition 1986,HEMISPHERE PUBLISHING CORP.
!     ISBN 0-899116-284-4
!
!     page 352 diagram 7-9 - alpha | Fs/Fc
!
      real*8 KBTAB(6,7),KSTAB(6,6)
      data ((KBTAB(i,j),j=1,7),i=1,6)
     &        /6.007d0 ,0.d0,15.d0,30.d0,45.d0,60.d0  ,90.d0  ,
     &           0.d0  ,0.d0, 0.d0, 0.d0, 0.d0, 0.d0  , 0.d0  ,
     &           0.1d0 ,0.d0, 0.d0, 0.d0, 0.d0, 0.d0  , 0.d0  ,
     &           0.2d0 ,0.d0, 0.d0, 0.d0, 0.d0, 0.d0  , 0.1d0 ,
     &           0.33d0,0.d0, 0.d0, 0.d0, 0.d0, 0.d0  , 0.2d0 ,
     &           0.5d0 ,0.d0, 0.d0, 0.d0, 0.d0, 0.1d0 , 0.25d0/
!
!     page 348-351 diagrams 7-5 to 7-8 - alpha | Fs/Fc
!
      data ((KSTAB(i,j),j=1,6),i=1,6)
     &        /6.006d0 ,0.d0,15.d0  ,30.d0  ,45.d0  , 60.d0  ,
     &           0.d0  ,0.d0, 0.d0  , 0.d0  , 0.d0  , 0.d0   , 
     &           0.1d0 ,0.d0, 0.d0  , 0.d0  , 0.05d0, 0.d0   ,
     &           0.2d0 ,0.d0, 0.d0  , 0.d0  , 0.14d0, 0.d0   ,
     &           0.33d0,0.d0, 0.14d0, 0.17d0, 0.14d0, 0.1d0  ,
     &           0.5d0 ,0.d0, 0.4d0 , 0.4d0 , 0.3d0 , 0.25d0/ 
!
!     page 352 diagram 7-9 R: zeta_c,st
!
      real*8 Z90TAB(6,13)
      data  ((Z90TAB(i,j),j=1,13),i=1,6)/
     & 6.013d0,0.d0,0.03d0,0.05d0,0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,
     & 0.7d0,0.8d0,1.0d0,
     & .06d0,.02d0,.05d0,.08d0,.08d0,.07d0,.01d0,-.15d0,1.d9,1.d9,
     & 1.d9,1.d9,1.d9,
     & .10d0,.04d0,.08d0,.10d0,.20d0,.26d0,.20d0,.05d0,-.13d0,
     & 1.d9,1.d9,1.d9,1.d9,
     & .20d0,.08d0,.12d0,.18d0,.25d0,.34d0,.32d0,.26d0,.16d0,
     & .02d0,-.14d0,1.d9,1.d9,
     & .33d0,.45d0,.50d0,.52d0,.59d0,.66d0,.64d0,.62d0,.58d0,
     & .44d0,.27d0,.08d0,-.34d0,
     & .50d0,1.00d0,1.04d0,1.06d0,1.16d0,1.25d0,1.25d0,1.22d0,1.10d0,
     & .88d0,.70d0,.45d0,0./
!
!     table to check the location of V2V0 in Z90TAB 
!
      real*8 Z90LIMX (5),Z90LIMY(5)
      data Z90LIMX    
     &        /0.06d0,0.1d0,0.2d0,0.33,0.5d0 /
!
      data Z90LIMY
     &   /0.1d0,0.1d0,0.3d0,0.5d0,0.7d0/
!
      data n0 /0/
      data n1 /1/
      data n2 /2/
      data n5 /5/
      data n6 /6/
      data n8 /8/
      data n10 /10/
      data n11 /11/
      data n12 /12/
      data n14 /14/
      data n15 /15/
      data n22 /22/
!     
      pi=4.d0*datan(1.d0)
      index=ielprop(nelem)
!     
      if((lakon(nelem)(2:5).eq.'REUS').or.
     &    (lakon(nelem)(2:5).eq.'LPUS')) then
!     
!     user defined zeta
!     
         zeta=prop(index+4)
!     
         return
!
      elseif((lakon(nelem)(2:5).eq.'REEN').or.
     &       (lakon(nelem)(2:5).eq.'LPEN')) then
!     
!     entrance 
!     
         zeta=prop(index+4)
!     
         return
!     
      elseif((lakon(nelem)(2:7).eq.'RELOID').or.
     &       (lakon(nelem)(2:7).eq.'LPLOID')) then
!     
!     THICK EDGED ORIFICE IN STRAIGHT CONDUIT (L/DH > 0.015)
!     I.E. IDEL'CHIK p175
!     
!     Input parameters
!     
!     Inlet/outlet sections
         a1=prop(index+1)
         a2=prop(index+2)
!     Hydraulic diameter
         dh=prop(index+3)
         if((dh.eq.0.d0).and.(A1.le.A2)) then
            dh=dsqrt(4d0*A1/Pi)
         elseif((dh.eq.0.d0).and.(A1.gt.A2)) then
            dh=dsqrt(4d0*A2/Pi)
         endif
!     Length        
         dl=prop(index+4)
!         
         lzd=dl/dh
         a2za1=min (a1/a2, 1.)
!         
         fa2za1=1.d0-a2za1
!
         iwrit1= 0
         if( lzd.gt.2.4d0) iwrit1= 1
!     
         ldumm=1.D0
         dhdumm=-1.D0
         ks=0.d0
         form_fact=1.d0
!
         call friction_coefficient(ldumm,dhdumm,ks,reynolds,
     &        form_fact,lambda)
!
         call onedint(XLZD,YTOR,n10,lzd,thau,n1,n1,n0,ier)
         zeta0=((0.5d0+thau*dsqrt(fa2za1))+fa2za1) * fa2za1
!
         if(reynolds.gt.1.d+05 ) then
            zeta=zeta0 + lambda * dabs(lzd)
         else
            call onedint(XRE,YERE,n14,reynolds,ereo,n1,n1,n0,ier)
!
            call twodint(zzeta,n15,n11,reynolds,
     &           a2za1,zetap,n1,IEXP,IER)
            zeta=zetap + ereo * zeta0 + lambda * dabs(lzd)
            IF( a2za1.gt.0.95 ) IWRIT1=1
         endif
!     
         if(dabs(lzd) .le. 0.015d0 )then 
            write(*,*) '*WARNING in zeta_calc: L/DH outside valid' 
            write(*,*) '         range ie less than 0.015 !'
         endif
!
         if( iwrit1.eq.1 ) then
            write(*,*) 
     &    'WARNING in zeta_calc: geometry data outside valid range' 
            write(*,*) 
     & '         l/dh greater than 2.4- extrapolated value(s) !'
         endif
!
      elseif((lakon(nelem)(2:7).eq.'REWAOR').or.
     &       (lakon(nelem)(2:7).eq.'LPWAOR'))then
!     
!     THICK-WALLED ORIFICE IN LARGE WALL (L/DH > 0.015)
!     I.E. IDL'CHIK page 174
!
!     Input parameters
!     
!     Inlet/outlet sections
         a1=prop(index+1)
         a2=prop(index+2)
!     Hydraulic diameter
         dh=prop(index+3)
         if((dh.eq.0.d0).and.(A1.le.A2)) then
            dh=dsqrt(4d0*A1/Pi)
         elseif((dh.eq.0.d0).and.(A1.gt.A2)) then
            dh=dsqrt(4d0*A2/Pi)
         endif
!     Length        
         dl=prop(index+4)
!     
         lzd=dl/dh
         ldumm=1.D0
         dhdumm=-1.D0
         ks=0.d0
         form_fact=1.d0
!     
         call friction_coefficient(ldumm,dhdumm,ks,reynolds,
     &        form_fact,lambda)
         call onedint (XLQD,YZETA1,n12,lzd,zeta01,n1,n1,n0,IER)
!     
         iwrit1=0
         if(lzd.gt.10.) iwrit1=1
!     
         if(reynolds.le.1.d+05) then
!     
            call onedint (XRE2,YZETA2,n14,reynolds,zeta02,n1,n1,n10,IER)
            call onedint (XRE2,YERE2,n14,reynolds,EREO,n1,n1,n0,IER)
!     
            zeta=zeta02+0.342*ereo*zeta01+lambda*lzd
!     
         elseif(reynolds.gt.1.d+05) then
            zeta=zeta01+lambda*lzd
         endif
         if(lzd.le.0.015d0) then
            write(*,*) '*WARNING in zeta_calc' 
            write(*,*) 
     &       '         l/dh outside valid range i.e. less than 0.015 !'
         endif
         if(iwrit1.eq.1) then
            write(*,*) '*WARNING in zeta_calc :extrapolated value(s)!'
         endif
!     
         return         
!     
      elseif((lakon(nelem)(2:5).eq.'REEL').or.
     &       (lakon(nelem)(2:5).eq.'LPEL')) then
!     
!     SUDDEN EXPANSION OF A STREAM WITH UNIFORM VELOCITY DISTRIBUTION
!     I.E. IDL'CHIK page 160      
! 
!     Input parameters
!    
!     Inlet/outlet sections
         a1=prop(index+1)
         a2=prop(index+2)
!     
         a2za1=a1/a2
         iwrit1=0
!     
         if(reynolds.LE.10.) then
            zeta=26.0/reynolds
         elseif(reynolds.gt.10.and.reynolds.le.3.5d+03) then
            call twodint(zzeta3,n14,n11,reynolds,a2za1,zeta,n1,IEXP3,
     &               IER)
            if(a2za1.lt.0.01d0.or.a2za1.gt.0.6d0) iwrit1=1
         else
            zeta=(1.-a2za1)**2
         endif
!     
         if(iwrit1.eq.1) then
            write(*,*) '*WARNING in zeta_calc: extrapolated value(s)!'
         endif
         return
!     
      elseif((lakon(nelem)(2:5).eq.'RECO').or.
     &       (lakon(nelem)(2:5).eq.'LPCO'))then
!     
!     SUDDEN CONTRACTION WITH & WITHOUT CONICAL BELLMOUTH ENTRY
!     I.E. IDL'CHIK p 168
! 
!     Input parameters
!    
!     Inlet/outlet sections
         a1=prop(index+1)
         a2=prop(index+2)
!     Hydraulic diameter
         dh=prop(index+3)
         if((dh.eq.0.d0).and.(A1.le.A2)) then
            dh=dsqrt(4d0*A1/Pi)
         elseif((dh.eq.0.d0).and.(A1.gt.A2)) then
            dh=dsqrt(4d0*A2/Pi)
         endif
!     Length
         dl=prop(index+4)
!     Angle
         alpha=prop(index+5)
!     
         a2za1=a2/a1
         iwrit1=0
         dl=abs(dl)
         lzd=dl/dh
!     
         if(dl.eq.0.d0) then
            if(reynolds.le.10.) then
               zeta=27.0/reynolds
            elseif(reynolds.gt.10.and.reynolds.le.1.d+04) then
              call twodint(ZZETA41,n14,n11,reynolds,a2za1,zeta,n1,IEXP,
     &              IER)
               if(a2za1.le.0.1d0.or.a2za1.gt.0.6d0) iwrit1=1
            elseif(reynolds.gt.1.d+04) then
               zeta=0.5d0*(1.-a2za1)
            endif
         elseif(dl.gt.0.d0) then
            call twodint(ZZETA42,n10,n0,alpha,lzd,zeta0,n1,IEXP,IER)
            zeta=zeta0*(1.-a2za1)
            if(lzd.lt.0.025d0 .or. lzd.gt.0.6d0) iwrit1=1
            if(reynolds  .le. 1.d+04) then
               write(*,*) '*WARNING in zeta_calc: reynolds outside valid
     & range i.e. < 10 000 !'
            endif   
         endif
!     
         if( iwrit1.eq.1 ) then
            WRITE(*,*) '*WARNING in zeta_calc: extrapolierte Werte!'
         endif
!     
         return
!
      elseif((lakon(nelem)(2:7).eq.'REBEID').or.
     &       (lakon(nelem)(2:7).eq.'LPBEID')) then
!
!
!        SHARP ELBOW (R/DH=0) AT 0 < DELTA < 180
!        I.E. IDL'CHIK page 294
!     
!        SHARP BENDS 0.5 < R/DH < 1.5 AND 0 < DELTA < 180
!        I.E. IDL'CHIK page 289-290
!
!        SMOOTH BENDS (R/DH > 1.5) AT 0 < DELTA < 180
!        I.E. IDL'CHIK page 289-290
!
!     Input parameters
!     
!     Inlet/outlet sections
         a1=prop(index+1)
         a2=prop(index+2)
!     Hydraulic diameter
         dh=prop(index+3)
         if((dh.eq.0.d0).and.(A1.le.A2)) then
            dh=dsqrt(4d0*A1/Pi)
         elseif((dh.eq.0.d0).and.(A1.gt.A2)) then
            dh=dsqrt(4d0*A2/Pi)
         endif
!     radius
         rad=prop(index+4)
!     angle
         delta=prop(index+5)
!     heigth/width (square section)
         a0=prop(index+6)
         b0=prop(index+7)
!
      iwrit1=0
      iwrit2=0
      rzdh=rad/dh
      if(a0.eq.0.d0)  azb=1.0d0
      if(a0.gt.0.d0) azb=a0/b0
!
      if(rzdh.le.0.5d0) then
         call onedint(XAQB,YC,n12,azb,C,n1,n1,n0,IER)
         zeta1=0.95*(SIN(delta*0.0087))**2+2.05*(SIN(delta*0.0087))**4
         call onedint(XDELTA,YA,n10,delta,A,n1,n1,n10,IER)
         zeta=c*a*zeta1
         if(azb.le.0.25d0.or.azb.gt.8.0d0) iwrit2=1
         if(reynolds.lt.4.d+04) then
            if(reynolds.le.3.d+03) iwrit1=1
            REI=MAX(2999.d0,reynolds)
            ldumm=1.D0
            dhdumm=-1.D0
            ks=0.d0
            form_fact=1.d0
            call friction_coefficient(ldumm,dhdumm,ks,REI,form_fact
     &           ,lambda)
            re_val=4.d+04
            call friction_coefficient(ldumm,dhdumm,ks,re_val,form_fact
     &           , lam)
            zeta=zeta*lambda/lam
         endif
!
      elseif(rzdh.gt.0.5d0.and.rzdh.lt.1.5d0) then
         call onedint(XDELTA,YA1,n10,delta,AI,n1,n1,n10,IER)
         call onedint(XRQDH,YB1,n8,rzdh,B1,n1,n1,n10,IER)
         call onedint(XAQB,YC1,n12,azb,C1,n1,n1,n10,IER)
         REI=MAX(2.d5,reynolds)
         ldumm=1.D0
         dhdumm=-1.D0
         ks=0.d0
         form_fact=1.d0
         call friction_coefficient(ldumm,dhdumm,ks,REI,form_fact
     &        , lambda)
         zeta=AI*B1*C1+0.0175d0*delta*rzdh*lambda
         if(azb.lt.0.25d0.or.azb.gt.8.0d0) iwrit2=1
         if(reynolds.lt.2.d+05) then
            IF(reynolds.lt.3.d+03) iwrit1=1
            REI=MAX(2999.d0,reynolds)
            call friction_coefficient(ldumm,dhdumm,ks,REI,form_fact
     &           ,lambda)
            re_val=2.d+05
            call friction_coefficient(ldumm,dhdumm,ks,re_val,form_fact
     &           , lam)
            zeta=zeta*lambda/lam
         endif
!
      elseif(rzdh.ge.1.5d0.and.rzdh.lt.50.d0) then
         call onedint(XDELTA,YA1,n10,delta,AI,n1,n1,n10,IER)
         call onedint(XAQB,YC2,n12,azb,C2,n1,n1,n10,IER)
         call onedint(XRZDH,YB2,n8,rzdh,B2,n1,n1,n0,IER)
         REI=MAX(2.d5,reynolds)
         ldumm=1.D0
         dhdumm=-1.D0
         ks=0.d0
         form_fact=1.d0
         call friction_coefficient(ldumm,dhdumm,ks,REI,form_fact
     &        ,lambda)
         zeta=AI*B2*C2+0.0175*delta*rzdh*lambda
         if(azb.lt.0.25d0.or.azb.gt.8.0d0) iwrit2=1
         if(reynolds.lt.2.d+05) then
            if(reynolds.lt.3.d+03) iwrit1=1
            REI=MAX(2999.d0,reynolds)
             call friction_coefficient(ldumm,dhdumm,ks,REI,form_fact
     &           ,lambda)
             re_val=2.d+05
            call friction_coefficient(ldumm,dhdumm,ks,re_val,form_fact
     &           , lam)
            zeta=zeta*lambda/lam
         endif
!
      elseif(rzdh.ge.50.d0) then
         zeta=0.0175d0*rzdh*delta*lambda
         if(reynolds.lt.2.d+04) then
             write (*,*)'Reynolds outside valid range i.e. < 20 000!'
         endif
      endif
!
      if(iwrit1.eq.1) then
!     
         write (*,*) 'Reynolds outside valid range i.e. < 3 000!'
      endif
!
      if(iwrit2.eq.1) then
         write(*,*) '*WARNING in zeta_calc: extrapolated value(s)!'
      endif
      return
!
      elseif((lakon(nelem)(2:7).eq.'REBEMI').or.
     &       (lakon(nelem)(2:7).eq.'LPBEMI')) then
!
!     SMOOTH BENDS B.H.R.A HANDBOOK
!
!     Input parameters
!
!     Inlet/outlet sections
         a1=prop(index+1)
         a2=prop(index+2)
!     Hydraulic diameter
         dh=prop(index+3)
!     Radius:
         rad=prop(index+4)
!     angle delta:
         delta=prop(index+5)
!     
         rzdh=Rad/DH
!     
         iwrit1=0
         if( delta.lt.10.d0 .or. delta.gt.180.d0  .or.
     &        rzdh .lt.0.5d0 .or. rzdh.  gt. 10.d0        ) iwrit1=1
!     
         call twodint(ZZETAO,n14,n11,rzdh,delta,zeta0,n1,IEXP6,IER)
         call twodint(KRE, n22,n11,reynolds,rzdh, k,n1,IEXP6,IER)
         zeta=zeta0 * k
!     
         if( reynolds.lt.1.d+3 .or. reynolds.gt.1.d+6 ) then 
            write (*,*)'Reynolds outside valid range <1.E+3 or >1.0E+6'
         endif
!     
         if( iwrit1.eq.1 ) then
            write (*,*)': geometry data outside valid range '
            write (*,*)' - extrapolated value(s)!'
         endif
         RETURN
!
      elseif((lakon(nelem)(2:7).eq.'REBEMA').or.
     &       (lakon(nelem)(2:7).eq.'LPBEMA')) then
!            
!     Own tables and formula to be included
!
         Write(*,*) '*WARNING in zeta_calc: ZETA implicitly equal 1'
         zeta=1.d0
           
      RETURN
!
      elseif((lakon(nelem)(2:5).eq.'REEX').or.
     &       (lakon(nelem)(2:5).eq.'LPEX')) then
!
!     EXIT LOSS COEFFICIENT FOR LAMINAR FLOWS DEPENDING ON THE
!     ACTUAL VELOCITY DISTRIBUTION AT THE EXIT
!
!     Input parameters
!     
!     Inlet/outlet sections
         a1=prop(index+1)
         a2=prop(index+2)
!     Hydraulic diameter
         dh=prop(index+3)
         if((dh.eq.0.d0).and.(A1.le.A2)) then
            dh=dsqrt(4d0*A1/Pi)
         elseif((dh.eq.0.d0).and.(A1.gt.A2)) then
            dh=dsqrt(4d0*A2/Pi)
         endif
!     Reference element
         nelem_ref=nint(prop(index+4))
!
         if(lakon(nelem_ref)(2:5).ne.'GAPF') then
            write(*,*) '*ERROR in zeta_calc :the reference element is no
     &t of type GASPIPE'
           call exit(201)
         endif
!
         if(lakon(nelem_ref)(2:6).eq.'GAPFI') then
            isothermal=.true.
         endif
!     Length of the previous pipe element
         dl=abs(prop(ielprop(nelem_ref)+3))
!    
         if(reynolds .le. 2300.) then
!     (LAMINAR FLOW)
            ldre=dl/dh/reynolds
            call onedint (XDRE,ZETAEX,n12,ldre,zeta,n1,n1,n0,IER)
         elseif((reynolds.gt.2300).and.(reynolds.lt.3000)) then
!     (TRANSITION LAMINAR-TURBULENT)
            ldre=dl/dh/2300.
            call onedint (XDRE,ZETAEX,n12,ldre,zetah,n1,n1,n0,IER)
            zeta=zetah-(zetah-1.)*((reynolds-2300.)/700.)
         else
!     (TURBULENT FLOW, RE.GT.3000)
            zeta=1.
       endif
!     
      RETURN
!
      elseif((lakon(nelem)(2:7).eq.'RELOLI').or.
     &       (lakon(nelem)(2:7).eq.'LPLOLI')) then 
!     
!     'METHOD OF LICHTAROWICZ'
!     "Discharge coeffcients for incompressible non-cavitating 
!     flow through long orifices"
!     A. Lichtarowicz, R.K duggins and E. Markland
!     Journal  Mechanical Engineering Science , vol 7, No. 2, 1965
!
!     TOTAL PRESSURE LOSS COEFFICIENT FOR LONG ORIFICES AND LOW REYNOLDS
!     NUMBERS ( RE < 2.E04 )
!
!     Input parameters
!     
!     Inlet/outlet sections
         a1=prop(index+1)
         a2=prop(index+2)
!     Hydraulic diameter
         dh=prop(index+3)
         if((dh.eq.0.d0).and.(A1.le.A2)) then
            dh=dsqrt(4d0*A1/Pi)
         elseif((dh.eq.0.d0).and.(A1.gt.A2)) then
            dh=dsqrt(4d0*A2/Pi)
         endif
!     Length
         dl=prop(index+4)
!     Isotermal
!
         lzd=dabs(dl)/dh
!     
         cdu=0.827-0.0085*lzd
         km=a1/a2
         call cd_lichtarowicz(cd,cdu,reynolds,km,lzd)
         if(reynolds.gt.2.d04) then
            write(*,*) 
     &        '*WARNING in zeta_calc: range of application exceeded !'
         endif
!     
         zeta=1./cd**2
!     
         return
!     
!     Branch
!     
      elseif((lakon(nelem)(2:5).eq.'REBR').or.
     &       (lakon(nelem)(2:5).eq.'LPBR')) then 
         nelem0=nint(prop(index+1))
         nelem1=nint(prop(index+2))
         nelem2=nint(prop(index+3))
         A0=prop(index+4)
         A1=prop(index+5)
         A2=prop(index+6)
         alpha1=prop(index+7)
         alpha2=prop(index+8)
!     
!     node definition
!     
         node10=kon(ipkon(nelem0)+1)
         node20=kon(ipkon(nelem0)+3)
         nodem0=kon(ipkon(nelem0)+2)
!     
         node11=kon(ipkon(nelem1)+1)
         node21=kon(ipkon(nelem1)+3)
         nodem1=kon(ipkon(nelem1)+2)
!     
         node12=kon(ipkon(nelem2)+1)
         node22=kon(ipkon(nelem2)+3)
         nodem2=kon(ipkon(nelem2)+2)
!     
!     determining the nodes which are not in common
!     
         if(node10.eq.node11) then
            node0=node10 
            node1=node21
            if(node11.eq.node12) then
               node2=node22
            elseif(node11.eq.node22) then
               node2=node12
            endif
         elseif(node10.eq.node21) then
            node0=node10
            node1=node11
            if(node21.eq.node12) then
               node0=node22
            elseif(node21.eq.node22) then
               node2=node12
            endif
         elseif(node20.eq.node11) then
            node0=node20
            node1=node21
            if(node11.eq.node12) then
               node2=node22
            elseif(node11.eq.node22) then
               node2=node12 
            endif
         elseif(node20.eq.node21) then
            node0=node20
            node1=node11
            if(node11.eq.node21) then
               node2=node22
           elseif(node21.eq.node22) then
               node2=node12
            endif
         endif
!     
!     density
!     
         if(lakon(nelem)(2:3).eq.'RE') then
!
!           for gases
!            
            qred_crit=dsqrt(kappa/R)*
     &           (1+0.5d0*(kappa-1))**(-0.5d0*(kappa+1)/(kappa-1))
!     
            icase=0
!     
            Tt0=v(0,node0)
            xflow0=iaxial*v(1,nodem0)
c            xflow0=v(1,nodem0)
            pt0=v(2,node0)
!     
            Qred_0=dabs(xflow0)*dsqrt(Tt0)/(A0*pt0)
            if(Qred_0.gt.qred_crit)
     &           then
               xflow0=qred_crit*(A0*pt0)/dsqrt(Tt0)
            endif
!     
            call ts_calc(xflow0,Tt0,pt0,kappa,r,a0,Ts0,icase)
            M0=dsqrt(2/(kappa-1)*(Tt0/Ts0-1))
!     
            rho0=pt0/(R*Tt0)*(Tt0/Ts0)**(-1/(kappa-1))
!     
            Tt1=v(0,node1)
            xflow1=iaxial*v(1,nodem1)
c            xflow1=v(1,nodem1)
            pt1=v(2,node0)
!     
            Qred_1=dabs(xflow1)*dsqrt(Tt1)/(A1*pt1)
            if(Qred_1.gt.qred_crit)
     &           then
               xflow1=qred_crit*(A1*pt1)/dsqrt(Tt1)
            endif
!     
            call ts_calc(xflow1,Tt1,pt1,kappa,r,a1,Ts1,icase)
            M1=dsqrt(2/(kappa-1)*(Tt1/Ts1-1))
!     
            rho1=pt1/(R*Tt1)*(Tt1/Ts1)**(-1/(kappa-1))
!     
            Tt2=v(0,node2)
            xflow2=iaxial*v(1,nodem2)
c            xflow2=v(1,nodem2)
            pt2=v(2,node0)
!     
            Qred_2=dabs(xflow2)*dsqrt(Tt2)/(A2*pt2)
            if(Qred_2.gt.qred_crit) then
               xflow2=qred_crit*(A2*pt2)/dsqrt(Tt2)
            endif
!     
            call ts_calc(xflow2,Tt2,pt2,kappa,r,a2,Ts2,icase)
            M2=dsqrt(2/(kappa-1)*(Tt2/Ts2-1))
            rho2=pt2/(R*Tt2)*(Tt2/Ts2)**(-1/(kappa-1))
         else
!
!           for liquids the density is supposed to be constant
!           across the element
!
            rho0=1.d0
            rho1=1.d0
            rho2=1.d0
         endif
!     
!     volumic flows (positive)
!     
         V0=dabs(xflow0/rho0)
         V1=dabs(xflow1/rho1)
         V2=dabs(xflow2/rho2)
!
         V1V0=V1/V0
         V2V0=V2/V0
! 
         a0a1=a0/a1
         a0a2=a0/a2
         a2a0=1/a0a2
!     
         W0W1=1/(V1V0*a0a1)
         W0W2=1/(V2V0*a0a2)
!     
!     Branch Joint Genium 
!     Branching Flow Part IV - TEES
!     Fluid Flow Division
!     Section 404.2 page 4 December 1986
!     Genium Publishing (see www.genium.com)
!     
         if((lakon(nelem)(2:7).eq.'REBRJG').or.
     &      (lakon(nelem)(2:7).eq.'LPBRJG')) then
!     
            ang1s=(1.41d0-0.00594*alpha1)*alpha1*pi/180
            ang2s=(1.41d0-0.00594*alpha2)*alpha2*pi/180
!     
            cang1s=dcos(ang1s)
            cang2s=dcos(ang2s)
!     
!     linear part
!     
            zetlin=2.d0*(V1V0**2*a0a1*cang1s+V2V0**2*a0a2*cang2s)
!
            if(nelem.eq.nelem1) then
               call onedint(XANG,YANG,n11,alpha1,lam10,n1,n2,n22,ier)
               zeta=lam10/64*(V1V0*a0a1)**2-zetlin+1d0
               zeta=zeta*(W0W1)**2
!     
            elseif(nelem.eq.nelem2) then
               call onedint(XANG,YANG,n11,alpha2,lam20,n1,n2,n22,ier)
               zeta=lam20/64*(V2V0*a0a2)**2-zetlin+1d0
               zeta=zeta*(W0W2)**2
            endif
            if(zeta.lt.0) then
               write(*,*) '*WARNING in zeta_calc: in Element',nelem,
     &               'TYPE= ',lakon(nelem)(2:5)
               write(*,*) '         Zeta value negative is set to 0.01'
               zeta=0.01d0
            endif
            return
!     
         elseif((lakon(nelem)(2:8).eq.'REBRJI1').or.
     &          (lakon(nelem)(2:8).eq.'LPBRJI1')) then
!
!     Branch Joint Idelchik 1
!     Diagrams of resistance coefficients p260-p266 section VII
!     I.E. IDEL'CHIK 'HANDBOOK OF HYDRAULIC RESISTANCE'
!     2nd edition 1986,HEMISPHERE PUBLISHING CORP.
!     ISBN 0-899116-284-4
!     
            a0a2=a0/a2
            if(alpha2.lt.60.d0) then
               if(nelem.eq.nelem1) then
                  zeta=1.d0-V1V0**2
     &                 -2.d0*a0a2*V2V0**2*dcos(alpha2*pi/180)
                   zeta=zeta*(W0W1)**2
               elseif(nelem.eq.nelem2) then
                  zeta=1.d0-V1V0**2
     &                 -2.d0*a0a2*V2V0**2*dcos(alpha2*pi/180)
     &                 +(a0a2*V2V0)**2-V1V0**2
                  zeta=zeta*(W0W2)**2
               endif
!     
            elseif(alpha2.eq.60.d0) then
!     
!     proceeding as for alpha2<60 with cos(alpha2)=0.5
!     
               if(nelem.eq.nelem1) then
                  zeta=1.d0-V1V0**2-a0a2*V2V0**2
                  zeta=zeta*(W0W1)**2
               elseif(nelem.eq.nelem2) then
                  zeta=1.d0-V1V0**2-a0a2*V2V0**2
     &                 +(a0a2*V2V0)**2-V1V0**2
                  zeta=zeta*(W0W2)**2
               endif
!     
            elseif(alpha2.lt.90.d0) then
!     
!     linear interpolation between alpha2=60 and alpha2=90
!     
               z1_60=1.d0-V1V0**2-a0a2*V2V0**2
               z1_90=(1.55d0-V2V0)*V2V0
               if(nelem.eq.nelem1) then
                  zeta=z1_60+(z1_90-z1_60)*(alpha2-60.d0)/30
                  zeta=zeta*(W0W1)**2
               elseif(nelem.eq.nelem2) then
                  z2_60=z1_60+(a0a2*V2V0)**2-V1V0**2
                  call onedint(TA2A0,TAFAKT,n12,a2a0,afakt,
     &                 n1,n1,n11,ier)
                  z2_90=afakt*(1.d0+(a0a2*V2V0)**2-2.d0*V1V0**2)
                  zeta=z2_60+(z2_90-z2_60)*(alpha2-60.d0)/30d0
                  zeta=zeta*(W0W2)**2
               endif
!     
            elseif(alpha2.eq.90.d0) then
               if(nelem.eq.nelem1) then
                  zeta=(1.55d0-V2V0)*V2V0
                  zeta=zeta*(W0W1)**2
               elseif(nelem.eq.nelem2) then
                  call onedint(TA2A0,TAFAKT,n12,a2a0,afakt,
     &                 n1,n1,n11,ier) 
                  zeta=afakt*(1.d0+(a0a2*V2V0)**2-2.d0*V1V0**2)
                  zeta=zeta*(W0W2)**2
               endif
            endif
            if(zeta.lt.0) then
               write(*,*) '*WARNING in zeta_calc: in Element',nelem,
     &               'TYPE= ',lakon(nelem)(2:5)
               write(*,*) '         Zeta value negative is set to 0.01'
               zeta=0.01d0
            endif
            return
!     
         elseif((lakon(nelem)(2:8).eq.'REBRJI2').or.
     &          (lakon(nelem)(2:8).eq.'LPBRJI2')) then
!     
!     Branch Joint Idelchik 2
!     Diagrams of resistance coefficients page 348-352 
!          I.E. IDEL'CHIK 'HANDBOOK OF HYDRAULIC RESISTANCE'
!     2nd edition 1986,HEMISPHERE PUBLISHING CORP.
!     ISBN 0-899116-284-4 page 348-352 
!     
            if(alpha2.lt.60.d0) then
               if(nelem.eq.nelem1) then
                  zeta=1+a0a1*V1V0**2*(a0a1-2.)
     &                 -2d0*a0a2*V2V0**2*dcos(alpha2*pi/180)
!     correction term
                  call twodint(KSTAB,n6,n11,a2a0,alpha2,ks2,n1
     &                 ,iexpbr1,ier)
                  zeta=zeta+ks2
                  zeta=zeta*(W0W1)**2
               elseif(nelem.eq.nelem2) then
                  zeta=1+a0a1*V1V0**2*(a0a1-2.)
     &                 -2d0*a0a2*V2V0**2*dcos(alpha2*pi/180)
     &                 -(a0a1*V1V0)**2+(a0a2*V2V0)**2
                  call twodint(KBTAB,n6,n11,a2a0,alpha2,kb,n1,
     &                 iexpbr1,ier)
                  zeta=zeta+kb
                  zeta=zeta*(W0W2)**2
               endif
!     
            elseif(alpha2.eq.60.d0) then
!     as for alpha2 < 60 , with dcos(alpha2)=0.5
               if(nelem.eq.nelem1) then
                  zeta=1+a0a1*V1V0**2*(a0a1-2.)-a0a2*V2V0**2
                  call twodint(KSTAB,n6,n11,a2a0,alpha2,ks2,n1,
     &                 iexpbr1,ier)
                  zeta=zeta+ks2
                  zeta=zeta*(W0W1)**2
               elseif(nelem.eq.nelem2) then
                  zeta=1+a0a1*V1V0**2*(a0a1-2.)-a0a2*V2V0**2
     &                 -(a0a1*V1V0)**2+(a0a2*V2V0)**2
                  call twodint(KBTAB,n6,n11,a2a0,alpha2,kb,n1,
     &                 iexpbr1,ier)
                  zeta=zeta+kb
                  zeta=zeta*(W0W2)**2
               endif
!     
            elseif(alpha2.lt.90.d0) then
!     linear interpolation between alpha2=60 and alpha2=90
               z1_60=1+a0a1*V1V0**2*(a0a1-2.)-a0a2*V2V0**2
!     correction term
               call twodint(KSTAB,n6,n11,a2a0,alpha2,ks2,n1,
     &              iexpbr1,ier)
               z1_60=z1_60+ks2
               if(nelem.eq.nelem1) then
                  call twodint(Z90TAB,n6,n11,a2a0,V2V0,z1_90,
     &                 n1,iexpbr1,ier)
                  zeta=z1_60+(z1_90-z1_60)*(alpha2-60)/30
                  zeta=zeta*(W0W1)**2
               elseif(nelem.eq.nelem2) then
                  z2_60=z1_60-(a0a1*V1V0)**2+(a0a2*v2v0)**2
                  call twodint(KBTAB,n6,n11,a2a0,alpha2,kb,n1,
     &                 iexpbr1,ier)
                  z2_60=z2_60+kb-ks2
                  z2_90=1.+(a0a2*V2V0)**2-2*a0a1*V1V0**2+kb
                  zeta=z2_60+(z2_90-z2_60)*(alpha2-60)/30
                  zeta=zeta*(W0W2)**2
               endif
            elseif(alpha2.eq.90.d0) then
               if(nelem.eq.nelem2) then
                  call twodint(KBTAB,n6,n11,a2a0,alpha2,kb,n1,
     &                 iexpbr1,ier)
                  zeta=1.+(a0a2*V2V0)**2-2*a0a1*V1V0**2+kb
                  zeta=zeta*(W0W2)**2
               elseif(nelem.eq.nelem1) then
!     table interpolation
                  call twodint(Z90TAB,n6,n11,a2a0,V2V0,zeta,
     &                 n1,iexpbr1,ier)
                  zeta=zeta*(W0W1)**2
!     cheching whether the table eveluation in the eptrapolated domain
!     (This procedure is guessed from the original table)
!     
                  Z90LIM11=Z90LIMX(1)
                  Z90LIM51=Z90LIMX(5)
                  if((a2a0.ge.Z90LIM11)
     &                 .and.(a2a0.le.Z90LIM51))then
                     call onedint(Z90LIMX,Z90LIMY,n5,A2A0,
     &                    V2V0L,n1,n1,n11,ier)
                     if(V2V0.gt.V2V0L) then
                        write(*,*) 'WARNING in zeta_calc: in element',
     &                                 nelem
                        write(*,*) 
     &            '        V2V0 in the extrapolated domain'
                        write(*,*) '        for zeta table (branch 1)'
                     endif
                  endif
               endif
            endif
            if(zeta.lt.0) then
               write(*,*) '*WARNING in zeta_calc: in Element',nelem,
     &               'TYPE= ',lakon(nelem)(2:5)
               write(*,*) '         Zeta value negative is set to 0.01'
               zeta=0.01d0
            endif
            return
!
         elseif((lakon(nelem)(2:7).eq.'REBRSG').or.
     &          (lakon(nelem)(2:7).eq.'LPBRSG')) then
!
!     Branch Split Genium 
!     Branching Flow Part IV - TEES
!     Fluid Flow Division
!     Section 404.2 page 3 December 1986
!     Genium Publishing (see www.genium.com)
!          
            if(nelem.eq.nelem1) then
!     
               ang1s=(1.41d0-0.00594*alpha1)*alpha1*pi/180
!     
               cang1s=dcos(ang1s)
!     
               if(alpha1.le.22.5) then
                  lam11=0.0712*alpha1**0.7041+0.37
                  lam12=0.0592*alpha1**0.7029+0.37
               else
                  lam11=1.d0
                  lam12=0.9d0
               endif
               zeta=lam11+(2.d0*lam12-lam11)*(V1V0*a0a1)**2
     &              -2d0*lam12*V1V0*a0a1*cang1s
               zeta=zeta*(W0W1)**2
!     
            elseif(nelem.eq.nelem2) then
!     
               ang2s=(1.41d0-0.00594*alpha2)*alpha2*pi/180
! 
               cang2s=dcos(ang2s)
!     
               if(alpha2.le.22.5) then
                  lam21=0.0712*alpha2**0.7041+0.37
                  lam22=0.0592*alpha2**0.7029+0.37
               else
                  lam21=1.d0
                  lam22=0.9d0
               endif
!     
               zeta=lam21+(2.d0*lam22-lam21)*(V2V0*a0a2)**2
     &              -2d0*lam22*V2V0*a0a2*cang2s
               zeta=zeta*(W0W2)**2
!     
            endif
            return 
!     
         elseif((lakon(nelem)(2:8).eq.'REBRSI1').or.
     &          (lakon(nelem)(2:8).eq.'LPBRSI1'))  then 
!
!     Branch Split Idelchik 1
!     Diagrams of resistance coefficients p280,p282 section VII
!     I.E. IDEL'CHIK 'HANDBOOK OF HYDRAULIC RESISTANCE'
!     2nd edition 1986,HEMISPHERE PUBLISHING CORP.
!     ISBN 0-899116-284-4
!          
            W1W0=V1V0*a0a1
            W2W0=V2V0*a0a2
!     
            if(nelem.eq.nelem1) then
               zeta=0.4d0*(1-W1W0)**2
               zeta=zeta*(W0W1)**2
!
            elseif(nelem.eq.nelem2) then
!
               dh0=dsqrt(A0*4d0/Pi)
               if(dh0.eq.0.d0) then
                  dh0=dsqrt(4d0*A0/Pi)
               endif
               dh2=dsqrt(A2*4d0/Pi)
               if(dh2.eq.0.d0) then
                  dh2=dsqrt(4d0*A2/Pi)
               endif
!
               hq=dh2/dh0
               if(alpha2.le.60.or.hq.le.2.d0/3.d0) then
                  zeta=0.95d0*((W2W0-2d0*dcos(alpha2*pi/180))
     &                 *W2W0+1.d0)
                  zeta=zeta*(W0W2)**2
               else
                  z2d390=0.95d0*((W2W0-2d0*dcos(90.d0*pi/180))
     &                 *W2W0+1.d0)
                  z1p090=0.95*(0.34d0+W2W0**2)
                  z90=z2d390+(3*hq-2.d0)*(z1p090-z2d390)
                  Z60=0.95d0*((W2W0-2d0*dcos(60.d0*pi/180))
     &                 *W2W0+1.d0)
                  zeta=z60+(alpha2/30.d0-2.d0)*(z90-z60)
                  zeta=zeta*(W0W2)**2
               endif
            endif
            return
!                 
         elseif((lakon(nelem)(2:8).eq.'REBRSI2').or.
     &          (lakon(nelem)(2:8).eq.'LPBRSI2')) then 
!
!     Branch Split Idelchik 2
!     Diagrams of resistance coefficients p289,section VII
!     I.E. IDEL'CHIK 'HANDBOOK OF HYDRAULIC RESISTANCE'
!     2nd edition 1986,HEMISPHERE PUBLISHING CORP.
!     ISBN 0-899116-284-4
!     
            if(nelem.eq.nelem1) then
               W1W0=V1V0*a0a1
               W0W1=1/W1W0
               zeta=1.d0+0.3d0*W1W0**2
               zeta=zeta*(W0W1)**2
            elseif(nelem.eq.nelem2) then
               W2W0=V2V0*a0a2
               W0W2=1/W2W0
               zeta=1.d0+0.3d0*W2W0**2
               zeta=zeta*(W0W2)**2
            endif
            return
         endif
         if(zeta.lt.0) then
            write(*,*) '*WARNING in zeta_calc: in Element',nelem,
     &           'TYPE= ',lakon(nelem)(2:5)
            write(*,*) '         Zeta value negative is set to 0.01'
            zeta=0.01d0
         endif
      endif
!   
      return
      end
      
      
