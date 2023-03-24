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
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine rcavi2(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,mi,ttime,time,
     &     iaxial,iplausi)
!     
!     rotating cavity element
!
!     author: Yannick Muller
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(*),idirf(*),index,kflag,mi(*),
     &     inv,ipkon(*),kon(*),kgas,nelem_in,nelem_out,iaxial,
     &     element0,node10,node20,node11,node21,node12,node22,node_cav,
     &     node_main,node_main2,node_in1,node_out1,node_in2,node_out2,
     &     iplausi
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),kappa,R,a,d,
     &     p1,p2,T1,T2,Aeff,C1,C2,C3,cd,cp,physcon(*),p2p1,km1,dvi,
     &     kp1,kdkm1,tdkp1,km1dk,x,y,ca1,cb1,ca2,cb2,dT1,alambda,
     &     reynolds,pi,xflow_oil,s,Tcav,pcav,pmin,pmax,ttime,time,
     &     Tref,Alpha1, Alpha2, Alpha3, GF,kf,MRTAP_ref_ein,
     &     MRTAP_ref_aus, m_ref_ein, m_ref_aus,maus_zu_mref,
     &     mein_zu_mref, A_aus, A_ein, A_ges,m_aus, m_ein, m_sperr
!
!
!
      pi=4.d0*datan(1.d0)   
!
      if (kflag.eq.0) then
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
!     
      elseif (kflag.eq.1) then      
         if(v(1,nodem).ne.0.d0) return
!     
         p1=v(2,node1)
         call rcavi_cp_lt(xflow)
         
         call rcavi_cp_nt(xflow)
      elseif (kflag.eq.2) then
!
      elseif (kflag.eq.3) then
!
      endif
      return
      end
