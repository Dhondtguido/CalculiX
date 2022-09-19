!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
!
!  regularization function for normal contact mortar
! see phd-thesis Sitzmann Chapter 3.2.1., semi-smooth Newton for normal contact 
!
!  [in]lambdap  contact pressure in normal direction
!  [in]divmode indicates whether function or derivate 
!                             should be called
!                    =0 function called
!                    =1 derivative called    
!  [in]regmode        selects regularization funtion
!                    =1 perturbed Lagrange
!                    =2 piece wise linear with given data points
!                    =3 expontential contact law
!  [out]gnc        result regularization function
!  [in] aninvloc        stiffness constant for perturbed Lagrange
!  [in] p0 parameter for exponential regularization
!  [in] beta parameter for exponential regularization
!
      subroutine regularization_gn_c(lambdap,divmode,regmode,
     &     gnc,aninvloc,p0,beta,elcon,nelcon,itie,ntmat_,
     &     plicon,nplicon,npmat_,ncmat_,tietol,scal)
!     
!     regularization function for normal contact
!     Author: Saskia Sitzmann
!     
      implicit none
!     
      integer divmode,regmode,i,ndata,kode,npmat_,ncmat_,
     &     itie,ntmat_,nplicon(0:ntmat_,*),nelcon(2,*),
     &     imat
!     
      real*8 lambdap,gnc,pn_d(40),gn_d(40),aninv1(40),t(40),
     &     beta,p0,aninvloc,elconloc(ncmat_),plconloc(802),t1l,
     &     elcon(0:ncmat_,ntmat_,*),plicon(0:2*npmat_,ntmat_,*),
     &     tietol(4,*),scal
!      
      kode=-51
      t1l=0.0
      imat=int(tietol(2,itie+1))
!     
      gnc=0.0
!
!     perturbed Lagrange
!
      if(regmode.eq.1)then
         if(divmode.eq.0)then
            gnc=aninvloc*lambdap
         elseif(divmode.eq.1)then
            gnc=aninvloc
         else
            write(*,*)'error in regularzation_gn_c.f!'
            call exit(201)
         endif
!     
!     multiple perturbed Lagrange
!
      else if(regmode.eq.2)then
         call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &        elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
         ndata=int(plconloc(801))
!        is this right?
c         ndata=int(plconloc(81))
         do i=1,ndata
            gn_d(i)=plconloc(2*i-1)*scal
            pn_d(i)=plconloc(2*i)*scal
         enddo 
         do i=1,ndata-1
            aninv1(i)=(gn_d(i+1)-gn_d(i))/(pn_d(i+1)-pn_d(i))
            t(i)=gn_d(i)-aninv1(i)*pn_d(i)
         enddo
         if(divmode.eq.0)then
            if(lambdap.lt.pn_d(1))then
               gnc=aninv1(1)*lambdap+t(1)
            endif
            do i=1,ndata-2
               if(pn_d(i+1).gt.lambdap.and.pn_d(i).le.lambdap)then
                  gnc=aninv1(i)*lambdap+t(i)
               endif
            enddo
            if(pn_d(ndata-1).le.lambdap)then
               gnc=aninv1(ndata-1)*lambdap+t(ndata-1)
            endif
         elseif(divmode.eq.1)then
            if(lambdap.lt.pn_d(1))then
               gnc=aninv1(1)
            endif
            do i=1,ndata-2
               if(pn_d(i+1).gt.lambdap.and.pn_d(i).le.lambdap)then
                  gnc=aninv1(i)
               endif
            enddo
            if(pn_d(ndata-1).le.lambdap)then
               gnc=aninv1(ndata-1)
            endif
         else
            write(*,*)'error in regularzation_gn_c.f!'
            call exit(201)
         endif
!     
!     exponetial perturbed Lagrange
!
      else if (regmode.eq.3)then
         if(divmode.eq.0)then
            if(lambdap.gt.0.0)then
               gnc=scal*beta*log(((lambdap/scal)+p0)/p0)
            else
               gnc=beta/(p0)*lambdap
            endif
         elseif(divmode.eq.1)then
            if(lambdap.gt.0.0)then
               gnc=beta/(lambdap+p0)
            else
               gnc=beta/(p0)
            endif
         else
            write(*,*)'error in regularzation_gn_c.f!'
            call exit(201)
         endif
      endif
c     call exit(201)
      return
      end
