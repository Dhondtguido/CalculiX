!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine us3_materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &  nalcon,imat,amat,iorien,pgauss,orab,ntmat_,elas,rho,iel,
     &  ithermal,alzero,mattyp,t0l,t1l,ihyper,istiff,elconloc,eth,
     &  kode,plicon,nplicon,plkcon,nplkcon,npmat_,plconloc,mi,dtime,
     &  iint,xstiff,ncmat_)
!
      implicit none
!
!     determines the material data for element iel
!     3-node shell element
!     author: Gil Rama
!
!     istiff=0: only interpolation of material data
!     istiff=1: copy the consistent tangent matrix from the field
!               xstiff and check for zero entries
!
      character*80 amat
!
      integer nelcon(2,*),nrhcon(*),nalcon(2,*),
     &  imat,iorien,ithermal(*),j,k,mattyp,kal(2,6),j1,j2,j3,j4,
     &  jj,ntmat_,istiff,nelconst,ihyper,kode,itemp,kin,nelas,
     &  iel,iint,mi(*),ncmat_,id,two,seven,
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_
!
      real*8 elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  alcon(0:6,ntmat_,*),eth(6),xstiff(27,mi(1),*),
     &  orab(7,*),elas(21),alph(6),alzero(*),rho,t0l,t1l,
     &  skl(3,3),xa(3,3),elconloc(*),emax,pgauss(3),
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  plconloc(802),dtime
!
!
!
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!
      two=2
      seven=7
!
!     nelconst: # constants read from file
!     nelas: # constants in the local tangent stiffness matrix
!
      if(istiff.eq.1) then

         nelas=nelcon(1,imat)
         if((nelas.lt.0).or.((nelas.ne.2).and.(iorien.ne.0))) nelas=21
!
!        calculating the density (needed for the mass matrix and
!        gravity or centrifugal loading)
!
         if(ithermal(1).eq.0) then
            rho=rhcon(1,1,imat)
         else
            call ident2(rhcon(0,1,imat),t1l,nrhcon(imat),two,id)
            if(nrhcon(imat).eq.0) then
               continue
            elseif(nrhcon(imat).eq.1) then
               rho=rhcon(1,1,imat)
            elseif(id.eq.0) then
               rho=rhcon(1,1,imat)
            elseif(id.eq.nrhcon(imat)) then
               rho=rhcon(1,id,imat)
            else
               rho=rhcon(1,id,imat)+
     &              (rhcon(1,id+1,imat)-rhcon(1,id,imat))*
     &              (t1l-rhcon(0,id,imat))/
     &              (rhcon(0,id+1,imat)-rhcon(0,id,imat))
            endif
         endif
!     
!        for nonlinear behavior (nonlinear geometric or
!        nonlinear material behavior): copy the stiffness matrix
!        from the last stress calculation
!     
         do j=1,21
            elas(j)=xstiff(j,iint,iel)
         enddo
!
!        check whether the fully anisotropic case can be
!        considered as orthotropic         
!
         if(nelas.eq.21) then
            emax=0.d0
            do j=1,9
               emax=max(emax,dabs(elas(j)))
            enddo
            do j=10,21
               if(dabs(elas(j)).gt.emax*1.d-10) then
                  emax=-1.d0
                  exit
               endif
            enddo
            if(emax.gt.0.d0) nelas=9
         endif
!
!        determining the type: isotropic, orthotropic or anisotropic         
!
         if(nelas.le.2) then
            mattyp=1
         elseif(nelas.le.9) then
            mattyp=2
         else
            mattyp=3
         endif
!
      else
!
         nelconst=nelcon(1,imat)
!
         if(nelconst.lt.0) then
!
!           inelastic material or user material
!
            if(nelconst.eq.-1) then
               nelconst=3
            elseif(nelconst.eq.-2) then
               nelconst=3
            elseif(nelconst.eq.-3) then
               nelconst=2
            elseif(nelconst.eq.-4) then
               nelconst=3
            elseif(nelconst.eq.-5) then
               nelconst=6
            elseif(nelconst.eq.-6) then
               nelconst=9
            elseif(nelconst.eq.-7) then
               nelconst=3
            elseif(nelconst.eq.-8) then
               nelconst=7
            elseif(nelconst.eq.-9) then
               nelconst=12
            elseif(nelconst.eq.-10) then
               nelconst=2
            elseif(nelconst.eq.-11) then
               nelconst=4
            elseif(nelconst.eq.-12) then
               nelconst=6
            elseif(nelconst.eq.-13) then
               nelconst=5
            elseif(nelconst.eq.-14) then
               nelconst=6
            elseif(nelconst.eq.-15) then
               nelconst=3
            elseif(nelconst.eq.-16) then
               nelconst=6
            elseif(nelconst.eq.-17) then
               nelconst=9
            elseif(nelconst.eq.-50) then
               nelconst=5
            elseif(nelconst.eq.-51) then
               nelconst=2
            elseif(nelconst.eq.-52) then
               nelconst=5
            elseif(nelconst.le.-100) then
               nelconst=-nelconst-100
            endif
!     
         endif
!
!        in case no initial temperatures are defined, the calculation
!        is assumed athermal, and the first available set material
!        constants are used
!
         if(ithermal(1).eq.0) then
            if(ihyper.ne.1) then
               do k=1,nelconst
                  elconloc(k)=elcon(k,1,imat)
               enddo
            else
               do k=1,nelconst
                  elconloc(k)=elcon(k,1,imat)
               enddo
!     
               itemp=1
!     
               if((kode.lt.-50).and.(kode.gt.-100)) then
                  plconloc(1)=0.d0
                  plconloc(2)=0.d0
                  plconloc(3)=0.d0
                  plconloc(801)=nplicon(1,imat)+0.5d0
                  plconloc(802)=nplkcon(1,imat)+0.5d0
!     
!     isotropic hardening
!     
                  if(nplicon(1,imat).ne.0) then
                     kin=0
                     call plcopy(plicon,nplicon,plconloc,npmat_,ntmat_,
     &                    imat,itemp,iel,kin)
                  endif
!     
!     kinematic hardening
!     
                  if(nplkcon(1,imat).ne.0) then
                     kin=1
                     call plcopy(plkcon,nplkcon,plconloc,npmat_,ntmat_,
     &                    imat,itemp,iel,kin)
                  endif
!     
               endif
!     
            endif
         else
!     
!     calculating the expansion coefficients
! 
!     
            
            call ident2(alcon(0,1,imat),t1l,nalcon(2,imat),seven,id)
            if(nalcon(2,imat).eq.0) then
               do k=1,6
                  alph(k)=0.d0
               enddo
               continue
            elseif(nalcon(2,imat).eq.1) then
               do k=1,nalcon(1,imat)
                  alph(k)=alcon(k,1,imat)!*(t1l-alzero(imat))
               enddo
            elseif(id.eq.0) then
               do k=1,nalcon(1,imat)
                  alph(k)=alcon(k,1,imat)!*(t1l-alzero(imat))
               enddo
            elseif(id.eq.nalcon(2,imat)) then
               do k=1,nalcon(1,imat)
                  alph(k)=alcon(k,id,imat)!*(t1l-alzero(imat))
               enddo
            else
               do k=1,nalcon(1,imat)
                  alph(k)=(alcon(k,id,imat)+
     &                 (alcon(k,id+1,imat)-alcon(k,id,imat))*
     &                 (t1l-alcon(0,id,imat))/
     &                 (alcon(0,id+1,imat)-alcon(0,id,imat)))
!     &                 *(t1l-alzero(imat))
               enddo
            endif
!     
!     subtracting the initial temperature influence       
!     
!             call ident2(alcon(0,1,imat),t0l,nalcon(2,imat),seven,id)
!             if(nalcon(2,imat).eq.0) then
!                continue
!             elseif(nalcon(2,imat).eq.1) then
!                do k=1,nalcon(1,imat)
!                   alph(k)=alph(k)-alcon(k,1,imat)*(t0l-alzero(imat))
!                enddo
!             elseif(id.eq.0) then
!                do k=1,nalcon(1,imat)
!                   alph(k)=alph(k)-alcon(k,1,imat)*(t0l-alzero(imat))
!                enddo
!             elseif(id.eq.nalcon(2,imat)) then
!                do k=1,nalcon(1,imat)
!                   alph(k)=alph(k)-alcon(k,id,imat)*(t0l-alzero(imat))
!                enddo
!             else
!                do k=1,nalcon(1,imat)
!                   alph(k)=alph(k)-(alcon(k,id,imat)+
!      &                 (alcon(k,id+1,imat)-alcon(k,id,imat))*
!      &                 (t0l-alcon(0,id,imat))/
!      &                 (alcon(0,id+1,imat)-alcon(0,id,imat)))
!      &                 *(t0l-alzero(imat))
!                enddo
!             endif
!     
!           storing the thermal strains
!     
            if(nalcon(1,imat).eq.1) then
               do k=1,3
                  eth(k)=alph(1)
               enddo
               do k=4,6
                  eth(k)=0.d0
               enddo
            elseif(nalcon(1,imat).eq.3) then
               do k=1,3
                  eth(k)=alph(k)
               enddo
               do k=4,6
                  eth(k)=0.d0
               enddo
            else
               do k=1,6
                  eth(k)=alph(k)
               enddo
            endif
!     
!     calculating the hardening coefficients
!     
!     for the calculation of the stresses, the whole curve
!     has to be stored:
!     plconloc(2*k-1), k=1...20: equivalent plastic strain values (iso)
!     plconloc(2*k),k=1...20:    corresponding stresses (iso)
!     plconloc(39+2*k),k=1...20: equivalent plastic strain values (kin)
!     plconloc(40+2*k),k=1...20: corresponding stresses (kin)
!     
!     initialization
!     
            if((kode.lt.-50).and.(kode.gt.-100)) then
               if(npmat_.eq.0) then
                  plconloc(801)=0.5d0
                  plconloc(802)=0.5d0
               else
                  plconloc(1)=0.d0
                  plconloc(2)=0.d0
                  plconloc(3)=0.d0
                  plconloc(801)=nplicon(1,imat)+0.5d0
                  plconloc(802)=nplkcon(1,imat)+0.5d0
!     
!     isotropic hardening
!     
                  if(nplicon(1,imat).ne.0) then
!     
                     if(nplicon(0,imat).eq.1) then
                        id=-1
                     else
                        call ident2(plicon(0,1,imat),t1l,
     &                       nplicon(0,imat),2*npmat_+1,id)
                     endif
!     
                     if(nplicon(0,imat).eq.0) then
                        continue
                     elseif((nplicon(0,imat).eq.1).or.(id.eq.0).or.
     &                       (id.eq.nplicon(0,imat))) then
                        if(id.le.0) then
                           itemp=1
                        else
                           itemp=id
                        endif
                        kin=0
                        call plcopy(plicon,nplicon,plconloc,npmat_,
     &                       ntmat_,imat,itemp,iel,kin)
                        if((id.eq.0).or.(id.eq.nplicon(0,imat))) then
                        endif
                     else
                        kin=0
                        call plmix(plicon,nplicon,plconloc,npmat_,
     &                       ntmat_,imat,id+1,t1l,iel,kin)
                     endif
                  endif
!     
!     kinematic hardening
!     
                  if(nplkcon(1,imat).ne.0) then
!     
                     if(nplkcon(0,imat).eq.1) then
                        id=-1
                     else
                        call ident2(plkcon(0,1,imat),t1l,
     &                       nplkcon(0,imat),2*npmat_+1,id)
                     endif
!     
                     if(nplkcon(0,imat).eq.0) then
                        continue
                     elseif((nplkcon(0,imat).eq.1).or.(id.eq.0).or.
     &                       (id.eq.nplkcon(0,imat))) then
                        if(id.le.0)then
                           itemp=1
                        else
                           itemp=id
                        endif
                        kin=1
                        call plcopy(plkcon,nplkcon,plconloc,npmat_,
     &                       ntmat_,imat,itemp,iel,kin)
                        if((id.eq.0).or.(id.eq.nplkcon(0,imat))) then
                        endif
                     else
                        kin=1
                        call plmix(plkcon,nplkcon,plconloc,npmat_,
     &                       ntmat_,imat,id+1,t1l,iel,kin)
                     endif
                  endif
               endif
            endif
!     
!     calculating the elastic constants
!     
            call ident2(elcon(0,1,imat),t1l,nelcon(2,imat),ncmat_+1,id)
            if(nelcon(2,imat).eq.0) then
               continue
            elseif(nelcon(2,imat).eq.1) then
                  do k=1,nelconst
                     elconloc(k)=elcon(k,1,imat)
                  enddo
            elseif(id.eq.0) then
                  do k=1,nelconst
                     elconloc(k)=elcon(k,1,imat)
                  enddo
            elseif(id.eq.nelcon(2,imat)) then
                  do k=1,nelconst
                     elconloc(k)=elcon(k,id,imat)
                  enddo
            else
                  do k=1,nelconst
                     elconloc(k)=elcon(k,id,imat)+
     &                    (elcon(k,id+1,imat)-elcon(k,id,imat))*
     &                    (t1l-elcon(0,id,imat))/
     &                    (elcon(0,id+1,imat)-elcon(0,id,imat))
                  enddo
            endif
!
!           modifying the thermal constants if anisotropic and
!           a transformation was defined
!
            if((iorien.ne.0).and.(nalcon(1,imat).gt.1)) then
!     
!              calculating the transformation matrix
!     
               call transformatrix(orab(1,iorien),pgauss,skl)
!     
!              transforming the thermal strain
!     
               xa(1,1)=eth(1)
               xa(1,2)=eth(4)
               xa(1,3)=eth(5)
               xa(2,1)=eth(4)
               xa(2,2)=eth(2)
               xa(2,3)=eth(6)
               xa(3,1)=eth(5)
               xa(3,2)=eth(6)
               xa(3,3)=eth(3)
!     
               do jj=1,6
                  eth(jj)=0.d0
                  j1=kal(1,jj)
                  j2=kal(2,jj)
                  do j3=1,3
                     do j4=1,3
                        eth(jj)=eth(jj)+
     &                       xa(j3,j4)*skl(j1,j3)*skl(j2,j4)
                     enddo
                  enddo
               enddo
            endif
         endif
      endif
!
      return
      end

      !
      SUBROUTINE us3_csys_cr(xg,tm,tmg)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)  :: xg(3,3)         ! element coordinates in global csys
      REAL*8, INTENT(OUT) :: tm(3,3),tmg(18,18)         ! transformation matrix (e1,e2,e3)T
      REAL*8 :: e1(3),e2(3),e3(3), dl
      !
      INTEGER :: j,i
      !
      tm(:,:) = 0.d0
      !    element frame (e1,e2,e3)
      !     e1 = 1 -> 2
      do j = 1,3
       e1(j) = xg(2,j)-xg(1,j)
      enddo
      ! norm it
      dl = dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
      do j = 1,3
         e1(j) = e1(j)/dl
      enddo 
      !
      !     e2 = 2 -> 3 (temp)
      !     
      do j = 1,3
        e2(j) = xg(3,j)-xg(1,j)
      enddo
      !
      !     e3 = e1 x e2
      !
      e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
      e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
      e3(3) = e1(1)*e2(2) - e1(2)*e2(1)
      ! norm it
      dl = dsqrt(e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3))
      do j = 1,3
         e3(j) = e3(j)/dl
      enddo 
      !
      !     e2 = e3 x e1
      !
      e2(:) = 0.d0
      e2(1) = e3(2)*e1(3) - e3(3)*e1(2)
      e2(2) = e3(3)*e1(1) - e3(1)*e1(3)
      e2(3) = e3(1)*e1(2) - e3(2)*e1(1)
      ! norm it
      dl = dsqrt(e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3))
      do j = 1,3
         e2(j) = e2(j)/dl
      enddo         
      !     transformation matrix from the global into the local system
      do j = 1,3
         tm(1,j) = e1(j)
         tm(2,j) = e2(j)
         tm(3,j) = e3(j)
      enddo
      !
      tmg(:,:) = 0.d0
      do i = 1,3
        do j = 1,3
          tmg(i,j)       = tm(i,j) ! f1
          tmg(i+3,j+3)   = tm(i,j) ! r1
          tmg(i+6,j+6)   = tm(i,j) ! f2
          tmg(i+9,j+9)   = tm(i,j) ! r2
          tmg(i+12,j+12) = tm(i,j) ! f3
          tmg(i+15,j+15) = tm(i,j) ! r3
        enddo
      enddo
      !
      RETURN
      END  
      !
      SUBROUTINE us3_csys(xg,tm,tmg) 
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)  :: xg(3,3)         ! element coordinates in global csys
      REAL*8, INTENT(OUT) :: tm(3,3),tmg(18,18)         ! transformation matrix (e1,e2,e3)T
      REAL*8 :: e1(3),e2(3),e3(3),dl,dd,xno(3)
      REAL*8 :: xi,et,xl(3,8),xs(3,7),p(3),shp(7,8)
      REAL*8 :: a(3,3)
      !
      INTEGER :: iflag,j,i,l,k
      !
      tm(:,:) = 0.d0
      !
      do i=1,3
        do j=1,3
      	 xl(i,j) = xg(j,i)
        enddo
      enddo
      xi = 0.d0
      et = 0.d0
      iflag = 2
      call shape3tri(xi,et,xl,xno,xs,shp,iflag)
      dd = dsqrt(xno(1)*xno(1)+xno(2)*xno(2)+xno(3)*xno(3))
      do l = 1,3
        xno(l)=xno(l)/dd
      enddo
      ! coordinates of the point at (xi,et)=(0.,0.)      
      do l = 1,3
        p(l) = 0.d0
        do k = 1,3
          p(l) = p(l) + shp(4,k)*xl(l,k)
        enddo
      enddo      
      ! unit matrix
      do k = 1,3
        do l = 1,3
          a(k,l) = 0.d0
        enddo
        a(k,k) = 1.d0
      enddo
      dd = a(1,1)*xno(1)+a(2,1)*xno(2)+a(3,1)*xno(3)
      if(dabs(dd).gt.0.999999999536d0) then
        ! project the z-axis
        dd = a(1,3)*xno(1)+a(2,3)*xno(2)+a(3,3)*xno(3)
         do l = 1,3
           e1(l) = a(l,3)-dd*xno(l)
         enddo
      else
        ! project the x-axis
        do l =1 ,3
          e1(l) = a(l,1)-dd*xno(l)
        enddo
      endif
      !
      dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
      do l = 1,3
        e1(l) = e1(l)/dd
      enddo
      ! e2 = n x e1
      e2(1) = xno(2)*e1(3) - e1(2)*xno(3)
      e2(2) = xno(3)*e1(1) - e1(3)*xno(1)
      e2(3) = xno(1)*e1(2) - e1(1)*xno(2)
      do l = 1,3
        e3(l) = xno(l)
      enddo
      ! put in tm and out
      do j = 1,3
         tm(1,j) = e1(j)
         tm(2,j) = e2(j)
         tm(3,j) = e3(j)
      enddo
      tmg(:,:) = 0.d0
      do i = 1,3
        do j = 1,3
          tmg(i,j)       = tm(i,j) ! f1
          tmg(i+3,j+3)   = tm(i,j) ! r1
          tmg(i+6,j+6)   = tm(i,j) ! f2
          tmg(i+9,j+9)   = tm(i,j) ! r2
          tmg(i+12,j+12) = tm(i,j) ! f3
          tmg(i+15,j+15) = tm(i,j) ! r3
        enddo
      enddo
      !      
      RETURN
      END
      !
      SUBROUTINE us3_Bs(X,Bs)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)  :: X(3,3)                 
      REAL*8, INTENT(OUT) :: Bs(2,18)   
      !
      REAL*8 :: x1,x2,x3,y1,y2,y3,x31,y31,y12,x21 
      REAL*8 :: Ae,x0,y0,X_1(3),Y_1(3),X_2(3),Y_2(3),Ai
      REAL*8 :: X_3(3),Y_3(3),bs1(2,6),bs2(2,6),bs3(2,6)
      REAL*8 :: B1(2,18),B2(2,18),B3(2,18),a3
      !
      Bs1(:,:) = 0.d0
      Bs2(:,:) = 0.d0
      Bs3(:,:) = 0.d0
      Bs(:,:)  = 0.d0
      !
      a3 = 1.d0/3.d0
      !
      x1 = X(1,1)
      x2 = X(2,1)
      x3 = X(3,1)
      y1 = X(1,2)
      y2 = X(2,2)
      y3 = X(3,2)
      x31 = x3-x1
      y31 = y3-y1
      y12 = y1-y2
      x21 = x2-x1
      !
      Ae  = 0.5d0*(x21*y31-x31*(-y12)) ! area
      !
      x0 = (x1+x2+x3)/3.d0
      y0 = (y1+y2+y3)/3.d0
      X_1 = (/X0,x1,x2/)
      Y_1 = (/Y0,y1,y2/) !  coords tri 1
      X_2 = (/X0,x2,x3/)
      Y_2 = (/Y0,y2,y3/) !  coords tri 2
      X_3 = (/X0,x3,x1/)
      Y_3 = (/Y0,y3,y1/) !  coords tri 3 
      !
      ! cell smoothing
      !
      call us3_CS(X_1,Y_1,bs1,bs2,bs3,Ai)
      B1(:,1:6)   = a3*bs1(:,:) + bs2(:,:)
      B1(:,7:12)  = a3*bs1(:,:) + bs3(:,:)
      B1(:,13:18) = a3*bs1(:,:)
      B1 = B1*Ai
      !
      call us3_CS(X_2,Y_2,bs1,bs2,bs3,Ai)
      B2(:,1:6)   = a3*bs1(:,:) 
      B2(:,7:12)  = a3*bs1(:,:) + bs2(:,:)
      B2(:,13:18) = a3*bs1(:,:) + bs3(:,:)
      B2 = B2*Ai      
      !
      call us3_CS(X_3,Y_3,bs1,bs2,bs3,Ai)
      B3(:,1:6)   = a3*bs1(:,:) + bs3(:,:)
      B3(:,7:12)  = a3*bs1(:,:) 
      B3(:,13:18) = a3*bs1(:,:) + bs2(:,:)
      B3 = B3*Ai       
      !
      Bs = (1.d0/Ae)*(B1(:,:)+B2(:,:)+B3(:,:))
      !
      RETURN
      END     
      !
      SUBROUTINE us3_Kp(X,Db,Ds,Kp)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)  :: X(3,3)         
      REAL*8, INTENT(IN)  :: Db(3,3)         
      REAL*8, INTENT(IN)  :: Ds(2,2)         
      REAL*8, INTENT(OUT) :: Kp(18,18)   
      !
      REAL*8 :: Bs(2,18),Bb(3,18),Ae
      !
       Ae  = 0.5d0*((X(2,1)-X(1,1))*(X(3,2)-X(1,2))
     & -(X(3,1)-X(1,1))*(-(X(1,2)-X(2,2)))) ! area
      !
      call us3_Bs(X,Bs)
      call us3_Bb(X(:,1),X(:,2),Bb)
      !
      Kp = (matmul(matmul(transpose(Bs),Ds),Bs) +
     & matmul(matmul(transpose(Bb),Db),Bb))*Ae
      !
      RETURN
      END
      !
      SUBROUTINE us3_CS(X,Y,bs1,bs2,bs3,Ae)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)  :: X(3),Y(3)               
      REAL*8, INTENT(OUT) :: bs1(2,6),bs2(2,6),bs3(2,6),Ae 
      !
      REAL*8 :: x13,x32,x21,y31,y12,y23 
      REAL*8 :: a1,a2,a3,a4
      !
      bs1(:,:) = 0.d0
      bs2(:,:) = 0.d0
      bs3(:,:) = 0.d0
      !
      x21 = X(2)-X(1)
      x13 = X(1)-X(3) 
      y31 = Y(3)-Y(1)
      y12 = Y(1)-Y(2)        
      x32 = X(3)-X(2)
      y23 = Y(2)-Y(3)
      ! triangle area - XY-tri
      Ae = 0.5d0*(x21*y31-x13*y12)
      a1 = 0.5d0*y12*x13
      a2 = 0.5d0*y31*x21  
      a3 = 0.5d0*x21*x13
      a4 = 0.5d0*y12*y31                               
      bs1(1,3) = +0.5d0*x32/Ae 
      bs1(1,4) = -0.5d0
      bs1(2,3) = +0.5d0*y23/Ae 
      bs1(2,5) = +0.5d0
      !
      bs2(1,3) = +0.5d0*x13/Ae 
      bs2(1,4) = +0.5d0*a1/Ae
      bs2(1,5) = +0.5d0*a3/Ae
      bs2(2,3) = +0.5d0*y31/Ae 
      bs2(2,4) = +0.5d0*a4/Ae
      bs2(2,5) = +0.5d0*a2/Ae
      !
      bs3(1,3) = +0.5d0*x21/Ae 
      bs3(1,4) = -0.5d0*a2/Ae
      bs3(1,5) = -0.5d0*a3/Ae
      bs3(2,3) = +0.5d0*y12/Ae 
      bs3(2,4) = -0.5d0*a4/Ae
      bs3(2,5) = -0.5d0*a1/Ae      
      !
      RETURN
      END
      !  
      SUBROUTINE us3_Bb(X,Y,bb)
      IMPLICIT NONE        
      REAL*8, INTENT(IN)  :: X(3),Y(3)               
      REAL*8, INTENT(OUT) :: bb(3,18)
      !
      REAL*8 :: x13,x31,x21,y31,y12
      REAL*8 :: y23,x32,Ae
      REAL*8 :: dNdx1,dNdx2,dNdx3
      REAL*8 :: dNdy1,dNdy2,dNdy3
      !
      bb(:,:) = 0.d0
      x21 = X(2)-X(1)
      x31 = X(3)-X(1)
      x32 = X(3)-X(2)
      y23 = Y(2)-Y(3)
      y31 = Y(3)-Y(1)
      y12 = Y(1)-Y(2)
      !
      Ae = 0.5d0*(x21*y31-x31*(-y12))
      !
      dNdx1 =  y23/(2.d0*Ae)
      dNdy1 =  x32/(2.d0*Ae)
      dNdx2 =  y31/(2.d0*Ae)
      dNdy2 = -x31/(2.d0*Ae)
      dNdx3 =  y12/(2.d0*Ae)
      dNdy3 =  x21/(2.d0*Ae)
      !
      bb(1,5)  = +dNdx1
      bb(1,11) = +dNdx2
      bb(1,17) = +dNdx3
      bb(2,4)  = -dNdy1
      bb(2,10) = -dNdy2
      bb(2,16) = -dNdy3
      bb(3,4)  = -dNdx1
      bb(3,5)  = +dNdy1
      bb(3,10) = -dNdx2
      bb(3,11) = +dNdy2
      bb(3,16) = -dNdx3
      bb(3,17) = +dNdy3
      !
      RETURN
      END
      ! 
      !  
      SUBROUTINE us3_Km(X,K,Qin,h)  
      IMPLICIT NONE  
      REAL*8, INTENT(IN)  :: X(3,3),Qin(3,3),h             
      REAL*8, INTENT(OUT) :: K(18,18)
      !
      REAL*8 :: alpha,ab,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9 
      REAL*8 :: x12,x23,x31,y12,y23,y31,A14
      REAL*8 :: x21,x32,x13,y21,y32,y13,T0(3,9),Te(3,3)
      REAL*8 :: Ae,A2,A4,LL21,LL32,LL13,h2,L(9,3),V 
      REAL*8 :: Q1(3,3),Q2(3,3),Q3(3,3),Q4(3,3),Q5(3,3),Q6(3,3)
      REAL*8 :: Enat(3,3),KO(3,3),Kb(9,9),Kh(9,9),Km(9,9)
      !
      INTEGER :: i
      !
      K(:,:)  = 0.d0
      Km(:,:) = 0.d0
      Kh(:,:) = 0.d0
      Kb(:,:) = 0.d0
      L(:,:)  = 0.d0
      T0(:,:) = 0.d0
      Te(:,:) = 0.d0
      Q1(:,:) = 0.d0
      Q2(:,:) = 0.d0
      Q3(:,:) = 0.d0 
      Q4(:,:) = 0.d0
      Q5(:,:) = 0.d0
      Q6(:,:) = 0.d0
      KO(:,:) = 0.d0
      !
      alpha = 1.d0/8.d0
      ab = alpha/6.d0
      b0 = (alpha**2)/4.d0 ! (4.0*(1.0+factor**2))
      b1 =  1.d0
      b2 =  2.d0
      b3 =  1.d0
      b4 =  0.d0
      b5 =  1.d0
      b6 = -1.d0
      b7 = -1.d0
      b8 = -1.d0
      b9 = -2.
      !
      x12 = X(1,1) - X(2,1)
      x23 = X(2,1) - X(3,1)
      x31 = X(3,1) - X(1,1)
      y12 = X(1,2) - X(2,2)
      y23 = X(2,2) - X(3,2)
      y31 = X(3,2) - X(1,2)
      x21 = -x12
      x32 = -x23
      x13 = -x31
      y21 = -y12
      y32 = -y23
      y13 = -y31
      !
      Ae = 0.5d0*(x21*y31-x31*(-y12))
      A2 = 2.d0*Ae
      A4 = 4.d0*Ae
      h2 = 0.5d0*h
      V  = Ae*h
      !
      LL21 = (x21**2) + (y21**2)
      LL32 = (x32**2) + (y32**2)
      LL13 = (x13**2) + (y13**2)
      !
      ! lumping matrix
      !
      L(1,1) = h2*y23
      L(1,3) = h2*x32
      L(2,2) = h2*x32
      L(2,3) = h2*y23
      L(3,1) = h2*y23*(y13-y21)*ab
      L(3,2) = h2*x32*(x31-x12)*ab
      L(3,3) = h2*(x31*y13-x12*y21)*2.d0*ab 
      !
      L(4,1) = h2*y31
      L(4,3) = h2*x13
      L(5,2) = h2*x13
      L(5,3) = h2*y31
      L(6,1) = h2*y31*(y21-y32)*ab
      L(6,2) = h2*x13*(x12-x23)*ab
      L(6,3) = h2*(x12*y21-x23*y32)*2.d0*ab
      !
      L(7,1) = h2*y12
      L(7,3) = h2*x21
      L(8,2) = h2*x21
      L(8,3) = h2*y12
      L(9,1) = h2*y12*(y32-y13)*ab
      L(9,2) = h2*x21*(x23-x31)*ab
      L(9,3) = h2*(x23*y32-x31*y13)*2.d0*ab 
      !
      ! basic stiffness
      !
      Kb = matmul(matmul(L,Qin),transpose(L))/V
      !
      ! trasformation hierachical rotations
      !
      T0(1,1) = x32/A4
      T0(1,2) = y32/A4      
      T0(1,3) = 1.d0
      T0(1,4) = x13/A4
      T0(1,5) = y13/A4      
      T0(1,7) = x21/A4
      T0(1,8) = y21/A4 
      !
      T0(2,1) = x32/A4
      T0(2,2) = y32/A4      
      T0(2,4) = x13/A4
      T0(2,5) = y13/A4
      T0(2,6) = 1.d0      
      T0(2,7) = x21/A4
      T0(2,8) = y21/A4 
      !
      T0(3,1) = x32/A4
      T0(3,2) = y32/A4      
      T0(3,4) = x13/A4
      T0(3,5) = y13/A4      
      T0(3,7) = x21/A4
      T0(3,8) = y21/A4
      T0(3,9) = 1.d0
      !
      !  transformation natural pattern
      !
      A14 = (1.d0/(Ae*A4))
      Te(1,1) = A14*y23*y13*LL21
      Te(1,2) = A14*y31*y21*LL32
      Te(1,3) = A14*y12*y32*LL13
      Te(2,1) = A14*x23*x13*LL21
      Te(2,2) = A14*x31*x21*LL32
      Te(2,3) = A14*x12*x32*LL13
      Te(3,1) = A14*(y23*x31+x32*y13)*LL21 
      Te(3,2) = A14*(y31*x12+x13*y21)*LL32 
      Te(3,3) = A14*(y12*x23+x21*y32)*LL13
      !
      A14 = (A2/3.d0)
      !
      !  nodal strain-displ.-matrix
      !
      Q1(1,1) = A14*b1/LL21
      Q1(1,2) = A14*b2/LL21
      Q1(1,3) = A14*b3/LL21
      Q1(2,1) = A14*b4/LL32
      Q1(2,2) = A14*b5/LL32
      Q1(2,3) = A14*b6/LL32
      Q1(3,1) = A14*b7/LL13
      Q1(3,2) = A14*b8/LL13
      Q1(3,3) = A14*b9/LL13
      !
      Q2(1,1) = A14*b9/LL21
      Q2(1,2) = A14*b7/LL21
      Q2(1,3) = A14*b8/LL21
      Q2(2,1) = A14*b3/LL32
      Q2(2,2) = A14*b1/LL32
      Q2(2,3) = A14*b2/LL32
      Q2(3,1) = A14*b6/LL13
      Q2(3,2) = A14*b4/LL13
      Q2(3,3) = A14*b5/LL13
      !
      Q3(1,1) = A14*b5/LL21
      Q3(1,2) = A14*b6/LL21
      Q3(1,3) = A14*b4/LL21
      Q3(2,1) = A14*b8/LL32
      Q3(2,2) = A14*b9/LL32
      Q3(2,3) = A14*b7/LL32
      Q3(3,1) = A14*b2/LL13
      Q3(3,2) = A14*b3/LL13
      Q3(3,3) = A14*b1/LL13      
      !
      Q4 = (Q1 + Q2)*0.5d0
      Q5 = (Q2 + Q3)*0.5d0
      Q6 = (Q3 + Q1)*0.5d0
      !
      Enat = matmul(matmul(transpose(Te),Qin),Te)
      !
      !  higher stiffness with respect to hier...rots
      !
      KO = (3.d0/4.d0)*b0*V*(matmul(matmul(transpose(Q4),Enat),Q4)
     & + matmul(matmul(transpose(Q5),Enat),Q5)
     & + matmul(matmul(transpose(Q6),Enat),Q6))
      !  higher stiffness [18x18)
      Kh = matmul(matmul(transpose(T0),KO),T0)
      Km = Kb + Kh
      !
      do i=1,3
        K(1+(i-1)*6,1:2)   = Km(1+(i-1)*3,1:2)   
        K(1+(i-1)*6,6)     = Km(1+(i-1)*3,3)
        K(1+(i-1)*6,7:8)   = Km(1+(i-1)*3,4:5)   
        K(1+(i-1)*6,12)    = Km(1+(i-1)*3,6)
        K(1+(i-1)*6,13:14) = Km(1+(i-1)*3,7:8)   
        K(1+(i-1)*6,18)    = Km(1+(i-1)*3,9)
        !zeile 2,8,14:
        K(2+(i-1)*6,1:2)   = Km(2+(i-1)*3,1:2)   
        K(2+(i-1)*6,6)     = Km(2+(i-1)*3,3)
        K(2+(i-1)*6,7:8)   = Km(2+(i-1)*3,4:5)   
        K(2+(i-1)*6,12)    = Km(2+(i-1)*3,6)
        K(2+(i-1)*6,13:14) = Km(2+(i-1)*3,7:8)   
        K(2+(i-1)*6,18)    = Km(2+(i-1)*3,9)
        !
        K(6+(i-1)*6,1:2)   = Km(3+(i-1)*3,1:2)   
        K(6+(i-1)*6,6)     = Km(3+(i-1)*3,3)
        K(6+(i-1)*6,7:8)   = Km(3+(i-1)*3,4:5)   
        K(6+(i-1)*6,12)    = Km(3+(i-1)*3,6)
        K(6+(i-1)*6,13:14) = Km(3+(i-1)*3,7:8)   
        K(6+(i-1)*6,18)    = Km(3+(i-1)*3,9)
      enddo
      !
      RETURN
      END
      !
      SUBROUTINE bm_ANDES(X,bm,Qin,h,r,s)
      IMPLICIT NONE  
      REAL*8, INTENT(IN)  :: X(3,3),Qin(3,3),h,r,s             
      REAL*8, INTENT(OUT) :: bm(3,18)
      !
      REAL*8 :: alpha,ab,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9 
      REAL*8 :: x12,x23,x31,y12,y23,y31,A14
      REAL*8 :: x21,x32,x13,y21,y32,y13,T0(3,9),Te(3,3)
      REAL*8 :: Ae,A2,A4,LL21,LL32,LL13,h2,L(9,3),V 
      REAL*8 :: Q1(3,3),Q2(3,3),Q3(3,3),Q4(3,3),Q5(3,3),Q6(3,3)
      REAL*8 :: Enat(3,3),matQ(3,3),bm2(3,9)
      !
      INTEGER :: i,i3,i6
      !
      L(:,:)  = 0.d0
      T0(:,:) = 0.d0
      Te(:,:) = 0.d0
      Q1(:,:) = 0.d0
      Q2(:,:) = 0.d0
      Q3(:,:) = 0.d0 
      Q4(:,:) = 0.d0
      Q5(:,:) = 0.d0
      Q6(:,:) = 0.d0
      !
      alpha = 1.d0/8.d0
      ab = alpha/6.d0
      b0 = (alpha**2)/4.d0 ! (4.0*(1.0+factor**2))
      b1 =  1.d0
      b2 =  2.d0
      b3 =  1.d0
      b4 =  0.d0
      b5 =  1.d0
      b6 = -1.d0
      b7 = -1.d0
      b8 = -1.d0
      b9 = -2.d0
      !
      x12 = X(1,1) - X(2,1)
      x23 = X(2,1) - X(3,1)
      x31 = X(3,1) - X(1,1)
      y12 = X(1,2) - X(2,2)
      y23 = X(2,2) - X(3,2)
      y31 = X(3,2) - X(1,2)
      x21 = -x12
      x32 = -x23
      x13 = -x31
      y21 = -y12
      y32 = -y23
      y13 = -y31
      !
      Ae = 0.5d0*(x21*y31-x31*(-y12))
      A2 = 2.d0*Ae
      A4 = 4.d0*Ae
      h2 = 0.5d0*h
      V  = Ae*h
      !
      LL21 = (x21**2) + (y21**2)
      LL32 = (x32**2) + (y32**2)
      LL13 = (x13**2) + (y13**2)
      !
      ! lumping matrix
      !
      L(1,1) = h2*y23
      L(1,3) = h2*x32
      L(2,2) = h2*x32
      L(2,3) = h2*y23
      L(3,1) = h2*y23*(y13-y21)*ab
      L(3,2) = h2*x32*(x31-x12)*ab
      L(3,3) = h2*(x31*y13-x12*y21)*2.d0*ab 
      !
      L(4,1) = h2*y31
      L(4,3) = h2*x13
      L(5,2) = h2*x13
      L(5,3) = h2*y31
      L(6,1) = h2*y31*(y21-y32)*ab
      L(6,2) = h2*x13*(x12-x23)*ab
      L(6,3) = h2*(x12*y21-x23*y32)*2.d0*ab
      !
      L(7,1) = h2*y12
      L(7,3) = h2*x21
      L(8,2) = h2*x21
      L(8,3) = h2*y12
      L(9,1) = h2*y12*(y32-y13)*ab
      L(9,2) = h2*x21*(x23-x31)*ab
      L(9,3) = h2*(x23*y32-x31*y13)*2.d0*ab 
      !
      !
      ! trasformation hierachical rotations
      !
      T0(1,1) = x32/A4
      T0(1,2) = y32/A4      
      T0(1,3) = 1.d0
      T0(1,4) = x13/A4
      T0(1,5) = y13/A4      
      T0(1,7) = x21/A4
      T0(1,8) = y21/A4 
      !
      T0(2,1) = x32/A4
      T0(2,2) = y32/A4      
      T0(2,4) = x13/A4
      T0(2,5) = y13/A4
      T0(2,6) = 1.d0      
      T0(2,7) = x21/A4
      T0(2,8) = y21/A4 
      !
      T0(3,1) = x32/A4
      T0(3,2) = y32/A4      
      T0(3,4) = x13/A4
      T0(3,5) = y13/A4      
      T0(3,7) = x21/A4
      T0(3,8) = y21/A4
      T0(3,9) = 1.d0
      !
      !  transformation natural pattern
      !
      A14 = (1.d0/(Ae*A4))
      Te(1,1) = A14*y23*y13*LL21
      Te(1,2) = A14*y31*y21*LL32
      Te(1,3) = A14*y12*y32*LL13
      Te(2,1) = A14*x23*x13*LL21
      Te(2,2) = A14*x31*x21*LL32
      Te(2,3) = A14*x12*x32*LL13
      Te(3,1) = A14*(y23*x31+x32*y13)*LL21 
      Te(3,2) = A14*(y31*x12+x13*y21)*LL32 
      Te(3,3) = A14*(y12*x23+x21*y32)*LL13
      !
      A14 = (A2/3.d0)
      !
      !  nodal strain-displ.-matrix
      !
                   
      Q1(1,1) = A14*b1/LL21
      Q1(1,2) = A14*b2/LL21
      Q1(1,3) = A14*b3/LL21
      Q1(2,1) = A14*b4/LL32
      Q1(2,2) = A14*b5/LL32
      Q1(2,3) = A14*b6/LL32
      Q1(3,1) = A14*b7/LL13
      Q1(3,2) = A14*b8/LL13
      Q1(3,3) = A14*b9/LL13
      !
      Q2(1,1) = A14*b9/LL21
      Q2(1,2) = A14*b7/LL21
      Q2(1,3) = A14*b8/LL21
      Q2(2,1) = A14*b3/LL32
      Q2(2,2) = A14*b1/LL32
      Q2(2,3) = A14*b2/LL32
      Q2(3,1) = A14*b6/LL13
      Q2(3,2) = A14*b4/LL13
      Q2(3,3) = A14*b5/LL13
      !
      Q3(1,1) = A14*b5/LL21
      Q3(1,2) = A14*b6/LL21
      Q3(1,3) = A14*b4/LL21
      Q3(2,1) = A14*b8/LL32
      Q3(2,2) = A14*b9/LL32
      Q3(2,3) = A14*b7/LL32
      Q3(3,1) = A14*b2/LL13
      Q3(3,2) = A14*b3/LL13
      Q3(3,3) = A14*b1/LL13      
      !
      Q4 = (Q1 + Q2)*0.5d0
      Q5 = (Q2 + Q3)*0.5d0
      Q6 = (Q3 + Q1)*0.5d0
      !
      Enat = matmul(matmul(transpose(Te),Qin),Te)     
      !     
      matQ = (1.d0-r-s)*Q1 + r*Q2 + s*Q3
      bm2  = (3.d0/2.d0)*(b0**0.5)*matmul(matmul(Te,matQ),T0) 
      bm2 = bm2 + transpose(L)/V
      !
      bm(:,:) = 0.d0
      !
      do i=1,3
        !
        i3 = (i-1)*3
        i6 = (i-1)*6
        !
        bm(1:3,1+i6:2+i6) = bm2(1:3,1+i3:2+i3)
        bm(1:3,6+i6) = bm2(1:3,3+i3) 
        !
      enddo
      !
      RETURN
      END
      !
      SUBROUTINE us3_matma(e,un,h,Dm,Db,Ds,Qin,Qs,x,iflag)  
      IMPLICIT NONE  
      REAL*8, INTENT(IN)  :: e,un,h,x(3,3) 
      INTEGER, INTENT(IN) :: iflag
      REAL*8, INTENT(OUT) :: Dm(3,3),Db(3,3),Ds(2,2),Qin(3,3),Qs(2,2)
      REAL*8 :: q1,kap,hei(3),he
      INTEGER :: k
      !
      ! Reduced material matrix (S33=0)+ integration in thichkness dirc.
      ! pre- in thickness dirc. integrated:
      if(iflag.eq.1) then
      Qin(:,:) = 0.d0
      q1 = e/(1.d0-un**2)
      Qin(1,1) = q1 
      Qin(1,2) = q1*un 
      Qin(2,1) = q1*un
      Qin(2,2) = q1 
      Qin(3,3) = q1*(1.d0-un)/2.d0
      ! membrane:
      Dm = Qin*h 
      ! bending:
      Db = Qin*h**3/12.d0
      ! shear
      
            !
      ! shear correction
      !
      ! l1 = 1->2
!       hei(1) = dsqrt((x(2,1)-x(1,1))**2+(x(2,2)-x(1,2))**2)
!       ! l1 = 2->3
!       hei(2) = dsqrt((x(3,1)-x(2,1))**2+(x(3,2)-x(2,2))**2)
!       ! l1 = 3->1
!       hei(3) = dsqrt((x(3,1)-x(1,1))**2+(x(3,2)-x(1,2))**2)
!       !
!       he = 0.d0      
!       do k = 1,3
!         if(he.LT.abs(hei(k))) then
!           he = abs(hei(k))
!         endif
!       enddo
      !
!       kap = 5.d0/6.d0
!       q1 = (h**3*kap/(h**2+0.1d0*he**2)*e/2.0/(1+un))
!       Ds(:,:)  = 0.0
!       Ds(1,1)  = q1
!       Ds(2,2)  = q1  
!       Qs(:,:) = 0.d0
!       q1 = (h**2*kap/(h**2+0.1d0*he**2)*e/2.0/(1+un))
!       Qs(1,1)  = q1
!       Qs(2,2)  = q1 

      Qs(:,:) = 0.d0
      kap = 5.d0/6.d0
      q1 = e/(2.d0*(1.d0+un))
      Qs(1,1)  = q1*kap
      Qs(2,2)  = q1*kap  
      Ds = Qs*h 
      else
      Qin(:,:) = 0.d0
      q1 = e/(1.d0-un**2)
      Qin(1,1) = q1 
      Qin(1,2) = q1*un 
      Qin(2,1) = q1*un
      Qin(2,2) = q1 
      Qin(3,3) = q1*(1.d0-un)/2.d0
      Qs(:,:) = 0.d0
      kap = 5.d0/6.d0
      q1 = e/(2.d0*(1.d0+un))
      Qs(1,1)  = q1*kap
      Qs(2,2)  = q1*kap  
      endif       
      RETURN
      END
      !
      SUBROUTINE us3_M(X,h,rho,M)  
      IMPLICIT NONE  
      REAL*8, INTENT(IN)  :: X(3,3),rho,h             
      REAL*8, INTENT(OUT) :: M(18,18) 
      REAL*8 :: MLST(12,12),Ae,Ai,TCH(12,9),alpha
      REAL*8 :: Mm(9,9),points3(3,2),w3,r,s,Nrs(3)
      REAL*8 :: N_u(6,18),q1,m_3t(6,6)
      REAL*8 :: x12,x23,x31,y12,y23,y31
      REAL*8 :: x21,x32,x13,y21,y32,y13
      INTEGER :: i,j
      !
      alpha = 1.d0/8.d0
      !
      points3(1,1) = 0.1666666666667
      points3(1,2) = 0.1666666666667
      points3(2,1) = 0.6666666666667
      points3(2,2) = 0.1666666666667
      points3(3,1) = 0.1666666666667
      points3(3,2) = 0.6666666666667
      w3 = 0.333333333333333
      !
      x12 = X(1,1) - X(2,1)
      x23 = X(2,1) - X(3,1)
      x31 = X(3,1) - X(1,1)
      y12 = X(1,2) - X(2,2)
      y23 = X(2,2) - X(3,2)
      y31 = X(3,2) - X(1,2)
      x21 = -x12
      x32 = -x23
      x13 = -x31
      y21 = -y12
      y32 = -y23
      y13 = -y31
!       !
      Ae = 0.5d0*(x21*y31-x31*(-y12))      
!       Ai = Ae*h*rho/180.d0
!       !
!       MLST(:,:) = 0.d0
!       MLST(1,1)=+Ai*6.d0
!       MLST(1,3)=-Ai*1.d0
!       MLST(1,5)=-Ai*1.d0
!       MLST(1,9)=-Ai*4.d0
!       MLST(2,2)=+Ai*6.d0
!       MLST(2,4)=-Ai*1.d0
!       MLST(2,6)=-Ai*1.d0
!       MLST(2,10)=-Ai*4.d0
!       MLST(3,5)=-Ai*4.d0
!       MLST(3,7)=+Ai*32.d0
!       MLST(3,9)=+Ai*16.d0
!       MLST(3,11)=+Ai*16.d0
!       MLST(4,6)=-Ai*4.d0
!       MLST(4,8)=+Ai*32.d0
!       MLST(4,10)=+Ai*16.d0
!       MLST(4,12)=+Ai*16.d0
!       MLST(5,1)=-Ai*1.d0
!       MLST(5,3) =+Ai*6.d0
!       MLST(5,5)=-Ai*1.d0
!       MLST(5,11)=-Ai*4.d0
!       MLST(6,2)=-Ai*1.d0
!       MLST(6,4) =+Ai*6.d0
!       MLST(6,6)=-Ai*1.d0
!       MLST(6,12)=-Ai*4.d0
!       MLST(7,1)=-Ai*4.d0
!       MLST(7,7) =+Ai*16.d0
!       MLST(7,9)=+Ai*32.d0
!       MLST(7,11)=+Ai*16.d0
!       MLST(8,2)=-Ai*4.d0
!       MLST(8,8) =+Ai*16.d0
!       MLST(8,10)=+Ai*32.d0
!       MLST(8,12)=+Ai*16.d0
!       MLST(9,1)=-Ai*1.d0
!       MLST(9,3) =-Ai*1.d0
!       MLST(9,5)=+Ai*6.d0
!       MLST(9,7) =-Ai*4.d0
!       MLST(10,2)=-Ai*1.d0
!       MLST(10,4) =-Ai*1.d0
!       MLST(10,6)=+Ai*6.d0
!       MLST(10,8) =-Ai*4.d0                                                        
!       MLST(11,3)=-Ai*4.d0
!       MLST(11,7)  =+Ai*16.d0
!       MLST(11,9)=+Ai*16.d0
!       MLST(11,11) =+Ai*32.d0
!       MLST(12,4)=-Ai*4.d0
!       MLST(12,8)  =+Ai*16.d0
!       MLST(12,10)=+Ai*16.d0
!       MLST(12,12) =+Ai*32.d0 
!       !
!       TCH(:,:) = 0.d0
!       TCH(1,1) = 1.d0
!       TCH(2,2) = 1.d0
!       TCH(3,1) = 0.5d0
!       TCH(3,3) = 0.125d0*alpha*y12 
!       TCH(3,6) = 0.125d0*alpha*y21  
!       TCH(3,7) = 0.5d0
!       TCH(4,2) = 0.5d0
!       TCH(4,3) = 0.125d0*alpha*x21 
!       TCH(4,6) = 0.125d0*alpha*x12  
!       TCH(4,8) = 0.5d0     
!       TCH(5,4) = 1.d0
!       TCH(6,5) = 1.d0 
!       TCH(7,4) = 0.5d0
!       TCH(7,6) = 0.125d0*alpha*y23
!       TCH(7,9) = 0.125d0*alpha*y32
!       TCH(7,7) = 0.5d0
!       TCH(8,5) = 0.5d0
!       TCH(8,6) = 0.125d0*alpha*x32
!       TCH(8,9) = 0.125d0*alpha*x23
!       TCH(8,8) = 0.5d0
!       TCH(9,7)  = 1.d0
!       TCH(10,8) = 1.d0
!       TCH(11,1) = 0.5d0
!       TCH(11,3) = 0.125d0*alpha*y13
!       TCH(11,9) = 0.125d0*alpha*y31
!       TCH(11,7) = 0.5d0
!       TCH(12,2) = 0.5d0
!       TCH(12,3) = 0.125d0*alpha*x31
!       TCH(12,9) = 0.125d0*alpha*x13
!       TCH(12,8) = 0.5d0
!       !
!       Mm = matmul(matmul(transpose(TCH),MLST),TCH)
!       !
!       do i=1,3
!         M(1+(i-1)*6,1:2)   = Mm(1+(i-1)*3,1:2)   
!         M(1+(i-1)*6,6)     = Mm(1+(i-1)*3,3)
!         M(1+(i-1)*6,7:8)   = Mm(1+(i-1)*3,4:5)   
!         M(1+(i-1)*6,12)    = Mm(1+(i-1)*3,6)
!         M(1+(i-1)*6,13:14) = Mm(1+(i-1)*3,7:8)   
!         M(1+(i-1)*6,18)    = Mm(1+(i-1)*3,9)
!         !zeile 2,8,14:
!         M(2+(i-1)*6,1:2)   = Mm(2+(i-1)*3,1:2)   
!         M(2+(i-1)*6,6)     = Mm(2+(i-1)*3,3)
!         M(2+(i-1)*6,7:8)   = Mm(2+(i-1)*3,4:5)   
!         M(2+(i-1)*6,12)    = Mm(2+(i-1)*3,6)
!         M(2+(i-1)*6,13:14) = Mm(2+(i-1)*3,7:8)   
!         M(2+(i-1)*6,18)    = Mm(2+(i-1)*3,9)
!         !
!         M(6+(i-1)*6,1:2)   = Mm(3+(i-1)*3,1:2)   
!         M(6+(i-1)*6,6)     = Mm(3+(i-1)*3,3)
!         M(6+(i-1)*6,7:8)   = Mm(3+(i-1)*3,4:5)   
!         M(6+(i-1)*6,12)    = Mm(3+(i-1)*3,6)
!         M(6+(i-1)*6,13:14) = Mm(3+(i-1)*3,7:8)   
!         M(6+(i-1)*6,18)    = Mm(3+(i-1)*3,9)
!       enddo
!       !
      M(:,:) = 0.d0
      !
      N_u(:,:)     = 0.d0 
      m_3t(:,:)    = 0.d0
      q1 = rho*h
      m_3t(1,1) = q1
      m_3t(2,2) = q1
      m_3t(3,3) = q1
      q1 = (rho*h**3)/12.d0
      m_3t(4,4) = q1 
      m_3t(5,5) = q1 
      !
      do i=1,3
        r = points3(i,1)
        s = points3(i,2)
        Nrs(1) = 1.d0-r-s
        Nrs(2) = r
        Nrs(3) = s
        !
        do j = 1,3
          N_u(1,1+(j-1)*6) = Nrs(j)
          N_u(2,2+(j-1)*6) = Nrs(j)
          N_u(3,3+(j-1)*6) = Nrs(j)
          N_u(4,4+(j-1)*6) = Nrs(j)
          N_u(5,5+(j-1)*6) = Nrs(j)
        enddo
        !
        M = M + matmul(matmul(transpose(N_u),m_3t),N_u)*Ae*w3
      enddo
      !
      RETURN
      END
      !       
      !  
      SUBROUTINE us3_ae(X,Ae)  
      IMPLICIT NONE  
      REAL*8, INTENT(IN)  :: X(3,3)            
      REAL*8, INTENT(OUT) :: Ae
      !
      Ae  = 0.5d0*((X(2,1)-X(1,1))*(X(3,2)-X(1,2))
     & -(X(3,1)-X(1,1))*(-(X(1,2)-X(2,2)))) ! area
      ! 
      RETURN
      END
      !  
      SUBROUTINE us3_xu(x,ushell,ueg,xg,tm,vl)   
      IMPLICIT NONE  
      REAL*8, INTENT(IN)  :: xg(3,3),tm(3,3),vl(6,3)            
      REAL*8, INTENT(OUT) :: x(3,3),ushell(18),ueg(18)
      !
      x(:,:) = 0.d0
      x(1,:) = matmul(tm,xg(1,:))
      x(2,:) = matmul(tm,xg(2,:))
      x(3,:) = matmul(tm,xg(3,:))
      !
      ushell(1:3)   = matmul(tm,vl(1:3,1))
      ushell(4:6)   = matmul(tm,vl(4:6,1))
      ushell(7:9)   = matmul(tm,vl(1:3,2))
      ushell(10:12) = matmul(tm,vl(4:6,2))
      ushell(13:15) = matmul(tm,vl(1:3,3))
      ushell(16:18) = matmul(tm,vl(4:6,3))
      !
      ueg(:) = 0.d0
      ueg(1:3)   = vl(1:3,1)
      ueg(4:6)   = vl(4:6,1)
      ueg(7:9)   = vl(1:3,2)
      ueg(10:12) = vl(4:6,2)
      ueg(13:15) = vl(1:3,3)
      ueg(16:18) = vl(4:6,3)
      !
      RETURN
      END
      !
      SUBROUTINE us3_linel_Qi(e,un,Qin,Qs) 
      IMPLICIT NONE  
      REAL*8, INTENT(IN)  :: e,un
      REAL*8, INTENT(OUT) :: Qin(3,3),Qs(2,2)
      REAL*8 :: q1,kap
      INTEGER :: k
      !
      ! Reduced material matrix (S33=0)
      !
      Qin(:,:) = 0.d0
      q1 = e/(1.d0-un**2)
      Qin(1,1) = q1 
      Qin(1,2) = q1*un 
      Qin(2,1) = q1*un
      Qin(2,2) = q1 
      Qin(3,3) = q1*(1.d0-un)/2.d0
      !
      Qs(:,:) = 0.d0
      kap = 5.d0/6.d0
      q1 = e/(2.d0*(1.d0+un))
      Qs(1,1)  = q1*kap
      Qs(2,2)  = q1*kap  
      RETURN
      END
      !      
