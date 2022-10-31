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
      subroutine materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &  imat,amat,iorien,pgauss,orab,ntmat_,elas,rho,iel,ithermal,
     &  alzero,mattyp,t0l,t1l,ihyper,istiff,elconloc,eth,kode,plicon,
     &  nplicon,plkcon,nplkcon,npmat_,plconloc,mi,dtime,iint,
     &  xstiff,ncmat_)
!
      implicit none
!
!     determines the material data for element iel
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
            call ident2(alcon(0,1,imat),t1l,nalcon(2,imat),seven,id)
            if(nalcon(2,imat).eq.0) then
               do k=1,6
                  alph(k)=0.d0
               enddo
               continue
            elseif(nalcon(2,imat).eq.1) then
               do k=1,nalcon(1,imat)
                  alph(k)=alcon(k,1,imat)*(t1l-alzero(imat))
               enddo
            elseif(id.eq.0) then
               do k=1,nalcon(1,imat)
                  alph(k)=alcon(k,1,imat)*(t1l-alzero(imat))
               enddo
            elseif(id.eq.nalcon(2,imat)) then
               do k=1,nalcon(1,imat)
                  alph(k)=alcon(k,id,imat)*(t1l-alzero(imat))
               enddo
            else
               do k=1,nalcon(1,imat)
                  alph(k)=(alcon(k,id,imat)+
     &                 (alcon(k,id+1,imat)-alcon(k,id,imat))*
     &                 (t1l-alcon(0,id,imat))/
     &                 (alcon(0,id+1,imat)-alcon(0,id,imat)))
     &                 *(t1l-alzero(imat))
               enddo
            endif
!     
!     subtracting the initial temperature influence       
!     
            call ident2(alcon(0,1,imat),t0l,nalcon(2,imat),seven,id)
            if(nalcon(2,imat).eq.0) then
               continue
            elseif(nalcon(2,imat).eq.1) then
               do k=1,nalcon(1,imat)
                  alph(k)=alph(k)-alcon(k,1,imat)*(t0l-alzero(imat))
               enddo
            elseif(id.eq.0) then
               do k=1,nalcon(1,imat)
                  alph(k)=alph(k)-alcon(k,1,imat)*(t0l-alzero(imat))
               enddo
            elseif(id.eq.nalcon(2,imat)) then
               do k=1,nalcon(1,imat)
                  alph(k)=alph(k)-alcon(k,id,imat)*(t0l-alzero(imat))
               enddo
            else
               do k=1,nalcon(1,imat)
                  alph(k)=alph(k)-(alcon(k,id,imat)+
     &                 (alcon(k,id+1,imat)-alcon(k,id,imat))*
     &                 (t0l-alcon(0,id,imat))/
     &                 (alcon(0,id+1,imat)-alcon(0,id,imat)))
     &                 *(t0l-alzero(imat))
               enddo
            endif
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
