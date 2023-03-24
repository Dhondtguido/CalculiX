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
!     function updating transformed u and inserting inactive dofs due to transformed mpcs,
!     used in stressmortar
!      
!     Author: Saskia Sitzmann
!
      subroutine resultsini_mortar(nk,v,ithermal,iperturb,
     &  nactdof,iout,vold,b,nodeboun,ndirboun,
     &  xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,
     &  bet,gam,dtime,mi)
!
!     initialization 
!
!     1. storing the calculated primary variables nodewise
!     2. inserting the boundary conditions nodewise (SPC's and MPC's)
!     3. determining which derived variables (strains, stresses,
!        internal forces...) have to be calculated
!
      implicit none
!
      character*20 labmpc(*)
!
      integer mi(*),nactdof(0:mi(2),*),nodeboun(*),ndirboun(*),
     &  ipompc(*),nodempc(3,*),mt,nk,ithermal(*),i,j,
     &  iperturb(*),iout,nboun,nmpc,nmethod,ist,ndir,node,index,
     &  incrementalmpc
!
      real*8 v(0:mi(2),*),cam(5),vold(0:mi(2),*),b(*),xboun(*),
     &  coefmpc(*),bet,gam,dtime,scal1,scal2,bnac,fixed_disp
!
      mt=mi(2)+1
!
      if((iout.ne.2).and.(iout.gt.-1)) then
!     
         if((nmethod.ne.4).or.(iperturb(1).le.1)) then
            if(ithermal(1).ne.2) then
               do i=1,nk
                  do j=1,mi(2)
                     if(nactdof(j,i).gt.0) then
                        bnac=b(nactdof(j,i))
                     else
                        cycle
                     endif
                     v(j,i)=v(j,i)+bnac
                     if((iperturb(1).ne.0).and.(abs(nmethod).eq.1)) then
                        if(dabs(bnac).gt.cam(1)) then
                           cam(1)=dabs(bnac)
                           cam(4)=nactdof(j,i)-0.5d0
                        endif
                     endif
                  enddo
               enddo
            endif
            if(ithermal(1).gt.1) then
               do i=1,nk
                  if(nactdof(0,i).gt.0) then
                     bnac=b(nactdof(0,i))
                  else
                     cycle
                  endif
                  v(0,i)=v(0,i)+bnac
                  if((iperturb(1).ne.0).and.(abs(nmethod).eq.1)) then
                     if(dabs(bnac).gt.cam(2)) then
                        cam(2)=dabs(bnac)
                        cam(5)=nactdof(0,i)-0.5d0
                     endif
                  endif
               enddo
            endif
!     
         else
!     
!     direct integration dynamic step
!     b contains the acceleration increment
!     
            if(ithermal(1).ne.2) then
               scal1=bet*dtime*dtime
               scal2=gam*dtime
               do i=1,nk
                  do j=1,mi(2)
                     if(nactdof(j,i).gt.0) then
                        bnac=b(nactdof(j,i))
                     else
                        cycle
                     endif
                     v(j,i)=v(j,i)+scal1*bnac
                  enddo
               enddo
            endif
            if(ithermal(1).gt.1) then
               do i=1,nk
                  if(nactdof(0,i).gt.0) then
                     bnac=b(nactdof(0,i))
                  else
                     cycle
                  endif
                  v(0,i)=v(0,i)+bnac
               enddo
            endif
         endif
!     
      endif
!     
!     SPC's and MPC's have to be taken into account for 
!     iout=0,1 and -1
!     
      if(abs(iout).lt.2) then
!     
!     inserting the boundary conditions
!     
         do i=1,nboun
            if(ndirboun(i).gt.mi(2)) cycle
            fixed_disp=xboun(i)
            v(ndirboun(i),nodeboun(i))=fixed_disp
     &           -vold(ndirboun(i),nodeboun(i))
         enddo
!     
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!     
         do i=1,nmpc
            if((labmpc(i)(1:20).eq.'                    ').or.
     &           (labmpc(i)(1:7).eq.'CONTACT').or.
     &           (labmpc(i)(1:6).eq.'CYCLIC').or.
     &           (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
               incrementalmpc=0
            else
               if((nmethod.eq.2).or.(nmethod.eq.3).or.
     &              ((iperturb(1).eq.0).and.(abs(nmethod).eq.1)))
     &              then
                  incrementalmpc=0
               else
                  incrementalmpc=1
               endif
            endif
            ist=ipompc(i)
            node=nodempc(1,ist)
            ndir=nodempc(2,ist)
            if(ndir.eq.0) then
               if(ithermal(1).lt.2) cycle
            elseif(ndir.gt.mi(2)) then
               cycle
            else
               if(ithermal(1).eq.2) cycle
            endif
            index=nodempc(3,ist)
            fixed_disp=0.d0
            if(index.ne.0) then
               do
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    v(nodempc(2,index),nodempc(1,index))
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
            fixed_disp=fixed_disp/coefmpc(ist)
            v(ndir,node)=fixed_disp
         enddo
      endif
!     
      return
      end
