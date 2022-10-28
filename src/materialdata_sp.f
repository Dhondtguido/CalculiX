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
      subroutine materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &  elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!
      implicit none
!
!     determines the material data for element i
!
      integer nelcon(2,*),imat,i,k,kin,ntmat_,nelconst,kode,
     &  itemp,ncmat_,id,nplicon(0:ntmat_,*),npmat_
!
      real*8 elcon(0:ncmat_,ntmat_,*),t1l,elconloc(*),
     &   plicon(0:2*npmat_,ntmat_,*),plconloc(802)
!
!
!
!     nelconst: # constants read from file
!     
!     calculating the hardening coefficients
!     
!     for the calculation of the spring stiffness, the whole curve
!     has to be stored:
!     plconloc(2*k-1), k=1...200: displacement
!     plconloc(2*k),k=1...200:    force
!     
      if(kode.lt.-50) then
         if(npmat_.eq.0) then
            plconloc(801)=0.5d0
            plconloc(802)=0.5d0
         else
            plconloc(1)=0.d0
            plconloc(2)=0.d0
            plconloc(3)=0.d0
            plconloc(801)=nplicon(1,imat)+0.5d0
            plconloc(802)=0.5d0
!     
!     nonlinear spring characteristic or gap conductance characteristic
!     
            if(nplicon(1,imat).ne.0) then
!     
               if(nplicon(0,imat).eq.1) then
                  id=-1
               else
                  call ident2(plicon(0,1,imat),t1l,nplicon(0,imat),
     &                 2*npmat_+1,id)
               endif
!     
               if(nplicon(0,imat).eq.0) then
                  continue
               elseif((nplicon(0,imat).eq.1).or.(id.eq.0).or.
     &                 (id.eq.nplicon(0,imat))) then
                  if(id.le.0) then
                     itemp=1
                  else
                     itemp=id
                  endif
                  kin=0
                  call plcopy(plicon,nplicon,plconloc,npmat_,ntmat_,
     &                 imat,itemp,i,kin)
                  if((id.eq.0).or.(id.eq.nplicon(0,imat))) then
                  endif
               else
                  kin=0
                  call plmix(plicon,nplicon,plconloc,npmat_,ntmat_,
     &                 imat,id+1,t1l,i,kin)
               endif
            endif
         endif
      else
!     
!     linear spring characteristic
!     
         nelconst=nelcon(1,imat)
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
     &              (elcon(k,id+1,imat)-elcon(k,id,imat))*
     &              (t1l-elcon(0,id,imat))/
     &              (elcon(0,id+1,imat)-elcon(0,id,imat))
            enddo
         endif
      endif
!     
      return
      end
      
