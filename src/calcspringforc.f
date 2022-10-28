!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine calcspringforc(imat,elcon,nelcon,ncmat_,ntmat_,t1l,
     &  kode,plicon,nplicon,npmat_,senergy,nener,fk,val)
!
!     calculates the spring forc and the spring energy (node-to-face penalty)
!
      implicit none
!
      integer i,imat,ncmat_,ntmat_,kode,niso,id,nplicon(0:ntmat_,*),
     &    npmat_,nelcon(2,*),nener
!
      real*8 t1l,elcon(0:ncmat_,ntmat_,*),elconloc(ncmat_),
     &  xk,fk,val,xiso(200),yiso(200),plicon(0:2*npmat_,ntmat_,*),
     &  senergy,plconloc(802)
!
!     interpolating the material data
!     
      call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &     elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!     
!     calculating the spring force and the spring constant
!     
      if(kode.eq.2)then
         xk=elconloc(1)
         fk=xk*val
         if(nener.eq.1) then
            senergy=fk*val/2.d0
         endif
      else
         niso=int(plconloc(801))
         do i=1,niso
            xiso(i)=plconloc(2*i-1)
            yiso(i)=plconloc(2*i)
         enddo
         call ident(xiso,val,niso,id)
         if(id.eq.0) then
            xk=0.d0
            fk=yiso(1)
            if(nener.eq.1) then
               senergy=fk*val
            endif
         elseif(id.eq.niso) then
            xk=0.d0
            fk=yiso(niso)
            if(nener.eq.1) then
               senergy=yiso(1)*xiso(1)
               do i=2,niso
                  senergy=senergy+(xiso(i)-xiso(i-1))*
     &                 (yiso(i)+yiso(i-1))/2.d0
               enddo
               senergy=senergy+(val-xiso(niso))*yiso(niso)
            endif
         else
            xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
            fk=yiso(id)+xk*(val-xiso(id))
            if(nener.eq.1) then
               senergy=yiso(1)*xiso(1)
               do i=2, id
                  senergy=senergy+(xiso(i)-xiso(i-1))*
     &                 (yiso(i)+yiso(i-1))/2.d0
               enddo
               senergy=senergy+(val-xiso(id))*(fk+yiso(id))/2.d0
            endif
         endif
      endif
!     
      return
      end

