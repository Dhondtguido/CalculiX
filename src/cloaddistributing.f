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
      subroutine cloaddistributing(noderef,iref,val,nodeforc,ndirforc,
     &     xforc,nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,
     &     co,ikforc,ilforc,isector,add,user,idefforc,ipompc,nodempc,
     &     nmpc,ikmpc,ilmpc,labmpc,edc,id,orab,coeffc,ier)
!     
!     distributed a distributing load among the nodes belonging to
!     the distributing surface
!     
      implicit none
!     
      logical add,user
!     
      character*20 labmpc(*)
!     
      integer nodeforc(2,*),ndirforc(*),noderef,iref,nforc,nforc_,j,
     &     iamforc(*),iamplitude,nam,ntrans,inotr(2,*),irefm3,ier,
     &     ikforc(*),ilforc(*),idof,id,isector,idefforc(*),ipompc(*),
     &     nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),istart,
     &     iend,iorien,node,idir
!     
      real*8 xforc(*),val,trab(7,*),co(3,*),orab(7,*),edc(12,*),
     &     coeffc(0:6,*),forcval,skl(3,3),xfor(6),e1(3),e2(3),e3(3)
!
!     applicable distributing coefficients in field coeffc
!
      istart=int(edc(1,id))
      iend=int(edc(2,id))
!
!     applicable orientation
!
      iorien=int(edc(3,id))
!
!     local coordinate system
!
      e1(1)=edc(4,id)
      e1(2)=edc(5,id)
      e1(3)=edc(6,id)
      e2(1)=edc(7,id)
      e2(2)=edc(8,id)
      e2(3)=edc(9,id)
      e3(1)=edc(10,id)
      e3(2)=edc(11,id)
      e3(3)=edc(12,id)
!
!     check for a transformation in the coupling node (= noderef)
!
      if(ntrans.gt.0) then
        if(inotr(1,noderef).ne.0) then
          write(*,*) '*ERROR in cloaddistributing: a transformation'
          write(*,*) '       is not allowed in the reference node'
          write(*,*) '       of a *distributing definition;'
          write(*,*) '       reference node: ',noderef
          write(*,*)
          ier=1
          return
        endif
      endif
!
      if(iorien.eq.0) then
!
!       each line in coeffc leads to a concentrated force
!
        do j=istart,iend
          idof=int(coeffc(0,j))
          node=int(idof/10.d0)
          if(ntrans.gt.0) then
            if(inotr(1,node).gt.0) then
              write(*,*)
     &'*ERROR in cloaddistributing: no transformation'
              write(*,*) '       is allowed in a node belonging to a'
              write(*,*) '       distributing coupling surface;'
              write(*,*) '       node at fault:',node
              write(*,*)
              ier=1
              cycle
            endif
          endif
          idir=int(idof)-10*node
          if(iref.le.3) then
!
!           force: local e1-e2-e3 surface coordinate system applies
!
            forcval=val*(e1(iref)*coeffc(1,j)+
     &                   e2(iref)*coeffc(2,j)+
     &                   e3(iref)*coeffc(3,j))
!          
            call forcadd(node,idir,forcval,nodeforc,ndirforc,xforc,
     &           nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,
     &           co,ikforc,ilforc,isector,add,user,idefforc,ipompc,
     &           nodempc,nmpc,ikmpc,ilmpc,labmpc)
          else
!
!           moment: local e1-e2-e3 surface coordinate system applies
!
            irefm3=iref-3
            forcval=val*(e1(irefm3)*coeffc(4,j)+
     &                   e2(irefm3)*coeffc(5,j)+
     &                   e3(irefm3)*coeffc(6,j))
!          
            call forcadd(node,idir,forcval,nodeforc,ndirforc,xforc,
     &           nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,
     &           co,ikforc,ilforc,isector,add,user,idefforc,ipompc,
     &           nodempc,nmpc,ikmpc,ilmpc,labmpc)
          endif
        enddo
      else
!
!       orientation is active
!
        call transformatrix(orab(1,iorien),co(1,noderef),skl)
!
        if(iref.le.3) then
!
!         force applied in the reference node
!         local e1-e2-e3 surface coordinate system applies
!
c          xfor(1)=val*skl(1,iref)
c          xfor(2)=val*skl(2,iref)
c          xfor(3)=val*skl(3,iref)
          xfor(1)=val*(e1(1)*skl(1,iref)+
     &                 e1(2)*skl(2,iref)+
     &                 e1(3)*skl(3,iref))
          xfor(2)=val*(e2(1)*skl(1,iref)+
     &                 e2(2)*skl(2,iref)+
     &                 e2(3)*skl(3,iref))
          xfor(3)=val*(e3(1)*skl(1,iref)+
     &                 e3(2)*skl(2,iref)+
     &                 e3(3)*skl(3,iref))
!     
          do j=istart,iend
            idof=int(coeffc(0,j))
            node=int(idof/10.d0)
            if(ntrans.gt.0) then
              if(inotr(1,node).gt.0) then
                write(*,*)
     &'*ERROR in cloaddistributing: no transformation'
                write(*,*) '       is allowed in a node belonging to a'
                write(*,*) '       distributing coupling surface;'
                write(*,*) '       node at fault:',node
                write(*,*)
                ier=1
                cycle
              endif
            endif
            idir=int(idof)-10*node
            forcval=xfor(1)*coeffc(1,j)
     &             +xfor(2)*coeffc(2,j)
     &             +xfor(3)*coeffc(3,j)
!     
            call forcadd(node,idir,forcval,nodeforc,ndirforc,xforc,
     &           nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,
     &           co,ikforc,ilforc,isector,add,user,idefforc,ipompc,
     &           nodempc,nmpc,ikmpc,ilmpc,labmpc)
          enddo
        else
!
!         moment applied in the reference node
!                local e1-e2-e3 surface coordinate system applies
!
          irefm3=iref-3
          xfor(4)=val*(e1(1)*skl(1,irefm3)+
     &                 e1(2)*skl(2,irefm3)+
     &                 e1(3)*skl(3,irefm3))
          xfor(5)=val*(e2(1)*skl(1,irefm3)+
     &                 e2(2)*skl(2,irefm3)+
     &                 e2(3)*skl(3,irefm3))
          xfor(6)=val*(e3(1)*skl(1,irefm3)+
     &                 e3(2)*skl(2,irefm3)+
     &                 e3(3)*skl(3,irefm3))
!     
          do j=istart,iend
            idof=int(coeffc(0,j))
            node=int(idof/10.d0)
            if(ntrans.gt.0) then
              if(inotr(1,node).gt.0) then
                write(*,*)
     &'*ERROR in cloaddistributing: no transformation'
                write(*,*) '       is allowed in a node belonging to a'
                write(*,*) '       distributing coupling surface;'
                write(*,*) '       node at fault:',node
                write(*,*)
                ier=1
                cycle
              endif
            endif
            idir=int(idof)-10*node
            forcval=xfor(4)*coeffc(4,j)
     &             +xfor(5)*coeffc(5,j)
     &             +xfor(6)*coeffc(6,j)
!     
            call forcadd(node,idir,forcval,nodeforc,ndirforc,xforc,
     &           nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,
     &           co,ikforc,ilforc,isector,add,user,idefforc,ipompc,
     &           nodempc,nmpc,ikmpc,ilmpc,labmpc)
          enddo
        endif
      endif
!     
      return
      end
      
