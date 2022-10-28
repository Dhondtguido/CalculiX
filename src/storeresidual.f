!
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
      subroutine storeresidual(nactdof,b,fn,filab,ithermal,nk,sti,stn,
     &  ipkon,inum,kon,lakon,ne,mi,orab,ielorien,co,
     &  itg,ntg,vold,ielmat,thicke,ielprop,prop)
!
!     This routine is called in case of divergence:
!     stores the residual forces in fn and changes the
!     file storage labels so that the independent
!     variables (displacements and/or temperatures) and
!     the corresponding residual forces are stored in the
!     frd file
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*)
      character*87 filab(*)
!
      integer mi(*),nactdof(0:mi(2),*),ithermal(*),i,j,nk,
     &  nfield,ndim,iorienglob,icfdout,ielmat(mi(3),*),
     &  ipkon(*),inum(*),kon(*),ielprop(*),
     &  ne,ielorien,itg(*),ntg,mt,nlabel
!
      real*8 b(*),fn(0:mi(2),*),sti(6,mi(1),*),stn(6,*),orab(7,*),
     &  co(3,*),vold(0:mi(2),*),thicke(mi(3),*),prop(*)
!
      mt=mi(2)+1
!
      nlabel=55
!
!     storing the residual forces in field fn
!
      do i=1,nk
         do j=0,mi(2)
            if(nactdof(j,i).gt.0) then
               fn(j,i)=b(nactdof(j,i))
            else
               fn(j,i)=0.d0
            endif
         enddo
      enddo
!
!     adapting the storage labels
!
      do i=1,nlabel
         filab(i)(1:4)='    '
      enddo
!
      if(ithermal(1).ne.2) then
         filab(1)(1:3)='U  '
         filab(5)(1:4)='RF  '
      else
         filab(1)(1:4)='    '
         filab(5)(1:4)='    '
      endif
!
      if(ithermal(1).gt.1) then
         filab(2)(1:4)='NT  '
         filab(10)(1:4)='RFL '
         filab(14)(1:4)='TT  '
         filab(15)(1:4)='MF  '
         filab(16)(1:4)='TP  '
         filab(17)(1:4)='ST  '
      else
         filab(2)(1:4)='    '
         filab(10)(1:4)='    '
         filab(14)(1:4)='    '
         filab(15)(1:4)='    '
         filab(16)(1:4)='    '
         filab(17)(1:4)='    '
      endif
!
!     calculating inum
!
      nfield=0
      ndim=0
      iorienglob=0
      cflag=filab(1)(5:5)
      icfdout=0
      force=.false.
      call extrapolate(sti,stn,ipkon,inum,kon,lakon,nfield,nk,
     &     ne,mi(1),ndim,orab,ielorien,co,iorienglob,cflag,
     &     vold,force,ielmat,thicke,ielprop,prop)
!
      if(ithermal(1).gt.1) then
         call networkextrapolate(vold,ipkon,inum,kon,lakon,ne,mi)
      endif
!
!     interpolation for 1d/2d elements
!
      if(filab(1)(5:5).eq.'I') then
         nfield=mt
         cflag=filab(1)(5:5)
         force=.false.
         call map3dto1d2d(vold,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,cflag,co,vold,force,mi,ielprop,prop)
      endif
!
!     marking gas nodes by multiplying inum by -1
!
      do i=1,ntg
         inum(itg(i))=-inum(itg(i))
      enddo
!
      return
      end


