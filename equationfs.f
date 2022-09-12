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
      subroutine equationfs(inpc,textpart,ipompc,nodempc,coefmpc,
     &  nmpc,nmpc_,mpcfree,co,trab,ntrans,ikmpc,ilmpc,
     &  labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &  lakon,ne,nload,sideload,ipkon,kon,nelemload,ier)
!
!     reading the input deck: *EQUATIONF
!
      implicit none
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*),label,sideload(*)
      character*132 textpart(16)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,istep,istat,
     &  n,i,j,ii,key,nterm,number,ntrans,ndir,indexe,
     &  mpcfreeold,ikmpc(*),ilmpc(*),id,idof,itr,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier,
     &  k,m,iface,ifacel,nelem,ifaceq(8,6),ifacet(6,4),ifacew(8,5),
     &  loadid,ne,nope,nopes,ipkon(*),nload,nelemload(2,*),kon(*)
!
      real*8 coefmpc(*),co(3,*),trab(7,*),a(3,3),x,cg(3)
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      do m=2,n
         write(*,*) 
     &        '*WARNING reading *EQUATIONF: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(m)(1:index(textpart(m),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*EQUATIONF%")
      enddo
!
      if(istep.gt.0) then
         write(*,*) 
     &     '*ERROR reading *EQUATIONF: *EQUATIONF should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         read(textpart(1)(1:10),'(i10)',iostat=istat) nterm
!
         nmpc=nmpc+1
         if(nmpc.gt.nmpc_) then
            write(*,*) '*ERROR reading *EQUATIONF: increase nmpc_'
            ier=1
            return
         endif
!
         labmpc(nmpc)='FLUID               '
         ipompc(nmpc)=mpcfree
         ii=0
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) '*ERROR reading *EQUATIONF: mpc definition ',
     &             nmpc
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*EQUATIONF%",ier)
               return
            endif
!
            do i=1,n/4
!
               read(textpart((i-1)*4+1)(1:10),'(i10)',iostat=istat)nelem
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*EQUATIONF%",ier)
                  return
               endif
               if((nelem.gt.ne).or.(nelem.le.0)) then
                  write(*,*) '*ERROR reading *EQUATIONF:'
                  write(*,*) '       element ',nelem,' is not defined'
                  ier=1
                  return
               endif
!
               read(textpart((i-1)*4+2)(2:2),'(i1)',iostat=istat) 
     &             ifacel
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*EQUATIONF%",ier)
                  return
               endif
               iface=10*nelem+ifacel
!
               read(textpart((i-1)*4+3)(1:10),'(i10)',iostat=istat) ndir
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*EQUATIONF%",ier)
                  return
               endif
               if(ndir.le.3) then
               elseif(ndir.eq.4) then
                  write(*,*) '*ERROR reading *EQUATIONF: an equation'
                  write(*,*) '       on DOF 4 is not allowed'
                  ier=1
                  return
               elseif(ndir.eq.5) then
                  write(*,*) '*ERROR reading *EQUATIONF: an equation'
                  write(*,*) '       on DOF 5 is not allowed'
                  ier=1
                  return
               elseif(ndir.eq.6) then
                  write(*,*) '*ERROR reading *EQUATIONF: an equation'
                  write(*,*) '       on DOF 6 is not allowed'
                  ier=1
                  return
               elseif(ndir.eq.8) then
                  ndir=4
               elseif(ndir.eq.11) then
                  ndir=0
               else
                  write(*,*) '*ERROR reading *EQUATIONF:'
                  write(*,*) '       direction',ndir,' is not defined'
                  ier=1
                  return
               endif
!
               read(textpart((i-1)*4+4)(1:20),'(f20.0)',iostat=istat) x
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*EQUATIONF%",ier)
                  return
               endif
!
!              check whether the face is transformed
!
               if(ntrans.le.0) then
                  itr=0
               else
                  label(1:20)='T                   '
                  write(label(2:2),'(i1)') ifacel
                  call identifytransform(nelem,label,nelemload,sideload,
     &                 nload,loadid)
                  if(loadid.eq.0) then
                     itr=0
                  else
                     itr=nelemload(2,loadid)
                  endif
               endif
!
               if((itr.eq.0).or.(ndir.eq.0).or.(ndir.eq.4)) then
                  nodempc(1,mpcfree)=iface
                  nodempc(2,mpcfree)=ndir
                  coefmpc(mpcfree)=x
!
!                 updating ikmpc and ilmpc
!
                  if(ii.eq.0) then
                     idof=-(8*(iface-1)+ndir)
                     call nident(ikmpc,idof,nmpc-1,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           write(*,100)
     &                   (ikmpc(id))/8+1,ikmpc(id)-8*((ikmpc(id))/8)
                           ier=1
                           return
                        endif
                     endif
                     do j=nmpc,id+2,-1
                        ikmpc(j)=ikmpc(j-1)
                        ilmpc(j)=ilmpc(j-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
                  endif
!
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                 '*ERROR reading *EQUATIONF: increase memmpc_'
                     ier=1
                     return
                  endif
               else
!
!        transformation applies
!
!        number of nodes belonging to the face
!
                  if(lakon(nelem)(4:4).eq.'8') then
                     nope=8
                     nopes=4
                  elseif(lakon(nelem)(4:4).eq.'4') then
                     nope=4
                     nopes=3
                  elseif(lakon(nelem)(4:4).eq.'6') then
                     nope=6
                     if(ifacel.le.2) then
                        nopes=3
                     else
                        nopes=4
                     endif
                  endif
!     
!     determining the center of gravity 
!     
                  do j=1,3
                     cg(j)=0.d0
                  enddo
!     
                  indexe=ipkon(nelem)
                  if(nope.eq.8) then
                     do k=1,nopes
                        do j=1,3
                           cg(j)=cg(j)+
     &                           co(j,kon(indexe+ifaceq(k,ifacel)))
                        enddo
                     enddo
                  elseif(nope.eq.4) then
                     do k=1,nopes
                        do j=1,3
                           cg(j)=cg(j)+
     &                           co(j,kon(indexe+ifacet(k,ifacel)))
                        enddo
                     enddo
                  else
                     do k=1,nopes
                        do j=1,3
                           cg(j)=cg(j)+
     &                           co(j,kon(indexe+ifacew(k,ifacel)))
                        enddo
                     enddo
                  endif
                  do j=1,3
                     cg(j)=cg(j)/nopes
                  enddo
!     
!     determining the transformation coefficients at the center
!     of gravity
!     
                  call transformatrix(trab(1,itr),
     &                 cg,a)
!
                  number=ndir-1
                  if(ii.eq.0) then
!
!                    determining which direction to use for the
!                    dependent side: should not occur on the dependent
!                    side in another MPC and should have a nonzero
!                    coefficient
!
                     do j=1,3
                        number=number+1
                        if(number.gt.3) number=1
                        idof=-(8*(iface-1)+number)
                        call nident(ikmpc,idof,nmpc-1,id)
                        if(id.gt.0) then
                           if(ikmpc(id).eq.idof) then
                              cycle
                           endif
                        endif
                        if(dabs(a(number,ndir)).lt.1.d-5) cycle
                        exit
                     enddo
                     if(j.gt.3) then
                        write(*,*) 
     &                       '*ERROR reading *EQUATIONF: MPC on face'
                        write(*,*) ifacel,' of element',nelem
                        write(*,*) ' in transformed coordinates'
                        write(*,*) ' cannot be converted in MPC: all'
                        write(*,*) ' DOFs in the node are used as'
                        write(*,*) ' dependent nodes in other MPCs'
                        ier=1
                        return
                     endif
                     number=number-1
!
!                    updating ikmpc and ilmpc
!
                     do j=nmpc,id+2,-1
                        ikmpc(j)=ikmpc(j-1)
                        ilmpc(j)=ilmpc(j-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
                  endif
!
                  do j=1,3
                     number=number+1
                     if(number.gt.3) number=1
                     if(dabs(a(number,ndir)).lt.1.d-5) cycle
                     nodempc(1,mpcfree)=iface
                     nodempc(2,mpcfree)=number
                     coefmpc(mpcfree)=x*a(number,ndir)
                     mpcfreeold=mpcfree
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                    '*ERROR reading *EQUATIONF: increase memmpc_'
                        ier=1
                        return
                     endif
                  enddo
               endif
!
               ii=ii+1
            enddo
!
            if(ii.eq.nterm) then
               nodempc(3,mpcfreeold)=0
               exit
            endif
         enddo
      enddo
!
 100  format(/,'*ERROR reading *EQUATIONF: the DOF corresponding to',
     &           /,'iface ',i1,' of element',i10,' in direction',
     &           i5,' is detected on',
     &           /,'the dependent side of two different MPC''s') 
      return
      end

