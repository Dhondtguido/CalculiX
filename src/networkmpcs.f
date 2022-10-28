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
      subroutine networkmpcs(inpc,textpart,ipompc,nodempc,coefmpc,
     &  nmpc,nmpc_,mpcfree,nk,ikmpc,ilmpc,
     &  labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *NETWORK MPC
!
      implicit none
!
      character*1 inpc(*)
      character*13 type
      character*20 labmpc(*)
      character*132 textpart(16)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,istep,istat,
     &  n,i,j,ii,key,nterm,nk,node,ndir,mpcfreeold,ikmpc(*),ilmpc(*),
     &  id,idof,iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),m,
     &  ier
!
      real*8 coefmpc(*),x
!
      do m=2,n
         if(textpart(m)(1:5).eq.'TYPE=') then
            type(1:13)=textpart(m)(6:18)
            if(textpart(m)(19:19).ne.' ') then
               write(*,*) '*ERROR reading *NETWORK MPC: type'
               write(*,*) '       of network mpc is too long'
               call inputwarning(inpc,ipoinpc,iline,
     &"*NETWORK MPC%")
             endif
         else
             write(*,*) 
     &        '*WARNING reading *NETWORK MPC: parameter not recognized:'
             write(*,*) '         ',
     &            textpart(m)(1:index(textpart(m),' ')-1)
             call inputwarning(inpc,ipoinpc,iline,
     &"*NETWORK MPC%")
          endif
       enddo
!
      if(istep.gt.0) then
         write(*,*) 
     &     '*ERROR reading *NETWORK MPC: *NETWORK MPC should be placed'
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
            write(*,*) '*ERROR reading *NETWORK MPC: increase nmpc_'
            ier=1
            return
         endif
!
         labmpc(nmpc)(1:7)='NETWORK'
         labmpc(nmpc)(8:20)=type(1:13)
         ipompc(nmpc)=mpcfree
         ii=0
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) 
     &           '*ERROR reading *NETWORK MPC: mpc definition ',nmpc
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*NETWORK MPC%",ier)
               return
            endif
!
            do i=1,n/3
!
               read(textpart((i-1)*3+1)(1:10),'(i10)',iostat=istat) node
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*NETWORK MPC%",ier)
                  return
               endif
               if((node.gt.nk).or.(node.le.0)) then
                  write(*,*) '*ERROR reading *NETWORK MPC:'
                  write(*,*) '       node ',node,' is not defined'
                  ier=1
                  return
               endif
!
               read(textpart((i-1)*3+2)(1:10),'(i10)',iostat=istat) ndir
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*NETWORK MPC%",ier)
                  return
               endif
               if(ndir.le.6) then
               elseif(ndir.eq.8) then
                  ndir=4
               elseif(ndir.eq.11) then
                  ndir=0
               else
                  write(*,*) '*ERROR reading *NETWORK MPC:'
                  write(*,*) '       direction',ndir,' is not defined'
                  ier=1
                  return
               endif
!
               read(textpart((i-1)*3+3)(1:20),'(f20.0)',iostat=istat) x
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*NETWORK MPC%",ier)
                  return
               endif
!
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=ndir
               coefmpc(mpcfree)=x
!
!                 updating ikmpc and ilmpc
!
               if(ii.eq.0) then
                  idof=8*(node-1)+ndir
                  call nident(ikmpc,idof,nmpc-1,id)
                  if(id.gt.0) then
                     if(ikmpc(id).eq.idof) then
                        write(*,100)
     &                       (ikmpc(id))/8+1,ikmpc(id)-8*((ikmpc(id))/8)
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
     &                  '*ERROR reading *NETWORK MPC: increase memmpc_'
                  ier=1
                  return
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
 100  format(/,'*ERROR reading *NETWORK MPC: the DOF corresponding to',
     &           /,'node ',i10,' in direction',i1,' is detected on',
     &           /,'the dependent side of two different MPC''s') 
      return
      end

