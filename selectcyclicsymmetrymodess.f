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
      subroutine selectcyclicsymmetrymodess(inpc,textpart,cs,ics,
     &  tieset,istartset,
     &  iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,ilmpc,
     &  mpcfree,mcs,set,nset,labmpc,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,nmethod,key,ipoinpc,ier)
!
!     reading the input deck: *SELECT CYCLIC SYMMETRY MODES
!
      implicit none
!
      character*1 inpc(*)
      character*20 labmpc(*)
      character*81 set(*),leftset,tieset(3,*)
      character*132 textpart(16)
!
      integer istep,istat,n,key,i,ns(5),ics(*),istartset(*),ier,
     &  iendset(*),ialset(*),id,ipompc(*),nodempc(3,*),nmpc,nmpc_,
     &  ikmpc(*),ilmpc(*),mpcfree,i1(2),i2(2),i3,i4,i5,j,k,
     &  mpcfreeold,idof,node,ileft,nset,irepeat,ipoinpc(0:*),
     &  mpc,iline,ipol,inl,ipoinp(2,*),inp(3,*),mcs,lprev,ij,nmethod
!
      real*8 coefmpc(*),csab(7),x1(2),x2(2),x3,x4,x5,dd,xn,yn,zn,
     &  cs(17,*)
!
!     irepeat indicates whether the step was preceded by another
!     cyclic symmetry step (irepeat=1) or not (irepeat=0)
!
      data irepeat /0/
      save irepeat
!
      if(istep.eq.0) then
         write(*,*)'*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
         write(*,*)'       *SELECT CYCLIC SYMMETRY MODES'
         write(*,*)'       should be placed within a step definition'
         ier=1
         return
      endif
!
!     check whether in case of cyclic symmetry the frequency procedure
!     is chosen
!
      if((nmethod.ne.2).and.(nmethod.ne.13)) then
         write(*,*) '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
         write(*,*) '       the only valid procedures'
         write(*,*) '       for cyclic symmetry calculations'
         write(*,*) '       with nodal diameters are *FREQUENCY'
         write(*,*) '       and *GREEN'
         ier=1
         return
      endif
!
      ns(2)=0
      ns(3)=0
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NMIN=') then
            read(textpart(i)(6:15),'(i10)',iostat=istat) ns(2)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*SELECT CYCLIC SYMMETRY MODES%",ier)
               return
            endif
         elseif(textpart(i)(1:5).eq.'NMAX=') then
            read(textpart(i)(6:15),'(i10)',iostat=istat) ns(3)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*SELECT CYCLIC SYMMETRY MODES%",ier)
               return
            endif
         else
            write(*,*) '*WARNING reading *SELECT CYCLIC SYMMETRY MODES:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SELECT CYCLIC SYMMETRY MODES%")
         endif
      enddo
!
!     check the input
!
      if(ns(2).lt.0) then
         ns(2)=0
         write(*,*) '*WARNING reading *SELECT CYCLIC SYMMETRY MODES:'
         write(*,*) '         minimum nodal'
         write(*,*) '         diameter must be nonnegative'
      endif
      if(ns(3).lt.ns(2)) then
         write(*,*) '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
         write(*,*) '       maximum nodal'
         write(*,*) '       diameter should not exceed minimal one'
         ier=1
         return
      endif
!
!     loop over all cyclic symmetry parts
!
      do ij=1,mcs
         ns(1)=int(cs(1,ij))
         ns(4)=int(cs(4,ij))
         leftset=tieset(2,int(cs(17,ij)))
         lprev=int(cs(14,ij))
         do i=1,7
            csab(i)=cs(5+i,ij)
         enddo
!
!     check whether cyclic symmetry axis is part of the structure
!
c         do i=1,nset
c            if(set(i).eq.leftset) exit
c         enddo
c         ileft=i
         call cident81(set,leftset,nset,id)
         ileft=nset+1
         if(id.gt.0) then
           if(leftset.eq.set(id)) then
             ileft=id
           endif
         endif
!     
!     if this step was preceded by a cyclic symmetry step:
!     check for MPC's for nodes on the cyclic symmetry axis
!     and delete them
!     
         if(irepeat.eq.1) then
            do i=1,ns(4)
               node=ics(lprev+i)
               if(node.lt.0) then
                  node=-node
                  do k=1,3
                     idof=8*(node-1)+k
                     call nident(ikmpc,idof,nmpc,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
c                           write(*,*) 'removing MPC',node,k
                           mpc=ilmpc(id)
                           call mpcrem(mpc,mpcfree,nodempc,nmpc,ikmpc,
     &                           ilmpc,labmpc,coefmpc,ipompc)
                        endif
                     endif
                  enddo
               endif
            enddo
         endif
!     
         do i=1,ns(4)
            node=ics(lprev+i)
            if(node.lt.0) then
               node=-node
               if(ns(2).ne.ns(3)) then
                  if((ns(2).eq.0).or.(ns(2).eq.1)) then
                     write(*,*) '*ERROR: axis of cyclic symmetry'
                     write(*,*) '        is part of the structure;'
                     write(*,*) '        nodal diameters 0, 1, and'
                     write(*,*) '        those above must be each in'
                     write(*,*) '        separate steps.'
                     ier=1
                     return
                  endif
               endif
!     
!     specifying special MPC's for nodes on the axis
!     
!     normal along the axis
!     
               xn=csab(4)-csab(1)
               yn=csab(5)-csab(2)
               zn=csab(6)-csab(3)
               dd=dsqrt(xn*xn+yn*yn+zn*zn)
               xn=xn/dd
               yn=yn/dd
               zn=zn/dd
!     
!     nodal diameter 0
!     
               if(ns(2).eq.0) then
                  if(dabs(xn).gt.1.d-10) then
                     i1(1)=2
                     i1(2)=3
                     i2(1)=1
                     i2(2)=1
                     x1(1)=xn
                     x1(2)=xn
                     x2(1)=-yn
                     x2(2)=-zn
                  elseif(dabs(yn).gt.1.d-10) then
                     i1(1)=1
                     i1(2)=3
                     i2(1)=2
                     i2(2)=2
                     x1(1)=yn
                     x1(2)=yn
                     x2(1)=-xn
                     x2(2)=-zn
                  elseif(dabs(zn).gt.1.d-10) then
                     i1(1)=1
                     i1(2)=2
                     i2(1)=3
                     i2(2)=3
                     x1(1)=zn
                     x1(2)=zn
                     x2(1)=-xn
                     x2(2)=-yn
                  endif
!     
!     generating two MPC's expressing that the nodes cannot
!     move in planes perpendicular to the cyclic symmetry
!     axis
!     
                  do k=1,2
                     idof=8*(node-1)+i1(k)
                     call nident(ikmpc,idof,nmpc,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           write(*,*) 
     &                 '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                           write(*,*) '       node',node,
     &                          ' on cyclic symmetry'
                           write(*,*) '       axis is used in other MPC'
                           ier=1
                           return
                        endif
                     endif
                     nmpc=nmpc+1
                     ipompc(nmpc)=mpcfree
                     labmpc(nmpc)='                    '
!     
!     updating ikmpc and ilmpc
!     
                     do j=nmpc,id+2,-1
                        ikmpc(j)=ikmpc(j-1)
                        ilmpc(j)=ilmpc(j-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
!     
                     nodempc(1,mpcfree)=node
                     nodempc(2,mpcfree)=i1(k)
                     coefmpc(mpcfree)=x1(k)
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                  '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                        write(*,*) '       increase memmpc_'
                        ier=1
                        return
                     endif
                     nodempc(1,mpcfree)=node
                     nodempc(2,mpcfree)=i2(k)
                     coefmpc(mpcfree)=x2(k)
                     mpcfreeold=mpcfree
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                  '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                        write(*,*) '       increase memmpc_'
                        ier=1
                        return
                     endif
                     nodempc(3,mpcfreeold)=0
                  enddo
               elseif(ns(2).eq.1) then
!     
!     nodal diameter 1
!     
                  if(dabs(xn).gt.1.d-10) then
                     i3=1
                     i4=2
                     i5=3
                     x3=xn
                     x4=yn
                     x5=zn
                  elseif(dabs(yn).gt.1.d-10) then
                     i3=2
                     i4=2
                     i5=3
                     x3=yn
                     x4=xn
                     x5=zn
                  else
                     i3=3
                     i4=1
                     i5=2
                     x3=zn
                     x4=xn
                     x5=yn
                  endif
!     
!     generating one MPC expressing that the nodes should
!     not move along the axis
!     
                  idof=8*(node-1)+i3
                  call nident(ikmpc,idof,nmpc,id)
                  if(id.gt.0) then
                     if(ikmpc(id).eq.idof) then
                        write(*,*) 
     &                  '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                        write(*,*) '       node',node,
     &                       ' on cyclic symmetry'
                        write(*,*) '       axis is used in other MPC'
                        ier=1
                        return
                     endif
                  endif
                  nmpc=nmpc+1
                  ipompc(nmpc)=mpcfree
                  labmpc(nmpc)='                    '
!     
!     updating ikmpc and ilmpc
!     
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
!     
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=i3
                  coefmpc(mpcfree)=x3
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &               '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                     write(*,*) '       increase memmpc_'
                     ier=1
                     return
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=i4
                  coefmpc(mpcfree)=x4
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &               '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                     write(*,*) '       increase memmpc_'
                     ier=1
                     return
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=i5
                  coefmpc(mpcfree)=x5
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &               '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                     write(*,*) '       increase memmpc_'
                     ier=1
                     return
                  endif
                  nodempc(3,mpcfreeold)=0
               else
                  do k=1,3
                     idof=8*(node-1)+k
                     call nident(ikmpc,idof,nmpc,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           write(*,*) 
     &                   '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                           write(*,*) '       node',node,
     &                          ' on cyclic symmetry'
                           write(*,*) '       axis is used in other MPC'
                           ier=1
                           return
                        endif
                     endif
                     nmpc=nmpc+1
                     ipompc(nmpc)=mpcfree
                     labmpc(nmpc)='                    '
!     
!     updating ikmpc and ilmpc
!     
                     do j=nmpc,id+2,-1
                        ikmpc(j)=ikmpc(j-1)
                        ilmpc(j)=ilmpc(j-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
!     
                     nodempc(1,mpcfree)=node
                     nodempc(2,mpcfree)=k
                     coefmpc(mpcfree)=1.d0
                     mpcfreeold=mpcfree
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                  '*ERROR reading *SELECT CYCLIC SYMMETRY MODES:'
                        write(*,*) '       increase memmpc_'
                        ier=1
                        return
                     endif
                     nodempc(3,mpcfreeold)=0
                  enddo
               endif
            endif
         enddo
!     
         cs(2,ij)=ns(2)+0.5d0
         cs(3,ij)=ns(3)+0.5d0
      enddo
!     
      if(irepeat.eq.0) irepeat=1
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!    
c      do j=1,nmpc
c         call writempc(ipompc,nodempc,coefmpc,labmpc,j)
c      enddo
!
      return
      end

