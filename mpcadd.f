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
      subroutine mpcadd(nodedep,is,ie,nboun,ipompc,nodempc,
     &  coefmpc,nmpc,nmpc_,mpcfree,orab,ikboun,ikmpc,ilmpc,co,labmpc,
     &  label,nodeind,iorientation)
!
!     generates an equality MPC between node "nodedep" (dependent) and
!     node "nodeind" (independent), both at the same position, 
!     in orientation system "iorientation" for local degrees of
!     freedom "is" up to "ie".
!
      implicit none
!
      character*20 labmpc(*),label
!
      integer nodedep,is,ie,nboun,i,j,ipompc(*),nodempc(3,*),nmpc,nmpc_,
     &  mpcfree,ikboun(*),ikmpc(*),ilmpc(*),idof,number,id,
     &  mpcfreeold,three,kflag,iy(3),inumber,nodeind,iorientation
!
      real*8 coefmpc(*),a(3,3),co(3,*),orab(7,*),dx(3),p(3)
!
      loop: do i=is,ie
!
         if(iorientation.eq.0) then
!
!        no transformation applies: simple SPC
!
            idof=8*(nodedep-1)+i
            call nident(ikboun,idof,nboun,id)
            if(id.gt.0) then
               if(ikboun(id).eq.idof) then
                  cycle loop
               endif
            endif
            call nident(ikmpc,idof,nmpc,id)
            if(id.ne.0) then
               if(ikmpc(id).eq.idof) cycle loop
            endif
!
!           new MPC
!
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) '*ERROR in mpcadd: increase nmpc_'
               call exit(201)
            endif
            labmpc(nmpc)=label
            ipompc(nmpc)=mpcfree
            do j=nmpc,id+2,-1
               ikmpc(j)=ikmpc(j-1)
               ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!
            nodempc(1,mpcfree)=nodedep
            nodempc(2,mpcfree)=i
            coefmpc(mpcfree)=1.d0
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*) '*ERROR in mpcadd: increase memmpc_'
               call exit(201)
            endif
!
            nodempc(1,mpcfree)=nodeind
            nodempc(2,mpcfree)=i
            coefmpc(mpcfree)=-1.d0
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*) '*ERROR in mpcadd: increase memmpc_'
               call exit(201)
            endif
            nodempc(3,mpcfreeold)=0
!
         else
!
!        transformation applies
!
            do j=1,3
               p(j)=co(j,nodedep)
            enddo
            call transformatrix(orab(1,iorientation),p,a)
!
!           new mpc
!
            iy(1)=1
            iy(2)=2
            iy(3)=3
            dx(1)=dabs(a(1,i))
            dx(2)=dabs(a(2,i))
            dx(3)=dabs(a(3,i))
            three=3
            kflag=-2
            call dsort(dx,iy,three,kflag)
            do inumber=1,3
               number=iy(inumber)
               idof=8*(nodedep-1)+number
               call nident(ikmpc,idof,nmpc,id)
               if(id.ne.0) then
                  if(ikmpc(id).eq.idof) cycle
               endif
               if(dabs(a(number,i)).lt.1.d-5) cycle
               nmpc=nmpc+1
               if(nmpc.gt.nmpc_) then
                  write(*,*) '*ERROR in mpcadd: increase nmpc_'
                  call exit(201)
               endif
               labmpc(nmpc)=label
               ipompc(nmpc)=mpcfree
               do j=nmpc,id+2,-1
                  ikmpc(j)=ikmpc(j-1)
                  ilmpc(j)=ilmpc(j-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
               exit
            enddo
!
!           check whether a dependent term was found; if none was
!           found this can be due to the fact that:
!           - all dofs were used by other MPC's
!           - the MPC coefficients were too small
!           - or a combination of both
!
            if(inumber.gt.3) cycle
!
            inumber=inumber-1
            do j=1,3
               inumber=inumber+1
               if(inumber.gt.3) inumber=1
               number=iy(inumber)
               if(dabs(a(number,i)).lt.1.d-30) cycle
!
               nodempc(1,mpcfree)=nodedep
               nodempc(2,mpcfree)=number
               coefmpc(mpcfree)=a(number,i)
               mpcfree=nodempc(3,mpcfree)
               if(mpcfree.eq.0) then
                  write(*,*) '*ERROR in mpcadd: increase memmpc_'
                  call exit(201)
               endif
!
               nodempc(1,mpcfree)=nodeind
               nodempc(2,mpcfree)=number
               coefmpc(mpcfree)=-a(number,i)
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
               if(mpcfree.eq.0) then
                  write(*,*) '*ERROR in mpcadd: increase memmpc_'
                  call exit(201)
               endif
            enddo
!
            nodempc(3,mpcfreeold)=0
         endif
      enddo loop
!
      return
      end
