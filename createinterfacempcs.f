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
      subroutine createinterfacempcs(imastnode,xmastnor,nmastnode,
     &  ikmpc,ilmpc,nmpc,ipompc,nodempc,coefmpc,labmpc,mpcfree,ikboun,
     &  nboun)
!
      character*20 labmpc(*)
!
      integer imastnode(*),nmastnode,ikmpc(*),ilmpc(*),nmpc,
     &  nodempc(3,*),kflag,i,j,node,lnor(3),three,k,id,id1,
     &  ikboun(*),nboun,mpcfree,m,ipompc(*),mpcfreeold
!
      real*8 xmastnor(3,*),coefmpc(*),xnor(3)
!
      kflag=-2
      three=3
!
      loop: do i=1,nmastnode
         node=imastnode(i)
!
!        sorting the components of the normal in the node
!
         do j=1,3
            xnor(j)=xmastnor(j,i)
            lnor(j)=j
         enddo
         call dsort(xnor,lnor,three,kflag)
!
         do k=1,3
            j=lnor(k)
            idof=8*(node-1)+j
            call nident(ikmpc,idof,nmpc,id)
            if(id.gt.0) then
               if(ikmpc(id).eq.idof)cycle
            endif
            call nident(ikboun,idof,nboun,id1)
            if(id1.gt.0) then
               if(ikboun(id1).eq.idof)cycle
            endif
            if(dabs(xnor(k)).lt.1.d-20)cycle
!
!           create a MPC corresponding to A.n=0
!            
            nmpc=nmpc+1
            labmpc(nmpc)='                    '
            ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
            do m=nmpc,id+2,-1
               ikmpc(m)=ikmpc(m-1)
               ilmpc(m)=ilmpc(m-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!     
            nodempc(1,mpcfree)=node
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=xmastnor(j,i)
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*)
     &            '*ERROR in createinterfacempcs: increase memmpc_'
               call exit(201)
            endif
!     
            j=j+1
            if(j.gt.3)j=1
            nodempc(1,mpcfree)=node
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=xmastnor(j,i)
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*)
     &            '*ERROR in createinterfacempcs: increase memmpc_'
               call exit(201)
            endif
!     
            j=j+1
            if(j.gt.3)j=1
            nodempc(1,mpcfree)=node
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=xmastnor(j,i)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*)
     &            '*ERROR in createinterfacempcs: increase memmpc_'
               call exit(201)
            endif
            nodempc(3,mpcfreeold)=0
!
            cycle loop
         enddo
!
         write(*,*) '*WARNING in createinterfacempcs: no A.n MPC'
         write(*,*) '         created for node ',node
!
      enddo loop
!     
      return
      end
