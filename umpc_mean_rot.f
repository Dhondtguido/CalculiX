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
      subroutine umpc_mean_rot(x,u,f,a,jdof,n,force,iit,idiscon,
     &          jnode,ikmpc,nmpc,ikboun,nboun,label)
!
!     updates the coefficients in a mean rotation mpc
!
!     INPUT:
!
!     x(3,1..n)          Carthesian coordinates of the nodes in the
!                        user mpc.
!     u(3,1..n)          Actual displacements of the nodes in the
!                        user mpc.     
!     jdof(1..n)         Actual degrees of freedom of the mpc terms
!     n                  number of terms in the user mpc
!     force              Actual value of the mpc force
!     iit                iteration number
!
!     OUTPUT:
!
!     f                  Actual value of the mpc. If the mpc is
!                        exactly satisfied, this value is zero
!     a(n)               coefficients of the linearized mpc
!     jdof(1..n)         Corrected degrees of freedom of the mpc terms
!     idiscon            0: no discontinuity
!                        1: discontinuity
!                        If a discontinuity arises the previous
!                        results are not extrapolated at the start of
!                        a new increment
!
      implicit none
!
      character*20 label
!
      integer jdof(*),n,nkn,i,j,k,imax,iit,idiscon,node,nodeold,
     &  nodemax,idirold,idirmax,jnode(*),idof,id,iold(3),inew(3),
     &  ikmpc(*),nmpc,ikboun(*),nboun,idir,nendnode
!
      real*8 x(3,*),u(3,*),f,a(*),aa(3),cgx(3),cgu(3),pi(3),b(3),
     &  xi(3),dd,al,a1,amax,c1,c2,c3,c4,c9,c10,force,xcoef,
     &  transcoef(3),stdev
!
      nkn=(n-1)/3
      if(3*nkn.ne.n-1) then
         write(*,*)
     &     '*ERROR in meanrotmpc: MPC has wrong number of terms'
         call exit(201)
      endif
!
!     normal along the rotation axis
!
      dd=0.d0
      do i=1,3
         aa(i)=x(i,n)
         dd=dd+aa(i)**2
      enddo
      dd=dsqrt(dd)
      if(dd.lt.1.d-10) then
         write(*,*) 
     &     '*ERROR in meanrotmpc: rotation vector has zero length'
         call exit(201)
      endif
      do i=1,3
         aa(i)=aa(i)/dd
      enddo
c      write(*,*) 'umpc_mean_rot aa',(aa(i),i=1,3)
!
!     finding the center of gravity of the position and the
!     displacements of the nodes involved in the MPC
!
      do i=1,3
         cgx(i)=0.d0
         cgu(i)=0.d0
      enddo
!
c      write(*,*) 'umpc_mean_rot, nkn',nkn
      do i=1,nkn
         do j=1,3
            cgx(j)=cgx(j)+x(j,3*i-2)
            cgu(j)=cgu(j)+u(j,3*i-2)
         enddo
c         write(*,*) 'umpc_mean_rot x',i,(x(j,3*i-2),j=1,3)
c         write(*,*) 'umpc_mean_rot u',i,(u(j,3*i-2),j=1,3)
      enddo
!
      do i=1,3
         cgx(i)=cgx(i)/nkn
         cgu(i)=cgu(i)/nkn
      enddo
c      write(*,*) 'umpc_mean_rot cgx',(cgx(i),i=1,3)
!
!     calculating a standard deviation; this quantity will
!     serve as a limit for checking the closeness of individual
!     nodes to the center of gravity
!
      stdev=0.d0
      do i=1,nkn
         do j=1,3
            stdev=stdev+(x(j,3*i-2)-cgx(j))**2
         enddo
      enddo
      stdev=stdev/nkn
c      write(*,*) 'umpc_mean_rot stdev',stdev
!
!     initializing a
!
      do i=1,n
         a(i)=0.d0
      enddo
!
!     calculating the partial derivatives and storing them in a
!
      f=0.d0
      do i=1,nkn
!
!        relative positions
!
         do j=1,3
            pi(j)=x(j,3*i-2)-cgx(j)
            xi(j)=u(j,3*i-2)-cgu(j)+pi(j)
         enddo
c         write(*,*) 'umpc_mean_rot pi',(pi(j),j=1,3)
!
!              projection on a plane orthogonal to the rotation vector
!
         c1=pi(1)*aa(1)+pi(2)*aa(2)+pi(3)*aa(3)
         do j=1,3
            pi(j)=pi(j)-c1*aa(j)
         enddo
c         write(*,*) 'umpc_mean_rot c1 pi',c1,(pi(j),j=1,3)
!
         c1=xi(1)*aa(1)+xi(2)*aa(2)+xi(3)*aa(3)
         do j=1,3
            xi(j)=xi(j)-c1*aa(j)
         enddo
!
         c1=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)
c         if(c1.lt.1.d-20) then
         if(c1.lt.stdev*1.d-10) then
            if(label(8:9).ne.'BS') then
               write(*,*) '*WARNING in meanrotmpc: node ',jnode(3*i-2)
               write(*,*) '         is very close to the '
               write(*,*) '         rotation axis through the'
               write(*,*) '         center of gravity of'
               write(*,*) '         the nodal cloud in a'
               write(*,*) '         mean rotation MPC.'
               write(*,*) '         This node is not taken'
               write(*,*) '         into account in the MPC'
            endif
            cycle
         endif
         c3=xi(1)*xi(1)+xi(2)*xi(2)+xi(3)*xi(3)
         c2=dsqrt(c1*c3)
!
         al=(aa(1)*pi(2)*xi(3)+aa(2)*pi(3)*xi(1)+aa(3)*pi(1)*xi(2)
     &     -aa(3)*pi(2)*xi(1)-aa(1)*pi(3)*xi(2)-aa(2)*pi(1)*xi(3))
     &     /c2
!
         f=f+dasin(al)
!
         do j=1,3
            if(j.eq.1) then
               c4=aa(2)*pi(3)-aa(3)*pi(2)
            elseif(j.eq.2) then
               c4=aa(3)*pi(1)-aa(1)*pi(3)
            else
               c4=aa(1)*pi(2)-aa(2)*pi(1)
            endif
            c9=(c4/c2-al*xi(j)/c3)/dsqrt(1.d0-al*al)
c            write(*,*) 'umpc_mean_rot c4,c9',j,c4,c9
!
            do k=1,nkn
               if(i.eq.k) then
                  c10=c9*(1.d0-1.d0/real(nkn))
               else
                  c10=-c9/real(nkn)
               endif
               a(k*3-3+j)=a(k*3-3+j)+c10
c               write(*,*) 'umpc_mean_rot c10',j,c10,a(k*3-3+j)
            enddo
         enddo
      enddo
c      do j=1,n
c         write(*,*) 'umpc_mean_rot a ',a(j)
c      enddo
      a(n)=-nkn
      f=f-nkn*u(1,n)
!
!     assigning the degrees of freedom
!     the first three dofs may not be in ascending order
!
      do i=1,3
         b(i)=a(i)
      enddo
      do i=1,3
         a(i)=b(jdof(i))
      enddo
!
      do i=4,nkn
         jdof(i*3-2)=1
         jdof(i*3-1)=2
         jdof(i*3)=3
      enddo
!
      jdof(n)=1
!
!     looking for the maximum tangent to decide which DOF should be
!     taken to be the dependent one
!
      if(dabs(a(1)).lt.1.d-5) then
!
!        special treatment of beam and shell mean rotation MPC's
!        cf. usermpc.f
!
         if(label(8:9).eq.'BS') then
            do i=1,3
               transcoef(i)=a(n-4+i)
            enddo
         endif
!
         imax=0
         amax=1.d-5
!
         if(label(8:9).eq.'BS') then
            nendnode=n-4
         else
            nendnode=n-1
         endif
         do i=2,nendnode
            if(dabs(a(i)).gt.amax) then
               idir=jdof(i)
!
               if(label(8:9).eq.'BS') then
                  if(a(i)*transcoef(idir).gt.0) cycle
               endif
!
               node=jnode(i)
               idof=8*(node-1)+idir
!
!              check whether dof is not used in another MPC
!
               call nident(ikmpc,idof,nmpc,id)
               if(id.gt.0) then
                  if(ikmpc(id).eq.idof) cycle
               endif
!
!              check whether dof is not used in a SPC
!
               call nident(ikboun,idof,nboun,id)
               if(id.gt.0) then
                  if(ikboun(id).eq.idof) cycle
               endif
!
               amax=dabs(a(i))
               imax=i
            endif
         enddo
!
         if(imax.eq.0) then
            write(*,*) '*ERROR in umpc_mean_rot'
            write(*,*) '       no mean rotation MPC can be'
            write(*,*) '       generated for the MPC containing'
            write(*,*) '       node ',jnode(1)
            call exit(201)
         endif
!
         nodeold=jnode(1)
         idirold=jdof(1)
         nodemax=jnode(imax)
         idirmax=jdof(imax)
!
!        exchange the node information, if needed
!
         if(nodemax.ne.nodeold) then
            do i=1,3
               iold(jdof(i))=i
            enddo
            do i=4,n-4
               if(jnode(i).eq.nodemax) then
                  if(jdof(i).eq.1) then
                     inew(1)=i
                  elseif(jdof(i).eq.2) then
                     inew(2)=i
                  else
                     inew(3)=i
                  endif
               endif
            enddo
!
            do j=1,3
               node=jnode(iold(j))
               idir=jdof(iold(j))
               xcoef=a(iold(j))
               jnode(iold(j))=jnode(inew(j))
               jdof(iold(j))=jdof(inew(j))
               a(iold(j))=a(inew(j))
               jnode(inew(j))=node
               jdof(inew(j))=idir
               a(inew(j))=xcoef
            enddo
         endif
!
!        exchange the direction information, if needed
!
         iold(1)=1
         do i=1,3
            if(jdof(i).eq.idirmax) then
               inew(1)=i
               exit
            endif
         enddo
!
         if(iold(1).ne.inew(1)) then
            do j=1,1
               idir=jdof(iold(j))
               xcoef=a(iold(j))
               jdof(iold(j))=jdof(inew(j))
               a(iold(j))=a(inew(j))
               jdof(inew(j))=idir
               a(inew(j))=xcoef
            enddo
         endif
      endif
!
      return
      end
