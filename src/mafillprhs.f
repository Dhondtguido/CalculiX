!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine mafillprhs(nk,kon,ipkon,lakon,ipompc,nodempc,
     &     coefmpc,nmpc,b1,nactdoh,mi,v,theta1,ne,dtimef,ipvar,var,
     &     compressible,num_cpus)
!     
!     filling the rhs b of the pressure equations (step 2)
!
      use omp_lib
!
      implicit none
!     
      character*8 lakon(*)
!     
      integer kon(*),ipompc(*),nodempc(3,*),nactdoh(*),compressible,
     &     mi(*),ipkon(*),ne,ipvar(*),nk,nmpc,i,j,
     &     id,ist,index,jdof1,node,indexe,nope,num_cpus,tid
!     
      real*8 coefmpc(*),b1(nk,0:mi(2)),v(nk,0:mi(2)),
     &     ff(8),theta1,var(*),dtimef
!
!     We use heap allocated b1_ where each thread owns a slice
!     instead of an OpenMP array reduction clause into b1 which
!     might exceed the default thread stack sizes for very large
!     models.
!
      real*8, allocatable :: b1_(:,:)
!
      allocate(b1_(nk,num_cpus))
!
      b1_=0
!
!$omp parallel private(tid) num_threads(num_cpus)
      tid = omp_get_thread_num() + 1
!$omp do private(j,index,indexe,nope,node,jdof1,id,ist,ff)
      do i=1,ne
!     
        indexe=ipkon(i)
        if(lakon(i)(4:4).eq.'8') then
          nope=8
        elseif(lakon(i)(4:4).eq.'4') then
          nope=4
        elseif(lakon(i)(4:4).eq.'6') then
          nope=6
        else
          cycle
        endif
!     
        call e_c3d_prhs(nk,kon(indexe+1),lakon(i),ff,i,v,dtimef,mi(1),
     &       theta1,ipvar,var)
!     
        if(compressible.eq.1) then
!
!         compressible fluid: dof = node number
!
          do j=1,nope
            node=kon(indexe+j)
            b1_(node,tid)=b1_(node,tid)+ff(j)
          enddo
        else
!
!     compressible fluid: dofs may have been removed due to
!     SPC's and/or MPC's
!
          do j=1,nope
!     
            node=kon(indexe+j)
            jdof1=nactdoh(node)
!     
!     inclusion of ff
!     
            if(jdof1.le.0) then
              if(nmpc.ne.0) then
                if(jdof1.ne.2*(jdof1/2)) then
                  id=(-jdof1+1)/2
                  ist=ipompc(id)
                  index=nodempc(3,ist)
                  if(index.eq.0) cycle
                  do
                    jdof1=nactdoh(nodempc(1,index))
                    if(jdof1.gt.0) then
                      b1_(jdof1,tid)=b1_(jdof1,tid)
     &                     -coefmpc(index)*ff(j)
     &                     /coefmpc(ist)
                    endif
                    index=nodempc(3,index)
                    if(index.eq.0) exit
                  enddo
                endif
              endif
              cycle
            endif
            b1_(jdof1,tid)=b1_(jdof1,tid)+ff(j)
!     
          enddo
        endif
      enddo
!$omp end do
!
      do j = 1, num_cpus
!$omp do
         do i=1,nk
            b1(i,4)=b1(i,4)+b1_(i,j)
         end do
!$omp end do
      end do
!$omp end parallel
!
      deallocate(b1_)
!
c      write(*,*) 'mafillprhs '
c      do i=1,nk
c        write(*,*) i,b1(i,4)
c      enddo
!     
      return
      end
      
