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
      subroutine createblock_struct(neq,ipointers,icolpardiso,aupardiso,
     &  nestart,num_cpus,ja,nz_num)
!
!     generates num_cpu blocks
!
      integer neq,ipointers(*),icolpardiso(*),nestart(*),num_cpus,ja(*),
     &  numd,i,j,k,m,icol,isubtract,nz_num
!
      real*8 aupardiso(*)
! 
!     nestart(i) points to the element before the block for which
!     cpu i is responsible
!
      numd=int(neq/num_cpus)+1
      nestart(1)=0
      do i=2,num_cpus
         nestart(i)=nestart(i-1)+numd
      enddo
      nestart(num_cpus+1)=neq
!
!     ipointers(i) points to the first entry of row i in icolpardiso
!     ja(i) points to the entry in icolpardiso before the start of row i
!
      j=0
!
      do k=1,num_cpus
         do i=nestart(k)+1,nestart(k+1)
            ja(i)=j
            do m=ipointers(i),ipointers(i+1)-1
               icol=icolpardiso(m)
               if((icol.gt.nestart(k)).and.
     &            (icol.le.nestart(k+1))) then
                  j=j+1
                  icolpardiso(j)=icol
                  aupardiso(j)=aupardiso(m)
               endif
             enddo
          enddo
       enddo
       ja(neq+1)=j
       nz_num=j
!
!     subtracting from iam the number of elements belonging
!     to the preceding blocks
!
      do k=2,num_cpus
         isubtract=nestart(k)
         do i=nestart(k)+1,nestart(k+1)
            do j=ja(i)+1,ja(i+1)
               icolpardiso(j)=icolpardiso(j)-isubtract
            enddo
         enddo
      enddo
!
      return
      end
