!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine writeinputdeck2(feasdir,nodedesi,ndesi,
     &     inoel,iponoel,xdesi,co,lakon,ipkon,kon,tinc,nk)
!     
!     storing the feasible directionvalues to the mesh modification file
!
      character*8 lakon(*)
!
      integer nodedesi(*),ndesi,inoel(2,*),iponoel(*),ipkon(*),
     &     kon(*),i,j,k,m,nsort(3),two,three,four
!
      real*8 feasdir(2,*),xdesi(3,*),co(3,*),tinc,xnorloc(3),sort(3),dd
!
      two=2
      three=3
      four=4
!
      write(20,107)
      do m=1,ndesi
        write(20,108) nk+m
      enddo
 107  format('*NODE')
 108  format(i10,',0.,0.,0.')
!
      do m=1,ndesi
        i=nodedesi(m)
!     
!     consideration of plain stress/strain 2d-elements
!     and axisymmetric elements        
!     
        ielem=inoel(1,iponoel(i))
        if((lakon(ielem)(7:7).eq.'A').or.
     &       (lakon(ielem)(7:7).eq.'S').or.
     &       (lakon(ielem)(7:7).eq.'E')) then
!     
          if(lakon(ielem)(4:5).eq.'20') then
            nope=20
          elseif (lakon(ielem)(4:4).eq.'8') then
            nope=8
          elseif (lakon(ielem)(4:5).eq.'15') then
            nope=15
          elseif (lakon(ielem)(4:4).eq.'6') then
            nope=6
          else
            cycle
          endif
!     
          indexe=ipkon(ielem)
          do j=1,nope
            if(i.eq.kon(indexe+j)) then
              exit
            endif
          enddo
!     
!     replace 3D node number by 2D node number     
!     the only 3D nodes to be replaced are the expanded node with   
!     the lowest node number (of the three expanded nodes; so nodes
!     with local numbers 1,2,3,4,9,10,11,12 for a C3D20(R) element). 
!     The other two are skipped by the statement nodedesiinv(.)=-1  
!     
          if(lakon(ielem)(4:5).eq.'20') then
            if(j.gt.4) j=j-4
            node=kon(indexe+nope+j)
          elseif(lakon(ielem)(4:5).eq.'15') then
            if(j.gt.3) j=j-3
            node=kon(indexe+nope+j)
          elseif(lakon(ielem)(4:4).eq.'8') then
            node=kon(indexe+nope+j)
          elseif(lakon(ielem)(4:5).eq.'6') then
            node=kon(indexe+nope+j)
          endif
!     
        elseif(lakon(ielem)(7:7).eq.'L') then
!     
!     no output for shell elements necessary
!     
          cycle
        else
!     
!     in case of a 3D model no change of node number     
          node=i
        endif
!
        dd=dsqrt(xdesi(1,m)**2+xdesi(2,m)**2+xdesi(3,m)**2)
        do k=1,3
          xnorloc(k)=xdesi(k,m)/dd
          sort(k)=dabs(xdesi(k,m)/dd)
          nsort(k)=k
        enddo
        call dsort(sort,nsort,three,two)
!
!       if the two largest values are equal, the entry with the
!       smallest index is taken as dependent entry
!
        if(dabs(sort(2)-sort(3)).lt.1.d-10) then
          if(nsort(2).lt.nsort(3)) then
            idummy=nsort(2)
            nsort(2)=nsort(3)
            nsort(3)=idummy
          endif
        endif
!
        write(20,101) four
        write(20,102) node,nsort(3),xnorloc(nsort(3)),
     &       node,nsort(2),xnorloc(nsort(2)),
     &       node,nsort(1),xnorloc(nsort(1)),
     &       nk+m
        write(20,103) nk+m,feasdir(2,i)*tinc
      enddo
      write(20,105)
!      
 101  format('*EQUATION',/,i1)
 102  format(3(i10,",",i1,",",e20.13,","),i10,',1,-1.')
 103  format('*BOUNDARY',/,i10,',1,1,',e20.13)
 105  format('*STEP',/,'*STATIC',/,'*NODE FILE',/,'U',/,'*END STEP')
      close(20)
!     
      return        
      end
