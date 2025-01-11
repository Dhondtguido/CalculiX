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
      subroutine identdesifaces(iregion,nsurfs,ipkonfa,lakonfa,konfa,
     &   ndesifaces,idesiface,nodedesiinv)
!
      implicit none
!
      character*8 lakonfa(*)
!
      integer i,m,nopedesi,index1,ipkonfa(*),nopem,nodeface(8),nnodes,
     &   ndesifaces,idesiface(*),nodedesiinv(*),konfa(*),nsurfs,
     &   iregion
!
!
!     lists which external faces have enough design variables on it
!     ndesifaces       number of surfaces containing design variables
!     idesiface(i)   pointer to surface in ipkonfa(*)
!
      ndesifaces=0
      do i=1,nsurfs
         if(lakonfa(i)(1:2).eq.'S3') then
            nopedesi=3
         elseif(lakonfa(i)(1:2).eq.'S4') then
            nopedesi=3
         elseif(lakonfa(i)(1:2).eq.'S6') then
            nopedesi=4
         elseif(lakonfa(i)(1:2).eq.'S8') then
            nopedesi=5
         else
            exit
         endif
         if(iregion.eq.0) nopedesi=1
!
         index1=ipkonfa(i)
         nopem=ipkonfa(i+1)-ipkonfa(i)
         do m=1,nopem
            nodeface(m)=konfa(index1+m)
         enddo
!
         nnodes=0
         do m=1,nopem
            if(nodedesiinv(nodeface(m)).ne.0) then
               nnodes=nnodes+1
            endif
         enddo
!     
         if(nnodes.ge.nopedesi) then
            ndesifaces=ndesifaces+1
            idesiface(ndesifaces)=i
         endif
      enddo
!
      return
      end
