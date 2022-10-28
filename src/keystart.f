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
      subroutine keystart(ifreeinp,ipoinp,inp,name,iline,ikey)
!
      implicit none
!
!     stores the order in which the input is to be read in fields
!     ipoinp and inp; before adding keywords to the list below
!     look at the CalculiX documentation, section "Program structure"
!     subsection "Allocation of the fields" subsubsection "readinput".
!
!     order:
!     1)  *RESTART,READ
!     2)  *NODE
!     3)  *USER ELEMENT
!     4)  *ELEMENT
!     5)  *NSET
!     6)  *ELSET
!     7)  *SURFACE
!     8)  *TRANSFORM
!     9)  *MATERIAL
!     10) *DISTRIBUTION
!     11) *ORIENTATION
!     12) *TIE
!     13) *SURFACE INTERACTION
!     14) *INITIAL CONDITIONS
!     15) *AMPLITUDE
!     16) *CONTACT PAIR
!     17) *COUPLING
!     18) everything else
!
      integer nentries
      parameter(nentries=18)
!
      character*20 name,nameref(nentries)
!
      integer ifreeinp,ipoinp(2,*),inp(3,*),namelen(nentries),i,ikey,
     &  iline
!
!     order in which the cards have to be read
!
      data nameref /'RESTART,READ','NODE','USERELEMENT','ELEMENT',
     &              'NSET',
     &              'ELSET','SURFACE','TRANSFORM','MATERIAL',
     &              'DISTRIBUTION',
     &              'ORIENTATION','TIE','INTERACTION',
     &              'INITIALCONDITIONS','AMPLITUDE',
     &              'CONTACTPAIR','COUPLING','REST'/
!
!     length of the names in field nameref
!
      data namelen /12,4,11,7,4,5,7,9,8,12,11,3,11,17,9,11,8,4/
!
      do i=1,nentries
         if(name(1:namelen(i)).eq.nameref(i)(1:namelen(i))) then
            if(ikey.eq.i) return
            if(ikey.gt.0) inp(2,ipoinp(2,ikey))=iline-1
            ikey=i
            if(ipoinp(1,i).eq.0) then
               ipoinp(1,i)=ifreeinp
            else
               inp(3,ipoinp(2,i))=ifreeinp
            endif
            ipoinp(2,i)=ifreeinp
            exit
         endif
      enddo
      inp(1,ifreeinp)=iline
      ifreeinp=ifreeinp+1
!
      return
      end
