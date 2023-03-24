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
      subroutine propertynet(ieg,nflow,prop,ielprop,lakon,iin,
     &       prop_store,ttime,time,nam,amname,namta,amta)
!
!     user subroutine propertynet
!
!
!     INPUT:
!
!     ieg(i)             global element number corresponding to
!                        network element i (i=1,...,nflow)
!     nflow              number of network elements
!     ielprop(i)         property to the position in fields prop and
!                        prop_store after which the properties for 
!                        element i start (prop(ielprop(i)+1),
!                        prop(ielprop(i)+2).....). The number is dictated
!                        by the type of element.
!     lakon(i)           label of element i
!     iin                gas network iteration number
!     prop_store         property values as specified in the
!                        input deck
!     ttime              total time
!     time               step time
!     nam                number of amplitudes
!     amname(i)          amplitude name of amplitude i
!     namta(1,i)         location of first (time,amplitude) pair in
!                        field amta
!     namta(2,i)         location of last (time,amplitude) pair in
!                        field amta
!     namta(3,i)         in absolute value the amplitude it refers to; if
!                        abs(namta(3,i))=i it refers to itself. If
!                        abs(namta(3,i))=j, amplitude i is a time delay
!                        of amplitude j; in that case the value of the
!                        time delay is stored in
!                        amta(1,namta(1,i)); in the latter case
!                        amta(2,namta(1,i)) is without meaning; if
!                        namta(3,i)>0 the time in amta for amplitude i is
!                        step time, else it is total time.
!     amta(1,i)          time of (time,amplitude)-pair i
!     amta(2,i)          amplitude of (time,amplitude)-pair i
!
!     OUTPUT:
!
!     prop               actual property values
!           
      implicit none
!
      character*8 lakon(*)
      character*80 amname(*)
!
      integer ieg(*),nflow,ielprop(*),iin,nam,namta(3,*)
!
      real*8 prop(*),prop_store(*),ttime,time,amta(2,*)
!
!
!
      return
      end

