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
!     S.W. Sloan, Adv.Eng.Software,1987,9(1),34-55.
!     Permission for use with the GPL license granted by Prof. Scott
!     Sloan on 17. Nov. 2013
!
      subroutine qsorti(n,list,key)
!     
      implicit none
!     
      integer list(*),key(*),n,ll,lr,lm,nl,nr,ltemp,stktop,maxstk,guess
!     
      parameter(maxstk=32)
!     
      integer lstack(maxstk),rstack(maxstk)
!     
      ll= 1
      lr=n
      stktop=0
 10   if(ll.lt.lr) then
         nl=ll
         nr=lr
         lm=(ll+lr)/2
         guess=key(list(lm))
 20      if (key(list(nl)).lt.guess) then
            nl=nl+1
            goto 20
         end if
 30      if (guess.lt.key(list(nr))) then
            nr=nr-1
            goto 30
         end if
         if(nl.lt.(nr-1)) then
            ltemp=list(nl)
            list(nl)=list(nr)
            list(nr)=ltemp
            nl=nl+1
            nr=nr-1
            goto 20
         end if
         if(nl.le.nr) then
            if(nl.lt.nr) then
               ltemp=list(nl)
               list(nl)=list(nr)
               list(nr)=ltemp
            end if
            nl=nl+1
            nr=nr-1
         end if
         stktop=stktop+1
         if(nr.lt.lm) then
            lstack(stktop)=nl
            rstack(stktop)=lr
            lr=nr
         else
            lstack(stktop)=ll
            rstack(stktop)=nr
            ll=nl
         end if
         goto 10
      end if
      if (stktop.ne.0) then
         ll=lstack(stktop)
         lr=rstack(stktop)
         stktop=stktop-1
         goto 10
      end if
      end
      
