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
      subroutine plmix(plcon,nplcon,plconloc,npmat_,ntmat_,
     &  imat,j,temp,nelem,kin)
!
!     interpolates the hardening data for material imat and temperature
!     j and j-1 to obtain data for temperature temp. The data is taken
!     from plcon and stored in plconloc.
!     The Von Mises stress is interpolated for a given equivalent
!     plastic strain. If the equivalent strain data points for
!     temperature j and j-1 do not coincide, the union of both is
!     taken. If this union exceeds 200 (ierror=1), the equivalent plastic
!     strain range is divided into 199 intervals yielding 200 new
!     equivalent strain data points, for which the Von Mises stress
!     is interpolated.
!     Attention: in plcon the odd storage spaces contain the Von
!                Mises stress, the even ones the equivalent plastic
!                strain. For plconloc, this order is reversed.
!
      implicit none
!
      integer imat,ndata,ntmat_,npmat_,nplcon(0:ntmat_,*),nelem,
     &  kin,k,j,k1,k2,ierror,ndata1,ndata2,itemp
!
      real*8 eplmin,eplmax,depl,epla,plcon(0:2*npmat_,ntmat_,*),
     &  plconloc(802),dummy,temp,ep1,ep2,t1,t2,s1,s2,ratio
!
!
!
      ndata=0
      ierror=0
!
      ndata1=nplcon(j-1,imat)
      ndata2=nplcon(j,imat)
      t1=plcon(0,j-1,imat)
      t2=plcon(0,j,imat)
      ratio=(temp-t1)/(t2-t1)
!
!     the interval on which the stress interpolation is performed
!     is the intersection of the domain of the two curves
!
      k1=1
      k2=1
      ep1=plcon(2,j-1,imat)
      ep2=plcon(2,j,imat)
      if(ep1.gt.ep2) then
         do k2=1,ndata2
            ep2=plcon(2*k2,j,imat)
            if(ep2.gt.ep1) exit
         enddo
         if(k2.gt.ndata2) then
            write(*,*) '*ERROR in plmix: there exist two temperatures'
            write(*,*) '       for which the hardening curves are'
            write(*,*) '       disjunct'
            call exit(201)
         endif
      elseif(ep2.gt.ep1) then
         do k1=1,ndata1
            ep1=plcon(2*k1,j-1,imat)
            if(ep1.gt.ep2) exit
         enddo
         if(k1.gt.ndata1) then
            write(*,*) '*ERROR in plmix: there exist two temperatures'
            write(*,*) '       for which the hardening curves are'
            write(*,*) '       disjunct'
            call exit(201)
         endif
      endif
!
      do
         s1=plcon(2*k1-1,j-1,imat)
         s2=plcon(2*k2-1,j,imat)
         ep1=plcon(2*k1,j-1,imat)
         ep2=plcon(2*k2,j,imat)
!
         if(dabs(ep1-ep2).lt.1.d-10) then
            if(k2.lt.ndata2) then
               k2=k2+1
            elseif(k1.lt.ndata1) then
               k1=k1+1
            else
               ndata=ndata+1
               if(ndata.gt.200) then
                  ierror=1
                  exit
               endif
               if(kin.eq.0) then
                  plconloc(2*ndata-1)=ep1+ratio*(ep2-ep1)
                  plconloc(2*ndata)=s1+ratio*(s2-s1)
               else
                  plconloc(399+2*ndata)=ep1+ratio*(ep2-ep1)
                  plconloc(400+2*ndata)=s1+ratio*(s2-s1)
               endif
               exit
            endif
            cycle
         endif
         if(ep1.lt.ep2) then
            ndata=ndata+1
            if(ndata.gt.200) then
               ierror=1
               exit
            endif
            call plinterpol(plcon,nplcon,j,s2,dummy,npmat_,ntmat_,
     &        imat,nelem,ep1)
            if(kin.eq.0) then
               plconloc(2*ndata-1)=ep1
               plconloc(2*ndata)=s1+ratio*(s2-s1)
            else
               plconloc(399+2*ndata)=ep1
               plconloc(400+2*ndata)=s1+ratio*(s2-s1)
            endif
            if(k1.lt.ndata1) then
               k1=k1+1
               cycle
            else
               exit
            endif
         else
            ndata=ndata+1
            if(ndata.gt.200) then
               ierror=1
               exit
            endif
            call plinterpol(plcon,nplcon,j-1,s1,dummy,npmat_,ntmat_,
     &        imat,nelem,ep2)
            if(kin.eq.0) then
               plconloc(2*ndata-1)=ep2
               plconloc(2*ndata)=s1+ratio*(s2-s1)
            else
               plconloc(399+2*ndata)=ep2
               plconloc(400+2*ndata)=s1+ratio*(s2-s1)
            endif
            if(k2.lt.ndata2) then
               k2=k2+1
               cycle
            else
               exit
            endif
         endif
      enddo
!
!     if more than 200 data points result, the interval is divided into
!     199 equidistant intervals
!
      if(ierror.eq.0) then
         if(kin.eq.0) then
            plconloc(801)=real(ndata)+0.5d0
         else
            plconloc(802)=real(ndata)+0.5d0
         endif
      else
         if(kin.eq.0) then
            eplmin=max(plcon(2,j-1,imat),plcon(2,j,imat))
            eplmax=min(plcon(2*ndata1,j-1,imat),plcon(2*ndata2,j,imat))
     &         -1.d-10
            depl=(eplmax-eplmin)/199.d0
            do k=1,200
               epla=eplmin+(k-1)*depl
               itemp=j-1
               call plinterpol(plcon,nplcon,itemp,s1,
     &              dummy,npmat_,ntmat_,imat,nelem,epla)
               itemp=j
               call plinterpol(plcon,nplcon,itemp,s2,
     &              dummy,npmat_,ntmat_,imat,nelem,epla)
               plconloc(2*k-1)=epla
               plconloc(2*k)=s1+ratio*(s2-s1)
            enddo
            plconloc(801)=200.5d0
         else
            eplmin=max(plcon(2,j-1,imat),plcon(2,j,imat))
            eplmax=min(plcon(2*ndata1,j-1,imat),plcon(2*ndata2,j,imat))
     &         -1.d-10
            depl=(eplmax-eplmin)/199.d0
            do k=1,200
               epla=eplmin+(k-1)*depl
               itemp=j-1
               call plinterpol(plcon,nplcon,itemp,s1,
     &              dummy,npmat_,ntmat_,imat,nelem,epla)
               itemp=j
               call plinterpol(plcon,nplcon,itemp,s2,
     &              dummy,npmat_,ntmat_,imat,nelem,epla)
               plconloc(399+2*k)=epla
               plconloc(400+2*k)=s1+ratio*(s2-s1)
            enddo
            plconloc(802)=200.5d0
         endif
      endif
!
      return
      end
