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
      subroutine frdnet(conet,nknet,konnet,nenet,vnet,time)
!     
!     stores the results in frd format
!     
      implicit none
!     
      character*1 c
      character*3 m1,m2,m3,m4,m5
      character*5 p0,p1,p2,p3,p4,p5,p6,p8,p10,p11,p12
      character*8 date,newclock,fmat
      character*10 clock
      character*20 newdate
      character*80 matname
      character*132 text
!     
      integer konnet(*),nknet,nenet,kode,i,j,indexe,one,null
!     
      real*8 conet(3,*),vnet(3,*),time
!     
      kode=1
!     
      c='C'
!     
      m1=' -1'
      m2=' -2'
      m3=' -3'
      m4=' -4'
      m5=' -5'
!     
      p0='    0'
      p1='    1'
      p2='    2'
      p3='    3'
      p4='    4'
      p5='    5'
      p6='    6'
      p8='    8'
      p10='   10'
      p11='   11'
      p12='   12'
!     
      if(time.le.0.d0) then
        fmat(1:8)='(e12.5) '
      elseif((dlog10(time).ge.0.d0).and.(dlog10(time).lt.11.d0)) then
        fmat(1:5)='(f12.'
        write(fmat(6:7),'(i2)') 11-int(dlog10(time)+1.d0)
        fmat(8:8)=')'
      else
        fmat(1:8)='(e12.5) '
      endif
!     
      null=0
      one=1
!     
      if(kode.eq.1) then
!     
        write(20,'(a5,a1)') p1,c
        call date_and_time(date,clock)
        newdate(1:20)='                    '
        newdate(1:2)=date(7:8)
        newdate(3:3)='.'
        if(date(5:6).eq.'01') then
          newdate(4:11)='january.'
          newdate(12:15)=date(1:4)
        elseif(date(5:6).eq.'02') then
          newdate(4:12)='february.'
          newdate(13:16)=date(1:4)
        elseif(date(5:6).eq.'03') then
          newdate(4:9)='march.'
          newdate(10:13)=date(1:4)
        elseif(date(5:6).eq.'04') then
          newdate(4:9)='april.'
          newdate(10:13)=date(1:4)
        elseif(date(5:6).eq.'05') then
          newdate(4:7)='may.'
          newdate(8:11)=date(1:4)
        elseif(date(5:6).eq.'06') then
          newdate(4:8)='june.'
          newdate(9:12)=date(1:4)
        elseif(date(5:6).eq.'07') then
          newdate(4:8)='july.'
          newdate(9:12)=date(1:4)
        elseif(date(5:6).eq.'08') then
          newdate(4:10)='august.'
          newdate(11:14)=date(1:4)
        elseif(date(5:6).eq.'09') then
          newdate(4:13)='september.'
          newdate(14:17)=date(1:4)
        elseif(date(5:6).eq.'10') then
          newdate(4:11)='october.'
          newdate(12:15)=date(1:4)
        elseif(date(5:6).eq.'11') then
          newdate(4:12)='november.'
          newdate(13:16)=date(1:4)
        elseif(date(5:6).eq.'12') then
          newdate(4:12)='december.'
          newdate(13:16)=date(1:4)
        endif
        newclock(1:2)=clock(1:2)
        newclock(3:3)=':'
        newclock(4:5)=clock(3:4)
        newclock(6:6)=':'
        newclock(7:8)=clock(5:6)
        write(20,'(a5,''UUSER'')') p1
        write(20,'(a5,''UDATE'',14x,a20)') p1,newdate
        write(20,'(a5,''UTIME'',14x,a8)') p1,newclock
        write(20,'(a5,''UHOST'')') p1
        write(20,'(a5,''UPGM               CalculiX'')') p1
        write(20,'(a5,''UDIR'')') p1
        write(20,'(a5,''UDBN'')') p1
!     
!     storing the coordinates of the nodes
!     
        write(20,'(a5,a1,67x,i1)') p2,c,one
!     
        do i=1,nknet
          write(20,100) m1,i,(conet(j,i),j=1,3)
        enddo
!     
        write(20,'(a3)') m3
!     
!     storing the element topology
!     
        write(20,'(a5,a1,67x,i1)') p3,c,one
!
        do i=1,5
          matname(i:i)=' '
        enddo
        do i=1,nenet
          indexe=8*(i-1)
          write(20,'(a3,i10,3a5)') m1,i,p1,p0,
     &         matname(1:5)
          write(20,'(a3,8i10)') m2,(konnet(indexe+j),j=1,8)
        enddo
!     
        write(20,'(a3)') m3
!     
      endif
!     
!     storing the fluid depth
!     
      do i=37,132
        text(i:i)=' '
      enddo
      text='    1PSTEP'
      write(text(25:36),'(i12)') kode
      write(20,'(a132)') text
!     
      text=
     & '  100CL       .00000E+00                                 3    1'
      text(75:75)='1'
      write(text(25:36),'(i12)') nknet
      write(text(8:12),'(i5)') 100+kode
      write(text(13:24),fmat) time
      write(text(59:63),'(i5)') kode
      write(20,'(a132)') text
      text=' -4  H           1    1'
      write(20,'(a132)') text
      text=' -5  H           1    1    0    0'
      write(20,'(a132)') text
!     
      do i=1,nknet
        write(20,100) m1,i,vnet(1,i)
      enddo
!     
      write(20,'(a3)') m3
!     
!     storing the Froude number
!     
      do i=37,132
        text(i:i)=' '
      enddo
      text='    1PSTEP'
      write(text(25:36),'(i12)') kode
      write(20,'(a132)') text
!     
      text=
     & '  100CL       .00000E+00                                 3    1'
      text(75:75)='1'
      write(text(25:36),'(i12)') nknet
      write(text(8:12),'(i5)') 100+kode
      write(text(13:24),fmat) time
      write(text(59:63),'(i5)') kode
      write(20,'(a132)') text
      text=' -4  Fr          1    1'
      write(20,'(a132)') text
      text=' -5  Fr          1    1    0    0'
      write(20,'(a132)') text
!     
      do i=1,nknet
        write(20,100) m1,i,vnet(2,i)
      enddo
!     
      write(20,'(a3)') m3
!     
!     storing the total head
!     
      do i=37,132
        text(i:i)=' '
      enddo
      text='    1PSTEP'
      write(text(25:36),'(i12)') kode
      write(20,'(a132)') text
!     
      text=
     & '  100CL       .00000E+00                                 3    1'
      text(75:75)='1'
      write(text(25:36),'(i12)') nknet
      write(text(8:12),'(i5)') 100+kode
      write(text(13:24),fmat) time
      write(text(59:63),'(i5)') kode
      write(20,'(a132)') text
      text=' -4  THEAD       1    1'
      write(20,'(a132)') text
      text=' -5  THEAD       1    1    0    0'
      write(20,'(a132)') text
!     
      do i=1,nknet
        write(20,100) m1,i,vnet(3,i)
      enddo
!     
      write(20,'(a3)') m3
!     
 100  format(a3,i10,1p,6e12.5)
!     
      return
      end
