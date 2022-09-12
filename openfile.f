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
      subroutine openfile(jobname)
!
      implicit none
!
      character*132 jobname,fndat,fnfrd,fnsta,fncvg,fncel,fnfcv
      integer i
!
!     opening the input and output file
!
      do i=1,132
         if(jobname(i:i).eq.' ') exit
      enddo
      i=i-1
      if(i.gt.128) then
         write(*,*) '*ERROR in openfile: input file name is too long:'
         write(*,'(a132)') jobname(1:132)
         write(*,*) '       exceeds 128 characters'
         call exit(201)
      endif
!
      fndat=jobname(1:i)//'.dat'
      open(5,file=fndat(1:i+4),status='unknown',err=51)
      close(5,status='delete',err=52)
      open(5,file=fndat(1:i+4),status='unknown',err=51)
!
!     delete the .frd file (it is reopened in frd.c)
!
      fnfrd=jobname(1:i)//'.frd'
      open(7,file=fnfrd(1:i+4),status='unknown',err=71)
      close(7,status='delete',err=72)
!
!     delete the .fcv file (it is reopened in compfluid.c)
!
      fnfcv=jobname(1:i)//'.fcv'
      open(12,file=fnfcv(1:i+4),status='unknown',err=74)
      close(12,status='delete',err=75)
!
!     .sta-file
!
      fnsta=jobname(1:i)//'.sta'
      open(8,file=fnsta(1:i+4),status='unknown',err=81)
      close(8,status='delete',err=82)
      open(8,file=fnsta(1:i+4),status='unknown',err=81)
      write(8,100)
      write(8,101)
 100  format('SUMMARY OF JOB INFORMATION')
 101  format('  STEP      INC     ATT  ITRS     TOT TIME     STEP TIME      
     &    INC TIME')
!
!     .cvg-file
!
      fncvg=jobname(1:i)//'.cvg'
      open(11,file=fncvg(1:i+4),status='unknown',err=91)
      close(11,status='delete',err=92)
      open(11,file=fncvg(1:i+4),status='unknown',err=91)
      write(11,102)
      write(11,103)
      write(11,104)
      write(11,105)
 102  format('SUMMARY OF C0NVERGENCE INFORMATION')
 103  format('  STEP   INC  ATT   ITER     CONT.   RESID.        CORR.   
     &    RESID.      CORR.'     )
 104  format('                              EL.    FORCE         DISP    
     &    FLUX        TEMP.'   )
 105  format('                              (#)     (%)           (%)    
     &     (%)         (%)')
!
!     contact elements
!
      fncel=jobname(1:i)//'.cel'
      open(27,file=fncel(1:i+4),status='unknown',err=93)
      close(27,status='delete',err=94)
!
      return
!
 51   write(*,*) '*ERROR in openfile: could not open file ',fndat(1:i+4)
      call exit(201)
 52   write(*,*) '*ERROR in openfile: could not delete file ',
     &  fndat(1:i+4)
      call exit(201)
 71   write(*,*) '*ERROR in openfile: could not open file ',fnfrd(1:i+4)
      call exit(201)
 72   write(*,*) '*ERROR in openfile: could not delete file ',
     &  fnfrd(1:i+4)
      call exit(201)
 74   write(*,*) '*ERROR in openfile: could not open file ',fnfcv(1:i+4)
      call exit(201)
 75   write(*,*) '*ERROR in openfile: could not delete file ',
     &  fnfcv(1:i+5)
      call exit(201)
 81   write(*,*) '*ERROR in openfile: could not open file ',fnsta(1:i+4)
      call exit(201)
 82   write(*,*) '*ERROR in openfile: could not delete file ',
     &  fnsta(1:i+4)
 91   write(*,*) '*ERROR in openfile: could not open file ',fncvg(1:i+4)
      call exit(201)
 92   write(*,*) '*ERROR in openfile: could not delete file ',
     &  fncvg(1:i+4)
 93   write(*,*) '*ERROR in openfile: could not open file ',fncel(1:19)
      call exit(201)
 94   write(*,*) '*ERROR in openfile: could not delete file ',
     &  fncel(1:19)
      call exit(201)
      end
