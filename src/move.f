!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, qUSA.
!

C> move
C> ====
C>
C> Move a file on disk from source path to target path.
C>
C> For the GNU and Intel compilers, the move operation is delegated to
C> their RENAME compiler extensions. For other compilers, a standard
C> Fortran implementation move_ moves the file by copying it first
C> to target and then deleting the source.
C>
C> Parameters
C> ----------
C> @param source Filename or path of the source file.
C> @param target Filename of path of the move target.
C> @param iostat Error code.
C>
      subroutine move(source, target, iostat)
!
!     moves a file on disk
!
#if defined(__INTEL_COMPILER)
      use ifport, only: rename
#endif
!
      implicit none
!
      character*132 source, target
      integer iostat
!
#if defined(__GFORTRAN__) || defined(__INTEL_COMPILER)
!     use GNU extension or IFPORT library function
      iostat = rename(source, target)
#else
!     fall back to copying and deleting
      call move_(source, target, iostat)
#endif
!
      end
!

C> move_
C> =====
C>
C> Move a file on disk by copying source to target and deleting source.
C>
C> This serves as a fallback for compilers which do not provide an
C> extension for renaming a file.
C>
C> Notes
C> -----
C> The copying is done by opening source and target as unformatted files
C> with stream access and reading/writing the data via a fixed size
C> buffer array of char type.
C>
C> Parameters
C> ----------
C> @param source Filename or path of the source file.
C> @param target Filename of path of the move target.
C> @param iostat Error code.
C>
      subroutine move_(source, target, iostat)
!
!     copy raw bits of source into target, then remove source
!
      use iso_fortran_env, only: file_storage_size, iostat_end, int64
!
      implicit none
!
      character(len=132) source, target
      integer source_unit, target_unit, iostat
      integer(int64) source_size, buffer_size, buffer_unit, remaining,
     &     i, iend
      integer, parameter :: buffer_units = 1024
!
      character(len=1), dimension(buffer_units) :: buffer
!
      buffer_unit = storage_size(buffer)
      buffer_size = buffer_unit*size(buffer)
!
!     open source file for reading
!
      open(newunit=source_unit, file=source, access='stream',
     &     form='unformatted', status='old', action='read',
     &     iostat=iostat)
      if (iostat .ne. 0) then
         write(*,*) "Error opening source file: ", iostat
         return
      end if
!
!     open destination file for writing
!
      open(newunit=target_unit, file=target, access='stream',
     &     form='unformatted', status='new', action='write',
     &     iostat=iostat)
      if (iostat .ne. 0) then
         write(*,*) "Error opening target file: ", iostat
         return
      end if
!
!     determine source file size in bits
!
      inquire(source_unit, size=source_size)
      source_size = file_storage_size*source_size
!
!     read from source file in buffer and write buffer to target file
!
      i = 1
      do while (i*buffer_size <= source_size)
         read(source_unit, iostat=iostat) buffer
         if (iostat .ne. 0) return
         write(target_unit) buffer
         i = i+1
      end do
      remaining = source_size - (i-1)*buffer_size
      iend = (remaining+buffer_unit-1)/buffer_unit
!
      read(source_unit, iostat=iostat) buffer
      if (iostat .ne. iostat_end) then
         return
      else
         iostat = 0
      end if
      write(target_unit) buffer(1:iend)
!
!     close both files and delete the source file
!
      close(source_unit, status='delete')
      close(target_unit)
!
      end