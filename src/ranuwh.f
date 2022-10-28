      REAL FUNCTION RANUWH ()
C
C Portable uniform random number generator, using the
C method of Wichmann and Hill. This file contains the program units
C RANUWH, INIRAN and RAWHIN. These are auxiliary functions for FMINSI,
C but may also be used separately.
C
C   FMINSI - Fortran subroutines for unconstrained function minimization
C   Copyright (C) 1992, 1995, 2001  Hugo Pfoertner
C
C   This library is free software; you can redistribute it and/or
C   modify it under the terms of the GNU Lesser General Public
C   License as published by the Free Software Foundation; either
C   version 2.1 of the License, or (at your option) any later version.
C
C   This library is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C   Lesser General Public License for more details.
C
C   You should have received a copy of the GNU Lesser General Public
C   License along with this library; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C   USA
C
C   Contact info: mailto:hugo@pfoertner.org
C   or use the information provided at http://www.pfoertner.org/
C
C Author: Hugo Pfoertner, Oberhaching, Germany
C 
C Uniform random numbers R in the range 0.0 <= R < 1.0 are created
C This function is written in portable Fortran 77 and does not depend
C on the implementation of overflow of integers. It should return
C the same sequence of random numbers on any machine, independent of
C the machine's architecture. Only standard integer arithmetic is used.
C
C
C Version History (left in German):
C 02.06.01 English translation of comments, LGPL header added
C 26.08.95 EXTERNAL RAEWIN ENGEFUEGT, UM UEBER BLOCKDATA-LINK
C          STARTBELEGUNG AUCH OHNE INIRAN-AUFRUF ZU ERZWINGEN
C 07.12.92 BASISVERSION
C
C Reference: WICHMANN AND HILL: APPL. STATIST. (JRSSC),
C                               (31) 188-190, (1982)
C
C Memory, implemented in Common block RWHRAN
C Before the first call the memory of the function has to be
C preset by calling INIRAN
C
      INTEGER IX, IY, IZ
      COMMON /RWHRAN/ IX, IY, IZ
C
C The block data program unit is referenced here to assure a
C proper initialization of the memory, even if the call to
C INIRAN is forgotten. The EXTERNAL statement forces the linker
C to include the initial values.
C
      EXTERNAL RAWHIN
C
C Advance memory
      IX = 171 * MOD ( IX, 177) -  2 * ( IX / 177 )
      IY = 172 * MOD ( IY, 176) - 35 * ( IY / 176 )
      IZ = 170 * MOD ( IZ, 178) - 63 * ( IZ / 178 )
C
C Bring into non-negative range
      IF ( IX .LT. 0 ) IX = IX + 30269
      IF ( IY .LT. 0 ) IY = IY + 30307
      IF ( IZ .LT. 0 ) IZ = IZ + 30323
C
C Truncate to 0..1 real number
      RANUWH = MOD ( REAL(IX) / 30269.0
     &             + REAL(IY) / 30307.0
     &             + REAL(IZ) / 30323.0,  1.0 )
C
      RETURN
C End of function RANUWH
      END
C *******************************************************************
      SUBROUTINE INIRAN
C
C This subroutine is part of the FMINSI program library.
C Purpose:
C Set starting values for uniform random number generator RANUWH
C See copyright notice given in function RANUWH
C
C Author: Hugo Pfoertner, Oberhaching
C
C Version History:
C 07.12.92 Initial version
C
C Reference: WICHMANN AND HILL: APPL. STATIST. (JRSSC),
C                               (31) 188-190, (1982)
C
C Memory:
      INTEGER IX, IY, IZ
      COMMON /RWHRAN/ IX, IY, IZ
C
C Set initial values
      IX = 1974
      IY = 235
      IZ = 337
C
      RETURN
C End of subroutine INIRAN
      END
C *******************************************************************
      BLOCKDATA RAWHIN
C
C This program unit is part of the FMINSI program library
C Purpose:
C Forced setting of initial values (e.g. if call to INIRAN is omitted)
C See copyright notice given in function RANUWH
C
C Author: Hugo Pfoertner, Oberhaching
C
C Version History:
C 26.08.95 Initial version
C
      INTEGER IX, IY, IZ
      COMMON /RWHRAN/ IX, IY, IZ
      DATA IX, IY, IZ / 1974, 235, 337 /
C End of blockdata RAWHIN
      END
