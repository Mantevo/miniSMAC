c copyright notice
c                      Copyright 2013 Sandia Corporation             
c                                                                    
c     Under the terms of Contract DE-AC04-94AL85000 with Sandia      
c     Corporation, the U.S. Government retains certain rights in this
c     software.                                                      
c                                                                    
c     Redistribution and use in source and binary forms, with or     
c     without modification, are permitted provided that the following
c     conditions are met:                                            
c                                                                    
c      1. Redistributions of source code must retain the above       
c         copyright notice, this list of conditions and the following
c         disclaimer.                                                
c                                                                    
c      2. Redistributions in binary form must reproduce the above    
c         copyright notice, this list of conditions and the following
c         disclaimer in the documentation and/or other materials     
c         provided with the distribution.                            
c                                                                    
c      3. Neither the name of the Corporation nor the names of the   
c         contributors may be used to endorse or promote products    
c         derived from this software without specific prior written  
c         permission.                                                
c                                                                    
c     THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
c     EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,  
c     THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A    
c     PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA    
c     CORPORATION OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,      
c     INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL     
c     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         
c     SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;   
c     OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
c     LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      
c     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF  
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.                                                   
c     
c
c************************************************************
c Subroutines in this file:
c  1. ilusol3
c  2. ilusol1 
c
c************************************************************
c
c
c************************************************************************
      subroutine ilusol3(jmax,kmax,dq,rhs,btc,bjm,bjp,bkm,bkp)
c************************************************************************
c  Scalar version of the routine.
c
c  Approximate solve of system M*dq = rhs, where M is
c  ilu decomposition of lhs of original system.
c
c  This routine assumes that btc contains the modified block array 'bprm'
c    computed by iluini3, and that these blocks are stored in lu
c    factored form, with the diagonals stored as reciprocals.
c
c  M = (btc + bjp + bkp)(btc^-1)(btc + bjm + bkm)
c------------------------------------------------------------------------
#include "precis.h"
#include "mpif.h"
#include "mpi_params.f"
      common/iparmi/nt,niter,ntmax,ntime,iflxo,ivis,iturb


      __REAL dq(jmax,kmax,3), rhs(jmax,kmax,3),
     &    btc(jmax,kmax,3,3), bjm(jmax,kmax,3,3), bjp(jmax,kmax,3,3),
     &    bkm(jmax,kmax,3,3), bkp(jmax,kmax,3,3)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Solve (btc + bjm + bkm)*rhs = rhs
c-----
      j = 1
      k = 1
      z1 =  rhs(j,k,1)*btc(j,k,1,1)
      z2 = (rhs(j,k,2) - btc(j,k,2,1)*z1)*btc(j,k,2,2)
      z3 = (rhs(j,k,3) - btc(j,k,3,1)*z1
     &                 - btc(j,k,3,2)*z2)*btc(j,k,3,3)

      rhs(j,k,3) = z3
      rhs(j,k,2) = z2 - btc(j,k,2,3)*rhs(j,k,3)
      rhs(j,k,1) = z1 - btc(j,k,1,3)*rhs(j,k,3)- btc(j,k,1,2)*rhs(j,k,2)

c      if (DEBUG.and.nodeid.eq.4) then
      if (DEBUG.and.nodeid.eq.0) then
       print*
c       print*,' DELETE15.15: nodeid,z1,z2,z3 = ',
c     &  nodeid,z1,z2,z3
c       print*,' DELETE15.15: nodeid, dq(1,1,(1,2,3)) = ', 
c     &  nodeid,dq(j,k,1),dq(j,k,2),dq(j,k,3)
c       print*,' DELETE15.15: nodeid,rhs(1,1,(1,2,3)) = ',
c     &  nodeid,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3)
c       print*,' DELETE15.15: nodeid,btc(j,k,1,1) = ',
c     &  nodeid,btc(j,k,1,1)
c       print*, ' DELETE15.15: nodeid,rhs(2,2,(1,2,3)) = ',
c     &  nodeid,rhs(2,2,1),rhs(2,2,2),rhs(2,2,3)
       print*, ' DELETE15.15: nodeid,nt,rhs(jmax,kmax,(1,2,3) = ',
     &  nodeid,nt,rhs(jmax,kmax,1),rhs(jmax,kmax,2),rhs(jmax,kmax,3)
       print*
       print*, 'BEFORE:'
       print*,'   nt   j   k   rhs1            rsh2            rsh3'
       do 6050 j=1,jmax
c       do 6050 k=1,kmax
       do 6050 k=1,1
        write(*,6055) nt,j,k,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3)
6055    format(3i5,1pe12.5,2x,1pe12.5,2x,1pe12.5)
6050   continue
      print*
      print*,'   nt   j   k   btc11      btc22        btc33'
      do 7060 j=1,jmax
      do 7060 k=1,1
       write(*,7065) nt,j,k,btc(j,k,1,1),btc(j,k,2,2),btc(j,k,3,3)
7065   format(3i5,1pe12.5,2x,1pe12.5,2x,1pe12.5)
7060  continue
      endif

c      stop 'stop: in ilusol.f/ilusol3 120'

c
      k=1
      do 20 j=2,jmax
         rr1 = rhs(j,k,1) - bjm(j,k,1,1)*rhs(j-1,k,1)
     &                    - bjm(j,k,1,2)*rhs(j-1,k,2)
     &                    - bjm(j,k,1,3)*rhs(j-1,k,3)
         rr2 = rhs(j,k,2) - bjm(j,k,2,1)*rhs(j-1,k,1)
     &                    - bjm(j,k,2,2)*rhs(j-1,k,2)
     &                    - bjm(j,k,2,3)*rhs(j-1,k,3)
         rr3 = rhs(j,k,3) - bjm(j,k,3,1)*rhs(j-1,k,1)
     &                    - bjm(j,k,3,2)*rhs(j-1,k,2)
     &                    - bjm(j,k,3,3)*rhs(j-1,k,3)
c
         z1 =  rr1*btc(j,k,1,1)
         z2 = (rr2 - btc(j,k,2,1)*z1)*btc(j,k,2,2)
         z3 = (rr3 - btc(j,k,3,1)*z1 - btc(j,k,3,2)*z2)*btc(j,k,3,3)
         rhs(j,k,3) = z3
         rhs(j,k,2) = z2 - btc(j,k,2,3)*rhs(j,k,3)
         rhs(j,k,1) = z1 - btc(j,k,1,3)*rhs(j,k,3)
     &                   - btc(j,k,1,2)*rhs(j,k,2)
20    continue

      if (DEBUG.and.nodeid.eq.4) then
       print*,' DELETE16.1: nodeid,rhs(jmax,kmax,(1,2,3)) = ',
     &  nodeid,rhs(jmax,kmax,1),rhs(jmax,kmax,2),rhs(jmax,kmax,3)
       print*
       print*, 'BETWEEN'
       print*,'   nt   j   k   rhs1            rsh2            rsh3'
       do 7050 j=1,jmax
       do 7050 k=1,kmax
        write(*,7055) nt,j,k,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3)
7055    format(3i5,1pe12.5,2x,1pe12.5,2x,1pe12.5)
7050   continue
c       stop 'stop: after DELETE16.1'
      endif

c   do 30 i=jmax+1,imax
      if (DEBUG.and.nodeid.eq.4) then
       print*
       print*, ' In do-loop 30:'
       print*, ' nodeid  j   k  btc11   btc12   btc13   btc21    btc22    
     &btc23   btc31   btc32    btc33'
      endif

      do 30 k=2,kmax
      do 30 j=1,jmax

       if (DEBUG.and.nodeid.eq.4) then
        write(*,4015) nodeid,j,k,btc(j,k,1,1),btc(j,k,1,2),btc(j,k,1,3),
     &    btc(j,k,2,1),btc(j,k,2,2),btc(j,k,2,3),btc(j,k,3,1),
     &    btc(j,k,3,2),btc(j,k,3,3)
4015   format(3i6,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,
     &   1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)
       endif

c++++++++++++++++++++++++++++++
        if (j.eq. 1) then
c++++++++++++++++++++++++++++++

         rr1 = rhs(j,k,1) - bjm(j,k,1,1)*rhs(jmax,k-1,1)
     &                    - bjm(j,k,1,2)*rhs(jmax,k-1,2)
     &                    - bjm(j,k,1,3)*rhs(jmax,k-1,3)
     &                    - bkm(j,k,1,1)*rhs(j,k-1,1)
     &                    - bkm(j,k,1,2)*rhs(j,k-1,2)
     &                    - bkm(j,k,1,3)*rhs(j,k-1,3)
         rr2 = rhs(j,k,2) - bjm(j,k,2,1)*rhs(jmax,k-1,1)
     &                    - bjm(j,k,2,2)*rhs(jmax,k-1,2)
     &                    - bjm(j,k,2,3)*rhs(jmax,k-1,3)
     &                    - bkm(j,k,2,1)*rhs(1,k-1,1)
     &                    - bkm(j,k,2,2)*rhs(1,k-1,2)
     &                    - bkm(j,k,2,3)*rhs(1,k-1,3)
         rr3 = rhs(j,k,3) - bjm(j,k,3,1)*rhs(jmax,k-1,1)
     &                    - bjm(j,k,3,2)*rhs(jmax,k-1,2)
     &                    - bjm(j,k,3,3)*rhs(jmax,k-1,3)
     &                    - bkm(j,k,3,1)*rhs(1,k-1,1)
     &                    - bkm(j,k,3,2)*rhs(1,k-1,2)
     &                    - bkm(j,k,3,3)*rhs(1,k-1,3)
      if (DEBUG.and.nodeid.eq.4) then
        write(*,4019) j,k,bjm(j,k,1,1),bjm(j,k,1,2),bjm(j,k,1,3),
     &   bkm(j,k,1,1),bkm(j,k,1,2),bkm(j,k,1,3)
4019  format('-----',2i6,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,
     &  1x,1pe12.5,1x,1pe12.5)
        write(*,4019) j,k,rhs(jmax,k-1,1),rhs(jmax,k-1,2),
     &   rhs(jmax,k-1,3),rhs(1,k-1,1),rhs(1,k-1,2),rhs(1,k-1,3)
        write(*,4020) rr1,rr2,rr2 
4020  format('----rr1,2,3----',1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)
      endif
c
         z1 =  rr1*btc(j,k,1,1)
         z2 = (rr2 - btc(j,k,2,1)*z1)*btc(j,k,2,2)
         z3 = (rr3 - btc(j,k,3,1)*z1 - btc(j,k,3,2)*z2)*btc(j,k,3,3)
      if (DEBUG.and.nodeid.eq.4) then
        write(*,4030) z1,z2,z3 
4030  format('----z1,2,3----',1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)
        write(*,4033) rr3,btc(j,k,3,1),z1,btc(j,k,3,2),z2,btc(j,k,3,3)
4033  format(' --z3 = ',1p6e13.5)
      endif
         
         rhs(j,k,3) = z3
         rhs(j,k,2) = z2 - btc(j,k,2,3)*rhs(j,k,3)
         rhs(j,k,1) = z1 - btc(j,k,1,3)*rhs(j,k,3) 
     &                   - btc(j,k,1,2)*rhs(j,k,2)
       if (DEBUG.and.nodeid.eq.4) then
       write(*,4017) j,k,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3)
4017  format('---new rhs---',2i6,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)
       print*
       print*
       endif

c++++++++++++++++++++++++++++++
        else
c++++++++++++++++++++++++++++++

         rr1 = rhs(j,k,1) - bjm(j,k,1,1)*rhs(j-1,k,1)
     &                    - bjm(j,k,1,2)*rhs(j-1,k,2)
     &                    - bjm(j,k,1,3)*rhs(j-1,k,3)
     &                    - bkm(j,k,1,1)*rhs(j,k-1,1)
     &                    - bkm(j,k,1,2)*rhs(j,k-1,2)
     &                    - bkm(j,k,1,3)*rhs(j,k-1,3)
         rr2 = rhs(j,k,2) - bjm(j,k,2,1)*rhs(j-1,k,1)
     &                    - bjm(j,k,2,2)*rhs(j-1,k,2)
     &                    - bjm(j,k,2,3)*rhs(j-1,k,3)
     &                    - bkm(j,k,2,1)*rhs(j,k-1,1)
     &                    - bkm(j,k,2,2)*rhs(j,k-1,2)
     &                    - bkm(j,k,2,3)*rhs(j,k-1,3)
         rr3 = rhs(j,k,3) - bjm(j,k,3,1)*rhs(j-1,k,1)
     &                    - bjm(j,k,3,2)*rhs(j-1,k,2)
     &                    - bjm(j,k,3,3)*rhs(j-1,k,3)
     &                    - bkm(j,k,3,1)*rhs(j,k-1,1)
     &                    - bkm(j,k,3,2)*rhs(j,k-1,2)
     &                    - bkm(j,k,3,3)*rhs(j,k-1,3)

      if (DEBUG.and.nodeid.eq.4) then
        write(*,4019) j,k,bjm(j,k,1,1),bjm(j,k,1,2),bjm(j,k,1,3),
     &   bkm(j,k,1,1),bkm(j,k,1,2),bkm(j,k,1,3)
        write(*,4019) j,k,rhs(jmax,k-1,1),rhs(jmax,k-1,2),
     &   rhs(jmax,k-1,3),rhs(1,k-1,1),rhs(1,k-1,2),rhs(1,k-1,3)
        write(*,4020) rr1,rr2,rr2 
      endif
c
         z1 =  rr1*btc(j,k,1,1)
         z2 = (rr2 - btc(j,k,2,1)*z1)*btc(j,k,2,2)
         z3 = (rr3 - btc(j,k,3,1)*z1 - btc(j,k,3,2)*z2)*btc(j,k,3,3)
      if (DEBUG.and.nodeid.eq.4) then
        write(*,4030) z1,z2,z3 
      endif
      
         rhs(j,k,3) = z3
         rhs(j,k,2) = z2 - btc(j,k,2,3)*rhs(j,k,3)
         rhs(j,k,1) = z1 - btc(j,k,1,3)*rhs(j,k,3) 
     &                   - btc(j,k,1,2)*rhs(j,k,2)
       if (DEBUG.and.nodeid.eq.4) then
        write(*,4017) j,k,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3)
        print*
        print*
       endif

c++++++++++++++++++++++++++++++
        endif
c++++++++++++++++++++++++++++++

30    continue

      if (DEBUG.and.nodeid.eq.4) then
       print*,' DELETE16.2: nodeid,rhs(jmax,kmax,(1,2,3)) = ',
     &  nodeid,rhs(jmax,kmax,1),rhs(jmax,kmax,2),rhs(jmax,kmax,3)
       print*
       print*, 'BETWEEN 2'
       print*,'   nt   j   k   rhs1            rhs2            rhs3'
       do 8050 j=1,jmax
       do 8050 k=1,kmax
        write(*,8055) nt,j,k,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3)
8055    format(3i5,1pe12.5,2x,1pe12.5,2x,1pe12.5)
8050   continue
c       stop 'stop: after DELETE16.2'
      endif

c-----
c  Solve btc^-1*rhs = rhs (Multiply by btc)
c-----
c      do 40 i=1,imax
      do 40 k=1,kmax
      do 40 j=1,jmax
         b11 = 1./btc(j,k,1,1)
         b12 = btc(j,k,1,2)*b11
         b13 = btc(j,k,1,3)*b11
         b21 = btc(j,k,2,1)
         b22 = 1./btc(j,k,2,2) + btc(j,k,2,1)*btc(j,k,1,2)
         b23 = btc(j,k,2,3)/btc(j,k,2,2) + btc(j,k,2,1)*btc(j,k,1,3)
         b31 = btc(j,k,3,1)
         b32 = btc(j,k,3,2) + btc(j,k,3,1)*btc(j,k,1,2)
         b33 = 1./btc(j,k,3,3) + btc(j,k,3,1)*btc(j,k,1,3)
     &                       + btc(j,k,3,2)*btc(j,k,2,3)
c         
         s1 = b11*rhs(j,k,1)
     &      + b12*rhs(j,k,2)
     &      + b13*rhs(j,k,3)
         s2 = b21*rhs(j,k,1)
     &      + b22*rhs(j,k,2)
     &      + b23*rhs(j,k,3)
         s3 = b31*rhs(j,k,1)
     &      + b32*rhs(j,k,2)
     &      + b33*rhs(j,k,3)
         rhs(j,k,1) = s1
         rhs(j,k,2) = s2
         rhs(j,k,3) = s3

c      if (j.eq.2.and.k.eq.2.and.nodeid.eq.4) then
      if (DEBUG.and.k.eq.2.and.nodeid.eq.4) then
c       print*,' DELETE15.0: j,b11,b12,b13 = ',j,b11,b12,b13
c       print*,' DELETE15.0: j,b21,b22,b23 = ',j,b21,b22,b23
c       print*,' DELETE15.0: j,b31,b32,b33 = ',j,b31,b32,b33
c       print*,' DELETE15.0: j,btc(j,2,2,2) = ',j,btc(j,2,2,2)
c       print*,' ------------------------------'
      endif 
40    continue

      if (DEBUG.and.nodeid.eq.4) then
       print*,' DELETE16.3: nodeid,rhs(jmax,kmax,(1,2,3)) = ',
     &  nodeid,rhs(jmax,kmax,1),rhs(jmax,kmax,2),rhs(jmax,kmax,3)
c       stop 'stop: after DELETE16.3'
      endif

c-----
c  Solve (btc + bjp + bkp)*dq = rhs
c-----
c      i = imax
      k = kmax
      j = jmax
      z1 =  rhs(j,k,1)*btc(j,k,1,1)
      z2 = (rhs(j,k,2) - btc(j,k,2,1)*z1)*btc(j,k,2,2)
      z3 = (rhs(j,k,3) - btc(j,k,3,1)*z1 - btc(j,k,3,2)*z2)*btc(j,k,3,3)
      dq(j,k,3) = z3
      dq(j,k,2) = z2 - btc(j,k,2,3)*dq(j,k,3)
      dq(j,k,1) = z1 - btc(j,k,1,3)*dq(j,k,3) - btc(j,k,1,2)*dq(j,k,2)

      if (DEBUG.and.nodeid.eq.4) then
       print*
       print*,' DELETE15.1: nt,nodeid, dq(jmax,kmax,(1,2,3)) = ', 
     &  nt,nodeid,dq(jmax,kmax,1),dq(jmax,kmax,2),dq(jmax,kmax,3)
       print*,' DELETE15.1: nt,nodeid,rhs(jmax,kmax,2),z2,z3 = ',
     &  nt,nodeid,rhs(jmax,kmax,2),z2,z3
       print*,' DELETE15.1: nt,btc(jmax,kmax,2,2) = ',
     &  nt,btc(jmax,kmax,2,2)

       print*
       print*, 'AFTER'
       print*,'   nt   j   k   rhs1            rsh2            rsh3'
       do 5050 j=1,jmax
       do 5050 k=1,kmax
        write(*,5055) nt,j,k,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3)
5055    format(3i5,1pe12.5,2x,1pe12.5,2x,1pe12.5)
5050   continue
      call flush(6)
c      stop 'stop: in ilusol.f/ilusol3 150'
      endif
c

c      do 120 i=imax-1,imax-jmax+1,-1
      k=kmax
      do 120 j=jmax-1,1,-1
         rr1 = rhs(j,k,1) - bjp(j,k,1,1)*dq(j+1,k,1)
     &                  - bjp(j,k,1,2)*dq(j+1,k,2)
     &                  - bjp(j,k,1,3)*dq(j+1,k,3)
         rr2 = rhs(j,k,2) - bjp(j,k,2,1)*dq(j+1,k,1)
     &                  - bjp(j,k,2,2)*dq(j+1,k,2)
     &                  - bjp(j,k,2,3)*dq(j+1,k,3)
         rr3 = rhs(j,k,3) - bjp(j,k,3,1)*dq(j+1,k,1)
     &                  - bjp(j,k,3,2)*dq(j+1,k,2)
     &                  - bjp(j,k,3,3)*dq(j+1,k,3)
c
         z1 =  rr1*btc(j,k,1,1)
         z2 = (rr2 - btc(j,k,2,1)*z1)*btc(j,k,2,2)
         z3 = (rr3 - btc(j,k,3,1)*z1 - btc(j,k,3,2)*z2)*btc(j,k,3,3)
         dq(j,k,3) = z3
         dq(j,k,2) = z2 - btc(j,k,2,3)*dq(j,k,3)
         dq(j,k,1) = z1 - btc(j,k,1,3)*dq(j,k,3)
     &                - btc(j,k,1,2)*dq(j,k,2)
120   continue

      if (DEBUG.and.nodeid.eq.4) then
       print*
       print*,' DELETE15.2: nodeid,dq(2,kmax,2),dq(2,kmax,3) = ', 
     &  nodeid,dq(2,kmax,2),dq(2,kmax,3)
      endif

c      do 130 i=imax-jmax,1,-1  ! kmax-1,1,-1; jmax,1,-1
      do 130 k=kmax-1,1,-1
      do 130 j=jmax,1,-1
        if (j.eq.jmax) then
         rr1 = rhs(j,k,1) - bkp(j,k,1,1)*dq(j,k+1,1)
     &                    - bkp(j,k,1,2)*dq(j,k+1,2)
     &                    - bkp(j,k,1,3)*dq(j,k+1,3)
     &                    - bjp(j,k,1,1)*dq(j,kmax,1)
     &                    - bjp(j,k,1,2)*dq(j,kmax,2)
     &                    - bjp(j,k,1,3)*dq(j,kmax,3)

         rr2 = rhs(j,k,2) - bkp(j,k,2,1)*dq(j,k+1,1)
     &                    - bkp(j,k,2,2)*dq(j,k+1,2)
     &                    - bkp(j,k,2,3)*dq(j,k+1,3)
     &                    - bjp(j,k,2,1)*dq(j,kmax,1)
     &                    - bjp(j,k,2,2)*dq(j,kmax,2)
     &                    - bjp(j,k,2,3)*dq(j,kmax,3)

         rr3 = rhs(j,k,3) - bkp(j,k,3,1)*dq(j,k+1,1) ! j,k+1
     &                    - bkp(j,k,3,2)*dq(j,k+1,2)
     &                    - bkp(j,k,3,3)*dq(j,k+1,3)
     &                    - bjp(j,k,3,1)*dq(1,kmax,1)   ! 1,kmax
     &                    - bjp(j,k,3,2)*dq(1,kmax,2)
     &                    - bjp(j,k,3,3)*dq(1,kmax,3)

       else

         rr1 = rhs(j,k,1) - bkp(j,k,1,1)*dq(j,k+1,1)
     &                    - bkp(j,k,1,2)*dq(j,k+1,2)
     &                    - bkp(j,k,1,3)*dq(j,k+1,3)
     &                    - bjp(j,k,1,1)*dq(j+1,k,1)
     &                    - bjp(j,k,1,2)*dq(j+1,k,2)
     &                    - bjp(j,k,1,3)*dq(j+1,k,3)

         rr2 = rhs(j,k,2) - bkp(j,k,2,1)*dq(j,k+1,1)
     &                    - bkp(j,k,2,2)*dq(j,k+1,2)
     &                    - bkp(j,k,2,3)*dq(j,k+1,3)
     &                    - bjp(j,k,2,1)*dq(j+1,k,1)
     &                    - bjp(j,k,2,2)*dq(j+1,k,2)
     &                    - bjp(j,k,2,3)*dq(j+1,k,3)

         rr3 = rhs(j,k,3) - bkp(j,k,3,1)*dq(j,k+1,1)
     &                    - bkp(j,k,3,2)*dq(j,k+1,2)
     &                    - bkp(j,k,3,3)*dq(j,k+1,3)
     &                    - bjp(j,k,3,1)*dq(j+1,k,1)
     &                    - bjp(j,k,3,2)*dq(j+1,k,2)
     &                    - bjp(j,k,3,3)*dq(j+1,k,3)

       endif
c
         z1 =  rr1*btc(j,k,1,1)
         z2 = (rr2 - btc(j,k,2,1)*z1)*btc(j,k,2,2)
         z3 = (rr3 - btc(j,k,3,1)*z1 - btc(j,k,3,2)*z2)*btc(j,k,3,3)
         dq(j,k,3) = z3
         dq(j,k,2) = z2 - btc(j,k,2,3)*dq(j,k,3)
         dq(j,k,1) = z1 - btc(j,k,1,3)*dq(j,k,3)
     &                - btc(j,k,1,2)*dq(j,k,2)
130   continue

      if (DEBUG.and.nodeid.eq.4) then
       print*
       print*,' DELETE15.3: nodeid, dq(2,2,2) = ', 
     &  nodeid,dq(2,2,2)
      endif

c-----
c   End of ilusol3
c-----
      return
      end

c
c
c************************************************************************
      subroutine ilusol1(jmax,kmax,dq,rhs,btc,bjm,bjp,bkm,bkp)
c************************************************************************
c  Scalar version of this routine.
c
c  Approximate solve of system M*dq = rhs, where M is
c  ilu decomposition of lhs of original system.
c
c  This routine assumes that btc contains the modified main diagional
c    computed by iluini1, and that these elements are stored as reciprocals.
c
c  M = (btc + bjp + bkp)(btc^-1)(btc + bjm + bkm)
c------------------------------------------------------------------------
#include "precis.h"
#include "mpi_params.f"
      __REAL dq(jmax,kmax), rhs(jmax,kmax),
     &    btc(jmax,kmax), bjm(jmax,kmax), bjp(jmax,kmax),
     &    bkm(jmax,kmax), bkp(jmax,kmax)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Solve (btc + bjm + bkm)*rhs = rhs
c-----
      i = 1
      rhs(i,i) = rhs(i,i)*btc(i,i)
c
      do 20 i=2,jmax
         rhs(i,1) = btc(i,1)*( rhs(i,1) - bjm(i,1)*rhs(i-1,1) )
20    continue

c      do 30 i=jmax+1,imax
      do 30 k=2,kmax
      do 30 j=1,jmax
       if (j.eq.1) then
         rhs(j,k) = btc(j,k)*( rhs(j,k) - bjm(j,k)*rhs(jmax,k-1)
     &                            - bkm(j,k)*rhs(j,k-1) )
       else
         rhs(j,k) = btc(j,k)*( rhs(j,k) - bjm(j,k)*rhs(j-1,k)
     &                            - bkm(j,k)*rhs(j,k-1) )
       endif
30    continue
c-----
c  Solve btc^-1*rhs = rhs (Multiply by btc)
c-----
      do 40 k=1,kmax
      do 40 j=1,jmax
         rhs(j,k) = rhs(j,k)/btc(j,k)
40    continue
c-----
c  Solve (btc + bjp + bkp)*dq = rhs
c-----
c      i = imax
      dq(jmax,kmax) = rhs(jmax,kmax)*btc(jmax,kmax)
c
c  j=jmax-1,1,-1, k=kmax
c
c      do 120 i=imax-1,imax-jmax+1,-1   ! jmax-1,1,-1, at k = kmax
      k=kmax
      do 120 j=jmax-1,1,-1
         dq(j,k) = btc(j,k)*( rhs(j,k) - bjp(j,k)*dq(j+1,k) )
120   continue
c
c  j=jmax,1,-1; k=kmax-1,1,-1
c
      if (DEBUG.and.nodeid.eq.4) then
       print*,'j  k    btc   rhs    bkp   bjp  dq(j1,k)   dq(j,k1)'
      endif
c      do 130 i=imax-jmax,1,-1
      do 130 k=kmax-1,1,-1
      do 130 j=jmax,1,-1
      if (DEBUG.and.nodeid.eq.4.and.j.ne.jmax) then
        write(*,131) j,k,btc(j,k),rhs(j,k),bkp(j,k),bjp(j,k),
     &   dq(j+1,k),dq(j,k+1)
131    format(i3,i3,1p6e12.5)
      endif

       if (j.eq.jmax) then
         dq(j,k) = btc(j,k)*( rhs(j,k) - bkp(j,k)*dq(j,k+1)
     &                           - bjp(j,k)*dq(1,k+1) )
       else
         dq(j,k) = btc(j,k)*( rhs(j,k) - bkp(j,k)*dq(j,k+1)
     &                           - bjp(j,k)*dq(j+1,k) )
       endif

130   continue

      if (DEBUG.and.nodeid.eq.4) then
      call flush(6)
      print*
c      print*, ' From ilusol/ilusol1.f with neqs =',neqs
c      print*, ' node    j    k   btc   bjm  rhs  dq    '
c      do 200 j=1,jmax
c      do 200 k=1,kmax
c       print*,  nodeid,j,k,btc(j,k),bjm(j,k),rhs(j,k),dq(j,k)
c200   continue
c      stop 'stop: at end of ilusol.f/ilusol1'
      endif
        

c-----
c   End of ilusol1
c-----
      return
      end
