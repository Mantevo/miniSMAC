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
c  1. lrelax
c  2. ljsweep
c  3. lksweep
c
c************************************************************
c
c
c************************************************************************
      subroutine lrelax(jmax,kmax,q,dq,s,
     &                  btc,bjm,bjp,bkm,bkp)
c************************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"

#ifdef _OPENMP
#include "omp_lib.h"
#endif

      __INTEGER stat(MPI_STATUS_SIZE)

      __REAL q(jmax,kmax,3),dq(jmax,kmax,3),s(jmax,kmax,3),
     &  btc(jmax,kmax,3,3),bjm(jmax,kmax,3,3),bjp(jmax,kmax,3,3),
     &  bkm(jmax,kmax,3,3), bkp(jmax,kmax,3,3)
      __INTEGER jmax, kmax
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Zero out dq array 
c-----

!$omp parallel do collapse(2) shared(kmax,jmax)
      do 9 i=1,3
      do 9 k=1,kmax
      do 9 j=1,jmax
         dq(j,k,i) = 0.0
9     continue
!$omp end parallel do

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, 'DWB10.5: BEFORE call to btlu3j'
       print*, 'DWB10.5: jmax,kmax,dq(2,1,1) = ',jmax,kmax,dq(2,1,1)
       print*, 'DWB10.5: btc(2,2,2,2) =',btc(2,2,2,2)
c       print*, 'DWB10.5: bjm(2,2,2,2), btc(2,2,2,2), bjp(2,2,2,2) = ',
c     &  bjm(2,2,2,2),btc(2,2,2,2),bjp(2,2,2,2)
       print*, 'DWB10.5: btc(2,2,1,1) = ',btc(2,2,1,1)
       print*
c       stop 'stop: in miniSMAC lrelax.f/lrelax 10'
      endif 

c-----
c  J-sweep: solve for lines of varying j and constant k.
c  First perform lu decomposition of each j-line for all grid
c  points in zones where j-sweeps will be done.
c-----
      njsmax = 0
c njsp = njswp = 2: set in initia.f
         njsmax = max( njsmax, njsp )
         if(njsp .gt. 0) then
           call btlu3j(jmax,kmax,bjm,btc,bjp,1,jmax,1,kmax)
         endif
      if (DEBUG.and.nodeid.eq.0) then
       print*, 'DWB10.6: AFTER call to btlu3j'
       print*, 'DWB10.6: jmax,kmax,dq(2,1,1) = ',jmax,kmax,dq(2,1,1)
       print*, 'DWB10.6: bjm(2,2,2,2), btc(2,2,2,2), bjp(2,2,2,2) = ',
     &  bjm(2,2,2,2),btc(2,2,2,2),bjp(2,2,2,2)
      endif 

c-----
c  Top of j-sweep loop
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*, ' DWB11.5: nodeid,njsmax,njsp = ',
     &    nodeid,njsmax,njsp
        print*, ' DWB11.5: nodeid,jmax,kmax,dq(2,1,1) = ',
     &    jmax,kmax,dq(2,1,1)
      endif

c default value of njsmax = 2 sweeps; njsp=njswp=2

      do 30 nj=1,njsmax
           if (DEBUG.and.nodeid.eq.0) then
            print*, 'DWB11: N5 - before ljsweep, nj,dq(2,1,1) = ',
     &            nj,dq(2,1,1)
            print*, 'DWB11: N5 - before ljsweep, nj,dq(10,1,2) = ',
     &            nj,dq(10,1,2)
           endif
            if(njsp .ge. nj) then
             call ljsweep(jmax,kmax,1,jmax,1,kmax,nodeid+1,nj,
     &        dq,s,btc,bjm,bjp,bkm,bkp
     &        )
              if (DEBUG.and.nodeid.eq.0) then
               print*, 'DWB9: after ljsweep, njsp,nj,dq(2,1,1) = ',
     &          njsp,nj,dq(2,1,1)
                print*, 'DWB9: after ljsweep, njsp,nj,dq(10,1,2) = ',
     &          njsp,nj,dq(10,1,2)
c               stop 'stop: after ljsweep in lrelax.f'
              endif

c-----
c  Pass boundary info between grids after each sweep
c    Note: default nbcimp = 1 in initia.f
c-----
               neqs = 3
               if( (nj/nbcimp)*nbcimp .eq. nj) then
               call bcimpds(jmax,kmax,neqs,q,dq,s,nodeid+1)
               if (DEBUG.and.nodeid.eq.0) then
                print*, 'DWB10: after bcimpds, nodeid,njsp,nj,dq(2,1,1), 
     &nt = ',
     &              nodeid,njsp,nj,dq(2,1,1),nt
                print*, 'DWB10: after bcimpds,nodeid,njsp,nj,dq(10,1,2),
     &nt = ',
     &              nodeid,njsp,nj,dq(10,1,2),nt
               endif
               endif
            endif
30    continue

      if (DEBUG.and.nt.ge.1) then 
c       stop 'stop: after bcimpds in lrelax.f'
      endif

c-----
c  Undo decomposition for j-sweeps
c-----
c njsp = 2 --> default defined in initia.f; defines number of 
c              sweeps in j or k direction
         if(njsp .gt. 0) then
          call btund3j(jmax,kmax,bjm,btc,bjp,
     &     1,jmax,1,kmax)
         endif

c-----
c  K-sweep: solve for lines of varying k and constant j.
c  First perform lu decomposition of each k-line for all grid
c  points in zones where k-sweeps will be done.
c-----
      nksmax = 0
c NOTE: nksp = 2, set in initia.f, identical to nkswp
         nksmax = max( nksmax, nksp )
         if(nksp .gt. 0) then
               call btlu3k(jmax,kmax,bkm,
     &            btc,bkp,1,jmax,1,kmax)
         endif

c-----
c  Top of k-sweep loop
c-----
      do 70 nk=1,nksmax
            if(nksp .ge. nk) then
                  call lksweep(jmax,kmax,1,jmax,1,kmax,
     &            nodeid+1,nk,dq,s,btc,bjm,bjp,bkm,bkp
     &            )
c-----
c  Pass boundary info between grids
c    Note: default nbcimp = 1 in initia.f
c-----
               neqs = 3
               if( (nk/nbcimp)*nbcimp .eq. nk)
     &         call bcimpds(jmax,kmax,neqs,q,dq,s,nodeid+1)

            endif
70    continue

c-----
c  Undo decomposition for k-sweeps
c-----
         if(nksp .gt. 0)
     &      call btund3k(jmax,kmax,bkm,btc,bkp,1,jmax,1,kmax)

c-----
c  End of lrelax
c-----
      return
      end
c
c
c****************************************************************
      subroutine ljsweep(jmax,kmax,jbeg,jend,kbeg,kend,nz,nj,dq,s,
     &                   btc,bjm,bjp,bkm,bkp)
c****************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __REAL btc(jmax,kmax,3,3), bjm(jmax,kmax,3,3),
     &     bjp(jmax,kmax,3,3), bkm(jmax,kmax,3,3),
     &     bkp(jmax,kmax,3,3), dq(jmax,kmax,3),
     &       s(jmax,kmax,3)
      common/btri/a(jkmax,3,3),b(jkmax,3,3),c(jkmax,3,3),rhs(jkmax,3),
     &            d(jkmax,3,3),e(jkmax,3,3)
c      double precision, parameter :: half=0.5
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Form forcing function by moving off-line implicit terms*delta q
c  to rhs and then solve system one line at a time.
c-----
      half = 0.5d0

      if( (nj/2)*2 .eq. nj) then
c sweep in -k direction
         kb = kend
         ke = kbeg
         ki = -1
         if (DEBUG) then
          if (nz-1.eq.8) then
            print*
            print*, ' DWB99: In lrelax.f/ljsweep: node,nj,kb,ke,ki = ',
     &        nz-1,nj,kb,ke,ki
            print*
          endif
         endif
      else
c sweep in +k direction
         kb = kbeg
         ke = kend
         ki = 1
         if (DEBUG) then
          if (nz-1.eq.8) then
            print*
            print*, ' DWB98: In lrelax.f/ljsweep: node,nj,kb,ke,ki = ',
     &        nz-1,nj,kb,ke,ki
            print*
          endif
         endif
      endif

c sweep over k
      do 20 kk=kb,ke,ki
         k = nint( mod( float(kk+kmax), float(kmax)+.1 ) )
         km = max(k-1, 1)
         kp = min(k+1, kmax)
         rkm = float(k-1-km)
         rkp = float(kp-k-1)
c NOTE: 'sign(a,b)' returns the value of 'a' with the sign of 'b'
         swkm = half + sign(half, rkm)
         swkp = half + sign(half, rkp)


         do 12 j=jbeg,jend
            rhs(j,1) = s(j,k,1)
     &      - ( bkm(j,k,1,1)*dq(j,km,1) + bkm(j,k,1,2)*dq(j,km,2)
     &        + bkm(j,k,1,3)*dq(j,km,3) )*swkm
     &      - ( bkp(j,k,1,1)*dq(j,kp,1) + bkp(j,k,1,2)*dq(j,kp,2)
     &        + bkp(j,k,1,3)*dq(j,kp,3) )*swkp
            rhs(j,2) = s(j,k,2)
     &      - ( bkm(j,k,2,1)*dq(j,km,1) + bkm(j,k,2,2)*dq(j,km,2)
     &        + bkm(j,k,2,3)*dq(j,km,3) )*swkm
     &      - ( bkp(j,k,2,1)*dq(j,kp,1) + bkp(j,k,2,2)*dq(j,kp,2)
     &        + bkp(j,k,2,3)*dq(j,kp,3) )*swkp
            rhs(j,3) = s(j,k,3)
     &      - ( bkm(j,k,3,1)*dq(j,km,1) + bkm(j,k,3,2)*dq(j,km,2)
     &        + bkm(j,k,3,3)*dq(j,km,3) )*swkm
     &      - ( bkp(j,k,3,1)*dq(j,kp,1) + bkp(j,k,3,2)*dq(j,kp,2)
     &        + bkp(j,k,3,3)*dq(j,kp,3) )*swkp
      
12       continue

!$omp parallel do
!$omp+ shared(jbeg,jend)
         do 15 m=1,3
         do 15 j=jbeg,jend
            a(j,m,1) = bjm(j,k,m,1)
            b(j,m,1) = btc(j,k,m,1)
            c(j,m,1) = bjp(j,k,m,1)
            a(j,m,2) = bjm(j,k,m,2)
            b(j,m,2) = btc(j,k,m,2)
            c(j,m,2) = bjp(j,k,m,2)
            a(j,m,3) = bjm(j,k,m,3)
            b(j,m,3) = btc(j,k,m,3)
            c(j,m,3) = bjp(j,k,m,3)
15       continue
!$omp end parallel do

c
c  Inlining of btso3j
c
c-----
c  Forward sweep
c-----

         j = jbeg
         d1 = b(j,1,1)* rhs(j,1)
         d2 = b(j,2,2)*(rhs(j,2) - b(j,2,1)*d1)
         d3 = b(j,3,3)*(rhs(j,3) - b(j,3,1)*d1 - b(j,3,2)*d2)
         rhs(j,3) = d3
         rhs(j,2) = d2 - b(j,2,3)*rhs(j,3)
         rhs(j,1) = d1 - b(j,1,3)*rhs(j,3) - b(j,1,2)*rhs(j,2)

         do 120 j=jbeg+1,jend
            do 110 n=1,3
               rhs(j,n) = rhs(j,n) - a(j,n,1)*rhs(j-1,1)
     &                             - a(j,n,2)*rhs(j-1,2)
     &                             - a(j,n,3)*rhs(j-1,3)
110         continue
            d1 = b(j,1,1)* rhs(j,1)
            d2 = b(j,2,2)*(rhs(j,2) - b(j,2,1)*d1)
            d3 = b(j,3,3)*(rhs(j,3) - b(j,3,1)*d1 - b(j,3,2)*d2)
            rhs(j,3) = d3
            rhs(j,2) = d2 - b(j,2,3)*rhs(j,3)
            rhs(j,1) = d1 - b(j,1,3)*rhs(j,3) - b(j,1,2)*rhs(j,2)
120      continue

c-----
c  Backward substitution
c-----

c!$omp single

      if (DEBUG.and.nodeid.eq.0.and.k.le.2) then
       print*
       print*, ' j  k  n  rhs_jn'
      endif

c!$omp end single


         do 140 j=jend-1,jbeg,-1
            do 130 n=1,3
               rhs(j,n) = rhs(j,n) - c(j,n,1)*rhs(j+1,1)
     &                             - c(j,n,2)*rhs(j+1,2)
     &                             - c(j,n,3)*rhs(j+1,3)
130         continue
140      continue

!$omp parallel do
!$omp+ shared(jbeg,jend,underr)
         do 10 j=jbeg,jend
            dq(j,k,1) = dq(j,k,1) + underr*(rhs(j,1) - dq(j,k,1))
            dq(j,k,2) = dq(j,k,2) + underr*(rhs(j,2) - dq(j,k,2))
            dq(j,k,3) = dq(j,k,3) + underr*(rhs(j,3) - dq(j,k,3))
10       continue
!$omp end parallel do

20    continue


      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' nodeid  j   k   dq(jk1)     dq(jk2)      dq(jk3)'
       if (kmax.gt.10) then
        kmaxx = 10
       else
        kmaxx = kmax
       endif
       if (jmax.gt.17) then
        jmaxx = 17
       else
        jmaxx = jmax
       endif
       do 1011 k=1,kmaxx
       do 1011 j=1,jmaxx
        write(*,1012) nodeid,j,k,dq(j,k,1),dq(j,k,2),dq(j,k,3)
1012    format(3i5,1p3e13.5)
1011   continue
       call flush(6)
       call flush(istdout)
c       stop 'stop: at end of lrelax.f/ljsweep'
      endif

c-----
c  End of ljsweep
c-----
      return
      end
c
c
c****************************************************************
      subroutine lksweep(jmax,kmax,jbeg,jend,kbeg,kend,nz,nk,dq,s,
     &                   btc,bjm,bjp,bkm,bkp)
c****************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __REAL btc(jmax,kmax,3,3), bjm(jmax,kmax,3,3),
     &     bjp(jmax,kmax,3,3), bkm(jmax,kmax,3,3),
     &     bkp(jmax,kmax,3,3), dq(jmax,kmax,3),
     &       s(jmax,kmax,3)
      common/btri/a(jkmax,3,3),b(jkmax,3,3),c(jkmax,3,3),rhs(jkmax,3),
     &            d(jkmax,3,3),e(jkmax,3,3)
c      double precision, parameter :: half=0.5
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c-----
c  Non-periodic
c-----
c-----
c  Form forcing function by moving off-line implicit terms*delta q
c  to rhs and then solve system one line at a time.
c-----
      half = 0.5d0

      if( (nk/2)*2 .eq. nk) then
         jb = jend
         je = jbeg
         ji = -1
       else
         jb = jbeg
         je = jend
         ji = 1
      endif


      do 40 j=jb,je,ji
         jm = max(j-1, 1)
         jp = min(j+1, jmax)
         rjm = float(j-1-jm)
         rjp = float(jp-j-1)
         swjm = half + sign(half, rjm)
         swjp = half + sign(half, rjp)

c         print *, ' for j = ',j,', kbeg,jb,je,ji = ',kbeg,jb,je,ji
c         print *, '    jm,jp,rjm,rjp,half,underr = ',jm,jp,rjm,rjp,
c     1     half,underr

         do 26 k=kbeg,kend
            rhs(k,1) = s(j,k,1)
     &      - ( bjm(j,k,1,1)*dq(jm,k,1) + bjm(j,k,1,2)*dq(jm,k,2)
     &        + bjm(j,k,1,3)*dq(jm,k,3) )*swjm
     &      - ( bjp(j,k,1,1)*dq(jp,k,1) + bjp(j,k,1,2)*dq(jp,k,2)
     &        + bjp(j,k,1,3)*dq(jp,k,3) )*swjp
            rhs(k,2) = s(j,k,2)
     &      - ( bjm(j,k,2,1)*dq(jm,k,1) + bjm(j,k,2,2)*dq(jm,k,2)
     &        + bjm(j,k,2,3)*dq(jm,k,3) )*swjm
     &      - ( bjp(j,k,2,1)*dq(jp,k,1) + bjp(j,k,2,2)*dq(jp,k,2)
     &        + bjp(j,k,2,3)*dq(jp,k,3) )*swjp
            rhs(k,3) = s(j,k,3)
     &      - ( bjm(j,k,3,1)*dq(jm,k,1) + bjm(j,k,3,2)*dq(jm,k,2)
     &        + bjm(j,k,3,3)*dq(jm,k,3) )*swjm
     &      - ( bjp(j,k,3,1)*dq(jp,k,1) + bjp(j,k,3,2)*dq(jp,k,2)
     &        + bjp(j,k,3,3)*dq(jp,k,3) )*swjp
26       continue

!$omp parallel do
!$omp+ shared(kbeg,kend)
         do 35 m=1,3
         do 35 k=kbeg,kend
            a(k,m,1) = bkm(j,k,m,1)
            b(k,m,1) = btc(j,k,m,1)
            c(k,m,1) = bkp(j,k,m,1)
            a(k,m,2) = bkm(j,k,m,2)
            b(k,m,2) = btc(j,k,m,2)
            c(k,m,2) = bkp(j,k,m,2)
            a(k,m,3) = bkm(j,k,m,3)
            b(k,m,3) = btc(j,k,m,3)
            c(k,m,3) = bkp(j,k,m,3)
35       continue
!$omp end parallel do

         k = kbeg
c         print *, ' for j = ',j,', k = kbeg,jb,je,ji = ',k,jb,je,ji
c         print *, '    jm,jp,rjm,rjp,half,underr = ',jm,jp,rjm,rjp,
c     1     half,underr

         d1 = b(k,1,1)* rhs(k,1)
         d2 = b(k,2,2)*(rhs(k,2) - b(k,2,1)*d1)
         d3 = b(k,3,3)*(rhs(k,3) - b(k,3,1)*d1 - b(k,3,2)*d2)
         rhs(k,3) = d3
         rhs(k,2) = d2 - b(k,2,3)*rhs(k,3)
         rhs(k,1) = d1 - b(k,1,3)*rhs(k,3) - b(k,1,2)*rhs(k,2)

         do 120 k=kbeg+1,kend
            do 110 n=1,3
               rhs(k,n) = rhs(k,n) - a(k,n,1)*rhs(k-1,1)
     &                             - a(k,n,2)*rhs(k-1,2)
     &                             - a(k,n,3)*rhs(k-1,3)
110          continue
            d1 = b(k,1,1)* rhs(k,1)
            d2 = b(k,2,2)*(rhs(k,2) - b(k,2,1)*d1)
            d3 = b(k,3,3)*(rhs(k,3) - b(k,3,1)*d1 - b(k,3,2)*d2)
            rhs(k,3) = d3
            rhs(k,2) = d2 - b(k,2,3)*rhs(k,3)
            rhs(k,1) = d1 - b(k,1,3)*rhs(k,3) - b(k,1,2)*rhs(k,2)
120       continue

c-----
c  Backward substitution
c-----

!$omp parallel do 
!$omp+ shared(kbeg,kend)
         do 140 k=kend-1,kbeg,-1
            do 130 n=1,3
               rhs(k,n) = rhs(k,n) - c(k,n,1)*rhs(k+1,1)
     &                             - c(k,n,2)*rhs(k+1,2)
     &                             - c(k,n,3)*rhs(k+1,3)
130          continue
140       continue
!$omp end parallel do

c
c End of btso3k inlining
c

!$omp parallel do
!$omp+ shared(kbeg,kend,underr)
         do 30 k=kbeg,kend
            dq(j,k,1) = dq(j,k,1) + underr*(rhs(k,1) - dq(j,k,1))
            dq(j,k,2) = dq(j,k,2) + underr*(rhs(k,2) - dq(j,k,2))
            dq(j,k,3) = dq(j,k,3) + underr*(rhs(k,3) - dq(j,k,3))
30       continue
!$omp end parallel do

40    continue

c-----
c  End of lksweep
c-----
      return
      end
