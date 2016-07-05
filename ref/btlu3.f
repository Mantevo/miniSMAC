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
c  1. btlu3j
c  2. btlu3k
c  3. btund3j
c  4. btund3k 
c
c************************************************************
c
c
c**************************************************************
      subroutine btlu3j(jmax,kmax,a,b,c,jbeg,jend,kbeg,kend)
c**************************************************************
c   LU decomposition of block tridiagonal system of equations
c   Upper blocks of upper blocks are overwritten in c.
c   LU decompostion of each diagonal block of lower matrix are
c   overwritten in b, and the diagonal entries of these blocks
c   are stored as reciprocals.
c
c   block size = (3x3)
c
c------------------------------------------------------------------
#include "precis.h"
#include "mpif.h"
#include "mpi_params.f"

#ifdef _OPENMP
#include "omp_lib.h"
#endif

      __REAL a(jmax,kmax,3,3), b(jmax,kmax,3,3), c(jmax,kmax,3,3)
c      integer ithread_num
c      integer function omp_get_thread_num()
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' At begin of btlu3j: b(2,1,2,2) = ',b(2,1,2,2)
c       stop 'stop: at begin of btlu3j'
      endif 
c-----
c  j=jbeg: Decompose b and then compute c
c-----
      if (DEBUG.and.nodeid.eq.0) then
c       print*, ' node 0:  ** before do 10'
c       print*, 'array equivalence:  a = bjm, b=btc, c=bjp'
c       print*, ' ----------'
c       print*, ' node 0:  ** a(1,1,1,1) = ',a(1,1,1,1)
c       print*, ' node 0:  ** b(1,1,1,1) = ',b(1,1,1,1)
c       print*, ' node 0:  ** c(1,1,1,1) = ',c(1,1,1,1)
c       print*, ' ----------'
       print*, ' node 0:  ** a(5,1,1,1) = ',a(5,1,1,1)
       print*, ' node 0:  ** b(5,1,1,1) = ',b(5,1,1,1)
       print*, ' node 0:  ** c(5,1,1,1) = ',c(5,1,1,1)
c       print*, ' ---------'
c       print*, ' node 0:  ** a(1,1,2,2) = ',a(1,1,2,2)
c       print*, ' node 0:  ** b(1,1,2,2) = ',b(1,1,2,2)
c       print*, ' node 0:  ** c(1,1,2,2) = ',c(1,1,2,2)
c       print*, ' ---------'
c       print*, ' node 0:  ** a(1,2,2,2) = ',a(1,2,2,2)
c       print*, ' node 0:  ** b(1,2,2,2) = ',b(1,2,2,2)
c       print*, ' node 0:  ** c(1,2,2,2) = ',c(1,2,2,2)
c       print*, ' ---------'
c       print*, ' node 0:  ** a(2,2,1,1) = ',a(2,2,1,1)
c       print*, ' node 0:  ** b(2,2,1,1) = ',b(2,2,1,1)
c       print*, ' node 0:  ** c(2,2,1,1) = ',c(2,2,1,1)
c       print*, ' ---------'
c       print*, ' node 0:  ** a(2,2,2,2) = ',a(2,2,2,2)
c       print*, ' node 0:  ** b(2,2,2,2) = ',b(2,2,2,2)
c       print*, ' node 0:  ** c(2,2,2,2) = ',c(2,2,2,2)
c       print*, ' ---------'
c       print*
c       print*, ' node 0: ** before do 10 in btlu3.f/btlu3'
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
       print*, 'nodeid  j  k   b(j,k,2,2) b(j,k,1,1) BEFORE 10 continue'
       do 2000 k=1,kmaxx
       do 2000 j=1,jmaxx
        write(*,2010) nodeid,j,k,b(j,k,2,2),b(j,k,1,1)
2010    format(3i5,1p2e13.5)
2000  continue
c       stop 'stop: in btlu3.f/btlu3j 10'
      endif


      j = jbeg

!$omp parallel do 
!$omp+ shared(j,jbeg,jend,kbeg,kend) 
!$omp+ private(d1,d2,d3) 

c      print*, ' jbeg,jend,kbeg,kend = ',jbeg,jend,kbeg,kend
c      print*, ' j = jbeg = ',j
c      stop 'stop before do 10'

      do 10 k=kbeg,kend
         b(j,k,1,1) = 1./b(j,k,1,1)
         b(j,k,1,2) = b(j,k,1,2)*b(j,k,1,1)
         b(j,k,1,3) = b(j,k,1,3)*b(j,k,1,1)
         b(j,k,2,2) = 1./(b(j,k,2,2) - b(j,k,2,1)*b(j,k,1,2) )
         b(j,k,2,3) =    (b(j,k,2,3) - b(j,k,2,1)*b(j,k,1,3))*b(j,k,2,2)
         b(j,k,3,2) =     b(j,k,3,2) - b(j,k,3,1)*b(j,k,1,2)
         b(j,k,3,3) = 1./(b(j,k,3,3) - b(j,k,3,1)*b(j,k,1,3)
     &                               - b(j,k,2,3)*b(j,k,3,2) )
c
         d1 = b(j,k,1,1)* c(j,k,1,1)
         d2 = b(j,k,2,2)*(c(j,k,2,1) - b(j,k,2,1)*d1)
         d3 = b(j,k,3,3)*(c(j,k,3,1) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
         c(j,k,3,1) = d3
         c(j,k,2,1) = d2 - b(j,k,2,3)*c(j,k,3,1)
         c(j,k,1,1) = d1 - b(j,k,1,3)*c(j,k,3,1)
     &                   - b(j,k,1,2)*c(j,k,2,1)
c
         d1 = b(j,k,1,1)* c(j,k,1,2)
         d2 = b(j,k,2,2)*(c(j,k,2,2) - b(j,k,2,1)*d1)
         d3 = b(j,k,3,3)*(c(j,k,3,2) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
         c(j,k,3,2) = d3
         c(j,k,2,2) = d2 - b(j,k,2,3)*c(j,k,3,2)
         c(j,k,1,2) = d1 - b(j,k,1,3)*c(j,k,3,2)
     &                   - b(j,k,1,2)*c(j,k,2,2)
c
         d1 = b(j,k,1,1)* c(j,k,1,3)
         d2 = b(j,k,2,2)*(c(j,k,2,3) - b(j,k,2,1)*d1)
         d3 = b(j,k,3,3)*(c(j,k,3,3) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
         c(j,k,3,3) = d3
         c(j,k,2,3) = d2 - b(j,k,2,3)*c(j,k,3,3)
         c(j,k,1,3) = d1 - b(j,k,1,3)*c(j,k,3,3)
     &                   - b(j,k,1,2)*c(j,k,2,3)
10    continue
!$omp end parallel do

      if (DEBUG) then
       print*, ' stop after do 10 in btlu3.f/btlu3j, nodeid =',nodeid
       stop 'stop after do 10 in btlu3.f/btlu3j'
      endif

      if (DEBUG.and.nodeid.eq.0) then
c       print*, ' node 0:  ***** after do 10'
c       print*, ' node 0:  ***** b(1,1,2,2) = ',b(1,2,2,2)
c       print*, ' node 0:  ***** a(1,1,2,2) = ',a(1,2,2,2)
c       print*, ' node 0:  ***** c(1,1,2,2) = ',c(1,2,2,2)
c       print*, ' node 0:  ***** b(1,2,2,2) = ',b(1,2,2,2)
c       print*, ' node 0:  ***** a(1,2,2,2) = ',a(1,2,2,2)
c       print*, ' node 0:  ***** c(1,2,2,2) = ',c(1,2,2,2)
c       print*, ' node 0:  ***** b(2,2,2,2) = ',b(1,2,2,2)
c       print*, ' node 0:  ***** a(2,2,2,2) = ',a(1,2,2,2)
c       print*, ' node 0:  ***** c(2,2,2,2) = ',c(1,2,2,2)
       print*, ' node 0:  ***** a(5,1,2,2) = ',a(5,1,2,2)
       print*, ' node 0:  ***** b(5,1,2,2) = ',b(5,1,2,2)
       print*, ' node 0:  ***** c(5,1,2,2) = ',c(5,1,2,2)
       print*, ' nodeid = ',nodeid
       print*
       print*, ' nodeid  j  k    b(j,k,2,2) AFTER 10 continue'
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
       do 3000 k=1,kmaxx
       do 3000 j=1,jmaxx
        write(*,2010) nodeid,j,k,b(j,k,2,2)
3000  continue
       call flush(6)
       call flush(istdout)
c       stop 'stop: in btlu3.f/btlu3j 20'
      endif

c-----
c  Interior points
c-----

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, '----------------------------'
       print*, 'Before do 50: jbeg+1, jend-1,kbeg,kend = ',
     &  jbeg+1,jend-1,kbeg,kend
      endif


      do 50 j=jbeg+1,jend-1

c-----
c  Form new b
c-----

!$omp parallel do
!$omp+ shared(kbeg,kend)

         do 20 k=kbeg,kend
            b(j,k,1,1) = b(j,k,1,1) - a(j,k,1,1)*c(j-1,k,1,1)
     &                              - a(j,k,1,2)*c(j-1,k,2,1)
     &                              - a(j,k,1,3)*c(j-1,k,3,1)
            b(j,k,2,1) = b(j,k,2,1) - a(j,k,2,1)*c(j-1,k,1,1)
     &                              - a(j,k,2,2)*c(j-1,k,2,1)
     &                              - a(j,k,2,3)*c(j-1,k,3,1)
            b(j,k,3,1) = b(j,k,3,1) - a(j,k,3,1)*c(j-1,k,1,1)
     &                              - a(j,k,3,2)*c(j-1,k,2,1)
     &                              - a(j,k,3,3)*c(j-1,k,3,1)
            b(j,k,1,2) = b(j,k,1,2) - a(j,k,1,1)*c(j-1,k,1,2)
     &                              - a(j,k,1,2)*c(j-1,k,2,2)
     &                              - a(j,k,1,3)*c(j-1,k,3,2)
            b(j,k,2,2) = b(j,k,2,2) - a(j,k,2,1)*c(j-1,k,1,2)
     &                              - a(j,k,2,2)*c(j-1,k,2,2)
     &                              - a(j,k,2,3)*c(j-1,k,3,2)
            b(j,k,3,2) = b(j,k,3,2) - a(j,k,3,1)*c(j-1,k,1,2)
     &                              - a(j,k,3,2)*c(j-1,k,2,2)
     &                              - a(j,k,3,3)*c(j-1,k,3,2)
            b(j,k,1,3) = b(j,k,1,3) - a(j,k,1,1)*c(j-1,k,1,3)
     &                              - a(j,k,1,2)*c(j-1,k,2,3)
     &                              - a(j,k,1,3)*c(j-1,k,3,3)
            b(j,k,2,3) = b(j,k,2,3) - a(j,k,2,1)*c(j-1,k,1,3)
     &                              - a(j,k,2,2)*c(j-1,k,2,3)
     &                              - a(j,k,2,3)*c(j-1,k,3,3)
            b(j,k,3,3) = b(j,k,3,3) - a(j,k,3,1)*c(j-1,k,1,3)
     &                              - a(j,k,3,2)*c(j-1,k,2,3)
     &                              - a(j,k,3,3)*c(j-1,k,3,3)
20       continue

!$omp end parallel do

c-----
c  lu decomposition of new b and form new c
c-----

!$omp parallel do
!$omp+ shared(kbeg,kend)
!$omp+ private(d1,d2,d3)
         do 30 k=kbeg,kend
            b(j,k,1,1) = 1./b(j,k,1,1)
            b(j,k,1,2) = b(j,k,1,2)*b(j,k,1,1)
            b(j,k,1,3) = b(j,k,1,3)*b(j,k,1,1)
            b(j,k,2,2) = 1./( b(j,k,2,2) - b(j,k,2,1)*b(j,k,1,2) )

            b(j,k,2,3) =(b(j,k,2,3) - b(j,k,2,1)*b(j,k,1,3))*b(j,k,2,2)
            b(j,k,3,2) = b(j,k,3,2) - b(j,k,3,1)*b(j,k,1,2)
            b(j,k,3,3) = 1./( b(j,k,3,3) - b(j,k,3,1)*b(j,k,1,3)
     &                                   - b(j,k,2,3)*b(j,k,3,2) )
c
            d1 = b(j,k,1,1)* c(j,k,1,1)
            d2 = b(j,k,2,2)*(c(j,k,2,1) - b(j,k,2,1)*d1)
            d3 = b(j,k,3,3)*(c(j,k,3,1) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
            c(j,k,3,1) = d3
            c(j,k,2,1) = d2 - b(j,k,2,3)*c(j,k,3,1)
            c(j,k,1,1) = d1 - b(j,k,1,3)*c(j,k,3,1)
     &                      - b(j,k,1,2)*c(j,k,2,1)
c
            d1 = b(j,k,1,1)* c(j,k,1,2)
            d2 = b(j,k,2,2)*(c(j,k,2,2) - b(j,k,2,1)*d1)
            d3 = b(j,k,3,3)*(c(j,k,3,2) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
            c(j,k,3,2) = d3
            c(j,k,2,2) = d2 - b(j,k,2,3)*c(j,k,3,2)
            c(j,k,1,2) = d1 - b(j,k,1,3)*c(j,k,3,2)
     &                      - b(j,k,1,2)*c(j,k,2,2)
c
            d1 = b(j,k,1,1)* c(j,k,1,3)
            d2 = b(j,k,2,2)*(c(j,k,2,3) - b(j,k,2,1)*d1)
            d3 = b(j,k,3,3)*(c(j,k,3,3) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
            c(j,k,3,3) = d3
            c(j,k,2,3) = d2 - b(j,k,2,3)*c(j,k,3,3)
            c(j,k,1,3) = d1 - b(j,k,1,3)*c(j,k,3,3)
     &                      - b(j,k,1,2)*c(j,k,2,3)
30       continue

!$omp end parallel do

50    continue

c-----
c  j=jend: Form new b
c-----

      j = jend

!$omp parallel do
!$omp+ shared(j,kbeg,kend)

      do 60 k=kbeg,kend
         b(j,k,1,1) = b(j,k,1,1) - a(j,k,1,1)*c(j-1,k,1,1)
     &                           - a(j,k,1,2)*c(j-1,k,2,1)
     &                           - a(j,k,1,3)*c(j-1,k,3,1)
         b(j,k,2,1) = b(j,k,2,1) - a(j,k,2,1)*c(j-1,k,1,1)
     &                           - a(j,k,2,2)*c(j-1,k,2,1)
     &                           - a(j,k,2,3)*c(j-1,k,3,1)
         b(j,k,3,1) = b(j,k,3,1) - a(j,k,3,1)*c(j-1,k,1,1)
     &                           - a(j,k,3,2)*c(j-1,k,2,1)
     &                           - a(j,k,3,3)*c(j-1,k,3,1)
         b(j,k,1,2) = b(j,k,1,2) - a(j,k,1,1)*c(j-1,k,1,2)
     &                           - a(j,k,1,2)*c(j-1,k,2,2)
     &                           - a(j,k,1,3)*c(j-1,k,3,2)
         b(j,k,2,2) = b(j,k,2,2) - a(j,k,2,1)*c(j-1,k,1,2)
     &                           - a(j,k,2,2)*c(j-1,k,2,2)
     &                           - a(j,k,2,3)*c(j-1,k,3,2)
         b(j,k,3,2) = b(j,k,3,2) - a(j,k,3,1)*c(j-1,k,1,2)
     &                           - a(j,k,3,2)*c(j-1,k,2,2)
     &                           - a(j,k,3,3)*c(j-1,k,3,2)
         b(j,k,1,3) = b(j,k,1,3) - a(j,k,1,1)*c(j-1,k,1,3)
     &                           - a(j,k,1,2)*c(j-1,k,2,3)
     &                           - a(j,k,1,3)*c(j-1,k,3,3)
         b(j,k,2,3) = b(j,k,2,3) - a(j,k,2,1)*c(j-1,k,1,3)
     &                           - a(j,k,2,2)*c(j-1,k,2,3)
     &                           - a(j,k,2,3)*c(j-1,k,3,3)
         b(j,k,3,3) = b(j,k,3,3) - a(j,k,3,1)*c(j-1,k,1,3)
     &                           - a(j,k,3,2)*c(j-1,k,2,3)
     &                           - a(j,k,3,3)*c(j-1,k,3,3)
60       continue

!$omp end parallel do


c-----
c  lu decomposition of new b
c-----

!$omp parallel do shared(kbeg,kend) 

      do 70 k=kbeg,kend
         b(j,k,1,1) = 1./b(j,k,1,1)
         b(j,k,1,2) = b(j,k,1,2)*b(j,k,1,1)
         b(j,k,1,3) = b(j,k,1,3)*b(j,k,1,1)
         b(j,k,2,2) = 1./(b(j,k,2,2) - b(j,k,2,1)*b(j,k,1,2) )
         b(j,k,2,3) =    (b(j,k,2,3) - b(j,k,2,1)*b(j,k,1,3))*b(j,k,2,2)
         b(j,k,3,2) =     b(j,k,3,2) - b(j,k,3,1)*b(j,k,1,2)
         b(j,k,3,3) = 1./(b(j,k,3,3) - b(j,k,3,1)*b(j,k,1,3)
     &                               - b(j,k,2,3)*b(j,k,3,2) )
70    continue
!$omp end parallel do


      if (DEBUG.and.nodeid.eq.0) then
       print*, ' node 0:  ***** after do 70'
       print*, ' node 0:  ***** b(2,1,2,2) = ',b(2,1,2,2)
       print*, ' node 0:  ***** b(3,1,2,2) = ',b(3,1,2,2)
       print*, ' node 0:  ***** b(4,1,2,2) = ',b(4,1,2,2)
       print*
c       stop 'stop: in btlu3.f/btlu3j 70'
      endif

c-----
c  End of btlu3j
c-----
      return
      end
c
c
c**************************************************************
      subroutine btlu3k(jmax,kmax,a,b,c,jbeg,jend,kbeg,kend)
c**************************************************************
c   LU decomposition of block tridiagonal system of equations
c   Upper blocks of upper blocks are overwritten in c.
c   LU decompostion of each diagonal block of lower matrix are
c   overwritten in b, and the diagonal entries of these blocks
c   are stored as reciprocals.
c
c   block size = (3x3)
c
c------------------------------------------------------------------
#include "precis.h"
      __REAL a(jmax,kmax,3,3), b(jmax,kmax,3,3), c(jmax,kmax,3,3)
c-----
c  k=kbeg: Decompose b and then compute c
c-----

      k = kbeg

!$omp parallel do shared(k,jbeg,jend,kbeg,kend) private(d1,d2,d3)

c      print*
c      print*, ' jbeg,jend,kbeg,kend = ',jbeg,jend,kbeg,kend
c      print*, ' k (which should be = kbeg) = ',k
c      print*
c      stop ' stop before do 10 in btlu3k'


      do 10 j=jbeg,jend
         b(j,k,1,1) = 1./b(j,k,1,1)
         b(j,k,1,2) = b(j,k,1,2)*b(j,k,1,1)
         b(j,k,1,3) = b(j,k,1,3)*b(j,k,1,1)
         b(j,k,2,2) = 1./(b(j,k,2,2) - b(j,k,2,1)*b(j,k,1,2) )
         b(j,k,2,3) =    (b(j,k,2,3) - b(j,k,2,1)*b(j,k,1,3))*b(j,k,2,2)
         b(j,k,3,2) =     b(j,k,3,2) - b(j,k,3,1)*b(j,k,1,2)
         b(j,k,3,3) = 1./(b(j,k,3,3) - b(j,k,3,1)*b(j,k,1,3)
     &                               - b(j,k,2,3)*b(j,k,3,2) )
c
         d1 = b(j,k,1,1)* c(j,k,1,1)
         d2 = b(j,k,2,2)*(c(j,k,2,1) - b(j,k,2,1)*d1)
         d3 = b(j,k,3,3)*(c(j,k,3,1) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
         c(j,k,3,1) = d3
         c(j,k,2,1) = d2 - b(j,k,2,3)*c(j,k,3,1)
         c(j,k,1,1) = d1 - b(j,k,1,3)*c(j,k,3,1)
     &                   - b(j,k,1,2)*c(j,k,2,1)
c
         d1 = b(j,k,1,1)* c(j,k,1,2)
         d2 = b(j,k,2,2)*(c(j,k,2,2) - b(j,k,2,1)*d1)
         d3 = b(j,k,3,3)*(c(j,k,3,2) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
         c(j,k,3,2) = d3
         c(j,k,2,2) = d2 - b(j,k,2,3)*c(j,k,3,2)
         c(j,k,1,2) = d1 - b(j,k,1,3)*c(j,k,3,2)
     &                   - b(j,k,1,2)*c(j,k,2,2)
c
         d1 = b(j,k,1,1)* c(j,k,1,3)
         d2 = b(j,k,2,2)*(c(j,k,2,3) - b(j,k,2,1)*d1)
         d3 = b(j,k,3,3)*(c(j,k,3,3) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
         c(j,k,3,3) = d3
         c(j,k,2,3) = d2 - b(j,k,2,3)*c(j,k,3,3)
         c(j,k,1,3) = d1 - b(j,k,1,3)*c(j,k,3,3)
     &                   - b(j,k,1,2)*c(j,k,2,3)
10    continue

!$omp end parallel do

c-----
c  Interior points
c-----

cdwb-!$omp parallel do shared(jbeg,jend,kbeg,kend) private(d1,d2,d3) 

      do 50 k=kbeg+1,kend-1
c-----
c  Form new b
c-----

!$omp parallel do
!$omp+ shared(jbeg,jend)

         do 20 j=jbeg,jend
            b(j,k,1,1) = b(j,k,1,1) - a(j,k,1,1)*c(j,k-1,1,1)
     &                              - a(j,k,1,2)*c(j,k-1,2,1)
     &                              - a(j,k,1,3)*c(j,k-1,3,1)
            b(j,k,2,1) = b(j,k,2,1) - a(j,k,2,1)*c(j,k-1,1,1)
     &                              - a(j,k,2,2)*c(j,k-1,2,1)
     &                              - a(j,k,2,3)*c(j,k-1,3,1)
            b(j,k,3,1) = b(j,k,3,1) - a(j,k,3,1)*c(j,k-1,1,1)
     &                              - a(j,k,3,2)*c(j,k-1,2,1)
     &                              - a(j,k,3,3)*c(j,k-1,3,1)
            b(j,k,1,2) = b(j,k,1,2) - a(j,k,1,1)*c(j,k-1,1,2)
     &                              - a(j,k,1,2)*c(j,k-1,2,2)
     &                              - a(j,k,1,3)*c(j,k-1,3,2)
            b(j,k,2,2) = b(j,k,2,2) - a(j,k,2,1)*c(j,k-1,1,2)
     &                              - a(j,k,2,2)*c(j,k-1,2,2)
     &                              - a(j,k,2,3)*c(j,k-1,3,2)
            b(j,k,3,2) = b(j,k,3,2) - a(j,k,3,1)*c(j,k-1,1,2)
     &                              - a(j,k,3,2)*c(j,k-1,2,2)
     &                              - a(j,k,3,3)*c(j,k-1,3,2)
            b(j,k,1,3) = b(j,k,1,3) - a(j,k,1,1)*c(j,k-1,1,3)
     &                              - a(j,k,1,2)*c(j,k-1,2,3)
     &                              - a(j,k,1,3)*c(j,k-1,3,3)
            b(j,k,2,3) = b(j,k,2,3) - a(j,k,2,1)*c(j,k-1,1,3)
     &                              - a(j,k,2,2)*c(j,k-1,2,3)
     &                              - a(j,k,2,3)*c(j,k-1,3,3)
            b(j,k,3,3) = b(j,k,3,3) - a(j,k,3,1)*c(j,k-1,1,3)
     &                              - a(j,k,3,2)*c(j,k-1,2,3)
     &                              - a(j,k,3,3)*c(j,k-1,3,3)
20       continue

c-----
c  lu decomposition of new b and form new c
c-----

!$omp parallel do 
!$omp+ shared(jbeg,jend)
!$omp+ private(d1,d2,d3)

         do 30 j=jbeg,jend
            b(j,k,1,1) = 1./b(j,k,1,1)
            b(j,k,1,2) = b(j,k,1,2)*b(j,k,1,1)
            b(j,k,1,3) = b(j,k,1,3)*b(j,k,1,1)
            b(j,k,2,2) = 1./( b(j,k,2,2) - b(j,k,2,1)*b(j,k,1,2) )
            b(j,k,2,3) =(b(j,k,2,3) - b(j,k,2,1)*b(j,k,1,3))*b(j,k,2,2)
            b(j,k,3,2) = b(j,k,3,2) - b(j,k,3,1)*b(j,k,1,2)
            b(j,k,3,3) = 1./( b(j,k,3,3) - b(j,k,3,1)*b(j,k,1,3)
     &                                   - b(j,k,2,3)*b(j,k,3,2) )
c
            d1 = b(j,k,1,1)* c(j,k,1,1)
            d2 = b(j,k,2,2)*(c(j,k,2,1) - b(j,k,2,1)*d1)
            d3 = b(j,k,3,3)*(c(j,k,3,1) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
            c(j,k,3,1) = d3
            c(j,k,2,1) = d2 - b(j,k,2,3)*c(j,k,3,1)
            c(j,k,1,1) = d1 - b(j,k,1,3)*c(j,k,3,1)
     &                      - b(j,k,1,2)*c(j,k,2,1)
c
            d1 = b(j,k,1,1)* c(j,k,1,2)
            d2 = b(j,k,2,2)*(c(j,k,2,2) - b(j,k,2,1)*d1)
            d3 = b(j,k,3,3)*(c(j,k,3,2) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
            c(j,k,3,2) = d3
            c(j,k,2,2) = d2 - b(j,k,2,3)*c(j,k,3,2)
            c(j,k,1,2) = d1 - b(j,k,1,3)*c(j,k,3,2)
     &                      - b(j,k,1,2)*c(j,k,2,2)
c
            d1 = b(j,k,1,1)* c(j,k,1,3)
            d2 = b(j,k,2,2)*(c(j,k,2,3) - b(j,k,2,1)*d1)
            d3 = b(j,k,3,3)*(c(j,k,3,3) - b(j,k,3,1)*d1 - b(j,k,3,2)*d2)
            c(j,k,3,3) = d3
            c(j,k,2,3) = d2 - b(j,k,2,3)*c(j,k,3,3)
            c(j,k,1,3) = d1 - b(j,k,1,3)*c(j,k,3,3)
     &                      - b(j,k,1,2)*c(j,k,2,3)
30       continue

!$omp end parallel do

50    continue

cdwb-!$omp end parallel do

c-----
c  k=kend: Form new b
c-----

      k = kend

!$omp parallel do 
!$omp+ shared(k,jbeg,jend) 

      do 60 j=jbeg,jend
         b(j,k,1,1) = b(j,k,1,1) - a(j,k,1,1)*c(j,k-1,1,1)
     &                           - a(j,k,1,2)*c(j,k-1,2,1)
     &                           - a(j,k,1,3)*c(j,k-1,3,1)
         b(j,k,2,1) = b(j,k,2,1) - a(j,k,2,1)*c(j,k-1,1,1)
     &                           - a(j,k,2,2)*c(j,k-1,2,1)
     &                           - a(j,k,2,3)*c(j,k-1,3,1)
         b(j,k,3,1) = b(j,k,3,1) - a(j,k,3,1)*c(j,k-1,1,1)
     &                           - a(j,k,3,2)*c(j,k-1,2,1)
     &                           - a(j,k,3,3)*c(j,k-1,3,1)
         b(j,k,1,2) = b(j,k,1,2) - a(j,k,1,1)*c(j,k-1,1,2)
     &                           - a(j,k,1,2)*c(j,k-1,2,2)
     &                           - a(j,k,1,3)*c(j,k-1,3,2)
         b(j,k,2,2) = b(j,k,2,2) - a(j,k,2,1)*c(j,k-1,1,2)
     &                           - a(j,k,2,2)*c(j,k-1,2,2)
     &                           - a(j,k,2,3)*c(j,k-1,3,2)
         b(j,k,3,2) = b(j,k,3,2) - a(j,k,3,1)*c(j,k-1,1,2)
     &                           - a(j,k,3,2)*c(j,k-1,2,2)
     &                           - a(j,k,3,3)*c(j,k-1,3,2)
         b(j,k,1,3) = b(j,k,1,3) - a(j,k,1,1)*c(j,k-1,1,3)
     &                           - a(j,k,1,2)*c(j,k-1,2,3)
     &                           - a(j,k,1,3)*c(j,k-1,3,3)
         b(j,k,2,3) = b(j,k,2,3) - a(j,k,2,1)*c(j,k-1,1,3)
     &                           - a(j,k,2,2)*c(j,k-1,2,3)
     &                           - a(j,k,2,3)*c(j,k-1,3,3)
         b(j,k,3,3) = b(j,k,3,3) - a(j,k,3,1)*c(j,k-1,1,3)
     &                           - a(j,k,3,2)*c(j,k-1,2,3)
     &                           - a(j,k,3,3)*c(j,k-1,3,3)
60       continue

!$omp end parallel do

c-----
c  lu decomposition of new b
c-----

!$omp parallel do shared(jbeg,jend)

      do 70 j=jbeg,jend
         b(j,k,1,1) = 1./b(j,k,1,1)
         b(j,k,1,2) = b(j,k,1,2)*b(j,k,1,1)
         b(j,k,1,3) = b(j,k,1,3)*b(j,k,1,1)
         b(j,k,2,2) = 1./(b(j,k,2,2) - b(j,k,2,1)*b(j,k,1,2) )
         b(j,k,2,3) =    (b(j,k,2,3) - b(j,k,2,1)*b(j,k,1,3))*b(j,k,2,2)
         b(j,k,3,2) =     b(j,k,3,2) - b(j,k,3,1)*b(j,k,1,2)
         b(j,k,3,3) = 1./(b(j,k,3,3) - b(j,k,3,1)*b(j,k,1,3)
     &                               - b(j,k,2,3)*b(j,k,3,2) )
70    continue

!$omp end parallel do

c-----
c  End of btlu3k
c-----
      return
      end
c
c
c****************************************************************
      subroutine btund3j(jmax,kmax,a,b,c,jbeg,jend,kbeg,kend)
c****************************************************************
c   Multiplies the lu decompostion of a tridiagonal system together
c   to undecompose the matricies.  This undoes everything done by
c   btlu3j.
c
c   block size = (3x3)
c
c------------------------------------------------------------------
#include "precis.h"
      __REAL a(jmax,kmax,3,3), b(jmax,kmax,3,3), c(jmax,kmax,3,3)
c-----
c  Undecompose each b block
c-----

!$omp parallel do shared(jbeg,jend,kbeg,kend)

      do 12 j=jbeg,jend

       do 10 k=kbeg,kend
         b(j,k,3,3) = 1./b(j,k,3,3) + b(j,k,3,1)*b(j,k,1,3)
     &                              + b(j,k,3,2)*b(j,k,2,3)
         b(j,k,3,2) = b(j,k,3,2) + b(j,k,3,1)*b(j,k,1,2)
         b(j,k,2,2) = 1./b(j,k,2,2)
         b(j,k,2,3) = b(j,k,2,3)*b(j,k,2,2) + b(j,k,2,1)*b(j,k,1,3)
         b(j,k,2,2) = b(j,k,2,2) + b(j,k,2,1)*b(j,k,1,2)
         b(j,k,1,1) = 1./b(j,k,1,1)
         b(j,k,1,3) = b(j,k,1,3)*b(j,k,1,1)
         b(j,k,1,2) = b(j,k,1,2)*b(j,k,1,1)
10     continue
12    continue

!$omp end parallel do

c-----
c  j = jend
c-----

      j = jend

!$omp parallel do 
!$omp+ shared(j,kbeg,kend) private(d1,d2,d3)

      do 20 k=kbeg,kend
         b(j,k,1,1) = b(j,k,1,1) + a(j,k,1,1)*c(j-1,k,1,1)
     &                           + a(j,k,1,2)*c(j-1,k,2,1)
     &                           + a(j,k,1,3)*c(j-1,k,3,1)
         b(j,k,2,1) = b(j,k,2,1) + a(j,k,2,1)*c(j-1,k,1,1)
     &                           + a(j,k,2,2)*c(j-1,k,2,1)
     &                           + a(j,k,2,3)*c(j-1,k,3,1)
         b(j,k,3,1) = b(j,k,3,1) + a(j,k,3,1)*c(j-1,k,1,1)
     &                           + a(j,k,3,2)*c(j-1,k,2,1)
     &                           + a(j,k,3,3)*c(j-1,k,3,1)
         b(j,k,1,2) = b(j,k,1,2) + a(j,k,1,1)*c(j-1,k,1,2)
     &                           + a(j,k,1,2)*c(j-1,k,2,2)
     &                           + a(j,k,1,3)*c(j-1,k,3,2)
         b(j,k,2,2) = b(j,k,2,2) + a(j,k,2,1)*c(j-1,k,1,2)
     &                           + a(j,k,2,2)*c(j-1,k,2,2)
     &                           + a(j,k,2,3)*c(j-1,k,3,2)
         b(j,k,3,2) = b(j,k,3,2) + a(j,k,3,1)*c(j-1,k,1,2)
     &                           + a(j,k,3,2)*c(j-1,k,2,2)
     &                           + a(j,k,3,3)*c(j-1,k,3,2)
         b(j,k,1,3) = b(j,k,1,3) + a(j,k,1,1)*c(j-1,k,1,3)
     &                           + a(j,k,1,2)*c(j-1,k,2,3)
     &                           + a(j,k,1,3)*c(j-1,k,3,3)
         b(j,k,2,3) = b(j,k,2,3) + a(j,k,2,1)*c(j-1,k,1,3)
     &                           + a(j,k,2,2)*c(j-1,k,2,3)
     &                           + a(j,k,2,3)*c(j-1,k,3,3)
         b(j,k,3,3) = b(j,k,3,3) + a(j,k,3,1)*c(j-1,k,1,3)
     &                           + a(j,k,3,2)*c(j-1,k,2,3)
     &                           + a(j,k,3,3)*c(j-1,k,3,3)
20    continue

!$omp end parallel do

c-----
c  Interior points
c------

      do 50 j=jend-1,jbeg+1,-1

!$omp parallel do shared(kbeg,kend) 
!$omp+   private(c11,c21,c12,c22,c13,c23) 

         do 30 k=kbeg,kend
            c11        = b(j,k,1,1)*c(j,k,1,1)
     &                 + b(j,k,1,2)*c(j,k,2,1)
     &                 + b(j,k,1,3)*c(j,k,3,1)
            c21        = b(j,k,2,1)*c(j,k,1,1)
     &                 + b(j,k,2,2)*c(j,k,2,1)
     &                 + b(j,k,2,3)*c(j,k,3,1)
            c(j,k,3,1) = b(j,k,3,1)*c(j,k,1,1)
     &                 + b(j,k,3,2)*c(j,k,2,1)
     &                 + b(j,k,3,3)*c(j,k,3,1)
            c12        = b(j,k,1,1)*c(j,k,1,2)
     &                 + b(j,k,1,2)*c(j,k,2,2)
     &                 + b(j,k,1,3)*c(j,k,3,2)
            c22        = b(j,k,2,1)*c(j,k,1,2)
     &                 + b(j,k,2,2)*c(j,k,2,2)
     &                 + b(j,k,2,3)*c(j,k,3,2)
            c(j,k,3,2) = b(j,k,3,1)*c(j,k,1,2)
     &                 + b(j,k,3,2)*c(j,k,2,2)
     &                 + b(j,k,3,3)*c(j,k,3,2)
            c13        = b(j,k,1,1)*c(j,k,1,3)
     &                 + b(j,k,1,2)*c(j,k,2,3)
     &                 + b(j,k,1,3)*c(j,k,3,3)
            c23        = b(j,k,2,1)*c(j,k,1,3)
     &                 + b(j,k,2,2)*c(j,k,2,3)
     &                 + b(j,k,2,3)*c(j,k,3,3)
            c(j,k,3,3) = b(j,k,3,1)*c(j,k,1,3)
     &                 + b(j,k,3,2)*c(j,k,2,3)
     &                 + b(j,k,3,3)*c(j,k,3,3)
            c(j,k,1,1) = c11
            c(j,k,2,1) = c21
            c(j,k,1,2) = c12
            c(j,k,2,2) = c22
            c(j,k,1,3) = c13
            c(j,k,2,3) = c23
30       continue

!$omp end parallel do

!$omp parallel do
!$omp+ shared(kbeg,kend,j)

         do 40 k=kbeg,kend
            b(j,k,1,1) = b(j,k,1,1) + a(j,k,1,1)*c(j-1,k,1,1)
     &                              + a(j,k,1,2)*c(j-1,k,2,1)
     &                              + a(j,k,1,3)*c(j-1,k,3,1)
            b(j,k,2,1) = b(j,k,2,1) + a(j,k,2,1)*c(j-1,k,1,1)
     &                              + a(j,k,2,2)*c(j-1,k,2,1)
     &                              + a(j,k,2,3)*c(j-1,k,3,1)
            b(j,k,3,1) = b(j,k,3,1) + a(j,k,3,1)*c(j-1,k,1,1)
     &                              + a(j,k,3,2)*c(j-1,k,2,1)
     &                              + a(j,k,3,3)*c(j-1,k,3,1)
            b(j,k,1,2) = b(j,k,1,2) + a(j,k,1,1)*c(j-1,k,1,2)
     &                              + a(j,k,1,2)*c(j-1,k,2,2)
     &                              + a(j,k,1,3)*c(j-1,k,3,2)
            b(j,k,2,2) = b(j,k,2,2) + a(j,k,2,1)*c(j-1,k,1,2)
     &                              + a(j,k,2,2)*c(j-1,k,2,2)
     &                              + a(j,k,2,3)*c(j-1,k,3,2)
            b(j,k,3,2) = b(j,k,3,2) + a(j,k,3,1)*c(j-1,k,1,2)
     &                              + a(j,k,3,2)*c(j-1,k,2,2)
     &                              + a(j,k,3,3)*c(j-1,k,3,2)
            b(j,k,1,3) = b(j,k,1,3) + a(j,k,1,1)*c(j-1,k,1,3)
     &                              + a(j,k,1,2)*c(j-1,k,2,3)
     &                              + a(j,k,1,3)*c(j-1,k,3,3)
            b(j,k,2,3) = b(j,k,2,3) + a(j,k,2,1)*c(j-1,k,1,3)
     &                              + a(j,k,2,2)*c(j-1,k,2,3)
     &                              + a(j,k,2,3)*c(j-1,k,3,3)
            b(j,k,3,3) = b(j,k,3,3) + a(j,k,3,1)*c(j-1,k,1,3)
     &                              + a(j,k,3,2)*c(j-1,k,2,3)
     &                              + a(j,k,3,3)*c(j-1,k,3,3)
40       continue

!$omp end parallel do

50    continue

c-----
c  j = jbeg
c-----

      j = jbeg

!$omp parallel do shared(j,kbeg,kend) 
!$omp+   private(c11,c21,c12,c22,c13,c23) 

      do 60 k=kbeg,kend
         c11        = b(j,k,1,1)*c(j,k,1,1)
     &              + b(j,k,1,2)*c(j,k,2,1)
     &              + b(j,k,1,3)*c(j,k,3,1)
         c21        = b(j,k,2,1)*c(j,k,1,1)
     &              + b(j,k,2,2)*c(j,k,2,1)
     &              + b(j,k,2,3)*c(j,k,3,1)
         c(j,k,3,1) = b(j,k,3,1)*c(j,k,1,1)
     &              + b(j,k,3,2)*c(j,k,2,1)
     &              + b(j,k,3,3)*c(j,k,3,1)
         c12        = b(j,k,1,1)*c(j,k,1,2)
     &              + b(j,k,1,2)*c(j,k,2,2)
     &              + b(j,k,1,3)*c(j,k,3,2)
         c22        = b(j,k,2,1)*c(j,k,1,2)
     &              + b(j,k,2,2)*c(j,k,2,2)
     &              + b(j,k,2,3)*c(j,k,3,2)
         c(j,k,3,2) = b(j,k,3,1)*c(j,k,1,2)
     &              + b(j,k,3,2)*c(j,k,2,2)
     &              + b(j,k,3,3)*c(j,k,3,2)
         c13        = b(j,k,1,1)*c(j,k,1,3)
     &              + b(j,k,1,2)*c(j,k,2,3)
     &              + b(j,k,1,3)*c(j,k,3,3)
         c23        = b(j,k,2,1)*c(j,k,1,3)
     &              + b(j,k,2,2)*c(j,k,2,3)
     &              + b(j,k,2,3)*c(j,k,3,3)
         c(j,k,3,3) = b(j,k,3,1)*c(j,k,1,3)
     &              + b(j,k,3,2)*c(j,k,2,3)
     &              + b(j,k,3,3)*c(j,k,3,3)
         c(j,k,1,1) = c11
         c(j,k,2,1) = c21
         c(j,k,1,2) = c12
         c(j,k,2,2) = c22
         c(j,k,1,3) = c13
         c(j,k,2,3) = c23
60    continue

!$omp end parallel do

c-----
c  End of btund3j
c-----
      return
      end
c
c
c****************************************************************
      subroutine btund3k(jmax,kmax,a,b,c,jbeg,jend,kbeg,kend)
c****************************************************************
c   Multiplies the lu decompostion of a tridiagonal system together
c   to undecompose the matricies.  This undoes everything done by
c   btlu3k.
c
c   block size = (3x3)
c
c------------------------------------------------------------------
#include "precis.h"
      __REAL a(jmax,kmax,3,3), b(jmax,kmax,3,3), c(jmax,kmax,3,3)
c-----
c  Undecompose each b block
c-----

!$omp parallel do shared(jbeg,jend,kbeg,kend)

      do 12 j=jbeg,jend
      do 10 k=kbeg,kend
         b(j,k,3,3) = 1./b(j,k,3,3) + b(j,k,3,1)*b(j,k,1,3)
     &                              + b(j,k,3,2)*b(j,k,2,3)
         b(j,k,3,2) = b(j,k,3,2) + b(j,k,3,1)*b(j,k,1,2)
         b(j,k,2,2) = 1./b(j,k,2,2)
         b(j,k,2,3) = b(j,k,2,3)*b(j,k,2,2) + b(j,k,2,1)*b(j,k,1,3)
         b(j,k,2,2) = b(j,k,2,2) + b(j,k,2,1)*b(j,k,1,2)
         b(j,k,1,1) = 1./b(j,k,1,1)
         b(j,k,1,3) = b(j,k,1,3)*b(j,k,1,1)
         b(j,k,1,2) = b(j,k,1,2)*b(j,k,1,1)
10    continue
12    continue

!$omp end parallel do

c-----
c  k = kend
c-----

      k = kend

!$omp parallel do shared(k,jbeg,jend) 

      do 20 j=jbeg,jend
         b(j,k,1,1) = b(j,k,1,1) + a(j,k,1,1)*c(j,k-1,1,1)
     &                           + a(j,k,1,2)*c(j,k-1,2,1)
     &                           + a(j,k,1,3)*c(j,k-1,3,1)
         b(j,k,2,1) = b(j,k,2,1) + a(j,k,2,1)*c(j,k-1,1,1)
     &                           + a(j,k,2,2)*c(j,k-1,2,1)
     &                           + a(j,k,2,3)*c(j,k-1,3,1)
         b(j,k,3,1) = b(j,k,3,1) + a(j,k,3,1)*c(j,k-1,1,1)
     &                           + a(j,k,3,2)*c(j,k-1,2,1)
     &                           + a(j,k,3,3)*c(j,k-1,3,1)
         b(j,k,1,2) = b(j,k,1,2) + a(j,k,1,1)*c(j,k-1,1,2)
     &                           + a(j,k,1,2)*c(j,k-1,2,2)
     &                           + a(j,k,1,3)*c(j,k-1,3,2)
         b(j,k,2,2) = b(j,k,2,2) + a(j,k,2,1)*c(j,k-1,1,2)
     &                           + a(j,k,2,2)*c(j,k-1,2,2)
     &                           + a(j,k,2,3)*c(j,k-1,3,2)
         b(j,k,3,2) = b(j,k,3,2) + a(j,k,3,1)*c(j,k-1,1,2)
     &                           + a(j,k,3,2)*c(j,k-1,2,2)
     &                           + a(j,k,3,3)*c(j,k-1,3,2)
         b(j,k,1,3) = b(j,k,1,3) + a(j,k,1,1)*c(j,k-1,1,3)
     &                           + a(j,k,1,2)*c(j,k-1,2,3)
     &                           + a(j,k,1,3)*c(j,k-1,3,3)
         b(j,k,2,3) = b(j,k,2,3) + a(j,k,2,1)*c(j,k-1,1,3)
     &                           + a(j,k,2,2)*c(j,k-1,2,3)
     &                           + a(j,k,2,3)*c(j,k-1,3,3)
         b(j,k,3,3) = b(j,k,3,3) + a(j,k,3,1)*c(j,k-1,1,3)
     &                           + a(j,k,3,2)*c(j,k-1,2,3)
     &                           + a(j,k,3,3)*c(j,k-1,3,3)
20    continue

!$omp end parallel do

c-----
c  Interior points
c------

      do 50 k=kend-1,kbeg+1,-1

!$omp parallel do 
!$omp+ shared(jbeg,jend)
!$omp+ private(c11,c21,c12,c13,c23)

         do 30 j=jbeg,jend

            c11        = b(j,k,1,1)*c(j,k,1,1)
     &                 + b(j,k,1,2)*c(j,k,2,1)
     &                 + b(j,k,1,3)*c(j,k,3,1)
            c21        = b(j,k,2,1)*c(j,k,1,1)
     &                 + b(j,k,2,2)*c(j,k,2,1)
     &                 + b(j,k,2,3)*c(j,k,3,1)
            c(j,k,3,1) = b(j,k,3,1)*c(j,k,1,1)
     &                 + b(j,k,3,2)*c(j,k,2,1)
     &                 + b(j,k,3,3)*c(j,k,3,1)
            c12        = b(j,k,1,1)*c(j,k,1,2)
     &                 + b(j,k,1,2)*c(j,k,2,2)
     &                 + b(j,k,1,3)*c(j,k,3,2)
            c22        = b(j,k,2,1)*c(j,k,1,2)
     &                 + b(j,k,2,2)*c(j,k,2,2)
     &                 + b(j,k,2,3)*c(j,k,3,2)
            c(j,k,3,2) = b(j,k,3,1)*c(j,k,1,2)
     &                 + b(j,k,3,2)*c(j,k,2,2)
     &                 + b(j,k,3,3)*c(j,k,3,2)
            c13        = b(j,k,1,1)*c(j,k,1,3)
     &                 + b(j,k,1,2)*c(j,k,2,3)
     &                 + b(j,k,1,3)*c(j,k,3,3)
            c23        = b(j,k,2,1)*c(j,k,1,3)
     &                 + b(j,k,2,2)*c(j,k,2,3)
     &                 + b(j,k,2,3)*c(j,k,3,3)
            c(j,k,3,3) = b(j,k,3,1)*c(j,k,1,3)
     &                 + b(j,k,3,2)*c(j,k,2,3)
     &                 + b(j,k,3,3)*c(j,k,3,3)
            c(j,k,1,1) = c11
            c(j,k,2,1) = c21
            c(j,k,1,2) = c12
            c(j,k,2,2) = c22
            c(j,k,1,3) = c13
            c(j,k,2,3) = c23
30       continue

!$omp end parallel do

!$omp parallel do
!$omp+ shared (jbeg,jend)

         do 40 j=jbeg,jend

            b(j,k,1,1) = b(j,k,1,1) + a(j,k,1,1)*c(j,k-1,1,1)
     &                              + a(j,k,1,2)*c(j,k-1,2,1)
     &                              + a(j,k,1,3)*c(j,k-1,3,1)
            b(j,k,2,1) = b(j,k,2,1) + a(j,k,2,1)*c(j,k-1,1,1)
     &                              + a(j,k,2,2)*c(j,k-1,2,1)
     &                              + a(j,k,2,3)*c(j,k-1,3,1)
            b(j,k,3,1) = b(j,k,3,1) + a(j,k,3,1)*c(j,k-1,1,1)
     &                              + a(j,k,3,2)*c(j,k-1,2,1)
     &                              + a(j,k,3,3)*c(j,k-1,3,1)
            b(j,k,1,2) = b(j,k,1,2) + a(j,k,1,1)*c(j,k-1,1,2)
     &                              + a(j,k,1,2)*c(j,k-1,2,2)
     &                              + a(j,k,1,3)*c(j,k-1,3,2)
            b(j,k,2,2) = b(j,k,2,2) + a(j,k,2,1)*c(j,k-1,1,2)
     &                              + a(j,k,2,2)*c(j,k-1,2,2)
     &                              + a(j,k,2,3)*c(j,k-1,3,2)
            b(j,k,3,2) = b(j,k,3,2) + a(j,k,3,1)*c(j,k-1,1,2)
     &                              + a(j,k,3,2)*c(j,k-1,2,2)
     &                              + a(j,k,3,3)*c(j,k-1,3,2)
            b(j,k,1,3) = b(j,k,1,3) + a(j,k,1,1)*c(j,k-1,1,3)
     &                              + a(j,k,1,2)*c(j,k-1,2,3)
     &                              + a(j,k,1,3)*c(j,k-1,3,3)
            b(j,k,2,3) = b(j,k,2,3) + a(j,k,2,1)*c(j,k-1,1,3)
     &                              + a(j,k,2,2)*c(j,k-1,2,3)
     &                              + a(j,k,2,3)*c(j,k-1,3,3)
            b(j,k,3,3) = b(j,k,3,3) + a(j,k,3,1)*c(j,k-1,1,3)
     &                              + a(j,k,3,2)*c(j,k-1,2,3)
     &                              + a(j,k,3,3)*c(j,k-1,3,3)

40       continue

!$omp end parallel do

50    continue

c-----
c  k = kbeg
c-----

      k = kbeg

!$omp parallel do shared(k,jbeg,jend) 
!$omp+   private(c11,c21,c12,c22,c13,c23) 

      do 60 j=jbeg,jend
         c11        = b(j,k,1,1)*c(j,k,1,1)
     &              + b(j,k,1,2)*c(j,k,2,1)
     &              + b(j,k,1,3)*c(j,k,3,1)
         c21        = b(j,k,2,1)*c(j,k,1,1)
     &              + b(j,k,2,2)*c(j,k,2,1)
     &              + b(j,k,2,3)*c(j,k,3,1)
         c(j,k,3,1) = b(j,k,3,1)*c(j,k,1,1)
     &              + b(j,k,3,2)*c(j,k,2,1)
     &              + b(j,k,3,3)*c(j,k,3,1)
         c12        = b(j,k,1,1)*c(j,k,1,2)
     &              + b(j,k,1,2)*c(j,k,2,2)
     &              + b(j,k,1,3)*c(j,k,3,2)
         c22        = b(j,k,2,1)*c(j,k,1,2)
     &              + b(j,k,2,2)*c(j,k,2,2)
     &              + b(j,k,2,3)*c(j,k,3,2)
         c(j,k,3,2) = b(j,k,3,1)*c(j,k,1,2)
     &              + b(j,k,3,2)*c(j,k,2,2)
     &              + b(j,k,3,3)*c(j,k,3,2)
         c13        = b(j,k,1,1)*c(j,k,1,3)
     &              + b(j,k,1,2)*c(j,k,2,3)
     &              + b(j,k,1,3)*c(j,k,3,3)
         c23        = b(j,k,2,1)*c(j,k,1,3)
     &              + b(j,k,2,2)*c(j,k,2,3)
     &              + b(j,k,2,3)*c(j,k,3,3)
         c(j,k,3,3) = b(j,k,3,1)*c(j,k,1,3)
     &              + b(j,k,3,2)*c(j,k,2,3)
     &              + b(j,k,3,3)*c(j,k,3,3)
         c(j,k,1,1) = c11
         c(j,k,2,1) = c21
         c(j,k,1,2) = c12
         c(j,k,2,2) = c22
         c(j,k,1,3) = c13
         c(j,k,2,3) = c23
60    continue


!$omp end parallel do

c-----
c  End of btund3k
c-----
      return
      end
