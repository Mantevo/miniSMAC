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
c  1. metric
c
c************************************************************
c
c
c******************************************************************
      subroutine metric(jmax,kmax,x,y,rtxy,dj)  !,nz,xj,yj,xk,yk)
c******************************************************************
c   Computes and stores all the metrics.
c   Time varying metrics are set in timemet.
c      rtxy(i,1,1) = d(xi)/d(tau)
c      rtxy(i,1,2) = d(xi)/d(x)
c      rtxy(i,1,3) = d(xi)/d(y)
c      rtxy(i,2,1) = d(eta)/d(tau)
c      rtxy(i,2,2) = d(eta)/d(x)
c      rtxy(i,2,3) = d(eta)/d(y)
c   Metrics are stored in a form which is NOT divided by the Jacobian
c
c------------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)

      __REAL x(jmax,kmax),y(jmax,kmax),
     &       rtxy(jmax,kmax,2,3),dj(jmax,kmax),
     &       xj(jmax,kmax), yj(jmax,kmax),
     &       xk(jmax,kmax), yk(jmax,kmax)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG) then
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif

      if (DEBUG.and.nodeid.eq.1) then
        print*
        print*, ' nodeid  j   k      x(jk)       y(jk)'
        do k=15,17
        do j=30,34
         write(*,910) nodeid,j,k,x(j,k),y(j,k)
910      format(i5,1x,i3,1x,i3,1p2e15.7)
        enddo 
        enddo
c       print*, 'stop: in metric.f/metric 10'
c       stop 'stop: in metric.f/metric 10'
      endif


      if (DEBUG) then
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif

c-----
c  j-direction differences
c-----
      do 10 j=2,jmax-1
      do 10 k=1,kmax
c second-order accurate, central difference
         xj(j,k) = 0.5*( x(j+1,k) - x(j-1,k) )
         yj(j,k) = 0.5*( y(j+1,k) - y(j-1,k) )
10    continue

      j = 1
      do 20 k=1,kmax
c second-order accurate, forward difference
         xj(j,k) = 0.5*( -x(j+2,k) + 4.*x(j+1,k) - 3.*x(j,k) )
         yj(j,k) = 0.5*( -y(j+2,k) + 4.*y(j+1,k) - 3.*y(j,k) )
20    continue

      j = jmax
      do 30 k=1,kmax
c second-order accurate, backward difference
         xj(j,k) = 0.5*( 3.*x(j,k) - 4.*x(j-1,k) + x(j-2,k) )
         yj(j,k) = 0.5*( 3.*y(j,k) - 4.*y(j-1,k) + y(j-2,k) )
30    continue

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' j  k    y_j         yj_jk'
       do 32 k=1,10
       do 32 j=1,17
        write(*,927) j,k,y(j,k),yj(j,k)
927     format(2i3,1p2e15.7)
32     continue
      endif

c-----
c  k-direction differences
c-----
      if (DEBUG.and.nodeid.eq.0) then
       print*
      endif

      do 40 k=2,kmax-1
      do 40 j=1,jmax
         xk(j,k) = 0.5*( x(j,k+1) - x(j,k-1) )
         yk(j,k) = 0.5*( y(j,k+1) - y(j,k-1) )
         if (DEBUG.and.nodeid.eq.0.and.k.eq.2.and.j.le.17) then
          write(*,921) nodeid,j,k,x(j,k+1),x(j,k-1),xk(j,k)
c          write(*,921) nodeid,j,k,y(j,k+1),y(j,k-1),yk(j,k)
921       format('In do 40: ',3i3,3f16.11)
         endif
40    continue

      k = 1
      do 70 j=1,jmax
         xk(j,k) = 0.5*( -x(j,k+2) + 4.*x(j,k+1) - 3.*x(j,k) )
         yk(j,k) = 0.5*( -y(j,k+2) + 4.*y(j,k+1) - 3.*y(j,k) )
70    continue

      k = kmax
      do 80 j=1,jmax
         xk(j,k) = 0.5*( 3.*x(j,k) - 4.*x(j,k-1) + x(j,k-2) )
         yk(j,k) = 0.5*( 3.*y(j,k) - 4.*y(j,k-1) + y(j,k-2) )
80    continue

c-----

c  Compute Jacobian and metrics
c-----
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' nodeid j k      xj      xk      yj      yk     x    y'
      endif

      djmax =-1.e10
      djmin = 1.e10
      do 100 k=1,kmax
      do 100 j=1,jmax
         rdj = xj(j,k)*yk(j,k) - xk(j,k)*yj(j,k)
c-----
c  Debugging the mesh
c-----
c ... make sure mesh jacobian is not zero
        if (rdj.eq.0) then
         write(istdout,*) ' Node ',nodeid,': Zero jacobian: j,k,xj,yk,
     &xk,yj:',
     &                         j,k,xj(j,k),yk(j,k),xk(j,k),yj(j,k)
c         stop ' stopping: zero jacobian in metric.f/metric 15' 
        endif

c ... print selected metrics
c        if (DEBUG.and.nodeid.eq.1.and.j.le.33.and.j.ge.31.and.k.eq.16) 
        if (DEBUG.and.nodeid.eq.0.and.j.le.17.and.k.le.10) 
     &  then
         write(*,951) nodeid,j,k,xj(j,k),xk(j,k),yj(j,k),yk(j,k),x(j,k),
     &     y(j,k) 
         write(*,952) rdj
951      format(3i3,1p6e15.7)
952      format('    rdj = ',e15.7)
        endif
c-----
         dj(j,k) = 1./rdj
         rtxy(j,k,1,2) = yk(j,k)
         rtxy(j,k,1,3) =-xk(j,k)
         rtxy(j,k,2,2) =-yj(j,k)
         rtxy(j,k,2,3) = xj(j,k)
         djmax = max( djmax, dj(j,k) )
         djmin = min( djmin, dj(j,k) )
100   continue

      if (DEBUG) then
      write(*,*) 'In metric.f/metric, Node ',nodeid
      write(*,*) ' Max Jacobian = ',djmax
      write(*,*) ' Min Jacobian = ',djmin
      endif

      if (DEBUG) then
       print*
c       print*, ' stop: after do 100 in metric.f/metric'
c       stop ' stop: after do 100 in metric.f/metric'
      endif

c-----
c  Output metric values to plot3d files
c-----
c      d1 = 1.
c      d2 = 1.
c      d3 = 1.
c      d4 = 1.
c      rewind(33)
c      write(33) jmax,kmax
c      write(33) d1,d2,d3,d4
c      write(33) ((rtxy(j,k,1,2),j=1,jmax),k=1,kmax),
c     &          ((rtxy(j,k,1,3),j=1,jmax),k=1,kmax),
c     &          ((rtxy(j,k,2,2),j=1,jmax),k=1,kmax),
c     &          ((rtxy(j,k,2,3),j=1,jmax),k=1,kmax)
c      close(33)
cc
c      rewind(34)
c      write(34) jmax,kmax
c      write(34) d1,d2,d3,d4
c      write(34) ((dj(j,k),j=1,jmax),k=1,kmax),
c     &          (((1.,j=1,jmax),k=1,kmax),n=1,3)
c      close(34)
c-----
c  Initialize time metrics to zero
c-----
      do 500 j=1,jmax
      do 500 k=1,kmax
         rtxy(j,k,1,1) = 0.0
         rtxy(j,k,2,1) = 0.0
500   continue
 
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' nodeid  j   k  rtxy(jk12)  rtxy(jk13)  rtxy(jk22)  rtxy
     &(jk23)'
       if (jmax.gt.17) then
        jmaxx = 17
       else
        jmaxx = jmax
       endif
       if (kmax.gt.10) then
        kmaxx = 10
       else
        kmaxx = kmax
       endif
       do j=1,jmaxx
       do k=1,kmaxx
        write(*,900) nodeid,j,k,rtxy(j,k,1,2),rtxy(j,k,1,3),
     &  rtxy(j,k,2,2),rtxy(j,k,2,3)
900    format(i5,2i5,1p4e13.5)
       enddo
       enddo
c       stop 'stop: in metric.f after writing metrics'
      endif
       
c-----
c  End of metric
c-----
      return
      end

