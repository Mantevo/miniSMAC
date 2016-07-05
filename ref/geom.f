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
c  1. geominit
c  2. ic
c  3. force
c
c************************************************************
c
c
c*****************************************************************
      subroutine geominit
c*****************************************************************
c  This routines is called once at the beginning of the
c  run and can be used for user-defined input, and to
c  initialize any user-defined common block variables.
c--------------------------------------------------------
#include "common.f"
      common/mulelr/chord,pinf,xrot,yrot,xmom,ymom
      logical DEBUG
      __REAL ivarr(6)
#include "mpif.h"
#include "mpi_params.f"
      INTEGER stat(MPI_STATUS_SIZE)

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Defaults values for namelist variables
c-----
c      alpha = 13.87d0
      chord = 1.0d0
      xrot = 0.5d0
      yrot = 0.0d0
      xmom = 0.25d0
      ymom = 0.0d0

      return

c-----
c  End of geominit
c-----
      return
      end
c
c
c*************************************************************
      subroutine ic(jmax,kmax,q)
c*************************************************************
c  Imposes the default initial conditions: freestream flow at
c  angle of attack, and sets flow on walls to zero.
c  This is called once for each subgrid at the beginning of
c  an initial start.
c-------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      INTEGER nodeid,stat(MPI_STATUS_SIZE)
c
      __REAL q(jmax,kmax,3)

c-----
c   Freestream velocity
c-----
      pi = 4.d0*datan( 1.d0 )
      ca = cos( alpha*pi/180. )
      sa = sin( alpha*pi/180. )
      do 10 j=1,jmax
      do 10 k=1,kmax
         q(j,k,1) = 1.0
         q(j,k,2) = ca
         q(j,k,3) = sa
10    continue

c-----
c  No-slip on all walls
c-----
      do 40 nw=1,nwall
         if(nzwall(nw) .eq. nodeid+1) then
            do 30 j=jwall1(nw),jwall2(nw)
            do 30 k=kwall1(nw),kwall2(nw)
               q(j,k,2) = 0.0 ! u velocity 
               q(j,k,3) = 0.0 ! v velocity
30          continue
         endif
40    continue
c-----
c  End of ic
c-----
      return
      end
c
c
c****************************************************************
      subroutine force(jmax,kmax,x,y,q,rtxy,dj,x0,y0,
     &                 fxp,fyp,fxs,fys,fmom,ioverlap)
c****************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __INTEGER stat(MPI_STATUS_SIZE)
      __REAL x(jmax,kmax),y(jmax,kmax),q(jmax,kmax,3),
     &  rtxy(jmax,kmax,2,3),dj(jmax,kmax)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG) then
       if (nodeid.eq.0) then
        print*
        print*, ' >>> In sub force'
        print*, ' nwall = ',nwall
        print*, '        nw        nwall        jkwall'
        do nw=1,nwall
         print*, nw,nzwall(nw),jkwall(nw)
        enddo
        print*
       endif
c       stop 'stop in geom.f/force' 
      endif

c-----
c  Loop through each of the no-slip walls as specified by the input
c  file bcmain.dat.
c----
      do 100 nw=1,nwall
         if (nodeid+1.eq.nzwall(nw)) then
            if (jkwall(nw).eq.1) then
c jkwall(nw)=1 implies j=constant surface (see rbc.f/rbcmain)
             if (DEBUG) then
              print*
              print*, ' call fmj: nodeid,nw,nzwall(nw),jkwall(nw) = ',
     &          nodeid,nw,nzwall(nw),jkwall(nw)
              print*
              print*, ' ... call fmj in sub force, nodeid = ',nodeid
              print*, ' nsw,k,jkinc,jbeg,jend = ',nw,kwall1(nw),
     &          jkinc(nw),jwall1(nw),jwall2(nw)
              print*, ' nodeid,q(2,1,1),2,3 = ',nodeid,q(2,1,1),
     &          q(2,1,2),q(2,1,3)
             endif
c call fmj in fm.f
             call fmj(jmax,kmax,x,y,q,rtxy,dj,
     &             jwall1(nw),jkinc(nw),kwall1(nw),kwall2(nw),
     &             x0,y0,fxp,fyp,fxs,fys,fmom,nw,ioverlap)

            elseif(jkwall(nw) .eq. 2) then
c jkwall(nw)=2 implies k=constant surface (see rbc.f/rbcmain)
             if (DEBUG) then
              print*
              print*, ' call fmk: nodeid,nw,nzwall(nw),jkwall(nw) = ',
     &          nodeid,nw,nzwall(nw),jkwall(nw)
              print*
              print*, ' ... call fmk in sub force, nodeid = ',nodeid
              print*, ' nw,k,jkinc,jbeg,jend = ',nw,kwall1(nw),
     &           jkinc(nw),jwall1(nw),jwall2(nw)
c              if (nodeid.eq.8) then
               print*,' nodeid, q(2,1,1),2,3 = ',nodeid,q(2,1,1),
     &           q(2,1,2),q(2,1,3)
c              endif
             endif
c call fmk in fm.f
             call fmk(jmax,kmax,x,y,q,rtxy,dj,
     &               kwall1(nw),jwall1(nw),jwall2(nw),
     &               x0,y0,fxp,fyp,fxs,fys,fmom,nw,ioverlap)
            endif
         endif
100   continue
c-----
c  End of force
c-----
      return
      end
