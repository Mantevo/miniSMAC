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
c  1. fluxj1
c  2. fluxk1
c  3. fluxj3
c  4. fluxk3
c  5. fluxj5
c  6. fluxk5 
c
c************************************************************
c
c
c************************************************************
      subroutine fluxj1(jmax,kmax,q,rtxy,dj,jbeg,jend,k,
     & at,am,ap)
c************************************************************
c   Computes the numerical flux at the mid-cell location
c   for a line and loads the values into the
c   array f(j,n), j=jbeg,jend, n=1,3
c   Warning: jend must be less than jmax
c--------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax)
      __REAL at(jmax,3,3),am(jmax,3,3),ap(jmax,3,3)
c-----
c  Flux vectors
c-----
      do 10 j=jbeg,jend+1
         cvel = rtxy(j,k,1,1) + rtxy(j,k,1,2)*q(j,k,2)
     &                        + rtxy(j,k,1,3)*q(j,k,3)
         fv(j,1) = rtxy(j,k,1,1)*q(j,k,1) + beta*(cvel - rtxy(j,k,1,1))
         fv(j,2) = q(j,k,2)*cvel + rtxy(j,k,1,2)*q(j,k,1)
         fv(j,3) = q(j,k,3)*cvel + rtxy(j,k,1,3)*q(j,k,1)
10    continue
c-----
c  Compute delta fluxes
c-----
      j1 = max( jbeg-1, 1 )
      j2 = min( jend+1, jmax-1 )
      do 20 j=j1,j2
         dq1 = q(j+1,k,1) - q(j,k,1)
         dq2 = q(j+1,k,2) - q(j,k,2)
         dq3 = q(j+1,k,3) - q(j,k,3)
         dfp(j,1) = ap(j,1,1)*dq1 + ap(j,1,2)*dq2 + ap(j,1,3)*dq3
         dfp(j,2) = ap(j,2,1)*dq1 + ap(j,2,2)*dq2 + ap(j,2,3)*dq3
         dfp(j,3) = ap(j,3,1)*dq1 + ap(j,3,2)*dq2 + ap(j,3,3)*dq3
         dfm(j,1) = am(j,1,1)*dq1 + am(j,1,2)*dq2 + am(j,1,3)*dq3
         dfm(j,2) = am(j,2,1)*dq1 + am(j,2,2)*dq2 + am(j,2,3)*dq3
         dfm(j,3) = am(j,3,1)*dq1 + am(j,3,2)*dq2 + am(j,3,3)*dq3
20    continue
c-----
c  First order fluxes
c-----
      do 40 j=jbeg,jend
         f(j,1) = 0.5*( fv(j+1,1) + fv(j,1) - dfp(j,1) + dfm(j,1) )
         f(j,2) = 0.5*( fv(j+1,2) + fv(j,2) - dfp(j,2) + dfm(j,2) )
         f(j,3) = 0.5*( fv(j+1,3) + fv(j,3) - dfp(j,3) + dfm(j,3) )
40    continue
c-----
c  End of fluxj1
c-----
      return
      end
c
c
c************************************************************
      subroutine fluxk1(jmax,kmax,q,rtxy,dj,kbeg,kend,j,
     & at,am,ap)
c************************************************************
c   Computes the numerical flux at the mid-cell location
c   for a line and loads the values into the
c   array f(j,n), j=jbeg,jend, n=1,3
c   Warning: jend must be less than jmax
c--------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax)
      __REAL at(kmax,3,3),am(kmax,3,3),ap(kmax,3,3)
c-----
c  Flux vectors
c-----
      do 10 k=kbeg,kend+1
         cvel = rtxy(j,k,2,1) + rtxy(j,k,2,2)*q(j,k,2)
     &                        + rtxy(j,k,2,3)*q(j,k,3)
         fv(k,1) = rtxy(j,k,2,1)*q(j,k,1) + beta*(cvel - rtxy(j,k,2,1))
         fv(k,2) = q(j,k,2)*cvel + rtxy(j,k,2,2)*q(j,k,1)
         fv(k,3) = q(j,k,3)*cvel + rtxy(j,k,2,3)*q(j,k,1)
10    continue
c-----
c  Compute delta fluxes
c-----
      k1 = max( kbeg-1, 1 )
      k2 = min( kend+1, kmax-1 )
      do 20 k=k1,k2
         dq1 = q(j,k+1,1) - q(j,k,1)
         dq2 = q(j,k+1,2) - q(j,k,2)
         dq3 = q(j,k+1,3) - q(j,k,3)
         dfp(k,1) = ap(k,1,1)*dq1 + ap(k,1,2)*dq2 + ap(k,1,3)*dq3
         dfp(k,2) = ap(k,2,1)*dq1 + ap(k,2,2)*dq2 + ap(k,2,3)*dq3
         dfp(k,3) = ap(k,3,1)*dq1 + ap(k,3,2)*dq2 + ap(k,3,3)*dq3
         dfm(k,1) = am(k,1,1)*dq1 + am(k,1,2)*dq2 + am(k,1,3)*dq3
         dfm(k,2) = am(k,2,1)*dq1 + am(k,2,2)*dq2 + am(k,2,3)*dq3
         dfm(k,3) = am(k,3,1)*dq1 + am(k,3,2)*dq2 + am(k,3,3)*dq3
20    continue
c-----
c  First order fluxes
c-----
      do 40 k=kbeg,kend
         f(k,1) = 0.5*( fv(k+1,1) + fv(k,1) - dfp(k,1) + dfm(k,1) )
         f(k,2) = 0.5*( fv(k+1,2) + fv(k,2) - dfp(k,2) + dfm(k,2) )
         f(k,3) = 0.5*( fv(k+1,3) + fv(k,3) - dfp(k,3) + dfm(k,3) )
40    continue
c-----
c  End of fluxk1
c-----
      return
      end
c
c
c************************************************************
      subroutine fluxj3(jmax,kmax,q,rtxy,dj,jbeg,jend,k,
     & at,am,ap)
c************************************************************
c   Computes the numerical flux at the mid-cell location
c   for a line and loads the values into the
c   array f(j,n), j=jbeg,jend, n=1,3
c   Warning: jend must be less than jmax
c--------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)

      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax)
      __REAL at(jmax,3,3),am(jmax,3,3),ap(jmax,3,3)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Flux vectors
c-----
c      do 10 j=jbeg,jend+1
      do 10 j=1,jmax
         cvel = rtxy(j,k,1,1) + rtxy(j,k,1,2)*q(j,k,2)
     &                        + rtxy(j,k,1,3)*q(j,k,3)
         fv(j,1) = rtxy(j,k,1,1)*q(j,k,1) + beta*(cvel - rtxy(j,k,1,1))
         fv(j,2) = q(j,k,2)*cvel + rtxy(j,k,1,2)*q(j,k,1)
         fv(j,3) = q(j,k,3)*cvel + rtxy(j,k,1,3)*q(j,k,1)
10    continue
c-----
c  Compute delta fluxes
c-----
c      j1 = max( jbeg-1, 1 )
c      j2 = min( jend+1, jmax-1 )
c      do 20 j=j1,j2
      do 20 j=1,jmax-1
         dq1 = q(j+1,k,1) - q(j,k,1)
         dq2 = q(j+1,k,2) - q(j,k,2)
         dq3 = q(j+1,k,3) - q(j,k,3)
         dfp(j,1) = ap(j,1,1)*dq1 + ap(j,1,2)*dq2 + ap(j,1,3)*dq3
         dfp(j,2) = ap(j,2,1)*dq1 + ap(j,2,2)*dq2 + ap(j,2,3)*dq3
         dfp(j,3) = ap(j,3,1)*dq1 + ap(j,3,2)*dq2 + ap(j,3,3)*dq3
         dfm(j,1) = am(j,1,1)*dq1 + am(j,1,2)*dq2 + am(j,1,3)*dq3
         dfm(j,2) = am(j,2,1)*dq1 + am(j,2,2)*dq2 + am(j,2,3)*dq3
         dfm(j,3) = am(j,3,1)*dq1 + am(j,3,2)*dq2 + am(j,3,3)*dq3
20    continue
c-----
c  Third order fluxes
c-----
      r6 = 1./6.
      endeps = 0.01
c      do 40 j=jbeg,jend
      do 40 j=1,jmax-1
         jm = max0( 1, j-1 )
         jp = min0( jmax, j+1 )
         jpp = min0( jmax, j+2 )
c-----
c  Third+Second order flux-difference at endpoint
c-----
         fm1 = (dfp(j,1) + dfm(j,1) - dfp(jp,1) - dfm(jp,1))*r6
         fm2 = (dfp(j,2) + dfm(j,2) - dfp(jp,2) - dfm(jp,2))*r6
         fm3 = (dfp(j,3) + dfm(j,3) - dfp(jp,3) - dfm(jp,3))*r6
         fp1 = (dfp(jm,1) + dfm(jm,1) - dfp(j,1) - dfm(j,1))*r6
         fp2 = (dfp(jm,2) + dfm(jm,2) - dfp(j,2) - dfm(j,2))*r6
         fp3 = (dfp(jm,3) + dfm(jm,3) - dfp(j,3) - dfm(j,3))*r6
c-----
c  Second order plus small amount of first order flux-difference at endpoint
c-----
c         fm1 = 0.5*endeps*( -dfp(j,1) + dfm(j,1) )
c         fm2 = 0.5*endeps*( -dfp(j,2) + dfm(j,2) )
c         fm3 = 0.5*endeps*( -dfp(j,3) + dfm(j,3) )
c         fp1 = 0.5*endeps*( -dfp(j,1) + dfm(j,1) )
c         fp2 = 0.5*endeps*( -dfp(j,2) + dfm(j,2) )
c         fp3 = 0.5*endeps*( -dfp(j,3) + dfm(j,3) )
c-----
         f31 = (dfp(jm,1)-dfp(j,1)+dfm(j,1)-dfm(jp,1))*r6
         f32 = (dfp(jm,2)-dfp(j,2)+dfm(j,2)-dfm(jp,2))*r6
         f33 = (dfp(jm,3)-dfp(j,3)+dfm(j,3)-dfm(jp,3))*r6
c
         cp = float(jpp - jp)
         cm = float(j   - jm)
      if (DEBUG.and.nodeid.eq.0) then
c       print*, ' flux.f: jpp,jp,j,jm,cp,cm = ',jpp,jp,j,jm,cp,cm
      endif
c
         f(j,1) = 0.5*( fv(jp,1) + fv(j,1) )
     &          + cdis*(cp*cm*f31 + cm*(1.-cp)*fp1 + cp*(1.-cm)*fm1)
         f(j,2) = 0.5*( fv(jp,2) + fv(j,2) )
     &          + cdis*(cp*cm*f32 + cm*(1.-cp)*fp2 + cp*(1.-cm)*fm2)
         f(j,3) = 0.5*( fv(jp,3) + fv(j,3) )
     &          + cdis*(cp*cm*f33 + cm*(1.-cp)*fp3 + cp*(1.-cm)*fm3)
      if (DEBUG.and.nodeid.eq.0) then
       print*, 'j,fj1,fj2,fj3 = ',j,f(j,1),f(j,2),f(j,3)
      endif
40    continue
 
      if (DEBUG) then
c       stop ' stop: in flux.f'
      endif

c-----
c  End of fluxj3
c-----
      return
      end
c
c
c************************************************************
      subroutine fluxk3(jmax,kmax,q,rtxy,dj,kbeg,kend,j,
     & at,am,ap)
c************************************************************
c   Computes the numerical flux at the mid-cell location
c   for a line and loads the values into the
c   array f(j,n), j=jbeg,jend, n=1,3
c   Warning: jend must be less than jmax
c--------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax)
      __REAL at(kmax,3,3),am(kmax,3,3),ap(kmax,3,3)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG.and.nodeid.eq.0) then
       print*, ' ---------- in fluxk3 ----------'
      endif

c-----
c  Flux vectors
c-----
      do 10 k=kbeg,kend+1
         cvel = rtxy(j,k,2,1) + rtxy(j,k,2,2)*q(j,k,2)
     &                        + rtxy(j,k,2,3)*q(j,k,3)
         fv(k,1) = rtxy(j,k,2,1)*q(j,k,1) + beta*(cvel - rtxy(j,k,2,1))
         fv(k,2) = q(j,k,2)*cvel + rtxy(j,k,2,2)*q(j,k,1)
         fv(k,3) = q(j,k,3)*cvel + rtxy(j,k,2,3)*q(j,k,1)
10    continue
c-----
c  Compute delta fluxes
c-----
      k1 = max( kbeg-1, 1 )
      k2 = min( kend+1, kmax-1 )
      do 20 k=k1,k2
         dq1 = q(j,k+1,1) - q(j,k,1)
         dq2 = q(j,k+1,2) - q(j,k,2)
         dq3 = q(j,k+1,3) - q(j,k,3)
         dfp(k,1) = ap(k,1,1)*dq1 + ap(k,1,2)*dq2 + ap(k,1,3)*dq3
         dfp(k,2) = ap(k,2,1)*dq1 + ap(k,2,2)*dq2 + ap(k,2,3)*dq3
         dfp(k,3) = ap(k,3,1)*dq1 + ap(k,3,2)*dq2 + ap(k,3,3)*dq3
         dfm(k,1) = am(k,1,1)*dq1 + am(k,1,2)*dq2 + am(k,1,3)*dq3
         dfm(k,2) = am(k,2,1)*dq1 + am(k,2,2)*dq2 + am(k,2,3)*dq3
         dfm(k,3) = am(k,3,1)*dq1 + am(k,3,2)*dq2 + am(k,3,3)*dq3
20    continue
c-----
c  Third order fluxes
c-----
      r6 = 1./6.
      endeps = 0.01
      do 40 k=kbeg,kend
         km = max0( 1, k-1 )
         kp = min0( kmax, k+1 )
         kpp = min0( kmax, k+2 )
c-----
c  Third+Second order flux-difference at endpoint
c-----
         fm1 = (dfp(k,1) + dfm(k,1) - dfp(kp,1) - dfm(kp,1))*r6
         fm2 = (dfp(k,2) + dfm(k,2) - dfp(kp,2) - dfm(kp,2))*r6
         fm3 = (dfp(k,3) + dfm(k,3) - dfp(kp,3) - dfm(kp,3))*r6
         fp1 = (dfp(km,1) + dfm(km,1) - dfp(k,1) - dfm(k,1))*r6
         fp2 = (dfp(km,2) + dfm(km,2) - dfp(k,2) - dfm(k,2))*r6
         fp3 = (dfp(km,3) + dfm(km,3) - dfp(k,3) - dfm(k,3))*r6
c-----
c  Second order plus small amount of first order flux-difference at endpoint
c-----
c         fm1 = 0.5*endeps*( -dfp(k,1) + dfm(k,1) )
c         fm2 = 0.5*endeps*( -dfp(k,2) + dfm(k,2) )
c         fm3 = 0.5*endeps*( -dfp(k,3) + dfm(k,3) )
c         fp1 = 0.5*endeps*( -dfp(k,1) + dfm(k,1) )
c         fp2 = 0.5*endeps*( -dfp(k,2) + dfm(k,2) )
c         fp3 = 0.5*endeps*( -dfp(k,3) + dfm(k,3) )
c-----
         f31 = (dfp(km,1)-dfp(k,1)+dfm(k,1)-dfm(kp,1))*r6
         f32 = (dfp(km,2)-dfp(k,2)+dfm(k,2)-dfm(kp,2))*r6
         f33 = (dfp(km,3)-dfp(k,3)+dfm(k,3)-dfm(kp,3))*r6
c
         cp = float(kpp - kp)
         cm = float(k   - km)
c
         f(k,1) = 0.5*( fv(kp,1) + fv(k,1) )
     &          + cdis*(cp*cm*f31 + cm*(1.-cp)*fp1 + cp*(1.-cm)*fm1)
         f(k,2) = 0.5*( fv(kp,2) + fv(k,2) )
     &          + cdis*(cp*cm*f32 + cm*(1.-cp)*fp2 + cp*(1.-cm)*fm2)
         f(k,3) = 0.5*( fv(kp,3) + fv(k,3) )
     &          + cdis*(cp*cm*f33 + cm*(1.-cp)*fp3 + cp*(1.-cm)*fm3)

      if (DEBUG.and.nodeid.eq.0) then
       print*, ' k,fk1,fk2,fk3 = ',k,f(k,1),f(k,2),f(k,3)
      endif

40    continue

      if (DEBUG.and.nodeid.eq.0) then
       stop ' stopping in flux.f/fluxk3'
      endif

c-----
c  End of fluxk3
c-----
      return
      end
c
c
c************************************************************
      subroutine fluxj5(jmax,kmax,q,rtxy,dj,jbeg,jend,k,
     & at,am,ap)
c************************************************************
c   Computes the numerical flux at the mid-cell location
c   for a line and loads the values into the
c   array f(j,n), j=jbeg,jend, n=1,3
c   Warning: jend must be less than jmax
c--------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax)
      __REAL at(jmax,3,3),am(jmax,3,3),ap(jmax,3,3)
c-----
c  Flux vectors
c-----
      do 10 j=jbeg,jend+1
         cvel = rtxy(j,k,1,1) + rtxy(j,k,1,2)*q(j,k,2)
     &                        + rtxy(j,k,1,3)*q(j,k,3)
         fv(j,1) = rtxy(j,k,1,1)*q(j,k,1) + beta*(cvel - rtxy(j,k,1,1))
         fv(j,2) = q(j,k,2)*cvel + rtxy(j,k,1,2)*q(j,k,1)
         fv(j,3) = q(j,k,3)*cvel + rtxy(j,k,1,3)*q(j,k,1)
10    continue
      do 15 j=jbeg,jend
         f(j,1) = 0.5*( fv(j+1,1) + fv(j,1) )
         f(j,2) = 0.5*( fv(j+1,2) + fv(j,2) )
         f(j,3) = 0.5*( fv(j+1,3) + fv(j,3) )
15    continue
c-----
c  Compute delta fluxes
c-----
      j1 = max( jbeg-1, 1 )
      j2 = min( jend+1, jmax-1 )
      do 20 j=j1,j2
         dq1 = q(j+1,k,1) - q(j,k,1)
         dq2 = q(j+1,k,2) - q(j,k,2)
         dq3 = q(j+1,k,3) - q(j,k,3)
         dfp(j,1) = ap(j,1,1)*dq1 + ap(j,1,2)*dq2 + ap(j,1,3)*dq3
         dfp(j,2) = ap(j,2,1)*dq1 + ap(j,2,2)*dq2 + ap(j,2,3)*dq3
         dfp(j,3) = ap(j,3,1)*dq1 + ap(j,3,2)*dq2 + ap(j,3,3)*dq3
         dfm(j,1) = am(j,1,1)*dq1 + am(j,1,2)*dq2 + am(j,1,3)*dq3
         dfm(j,2) = am(j,2,1)*dq1 + am(j,2,2)*dq2 + am(j,2,3)*dq3
         dfm(j,3) = am(j,3,1)*dq1 + am(j,3,2)*dq2 + am(j,3,3)*dq3
20    continue
c-----
c  Fifth order fluxes
c-----
      r6 = 1./6.
      r60 = 1./60.
      do 40 j=jbeg,jend
         jmm = max0( 1, j-2 )
         jm = max0( 1, j-1 )
         jp = min0( jmax, j+1 )
         jpp = min0( jmax, j+2 )
         jppp = min0( jmax, j+3 )
c-----
         cpp= float(jppp- jpp)
         cp = float(jpp -  jp)
         cm = float(j   -  jm)
         cmm= float(jm  - jmm)
         c3 = cp*cm
         c5 = cpp*cp*cm*cmm
c-----
         f51 = r60*(
     &       - 2.*dfp(jmm,1) +11.*dfp(jm,1) -6.*dfp(j,1) -3.*dfp(jp,1)
     &       + 2.*dfm(jpp,1) -11.*dfm(jp,1) +6.*dfm(j,1) +3.*dfm(jm,1))
         f52 = r60*(
     &       - 2.*dfp(jmm,2) +11.*dfp(jm,2) -6.*dfp(j,2) -3.*dfp(jp,2)
     &       + 2.*dfm(jpp,2) -11.*dfm(jp,2) +6.*dfm(j,2) +3.*dfm(jm,2))
         f53 = r60*(
     &       - 2.*dfp(jmm,3) +11.*dfp(jm,3) -6.*dfp(j,3) -3.*dfp(jp,3)
     &       + 2.*dfm(jpp,3) -11.*dfm(jp,3) +6.*dfm(j,3) +3.*dfm(jm,3))
c-----
         f31 = (dfp(jm,1) - dfp(j ,1) + dfm(j ,1) - dfm(jp,1))*r6
         f32 = (dfp(jm,2) - dfp(j ,2) + dfm(j ,2) - dfm(jp,2))*r6
         f33 = (dfp(jm,3) - dfp(j ,3) + dfm(j ,3) - dfm(jp,3))*r6
         fm1 = (dfp(j ,1) - dfp(jp,1) + dfm(j ,1) - dfm(jp,1))*r6
         fm2 = (dfp(j ,2) - dfp(jp,2) + dfm(j ,2) - dfm(jp,2))*r6
         fm3 = (dfp(j ,3) - dfp(jp,3) + dfm(j ,3) - dfm(jp,3))*r6
         fp1 = (dfp(jm,1) - dfp(j ,1) + dfm(jm,1) - dfm(j ,1))*r6
         fp2 = (dfp(jm,2) - dfp(j ,2) + dfm(jm,2) - dfm(j ,2))*r6
         fp3 = (dfp(jm,3) - dfp(j ,3) + dfm(jm,3) - dfm(j ,3))*r6
c-----
         f(j,1) = f(j,1) + cdis*( c5*f51 + cp*cm*(1.-c5)*f31 
     &          + cm*(1.-cp)*(1.-c5)*fp1 + cp*(1.-cm)*(1.-c5)*fm1 ) 
         f(j,2) = f(j,2) + cdis*( c5*f52  + cp*cm*(1.-c5)*f32 
     &          + cm*(1.-cp)*(1.-c5)*fp2 + cp*(1.-cm)*(1.-c5)*fm2 ) 
         f(j,3) = f(j,3) + cdis*( c5*f53 + cp*cm*(1.-c5)*f33 
     &          + cm*(1.-cp)*(1.-c5)*fp3 + cp*(1.-cm)*(1.-c5)*fm3 ) 
c-----
40    continue
c-----
c  End of fluxj5
c-----
      return
      end
c
c
c************************************************************
      subroutine fluxk5(jmax,kmax,q,rtxy,dj,kbeg,kend,j,
     &  at,am,ap)
c************************************************************
c   Computes the numerical flux at the mid-cell location
c   for a line and loads the values into the
c   array f(j,n), j=jbeg,jend, n=1,3
c   Warning: jend must be less than jmax
c--------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax)
      __REAL at(kmax,3,3),am(kmax,3,3),ap(kmax,3,3)
c-----
c  Flux vectors
c-----
      do 10 k=kbeg,kend+1
         cvel = rtxy(j,k,2,1) + rtxy(j,k,2,2)*q(j,k,2)
     &                        + rtxy(j,k,2,3)*q(j,k,3)
         fv(k,1) = rtxy(j,k,2,1)*q(j,k,1) + beta*(cvel - rtxy(j,k,2,1))
         fv(k,2) = q(j,k,2)*cvel + rtxy(j,k,2,2)*q(j,k,1)
         fv(k,3) = q(j,k,3)*cvel + rtxy(j,k,2,3)*q(j,k,1)
10    continue
      do 15 k=kbeg,kend
         f(k,1) = 0.5*( fv(k+1,1) + fv(k,1) )
         f(k,2) = 0.5*( fv(k+1,2) + fv(k,2) )
         f(k,3) = 0.5*( fv(k+1,3) + fv(k,3) )
15    continue
c-----
c  Compute delta fluxes
c-----
      k1 = max( kbeg-1, 1 )
      k2 = min( kend+1, kmax-1 )
      do 20 k=k1,k2
         dq1 = q(j,k+1,1) - q(j,k,1)
         dq2 = q(j,k+1,2) - q(j,k,2)
         dq3 = q(j,k+1,3) - q(j,k,3)
         dfp(k,1) = ap(k,1,1)*dq1 + ap(k,1,2)*dq2 + ap(k,1,3)*dq3
         dfp(k,2) = ap(k,2,1)*dq1 + ap(k,2,2)*dq2 + ap(k,2,3)*dq3
         dfp(k,3) = ap(k,3,1)*dq1 + ap(k,3,2)*dq2 + ap(k,3,3)*dq3
         dfm(k,1) = am(k,1,1)*dq1 + am(k,1,2)*dq2 + am(k,1,3)*dq3
         dfm(k,2) = am(k,2,1)*dq1 + am(k,2,2)*dq2 + am(k,2,3)*dq3
         dfm(k,3) = am(k,3,1)*dq1 + am(k,3,2)*dq2 + am(k,3,3)*dq3
20    continue
c-----
c  Fifth order fluxes
c-----
      r6 = 1./6.
      r60 = 1./60.
      do 40 k=kbeg,kend
         kmm = max0( 1, k-2 )
         km = max0( 1, k-1 )
         kp = min0( kmax, k+1 )
         kpp = min0( kmax, k+2 )
         kppp = min0( kmax, k+3 )
c
         cpp= float(kppp- kpp)
         cp = float(kpp -  kp)
         cm = float(k   -  km)
         cmm= float(km  - kmm)
         c3 = cp*cm
         c5 = cpp*cp*cm*cmm
c-----
         f51 = r60*(
     &       - 2.*dfp(kmm,1) +11.*dfp(km,1) -6.*dfp(k,1) -3.*dfp(kp,1)
     &       + 2.*dfm(kpp,1) -11.*dfm(kp,1) +6.*dfm(k,1) +3.*dfm(km,1))
         f52 = r60*(
     &       - 2.*dfp(kmm,2) +11.*dfp(km,2) -6.*dfp(k,2) -3.*dfp(kp,2)
     &       + 2.*dfm(kpp,2) -11.*dfm(kp,2) +6.*dfm(k,2) +3.*dfm(km,2))
         f53 = r60*(
     &       - 2.*dfp(kmm,3) +11.*dfp(km,3) -6.*dfp(k,3) -3.*dfp(kp,3)
     &       + 2.*dfm(kpp,3) -11.*dfm(kp,3) +6.*dfm(k,3) +3.*dfm(km,3))
c-----
         f31 = (dfp(km,1) - dfp(k ,1) + dfm(k ,1) - dfm(kp,1))*r6
         f32 = (dfp(km,2) - dfp(k ,2) + dfm(k ,2) - dfm(kp,2))*r6
         f33 = (dfp(km,3) - dfp(k ,3) + dfm(k ,3) - dfm(kp,3))*r6
         fm1 = (dfp(k ,1) - dfp(kp,1) + dfm(k ,1) - dfm(kp,1))*r6
         fm2 = (dfp(k ,2) - dfp(kp,2) + dfm(k ,2) - dfm(kp,2))*r6
         fm3 = (dfp(k ,3) - dfp(kp,3) + dfm(k ,3) - dfm(kp,3))*r6
         fp1 = (dfp(km,1) - dfp(k ,1) + dfm(km,1) - dfm(k ,1))*r6
         fp2 = (dfp(km,2) - dfp(k ,2) + dfm(km,2) - dfm(k ,2))*r6
         fp3 = (dfp(km,3) - dfp(k ,3) + dfm(km,3) - dfm(k ,3))*r6
c-----
         f(k,1) = f(k,1) + cdis*( c5*f51 + cp*cm*(1.-c5)*f31 
     &          + cm*(1.-cp)*(1.-c5)*fp1 + cp*(1.-cm)*(1.-c5)*fm1 ) 
         f(k,2) = f(k,2) + cdis*( c5*f52  + cp*cm*(1.-c5)*f32 
     &          + cm*(1.-cp)*(1.-c5)*fp2 + cp*(1.-cm)*(1.-c5)*fm2 ) 
         f(k,3) = f(k,3) + cdis*( c5*f53 + cp*cm*(1.-c5)*f33 
     &          + cm*(1.-cp)*(1.-c5)*fp3 + cp*(1.-cm)*(1.-c5)*fm3 ) 
c-----
40    continue
c-----
c  End of fluxk5
c-----
      return
      end
