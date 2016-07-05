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
c  1. conup 
c
c************************************************************
c
c
c*************************************************************
      subroutine conup(jmax,kmax,q,rtxy,dj,s,
     &             btc,bjm,bjp,bkm,bkp)
c*************************************************************
c
c   Computes the implicit and explicit convective terms using
c   upwind differencing.
c   Uses 1st order differencing in forming the left-hand side terms.
c-------------------------------------------------------------
#include "common.f"
#include "mpi_params.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax),
     &  s(jmax,kmax,3), btc(jmax,kmax,3,3),
     &  bjm(jmax,kmax,3,3), bjp(jmax,kmax,3,3),
     &  bkm(jmax,kmax,3,3), bkp(jmax,kmax,3,3)
      __REAL, allocatable:: at(:,:,:),am(:,:,:),ap(:,:,:)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' --------- in conup.f/conup -----------'
       print*, ' s(17,2,1),f(2,1) = ',s(17,2,1),f(2,1)
c       stop 'stop: AT START OF conup.f/conup 1'
      endif

c-----
c  J-sweep
c-----
      allocate(at(jmax,3,3),am(jmax,3,3),ap(jmax,3,3))

      do 100 k=2,kmax-1

       if (DEBUG.and.nodeid.eq.0.and.k.eq.2) then
        print*
        print*, ' beta,cvel = ',beta,cvel
        print*, '  j   k    atj12       atj22       atj32'
       endif

         do 50 j=1,jmax
            cvel = rtxy(j,k,1,1) + rtxy(j,k,1,2)*q(j,k,2)
     &           + rtxy(j,k,1,3)*q(j,k,3)
            at(j,1,1) = rtxy(j,k,1,1)
            at(j,2,1) = rtxy(j,k,1,2)
            at(j,3,1) = rtxy(j,k,1,3)
            at(j,1,2) = rtxy(j,k,1,2)*beta
            at(j,2,2) = rtxy(j,k,1,2)*q(j,k,2) + cvel
            at(j,3,2) = rtxy(j,k,1,2)*q(j,k,3)
            at(j,1,3) = rtxy(j,k,1,3)*beta
            at(j,2,3) = rtxy(j,k,1,3)*q(j,k,2)
            at(j,3,3) = rtxy(j,k,1,3)*q(j,k,3) + cvel
            if (DEBUG.and.nodeid.eq.0.and.j.le.17.and.k.eq.2) then
             write(*,900) j,k,at(j,1,2),at(j,2,2),at(j,3,2)
900          format(2i4,1p3e13.5)
            endif
50       continue

      if (DEBUG.and.nodeid.eq.0) then
       call flush(6)
       call flush(istdout)
c       print*, 'stop: in conup.f/conup 50'
c       stop 'stop: in conup.f/conup 50'
      endif

c compute plus/minus Jacobian matrices at mid-cell location;
c  form the A- matrix by subtracting A+ from A matrix
c 
c         call amatpmj(jmax,kmax,q,rtxy,dj,1,jmax-1,k)
c  sub amatpmj is coded inline below in the 'do 30' loop
c-----
c  Compute matrix
c-----
      if (DEBUG.and.nodeid.eq.0.and.k.eq.2) then
       print*
       print*, ' j  k    ap_11    ap_22   ap_33   am_11   am_22   am_33'
      endif

      do 30 j=1,jmax-1
         rt = 0.5*( rtxy(j+1,k,1,1) + rtxy(j,k,1,1) )
         rx = 0.5*( rtxy(j+1,k,1,2) + rtxy(j,k,1,2) )
         ry = 0.5*( rtxy(j+1,k,1,3) + rtxy(j,k,1,3) )
         q2 = 0.5*( q(j+1,k,2) + q(j,k,2) )
         q3 = 0.5*( q(j+1,k,3) + q(j,k,3) )
         cvel = rx*q2 + ry*q3
         cc = sqrt( (cvel+0.5*rt)**2 + beta*(rx*rx + ry*ry) )
         eig = cvel + rt
         eig1 = 0.5*(eig + abs(eig))
         eig2 = cvel + rt*0.5 + cc
c-----
c  Load A+ matrix
c-----
         x21 =-ry*eig1/(cc*cc - 0.25*rt*rt)
         x31 = rx*eig1/(cc*cc - 0.25*rt*rt)
         x12 = 0.5*eig2/cc
         x22 = 0.5*eig2*(q2*eig2 + beta*rx)/(beta*cc*(cc-0.5*rt))
         x32 = 0.5*eig2*(q3*eig2 + beta*ry)/(beta*cc*(cc-0.5*rt))
         xi11 = ry*q2 - rx*q3
         xi12 =-q3*(cvel + rt) - ry*beta
         xi13 = q2*(cvel + rt) + rx*beta
         xi21 = cc - cvel - 0.5*rt
         xi22 = rx*beta
         xi23 = ry*beta
c-----
         ap(j,1,1) =            x12*xi21
         ap(j,2,1) = x21*xi11 + x22*xi21
         ap(j,3,1) = x31*xi11 + x32*xi21
         ap(j,1,2) =            x12*xi22
         ap(j,2,2) = x21*xi12 + x22*xi22
         ap(j,3,2) = x31*xi12 + x32*xi22
         ap(j,1,3) =            x12*xi23
         ap(j,2,3) = x21*xi13 + x22*xi23
         ap(j,3,3) = x31*xi13 + x32*xi23
c-----
c  Load A- matrix
c-----
         am(j,1,1) = -ap(j,1,1)
         am(j,2,1) = -ap(j,2,1) + rx
         am(j,3,1) = -ap(j,3,1) + ry
c
         am(j,1,2) = -ap(j,1,2) + rx*beta
         am(j,2,2) = -ap(j,2,2) + rt + 2.*rx*q2 + ry*q3
         am(j,3,2) = -ap(j,3,2) + rx*q3
c
         am(j,1,3) = -ap(j,1,3) + ry*beta
         am(j,2,3) = -ap(j,2,3) + ry*q2
         am(j,3,3) = -ap(j,3,3) + rt + rx*q2 + 2.*ry*q3
      if (DEBUG.and.nodeid.eq.0.and.k.eq.2.and.j.le.17) then
       write(*,913) j,k,ap(j,1,1),ap(j,2,2),ap(j,3,3),am(j,1,1),
     &   am(j,2,2),am(j,3,3)
913    format(2i3,1p6e13.5)
       write(*,914) rt,rx,ry,q2,q3
914    format(5x,1p5e13.5)
       write(*,915) rtxy(j+1,k,1,3),rtxy(j,k,1,3)
915    format(5x,1p2e13.5)
      endif
30    continue

c compute first-order, 3rd-order, or 5th-order flux differencing
c  splitting for the convective terms.
c  Default: iflxo = 3
         if(iflxo .eq. 1) then
            call fluxj1(jmax,kmax,q,rtxy,dj,1,jmax-1,k,
     &                  at,am,ap)
         elseif(iflxo .eq. 3) then
c DWB: fluxj3 checks out with original code
            call fluxj3(jmax,kmax,q,rtxy,dj,1,jmax-1,k,
     &                  at,am,ap)
         elseif(iflxo .eq. 5) then
            call fluxj5(jmax,kmax,q,rtxy,dj,1,jmax-1,k,
     &                  at,am,ap)
         else
            print*
            print*, ' ERROR: invalid value for iflxo in conup.f'
            print*, '   Check input parameters.'
            print*, '   iflxo must be 1, 3, or 5.'
            print*, '   Current value iflxo =',iflxo
            print*, ' Program terminated.'
            print*
            stop ' stop: error iflxo value, stopping in conup.f'
         endif

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' j  k   btc11    btc21   btc31    btc12   btc22   btc32'
      endif

          do 80 j=2,jmax-1
            bjm(j,k,1,1) = 0.5*(-at(j-1,1,1)-ap(j-1,1,1)+am(j-1,1,1))
            bjp(j,k,1,1) = 0.5*( at(j+1,1,1)-ap(j,1,1) + am(j,1,1) )
            btc(j,k,1,1) = 0.5*( ap(j,1,1) + ap(j-1,1,1)
     &                          - am(j,1,1) - am(j-1,1,1) )
            bjm(j,k,2,1) = 0.5*(-at(j-1,2,1)-ap(j-1,2,1)+am(j-1,2,1))
            bjp(j,k,2,1) = 0.5*( at(j+1,2,1)-ap(j,2,1) + am(j,2,1) )
            btc(j,k,2,1) = 0.5*( ap(j,2,1) + ap(j-1,2,1)
     &                          - am(j,2,1) - am(j-1,2,1) )
            bjm(j,k,3,1) = 0.5*(-at(j-1,3,1)-ap(j-1,3,1)+am(j-1,3,1))
            bjp(j,k,3,1) = 0.5*( at(j+1,3,1)-ap(j,3,1) + am(j,3,1) )
            btc(j,k,3,1) = 0.5*( ap(j,3,1) + ap(j-1,3,1)
     &                          - am(j,3,1) - am(j-1,3,1) )
            bjm(j,k,1,2) = 0.5*(-at(j-1,1,2)-ap(j-1,1,2)+am(j-1,1,2))
            bjp(j,k,1,2) = 0.5*( at(j+1,1,2)-ap(j,1,2) + am(j,1,2) )
            btc(j,k,1,2) = 0.5*( ap(j,1,2) + ap(j-1,1,2)
     &                          - am(j,1,2) - am(j-1,1,2) )
            bjm(j,k,2,2) = 0.5*(-at(j-1,2,2)-ap(j-1,2,2)+am(j-1,2,2))
            bjp(j,k,2,2) = 0.5*( at(j+1,2,2)-ap(j,2,2) + am(j,2,2) )
            btc(j,k,2,2) = 0.5*( ap(j,2,2) + ap(j-1,2,2)
     &                          - am(j,2,2) - am(j-1,2,2) )
            bjm(j,k,3,2) = 0.5*(-at(j-1,3,2)-ap(j-1,3,2)+am(j-1,3,2))
            bjp(j,k,3,2) = 0.5*( at(j+1,3,2)-ap(j,3,2) + am(j,3,2) )
            btc(j,k,3,2) = 0.5*( ap(j,3,2) + ap(j-1,3,2)
     &                          - am(j,3,2) - am(j-1,3,2) )
            bjm(j,k,1,3) = 0.5*(-at(j-1,1,3)-ap(j-1,1,3)+am(j-1,1,3))
            bjp(j,k,1,3) = 0.5*( at(j+1,1,3)-ap(j,1,3) + am(j,1,3) )
            btc(j,k,1,3) = 0.5*( ap(j,1,3) + ap(j-1,1,3)
     &                          - am(j,1,3) - am(j-1,1,3) )
            bjm(j,k,2,3) = 0.5*(-at(j-1,2,3)-ap(j-1,2,3)+am(j-1,2,3))
            bjp(j,k,2,3) = 0.5*( at(j+1,2,3)-ap(j,2,3) + am(j,2,3) )
            btc(j,k,2,3) = 0.5*( ap(j,2,3) + ap(j-1,2,3)
     &                          - am(j,2,3) - am(j-1,2,3) )
            bjm(j,k,3,3) = 0.5*(-at(j-1,3,3)-ap(j-1,3,3)+am(j-1,3,3))
            bjp(j,k,3,3) = 0.5*( at(j+1,3,3)-ap(j,3,3) + am(j,3,3) )
            btc(j,k,3,3) = 0.5*( ap(j,3,3) + ap(j-1,3,3)
     &                          - am(j,3,3) - am(j-1,3,3) )
      if (DEBUG.and.k.eq.2.and.j.le.17.and.nodeid.eq.0) then
       write(*,907) j,k,btc(j,k,1,1),btc(j,k,2,1),btc(j,k,3,1),
     &  btc(j,k,1,2),btc(j,k,2,2),btc(j,k,3,2)
907    format(2i4,1p6e13.5)
      endif

80       continue
 
      if (DEBUG) then
       call flush(6)
       call flush(istdout)
c       stop 'stop: in conup.f/conup 60'
      endif

         do 90 j=2,jmax-1
            s(j,k,1) = - ( f(j,1) - f(j-1,1) )
            s(j,k,2) = - ( f(j,2) - f(j-1,2) )
            s(j,k,3) = - ( f(j,3) - f(j-1,3) )
90       continue

100   continue

      if (DEBUG.and.nodeid.eq.0) then
c       stop 'stop: in conup.f/conup 70'
      endif
      
      deallocate(at,am,ap)

c-----
c  K-sweep non-periodic
c-----
      allocate(at(kmax,3,3),am(kmax,3,3),ap(kmax,3,3))

         do 300 j=2,jmax-1

            do 279 k=1,kmax
               cvel = rtxy(j,k,2,1) + rtxy(j,k,2,2)*q(j,k,2)
     &              + rtxy(j,k,2,3)*q(j,k,3)
               at(k,1,1) = rtxy(j,k,2,1)
               at(k,2,1) = rtxy(j,k,2,2)
               at(k,3,1) = rtxy(j,k,2,3)
               at(k,1,2) = rtxy(j,k,2,2)*beta
               at(k,2,2) = rtxy(j,k,2,2)*q(j,k,2) + cvel
               at(k,3,2) = rtxy(j,k,2,2)*q(j,k,3)
               at(k,1,3) = rtxy(j,k,2,3)*beta
               at(k,2,3) = rtxy(j,k,2,3)*q(j,k,2)
               at(k,3,3) = rtxy(j,k,2,3)*q(j,k,3) + cvel
279         continue
c
c            call amatpmk(jmax,kmax,q,rtxy,dj,1,kmax-1,j)
c-----
c  Compute matrix
c-----
      do 35 k=1,kmax-1
         rt = 0.5*( rtxy(j,k+1,2,1) + rtxy(j,k,2,1) )
         rx = 0.5*( rtxy(j,k+1,2,2) + rtxy(j,k,2,2) )
         ry = 0.5*( rtxy(j,k+1,2,3) + rtxy(j,k,2,3) )
         q2 = 0.5*( q(j,k+1,2) + q(j,k,2) )
         q3 = 0.5*( q(j,k+1,3) + q(j,k,3) )
         cvel = rx*q2 + ry*q3
         cc = sqrt( (cvel+0.5*rt)**2 + beta*(rx*rx + ry*ry) )
         eig = cvel + rt
         eig1 = 0.5*(eig + abs(eig))
         eig2 = cvel + rt*0.5 + cc
c-----
c  Load A+ matrix
c-----
         x21 =-ry*eig1/(cc*cc - 0.25*rt*rt)
         x31 = rx*eig1/(cc*cc - 0.25*rt*rt)
         x12 = 0.5*eig2/cc
         x22 = 0.5*eig2*(q2*eig2 + beta*rx)/(beta*cc*(cc-0.5*rt))
         x32 = 0.5*eig2*(q3*eig2 + beta*ry)/(beta*cc*(cc-0.5*rt))
         xi11 = ry*q2 - rx*q3
         xi12 =-q3*(cvel + rt) - ry*beta
         xi13 = q2*(cvel + rt) + rx*beta
         xi21 = cc - cvel - 0.5*rt
         xi22 = rx*beta
         xi23 = ry*beta
c-----
         ap(k,1,1) =            x12*xi21
         ap(k,2,1) = x21*xi11 + x22*xi21
         ap(k,3,1) = x31*xi11 + x32*xi21
         ap(k,1,2) =            x12*xi22
         ap(k,2,2) = x21*xi12 + x22*xi22
         ap(k,3,2) = x31*xi12 + x32*xi22
         ap(k,1,3) =            x12*xi23
         ap(k,2,3) = x21*xi13 + x22*xi23
         ap(k,3,3) = x31*xi13 + x32*xi23
c-----
c  Load A- matrix
c-----
         am(k,1,1) = -ap(k,1,1)
         am(k,2,1) = -ap(k,2,1) + rx
         am(k,3,1) = -ap(k,3,1) + ry
c
         am(k,1,2) = -ap(k,1,2) + rx*beta
         am(k,2,2) = -ap(k,2,2) + rt + 2.*rx*q2 + ry*q3
         am(k,3,2) = -ap(k,3,2) + rx*q3
c
         am(k,1,3) = -ap(k,1,3) + ry*beta
         am(k,2,3) = -ap(k,2,3) + ry*q2
         am(k,3,3) = -ap(k,3,3) + rt + rx*q2 + 2.*ry*q3
35    continue
c

c  Default: iflxo = 3
            if(iflxo .eq. 1) then
               call fluxk1(jmax,kmax,q,rtxy,dj,1,kmax-1,j,
     &                     at,am,ap)
            elseif(iflxo .eq. 3) then
               call fluxk3(jmax,kmax,q,rtxy,dj,1,kmax-1,j,
     &                     at,am,ap)
            elseif(iflxo .eq. 5) then
               call fluxk5(jmax,kmax,q,rtxy,dj,1,kmax-1,j,
     &                     at,am,ap)
            endif
c
      if (DEBUG.and.nodeid.eq.0.and.j.le.17) then
       print*
       print*
       print*
       print*,'-j   k   btc11    btc21    btc31   btc12   btc22   btc32'
      endif
c
            do 280 k=2,kmax-1
               bkm(j,k,1,1)=0.5*(-at(k-1,1,1) -ap(k-1,1,1)+am(k-1,1,1))
               bkp(j,k,1,1)=0.5*( at(k+1,1,1) -ap(k,1,1) + am(k,1,1) )
               btc(j,k,1,1)=btc(j,k,1,1) +.5*(ap(k,1,1) + ap(k-1,1,1)
     &                                        - am(k,1,1) - am(k-1,1,1))
               bkm(j,k,2,1)=0.5*(-at(k-1,2,1) -ap(k-1,2,1)+am(k-1,2,1))
               bkp(j,k,2,1)=0.5*( at(k+1,2,1) -ap(k,2,1) + am(k,2,1) )
               btc(j,k,2,1)=btc(j,k,2,1) +.5*(ap(k,2,1) + ap(k-1,2,1)
     &                                        - am(k,2,1) - am(k-1,2,1))
               bkm(j,k,3,1)=0.5*(-at(k-1,3,1) -ap(k-1,3,1)+am(k-1,3,1))
               bkp(j,k,3,1)=0.5*( at(k+1,3,1) -ap(k,3,1) + am(k,3,1) )
               btc(j,k,3,1)=btc(j,k,3,1) +.5*(ap(k,3,1) + ap(k-1,3,1)
     &                                        - am(k,3,1) - am(k-1,3,1))
               bkm(j,k,1,2)=0.5*(-at(k-1,1,2) -ap(k-1,1,2)+am(k-1,1,2))
               bkp(j,k,1,2)=0.5*( at(k+1,1,2) -ap(k,1,2) + am(k,1,2) )
               btc(j,k,1,2)=btc(j,k,1,2) +.5*(ap(k,1,2) + ap(k-1,1,2)
     &                                        - am(k,1,2) - am(k-1,1,2))
               bkm(j,k,2,2)=0.5*(-at(k-1,2,2) -ap(k-1,2,2)+am(k-1,2,2))
               bkp(j,k,2,2)=0.5*( at(k+1,2,2) -ap(k,2,2) + am(k,2,2) )
               btc(j,k,2,2)=btc(j,k,2,2) +.5*(ap(k,2,2) + ap(k-1,2,2)
     &                                        - am(k,2,2) - am(k-1,2,2))
               bkm(j,k,3,2)=0.5*(-at(k-1,3,2) -ap(k-1,3,2)+am(k-1,3,2))
               bkp(j,k,3,2)=0.5*( at(k+1,3,2) -ap(k,3,2) + am(k,3,2) )
               btc(j,k,3,2)=btc(j,k,3,2) +.5*(ap(k,3,2) + ap(k-1,3,2)
     &                                        - am(k,3,2) - am(k-1,3,2))
               bkm(j,k,1,3)=0.5*(-at(k-1,1,3) -ap(k-1,1,3)+am(k-1,1,3))
               bkp(j,k,1,3)=0.5*( at(k+1,1,3) -ap(k,1,3) + am(k,1,3) )
               btc(j,k,1,3)=btc(j,k,1,3) +.5*(ap(k,1,3) + ap(k-1,1,3)
     &                                        - am(k,1,3) - am(k-1,1,3))
               bkm(j,k,2,3)=0.5*(-at(k-1,2,3) -ap(k-1,2,3)+am(k-1,2,3))
               bkp(j,k,2,3)=0.5*( at(k+1,2,3) -ap(k,2,3) + am(k,2,3) )
               btc(j,k,2,3)=btc(j,k,2,3) +.5*(ap(k,2,3) + ap(k-1,2,3)
     &                                        - am(k,2,3) - am(k-1,2,3))
               bkm(j,k,3,3)=0.5*(-at(k-1,3,3) -ap(k-1,3,3)+am(k-1,3,3))
               bkp(j,k,3,3)=0.5*( at(k+1,3,3) -ap(k,3,3) + am(k,3,3) )
               btc(j,k,3,3)=btc(j,k,3,3) +.5*(ap(k,3,3) + ap(k-1,3,3)
     &                                        - am(k,3,3) - am(k-1,3,3))

      if (DEBUG.and.nodeid.eq.0.and.k.le.10.and.j.le.17) then
       write(*,907) j,k,btc(j,k,1,1),btc(j,k,2,1),btc(j,k,3,1),
     &  btc(j,k,1,2),btc(j,k,2,2),btc(j,k,3,2)
      endif

280         continue

c print initial values of s
      if (DEBUG.and.nodeid.eq.0.and.j.le.17) then
        print*
        print*, ' Initial values of s()'
        print*, ' nodeid j  k  sjk1      sjk2      sjk3     fk1     fk2 
     &      fk3'
c        do 289 k=2,kmax-1
        do 289 k=2,10
         write(*,288) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3),f(k,1),
     &    f(k,2),f(k,3)
288      format(i4,2i4,1p6e13.5)
289    continue
      endif

      if (DEBUG.and.nodeid.eq.0.and.j.le.17) then
       print*
       print*
      endif

            do 290 k=2,kmax-1
      if (DEBUG.and.nodeid.eq.0.and.j.le.17) then
       print*,' j,k = ',j,k
       print*, ' in 290, before: s1,f1,fm1 = ',s(j,k,1),f(k,1),f(k-1,1)
       print*, '         before: s2,f2,fm1 = ',s(j,k,2),f(k,2),f(k-1,2)
       print*, '         before: s3,f3,fm1 = ',s(j,k,3),f(k,3),f(k-1,3)
      endif 
               s(j,k,1) = s(j,k,1) - ( f(k,1) - f(k-1,1) )
               s(j,k,2) = s(j,k,2) - ( f(k,2) - f(k-1,2) )
               s(j,k,3) = s(j,k,3) - ( f(k,3) - f(k-1,3) )

      if (DEBUG.and.nodeid.eq.0.and.j.le.17) then
       print*, '        after: s1,f1,fm1 = ',s(j,k,1),f(k,1),f(k-1,1)
       print*, '        after: s2,f2,fm1 = ',s(j,k,2),f(k,2),f(k-1,2)
       print*, '        after: s3,f3,fm1 = ',s(j,k,3),f(k,3),f(k-1,3)
      endif
290         continue

300      continue

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' In conup.f/conup at end:'
       print*, ' btc(2,2,1,1),btc(2,2,2,2) =',btc(2,2,1,1),btc(2,2,2,2)
       print*
       call flush(6)
       call flush(istdout)
c       print*, ' stop: AT END OF conup.f/conup'
c       stop 'stop: conup.f/conup 999'
      endif

      deallocate(at,am,ap)

c-----
c  End of conup
c-----
      return
      end
