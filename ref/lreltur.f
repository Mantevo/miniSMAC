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
c  1. lreltur
c  2. turswpj
c  3. turswpk
c  4. scatri
c
c************************************************************
c
c
c************************************************************************
      subroutine lreltur(jmax,kmax,turvar,ds,rhs,
     &                   sjm,sjp,skm,skp,stc,nodeidp1)
c************************************************************************
c
c Purpose: solve system of equations with turbulence using line relaxation
c
c------------------------------------------------------------------------
#include "common.f"
#include "mpi_params.f"
      __REAL turvar(jmax,kmax), ds(jmax,kmax), rhs(jmax,kmax), 
     &     sjm(jmax,kmax),
     &     sjp(jmax,kmax), skm(jmax,kmax), skp(jmax,kmax), 
     &     stc(jmax,kmax)
      __INTEGER jmax, kmax
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG) then
       print*
       print*, ' >> Entering lreltur.f/lreltur, nodeid =',nodeid
       print*, '      nodeid,jmax,kmax = ',nodeid,jmax,kmax
       print*, '      ntjsp,ntksp = ',ntjsp,ntksp
       print*
      endif
c-----
c  Solve system of equations using line-relaxation
c  call bcimpds for implicit injection of boundary
c  info during sweep process
c-----
      do 20 k=1,kmax
      do 20 j=1,jmax
       ds(j,k) = 0.0
20    continue

      nsmax = 0
      nsmax = max(nsmax, ntjsp, ntksp )

      if (DEBUG) then
       print*
       print*, ' nodeid, nsmax = ',nodeid,nsmax
       print*
c       print*, ' stop: 10 in lreltur.f/lreltur'
c       stop 'stop: 10 in lreltur.f/lreltur'
      endif

      neqs = 1
      if (DEBUG) then
       print*, ' start: do 70 -- ntjsp,nsmax,neqs = ',ntjsp,nsmax,
     1    neqs
      endif

       do 70 njk=1,nsmax
c         do 40 nz=1,nzone
c-----
c  j-sweeps
c-----
        if (DEBUG) then
         print*, ' >> For j sweeps: njk = ',njk
         print*
        endif
            if(ntjsp .ge. njk) then
             if (DEBUG) then
              print*, ' -- call turswpj'
             endif
               call turswpj(jmax,kmax,sjm(1,1),
     &            sjp(1,1),skm(1,1),skp(1,1),stc(1,1),ds(1,1),
     &            rhs(1,1),njk)
               if (DEBUG) then
                print*, ' -- call bcimpds for j, nodeidp1 = ',nodeidp1
               endif
               call bcimpds(jmax,kmax,neqs,turvar(1,1),ds(1,1),
     &            rhs(1,1),nodeidp1)
c-----
            endif
40       continue

c         do 50 nz=1,nzone
c-----
c  k-sweeps
c-----
         if (DEBUG) then
          print*, ' >> For k sweeps: njk = ',njk
          print*
         endif
            if(ntksp .ge. njk) then
               call turswpk(jmax,kmax,sjm(1,1),
     &            sjp(1,1),skm(1,1),skp(1,1),stc(1,1),ds(1,1),
     &            rhs(1,1),njk)
               if (DEBUG) then
                print*, ' -- call bcimpds for k, nodeidp1 = ',nodeidp1
               endif
               call bcimpds(jmax,kmax,neqs,turvar(1,1),ds(1,1),
     &            rhs(1,1),nodeidp1)
            endif
50       continue
70    continue
c-----
c  End of lreltur
c-----
      if (DEBUG) then
       print*
       print*, ' stop: 20 in lreltur.f/lreltur at end'
       stop ' stop: 20 in lrelture.f/lreltur at end'
      endif

      return
      end
c
c
c************************************************************************
      subroutine turswpj(jmax,kmax,sjm,sjp,skm,skp,stc,ds,rhs,nj)
c************************************************************************
#include "common.f"
      __REAL sjm(jmax,kmax),sjp(jmax,kmax),
     &       skm(jmax,kmax),skp(jmax,kmax),
     &       stc(jmax,kmax), ds(jmax,kmax), rhs(jmax,kmax)
c
      __REAL aaa(jkmax), bbb(jkmax), ccc(jkmax), ff(jkmax)
c-----
c  Solution by line-relaxation: implicit in j-direction
c-----
c      if( (nj/2)*2 .eq. nj) then
c         kbeg = kmax
c         kend = 1
c         kinc = -1
c      else
c         kbeg = 1
c         kend = kmax
c         kinc = 1
c      endif
c-----
      kbeg = 1
      kend = kmax
      kinc = 1
      do 110 k=kbeg,kend,kinc
         if(k .eq. 1) then
               do 103 j=1,jmax
                  aaa(j) = sjm(j,k)
                  bbb(j) = stc(j,k)
                  ccc(j) = sjp(j,k)
                  ff(j) = rhs(j,k) - skp(j,k)*ds(j,k+1)
103            continue
         elseif(k .eq. kmax) then
               do 105 j=1,jmax
                  aaa(j) = sjm(j,k)
                  bbb(j) = stc(j,k)
                  ccc(j) = sjp(j,k)
                  ff(j) = rhs(j,k) - skm(j,k)*ds(j,k-1)
105            continue
         else
            do 106 j=1,jmax
               aaa(j) = sjm(j,k)
               bbb(j) = stc(j,k)
               ccc(j) = sjp(j,k)
               ff(j) = rhs(j,k) -skm(j,k)*ds(j,k-1) -skp(j,k)*ds(j,k+1)
106         continue
         endif
         call scatri(jmax,1,jmax,aaa,bbb,ccc,ff)
         do 108 j=1,jmax
            ds(j,k) = ff(j)
108      continue
110   continue

c-----
c  End point
c-----
      k = kbeg
      if(k .eq. 1) then
            do 203 j=1,jmax
               aaa(j) = sjm(j,k)
               bbb(j) = stc(j,k)
               ccc(j) = sjp(j,k)
               ff(j) = rhs(j,k) - skp(j,k)*ds(j,k+1)
203         continue
      elseif(k .eq. kmax) then
            do 205 j=1,jmax
               aaa(j) = sjm(j,k)
               bbb(j) = stc(j,k)
               ccc(j) = sjp(j,k)
               ff(j) = rhs(j,k) - skm(j,k)*ds(j,k-1)
205         continue
      else
         do 206 j=1,jmax
            aaa(j) = sjm(j,k)
            bbb(j) = stc(j,k)
            ccc(j) = sjp(j,k)
            ff(j) = rhs(j,k) - skm(j,k)*ds(j,k-1) - skp(j,k)*ds(j,k+1)
206      continue
      endif
      call scatri(jmax,1,jmax,aaa,bbb,ccc,ff)
      do 208 j=1,jmax
         ds(j,k) = ff(j)
208   continue

c-----
c  End of turswpj
c-----
      return
      end
c
c
c************************************************************************
      subroutine turswpk(jmax,kmax,sjm,sjp,skm,skp,stc,ds,rhs,nk)
c************************************************************************
#include "common.f"
      __REAL sjm(jmax,kmax),sjp(jmax,kmax),
     &       skm(jmax,kmax),skp(jmax,kmax),
     &       stc(jmax,kmax), ds(jmax,kmax), rhs(jmax,kmax)
c
      __REAL aaa(jkmax), bbb(jkmax), ccc(jkmax), ff(jkmax)
c-----
c  Solution by line-relaxation: implicit in k-direction
      jbeg = 1
      jend = jmax
      jinc = 1
      do 110 j=jbeg,jend,jinc
         if(j .eq. 1) then
            do 102 k=1,kmax
               aaa(k) = skm(j,k)
               bbb(k) = stc(j,k)
               ccc(k) = skp(j,k)
               ff(k) = rhs(j,k) - sjp(j,k)*ds(j+1,k)
102         continue
         elseif(j .eq. jmax) then
            do 104 k=1,kmax
               aaa(k) = skm(j,k)
               bbb(k) = stc(j,k)
               ccc(k) = skp(j,k)
               ff(k) = rhs(j,k) - sjm(j,k)*ds(j-1,k)
104         continue
         else
            do 106 k=1,kmax
               aaa(k) = skm(j,k)
               bbb(k) = stc(j,k)
               ccc(k) = skp(j,k)
               ff(k) =rhs(j,k) -sjm(j,k)*ds(j-1,k) -sjp(j,k)*ds(j+1,k)
106         continue
         endif
         call scatri(kmax,1,kmax,aaa,bbb,ccc,ff)
         do 108 k=1,kmax
            ds(j,k) = ff(k)
108      continue
110   continue

c-----
c  End point
c-----
      j = jbeg
      if(j .eq. 1) then
         do 202 k=1,kmax
            aaa(k) = skm(j,k)
            bbb(k) = stc(j,k)
            ccc(k) = skp(j,k)
            ff(k) = rhs(j,k) - sjp(j,k)*ds(j+1,k)
202      continue
      elseif(j .eq. jmax) then
         do 204 k=1,kmax
            aaa(k) = skm(j,k)
            bbb(k) = stc(j,k)
            ccc(k) = skp(j,k)
            ff(k) = rhs(j,k) - sjm(j,k)*ds(j-1,k)
204      continue
      else
         do 206 k=1,kmax
            aaa(k) = skm(j,k)
            bbb(k) = stc(j,k)
            ccc(k) = skp(j,k)
            ff(k) =rhs(j,k) -sjm(j,k)*ds(j-1,k) -sjp(j,k)*ds(j+1,k)
206      continue
      endif
      call scatri(kmax,1,kmax,aaa,bbb,ccc,ff)
      do 208 k=1,kmax
         ds(j,k) = ff(k)
208   continue

c-----
c  End of turswpk
c-----
      return
      end
c
c
c***********************************************************************
      subroutine scatri(imax,ibeg,iend,a,b,c,f)
c***********************************************************************
#include "precis.h"
      __REAL a(imax), b(imax), c(imax), f(imax)
c
      i = ibeg
      f(i) = f(i)/b(i)
      c(i) = c(i)/b(i)
      do 10 i=ibeg+1,iend-1
         z = 1./(b(i) - a(i)*c(i-1))
         c(i) = c(i)*z
         f(i) = (f(i) - a(i)*f(i-1))*z
10    continue
      i = iend
      f(i) = (f(i) - a(i)*f(i-1))/(b(i) - a(i)*c(i-1))
      do 20 i=iend-1,ibeg,-1
         f(i) = f(i) - c(i)*f(i+1)
20    continue
c-----
c  End of scatri
c-----
      return
      end
