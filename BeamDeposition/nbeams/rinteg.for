c//////////////////////////////////////////////////////////////////////////
c/
c/ NAME
c/    rinteg
c/
c/ DESCRIPTION
c/    This routine calculates the part of the integral along the
c/    major radius R at a given plane Z=zb, for the integration
c/    on the face of the beam.
c/    Written by John Mandrekas, GIT, for the SuperCode
c/    Revisions:
c/    04-14-92, John Mandrekas:
c/       Modified to handle beams with three energy components more
c/       efficiently
c/    02-10-92, John Mandrekas:
c/       Reversed integration limits when calling SINTEG
c/       Defined -nstart- for better interaction with GETRHO
c/       Allowed for uniform power distribution in case of circular
c/        beams (rectangular beams already had this capability)
c/
c/ PARAMETERS
c/    None
c/
c/ RETURNS
c/    Nothing
c/
c/////////////////////////////////////////////////////////////////////////

      subroutine rinteg

      implicit none

      include 'nbparams.inc'
      include 'BeamsLoc.inc'
      include 'nbplasma.inc'
      
c/    Local variable declarations:
c/    ---------------------------
      integer ie, ir
      real rterm, rbterm, rwterm, rbin, rbout, rwin, rwout, rcut,
     .     rl, ru, ruprl, rumrl, sumr1, sumr2, sumr3, dblpass, rsml2,
     .     rtest, rtoler, rcoeff, pwrfact, pwrterm, pwzterm1, pwzterm2, 
     .     ksibeam, sumksi1, sumksi2, sumksi3, sumterm

      real  dterm(3)
      data rtoler /1.0e-05/
            
c/    Intersection of beam with zb plane (rbin, rbout):
c/    ------------------------------------------------
c/    shape=0: circular cross section
c/    shape=1: rectangular cross section

      if (shape.eq.0) then
         rbterm = Sqrt(rbeam**2 -zb**2)
      else if (shape.eq.1) then
         rbterm = 0.5*width
      endif

      rterm  = Sqrt(rhoi**2 - (zb/elongr)**2)
      rbin   = rtan - rbterm
      rbout  = rtan + rbterm
      rwterm = Sqrt(rhoSurf**2 - (zb/elonga)**2)
      rwin   = rmajor - rwterm 
      rwout  = rmajor + rwterm

c/    Intersection of flux surface -rho- with zb plane (rflux):

      rflux = rmajor + shft + sgn*rterm
      
c/    If rflux < rbin, no intersection of beamline with rho:

      if (rflux.le.rbin) then
         do ie = 1, 3
            alongr(ie) = 0.0
            alongksi(ie)= 0.0
         enddo
         return
      endif
      
c/    Calculate limits of -rb- integration:

      rl = rbin
      ru = min(rflux, rbout)
      
c/    Start coordinate transformations for Gaussian integration:

      ruprl = ru + rl
      rumrl = ru - rl
      
c/    Avoid very small integration intervals

      rtest = rumrl/ruprl
      if (rtest.lt.rtoler) then
         do ie = 1, 3
            alongr(ie) = 0.0
            alongksi(ie)= 0.0
         enddo
         return
      endif

c/    Start integration along R:

      sumr1 = 0.0
      sumr2 = 0.0
      sumr3 = 0.0
      sumksi1 = 0.0
      sumksi2 = 0.0
      sumksi3 = 0.0
      do ir = 1, nr
         rb = 0.5*(rumrl*sr(ir) + ruprl)
         ksibeam = rb/rflux
         
c/    Check for double pass, and calculate integral along beam path
         
         rcut = rflux   
         nstart = 0     ! needed by GETRHO
         call sinteg (rwout, rcut, capd0)
         if (rb.lt.rwin) then
            do ie = 1, 3
               capd1(ie) = 0.0
            enddo
            dblpass = 0.0
         else  
            call sinteg (rcut, rb, capd1)
            dblpass = 1.0
         endif
         do ie = 1, 3
            if (fbpwr(ie,ibeam).ne.0.0) then
               dterm(ie) = exp(capd0(ie)) + 
     .           dblpass*exp(capd0(ie) + 2.0*capd1(ie))
            endif
         enddo

c/    Power distribution on the surface of the beam:
c/    ---------------------------------------------
c/       ptype = 0: uniform distribution
c/       ptype = 1: Gaussian or bi-Gaussian

         if (shape.eq.0) then       
            if (ptype.eq.0) pwrfact = 1.0
            if (ptype.eq.1) then
               rsml2 = zb**2 + (rtan - rb)**2
               pwrfact = exp(-rsml2/(gaussR**2))
            endif
         else if (shape.eq.1) then
            if (ptype.eq.0) pwrfact = 1.0
            if (ptype.eq.1) then
               pwrterm  = (rb-rtan)/gaussR
               pwzterm1 = (zb-zpos)/gaussZ
               pwzterm2 = (zb+zpos)/gaussZ
               pwrfact  = exp(-pwrterm**2)*(
     .           exp(-pwzterm1**2) + exp(-pwzterm2**2))
            endif
         endif

         sumterm = rflux*pwrfact*Sqrt((1.0 - sr(ir)**2)/
     .             (rflux**2 - rb**2))
         sumr1 = sumr1 + sumterm*dterm(1)
         sumr2 = sumr2 + sumterm*dterm(2)
         sumr3 = sumr3 + sumterm*dterm(3)
         sumksi1 = sumksi1 + sumterm*dterm(1)*ksibeam
         sumksi2 = sumksi2 + sumterm*dterm(2)*ksibeam
         sumksi3 = sumksi3 + sumterm*dterm(3)*ksibeam
      enddo
      
      rcoeff = 0.5*rumrl*wr
      alongr(1) = rcoeff*sumr1
      alongr(2) = rcoeff*sumr2
      alongr(3) = rcoeff*sumr3
      alongksi(1) = rcoeff*sumksi1
      alongksi(2) = rcoeff*sumksi2
      alongksi(3) = rcoeff*sumksi3

      return
      end                        
