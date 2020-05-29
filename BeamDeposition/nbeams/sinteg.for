c//////////////////////////////////////////////////////////////////////////
c/
c/ NAME
c/    sinteg
c/
c/ DESCRIPTION
c/    This routine calculates the attenuation of a beamlet in the 
c/    plasma along its path, from R=r1 to R=r2. The integral is 
c/    calculated using Gaussian integration formulas.
c/    Written by John Mandrekas, GIT, for the SuperCode
c/    Revisions:
c/    04-14-92, John Mandrekas:
c/       Modified to handle beams with three energy components more
c/       efficiently.
c/    02-07-92, John Mandrekas:
c/       Modified to work with new and faster GETRHO
c/
c/ PARAMETERS
c/    r1, r2:  Lower and upper limits of integration
c/    alongs:  The result of the integration along the chord for
c/             each of the energy components
c/
c////////////////////////////////////////////////////////////////////////////

      subroutine sinteg (r1, r2, alongs)

      implicit none
      include 'nbparams.inc'
      include 'BeamsLoc.inc'
      include 'nbplasma.inc'

c/    Local Variable Declarations:
c/    ---------------------------      
      real r1, r2, alongs(3)
      integer ie, is, jlo 
      real u1, u2, u2pu1, u2mu1, u, sums1, sums2, sums3, 
     .     sterm, rprime, rhotest, stest, stoler, usqr
      real ovlambda(3)

      data stoler /1.0e-5/
           
c/    Coordinate transformation, u = rb/R:

      u1 = rb/r1
      u2 = rb/r2
      u2pu1 = u2 + u1
      u2mu1 = u2 - u1

c/    Avoid very small integration intervals:

      stest = abs(u2mu1)
      if (stest.lt.stoler) then
         do ie = 1, 3
            alongs(ie) = 0.0
         enddo
         return
      endif

      sums1 = 0.0
      sums2 = 0.0
      sums3 = 0.0
      do is = 1, ns
         u = 0.5*(u2mu1*ss(is) + u2pu1)
         usqr = u**2
         rprime = rb/u
         sterm = (1.0/usqr)*Sqrt((1.0 - ss(is)**2)/
     .                 (1.0 - usqr))
         call getrho (rprime, rhotest, jlo)
         do ie = 1, 3
            if (fbpwr(ie,ibeam).ne.0.0) then
               ovlambda(ie) = omfp(jlo,ie,ibeam) +
     .         (omfp(jlo+1,ie,ibeam) - omfp(jlo,ie,ibeam))*
     .         (rhotest - rho(jlo))/(rho(jlo+1) - rho(jlo))
            endif
         enddo
         sums1 = sums1 + ovlambda(1)*sterm
         sums2 = sums2 + ovlambda(2)*sterm
         sums3 = sums3 + ovlambda(3)*sterm
      enddo
      
      alongs(1) = -0.5*rb*u2mu1*ws*sums1
      alongs(2) = -0.5*rb*u2mu1*ws*sums2
      alongs(3) = -0.5*rb*u2mu1*ws*sums3
      
      return
      end
