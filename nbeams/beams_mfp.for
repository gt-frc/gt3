c/////////////////////////////////////////////////////////////////////////
c/
c/ NAME
c/    Beams_mfp
c/
c/ DESCRIPTION
c/    This subroutine calculates the inverse mean free path of
c/    hydrogenic neutrals in a plasma using the beam stopping cross
c/    section parametrization by C. Boley (Nucl. Fusion 29, 2125, 1989)
c/    Version 0.00, John Mandrekas, GIT, 06-25-91
c/    Updated fits (SIGFIT91), 03/27/92
c/
c/ PARAMETERS
c/    None
c/ 
c/ RETURNS
c/    Nothing
c/
c/////////////////////////////////////////////////////////////////////////

      subroutine Beams_mfp
      
      implicit none
      
      include 'nbparams.inc'
      include 'BeamsLoc.inc'
      include 'nbplasma.inc'

      real sig_olson

c/    Local variable declarations:
c/    ---------------------------      
      integer i, ie, iz, izq 
      real te, dene, eb_amu, sig, sigzol, zq


c/    Get ready to call the subroutine SIGFIT:

      do i = 1, nrho
         te = elecTemp(i)
         dene = 1.0e14*elecDensity(i)
	 if (ihe4.ne.0) then
	    denz(1) = 1.0e14*dni(i,ihe4)
         else
	    denz(1) = 0.0
         endif
	 if (icarb.ne.0) then
	    denz(2) = 1.0e14*dni(i,icarb)
         else
	    denz(2) = 0.0
         endif
	 if (ioxy.ne.0) then
	    denz(3) = 1.0e14*dni(i,ioxy)
         else
	    denz(3) = 0.0
         endif
	 if (iFe.ne.0) then
	    denz(4) = 1.0e14*dni(i,iFe)
         else
	    denz(4) = 0.0
         endif

c/    Calculate inverse mean free path for each energy group if required:

c/    For	 the impurity species NOT included in Boley's routine, 
c/    calculate stopping cross section using Olson's formula and
c/    add to the total inverse mean free path:

	 do ie = 1, 3
	    if (fbpwr(ie,ibeam).ne.0.0) then
	       eb_amu = ebkev/(ie*ab)
    	       call sigfit (eb_amu, dene, te, denz, sig)
	       omfp(i,ie,ibeam) = 1.0e16 * elecDensity(i) * sig

	       if (inonb.NE.0) then
	          sigzol = 0.0
		  do iz = 1, inonb
		     izq = nonbindx(iz)
		     zq = znum(izq)
		     sigzol = sigzol + dni(i,izq)*sig_olson(eb_amu,zq)
		  enddo
		  omfp(i,ie,ibeam) = omfp(i,ie,ibeam) + 1.0e16*sigzol
	       endif

            endif

         enddo   ! End of DO loop over energy components

      enddo      ! End of DO loop over radial points

      return
      end
