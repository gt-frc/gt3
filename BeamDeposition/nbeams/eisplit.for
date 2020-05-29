c///////////////////////////////////////////////////////////////////////
c/
c/ NAME
c/    eiSplit
c/
c/ DESCRIPTION
c/    This routine calculates the fraction of energy deposited to
c/    plasma ions and electrons during the slowing down of a fast 
c/    ion in a plasma, assuming classical, instantaneous slowing
c/    down.
c/
c/ PARAMETERS
c/    fastEnergy : initial energy of fast ions in keV
c/    fastMass   : atomic mass of fast ions
c/
c/ RETURNS
c/    fion(i)    : vector with the fraction of the energy going
c/                 to plasma ions at zone -i-
c/    felec(i)   : vector with the fraction of the energy going
c/                 to plasma electrons.
c/
c/ CHANGES
c/    16-aug-00  : Changes in preparation for NTCC submission
c/
c///////////////////////////////////////////////////////////////////////

      subroutine eiSplit (fastEnergy, fastMass, fion, felec)

  
      implicit none
      
      include 'nbparams.inc'
      include 'nbconsts.inc'
      include 'nbplasma.inc'
      

c/    Local declarations:
c/    ------------------
      real fastEnergy, fastMass, fion(mxrho), felec(mxrho)
      real coulLoge, coulLogi, rterm, root3, sumbrack, term1, term2,
     .   twoThirds, ecrit, xcr, zbracket
      integer i, j
      
      root3 = Sqrt(3.0)
      twoThirds = 2.0/3.0
      do i = 1, nrho
         coulLoge = 37.8 + Log(elecTemp(i)/Sqrt(1.0e20*elecDensity(i)))
         rterm = Sqrt(elecTemp(i) * fastEnergy * fastMass /
     .      elecDensity(i))
         sumbrack = 0.0
	 do j = 1, nions
            coulLogi =  19.1 + Log((atw(j) / (fastMass+atw(j)))*rterm)
            sumbrack = sumbrack + dni(i,j) * znum(j)**2 
     .         * coulLogi / atw(j)
         enddo
         zbracket = sumbrack / (coulLoge * elecDensity(i))

c/    Critical energy for fast ion slowing down:
c/    -----------------------------------------
         ecrit = 14.8 * (zbracket**twoThirds) * fastMass * elecTemp(i)
         xcr = Sqrt(fastEnergy/ecrit)
         term1 = Log((xcr**2+2.*xcr+1.)/(xcr**2-xcr+1))/6.0
         term2 = (atan((2.*xcr-1)/root3) + pi/6.)/root3
         fion(i) = 2.0 * (term2 - term1) / (xcr**2)
         felec(i) = 1.0 - fion(i)
      enddo

      return
      end
