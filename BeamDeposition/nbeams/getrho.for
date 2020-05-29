c///////////////////////////////////////////////////////////////////////
c/
c/ NAME
c/    Beams_getrho
c/
c/ DESCRIPTION
c/   This subroutine finds the value of rho (rhotest) for given
c/   R (rprime). It is based on an iteration scheme, using the
c/   routine HUNT.FOR.
c/   Written by John Mandrekas, GIT, 02-03-92, for the SuperCode
c/
c/ PARAMETERS
c/    rprime  :  R for which we seek rho
c/
c/ RETURNS
c/    rhotest : The calculated value of rho (normalized)
c/    jlo     : rhotest is between rho(jlo) and rho(jlo+1)
c/
c///////////////////////////////////////////////////////////////////////
    

      subroutine getrho (rprime,  rhotest, jlo)
      
      implicit none
      include 'nbparams.inc'
      include 'nbplasma.inc'
      include 'BeamsLoc.inc'

      integer i, jtest, jlo
      real rprime, rhotest


c/    nstart  : if nstart = 0, use the outer flux surface
c/              for a first guess. Otherwise, use the previous
c/              value of jtest (defined in subroutine rinteg).

      if (nstart.eq.0) then
         jtest = nrho
         nstart = 1
	 jlo = 0.0
      endif

      do i = 1, maxit
         rhotest = Sqrt((rprime - rmajor - 
     .     rminor*shift(jtest))**2 + (zb/elong(jtest))**2)/rminor
         call hunt (rho, nrho, rhotest, jlo)
	 if (jlo.le.0) jlo = 1
	 if (jlo.ge.nrho) jlo = nrho -1
         jtest = jlo + 1
      enddo

      return
      end
    

