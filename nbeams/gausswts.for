c///////////////////////////////////////////////////////////////////////
c/
c/ NAME 
c/   Beams_gaussWts
c/
c/ DESCRIPTION
c/    calculates Gaussian weights and abscissae for the integrations
c/    needed in the neutral beam calculations
c/
c/ PARAMETERS
c/    n:    number of nodes
c/ 
c/ RETURNS
c/    w:    weight factor
c/    s:    vector with the weight functions
c/
c///////////////////////////////////////////////////////////////////////

      subroutine gaussWts (n, w, sgauss)

c/    Local Variable Declarations:
c/    ---------------------------
      integer i, n
      real w, sgauss(n), pi
      
      data pi /3.14159265/
 
      w = pi/float(n)
      do i = 1, n 
         sgauss(i) = cos(pi*(2.0*float(i) - 1.0)/
     .             (2.0*float(n)))
      enddo
 
      return
      end
