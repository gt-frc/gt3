
      real function sig_olson (eb, zq)

c///////////////////////////////////////////////////////////////////////
c/
c/    This function calculates the total electron loss cross section 
c/    (charge exchange and ionization)  between hydrogenic neutrals
c/    and heavy impurities.
c/    Written by John Mandrekas, GIT, 08/31/99
c/
c/    Reference:
c/    ---------
c/    R.E. Olson, et. al, Phys. Rev. Letters, 41 (1978) 163.
      
c/    Input:
c/    ------
c/    eb  : neutral energy in keV/amu
c/    zq  : ion charge state

c/    Output:
c/    ------
c/    sig_olson : Stopping cross section in cm^2
c/
c///////////////////////////////////////////////////////////////////////
      
      implicit none
      real eb, zq, sig_loss, factr

      factr = 32.0 * zq / eb
      sig_loss = 4.6e-16 * zq * factr * (1.0 - exp(-1.0/factr))

      sig_olson = sig_loss

      return
      end
