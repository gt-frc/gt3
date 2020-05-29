      subroutine coulomb(clog, dne, te, ami, amf, zi, zf, efast, isp)

c/    This subroutine calculates the Coulomb logarithms due to
c/    collisions between fast ions and plasma ions and electrons
c/    Reference: 
c/    D.R. Mikkelsen, C.E. Singer, Nuclear Technol./Fusion 4, 237 (1983)

c/    Input variables:
c/    dne    : electron density in units of 10**14 cm-3
c/    te     : electron temperature in keV
c/    ami    : plasma ion mass number (amu)
c/    amf    : fast ion mass number
c/    zi     : plasma ion atomic number
c/    zf     : fast ion atomic number
c/    efast  : fast ion energy (keV)
c/    isp    : isp=1 -> electron coulomb logarithm
c/             isp=2 -> beam-plasma ion coulomb logarithm
c/    Output :
c/    clog   : coulomb logarithm
c/    
c/    07/07/98: Fixed bug in calculation of electron logarithm
      
      implicit none
      integer isp
      real ami, amf, clog, dne, efast, elim, te, te10, 
     .	   zi, zf
      te10 = te/10.0
      elim = 100.*amf*(zf*zi)**2

      if (isp.eq.1) then
         clog = 17.0 + alog10((te10)/sqrt(dne))
      else
         if (efast.le.elim) then
            clog = 25.4 + alog10((efast*1.e-3*ami/(zf*zi*(ami+amf)))*
     x        sqrt(te10/dne))
         else if (efast.gt.elim) then
            clog = 23.7 + alog10((ami/(ami+amf))*sqrt(te10*amf*efast
     x        *1.e-3/dne))
         endif
      endif

      return
      end
