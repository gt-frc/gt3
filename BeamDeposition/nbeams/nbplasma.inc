c/    File nbplasma.inc with background plasma information needed by the
c/    Beams module:

      integer nonbindx(maxIons)
      integer ihydr, ideut, itrit, ihe3, ihe4, ioxy, icarb, iFe, inonb
      
      common /nbspec/ nonbindx, inonb, ihydr, ideut, itrit, ihe3, ihe4, 
     .   ioxy, icarb, iFe
     
      integer nions
      real dni(mxrho, maxIons), znum(maxIons), atw(maxIons)
      common /nbat/ dni, znum, atw, nions
      
      real elecDensity(mxrho), elecTemp(mxrho), ionTemp(mxrho),
     .   zef(mxrho)
      common /nteti/ elecDensity, elecTemp, ionTemp, zef
      
c/    Grid and MHD related parameters:
c/    -------------------------------
      integer nrho
      real rmajor, rminor, btor, volume
      real rho(mxrho), vpr(mxrho), dv(mxrho), elong(mxrho), 
     .     elong_rho(mxrho), shift(mxrho), shift_rho(mxrho) 
      
      common /nbplas/ nrho, rmajor, rminor, btor, volume, rho, 
     .   vpr, dv, elong, elong_rho, shift, shift_rho 
