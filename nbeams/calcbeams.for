      subroutine calcBeams(nbeams, amb, zbeam, ebeam, pbeam, inbfus,
     .   rtang, nbshape, bwidth, bheigh, nbptype, bgaussR, bgaussZ,
     .   bzpos, pwrfrac, maxiter, nion, aion, zion, ne20, ni20, tekev,
     .   tikev, zeff, r0, a, b0, volp, n, rnorm, vprime, dvol, darea,
     .   kappa, dkappa, shafr, dshafr, hofr, shinethru, jnbTot, pnbe,
     .   pnbi, beamDens, beamPress, beamFus, pbfuse, pbfusi, snBeamDD,
     .   snBeamDT, nbcur, etanb, gammanb, pNBAbsorb, pNBLoss, nbcurTot,
     .   etanbTot, beamBeta, pNBAbsorbTot, pNBLossTot, beamFusTot,
     .   beamFusChTot, snDTTotal, snDDTotal, iflag,taus)

c/	 use MSIMSL

c///////////////////////////////////////////////////////////////////////
c/
c/    This is the main calling routine of the Neutral Beam package.
c/    It calls various routines to calculate the most common NB
c/    deposition, heating and fast ion quantities including NB current
c/    drive and beam-target fusion parameters.
c/
c/    This package was originally developed for the SuperCode project
c/    (~1992) by John Mandrekas. It has been modified by the author for
c/    submission to the National Transport Code Collaboration (NTCC).
c/
c/    Latest release: 08/22/2000
c/
c/    Questions, comments and bug reports should be addressed to:
c/
c/                 John Mandrekas
c/                 Georgia Institute of Technology
c/                 Fusion Research Center
c/                 Atlanta, GA 30332-0225
c/                 Tel: (404) 894-7730
c/                 E-mail: john.mandrekas@me.gatech.edu
c/
c/    DESCRIPTION OF VARIABLES IN THE SUBROUTINE CALLING LIST
c/
c/    Neutral beam related input parameters:
c/    -------------------------------------
c/    Notice: We can have as many beamlines as desired, with the
c/    following limitations:
c/
c/    1) The beam ion species should be the same for all beamlines
c/    2) If beam-target fusion calculations are desired, then the
c/       beam ions should be either Deuterium or He3. For NB deposition
c/       and current drive calculations, any H species is OK.
c/
c/    nbeams    : Number of beamlines
c/    amb       : Atomic weight of beam ions (1 for H, 2 for D, etc.)
c/    zbeam     : Atomic number of beam ions
c/    ebeam(ib) : Beam energy of each beamline (keV)
c/    pbeam(ib) : Beam power of each beamline  (MW)
c/
c/    pwrfr(ie,ib) : The fraction of power in the energy component
c/                   ie (ie = 1,2,3 corresponding to full, half and
c/                   1/3 energy components) for each beamline ib
c/
c/    inbfus      : Flag for the beam-target fusion calculations:
c/                  0 - Do not calculate
c/                  1 - Calculate
c/    rtang(ib)   : Tangency radius of each beamline (m)
c/    nbshape(ib) : Flag determining the cross section shape of each
c/                  beamline. Available options are 0 for circular
c/                  and 1 for rectangular.
c/    bwidth(ib)  : Width (or diameter for circular) of each beamline (m)
c/    bheigh(ib)  : Height of each beamline (m)
c/    nbptype(ib) : Flag for power distribution in the beam cross section.
c/                  0 : uniform distribution
c/                  1 : Gaussian (circular) or bi-Gaussian (rectangular)
c/                  p ~ exp(-(r/gaussR)**2) - circular
c/                  p ~ exp(-(R/gaussR)**2)*(exp(-((Z-bzpos)/gaussZ)**2 +
c/                                          exp(-((Z+bzpos)/gaussZ)**2))
c/
c/    bgaussR(ib) : Gaussian width (radial direction for circular beams,
c/                  horizontal direction for rectangular beams). (m)
c/    bgaussZ(ib) : Gaussian width in vertical direction for rectangular
c/                  beams (m)
c/    bzpos(ib)   : Position of peak in Z direction for bi-Gaussian
c/                  power distribution (rectangular shapes only) (m)
c/
c/    maxiter     : Maximum number of iterations in the GETRHO routine
c/                  which finds the rho corresponding to a given R.
c/                  Recommended value: 2
c/
c/    Background Plasma Input Quantities:
c/    ----------------------------------
c/    nion        : Number of ion species (includes main ions and impurities)
c/    aion(j)     : Atomic Mass of j-th ion species (j = 1,nion)
c/    zion(j)     : Atomic Number of j-th ion species
c/    ne20(i)     : Electron density at grid point -i- (10^20 /m^3)
c/    ni20(i,j)   : Ion density at -i- for ion species -j- (10^20 /m^3)
c/    tekev(i)    : Electron temperature at -i- (keV)
c/    tikev(i)    : Ion temperature at -i- (keV)
c/    zeff(i)     : Plasma Zeff at grid point -i-
c/
c/    Geometry and MHD Input Variables:
c/    --------------------------------
c/    r0          : Major radius (m)
c/    a           : Minor radius (m)
c/    b0          : Toroidal magnetic field (T)
c/    volp        : Plasma volume (m^3)
c/    n           : Number of radial grid points
c/    rnorm(i)    : Normalized radial coordinate (0 at center, 1 at edge)
c/    vprime(i)   : dVol/dr (m^2)
c/    dvol(i)     : Volume of cell -i- (m^3)
c/    darea(i)    : Cross sectional area of cell (i) (m^2)
c/    kappa(i)    : Elongation at grid point i
c/    dkappa(i)   : d_kappa / drho at -i- (rho is normalized here)
c/    shafr(i)    : Shafranov shift (normalized to minor radius)
c/    dshafr(i)   : d_shift/d_rho (both shift and rho normalized)
c/
c/    Calculated (output) Variables:
c/    -----------------------------
c/    hofr(i,ie,ib) : Neutral beam deposition profile for beam ib,
c/                    energy component ie, at grid point i. Notice that
c/                    the deposition profile is normalized to 1!
c/
c/    shinethru(ie,ib) : Beam shinethrough for beam ib and energy component
c/                       ie.
c/
c/    jnbTot(i)     : Total neutral beam driven current density from all
c/                    beam lines and energy groups at grid point i (A/m^2)
c/    pnbe(i)       : Beam heating power to electrons at i (MW/m^3)
c/    pnbi(i)       : Beam heating power to ions at i (MW/m^3)
c/    beamDens(i)   : Fast beam ion density at grid point i (/m^3)
c/    beamPress(i)  : Fast beam ion pressure at i (Pa)
c/    beamFus(i)    : Beam-target total fusion power density at grid
c/                    point i (this includes neutron AND charge particle
c/                    contributions) (MW/m^3).
c/    pNBLoss(ib)   : Neutral beam power lost per beamline (MW).
c/    pNBAbsorb(ib) : Neutral beam power absorbed per beamline (MW).
c/    pbfuse(i)     : Beam-target fusion power to electrons (MW/m^3)
c/    pbfusi(i)     : Beam-target fusion power to ions (MW/m^3)
c/    snBeamDD(i)   : Volumetric DD neutron rate production (#n/s/m^3)
c/    snBeamDT(i)   : Volumetric DT neutron rate production (#n/s/m^3)
c/    nbur(ib)      : NB driven current from beam ib (A)
c/    etanb(ib)     : NB current drive efficiency for beam ib (A/W)
c/    gammanb(ib)   : NBCD figure of merit for beam ib  (10^20 A/W/m^2)
c/    pNBAbsorbTot  : Total NB power absorbed in the plasma.
c/    pNBLossTot    : Total NB power lost.
c/    nbcurTot      : Total NB driven current (A).
c/    etanbTot      : Total current drive efficiency (A/W).
c/    beamBeta      : Total beta of the fast ions
c/    beamFusTot    : Total beam-target fusion power (MW).
c/    beamFusChTot  : Total charged particle beam fusion power (MW).
c/    snDTTotal     : Total beam-target DT neutron rate production (#/s)
c/    snDDTotal     : Total beam-target DD neutron rate production (#/s)
c/
c/    iflag         : Error flag. If not equal to zero, then a problem
c/                    was encountered.
c/
c/    ADDITIONAL OUTPUT VARIABLES:
c/
c/    The above output variable list represents only a small selection of
c/    the various calculated quantities in the module. Some of the other
c/    calculated quantities that may be of interest are included below.
c/    If they are needed, then they should be placed in the calling
c/    statement of this subroutine.
c/
c/    beamFusElecTot : Electron beam-target fusion power (MW).
c/    beamFusIonTot  : Ion beam-target fusion power (MW).
c/
c/    pitchangl(i,ie,ib) : Cosine of the average fast ion angle for
c/                         the ie energy component of beamline ib at grid
c/                         point i.
c/
c/    jnbie(i,ie,ib)  : Neutral beam driven current density by the
c/                      ie component of beam ib at grid point i (A/m^2).
c/    jnb(i,ib)       : Neutral beam driven current density due to
c/                      beamline ib at grid point i (A/m^2).
c/    beamDDFus(i)    : Beam-target fusion power density from DD
c/                      reactions at grid point i (MW/m^3).
c/    beamDTFus(i)    : Beam-target fusion power density
c/                      from DT reactions (MW/m^3).
c/    beamFusDTHe4(i) : Beam-target fusion power density to the
c/		        He4 particles of the DT reaction (MW/m^3).
c/    beamFusDDHe3(i) : Beam-target fusion power density to the
c/                      He3 particles of the DD reaction (MW/m^3).
c/    beamFusDDp(i)   : Beam-target fusion power density
c/  		        in the protons of the DD reaction (MW/m^3).
c/    beamFusDDt(i)   : Beam-target fusion power density in the
c/ 		        tritons of the DD reaction (MW/m^3).
c/
c/    About the radial grid:
c/    ---------------------
c/    The regular grid consists of n grid points (uniform or nonuniform)
c/    from the magnetic axis (rnorm(1) = 0) to the last closed flux surface
c/    (rnorm(n) = 1.0). All output profile quantities and input radial plasma
c/    profiles (ne, ni, Te, Ti, etc.) are supposed to be defined at these
c/    grid points. However, dvol, darea, and vprime are defined at the
c/    boundaries of a cell centered around each grid point, which can
c/    be considered as a "half-grid", i.e. r(i+1/2) = 0.5*(r(i) + r(i+1)).
c/    So, dvol(i) is the volume enclosed between r(i+1/2) and r(i-1/2).
c/    The calculations are carried out from i = 1, to i = n-1, so there is
c/    no need to define dvol(n) or darea(n). However, vprime(n) should be
c/    defined.
c/
c/    Historical Note (J. Mandrekas, 08/23/00):
c/    ----------------------------------------
c/    This module was originally developed by the author for the SuperCode
c/    (S.W. Haney, L.J. Perkins, J. Galambos and J. Mandrekas) a mixed
c/    C++/Fortran code running under a shell. In this code, variables were
c/    either Public or Private and were included in 'Module Descriptor Files'
c/    which were similar to Fortran Common Blocks. I have moved the private
c/    variables into the BeamsLoc.inc common block. Most of the 'public'
c/    variables are now included in the subroutine list, while a lot of
c/    background plasma and geometry related variables are in the common
c/    block 'nbplasma.inc'. To avoid rewriting large parts of the code,
c/    some global variables are copied to local versions which are then
c/    stored in one of the common blocks.

      implicit none

      include 'nbparams.inc'
      include 'BeamsLoc.inc'
      include 'nbplasma.inc'
      include 'nbconsts.inc'

c/    Subroutine Variable Declarations:
c/    --------------------------------

      integer iflag, nbeams, inbfus, n, nion, maxiter
      integer nbshape(maxBeams), nbptype(maxBeams)

      real amb, zbeam, r0, a, b0, volp, nbcurTot, etanbTot, beamBeta,
     .     pNBAbsorbTot, pNBLossTot, beamFusTot, beamFusChTot,
     .     snDDTotal, snDTTotal, taus

      real ebeam(maxBeams), pbeam(maxBeams), rtang(maxBeams),
     .     bwidth(maxBeams), bheigh(maxBeams), bgaussR(maxBeams),
     .     bgaussZ(maxBeams), bzpos(maxBeams), pwrfrac(3, maxBeams),
     .     aion(maxIons), zion(maxIons), ne20(mxrho),
     .     ni20(mxrho,maxIons), tekev(mxrho), tikev(mxrho), zeff(mxrho),
     .     rnorm(mxrho), vprime(mxrho), dvol(mxrho), darea(mxrho),
     .     kappa(mxrho), dkappa(mxrho), shafr(mxrho), dshafr(mxrho),
     .     hofr(mxrho,3,maxBeams), shinethru(3,maxBeams), jnbTot(mxrho),
     .     pnbe(mxrho), pnbi(mxrho), beamDens(mxrho), beamPress(mxrho),
     .     beamFus(mxrho), pbfuse(mxrho), pbfusi(mxrho),snBeamDD(mxrho),
     .     snBeamDT(mxrho), nbcur(maxBeams), etanb(maxBeams),
     .     gammanb(maxBeams), pNBAbsorb(maxBeams), pNBLoss(maxBeams)

c/    Local Variable Declarations:
c/    ---------------------------


c/    Variables in the FASTIONS subroutine list, but not exported by this
c/    routine:

      real jnbie(mxrho,3,maxBeams), jnb(mxrho,maxBeams),
     .     beamDDFus(mxrho), beamDTFus(mxrho), beamFusDTHe4(mxrho),
     .     beamFusDDHe3(mxrho), beamFusDDp(mxrho), beamFusDDt(mxrho)

c/    Other local variables:

      real beamFusElecTot, beamFusIonTot
      real beamFusIon(mxrho), beamFusElec(mxrho)

      real eDDHe3, eDTHe4, eDDp, eDDt, massHe3, massHe4, massp, masst,
     .   elecVolAvg, sumjnb, sumpnbi, sumpnbe

      real sdot
      integer i, j, ib, ie, isp, nrhom1

      data eDDHe3/820.00/, eDDp/3020.0/, eDDt/1010.0/, eDTHe4/3500.00/,
     .     massHe3/3.0/, massHe4/4.0/, massp/1.0/, masst/3.0/

      iflag = 0

c/    Check dimensions and return with error message if not right:

      if (nbeams.GT.maxBeams) then
	 iflag = 1
	 write (6,*) 'ERROR: Maximum beams dimension exceeded'
	 return
      endif

      if (nion.GT.maxIons) then
	 iflag = 2
	 write (6,*) 'ERROR: Maximum number of ions dimension exceeded'
	 return
      endif

      if (n.GT.mxrho) then
	 iflag = 3
	 write (6,*) 'ERROR: Radial grid points dimension exceeded'
	 return
      endif

c/    Assign variables in the subroutine statement to local variables:
c/    ---------------------------------------------------------------
      ab = amb
      maxit = maxiter

      do i = 1, nbeams
	     do ie = 1, 3
            fbpwr(ie,i) = pwrfrac(ie,i)
         enddo
      enddo

      nrho = n
      nrhom1 = nrho - 1
      rmajor = r0
      rminor = a
      volume = volp
      nions = nion

      do j = 1, nions
	 atw(j) = aion(j)
	 znum(j) = zion(j)
      enddo

c/    Notice: Below we use the BLAS routine scopy which copies one
c/    vector into another. If BLAS is not available (unlikely) then
c/    this part can be replaced by a simple do loop.
c/
      call scopy (n, rnorm, 1, rho, 1)
      call scopy (n, vprime, 1, vpr, 1)
      call scopy (n, dvol, 1, dv, 1)
      call scopy (n, kappa, 1, elong, 1)
      call scopy (n, dkappa, 1, elong_rho, 1)
      call scopy (n, shafr, 1, shift, 1)
      call scopy (n, dshafr, 1, shift_rho, 1)
      call scopy (n, ne20, 1, elecDensity, 1)
      call scopy (n, tekev, 1, elecTemp, 1)
      call scopy (n, tikev, 1, ionTemp, 1)
      call scopy (n, zeff, 1, zef, 1)
      do j = 1, nions
        call scopy(n, ni20(1,j), 1, dni(1,j), 1)
      enddo

c/    The following part of the code is needed to separate the ion
c/    species that are included in Boley's beam stopping cross section
c/    routine from those that are not and will be treated using
c/    Olson's formulation:

      do isp = 1, nions
	 inonb = 0
	 if(atw(isp).EQ.1.0.AND.znum(isp).EQ.1.) then
	    ihydr = isp
         else if (atw(isp).EQ.2.0.AND.znum(isp).EQ.1.) then
	    ideut = isp
         else if (atw(isp).EQ.3.0.AND.znum(isp).EQ.1.) then
	    itrit = isp
         else if (atw(isp).EQ.3.0.AND.znum(isp).EQ.2.) then
	    ihe3 = isp
	    inonb = inonb + 1
	    nonbindx(inonb) = isp
         else if (atw(isp).EQ.4.0.AND.znum(isp).EQ.2.) then
	    ihe4 = isp
         else if (atw(isp).EQ.12.0) then
	    icarb = isp
         else if (atw(isp).EQ.16.0) then
	    ioxy = isp
         else if (atw(isp).EQ.56.0) then
	    iFe = isp
         else
	    inonb = inonb + 1
	    nonbindx(inonb) = isp
         endif
      enddo

c/    Calculate fast ion deposition:
c/    -----------------------------
      call hofrho (nbeams, ebeam, pbeam, rtang, nbshape, bwidth,
     .   bheigh, nbptype, bgaussR, bgaussZ, bzpos, hofr, shinethru,
     .   pNBAbsorb, pNBLoss, pNBAbsorbTot, pNBLossTot)

c/    Calculate fast ion related parameters:
c/    -------------------------------------
      call fastIons (nbeams, amb, zbeam, ebeam, pbeam, inbfus,
     .   hofr, shinethru, jnbie, jnb, jnbTot, beamDens, beamPress,
     .   beamFus, beamDTFus, beamDDFus, beamFusDTHe4, beamFusDDHe3,
     .   beamFusDDp, beamFusDDt, snBeamDT, snBeamDD, taus)

c/    Calculate total beam heating power density (in MW/m^3):
c/    (summed over all beams and energy components)
c/    Recall that now hofr is normalized, so we must multiply
c/    by 1.0 - shinethru(ie,ib):
c/
      do i = 1, nrho
         sumpnbi = 0.0
         sumpnbe = 0.0
         do ib = 1, nbeams
            do ie = 1, 3
              sumpnbi = sumpnbi + pbeam(ib)*fbpwr(ie,ib)*
     .         (1.0 - shinethru(ie,ib)) * hofr(i,ie,ib)*fNBion(i,ie,ib)
              sumpnbe = sumpnbe + pbeam(ib)*fbpwr(ie,ib)*
     .         (1.0 - shinethru(ie,ib)) * hofr(i,ie,ib)*fNBelec(i,ie,ib)
            enddo
         enddo
         pnbi(i) = sumpnbi/volume
         pnbe(i) = sumpnbe/volume
      enddo

c/    Total driven current per beamline (A) and current drive
c/    efficiencies (summed over energy components);

      elecVolAvg = sdot(nrhom1,elecDensity,1,dvol,1) / volume

      nbcurTot = 0.0
      do ib = 1, nbeams
         sumjnb = 0.0
         do i = 1, nrhom1
            sumjnb = sumjnb + jnb(i,ib) * darea(i)
         enddo
         nbcur(ib) = sumjnb
         if (pNBAbsorb(ib).ne.0.0)
     .      etanb(ib) = 1.0e-06*nbcur(ib) / pNBAbsorb(ib)
         gammanb(ib) =  elecVolAvg*rmajor*etanb(ib)
         nbcurTot = nbcurTot + nbcur(ib)
      enddo

      if (pNBAbsorbTot.ne.0.0) then
         etanbTot = 1.0e-06*nbcurTot / pNBAbsorbTot
      endif

c/    Calculate the beam beta:
c/    -----------------------
      beamBeta = sdot(nrhom1, beamPress, 1,dvol, 1)
      beamBeta = 2.0*muZero*beamBeta / (volume*b0**2)

c/    Calculate beam-target related quantities (if inbfus=1):
c/    ----------------------------------------------------
      beamFusTot = 0.0
      beamFusElecTot = 0.0
      beamFusIonTot = 0.0
      beamFusChTot = 0.0
      snDDTotal = 0.0
      snDTTotal = 0.0

      if (inbfus.ne.0) then

c/    Calculate division of power to electrons and ions for charged
c/    particle products of the beam-plasma interactions:

         call eiSplit (eDDHe3, massHe3, fDDHe3ion, fDDHe3elec)
         call eiSplit (eDDp, massp, fDDpion, fDDpelec)
         call eiSplit (eDDt, masst, fDDtion,fDDtelec)
         call eiSplit (eDTHe4, massHe4, fDTHe4ion, fDTHe4elec)

         do i = 1, nrho
            beamFusIon(i) = beamFusDDHe3(i) * fDDHe3ion(i) +
     .                      beamFusDDp(i) * fDDpion(i) +
     .                      beamFusDDt(i) * fDDtion(i) +
     .                      beamFusDTHe4(i) * fDTHe4ion(i)
            beamFusElec(i) = beamFusDDHe3(i) * fDDHe3elec(i) +
     .                       beamFusDDp(i) * fDDpelec(i) +
     .                       beamFusDDt(i) * fDDtelec(i) +
     .                       beamFusDTHe4(i) * fDTHe4elec(i)
         enddo

c/   Calculate totals of selected quantities:

         beamFusTot = sdot(nrhom1,beamFus,1,dvol,1)
         beamFusElecTot = sdot(nrhom1,beamFusElec,1,dvol,1)
         beamFusIonTot = sdot(nrhom1,beamFusIon,1,dvol,1)
         beamFusChTot = beamFusElecTot + beamFusIonTot
         snDDTotal = sdot(nrhom1,snBeamDD,1,dvol,1)
         snDTTotal = sdot(nrhom1,snBeamDT,1,dvol,1)

      endif

      return
      end
