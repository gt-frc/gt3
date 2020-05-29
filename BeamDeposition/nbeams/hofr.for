      subroutine hofrho (nbeams, ebeam, pbeam, rtang, nbshape, bwidth,
     .   bheigh, nbptype, bgaussR, bgaussZ, bzpos, hofr, shinethru,
     .   pNBAbsorb, pNBLoss, pNBAbsorbTot, pNBLossTot)

c/////////////////////////////////////////////////////////////////////////
c/
c/ NAME
c/    hofrho
c/
c/ DESCRIPTION
c/    This is the main routine for the calculation of the NB
c/    deposition profile using the diffuse beam model. This routine
c/    calls the appropriate routines for the R-Z integration on the
c/    face of the beam.
c/    Written by John Mandrekas, GIT, for the SuperCode
c/    Version 0.00, 08/08/91
c/    Modified for NTCC, 08/29/00, John Mandrekas
c/    References:
c/
c/    J.A. Rome, J.D. Callen, J.F. Clarke, Nucl. Fusion, 14, 141 (1974)
c/    H.C. Howe, 'Physics Models in the Tokamak Transport Code PROCTR',
c/                ORNL report, ORNL/TM-9537, 1985
c/    J. Mandrekas, 'Physics Models and User's Guide for the Neutral
c/                   Beam Module of the SuperCode,' GTFR-102, August 1992
c/
c/////////////////////////////////////////////////////////////////////////
c/
c/    Input parameters (notice, that other parameters needed by this 
c/    routine are passed by the BeamsLoc.inc and nbplasma.inc common
c/    blocks.
c/
c/    nbeams      : Number of beamlines
c/    ebeam(ib)   : Beam energy of each beamline (keV)
c/    pbeam(ib)   : Beam power of each beamline  (MW)
c/    rtang(ib)   : Tangency radius of each beamline (m)
c/    nbshape(ib) : Flag determining the cross section shape of each
c/                  beamline. Available optios are 0 for circlular
c/                  and 1 for rectangular.
c/    bwidth(ib)  : Width (or diameter for circular) of each beamline (m)
c/    bheigh(ib)  : Height of each beamline (m)
c/    nbptype(ib) : Flag for power distribution in the beam cross section.
c/                  0 : uniform distribution
c/                  1 : Gaussian (circular) or bi-Gaussian (rectangular)
c/    bgaussR(ib) : Gaussian width (radial direction for circular beams,
c/                  horizontal direction for rectangular beams). (m)
c/    bgaussZ(ib) : Gaussian width in vertical direction for rectangular
c/                  beams (m)
c/    bzpos(ib)   : Position of peak in Z direction for bi-Gaussian
c/                  power distribution (rectangular shapes only) (m)
c/ 
c/    Calculated Quantities:
c/    ---------------------      
c/    hofr(i,ie,ib)    : Neutral beam deposition profile for beam ib,
c/                       energy component ie, at grid point i. Notice
c/                       that the deposition profile is normalized to 1!
c/    shinethru(ie,ib) : Beam shinethrough for beam ib and energy
c/                       component ie.
c/    pNBLoss(ib)      : Neutral beam power lost per beamline (MW).
c/    pNBAbsorb(ib)    : Neutral beam power absorbed per beamline (MW).
c/    pNBAbsorbTot     : Total NB power absorbed in the plasma.
c/    pNBLossTot       : Total NB power lost.


      implicit none
      include 'nbparams.inc'
      include 'nbplasma.inc'
      include 'BeamsLoc.inc'
      include 'nbconsts.inc'

      integer nbeams
      integer nbshape(maxBeams), nbptype(maxBeams)
      real pNBAbsorbTot, pNBLossTot
      real ebeam(maxBeams), pbeam(maxBeams), rtang(maxBeams), 
     .     bwidth(maxBeams), bheigh(maxBeams), bgaussR(maxBeams), 
     .     bgaussZ(maxBeams), bzpos(maxBeams), pNBAbsorb(maxBeams),
     .     pNBLoss(maxBeams)
      real hofr(mxrho, 3, maxBeams), shinethru(3, maxBeams)

c/    Local variable declarations:
c/    ---------------------------
      integer i, nrhom1, ie, ib
      real erf
      real cnorm, fact, fact2, pi2, rtest1, rtest2, rlimit, rhoMax,
     .     prterm, pzterm1, pzterm2, pnbtotal, sumhofr, sumlossib
      real factr(3), hplus(3), hmnus(3), ksiplus(3), ksiminus(3)

      save rhoMax
      
      data rhoMax /0.99/

      pi2 = pi**2
      
c/    Calculate parameters for Gaussian integrations:
c/    ----------------------------------------------
      call gaussWts (nz, wz, sz)
      call gaussWts (nr, wr, sr)
      call gaussWts (ns, ws, ss)
     
      elong0 = elong(1) 
      elonga = elong(nrho)
      nrhom1 = nrho - 1
      rhoSurf = rminor*rhoMax

c/    Start do loop over beamlines:
c/    ----------------------------
      do ib = 1, nbeams
         ibeam = ib
         ebkev = ebeam(ib)
         rtan  = rtang(ib)
         shape = nbshape(ib)
         width = bwidth(ib)
         rbeam = 0.5 * width
         height= bheigh(ib)
         gaussR= bgaussR(ib)
         gaussZ= bgaussZ(ib)
         zpos  = bzpos(ib)
         ptype = nbptype(ib)
	
c/    Calculate beam stopping cross sections:
c/    --------------------------------------
         call Beams_mfp
      
c/    Normalization factor for current distribution:
c/    ---------------------------------------------
c/    Case for circular beam:
	
         if (shape.eq.0) then
            if (ptype.eq.0) cnorm = 1.0/(pi*rbeam**2)
            if (ptype.eq.1) then
               fact  = rbeam/gaussR
               fact2 = fact**2      
               cnorm = pi*gaussR**2*(1.0 - exp(-fact2))
               cnorm = 1.0/cnorm
             endif
c/    Case for rectangular beam:

         else if (shape.eq.1) then
            if (ptype.eq.0) cnorm = 1.0/(height*width)
            if (ptype.eq.1) then
               prterm = 0.5*width/gaussR
               pzterm1= 0.5*(height-2.0*zpos)/gaussZ
               pzterm2= 0.5*(height+2.0*zpos)/gaussZ
               cnorm = pi*gaussR*gaussZ*erf(prterm)*(
     .              erf(pzterm1) + erf(pzterm2))
               cnorm = 1.0/cnorm
            endif
         endif
      
c/    Getting ready for the big DO loop:
c/    ---------------------------------
         rtest1 = rtan - rbeam
         do iflux = 1, nrhom1
            rhoi  = rminor*rho(iflux)
            elongr = elong(iflux)
            delong = elong_rho(iflux)/rminor
            shft   = shift(iflux)*rminor
            dshft  = shift_rho(iflux)
            rtest2 = rmajor + shft + rhoi
           
c/    Check to see if there is intersection with this flux surface:

            if (rtest1.ge.rtest2) then
               do ie = 1, 3
                  hofr (iflux,ie,ib) = 0.0
                  pitchangl(iflux,ie,ib) = 0.0
               enddo
            else
            
c/    Limit of rho/vprime at the center:                     

               if (rhoi.eq.0.0) then
                  rlimit = 4.0*pi2*elong0*(rmajor + rminor*shift(1))
                  do ie = 1, 3
                     factr(ie) = 2.0*volume*omfp(1,ie,ib)/rlimit
                  enddo
               else
                  do ie = 1, 3
                     factr(ie) = (2.0*rhoi*volume / vpr(iflux)) *
     .                     omfp(iflux,ie,ib)
                  enddo
               endif
           
c/    Calculate contributions from outside intersections:

               sgn = 1.0
               call zinteg
               do ie = 1, 3
                  hplus(ie) = alongz(ie)
                  ksiplus(ie) = avksi(ie)
               enddo
	
c/    Calculate contributions from inside intersections:
               
               sgn = -1.0
               call zinteg
               do ie = 1, 3
                  hmnus(ie) = alongz(ie)
                  ksiminus(ie) = avksi(ie)
               enddo
 
c/    Total HOFR and average pitch angle at flux zone -iflux- :
	
               do ie = 1, 3
                  if (fbpwr(ie,ib).ne.0.0) then
                     hofr(iflux,ie,ib) = cnorm*factr(ie)*
     .                  (hplus(ie) + hmnus(ie))
                     pitchangl(iflux,ie,ib) = (ksiplus(ie) +
     .                 ksiminus(ie))/ (hplus(ie) + hmnus(ie))
                  else
                     hofr(iflux,ie,ib) = 0.0
                     pitchangl(iflux,ie,ib) = 0.0
                  endif
               enddo
            endif
         enddo

c/    Set deposition profile at the edge = 0
	
         do ie = 1, 3
            hofr(nrho,ie,ib) = 0.0
            pitchangl(nrho,ie,ib) = 0.0
         enddo
      enddo
	
c/    Calculate beam shinethrough, and other quantities of interest:
c/    -------------------------------------------------------------
      do ib = 1, nbeams
         do ie = 1, 3
            sumhofr = 0.0
            do i = 1, nrho
               sumhofr = sumhofr + hofr(i,ie,ib)*dv(i)
            enddo
            shinethru(ie,ib) = 1.0 - sumhofr/volume
         enddo
      enddo

c/    Renormalize the deposition profile, so that its volume integral
c/    divided by the total plasma volume is equal to one:
c/    ---------------------------------------------------
      do ib = 1, nbeams
         do ie = 1, 3
            if (shinethru(ie,ib).ne.1.0) then
              do i = 1, nrho
                hofr(i,ie,ib) = hofr(i,ie,ib) / (1.0 - shinethru(ie,ib))
              enddo
            endif
            if (shinethru(ie,ib).lt.0.0) shinethru(ie,ib) = 0.0
         enddo
      enddo
	
c/    Calculate the total NB power loss (shinethrough or missed
c/    the plasma, which is also part of shinethrough in this model).
c/    -------------------------------------------------------------
      pNBLossTot = 0.0
      pnbtotal =  0.0
      do ib = 1, nbeams
         sumlossib = 0.0
         pnbtotal = pnbtotal + pbeam(ib) 
         do ie = 1, 3
            sumlossib = sumlossib + pbeam(ib) * fbpwr(ie,ib) 
     .        * shinethru(ie,ib)
         enddo
         pNBLoss(ib) = sumlossib
         pNBAbsorb(ib) = pbeam(ib) - pNBLoss(ib)
         pNBLossTot = pNBLossTot + pNBLoss(ib)
      enddo
      pNBAbsorbTot = pnbtotal - pNBLossTot

      return
      end                
