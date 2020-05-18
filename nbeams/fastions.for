      subroutine fastIons (nbeams, amb, zbeam, ebeam, pbeam, inbfus, 
     .   hofr, shinethru, jnbie, jnb, jnbTot, beamDens, beamPress, 
     .   beamFus, beamDTFus, beamDDFus, beamFusDTHe4, beamFusDDHe3, 
     .   beamFusDDp, beamFusDDt, snBeamDT, snBeamDD, taus)

c///////////////////////////////////////////////////////////////////////
c/
c/ NAME
c/    fastIons
c/
c/ DESCRIPTION
c/
c/    This subroutine calculates the current density driven by the
c/    neutral beams as a function of the flux surface label -rho-
c/    taking into account the  electron return current and the
c/    trapped electron correction, based on the uniform field solution
c/    of the Fokker-Planck equation for the fast ions.
c/    Written for the SUPERCODE by John Mandrekas, GIT
c/    Revisions:
c/    08-16-00, John Mandrekas
c/       Removed obsolete variables (deutPresent, etc.)
c/    06-05-92, John Mandrekas
c/       Updated it, to work with normalized deposition profiles.
c/    05-29-92, John Mandrekas
c/       Changed name to fastIons, to better reflect what the routine
c/       does.
c/    05-29-92, John Mandrekas
c/       Modified to calculate division of power between electrons and
c/       ions for the two-temperature implementation. Also, more
c/       resolution in the beam-target fusion reactions.
c/    04-15-92, John Mandrekas
c/       Modified to handle beams with three energy components.
c/       Changed definition of epsilon, to include shift.
c/    01-08-92, John Mandrekas
c/       Modified to calculate beam-target reaction rates and fast
c/       beta. 
c/    10-29-91, Original Release
c/
c/    References:
c/    ----------
c/    For the fast ion distribution function:
c/    J.D. Gaffey, J. Plasma Physics 16, 149 (1976)
c/    For the fit to the trapped electron correction:
c/    D.R. Mikkelsen, C.E. Singer, Nuclear Technol./Fusion 4, 237 (1983)
c/
c/    Calculated quantities:
c/    ---------------------
c/    jnb(i,ie,ib) : neutral beam driven current density from the 
c/                   ie component of the ib beam at grid point i (A/m^2)
c/    jnb(i,ib)    : neutral beam driven current density from all energy
c/                   components of beam ib at i (A/m^2)
c/    jnbTot(i)    : total neutral beam driven current density at i
c/                   (summed over all beamlines and energy groups) (A/m^2)
c/    beamFus(i)   : beam-target fusion power density (MW/m^3)
c/                   (this includes charged particles AND neutrons)
c/                   from all possible reactions and all beamlines and
c/                   energy groups.
c/    beamDTFus(i) : beam-target fusion power density from DT beam-target
c/                   interactions (MW/m^3)
c/    beamDDFus(i) : beam-target fusion power density from DD beam-target
c/                   interactions (MW/m^3)
c/
c/    beamFusDTHe4(i) : charged particle (He4) beam-target fusion power
c/                      density at i (MW/m^3)
c/    beamFusDDHe3(i) : charged particle (He3) beam-target fusion power
c/                      density at i  (MW/m^3)
c/    beamFusDDp(i)   : charged particle (p) beam-target fusion power
c/                      density at i (MW/m^3)
c/    beamFusDDt(i)   : charged particle (t) beam-target fusion power
c/                      density at i (MW/m^3)
c/    beamPress(i)    : fast ion pressure from the beam ions (Pa)
c/    snBeamDD(i)     : rate of production of DD neutrons (#/s/m^3)
c/    snBeamDT(i)     : rate of production of DT neutrons (#/s/m^3)
c/    beamDens(i)     : fast ion density (#/m^3)
c/
c///////////////////////////////////////////////////////////////////////

     
      implicit none
      
      include 'nbparams.inc'
      include 'BeamsLoc.inc'
      include 'nbplasma.inc'
      include 'nbconsts.inc'

      integer nbeams, inbfus
      real amb, zbeam, taus
      real ebeam(maxBeams), pbeam(maxBeams)
      real hofr(mxrho,3,maxBeams), jnbie(mxrho,3,maxBeams), 
     .     jnb(mxrho,maxBeams), shinethru(3,maxBeams), jnbTot(mxrho),
     .     beamDens(mxrho), beamPress(mxrho), beamFus(mxrho), 
     .     beamDDFus(mxrho), beamDTFus(mxrho), beamFusDTHe4(mxrho), 
     .     beamFusDDHe3(mxrho), beamFusDDp(mxrho), beamFusDDt(mxrho), 
     .     snBeamDD(mxrho), snBeamDT(mxrho) 
     
c/    Local variable declarations:
c/    ---------------------------
      real srcfast(mxrho,3,maxBeams), bmcur(3,maxBeams), 
     .     vbeam(3,maxBeams)
     
      real ami, bcoef, capke, cloge, clogi, comfact, dne20, energy,
     .     ecrit, eddn, eddp, edhe3, edtn, epsilon, factor,
     .     eddnhe3, eddpp, eddpt, edtnhe4,
     .     fb, frate, jfast, rcoef, reslt, rint, rrddn, rrddp,
     .     rrdhe3, rrdt, root3, sumbrack, sumhat, tauf, term,
     .     sumjnbib, sumsDDib, sumsDTib, sumFusib, sumdnbib, sumPreib,
     .     sumDDFusib, sumDTFusib, sumFusDTHe4, sumFusDDHe3, sumFusDDp,
     .     sumFusDDt, sumjnbie, term1, term2, trapfact, xcr, yc, zbr,
     .     zhat, zte, zti, zz
      integer i, ie, isp, nfrst, nrhom1, jhy, jdeut, jtrit, jhe3, jhe4,
     .        ib, ibion, ireact
      external fb, frate
      common/intg/yc, zhat
      common/sigma/ecrit, comfact, ireact
      save nfrst, jhy, jdeut, jtrit, jhe3, jhe4, ibion

      data nfrst/0/, jhy/0/, jdeut/0/, jtrit/0/, jhe3/0/, jhe4/0/,
     .     ibion/0/
     
c/    Fusion reaction energies in MeV:

      data eddn/3.27/, eddp/4.03/, edhe3/18.3/, edtn/17.6/, 
     .   eddnhe3/0.82/, eddpp/3.02/, eddpt/1.01/, edtnhe4/3.50/

c/    Identify ion species for the fusion reactivity calculations:
c/    Current formulation allows injection of D or He3 beams.
c/    However, the beam species should also be a part of the plasma
c/    ion population:

      if (nfrst.eq.0.and.inbfus.eq.1) then
	 if (ideut.ne.0) jdeut = ideut
	 if (itrit.ne.0) jtrit = itrit
	 if (ihe3.ne.0) jhe3  = ihe3
	 if (ihe4.ne.0) jhe4  = ihe4
         if (amb.eq.2.0.and.zbeam.eq.1.0) ibion = jdeut
         if (amb.eq.3.0.and.zbeam.eq.2.0) ibion = jhe3
         if (ibion.eq.0) then ! cancel fusion calculations
            inbfus = 0
         endif
         nfrst = 1
      endif

      root3 = Sqrt(3.0)

c/    Calculate the equivalent atomic current of each beamline in A:
c/    (Recall that pbeam is in MW and ebeam in keV!):

      do  ib = 1, nbeams
	 do ie = 1, 3
            bmcur(ie,ib) = 1.0e3*fbpwr(ie,ib)*float(ie)*
     .	       pbeam(ib)/ebeam(ib)
            vbeam(ie,ib) = 4.3770e5*Sqrt(ebeam(ib)/(float(ie)*amb))
         enddo
      enddo

c/    Initialize arrays to avoid problems with some compilers:

      do i = 1, nrho
         snBeamDD(i) = 0.0
         snBeamDT(i) = 0.0
         jnbTot(i) = 0.0
	 beamFus(i) = 0.0
	 beamDDFus(i) = 0.0
	 beamDTFus(i) = 0.0
	 beamFusDTHe4(i) = 0.0
	 beamFusDDHe3(i) = 0.0
	 beamFusDDp(i) = 0.0
	 beamFusDDt(i) = 0.0
	 beamDens(i) = 0.0
	 beamPress(i) = 0.0
	 do ib = 1, nbeams
	    jnb(i,ib) = 0.0
	    do ie = 1, 3
	       jnbie(i,ie,ib) = 0.0
	       srcfast(i,ie,ib) = 0.0
	       fNBion(i,ie,ib) = 0.0
	       fNBelec(i,ie,ib) = 0.0
            enddo
         enddo
      enddo

      nrhom1 = nrho - 1
      do  i = 1, nrhom1
	 epsilon = rho(i)*rminor/(rmajor + rminor*shift(i))
         dne20 = elecDensity(i)
	 zte = elecTemp(i)
	 zti = ionTemp(i)
         call coulomb (cloge,dne20,zte,0.0,amb,0.0,0.0,0.0,1)
         taus = 0.2*amb*zte**1.5/(cloge*dne20*zbeam**2)

c/    Calculate trapped electron correction factor:

         trapfact = (1.55 + 0.85/zef(i))*Sqrt(epsilon) -
     x              (0.20 + 1.55/zef(i))*epsilon
         factor  = 1.0 - (1.0 - trapfact)*zbeam/zef(i)
	 sumjnbib = 0.0
	 sumsDDib = 0.0
	 sumsDTib = 0.0
	 sumFusib = 0.0
	 sumDDFusib = 0.0
	 sumDTFusib = 0.0
	 sumFusDTHe4 = 0.0
	 sumFusDDHe3 = 0.0
	 sumFusDDp = 0.0
	 sumFusDDt = 0.0
	 sumdnbib = 0.0
	 sumPreib = 0.0
         do ib = 1,  nbeams
	    sumjnbie = 0.0
	    do ie = 1, 3
	       if (fbpwr(ie,ib).ne.0.0) then
                  energy = ebeam(ib)/float(ie)
                  sumbrack = 0.0
                  sumhat = 0.0
                  do  isp = 1, nions
                     ami = atw(isp)
                     zz  = znum(isp)
                     call coulomb (clogi,dne20,zte,ami,amb,zz,zbeam,
     .    	       energy,2)
                     sumbrack = sumbrack + dni(i,isp)*zz**2*clogi/ami
                     sumhat = sumhat + dni(i,isp)*zz*zz*clogi/amb
                  enddo
                  zbr = sumbrack/(cloge*elecDensity(i))
                  zhat = sumhat/sumbrack

c/    Critical energy for fast ion slowing down:

                  ecrit = 14.8*(zbr**(2./3.))*amb*zte
                  yc = Sqrt(ecrit/energy)
                  call qsimp(fb, 0., 1., reslt)
                  capke = reslt*(1. + yc**3)**(zhat/3.)
                  tauf = taus*Log((1.0+(energy/ecrit)**1.5)/
     x                   (1.0+(zti/ecrit)**1.5))/3.0

c/    Calculate fast and net current densities (A/m2):

                  srcfast(i,ie,ib) = bmcur(ie,ib) *
     x		    (1.0 - shinethru(ie,ib)) * hofr(i,ie,ib)/volume
                  jfast = zbeam*srcfast(i,ie,ib)*taus*vbeam(ie,ib)*
     x                    capke*pitchangl(i,ie,ib)
	          jnbie(i,ie,ib) = factor*jfast
                  sumjnbie = sumjnbie + jnbie(i,ie,ib) 

                  if (inbfus.eq.1) then
                  
c/    Calculate rates of beam-plasma interactions:

                     rcoef = 1.36593e13*srcfast(i,ie,ib)*taus/
     x                       Sqrt(amb)
                     if (ibion.eq.jdeut) then              ! D injection
                        if (jdeut.ne.0) then
                           comfact = 0.50
                           ireact = 1                      ! D(d,n)3He
                           call qsimp (frate,zti,energy,rint)
                           rrddn = rcoef*dni(i,jdeut)*rint
                           ireact = 2                      ! D(d,p)T
                           call qsimp (frate,zti,energy,rint)
                           rrddp = rcoef*dni(i,jdeut)*rint
                        endif
                        if (jhe3.ne.0) then                 ! D(3He,p)4He
                           comfact = 3.0/5.0
                           ireact = 3
                           call qsimp (frate,zti,energy,rint)
                           rrdhe3 = dni(i,jhe3)*rcoef*rint
                        endif
                        if (jtrit.ne.0) then                ! D(t,n)4He
                           comfact = 3.0/5.0
                           ireact = 4
                           call qsimp (frate,zti,energy,rint)
                           rrdt = dni(i,jtrit)*rcoef*rint
                        endif
                     else if (ibion.eq.jhe3) then          ! He3 injection
                        comfact = 2.0/5.0
                        ireact = 3
                        call qsimp (frate,zti,energy,rint)
                        rrdhe3 = dni(i,jdeut)*rcoef*rint
                     endif
                     sumsDDib = sumsDDib + rrddn
                     sumsDTib = sumsDTib + rrdt
		     sumDDFusib = sumDDFusib + 1.6022e-19 * 
     .                   (eddn*rrddn + eddp*rrddp)
		     sumDTFusib = sumDTFusib + 1.6022e-19 * 
     .                   edtn*rrdt
                     sumFusib = sumFusib + 1.6022e-19*(
     .                      eddn*rrddn + eddp*rrddp + edhe3*rrdhe3 +
     .                      edtn*rrdt)
		     sumFusDTHe4 = sumFusDTHe4 + 1.6022e-19*
     .                  rrdt*edtnhe4
		     sumFusDDHe3 = sumFusDDHe3 + 1.6022e-19*
     .                  rrddn*eddnhe3
		     sumFusDDp = sumFusDDp + 1.6022e-19*
     .                  rrddp * eddpp
		     sumFusDDt = sumFusDDt + 1.6022e-19*
     .                  rrddp * eddpt
		     
                  endif

c/    Calculate fast ion density (#/m^3):
                  sumdnbib = sumdnbib + srcfast(i,ie,ib)*tauf/eCharge

c/    Calculate fast ion pressure (in Pa):

                  xcr = Sqrt(energy/ecrit)
                  term1 = Log((xcr**2+2.*xcr+1.)/(xcr**2-xcr+1))/6.0
                  term2 = (atan((2.*xcr-1)/root3) + pi/6.)/root3
                  term = xcr**2/2. + term1 - term2
                  bcoef = 2.*ecrit*1.6022e-16*taus/3.
                  sumPreib = sumPreib +
     .                     srcfast(i,ie,ib)*bcoef*term/eCharge

c/    Calculate power division between electrons and ions:

		  fNBion(i,ie,ib) = 2.0*(term2 - term1) / (xcr**2)
		  fNBelec(i,ie,ib) = 1.0 - fNBion(i,ie,ib)
               endif
	    enddo

	    jnb(i,ib) = sumjnbie
	    sumjnbib = sumjnbib + sumjnbie
         enddo
	 snBeamDD(i) = sumsDDib
	 snBeamDT(i) = sumsDTib
	 jnbTot(i) = sumjnbib
	 beamFus(i)  = sumFusib
	 beamDTFus(i) = sumDTFusib
	 beamDDFus(i) = sumDDFusib
	 beamFusDTHe4(i) = sumFusDTHe4
	 beamFusDDHe3(i) = sumFusDDHe3
	 beamFusDDp(i) = sumFusDDp
	 beamFusDDt(i) = sumFusDDt
	 beamDens(i) = sumdnbib
	 beamPress(i) = sumPreib
      enddo

      return
      end


c/////////////////////////////////////////////////////////////////////////
c/
c/ NAME
c/    fb
c/
c/ DESCRIPTION
c/    Calculates the distribution function for the fast ions
c/
c/ PARAMETERS
c/    y   : normalized fast ion energy
c/
c/ RETURNS
c/    fb  : value of distribution function
c/
c/
c/////////////////////////////////////////////////////////////////////////

      real function fb(y)

      real y, yc, zhat
      common/intg/yc,zhat
      
      fb=(y**3 / (y**3 + yc**3))**(zhat / 3.0 + 1.0)
      
      return
      end


