c/    This program is a driver to test the NB module nbeams. It is not
c/    part of the nbeams package. Its purpose is to create a background
c/    plasma (geometry, MHD and plasma profiles) and then call the NB
c/    routines to calculate various neutral beam heating and current drive
c/    parameters.
c/
c/    The program reads the input file 'inbeams.dat' which contains a
c/    namelist with various input variables, calls the main routine of
c/    the NB module (calcBeams) and prints out some information to the
c/    file 'outbeams.dat'
c/
c/    Written by John Mandrekas, GIT, for the NTCC 08/29/00
c/
c/    Description of some of the input variables in the namelist
c/    (Only the local variables are included here. For a description
c/    of the variables in the calcBeams list see the extensive comments
c/    in the beginning of the calcbeams.f routine)
c/
c/    e0         - elongation at center
c/    ea         - elongation at edge
c/    shft0      - shift of magnetic axis (normalized to minor radius)
c/    teavg      - Density-weighted volume average electron temperature (keV)
c/    tea        - Electron temperature at the edge (separatrix) (keV)
c/    tiavg      - Density-weighted volume average ion temperature (keV)
c/    tia        - Ion temperature at the edge (separatrix) (keV)
c/    alfat      - Temperature profile exponent
c/    denav      - Average electron density (10^20 /m^3)
c/    edgavr     - edge-to-average ratio, ne(a)/<ne>
c/    alfan      - Density profile exponent
c/    fion       - Concentration of each ion species (ni/ne)

c/
        implicit none
        include 'nbparams.inc'
        include 'BeamsLoc.inc'

c/      Variables in the calcBeams variable list:
c/      ----------------------------------------

        integer iflag, ie, ib, nbeams, inbfus, n, nion, maxiter
        integer nbshape(maxBeams), nbptype(maxBeams)

        real amb, zbeam, r0, a, b0, volp, nbcurTot, etanbTot, beamBeta,
     .     pNBAbsorbTot, pNBLossTot, beamFusTot, beamFusChTot,
     .     snDDTotal, snDTTotal, taus

        real ebeam(maxBeams), pbeam(maxBeams), rtang(maxBeams),
     .     bwidth(maxBeams), bheigh(maxBeams), bgaussR(maxBeams),
     .     bgaussZ(maxBeams), bzpos(maxBeams), pwrfrac(3,maxBeams),
     .     aion(maxIons), zion(maxIons), ne20(mxrho),
     .     ni20(mxrho,maxIons), tekev(mxrho), tikev(mxrho), zeff(mxrho),
     .     rnorm(mxrho), vprime(mxrho), dvol(mxrho), darea(mxrho),
     .     kappa(mxrho), dkappa(mxrho), shafr(mxrho), dshafr(mxrho),
     .     hofr(mxrho,3,maxBeams), shinethru(3,maxBeams), jnbTot(mxrho),
     .     pnbe(mxrho), pnbi(mxrho), beamDens(mxrho), beamPress(mxrho),
     .     pbfuse(mxrho), pbfusi(mxrho), snBeamDD(mxrho),beamFus(mxrho),
     .     snBeamDT(mxrho), nbcur(maxBeams), etanb(maxBeams),
     .     gammanb(maxBeams), pNBAbsorb(maxBeams), pNBLoss(maxBeams)

c/      Local variables:

        integer i, j, ip1, nm1, nbinp, nbout
        real e0, te0, ti0, ea, shft0, dlr, tt, pi, den0, L_tor, sumzef, rn2, ssum
        CHARACTER(len=100) argi, argo
        real r(mxrho), rmid(mxrho)

c       namelist /nbin/ nbeams, inbfus, n, nion, maxiter, nbshape,
c     .    nbptype, amb, zbeam, ebeam, pbeam, rtang, bwidth, bheigh,
c     .    bgaussR, bgaussZ, bzpos, pwrfrac, aion, zion, r0, a, b0,
c     .    e0, ea, shft0, teavg, tiavg, tea, tia, alfat, denav, edgavr,
c     .    alfan, fion

        namelist /nbin/ nbeams, inbfus, n, nion, maxiter, nbshape,
     .    nbptype, amb, zbeam, ebeam, pbeam, rtang, bwidth, bheigh,
     .    bgaussR, bgaussZ, bzpos, pwrfrac, aion, zion, r0, a, b0,
     .    e0, ea, shft0, ne20, ni20, tekev, tikev

        pi = acos(-1.0)

c/      Assign I/O unit numbers:
c/      -----------------------
        nbinp = 10
        nbout = 20

c/      Open files:
c/      ----------
        CALL getarg(1, argi)
        open (unit=nbinp,file =TRIM(argi),status ='old')
        read (nbinp,nbin)
        close(unit=nbinp,status='keep')
        
        CALL getarg(2,argo)
        open (nbout,file =TRIM(argo),status ='unknown')

c/      BACKGROUND PLASMA:

c/      Grid and other geometric quantities:
c/      -----------------------------------
        nm1 = n - 1
        dlr = a / float(nm1)
        r(1) = 0.0
        r(n) = a
c/        open(unit=30,file='bugging.txt',status='old')
        do i = 2, nm1
            r(i) = r(i-1) + dlr
c/            write(30,*) r(i)
        enddo
c/        close(30,status='keep')

c/      Interface grid:
c/      --------------
c/      Notice that i+1/2 corresponds to i+1 in the interface grid!
        rmid(1) = 0.0
        do i = 2, n
            rmid(i) = 0.5 * (r(i-1) + r(i))
        enddo

c/      Neutral Beam normalized grid:

        do i = 1, n
            rnorm(i) = r(i) / r(n)
            rn2 = rnorm(i)**2
            tt = 1.0 - rn2
            kappa(i)  = e0 + (ea - e0)*rn2
            dkappa(i) = 2.0*(ea - e0) * rnorm(i)
            shafr(i)  = shft0 * tt
            dshafr(i) = -2.0 * shft0 * rnorm(i)
        enddo

c/      Cylindrical limit of  important metric  quantities:
c/      --------------------------------------------------

        L_tor = 2.0 * pi * r0
        vprime(1) = 0.0
        do i = 1, nm1
            ip1 = i + 1
            dvol(i) = L_tor * pi * kappa(i) * (rmid(ip1)**2 - rmid(i)**2)
            vprime(ip1) = L_tor * 2.0 * pi * kappa(i) * rmid(ip1)
            darea(i) = dvol(i) / L_tor
        enddo

        dvol(n) = 0.0
        darea(n) = 0.0

c/      Plasma volumes:
c/      --------------
        volp = ssum(nm1, dvol, 1)
c/      Parabolic-like plasma profiles:
c/      ------------------------------
c       te0 = ((1.0 + alfat+alfan) / (1.0 + alfan - alfan*alfat*edgavr / (1.0 + alfat))) * teavg
c       ti0 = te0 * tiavg / teavg
c       den0 = (1.0 + alfan - alfan*edgavr) * denav

        te0 = tekev(1)
        ti0 = tikev(1)
        den0 = ne20(1)

        do  i = 1, n
c            tt = 1.0 - (r(i) / a)**2
c	        tekev(i) = (te0 - tea) * tt**alfat + tea
c           tikev(i) = (ti0 - tia) * tt**alfat + tia
c           ne20(i) = (den0 - edgavr*denav) * tt**alfan + edgavr*denav
            do j = 2, nion
                ni20(i,j) = (ne20(i) - ni20(i,j-1)*zion(j-1))/zion(j)
            enddo
        enddo

c/      Calculation of the plasma Zeff:
c/      ------------------------------

        do i = 1, n
            sumzef = 0.0
            do j = 1, nion
                sumzef = sumzef + ni20(i,j) * zion(j)**2 / ne20(i)
            enddo
            zeff(i) = sumzef
        enddo

c/      Call the beams package:
        print *,'nbeams: ',nbeams
        print *,'amb: ',amb
        print *,'zbeam: ',zbeam
        print *,'ebeam: ',ebeam
        print *,'pbeam: ',pbeam
        print *,'inbfus: ',inbfus
        print *,'rtang: ',rtang
        print *,'nbshape: ',nbshape
        print *,'bwidth: ',bwidth
        print *,'bheigh: ',bheigh
        print *,'nbptype: ',nbptype
        print *,'bgaussR: ', bgaussR
        print *,'bgaussZ: ',bgaussZ
        print *,'bzpos: ', bzpos
        print *,'pwrfrac: ', pwrfrac
        print *,'maxiter: ', maxiter
        print *,'nion: ', nion
        print *,'aion: ', aion
        print *,'zion: ', zion
        print *,'ne20: ',ne20
        print *,'ni20: ', ni20
        print *,'tekev: ',   tekev
        print *,'tikev: ', tikev
        print *,'zeff: ', zeff
        print *,'r0: ', r0
        print *,'a: ', a
        print *,'b0: ', b0
        print *,'volp: ', volp
        print *,'n: ', n
        print *,'rnorm: ',rnorm
        print *,'vprime: ',vprime
        print *,'dvol: ', dvol
        print *,'darea: ',darea
        print *,'kappa: ', kappa
        print *,'dkappa: ', dkappa
        print *,'shafr: ', shafr
        print *,'dshafr: ', dshafr



        call calcBeams(nbeams, amb, zbeam, ebeam, pbeam, inbfus,
     .   rtang, nbshape, bwidth, bheigh, nbptype, bgaussR, bgaussZ,
     .   bzpos, pwrfrac, maxiter, nion, aion, zion, ne20, ni20, tekev,
     .   tikev, zeff, r0, a, b0, volp, n, rnorm, vprime, dvol, darea,
     .   kappa, dkappa, shafr, dshafr, hofr, shinethru, jnbTot, pnbe,
     .   pnbi, beamDens, beamPress, beamFus, pbfuse, pbfusi, snBeamDD,
     .   snBeamDT, nbcur, etanb, gammanb, pNBAbsorb, pNBLoss, nbcurTot,
     .   etanbTot, beamBeta, pNBAbsorbTot, pNBLossTot, beamFusTot,
     .   beamFusChTot, snDTTotal, snDDTotal, iflag,taus)

        print *,'hofr: ', hofr
        print *,'shinethru: ', shinethru
        print *,'jnbTot: ', jnbTot
        print *,'pnbe: ',pnbe
        print *,'pnbi: ',pnbi
        print *,'beamDens: ', beamDens
        print *,'beamPress: ', beamPress
        print *,'beamFus: ', beamFus
        print *,'pbfuse: ', pbfuse
        print *,'pbfusi: ', pbfusi
        print *,'snBeamDD: ',snBeamDD
        print *,'snBeamDT: ', snBeamDT
        print *,'nbcur: ', nbcur
        print *,'etanb: ', etanb
        print *,'gammanb: ', gammanb
        print *,'pNBAbsorb: ', pNBAbsorb
        print *,'pNBLoss: ', pNBLoss
        print *,'nbcurTot: ',nbcurTot
        print *,'etanbTot: ', etanbTot
        print *,'beamBeta: ',beamBeta
        print *,'pNBAbsorbTot: ', pNBAbsorbTot
        print *,'pNBLossTot: ', pNBLossTot
        print *,'beamFusTot: ',beamFusTot
        print *,'beamFusChTot: ', beamFusChTot
        print *,'snDTTotal: ', snDTTotal
        print *,'snDDTotal: ', snDDTotal
        print *,'iflag: ',iflag
        print *,'taus: ',taus
c/      Write output quantities to file nbout:

c/      Global parameters
c/      -----------------
        write (nbout, 2000)
        write (nbout, 2100)
        write (nbout, 2200) pNBAbsorbTot, pNBLossTot, nbcurTot, etanbTot, beamBeta, taus, volp
        write (nbout, '(1x)')

        write (nbout, 2300)
        write (nbout, 2400)
        write (nbout, 2410) (pNBAbsorb(ib), ib = 1, nbeams)
        write (nbout, 2415) (pNBLoss(ib), ib = 1, nbeams)
        write (nbout, 2420) (nbcur(ib), ib = 1, nbeams)
        write (nbout, 2430) (etanb(ib), ib = 1, nbeams)
        write (nbout, 2440) (gammanb(ib), ib = 1, nbeams)
        write (nbout, '(1x)')
        write (nbout, 2450)
        write (nbout, 2460) (shinethru(1,ib), ib = 1, nbeams)
        write (nbout, 2470) (shinethru(2,ib), ib = 1, nbeams)
        write (nbout, 2480) (shinethru(3,ib), ib = 1, nbeams)
        write (nbout, '(1x)')

        if (inbfus.NE.0) then
            write (nbout, 2500)
            write (nbout, 2600) beamFusTot, beamFusChTot, snDTTotal, snDDTotal
        endif

        write (nbout, '(1x)')

c/      Write the deposition profile for each beamline and energy group:
c/      ---------------------------------------------------------------
        write (nbout, *) 'Neutral Beam Deposition Profiles'
        write (nbout, *) '--------------------------------'
        do ib = 1, nbeams
            write (nbout, 1000) ib
            do i = 1, n
                write (nbout,1100) rnorm(i), (hofr(i,ie,ib),ie=1,3)
            enddo
            write (nbout, '(1x)')
        enddo

c/      Write the pitch angle profile for each beamline and energy group:
c/      ---------------------------------------------------------------
        write (nbout, *) 'Pitch Angle at Birth'
        write (nbout, *) '--------------------------------'
        do ib = 1, nbeams
            write (nbout, 1200) ib
            do i = 1, n
                write (nbout,1300) rnorm(i), (pitchangl(i,ie,ib),ie=1,3)
            enddo
            write (nbout, '(1x)')
        enddo

c/      Write heating and current drive profile information:
c/      ----------------------------------------------------
        write (nbout, 3000)
        do i = 1, n
            write (nbout, 3100) rnorm(i), jnbTot(i), pnbe(i), pnbi(i), beamDens(i), beamPress(i), beamFus(i), dvol(i),darea(i)
        enddo


 1000 format (1x, 'Beam ', i2/
     .   5x, 'rho', 6x, 'hofr_1', 7x, 'hofr_2', 7x, 'hofr_3')
 1100 format (1x, f8.5, 3(1x,e12.5))
 1200 format (1x, 'Angle ', i2/
     .   5x, 'rho', 6x, 'zeta_1', 7x, 'zeta_2', 7x, 'zeta_3')
 1300 format (1x, f8.5, 3(1x,e12.5))

 2000 format (1x, 'Global NB Heating and Current Drive Parameters',/
     .        1x, '----------------------------------------------')
 2100 format (1x, 'TOTALS:')
 2200 format (1x, 'Total Absorbed Power      = ', f9.4, ' MW'/
     .        1x, 'Total Lost Power          = ', f9.4, ' MW'/
     .        1x, 'Total NB Driven Current   = ', e12.5,' A'/
     .        1x, 'Total NBCD Efficiency     = ', e12.5,' A/W'/
     .        1x, 'Total Beam Beta           = ', e12.5, /
     .        1x, 'Taus                      = ',  e12.5, ' s'/
     .        1x, 'Volp                      = ',  e12.5, ' m3' )

 2300 format (1x, 'Information per Beamline and Energy Group:'/
     .        1x, '-----------------------------------------')
 2400 format (23x,'Beam 1', 10x, 'Beam 2', 10x, 'Beam 3'/
     .        23x,'------', 10x, '------', 10x, '------')
 2410 format (1x, 'Absorbed Power   ', 1x, e12.5, 2(4x, e12.5))
 2415 format (1x, 'Lost Power       ', 1x, e12.5, 2(4x, e12.5))
 2420 format (1x, 'NB driven current', 1x, e12.5, 2(4x, e12.5))
 2430 format (1x, 'NBCD efficiency  ', 1x, e12.5, 2(4x, e12.5))
 2440 format (1x, 'NBCD gamma       ', 1x, e12.5, 2(4x, e12.5))
 2450 format (1x, 'Beam Shinethrough')
 2460 format (1x, '   energy group 1', 1x, e12.5, 2(4x, e12.5))
 2470 format (1x, '   energy group 2', 1x, e12.5, 2(4x, e12.5))
 2480 format (1x, '   energy group 3', 1x, e12.5, 2(4x, e12.5))

 2500 format (1x, 'Beam-Target Parameters'/
     .        1x, '----------------------')
 2600 format (1x, 'Total Beam-Target Fusion Power   = ', e12.5, ' MW'/
     .        1x, 'Total Power to Charged Particles = ', e12.5, ' MW'/
     .        1x, 'Total DT Neutron Rate            = ', e12.5, ' n/s'/
     .        1x, 'Total DD Neutron Rate            = ', e12.5, ' n/s')

 3000 format (3x, 'rho', 5x, 'jnbtot', 6x, 'pNBe', 8x, 'pNBi', 8x,
     .   'nbfast', 6x, 'pressb', 6x, 'pfusb', 6x, 'dvol',6x,'dA')
 3100 format (1x, f6.4, 8(1x, e11.4))

c/      close(nbout,status='unknown')
      close(nbout)
      stop
      end
