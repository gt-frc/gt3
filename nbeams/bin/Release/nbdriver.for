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

c/    Variables in the calcBeams variable list:
c/    ----------------------------------------
      
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
	 
	

c/    Local variables:
      
      integer i, j, ip1, nm1, nbinp, nbout
      real e0, te0, ti0, ea, shft0, dlr, tt, pi,
	.   den0, L_tor, sumzef, rn2, ssum

      real r(mxrho), rmid(mxrho)

c      namelist /nbin/ nbeams, inbfus, n, nion, maxiter, nbshape, 
c     .    nbptype, amb, zbeam, ebeam, pbeam, rtang, bwidth, bheigh,
c     .    bgaussR, bgaussZ, bzpos, pwrfrac, aion, zion, r0, a, b0,
c     .    e0, ea, shft0, teavg, tiavg, tea, tia, alfat, denav, edgavr, 
c     .    alfan, fion

      namelist /nbin/ nbeams, inbfus, n, nion, maxiter, nbshape, 
     .    nbptype, amb, zbeam, ebeam, pbeam, rtang, bwidth, bheigh,
     .    bgaussR, bgaussZ, bzpos, pwrfrac, aion, zion, r0, a, b0,
     .    e0, ea, shft0

      pi = acos(-1.0)

c/    Assign I/O unit numbers:
c/    -----------------------
      nbinp = 10
      nbout = 20

c/    Open files:
c/    ----------
      open (unit=nbinp,file ='inbeams.dat',status ='old')
      read (nbinp,nbin)
	close(unit=nbinp,status='keep')

	open (nbout,file ='outbeams.dat',status ='unknown')
      
c/    BACKGROUND PLASMA:

c/    Grid and other geometric quantities:
c/    -----------------------------------
      nm1 = n - 1
      dlr = a / float(nm1)
      r(1) = 0.0
	r(n) = a
	open(unit=30,file='bugging.txt',status='old')
      do i = 2, nm1
          r(i) = r(i-1) + dlr
	    write(30,*) r(i)
	enddo
	close(30,status='keep')
	
c/    Interface grid:
c/    --------------
c/    Notice that i+1/2 corresponds to i+1 in the interface grid!
      rmid(1) = 0.0
      do i = 2, n
        rmid(i) = 0.5 * (r(i-1) + r(i))
      enddo

c/    Neutral Beam normalized grid:

      do i = 1, n
         rnorm(i) = r(i) / r(n)
         rn2 = rnorm(i)**2
         tt = 1.0 - rn2
         kappa(i)  = e0 + (ea - e0)*rn2                                         
         dkappa(i) = 2.0*(ea - e0) * rnorm(i)                                
         shafr(i)  = shft0 * tt
         dshafr(i) = -2.0 * shft0 * rnorm(i)
      enddo
         
c/    Cylindrical limit of  important metric  quantities:
c/    --------------------------------------------------

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
	
c/    Plasma volumes:
c/    --------------
      volp = ssum(nm1, dvol, 1)
c/    Parabolic-like plasma profiles:
c/    ------------------------------
c      te0 = ((1.0 + alfat+alfan) / (1.0 + alfan - alfan*alfat*edgavr /
c     .   (1.0 + alfat))) * teavg
c      ti0 = te0 * tiavg / teavg
c      den0 = (1.0 + alfan - alfan*edgavr) * denav
     tekev(1)  =1.35834934112
     tekev(2)  =1.35156016814
     tekev(3)  =1.34281705212
     tekev(4)  =1.33219322804
     tekev(5)  =1.31976193091
     tekev(6)  =1.30559639572
     tekev(7)  =1.28976985748
     tekev(8)  =1.27235555117
     tekev(9)  =1.2534267118
     tekev(10)  =1.23305657437
     tekev(11)  =1.21131837388
     tekev(12)  =1.18828534531
     tekev(13)  =1.16403072368
     tekev(14)  =1.13862774397
     tekev(15)  =1.11214964119
     tekev(16)  =1.08466965034
     tekev(17)  =1.05626100641
     tekev(18)  =1.0269969444
     tekev(19)  =0.996950699308
     tekev(20)  =0.966195506134
     tekev(21)  =0.934804599874
     tekev(22)  =0.902851215528
     tekev(23)  =0.870408588093
     tekev(24)  =0.837549952567
     tekev(25)  =0.804348543949
     tekev(26)  =0.770877597235
     tekev(27)  =0.737210347426
     tekev(28)  =0.703420029518
     tekev(29)  =0.669579878509
     tekev(30)  =0.635763129399
     tekev(31)  =0.602043017184
     tekev(32)  =0.568492776862
     tekev(33)  =0.535185643433
     tekev(34)  =0.502194851894
     tekev(35)  =0.469593637243
     tekev(36)  =0.437455234478
     tekev(37)  =0.405852878597
     tekev(38)  =0.374859804598
     tekev(39)  =0.34454924748
     tekev(40)  =0.31499444224
     tekev(41)  =0.286268623877
     tekev(42)  =0.258445027388
     tekev(43)  =0.231596887772
     tekev(44)  =0.205797440027
     tekev(45)  =0.18111991915
     tekev(46)  =0.15763756014
     tekev(47)  =0.135423597996
     tekev(48)  =0.114551267714
     tekev(49)  =0.0950938042933
     tekev(50)  =0.0771244427318
     tekev(51)  =0.0607164180275
     tikev(1)  =1.42084741506
     tikev(2)  =1.39972456608
     tikev(3)  =1.37749530755
     tikev(4)  =1.35422269888
     tikev(5)  =1.32996979948
     tikev(6)  =1.30479966873
     tikev(7)  =1.27877536605
     tikev(8)  =1.25195995085
     tikev(9)  =1.22441648252
     tikev(10)  =1.19620802046
     tikev(11)  =1.16739762408
     tikev(12)  =1.13804835279
     tikev(13)  =1.10822326598
     tikev(14)  =1.07798542307
     tikev(15)  =1.04739788344
     tikev(16)  =1.01652370651
     tikev(17)  =0.98542595168
     tikev(18)  =0.954167678351
     tikev(19)  =0.922811945926
     tikev(20)  =0.891421813808
     tikev(21)  =0.860060341402
     tikev(22)  =0.82879058811
     tikev(23)  =0.797675613335
     tikev(24)  =0.766778476481
     tikev(25)  =0.73616223695
     tikev(26)  =0.705889954145
     tikev(27)  =0.676024687471
     tikev(28)  =0.646629496328
     tikev(29)  =0.617767440122
     tikev(30)  =0.589501578255
     tikev(31)  =0.56189497013
     tikev(32)  =0.53501067515
     tikev(33)  =0.508911752718
     tikev(34)  =0.483661262238
     tikev(35)  =0.459322263113
     tikev(36)  =0.435957814745
     tikev(37)  =0.413630976537
     tikev(38)  =0.392404807894
     tikev(39)  =0.372342368218
     tikev(40)  =0.353506716912
     tikev(41)  =0.335960913379
     tikev(42)  =0.319768017022
     tikev(43)  =0.304991087245
     tikev(44)  =0.29169318345
     tikev(45)  =0.279937365041
     tikev(46)  =0.269786691421
     tikev(47)  =0.261304221992
     tikev(48)  =0.254553016159
     tikev(49)  =0.249596133324
     tikev(50)  =0.24649663289
     tikev(51)  =0.24531757426
     ne20(1)  =0.468303825572
     ne20(2)  =0.459208458359
     ne20(3)  =0.45041915764
     ne20(4)  =0.441922421074
     ne20(5)  =0.43370474632
     ne20(6)  =0.425752631036
     ne20(7)  =0.418052572882
     ne20(8)  =0.410591069515
     ne20(9)  =0.403354618594
     ne20(10)  =0.396329717778
     ne20(11)  =0.389502864726
     ne20(12)  =0.382860557096
     ne20(13)  =0.376389292547
     ne20(14)  =0.370075568737
     ne20(15)  =0.363905883325
     ne20(16)  =0.357866733971
     ne20(17)  =0.351944618331
     ne20(18)  =0.346126034066
     ne20(19)  =0.340397478833
     ne20(20)  =0.334745450292
     ne20(21)  =0.329156446101
     ne20(22)  =0.323616963919
     ne20(23)  =0.318113501403
     ne20(24)  =0.312632556214
     ne20(25)  =0.307160626009
     ne20(26)  =0.301684208448
     ne20(27)  =0.296189801188
     ne20(28)  =0.290663901889
     ne20(29)  =0.285093008209
     ne20(30)  =0.279463617806
     ne20(31)  =0.27376222834
     ne20(32)  =0.26797533747
     ne20(33)  =0.262089442853
     ne20(34)  =0.256091042148
     ne20(35)  =0.249966633014
     ne20(36)  =0.24370271311
     ne20(37)  =0.237285780094
     ne20(38)  =0.230702331625
     ne20(39)  =0.223938865362
     ne20(40)  =0.216981878963
     ne20(41)  =0.209817870087
     ne20(42)  =0.202433336392
     ne20(43)  =0.194814775537
     ne20(44)  =0.186948685181
     ne20(45)  =0.178821562983
     ne20(46)  =0.170419906601
     ne20(47)  =0.161730213693
     ne20(48)  =0.152738981919
     ne20(49)  =0.143432708936
     ne20(50)  =0.133797892405
     ne20(51)  =0.123821029982
     ni20(1,1)  =0.462313523153
     ni20(2,1)  =0.450387339946
     ni20(3,1)  =0.439060846317
     ni20(4,1)  =0.428308912692
     ni20(5,1)  =0.418106409496
     ni20(6,1)  =0.408428207155
     ni20(7,1)  =0.399249176093
     ni20(8,1)  =0.390544186737
     ni20(9,1)  =0.382288109511
     ni20(10,1)  =0.374455814842
     ni20(11,1)  =0.367022173154
     ni20(12,1)  =0.359962054873
     ni20(13,1)  =0.353250330425
     ni20(14,1)  =0.346861870235
     ni20(15,1)  =0.340771544728
     ni20(16,1)  =0.33495422433
     ni20(17,1)  =0.329384779467
     ni20(18,1)  =0.324038080563
     ni20(19,1)  =0.318888998044
     ni20(20,1)  =0.313912402336
     ni20(21,1)  =0.309083163864
     ni20(22,1)  =0.304376153053
     ni20(23,1)  =0.299766240329
     ni20(24,1)  =0.295228296118
     ni20(25,1)  =0.290737190844
     ni20(26,1)  =0.286267794934
     ni20(27,1)  =0.281794978812
     ni20(28,1)  =0.277293612904
     ni20(29,1)  =0.272738567636
     ni20(30,1)  =0.268104713433
     ni20(31,1)  =0.26336692072
     ni20(32,1)  =0.258500059924
     ni20(33,1)  =0.253479001468
     ni20(34,1)  =0.248278615779
     ni20(35,1)  =0.242873773282
     ni20(36,1)  =0.237239344403
     ni20(37,1)  =0.231350199567
     ni20(38,1)  =0.2251812092
     ni20(39,1)  =0.218707243726
     ni20(40,1)  =0.211903173572
     ni20(41,1)  =0.204743869162
     ni20(42,1)  =0.197204200923
     ni20(43,1)  =0.189259039279
     ni20(44,1)  =0.180883254657
     ni20(45,1)  =0.172051717481
     ni20(46,1)  =0.162739298177
     ni20(47,1)  =0.15292086717
     ni20(48,1)  =0.142571294886
     ni20(49,1)  =0.131665451751
     ni20(50,1)  =0.120178208189
     ni20(51,1)  =0.108084434627
	  te0 = tekev(1)
	  ti0 = tikev(1)
	  den0 = ne20(1)

      do  i = 1, n
	  tt = 1.0 - (r(i) / a)**2
c	  tekev(i) = (te0 - tea) * tt**alfat + tea
c        tikev(i) = (ti0 - tia) * tt**alfat + tia
c        ne20(i) = (den0 - edgavr*denav) * tt**alfan + 
c     .     edgavr*denav
	do j = 2, nion
	   ni20(i,j) = (ne20(i) - ni20(i,j-1)*zion(j-1))/zion(j)
	enddo
      enddo
	
c/    Calculation of the plasma Zeff:
c/    ------------------------------
      
      do i = 1, n
         sumzef = 0.0
         do j = 1, nion
            sumzef = sumzef + ni20(i,j) * zion(j)**2 / ne20(i)
         enddo
         zeff(i) = sumzef
      enddo

c/    Call the beams package:

      call calcBeams(nbeams, amb, zbeam, ebeam, pbeam, inbfus,      
     .   rtang, nbshape, bwidth, bheigh, nbptype, bgaussR, bgaussZ,       
     .   bzpos, pwrfrac, maxiter, nion, aion, zion, ne20, ni20, tekev,    
     .   tikev, zeff, r0, a, b0, volp, n, rnorm, vprime, dvol, darea,     
     .   kappa, dkappa, shafr, dshafr, hofr, shinethru, jnbTot, pnbe,     
     .   pnbi, beamDens, beamPress, beamFus, pbfuse, pbfusi, snBeamDD,    
     .   snBeamDT, nbcur, etanb, gammanb, pNBAbsorb, pNBLoss, nbcurTot,   
     .   etanbTot, beamBeta, pNBAbsorbTot, pNBLossTot, beamFusTot,        
     .   beamFusChTot, snDTTotal, snDDTotal, iflag,taus)                       
                                                                          
c/    Write output quantities to file nbout:
	
 
c/    Global parameters
c/    -----------------
      write (nbout, 2000)
      write (nbout, 2100)
      write (nbout, 2200) pNBAbsorbTot, pNBLossTot, nbcurTot,
     .   etanbTot, beamBeta, taus, volp
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
         write (nbout, 2600) beamFusTot, beamFusChTot, snDTTotal, 
     .      snDDTotal
      endif
 
      write (nbout, '(1x)')

c/    Write the deposition profile for each beamline and energy group:
c/    ---------------------------------------------------------------
      write (nbout, *) 'Neutral Beam Deposition Profiles'
      write (nbout, *) '--------------------------------'
      do ib = 1, nbeams
         write (nbout, 1000) ib
         do i = 1, n
            write (nbout,1100) rnorm(i), (hofr(i,ie,ib),ie=1,3)
         enddo  
         write (nbout, '(1x)')
      enddo

c/    Write the pitch angle profile for each beamline and energy group:
c/    ---------------------------------------------------------------
      write (nbout, *) 'Pitch Angle at Birth'
      write (nbout, *) '--------------------------------'
      do ib = 1, nbeams
         write (nbout, 1200) ib
         do i = 1, n
            write (nbout,1300) rnorm(i), (pitchangl(i,ie,ib),ie=1,3)
         enddo
         write (nbout, '(1x)')
      enddo

c/    Write heating and current drive profile information:
c/    ----------------------------------------------------
      write (nbout, 3000)
      do i = 1, n
         write (nbout, 3100) rnorm(i), jnbTot(i), pnbe(i), pnbi(i), 
     .	    beamDens(i), beamPress(i), beamFus(i), dvol(i),darea(i)
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
     
	 close(nbout,status='unknown')
      stop 
      end