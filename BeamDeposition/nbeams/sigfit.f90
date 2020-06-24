subroutine sigfit(ebeam, dene, te, denz, sig)
    !-----------------------------------------------------------------------
    !.....author:
    !       Charles D. Boley
    !       Lawrence Livermore National Laboratory
    !       L-574
    !       Livermore, CA 94550
    !       (415) 423-7365

    !     version of 3/11/92

    !.....This subroutine computes the neutral beam stopping cross section,
    !     for a hydrogenic beam, as a function of beam energy, electron
    !     density, electron temperature, and impurity density.  Four
    !     impurities (He, C, O, Fe) are allowed.

    !     The code employs a fit based on results of the author's beam
    !     penetration code.  The latter code contains the detailed atomic
    !     cross sections recommended at the IAEA Specialists' Meeting on
    !     Required Atomic Data Base for Neutral Beam Penetration in Large
    !     Tokamaks (Vienna, 4/89).  The fit applies for beam energies from
    !     10 keV/amu to 10 MeV/amu and for all densities and electron
    !     temperatures of interest.

    !     The fit is independent of the ion temperature and hence does not
    !     distinguish among hydrogen ion isotopes in the plasma.  (The actual
    !     cross sections were evaluated at Ti = 15 keV.)  The ion temperature
    !     has little effect provided that it is small compared to the beam
    !     energy per amu.  (This condition normally is a triviality, except
    !     perhaps for the 1/3 energy component of a deuterium beam.)  The
    !     following table gives an idea of the variation with ion temperature
    !     at low beam energy, for injection into a pure H plasma with
    !     dene=1.e14 and Te=10 keV (energies in keV; cross sections in this
    !     table in units of 10**-16 cm**2):

    !     ebeam  sig(Ti=1)  sig(Ti=15)  sig(Ti=30)  sig(fit)
    !       10     11.43       9.96        8.91      11.18
    !       30      6.06       5.12        4.71       5.35
    !       50      3.57       3.56        3.39       3.75
    !      100      2.09       2.11        2.09       2.28

    !     The fit was evaluated with B(perpendicular component) = 5 T.
    !     Variations of the magnetic field have an insignificant effect on
    !     the stopping cross section.

    !     Accuracy of fit: rms error about 2.5%
    !                      max error about 12%

    !     Input parameters --
    !       ebeam: beam energy per amu (keV/amu) -- 10. .le. ebeam .le. 1.e4
    !       dene:  electron density (cm**-3) -- 1.e12 .le. dene .le. 1.e15
    !       te:    electron temperature (keV) -- te .ge. 1.
    !       denz:  impurity densities (cm**-3) -- array of dimension 4
    !         denz(1): He
    !         denz(2): C
    !         denz(3): O
    !         denz(4): Fe
    !         Note: denz(i)=0. is permissible
    !     Output parameter --
    !       sig:   beam stopping cross section (cm**2)
    !
    !     Modified by John Mandrekas (GIT) for the SUPERCODE
    !     New version with fits valid from lower beam energies (10 keV)
    !     03/27/92, John Mandrekas, GIT

    implicit none
    integer mz, mt1, mt2, mt3
    parameter(mz = 4, mt1 = 4, mt2 = 3, mt3 = 2)
    integer i, nth1, nth2, nth3, ntz1(mz), ntz2(mz), ntz3(mz), &
            iinit
    real A1(mt1, mt2, mt3), A2(mt1, mt2, mt3, mz)
    common/cfit/A1, nth1, nth2, nth3, A2, ntz1, ntz2, ntz3
    real s1, ebeam, te, sig, sum, dene, poly3f, &
            denz(mz), s2(mz), az(mz), xx(3)
    data az/2., 6., 8., 26./
    data iinit/0/
    if(iinit.eq.0) then
        iinit = 1
        call initfit
    endif

    xx(1) = alog(ebeam)
    xx(2) = alog(dene / 1.e13)
    xx(3) = alog(te)
    s1 = poly3f(A1, mt1, mt2, mt3, xx(1), nth1, xx(2), nth2, xx(3), nth3)
    do i = 1, 4
        if(denz(i).gt.0.) then
            s2(i) = poly3f(A2(1, 1, 1, i), mt1, mt2, mt3, xx(1), ntz1(i), &
                    xx(2), ntz2(i), xx(3), ntz3(i))
        else
            s2(i) = 0.
        endif
    enddo

    sum = 0.0
    do i = 1, 4
        sum = sum + (denz(i) * az(i) * (az(i) - 1.) / dene) * s2(i)
    enddo

    sig = (1.e-16 / ebeam) * exp(s1) * (1. + sum)

    return
end

subroutine initfit
    !-----------------------------------------------------------------------
    !-----New version, with different fits valid from lower beam energies
    !-----Received from C. Boley, LLNL
    !-----Modified for  use by the SuperCode, John Mandrekas, 03/27/92

    implicit none
    integer mz, mt1, mt2, mt3
    parameter(mz = 4, mt1 = 4, mt2 = 3, mt3 = 2)
    integer i, nth1, nth2, nth3, ntz1(mz), ntz2(mz), ntz3(mz)
    real A1(mt1, mt2, mt3), A2(mt1, mt2, mt3, mz)
    common/cfit/A1, nth1, nth2, nth3, A2, ntz1, ntz2, ntz3

    !..   pure plasma
    nth1 = 3
    nth2 = 3
    nth3 = 2
    A1(1, 1, 1) = 3.95e+00
    A1(1, 1, 2) = 1.60e-02
    A1(1, 2, 1) = -3.84e-02
    A1(1, 2, 2) = -5.98e-03
    A1(1, 3, 1) = -3.10e-03
    A1(1, 3, 2) = -1.09e-03
    A1(2, 1, 1) = 3.67e-01
    A1(2, 1, 2) = -2.15e-02
    A1(2, 2, 1) = 3.07e-02
    A1(2, 2, 2) = 1.78e-03
    A1(2, 3, 1) = 3.16e-03
    A1(2, 3, 2) = 3.47e-04
    A1(3, 1, 1) = -9.95e-03
    A1(3, 1, 2) = 6.19e-04
    A1(3, 2, 1) = -2.36e-03
    A1(3, 2, 2) = -1.67e-04
    A1(3, 3, 1) = -1.31e-04
    A1(3, 3, 2) = -2.28e-05

    do i = 1, 4
        ntz1(i) = 4
        ntz2(i) = 2
        ntz3(i) = 2
    enddo

    !..   He
    i = 1
    A2(1, 1, 1, i) = -1.76e+00
    A2(1, 1, 2, i) = -2.90e-01
    A2(1, 2, 1, i) = 5.43e-02
    A2(1, 2, 2, i) = 1.04e-02
    A2(2, 1, 1, i) = 7.49e-01
    A2(2, 1, 2, i) = 1.57e-01
    A2(2, 2, 1, i) = -3.74e-02
    A2(2, 2, 2, i) = -6.70e-03
    A2(3, 1, 1, i) = -7.28e-02
    A2(3, 1, 2, i) = -2.45e-02
    A2(3, 2, 1, i) = 6.11e-03
    A2(3, 2, 2, i) = 1.14e-03
    A2(4, 1, 1, i) = 2.05e-03
    A2(4, 1, 2, i) = 1.34e-03
    A2(4, 2, 1, i) = -2.82e-04
    A2(4, 2, 2, i) = -5.42e-05

    !..   C
    i = 2
    A2(1, 1, 1, i) = -1.89e-01
    A2(1, 1, 2, i) = -3.22e-02
    A2(1, 2, 1, i) = 5.43e-02
    A2(1, 2, 2, i) = 3.83e-03
    A2(2, 1, 1, i) = -1.41e-03
    A2(2, 1, 2, i) = 8.98e-03
    A2(2, 2, 1, i) = -3.34e-02
    A2(2, 2, 2, i) = -1.97e-03
    A2(3, 1, 1, i) = 3.10e-02
    A2(3, 1, 2, i) = 7.43e-04
    A2(3, 2, 1, i) = 5.08e-03
    A2(3, 2, 2, i) = 1.96e-04
    A2(4, 1, 1, i) = -2.54e-03
    A2(4, 1, 2, i) = -4.21e-05
    A2(4, 2, 1, i) = -2.20e-04
    A2(4, 2, 2, i) = 8.32e-07

    !..   O
    i = 3
    A2(1, 1, 1, i) = -1.07e-01
    A2(1, 1, 2, i) = -3.36e-02
    A2(1, 2, 1, i) = 4.90e-02
    A2(1, 2, 2, i) = 3.77e-03
    A2(2, 1, 1, i) = -3.72e-02
    A2(2, 1, 2, i) = 1.11e-02
    A2(2, 2, 1, i) = -2.89e-02
    A2(2, 2, 2, i) = -1.84e-03
    A2(3, 1, 1, i) = 3.34e-02
    A2(3, 1, 2, i) = 7.20e-06
    A2(3, 2, 1, i) = 4.14e-03
    A2(3, 2, 2, i) = 1.64e-04
    A2(4, 1, 1, i) = -2.50e-03
    A2(4, 1, 2, i) = 1.08e-05
    A2(4, 2, 1, i) = -1.65e-04
    A2(4, 2, 2, i) = 2.45e-06

    !..   Fe
    i = 4
    A2(1, 1, 1, i) = 5.46e-02
    A2(1, 1, 2, i) = -3.89e-02
    A2(1, 2, 1, i) = 1.71e-02
    A2(1, 2, 2, i) = -7.35e-04
    A2(2, 1, 1, i) = -8.97e-02
    A2(2, 1, 2, i) = 2.36e-02
    A2(2, 2, 1, i) = -7.46e-03
    A2(2, 2, 2, i) = 9.94e-04
    A2(3, 1, 1, i) = 2.96e-02
    A2(3, 1, 2, i) = -4.71e-03
    A2(3, 2, 1, i) = 2.52e-04
    A2(3, 2, 2, i) = -3.05e-04
    A2(4, 1, 1, i) = -1.75e-03
    A2(4, 1, 2, i) = 3.63e-04
    A2(4, 2, 1, i) = 4.10e-05
    A2(4, 2, 2, i) = 2.37e-05

    return
end

real function poly3f(a, ndim1, ndim2, ndim3, &
        x1, n1, x2, n2, x3, n3)

    !-----Modified by j. mandrekas, for the SuperCode
    !-----New version with new fits, 03/27/92

    implicit none

    integer ndim1, ndim2, ndim3, n1, n2, n3, i, i1, i2, i3
    real a(ndim1, ndim2, ndim3)
    real x1, x2, x3, y1(4), y2(4), y3(4)

    y1(1) = 1.0
    do i = 2, n1
        y1(i) = y1(i - 1) * x1
    enddo

    y2(1) = 1.0
    do i = 2, n2
        y2(i) = y2(i - 1) * x2
    enddo

    y3(1) = 1.0
    do i = 2, n3
        y3(i) = y3(i - 1) * x3
    enddo

    poly3f = 0.0
    do i1 = 1, n1
        do i2 = 1, n2
            do i3 = 1, n3
                poly3f = poly3f + a(i1, i2, i3) * y1(i1) * y2(i2) * y3(i3)
            enddo
        enddo
    enddo

    return
end