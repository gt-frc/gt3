c/////////////////////////////////////////////////////////////////////////
c/    
c/     NAME
c/        frate
c/    
c/     DESCRIPTION
c/    
c/        This function calculates the integrand needed in the calculation
c/    of the beam-plasma fusion reaction rates. The fusion cross sections
c/    are taken from the recent work of H. Bosch, IPP I/252, 1990 and
c/    are expressed in mbarns. The cross sections are given as functions
c/    of the COM energy (in keV). The Lab energy -e- is converted to
c/    the COM energy, -ecm- by multiplying with the factor -comfact-
c/    which is evaluated in the calling routine.
c/    Created by Andrew Renalds, GIT, 01-23-91
c/    Modified for use in the SuperCode by John Mandrekas, GIT, 01-07-92.
c/    Information regarding the different coefficients for the different
c/    energy regions of the D(3He,p)4He, and D(t,n)4He reactions moved
c/    from the calling routine (fastIons), here to hide unecessary details
c/    and fix problems at low energies, JM, 01-11-92
c/
c/ PARAMETERS
c/    e  : LAB energy in keV
c/
c/ RETURNS
c/    frate: the integrand for the velocity space integration
c/
c/////////////////////////////////////////////////////////////////////////

      real function frate(e)
      
      implicit none

      real a1(6), a2(6), a3(6), a4(6), a5(6), bg(6)
      real b1(6), b2(6), b3(6), b4(6)
      real ecrit, comfact, e, ecm, edtmax, sn, sd, s, sigv, dfterm,
     .     dhe3lim, dtlim, fresult
      integer ii, ireact
      common/sigma/ecrit, comfact, ireact


c/    e     : energy in the lab frame, in keV
c/    ecm   : energy in the COM frame, in keV
c/    ecrit : critical energy, in keV
c/    ii    : reaction index as follows:
c/      1 = coefficients for the D(d,n)3He   reaction 0.5-4900keV
c/      2 = coefficients for the D(d,p)T     reaction 0.5-5000keV
c/      3 = coefficients for the D(3He,p)4He reaction 0.3-900keV
c/      4 = coefficients for the D(3He,p)4He reaction 900-4800keV
c/      5 = coefficients for the D(t,n)4He   reaction 0.5-550keV
c/      6 = coefficients for the D(t,n)4He   reaction 550-4700keV

      data a1(1),a2(1),a3(1) / 5.3701e4,3.3027e2,-1.2706e-1 /
      data a4(1),a5(1) / 2.9327e-5,-2.5151e-9 /
      data b1(1),b2(1) / 0.0, 0.0 /
      data b3(1),b4(1) / 0.0, 0.0 /
      data bg(1) / 31.3970 /

      data a1(2),a2(2),a3(2) / 5.5576e4,2.1054e2,-3.2638e-2 /
      data a4(2),a5(2) / 1.4987e-6,1.8181e-10 /
      data b1(2),b2(2) / 0.0, 0.0 /
      data b3(2),b4(2) / 0.0, 0.0 /
      data bg(2) / 31.3970 /

      data a1(3),a2(3),a3(3) / 5.7501e6,2.5226e3,4.5566e1 /
      data a4(3),a5(3) / 0.0, 0.0 /
      data b1(3),b2(3) / -3.1995e-3,-8.5530e-6 /
      data b3(3),b4(3) / 5.9014e-8, 0.0 /
      data bg(3) / 68.7508 /

      data a1(4),a2(4),a3(4) / -8.3993e5, 0.0, 0.0 /
      data a4(4),a5(4) / 0.0, 0.0 /
      data b1(4),b2(4) / -2.6830e-3,1.1633e-6 /
      data b3(4),b4(4) / -2.1332e-10,1.4250e-14 /
      data bg(4) / 68.7508 /

      data a1(5),a2(5),a3(5) / 6.927e4,7.454e8,2.05e6 /
      data a4(5),a5(5) / 5.2002e4, 0.0 /
      data b1(5),b2(5) / 63.8,-0.995 /
      data b3(5),b4(5) / 6.981e-5,1.728e-4 /
      data bg(5) / 34.3827 /

      data a1(6),a2(6),a3(6) / -1.4714e6, 0.0, 0.0 /
      data a4(6),a5(6) / 0.0, 0.0 /
      data b1(6),b2(6) / -8.4127e-3,4.7983e-6 /
      data b3(6),b4(6) / -1.0748e-9,8.5184e-14 /
      data bg(6) / 34.3827 /

      data edtmax /4700.00/, dhe3lim/900.0/, dtlim/550.0/

      ecm = comfact*e                ! convert to COM energy
c-----check energy regions:
      if (ireact.eq.1) ii = 1        ! D(d,n)he3
      if (ireact.eq.2) ii = 2        ! D(d,p)T
      if (ireact.eq.3) then          ! D(he3,p)he4
         if (ecm.le.dhe3lim) ii = 3
         if (ecm.gt.dhe3lim) ii = 4
      endif
      if (ireact.eq.4) then          ! D(t,n)he4
         if (ecm.le.dtlim) ii = 5
         if (ecm.gt.dtlim) ii = 6
         if(ecm.gt.edtmax) ecm = edtmax
      endif

c-----Calculate cross section and integrand:
      sn = a1(ii)+ecm*(a2(ii)+ecm*(a3(ii)+ecm*(a4(ii)+ecm*a5(ii))))
      sd = 1.0+ecm*(b1(ii)+ecm*(b2(ii)+ecm*(b3(ii)+ecm*b4(ii))))
      s = sn/sd
      sigv = s/(ecm*exp(bg(ii)/Sqrt(ecm)))         ! in mbarns!
      dfterm = 1.0/(Sqrt(e)*(1.0+(ecrit/e)**1.5))
      fresult = sigv*dfterm

      frate = fresult
      
      return
      end
