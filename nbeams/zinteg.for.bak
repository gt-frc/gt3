c/////////////////////////////////////////////////////////////////////////
c/
c/ NAME 
c/    Beams_zinteg
c/
c/ DESCRIPTION
c/    This subroutine calculates the Z-part of the 2D integral
c/    for the calculation of the beam deposition profile.
c/    Written by John Mandrekas, GIT, for the SuperCode
c/    Revisions:
c/    04-14-92, John Mandrekas:
c/       Updated to handle three beams with three energy components 
c/    10-25-91, John Mandrekas:
c/       Original release
c/
c/ PARAMETERS
c/    None
c/
c/ RETURNS
c/    Nothing
c/
c/////////////////////////////////////////////////////////////////////////

      subroutine zinteg
      
      implicit none
      include 'nbparams.inc'
      include 'BeamsLoc.inc'
      include 'nbconsts.inc'
      

c/    Local variable declarations:
c/    ---------------------------
      integer ie, iz
      real  dterm, eterm, el2, el3, sumz1, sumz2, sumz3, zb2, zu, 
     .      zterm1, zterm, ztest, resltz, sumzksi1, sumzksi2, sumzksi3

c/    Limiting case for rho -> 0.0:

      if (rhoi.le.0.0) then
         zb = 0.0
         call rinteg
         zterm  = 0.5*pi*elong0
         do ie = 1, 3
            alongz(ie) = alongr(ie)*zterm
            avksi(ie)  = alongksi(ie)*zterm
         enddo
         return
      endif

c/    Integration limits in the z direction:

      ztest = rhoi*elongr
      if (shape.eq.0) zu = min (ztest, rbeam)
      if (shape.eq.1) zu = min (ztest, 0.5*height)
      
      dterm = dshft/rhoi
      eterm = delong/rhoi
      el2   = elongr**2
      el3   = elongr**3
      
c/    Start integration:      

      sumz1 = 0.0
      sumz2 = 0.0
      sumz3 = 0.0
      sumzksi1 = 0.0
      sumzksi2 = 0.0
      sumzksi3 = 0.0
      do iz  = 1, nz
         zb  = 0.5*zu*(1.0 + sz(iz))
         zb2 = zb**2
         call rinteg
         zterm1 = (1.0 + eterm*zb2/el3)/
     .           Sqrt(rhoi**2 - zb2/el2) + sgn*dterm
         zterm = zterm1*Sqrt(1.0 - sz(iz)**2)
         sumz1 = sumz1 + alongr(1)*zterm
         sumz2 = sumz2 + alongr(2)*zterm
         sumz3 = sumz3 + alongr(3)*zterm
         sumzksi1 = sumzksi1 + alongksi(1)*zterm
         sumzksi2 = sumzksi2 + alongksi(2)*zterm
         sumzksi3 = sumzksi3 + alongksi(3)*zterm
      enddo
     
      resltz = 0.5*zu*wz
      alongz(1) = resltz*sumz1
      alongz(2) = resltz*sumz2
      alongz(3) = resltz*sumz3
      avksi(1)  = resltz*sumzksi1
      avksi(2)  = resltz*sumzksi2
      avksi(3)  = resltz*sumzksi3
      
      return
      end
