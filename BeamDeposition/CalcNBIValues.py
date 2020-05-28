#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from shapely.geometry import Point, LineString
from scipy.interpolate import UnivariateSpline
import numpy as np

class calc_nbi_vals:
    def __init__(self, inp, core,):
        R0 = core.pts.axis.mag[0]


        R_tan = inp.rtang

        # center of central solenoid is at the origin
        # define vertical line at x=R_tan going from -y to +y, where
        # y = sqrt(R_tan^2 - (R_0(a) + a)^2)

        origin = Point(0,0)
        pt1 = Point(R_tan, -1*sqrt((core.R0_a + core.a)**2 - R_tan**2))
        pt2 = Point(R_tan, sqrt((core.R0_a + core.a)**2 - R_tan**2))

        beamline = LineString((pt1, pt2))
        # interpolate along half of the beam path (the second half will give identical values as the first half)
        # use a hundred points or so, just to get a smooth interpolation
        rho_beam = np.zeros(100)
        zeta_beam = np.zeros(100)
        for i,v in enumerate(np.linspace(0,0.5,100)):
            pt = beamline.interpolate(v, normalized=True)
            R = pt.distance(origin)
            rho_beam[i] = abs(R - R0) / core.a
            zeta_beam[i] = R_tan / R

        rho_zeta = np.column_stack((rho_beam,zeta_beam))
        rho_zeta_sorted = rho_zeta[rho_zeta[:,0].argsort()]

        #
        # A couple things need to be noted here:
        #   1.  If the radius of tangency of a beam is less than the R of the magnetic axis,
        #       then those inner flux surfaces that the beam hits four times will have two
        #       distinct zeta values for the IOL calculation. It hits the inboard side of the
        #       flux surface with a different direction cosine than when it hits the outbard
        #       side of the flux surface. This results in a bifurcated plot of zeta vs rho.
        #       You can see it with the following code:
        #
        #       plt.scatter(rho_zeta_sorted[:,0], rho_zeta_sorted[:,1])
        #
        #   2.  Although in principle, you could include the correct zeta values for each rho
        #       in the beam IOL calculation, the resulting F_orb, M_orb, and E_orb functions are
        #       not going to be terribly different from what they would be if used a single,
        #       appropriately averaged zeta value for the entire beam. This is partially because
        #       all of the zeta values are likely to be fairly high (i.e. ~> 0.8).
        #
        #   3.  If someone down the road wants to treat this more carefully, knock yourself out.
        #       At present, a single, deposition profile-weighted zeta will be calculated and used
        #       for the beam IOL calculation. If you change this, be sure to update these notes.
        #

        def moving_average_sortof(x, y, deltax):
            xmin = np.amin(x)
            xmax = np.amax(x)

            xnew = np.linspace(xmin,xmax,100)
            ynew = np.zeros(100)
            for i,xval in enumerate(xnew):
                ymax = np.amax(y[np.logical_and(x < xval+deltax, x > xval-deltax)])
                ymin = np.amin(y[np.logical_and(x < xval+deltax, x > xval-deltax)])
                ynew[i] = (ymax + ymin)/2
            return xnew, ynew

        rho_av_prof, zeta_av_prof = moving_average_sortof(rho_zeta_sorted[:,0], rho_zeta_sorted[:,1], 0.01)
        dPdr_norm1_interp = UnivariateSpline(inp.dPdr_norm1_data[:,0], inp.dPdr_norm1_data[:,1], k=3, s=0)

        zeta_av_val = np.sum(zeta_av_prof * dPdr_norm1_interp(rho_av_prof))

        # calculat dVdr, which is used in the calculation of dPdV
        dVdr = core.rho2vol.derivative()(core.rho)
        dVdr_1D = core.rho2vol.derivative()(core.rho[:,0])

        #Set final values for beam 1
        self.zeta_1 = zeta_av_val

        # assume dPdr_norm1 given with rho, like most other inputs
        self.dPdr_1 = inp.pbeam * dPdr_norm1_interp(core.rho)
        self.dPdr_1_1D = inp.pbeam * dPdr_norm1_interp(core.rho[:,0])

        self.dPdV_1 = self.dPdr_1 / dVdr
        self.dPdV_1_1D = self.dPdr_1_1D / dVdr_1D

        # TODO: add ability to have more than one beam when specifying deposition profiles as inputs.
        self.zeta_2 = zeta_av_val

        self.dPdr_2 = np.zeros(core.rho.shape)
        self.dPdr_2_1D = np.zeros(core.rho[:, 0].shape)

        self.dPdV_2 = np.zeros(core.rho.shape)
        self.dPdV_2_1D = np.zeros(core.rho[:, 0].shape)

        self.zeta_3 = zeta_av_val

        self.dPdr_3 = np.zeros(core.rho.shape)
        self.dPdr_3_1D = np.zeros(core.rho[:, 0].shape)

        self.dPdV_3 = np.zeros(core.rho.shape)
        self.dPdV_3_1D = np.zeros(core.rho[:, 0].shape)

        pwr_frac = np.array([1, 0, 0])
        self.beam_pwr_1 = inp.pbeam * pwr_frac[0]  # atomic deuterium
        self.beam_pwr_2 = inp.pbeam * pwr_frac[1]  # molecular deuterium D2
        self.beam_pwr_3 = inp.pbeam * pwr_frac[2]  # molecular deuterium D3

        self.beam_en_1 = inp.ebeam
        self.beam_en_2 = inp.ebeam / 2
        self.beam_en_3 = inp.ebeam / 3