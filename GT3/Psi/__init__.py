#!/usr/bin/python
from scipy.interpolate import griddata, interp1d

from GT3.utilities.PlotBase import PlotBase
from GT3.utilities.GT3LineString import GT3LineString
from shapely.ops import Point, MultiPoint, LineString
from shapely.geometry.polygon import LinearRing
import matplotlib.pyplot as plt
from matplotlib.path import Path
import numpy as np
from GT3.Core.Functions.Cut import cut
import GT3.constants as constants

u_0 = constants.mu_0


class Psi(PlotBase):
    def __init__(self, R, Z, psi, wall, sep_val=1.0):
        """
        The Psi class takes in exp R,Z,psi data and attempts to find important plasma lines and points such as the
        separatrix and x-points.

        :param psi:
        :type psi:
        :param R:
        :type R:
        :param Z:
        :type Z:
        """
        super().__init__()
        self.R = R
        self.Z = Z
        self.psi_exp = psi
        self.wall = wall
        self.sep_val = sep_val
        self._find_xpt_mag_axis()
        self.xpt = [self.xpt_l, self.xpt_u]
        self._calc_psi_norm()

        self.raw_dpsidR = np.abs(np.gradient(self.psi_norm_exp, self.R[0, :], axis=1))
        self.raw_1_R_dpsidR = self.raw_dpsidR / self.R
        self.raw_d_dR_1_R_dpsidR = np.abs(np.gradient(self.raw_1_R_dpsidR, self.R[0, :], axis=1))
        self.raw_dpsidZ = np.abs(np.gradient(self.psi_norm_exp, self.Z[:, 0], axis=0))
        self.raw_d2psidZ2 = np.abs(np.gradient(self.raw_dpsidZ, self.Z[:, 0], axis=0))
        self.raw_dpsidr = self.raw_dpsidR + self.raw_dpsidZ
        self.raw_j = -(self.raw_d_dR_1_R_dpsidR * self.R + self.raw_d2psidZ2) / (self.R * u_0)

        self._calc_pts_lines()
        self._calc_rho2psi_interp()

        # Calculate plasma radius considering potential for wall limiter
        # Inboard limiter
        if self.ibmp_pt[0] < self.wall.bounds[0]:
            self.a = np.average((self.obmp_pt[0] - self.mag_axis[0], self.mag_axis[0] - self.wall.bounds[0]))
        # Outboard limiter
        elif self.obmp_pt[0] >= self.wall.bounds[2]:
            self.a = np.average((self.wall.bounds[2] - self.mag_axis[0], self.mag_axis[0] - self.ibmp_pt[0]))
        else:
            self.a = np.average((self.obmp_pt[0] - self.mag_axis[0], self.mag_axis[0] - self.ibmp_pt[0]))

    def set_rho_mesh(self, rho):
        self.rho = rho
        self.psi = self.rho2psi(self.rho)
        self.psi_norm = self.rho2psinorm(self.rho)
    def get_xpt_lower(self):
        if self.xpt_l is not None:
            return self.xpt_l
        else:
            print("No Upper X-point")
            return None

    def get_xpt_upper(self):
        if self.xpt_u is not None:
            return self.xpt_u
        else:
            print("No Upper X-point")
            return None

    def get_num_xpts(self):
        if self.xpt_u is not None and self.xpt_l is not None:
            return 2
        else:
            return 1

    def get_mag_axis(self):
        return self.mag_axis

    def get_R(self):
        return self.R

    def get_Z(self):
        return self.Z

    def get_psi_exp(self):
        return self.psi_exp

    def get_norm_psi(self):
        return self.psi_norm

    def get_norm_xpt(self):
        return self.norm_xpt

    def get_inboard_midplane(self):
        return self.ibmp_pt

    def get_outboard_midplane(self):
        return self.obmp_pt

    def get_top(self):
        return self.top

    def get_bottom(self):
        return self.bottom

    def get_geom_axis(self):
        return self.geo_axis

    def get_inboard_strike(self):
        return self.ib_strike

    def get_outboard_strike(self):
        return self.ob_strike

    def get_inboard_div_line_short(self, xpt=True):
        if xpt:
            return self.ib_line_cut
        else:
            return self._remove_xpts(self.ib_line_cut)

    def get_inboard_div_line_full(self, xpt=True):
        if xpt:
            return self.ib_line
        else:
            return self._remove_xpts(self.ib_line)

    def get_outboard_div_line_short(self, xpt=True):
        if xpt:
            return self.ob_line_cut
        else:
            return self._remove_xpts(self.ob_line_cut)

    def get_outboard_div_line_full(self, xpt=True):
        if xpt:
            return self.ob_line
        else:
            return self._remove_xpts(self.ob_line)

    def _remove_xpts(self, ls):
        if self.xpt[0] is not None:
            if np.all(np.array(ls.coords[0]) == self.xpt[0]):
                return GT3LineString(ls.coords[1:])
            if np.all(np.array(ls.coords[-1]) == self.xpt[0]):
                return GT3LineString(ls.coords[:-1])
        if self.xpt[1] is not None:
            if np.all(np.array(ls.coords[0]) == self.xpt[1]):
                return GT3LineString(ls.coords[1:])
            if np.all(np.array(ls.coords[-1]) == self.xpt[1]):
                return GT3LineString(ls.coords[:-1])

    def get_lcfs(self):
        return self.lcfs_line

    def get_lcfs_closed(self):
        return self.lcfs_line_closed

    def get_lcfs_ib_2_ob(self):
        return self.entire_sep


    def _find_xpt_mag_axis(self):
        """

        """
        R = self.R
        Z = self.Z
        psi = self.psi_exp
        wall = self.wall
        xpt_u = None
        xpt_l = None
        m_axis = None
        # find x-point location
        dpsidR = np.gradient(psi, R[0, :], axis=1)
        dpsidZ = np.gradient(psi, Z[:, 0], axis=0)
        d2psidR2 = np.gradient(dpsidR, R[0, :], axis=1)
        d2psidZ2 = np.gradient(dpsidZ, Z[:, 0], axis=0)

        # find line(s) where dpsidR=0
        csR = plt.contour(R, Z, dpsidR, [0])
        csZ = plt.contour(R, Z, dpsidZ, [0])

        dpsidR_0 = csR.collections[0].get_paths()
        # dpsidR_0 = cntr.contour(R, Z, dpsidR).trace(0.0)

        # find line(s) where dpsidZ=0
        # dpsidZ_0 = cntr.contour(R, Z, dpsidZ).trace(0.0)
        dpsidZ_0 = csZ.collections[0].get_paths()
        for i, path1 in enumerate(dpsidR_0):
            for j, path2 in enumerate(dpsidZ_0):
                try:
                    # find intersection points between curves for dpsidR=0 and dpsidZ=0
                    ints = GT3LineString(path1.vertices).intersection(GT3LineString(path2.vertices))
                    # if there is only one intersection ('Point'), then we're probably not
                    # dealing with irrelevant noise in psi
                    if isinstance(ints, Point):
                        # check if local maximum or minimum
                        # If the point is within the walls and less than 0.5m from the vessel centroid, this is
                        # very likely the magnetic axis
                        if wall.convex_hull.contains(ints) and wall.centroid.distance(ints) < 0.5:
                            # we've found the magnetic axis
                            m_axis = np.array([ints.x, ints.y])
                        # If the point is not inside the walls, it's a magnet
                        elif not wall.convex_hull.contains(ints):
                            # we've found a magnet. Do nothing.
                            pass
                        elif wall.convex_hull.contains(ints) and ints.y < 0:
                            # The point is within the walls and not near the magnetic axis. This is likely the
                            # lower x-point.
                            # TODO: Provides lower x-point here, but upper x-point can be built out
                            xpt = np.array([ints.x, ints.y])
                            xpt_l = xpt
                        elif wall.convex_hull.contains(ints) and ints.y > 0:
                            # This is an upper x-point. This functionality needs to be built in but is here for
                            # later
                            xpt_u = np.array([ints.x, ints.y])

                        # uncomment this line when debugging
                        # print list(ints.coords), d2psidR2(ints.x, ints.y), d2psidZ2(ints.x, ints.y)

                    # If multiple points are found, the flux surfaces may have come back around onto each other.
                    if isinstance(ints, MultiPoint):
                        for point in ints:
                            # check if local maximum or minimum
                            # If the point is within the walls and less than 0.5m from the vessel centroid, this is
                            # very likely the magnetic axis
                            if wall.convex_hull.contains(point) and wall.centroid.distance(
                                    point) < 0.5:
                                # we've found the magnetic axis
                                m_axis = np.array([point.x, point.y])
                            # If the point is not inside the walls, it's a magnet
                            elif not wall.convex_hull.contains(point):
                                # we've found a magnet. Do nothing.
                                pass
                            elif wall.convex_hull.contains(point) and point.y < 0:
                                # The point is within the walls and not near the magnetic axis. This is likely the
                                # lower x-point.
                                # TODO: Provides lower x-point here, but upper x-point can be built out
                                xpt = np.array([point.x, point.y])
                                xpt_l = xpt
                            elif wall.convex_hull.contains(point) and point.y > 0:
                                # This is an upper x-point. This functionality needs to be built in but is here for
                                # later
                                xpt_u = np.array([point.x, point.y])

                except:
                    pass
        self.xpt_l = xpt_l
        self.xpt_u = xpt_u
        self.mag_axis = m_axis

    def _calc_psi_norm(self, **kwargs):
        # normalize psi
        # psi_interp = Rbf(R, Z, psi)
        # psi_min = psi_interp(axis_mag[0], axis_mag[1])
        #
        # psi_shifted = psi - psi_min  # set center to zero
        # psi_shifted_interp = Rbf(R, Z, psi_shifted)
        # psi_shifted_xpt = psi_shifted_interp(xpt[0], xpt[1])
        #
        # psi_norm = psi_shifted / psi_shifted_xpt
        R = self.R
        Z = self.Z
        psi = self.psi_exp
        axis_mag = self.mag_axis
        psi_min = griddata(np.column_stack((R.flatten(), Z.flatten())),
                           psi.flatten(),
                           [axis_mag[0], axis_mag[1]],
                           method='cubic')

        psi_shifted = psi - psi_min  # set center to zero
        psi_shifted_xpt_l, psi_shifted_xpt_u = None, None

        if self.xpt[0] is not None:

            psi_shifted_xpt_l = griddata(np.column_stack((R.flatten(), Z.flatten())),
                                   psi_shifted.flatten(),
                                   [self.xpt[0][0], self.xpt[0][1]],
                                   method='cubic')
        if self.xpt[1] is not None:
            psi_shifted_xpt_u = griddata(np.column_stack((R.flatten(), Z.flatten())),
                                   psi_shifted.flatten(),
                                   [self.xpt[1][0], self.xpt[1][1]],
                                   method='cubic')
        psi_shifted_xpt = [psi_shifted_xpt_l, psi_shifted_xpt_u]
        if self.xpt[1] is None:
            psi_norm = psi_shifted / np.average(psi_shifted_xpt_l)
            num_xpts = 1
            norm_xpt = self.xpt[0]
            if kwargs.get("debug"):
                print("DEBUG: 1 XPT - psi normalized to lower X-point  " + str(self.xpt[0]))
        elif self.xpt[0] is None:
            num_xpts = 1
            norm_xpt = self.xpt[1]
            if kwargs.get("debug"):
                print("DEBUG: 1 XPT - psi normalized to upper X-point:  " + str(self.xpt[1]))
            psi_norm = psi_shifted / np.average(psi_shifted_xpt_u)
        else:
            num_xpts = 2
            norm_xpt = self.xpt[1]
            if kwargs.get("debug"):
                print("DEBUG: 2 XPT -  psi normalized to upper X-point: " + str(self.xpt[1]))
            # psi_norm = psi_shifted / np.average(psi_shifted_xpt)
            psi_norm = psi_shifted / np.average(psi_shifted_xpt_u)

        self.psi_norm_exp = psi_norm
        self.norm_xpt = norm_xpt
        self.num_xpts = num_xpts

    def _calc_pts_lines(self):

        """
        Calculates points and lines for the flux surfaces

        :param debug: Activate debugging function and plots
        :type debug: bool
        :param xpt:
        :param wall:
        :type wall: LineString
        :param mag_axis:
        :return:
        :rtype: object
        """

        xpt = self.xpt
        wall = self.wall
        mag_axis = self.mag_axis
        sep_val = self.sep_val
        norm_xpt = self.norm_xpt

        if sep_val < 1.0:
            # We have manually overridden the separatrix to exclude the X-points and SOL. We're only able to get a few
            # parameters. We will assume only 1 contour will have this value
            contours = plt.contour(self.R[0], self.Z[:, 0], self.psi_norm_exp, [sep_val]).collections[0].get_paths()
            for contour in contours:
                ls = GT3LineString(contour.vertices)
                if wall.convex_hull.contains(ls):
                    ibmp = [ls.bounds[0], mag_axis[1]]
                    obmp = [ls.bounds[2], mag_axis[1]]
                    geo = [ls.centroid.xy[0][0], ls.centroid.xy[1][0]]
                    top = [ls.centroid.xy[0][0], ls.bounds[3]]
                    bottom = [ls.centroid.xy[0][0], ls.bounds[1]]
                    mag = mag_axis

                    self.ibmp_pt = ibmp
                    self.obmp_pt = obmp
                    self.top = top
                    self.bottom = bottom
                    self.geo_axis = geo
                    self.ib_strike = None
                    self.ob_strike = None
                    self.ib_line_cut = None
                    self.ob_line_cut = None
                    self.ib_line = None
                    self.ob_line = None

                    self.lcfs_line = ls
                    self.lcfs_line_closed = ls
                    self.entire_sep = ls
            raise Exception("Could not find an LCFS with overwritten sep_val")
        # We find the potential sepatartrix lines first
        contours_paths = plt.contour(self.R[0], self.Z[:, 0], self.psi_norm_exp, [sep_val]).collections[0].get_paths()

        if xpt[0] is not None and xpt[1] is not None:
            # We have 2 X-points. We'll basically ignore one, and SOL will be disabled in the rest of GT3. Any neutrals
            # calculations will have to be omitted. We also know that this line (possibly by definition) will contain
            # an x-point and have IB/OB strike points. We'll overwrite contours_paths to say that this is the only
            # contour worth caring about in the rest of the calculation.

            temp_contours = plt.contour(self.R[0], self.Z[:, 0], self.psi_norm_exp, [sep_val]).collections[
                0].get_paths()
            if len(temp_contours) == 1:
                ls = GT3LineString(temp_contours.vertices)
                if ls.convex_hull.contains(Point(xpt[0])):
                    # The upper X-point will define our separatrix
                    xpt_sep = xpt[0]
                elif ls.convex_hull.contains(Point(xpt[1])):
                    # The lower X-point will define our separatrix
                    xpt_sep = xpt[1]
                else:
                    raise Exception("There are 2 x-points, but neither are on the separatrix as found in CalcPtsLines.")
            else:
                for path in temp_contours:
                    ls = GT3LineString(path.vertices, wall=wall)
                    if ls.convex_hull.contains(Point(xpt[0])):
                        # The lower X-point will define our separatrix
                        xpt_sep = xpt[0]
                        xpt_temp = xpt_sep
                        if ls.contains(xpt_sep):
                            contours_paths = path
                        else:
                            print("FUck")
                        break
                    elif ls.convex_hull.contains(Point(xpt[1])):
                        # The upper X-point will define our separatrix
                        xpt_sep = xpt[1]
                        xpt_temp = xpt_sep
                        contours_paths = path
                        break
        else:
            # Grab the actual x-point
            if xpt[0] is not None:
                xpt_temp = xpt[0]
            else:
                xpt_temp = xpt[1]

        # This line might render code above unnecessary?

        xpt_temp = norm_xpt

        # create lines for seperatrix and divertor legs of seperatrix

        # Convert back to arrays since I don't want to deal with changing everything to shapely since this already works
        if isinstance(contours_paths, Path):
            contours = [contours_paths.vertices]
        else:
            contours = [a.vertices for a in contours_paths]
        # Replace points in the vic. of the xpt with the xpt
        contour_new = []
        for contour in contours:
            dist = np.sqrt((contour[:, 0] - xpt_temp[0]) ** 2 + (contour[:, 1] - xpt_temp[1]) ** 2)
            # TODO: Put this distance in the input file or make distance selection smarter
            contour[dist < 0.03] = xpt_temp
            contour_new.append(contour[np.where((contour != np.roll(contour, -1, axis=0)).all(axis=1))])

        # if there is only one contour, then it's drawing one continuous line from ib strike point to ob strike point
        # if there are two contours, then it's drawing separate seperatrix and divertor legs

        if len(contour_new) == 1:
            contour = contour_new[0]

            # make sure the line goes from inboard to outboard

            if contour[-1, 0] < contour[0, 0]:
                contour = np.flipud(contour)

            xpt_loc = np.where((contour == xpt_temp).all(axis=1))[0]
            try:
                ib = np.flipud(contour[:xpt_loc[0] + 1])
            except IndexError:
                # The separatrix did not converge to the x-point. We should cut and project 2 times to create the LCFS.
                raise Exception("The separatrix does not include the x-point.")
                ls = GT3LineString(contour)
                pt = Point(xpt_temp)
                nearest = nearest_points(ls, pt)[0]
                cuts = split(ls, nearest)
                temp_ls = GT3LineString(cuts[0].coords[:] + pt.coords[:] + cuts[1].coords[:])
                _plot_contours(temp_ls, title="Temporary LS from cuts")

            if len(xpt_loc) == 1:
                second_point = _find_ob(contour, xpt_sep)
                xpt_loc = np.append(xpt_loc, second_point)
                # raise Exception("Only 1 xpoint index found. Should be 2. Fix me")

            lcfs = contour[xpt_loc[0]:xpt_loc[1]]
            ob = contour[xpt_loc[1]:]

        elif len(contour_new) > 1:
            for contour in contour_new:
                # determine if contour is seperatrix, divertor legs, a main ib2ob line, or something else
                contour = np.array(contour)
                if not GT3LineString(contour).intersects(wall) and wall.convex_hull.contains(GT3LineString(contour)):
                    # Found the Seperatrix

                    lcfs = contour
                    continue
                # Meaningless noise lines will not pass near the x-point, so check that first to elliminate them
                if LineString(contour).distance(Point(xpt_temp)) < 0.01:
                    # if the largest y value is larger than the y value of the magnetic axis AND the line
                    # intersects the wall twice, then it's an ib2ob line and there are more than one psi_norm=1
                    # contours because the others are noise. Treat this one the same way we would if there were
                    # only one contour.

                    # count number of intersections with the wall
                    wall_ints = len(GT3LineString(contour).intersection(wall))
                    contour_ls = GT3LineString(contour)

                    if wall_ints >= 2 \
                            and np.amax(contour[:, 1]) > mag_axis[1]\
                            and not (contour >= mag_axis).all() \
                            and contour_ls.convex_hull.contains(Point(mag_axis)):

                        # This could be a wall-limited LCFS.
                        # TODO: Check and eliminate

                        # then we probably have an ib2ob line. Treat the same way as above
                        # We also are testing to make sure we don't have an upper null

                        # make sure the line goes from inboard to outboard
                        if contour[-1, 0] < contour[0, 0]:
                            contour = np.flipud(contour)

                        xpt_loc = np.where((contour == xpt_temp).all(axis=1))[0]

                        ib = np.flipud(contour[:xpt_loc[0] + 1])
                        lcfs = contour[xpt_loc[0]:xpt_loc[1]]
                        ob = contour[xpt_loc[1]:]

                    # if the largest y value is larger than the y value of the magnetic axis
                    elif np.amax(contour[:, 1]) > mag_axis[1]:
                        # we have the main seperatrix
                        lcfs = contour

                    else:  # we have the divertor legs
                        # make sure the line goes from inboard to outboard
                        if contour[-1, 0] < contour[0, 0]:
                            contour = np.flipud(contour)

                        # create ib and ob divertor legs, each starting at the x-point
                        xpt_loc = np.where((contour == xpt_temp).all(axis=1))[0][0]
                        ob = contour[xpt_loc:]
                        ib = np.flipud(contour[:xpt_loc + 1])

        # create lcfs lines
        lcfs_line = GT3LineString(lcfs)
        lcfs_line_closed = LinearRing(lcfs)

        # create ib and ob linestrings and truncate at the wall
        ib_line = GT3LineString(ib)
        ob_line = GT3LineString(ob)

        ib_strike = ib_line.intersection(wall)
        ob_strike = ob_line.intersection(wall)

        ib_line_cut = cut(ib_line, ib_line.project(ib_strike, normalized=True))[0]
        ob_line_cut = cut(ob_line, ob_line.project(ob_strike, normalized=True))[0]

        entire_sep = np.vstack((np.flipud(ib), lcfs[1:], ob))
        entire_sep_line = GT3LineString(entire_sep)

        # create points object
        obmp_pt = lcfs[np.argmax(lcfs, axis=0)[0]]
        ibmp_pt = lcfs[np.argmin(lcfs, axis=0)[0]]
        top_pt = lcfs[np.argmax(lcfs, axis=0)[1]]
        bot_pt = lcfs[np.argmin(lcfs, axis=0)[1]]
        axis_geo = [(obmp_pt[0] + ibmp_pt[0]) / 2, (obmp_pt[1] + ibmp_pt[1]) / 2]
        self.ibmp_pt = ibmp_pt
        self.obmp_pt = obmp_pt
        self.top = top_pt
        self.bottom = bot_pt
        self.geo_axis = axis_geo
        self.ib_strike = ib_strike
        self.ob_strike = ob_strike
        self.ib_line_cut = ib_line_cut
        self.ob_line_cut = ob_line_cut
        self.ib_line = ib_line
        self.ob_line = ob_line

        self.lcfs_line = lcfs_line
        self.lcfs_line_closed = lcfs_line_closed
        self.entire_sep = entire_sep_line

        # Create separatrix that is cut at the wall

        ib_line_cut_nox = self.get_inboard_div_line_short(xpt=False)
        ob_line_cut_nox = self.get_outboard_div_line_short(xpt=False)

        entire_sep_cut = np.vstack((
            np.array([ib_line_cut_nox.coords.xy[0], ib_line_cut_nox.coords.xy[1]]).T,
            lcfs[1:],
            np.array([ob_line_cut_nox.coords.xy[0], ob_line_cut_nox.coords.xy[1]]).T
        ))
        entire_sep_cut_line = GT3LineString(entire_sep_cut)

        self.entire_sep_cut = entire_sep_cut_line

    def _calc_rho2psi_interp(self):
        rho_vals = np.linspace(0, self.sep_val, 100)

        ptsRZ = np.zeros((len(rho_vals), 2))

        obmp_line = GT3LineString([Point(self.mag_axis), Point(self.obmp_pt)])
        for i, rho in enumerate(rho_vals):
            ptsRZ[i] = np.asarray(obmp_line.interpolate(rho, normalized=True).coords)

        psi_vals = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())),
                            self.psi_exp.flatten(),
                            ptsRZ,
                            method='cubic')
        if psi_vals[0] > psi_vals[-1]:
            print("Psi values are decreasing in Psi calculation. Flipping")
            psi_vals = np.flip(psi_vals)

        psinorm_vals = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())),
                                self.psi_norm_exp.flatten(),
                                ptsRZ,
                                method='cubic')

        psinorm_vals[0] = 0
        psinorm_vals[-1] = self.sep_val

        # For some reason, psi_vals will sometimes not be monotonically increasing, especially near the magnetic axis.
        # This prevents UnivariateSpline from working. To prevent this, we're just going to delete any points that are
        # lower than the previous point, prior to doing the fit. As long as the psi data isn't too wonky, this should be
        # fine.
        psivals_mi = []  # psivals that are monotonically increasing
        psinormvals_mi = []  # psinormvals corresponding to monotonically increasing psi
        rhovals_mi = []  # rhovals corresponding to monotonically increasing psi
        for i, psi_val in enumerate(psi_vals):
            if i == 0:
                rhovals_mi.append(0)
                psivals_mi.append(psi_vals[0])  # this is probably supposed to be zero as well, but we'll leave it for now
                psinormvals_mi.append(0)
            elif psi_val > psivals_mi[-1]:
                rhovals_mi.append(rho_vals[i])
                psivals_mi.append(psi_vals[i])
                psinormvals_mi.append(psinorm_vals[i])

        rhovals_mi = np.asarray(rhovals_mi)
        psivals_mi = np.asarray(psivals_mi)
        psinormvals_mi = np.asarray(psinormvals_mi)

        self.rho2psi = interp1d(rhovals_mi, psivals_mi, fill_value='extrapolate')
        self.rho2psinorm = interp1d(rhovals_mi, psinormvals_mi, fill_value='extrapolate')
        self.psi2rho = interp1d(psivals_mi, rhovals_mi, fill_value='extrapolate')
        self.psi2psinorm = interp1d(psivals_mi, psinormvals_mi, fill_value='extrapolate')
        self.psinorm2rho = interp1d(psinormvals_mi, rhovals_mi, fill_value='extrapolate')
        self.psinorm2psi = interp1d(psinormvals_mi, psivals_mi, fill_value='extrapolate')


    def _find_ob(self, contour, xpt):
        xpt_Pt = Point(xpt)
        dist = 10
        res_pt = None
        for coord in contour:
            if coord[0] > xpt[0]:
                if Point(coord).distance(xpt_Pt) < dist:
                    res_pt = coord
                    dist = Point(coord).distance(xpt_Pt)

        if not res_pt.all():
            raise Exception("Failed to find the second splice point for separatrix splicing.")
        else:
            return np.where((contour==res_pt).all(axis=1))[0][0]


    def _plot_exp_psi_raw(self, res=50):
        try:
            ax = self._plot_with_wall()
        except:
            ax = self._plot_without_wall()
        try:
            cs = ax.contour(self.R, self.Z, self.psi_exp, res)
            plt.clabel(cs, inline=1, fontsize=10)
            return cs
        except NameError:
            print("Psi not defined")
            pass

    def plot_exp_psi(self, res=50):
        try:
            ax = self._plot_with_wall()
        except:
            return
        try:
            ax.contour(self.R, self.Z, self.psi_norm_exp, res)
            return ax
        except NameError:
            print("Psi not defined")
            pass