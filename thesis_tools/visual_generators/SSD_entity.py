import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import shapely
from sympy.solvers import solve
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.patches import Arrow

n_points_v_ring = 360 # Number of points in speed rings
max_v = 100. # Scale diagrams to max v [m/s]
R_PZ = 50. # Radius of protected zone [m]
v_range = 25. # Range on the axes of a plot
vec_x = np.array([[1], [0]])
vec_y = np.array([[0], [1]])
unit_dartipp = Polygon([(0., -0.5), (1., -1.), (0., 1.), (-1., -1.)]) # unit darttip
darttip_scale = 1.5 # default scaling of darttip

unit_target_triangle = Polygon([(0., 0.), (1., 1.), (-1., 1.)])
target_triangle_scale = 1.0 # default scaling of target symbol

# Scenario specific constants
half_geofence_width = 200. # [m]

class SSD(object):
    def __init__(self, wind_spd = 0., wind_dir = 0.):
        self.wind_spd = wind_spd # [m/s]
        self.wind_dir = np.deg2rad(wind_dir) # [rad]

        self.v_wind_x = -self.wind_spd * np.sin(self.wind_dir) # wind blowing to east [m/s]
        self.v_wind_y = -self.wind_spd * np.cos(self.wind_dir) # wind blowing to north [m/s]

        self.geofence = ConvexGeofence()

        self.SSD_entities = {}

    def add_SSD_entity(self, callsign, pos_east, pos_north, trk, airspeed, v_max, v_min = None):
        self.SSD_entities[callsign] = SSD_entity(self, callsign, pos_east, pos_north, trk, airspeed, v_max, v_min)

class SSD_entity(object):
    def __init__(self, SSD, callsign, pos_east, pos_north, trk, airspeed, v_max, v_min = None):
        self.SSD = SSD
        self.callsign = callsign # string
        self.x = pos_east # [m]
        self.y = pos_north # [m]
        self.trk = np.deg2rad(trk) # [rad]
        self.hdg = np.deg2rad(trk_to_hdg(trk, airspeed, SSD.v_wind_y, SSD.v_wind_x))
        self.airspeed = airspeed # [m/s]
        self.v_max = v_max # [m/s]
        self.v_min = v_min # [None] (if RC) / [m/s] if FW

        self.gs, self.hdg = airspeed_to_groundsspeed(self.airspeed, self.trk, self.SSD.v_wind_x, self.SSD.v_wind_y)
        self.gs_x = self.gs * np.sin(self.trk)
        self.gs_y = self.gs * np.cos(self.trk)

    def construct_RV_polygon(self):
        # construct inner and outer speed ring coordinates
        angles = np.linspace(0, 2. * np.pi, n_points_v_ring)
        x_unit_circle = np.sin(angles) 
        y_unit_circle = np.cos(angles)

        # outer circle
        x_outer_v = x_unit_circle * self.v_max + SSD.v_wind_x
        y_outer_v = y_unit_circle * self.v_max + SSD.v_wind_y

        outer_v_polygon = Polygon(list(map(tuple, np.concatenate((x_outer_v, y_outer_v)).reshape((2, n_points_v_ring)).T)))

        # inner circle and construct RV
        if self.v_min is not None:
            x_inner_v = x_unit_circle * self.v_min + SSD.v_wind_x
            y_inner_v = y_unit_circle * self.v_min + SSD.v_wind_y

            inner_v_polygon = Polygon(list(map(tuple, np.concatenate((x_inner_v, y_inner_v)).reshape((2, n_points_v_ring)).T)))
            RV_polygon = outer_v_polygon.difference(inner_v_polygon)
        else:
            RV_polygon = outer_v_polygon

        return RV_polygon

    def construct_RV_polygon_HDG(self):
        # construct inner and outer speed ring coordinates
        angles = np.linspace(0, 2. * np.pi, n_points_v_ring)
        x_unit_circle = np.sin(angles) 
        y_unit_circle = np.cos(angles)

        # outer circle
        x_outer_v = x_unit_circle * (self.airspeed + 0.001) + SSD.v_wind_x
        y_outer_v = y_unit_circle * (self.airspeed + 0.001) + SSD.v_wind_y

        outer_v_polygon = Polygon(list(map(tuple, np.concatenate((x_outer_v, y_outer_v)).reshape((2, n_points_v_ring)).T)))

        # inner circle and construct RV
        if self.v_min is not None:
            x_inner_v = x_unit_circle * (self.airspeed - 0.001) + SSD.v_wind_x
            y_inner_v = y_unit_circle * (self.airspeed - 0.001) + SSD.v_wind_y

            inner_v_polygon = Polygon(list(map(tuple, np.concatenate((x_inner_v, y_inner_v)).reshape((2, n_points_v_ring)).T)))
            RV_polygon = outer_v_polygon.difference(inner_v_polygon)
        else:
            RV_polygon = outer_v_polygon

        return RV_polygon

    def construct_VO_polygon(self, intruder):
        # First construct CC and move points afterwards
        dx = self.x - intruder.x
        dy = self.y - intruder.y
        d_rel = np.sqrt(dx**2 + dy**2)
        alpha = np.arcsin(R_PZ / d_rel)
        d_psi = np.arctan2(-dx, -dy)

        # generate unit vector of 2 vertices and rotate by d_psi
        R = np.matrix([[ np.cos(d_psi), np.sin(d_psi)],\
                         [-np.sin(d_psi), np.cos(d_psi)]])

        n_1 = R * np.matrix([[np.sin(alpha)],[np.cos(alpha)]])
        n_2 = R * np.matrix([[-np.sin(alpha)], [np.cos(alpha)]])

        # compute polygon points of VO
        p_0 = (intruder.gs_x, intruder.gs_y) # (x, y)
        p_1 = (np.array(n_1)[0][0] * 2. * max_v + intruder.gs_x, np.array(n_1)[1][0] * 2. * max_v + intruder.gs_y) 
        p_2 = (np.array(n_2)[0][0] * 2. * max_v + intruder.gs_x, np.array(n_2)[1][0] * 2. * max_v + intruder.gs_y) 

        VO_polygon = Polygon([p_0, p_1, p_2])

        return VO_polygon

    # Constructs the forbidden velocity space of VOs
    def construct_VOs_polygon(self):
        VO_polygons = []
        for callsign, intruder in self.SSD.SSD_entities.items():
            if callsign == self.callsign:
                continue
            else:
                VO_polygons.append(self.construct_VO_polygon(intruder))

        VOs_polygon = VO_polygons[0]
        if len(VO_polygons) > 1:
            for i in np.arange(1, len(VO_polygons)):
                VOs_polygon = VOs_polygon.union(VO_polygons[i])
        return VOs_polygon
    
    # Construct a 'virtual' VO of a geofence segment
    def construct_VO_geofence_segment_polygon(self, intruder_id, segment_idx):
        # First load intruder and geofence segment data
        intruder = self.SSD.SSD_entities[intruder_id]

        d_int = np.array([[intruder.x - self.x], \
                          [intruder.y - self.y]])
        v_int = np.array([[intruder.gs_x], \
                          [intruder.gs_y]])

        segment = self.SSD.geofence.segments[segment_idx]

        # Check if intruder is converging towards segment, else return None
        v_int_dot_y_p = v_int[0][0] * segment.y_p[0][0] + v_int[1][0] * segment.y_p[1][0]
        if (v_int_dot_y_p >= 0):
            return None

        d_int_dot_x_p = d_int[0][0] * segment.x_p[0][0] + d_int[1][0] * segment.x_p[1][0]
        d_int_dot_y_p = d_int[0][0] * segment.y_p[0][0] + d_int[1][0] * segment.y_p[1][0]

        a_geo = np.array([[segment.a[0] - self.x], \
                          [segment.a[1] - self.y]])
        d_geo = -1. * (a_geo[0][0] * segment.y_p[0][0] + a_geo[1][0] * segment.y_p[1][0]) 

        # Secondly define double rotated reference system
        phi_p = 0.5 * np.arctan2(-d_int_dot_x_p, d_int_dot_y_p)
        phi_total = segment.phi + phi_p
        
        cos_phi_p = np.cos(phi_p)
        sin_phi_p = np.sin(phi_p)
        cos_phi_total = np.cos(phi_total)
        sin_phi_total = np.sin(phi_total)

        R = np.array([[cos_phi_total, -sin_phi_total], \
                      [sin_phi_total, cos_phi_total]])
        
        x_pp = np.linalg.solve(R.T, vec_x)
        y_pp = np.linalg.solve(R.T, vec_y)

        # dot products
        d_int_dot_x_pp = d_int[0][0] * x_pp[0][0] + d_int[1][0] * x_pp[1][0]
        d_int_dot_y_pp = d_int[0][0] * y_pp[0][0] + d_int[1][0] * y_pp[1][0]
        v_int_dot_x_pp = v_int[0][0] * x_pp[0][0] + v_int[1][0] * x_pp[1][0]
        v_int_dot_y_pp = v_int[0][0] * y_pp[0][0] + v_int[1][0] * y_pp[1][0]
        d_int_dot_v_int = d_int[0][0] * v_int[0][0] + d_int[1][0] * v_int[1][0]

        # Constants needed to calculate 'virtual' VO geometries
        C1 = 1. + sin_phi_p * d_int_dot_x_pp / d_geo
        C2 = 1. + cos_phi_p * d_int_dot_y_pp / d_geo
        C3 = -2. * v_int_dot_x_pp - sin_phi_p * d_int_dot_v_int / d_geo
        C4 = -2. * v_int_dot_y_pp - cos_phi_p * d_int_dot_v_int / d_geo
        Cx_pp = - C3 / (2. * C1)
        Cy_pp = - C4 / (2. * C2)
        a2 = (-intruder.gs**2 + C2 * Cy_pp**2) / C1 + Cx_pp**2 
        b2 = (-intruder.gs**2 + C1 * Cx_pp**2) / C2 + Cy_pp**2 

        # now construct polygon based on sign of b2
        VO = None
        # If ownshup outside geofence
        if a2 <= 0:
            pass
        elif b2 > 0:
            ellipse_angles = np.linspace(0., 2. * np.pi, n_points_v_ring)
            rotated_xs = np.sqrt(a2) * np.cos(ellipse_angles) + Cx_pp
            rotated_ys = np.sqrt(b2) * np.sin(ellipse_angles) + Cy_pp
        else:
            tmax = np.log((20. * self.v_max + np.sqrt(20.**2 * self.v_max**2 + a2)) / np.sqrt(a2))
            tmin = -tmax
            if (phi_p > 0):
                t = np.linspace(tmin, tmax, n_points_v_ring)
                rotated_xs = -np.sqrt(a2) * np.cosh(t) + Cx_pp
            else:
                t = np.linspace(tmax, tmin, n_points_v_ring)
                rotated_xs = np.sqrt(a2) * np.cosh(t) + Cx_pp
            rotated_ys = np.sqrt(-b2) * np.sinh(t) * Cy_pp

        if a2 > 0:
            non_rotated_xs = np.array([rotated_xs * cos_phi_total - rotated_ys * sin_phi_total])
            non_rotated_ys = np.array([rotated_xs * sin_phi_total + rotated_ys * cos_phi_total])
            xy_tuple_array = np.concatenate((non_rotated_xs, non_rotated_ys), axis=0).transpose()
            VO = Polygon(xy_tuple_array)
        return VO

    def construct_VOs_geofence_polygon(self):
        VO_polygons = []
        for callsign, intruder in self.SSD.SSD_entities.items():
            if callsign == self.callsign:
                continue
            else:
                for i in range(len(self.SSD.geofence.segments)):
                    polygon = self.construct_VO_geofence_segment_polygon(callsign, i)
                    if polygon is not None:
                        VO_polygons.append(polygon)

        # Return None if no geofence VOs
        if len(VO_polygons) == 0:
            return None

        # |Compute union of Vos    
        VOs_polygon = VO_polygons[0]
        if len(VO_polygons) > 1:
            for i in np.arange(1, len(VO_polygons)):
                VOs_polygon = VOs_polygon.union(VO_polygons[i])

        return VOs_polygon
    
    ############
    # PLOTTERS #
    ############

    def plot_RV(self):
        RV_polygon = self.construct_RV_polygon()

        fig, ax = plt.subplots()
        path = pathify(RV_polygon)
        patch = PathPatch(path, facecolor='lightgrey', edgecolor='black')
        ax.add_patch(patch)

        # draw own speed vector
        speed_vector = Arrow(self.SSD.v_wind_x, self.SSD.v_wind_y, self.gs_x, self.gs_y, facecolor='black', edgecolor='black')
        ax.add_patch(speed_vector)

        # draw dartip
        darttip = shapely.affinity.scale(unit_dartipp, xfact=darttip_scale, yfact=darttip_scale, origin=(0,0))
        darttip = shapely.affinity.rotate(darttip, -np.rad2deg(self.trk), origin=(0,0))
        darttip_path = pathify(darttip)
        darttip_patch = PathPatch(darttip_path, facecolor='black', edgecolor='black')
        ax.add_patch(darttip_patch)

        plt.xlim(-v_range, v_range)
        plt.ylim(-v_range, v_range)
        plt.axis('off')
        ax.set_aspect('equal')
        plt.show()

    def plot_VOs(self):
        VOs_polygon = self.construct_VOs_polygon()

        fig, ax = plt.subplots()

        if VOs_polygon.geom_type == 'Polygon':
            path = pathify(VOs_polygon)
            patch = PathPatch(path, facecolor='red', edgecolor='black')
            ax.add_patch(patch)
        elif VOs_polygon.geom_type == 'MultiPolygon':
            for geom in VOs_polygon.geoms:
                path = pathify(geom)
                patch = PathPatch(path, facecolor='red', edgecolor='black')
                ax.add_patch(patch)

        # draw own speed vector
        speed_vector = Arrow(self.SSD.v_wind_x, self.SSD.v_wind_y, self.gs_x, self.gs_y, facecolor='black', edgecolor='black')
        ax.add_patch(speed_vector)

        # draw dartip
        darttip = shapely.affinity.scale(unit_dartipp, xfact=darttip_scale, yfact=darttip_scale, origin=(0,0))
        darttip = shapely.affinity.rotate(darttip, -np.rad2deg(self.trk), origin=(0,0))
        darttip_path = pathify(darttip)
        darttip_patch = PathPatch(darttip_path, facecolor='black', edgecolor='black')
        ax.add_patch(darttip_patch)

        plt.xlim(-v_range, v_range)
        plt.ylim(-v_range, v_range)
        plt.axis('off')
        ax.set_aspect('equal')
        plt.show()
    
    # plots the 'virtual' VOs for geofence segments
    def plot_VOs_geofences(self):
        VOs_polygon = self.construct_VOs_geofence_polygon()

        fig, ax = plt.subplots()

        if VOs_polygon.geom_type == 'Polygon':
            path = pathify(VOs_polygon)
            patch = PathPatch(path, facecolor='red', edgecolor='black')
            ax.add_patch(patch)
        elif VOs_polygon.geom_type == 'MultiPolygon':
            for geom in VOs_polygon.geoms:
                path = pathify(geom)
                patch = PathPatch(path, facecolor='red', edgecolor='black')
                ax.add_patch(patch)

        # draw own speed vector
        speed_vector = Arrow(self.SSD.v_wind_x, self.SSD.v_wind_y, self.gs_x, self.gs_y, facecolor='black', edgecolor='black')
        ax.add_patch(speed_vector)

        # draw dartip
        darttip = shapely.affinity.scale(unit_dartipp, xfact=darttip_scale, yfact=darttip_scale, origin=(0,0))
        darttip = shapely.affinity.rotate(darttip, -np.rad2deg(self.trk), origin=(0,0))
        darttip_path = pathify(darttip)
        darttip_patch = PathPatch(darttip_path, facecolor='black', edgecolor='black')
        ax.add_patch(darttip_patch)

        plt.xlim(-v_range, v_range)
        plt.ylim(-v_range, v_range)
        plt.axis('off')
        ax.set_aspect('equal')
        plt.show()

    def plot_SSD(self):
        fig, ax = plt.subplots()

        RV_polygon = self.construct_RV_polygon()
        VOs_polygon = self.construct_VOs_polygon()
        VOs_gf_polygon = self.construct_VOs_geofence_polygon()

        if VOs_gf_polygon is not None:
            VOs_polygon = VOs_polygon.union(VOs_gf_polygon)

        ARV = RV_polygon.difference(VOs_polygon)
        FRV = RV_polygon.intersection(VOs_polygon)

        if FRV.geom_type == 'Polygon':
            path = pathify(FRV)
            patch = PathPatch(path, facecolor='red', edgecolor='black')
            ax.add_patch(patch)
        elif FRV.geom_type == 'MultiPolygon':
            for geom in FRV.geoms:
                path = pathify(geom)
                patch = PathPatch(path, facecolor='red', edgecolor='black')
                ax.add_patch(patch)

        if ARV.geom_type == 'Polygon':
            path = pathify(ARV)
            patch = PathPatch(path, facecolor='lightgrey', edgecolor='black')
            ax.add_patch(patch)
        elif ARV.geom_type == 'MultiPolygon':
            for geom in ARV.geoms:
                path = pathify(geom)
                patch = PathPatch(path, facecolor='lightgrey', edgecolor='black')
                ax.add_patch(patch)

        # draw own speed vector
        speed_vector = Arrow(self.SSD.v_wind_x, self.SSD.v_wind_y, self.gs_x, self.gs_y, facecolor='black', edgecolor='black')
        ax.add_patch(speed_vector)

        # draw dartip
        darttip = shapely.affinity.scale(unit_dartipp, xfact=darttip_scale, yfact=darttip_scale, origin=(0,0))
        darttip = shapely.affinity.rotate(darttip, -np.rad2deg(self.trk), origin=(0,0))
        darttip_path = pathify(darttip)
        darttip_patch = PathPatch(darttip_path, facecolor='black', edgecolor='black')
        ax.add_patch(darttip_patch)

        plt.xlim(-v_range, v_range)
        plt.ylim(-v_range, v_range)
        plt.axis('off')
        ax.set_aspect('equal')
        plt.show()

    def plot_SSD_solution_points(self, target_trk = None):
        # Define track in radians for DEST resolution
        if target_trk == None:
            target_trk = self.trk
        else:
            target_trk = np.deg2rad(target_trk)

        fig, ax = plt.subplots()

        RV_polygon = self.construct_RV_polygon()
        RV_HDG_polygon = self.construct_RV_polygon_HDG()
        VOs_polygon = self.construct_VOs_polygon()
        VOs_gf_polygon = self.construct_VOs_geofence_polygon()

        if VOs_gf_polygon is not None:
            VOs_polygon = VOs_polygon.union(VOs_gf_polygon)

        ARV = RV_polygon.difference(VOs_polygon)
        FRV = RV_polygon.intersection(VOs_polygon)

        if FRV.geom_type == 'Polygon':
            path = pathify(FRV)
            patch = PathPatch(path, facecolor='red', edgecolor='black')
            ax.add_patch(patch)
        elif FRV.geom_type == 'MultiPolygon':
            for geom in FRV.geoms:
                path = pathify(geom)
                patch = PathPatch(path, facecolor='red', edgecolor='black')
                ax.add_patch(patch)

        if ARV.geom_type == 'Polygon':
            path = pathify(ARV)
            patch = PathPatch(path, facecolor='lightgrey', edgecolor='black')
            ax.add_patch(patch)
        elif ARV.geom_type == 'MultiPolygon':
            for geom in ARV.geoms:
                path = pathify(geom)
                patch = PathPatch(path, facecolor='lightgrey', edgecolor='black')
                ax.add_patch(patch)

        # Calculate resolutions
        # First derive ARV vertices for OPT and DEST
        # Loop through all exteriors and append. Afterwards concatenate
        p = []
        q = []
        for geom in ARV.geoms:
            p.append(np.array(geom.exterior.coords))
            q.append(np.diff(np.row_stack((p[-1], p[-1][0])), axis=0))
            try:
                p.append(np.array(geom.interior.coords))
                q.append(np.diff(np.row_stack((p[-1], p[-1][0])), axis=0))
            except:
                pass
        p = np.concatenate(p)
        q = np.concatenate(q)
        # Calculate squared distance between edges
        l2 = np.sum(q ** 2, axis=1)
        # Catch l2 == 0 (exception)
        same = l2 < 1e-8
        l2[same] = 1.

        # Derive ARV vertices for HDG resolution strategy
        # Loop through all exteriors and append. Afterwards concatenate
        p_HDG = []
        q_HDG = []
        # Create ARV for HDG resolution strategy
        ARV_HDG = ARV.intersection(RV_HDG_polygon)
        for geom in ARV_HDG.geoms:
            p_HDG.append(np.array(geom.exterior.coords))
            q_HDG.append(np.diff(np.row_stack((p_HDG[-1], p_HDG[-1][0])), axis=0))
            try:
                p_HDG.append(np.array(geom.interior.coords))
                q_HDG.append(np.diff(np.row_stack((p_HDG[-1], p_HDG[-1][0])), axis=0))
            except:
                pass
        p_HDG = np.concatenate(p_HDG)
        q_HDG = np.concatenate(q_HDG)
        # Calculate squared distance between edges
        l2_HDG = np.sum(q_HDG ** 2, axis=1)
        # Catch l2 == 0 (exception)
        same_HDG = l2_HDG < 1e-8
        l2_HDG[same_HDG] = 1.

        # calc gs to target
        target_hdg = np.deg2rad(trk_to_hdg(np.rad2deg(target_trk), self.airspeed, self.SSD.v_wind_y, self.SSD.v_wind_x))
        gs_x_DEST = np.sin(target_hdg) * self.airspeed + self.SSD.v_wind_x
        gs_y_DEST = np.cos(target_hdg) * self.airspeed + self.SSD.v_wind_y

        # calc resolution for each ruleset
        # for OPT
        # Calc t
        t_OPT = np.sum((np.array([self.gs_x, self.gs_y]) - p) * q, axis=1) / l2
        # Speed of boolean indices only slightly faster (negligible)
        # t must be limited between 0 and 1
        t_OPT = np.clip(t_OPT, 0., 1.)
        t_OPT[same] = 0.
        # Calculate closest point to each edge
        x1_OPT = p[:, 0] + t_OPT * q[:, 0]
        y1_OPT = p[:, 1] + t_OPT * q[:, 1]
        # Get distance squared
        d2_OPT = (x1_OPT - self.gs_x) ** 2 + (y1_OPT - self.gs_y) ** 2
        # Sort distance
        ind_OPT = np.argsort(d2_OPT)
        x1_OPT = x1_OPT[ind_OPT]
        y1_OPT = y1_OPT[ind_OPT]

        # for DEST
        # Calc t
        t_DEST = np.sum((np.array([gs_x_DEST, gs_y_DEST]) - p) * q, axis=1) / l2
        # Speed of boolean indices only slightly faster (negligible)
        # t must be limited between 0 and 1
        t_DEST = np.clip(t_DEST, 0., 1.)
        t_DEST[same] = 0.
        # Calculate closest point to each edge
        x1_DEST = p[:, 0] + t_DEST * q[:, 0]
        y1_DEST = p[:, 1] + t_DEST * q[:, 1]
        # Get distance squared
        d2_DEST = (x1_DEST - gs_x_DEST) ** 2 + (y1_DEST - gs_y_DEST) ** 2
        # Sort distance
        ind_DEST = np.argsort(d2_DEST)
        x1_DEST = x1_DEST[ind_DEST]
        y1_DEST = y1_DEST[ind_DEST]

        # for HDG
        # Calc t
        t_HDG = np.sum((np.array([self.gs_x, self.gs_y]) - p_HDG) * q_HDG, axis=1) / l2_HDG
        # Speed of boolean indices only slightly faster (negligible)
        # t must be limited between 0 and 1
        t_HDG = np.clip(t_HDG, 0., 1.)
        t_HDG[same_HDG] = 0.
        # Calculate closest point to each edge
        x1_HDG = p_HDG[:, 0] + t_HDG * q_HDG[:, 0]
        y1_HDG = p_HDG[:, 1] + t_HDG * q_HDG[:, 1]
        # Get distance squared
        d2_HDG = (x1_HDG - self.gs_x) ** 2 + (y1_HDG - self.gs_y) ** 2
        # Sort distance
        ind_HDG = np.argsort(d2_HDG)
        x1_HDG = x1_HDG[ind_HDG]
        y1_HDG = y1_HDG[ind_HDG]
        
        # draw own speed vector
        speed_vector = Arrow(0., 0., float(self.gs_x), float(self.gs_y), facecolor='black', edgecolor='black')
        ax.add_patch(speed_vector)

        # draw dartip
        darttip = shapely.affinity.scale(unit_dartipp, xfact=darttip_scale, yfact=darttip_scale, origin=(0,0))
        darttip = shapely.affinity.rotate(darttip, -np.rad2deg(self.trk), origin=(0,0))
        darttip_path = pathify(darttip)
        darttip_patch = PathPatch(darttip_path, facecolor='black', edgecolor='black')
        ax.add_patch(darttip_patch)

        # draw target symbol
        target_symbol = shapely.affinity.scale(unit_target_triangle, xfact=target_triangle_scale, yfact=target_triangle_scale, origin=(0,0))
        target_symbol = shapely.affinity.rotate(target_symbol, -np.rad2deg(target_trk), origin=(0,0))
        dx = self.v_max * np.sin(target_trk) + self.SSD.v_wind_x
        dy = self.v_max * np.cos(target_trk) + self.SSD.v_wind_y
        target_symbol = shapely.affinity.translate(target_symbol, xoff=dx, yoff=dy)
        target_symbol_path = pathify(target_symbol)
        target_symbol_patch = PathPatch(target_symbol_path, facecolor='blue', edgecolor='black')
        ax.add_patch(target_symbol_patch)

        plt.xlim(-v_range, v_range)
        plt.ylim(-v_range, v_range)
        #plt.axis('off')
        ax.set_aspect('equal')

        # plot solutions
        plt.plot(x1_OPT[0], y1_OPT[0], 'go')
        plt.plot(x1_DEST[0], y1_DEST[0], 'go')
        plt.plot(x1_HDG[0], y1_HDG[0], 'go')

        plt.show()

def pathify(polygon):
    ''' Convert coordinates to path vertices. Objects produced by Shapely's
        analytic methods have the proper coordinate order, no need to sort.

        The codes will be all "LINETO" commands, except for "MOVETO"s at the
        beginning of each subpath
    '''
    vertices = list(polygon.exterior.coords)
    codes = [Path.MOVETO if i == 0 else Path.LINETO
             for i in range(len(polygon.exterior.coords))]

    for interior in polygon.interiors:
        vertices += list(interior.coords)
        codes += [Path.MOVETO if i == 0 else Path.LINETO
                  for i in range(len(interior.coords))]

    return Path(vertices, codes)

class GeofenceSegment(object):
    def __init__(self, a, b):
        self.a = a
        self.b = b

        self.dx = self.b[0] - self.a[0]
        self.dy = self.b[1] - self.a[1]

        self.phi = np.arctan2(self.dy, self.dx)
        self.cos_phi = np.cos(self.phi)
        self.sin_phi = np.sin(self.phi)

        self.R = np.array([[self.cos_phi, -self.sin_phi], \
                           [self.sin_phi, self.cos_phi]])

        self.x_p = np.linalg.solve(self.R.T, vec_x)
        self.y_p = np.linalg.solve(self.R.T, vec_y)

class ConvexGeofence(object):
    def __init__(self):
        self.segments = []

    def add_segment(self, a, b):
        segment = GeofenceSegment(a, b)
        self.segments.append(segment)

def airspeed_to_groundsspeed(airspeed, trk, v_wind_x, v_wind_y):
    if (v_wind_x == 0) and (v_wind_y == 0):
        return airspeed, trk

    hdg = sp.Symbol('hdg')
    gs = sp.Symbol('gs')
    solution = solve([airspeed * sp.sin(hdg) + v_wind_x - gs * sp.sin(trk), airspeed * sp.cos(hdg) + v_wind_y - gs * sp.cos(trk)], dict=True)
    for i in range(len(solution)):
        if solution[i][gs] >= 0:
            gs_sol = solution[i][gs]
            hdg_sol = solution[i][hdg]
    return gs_sol, hdg_sol

def trk_to_hdg(trk, airspeed, vn_wind, ve_wind):
    trk = np.deg2rad(trk)
    # Calculate decrab angle (right positivive)
    u_right = np.array([np.cos(trk + 0.5 * np.pi), np.sin(trk + 0.5 * np.pi)]) # (north, east) Unit vector right perpenidcular to track
    wind_right = vn_wind * u_right[0] + ve_wind * u_right[1] # Wind component along u_right
    decrab_angle = np.arcsin(-wind_right / airspeed)# to right positive

    # so the hdg is the trk + decrab_angle
    hdg = np.rad2deg((trk + decrab_angle) % (2. * np.pi))

    return hdg

SSD = SSD(0., 0.)
SSD.add_SSD_entity('UAV1', -100, 0, 0, 10, 20, 5)
SSD.add_SSD_entity('UAV2', 0, 100, 270, 10, 20, 5)
SSD.add_SSD_entity('UAV3', 50, -50, 300, 7.5, 20, 5)
SSD.add_SSD_entity('UAV4', -50, -150, 320, 12.5, 20, 5)
#SSD.geofence.add_segment((-half_geofence_width, half_geofence_width), (-half_geofence_width, -half_geofence_width))
#SSD.geofence.add_segment((-half_geofence_width, -half_geofence_width), (half_geofence_width, -half_geofence_width))
#SSD.geofence.add_segment((half_geofence_width, -half_geofence_width), (half_geofence_width, half_geofence_width))
#SSD.geofence.add_segment((half_geofence_width, half_geofence_width), (-half_geofence_width, half_geofence_width))
#SSD.SSD_entities['UAV1'].plot_RV()
#SSD.SSD_entities['UAV1'].plot_VOs()
#SSD.SSD_entities['UAV1'].plot_VOs_geofences()
#SSD.SSD_entities['UAV1'].plot_SSD()
SSD.SSD_entities['UAV1'].plot_SSD_solution_points(target_trk = 90.)