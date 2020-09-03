import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from sympy.solvers import solve
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

n_points_v_ring = 360 # Number of points in speed rings
max_v = 100. # Scale diagrams to max v [m/s]
R_PZ = 50. # Radius of protected zone [m]

class SSD(object):
    def __init__(self, wind_spd = 0., wind_dir = 0.):
        self.wind_spd = wind_spd # [m/s]
        self.wind_dir = np.deg2rad(wind_dir) # [rad]

        self.v_wind_x = -self.wind_spd * np.sin(self.wind_dir) # wind blowing to east [m/s]
        self.v_wind_y = -self.wind_spd * np.cos(self.wind_dir) # wind blowing to north [m/s]

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
        self.airspeed = airspeed # [m/s]
        self.v_max = v_max # [m/s]
        self.v_min = v_min # [None] (if RC) / [m/s] if FW

        self.gs, self.hdg = airspeed_to_groundsspeed(self.airspeed, self.trk, self.SSD.v_wind_x, self.SSD.v_wind_y)
        self.gs_x = self.gs * np.sin(trk)
        self.gs_y = self.gs * np.cos(trk)

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

        # Now loop over all intruders (except self)
        for callsign, intruder in self.SSD.SSD_entities.items():
            # Do not add self
            if callsign == self.callsign:
                continue

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

        print(np.array(n_1)[0][0] * 2. * max_v)

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
    
    ############
    # PLOTTERS #
    ############

    def plot_VOs(self):
        VOs_polygon = self.construct_VOs_polygon()

        fig, ax = plt.subplots()

        if VOs_polygon.geom_type == 'Polygon':
            path = pathify(VOs_polygon)
            patch = PathPatch(path, facecolor='red', edgecolor='red')
            ax.add_patch(patch)
        elif VOs_polygon.geom_type == 'MultiPolygon':
            for geom in VOs_polygon.geoms:
                path = pathify(geom)
                patch = PathPatch(path, facecolor='red', edgecolor='red')
                ax.add_patch(patch)
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

def airspeed_to_groundsspeed(airspeed, trk, v_wind_x, v_wind_y):
    if (v_wind_x == 0) and (v_wind_y == 0):
        return airspeed, trk

    hdg = sp.Symbol('hdg')
    gs = sp.Symbol('gs')
    solution = solve([airspeed * sp.sin(hdg) + v_wind_x - gs * sp.sin(trk), airspeed * sp.cos(hdg) + v_wind_y - gs * sp.cos(trk)], dict=True)
    print(solution)
    for i in range(len(solution)):
        if solution[i][gs] >= 0:
            gs_sol = solution[i][gs]
            hdg_sol = solution[i][hdg]
    return gs_sol, hdg_sol


SSD = SSD(0., 0.)
SSD.add_SSD_entity('UAV1', 0, 0, 0, 10, 20, 5)
SSD.add_SSD_entity('UAV2', 500, 0, 270, 10, 20, 5)
SSD.add_SSD_entity('UAV3', 250, 250, 225, 10, 20, 5)
SSD.add_SSD_entity('UAV4', 0, -300, 30, 10, 20, 5)
SSD.SSD_entities['UAV1'].plot_VOs()