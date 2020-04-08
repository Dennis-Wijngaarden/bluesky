import numpy as np
import parameters
import json
from scipy.spatial import ConvexHull
import random
import os

class RouteSegment():
    def __init__(self, p0, p1):
        dx = p1[0] - p0[0] # delta in East [m]
        dy = p1[1] - p0[1] # delta in North [m]
        self.abs = np.sqrt(dx**2 + dy**2) # length of vector [m]
        self.u = np.array([dx, dy]) / self.abs # unit vector of segment
        self.n_right = np.array([self.u[1], -self.u[0]]) # normal vector pointing to the right 
        self.n_left = np.array([self.u[1], self.u[0]]) # Normal vector pointing to the left
        self.trk = np.arctan2(dx, dy) # radians [-pi, pi]

def intersect_angle_ranges(angles1, angles2):
    # Check if angles 2 range must be adjusted by 2pi
    diff_interval = angles2[0] - angles1[0] # diff wrt angle 1
    if (diff_interval > np.pi):
        angles2 = angles2 - np.pi
    if (diff_interval < -np.pi):
        angles2 = angles2 + np.pi

    min_angle = max(angles1[0], angles2[0])
    max_angle = min(angles1[1], angles2[1])
    return np.array([min_angle, max_angle])

def get_trk_range(trk0, trk1):
    # First assume trk0 is minimum
    trk_0 = trk0
    trk_1 = trk1
    if (trk_1 < trk_0):
        trk_1 += 2. * np.pi
    if ((trk_1 - trk_0) < np.pi):
        angle_min = trk_0
        angle_max = trk_1
        return angle_min, angle_max
    
    trk_0 = trk0
    trk_1 = trk1
    if (trk_0 < trk_1):
        trk_0 += 2. * np.pi
    angle_min = trk_1
    angle_max = trk_0
    return angle_min, angle_max

def generate_geofence(loc_route, loc_output):

    # Check if there is a route.json, otherwise stop
    if (not os.path.isfile("thesis_tools/data/route.json")):
        print("first generate route using route_generator.py")
        exit()

    # initialize randomizeer
    random.seed()

    # Create geofence dictionary (index corresponds to mission number)
    geofence_data = []

    # load route data from route.json
    route_json = open(loc_route, "r")
    route_data = json.load(route_json)
    route_json.close()

    for i in range(parameters.N_missions):
        # Add list to geofence data to store geofence points for mission
        geofence_data.append([])
        for j in range(parameters.N_vehicles):
            # Create for every vehicle a dictionary to store data
            data_entry = {}
            # add geofence points, qdr and dist to data_entry
            data_entry['points'] = []
            data_entry['qdr'] = []
            data_entry['dist'] = []

            # Determine outer points of routes CCW
            route_hull = ConvexHull(points = np.array(route_data[i][j]['points']))
            points = np.array(route_data[i][j]['points'])[route_hull.vertices].tolist()

            for k in range(len(points)):
                segments = [None, None] # Route segments of adjacent points

                # Look at adjaccent points

                # if first point
                if (k == 0):
                    segments[0] = RouteSegment(points[-1], points[k])
                    segments[1] = RouteSegment(points[1],  points[k])
                # if last point
                if (k == len(points) - 1):
                    segments[0] = RouteSegment(points[-2], points[k])
                    segments[1] = RouteSegment(points[0],  points[k])
                # else
                else:
                    segments[0] = RouteSegment(points[k - 1], points[k])
                    segments[1] = RouteSegment(points[k + 1], points[k])

                # Calc angle range
                angle_min, angle_max = get_trk_range(segments[0].trk, segments[1].trk)
                
                angle = random.uniform(angle_min, angle_max)
                radius = random.uniform(0., parameters.max_gf_dist)
                dx = radius * np.sin(angle) # dx from points[k]
                dy = radius * np.cos(angle) # dy from points[k]
                point = [points[k][0] + dx, points[k][1] + dy]
                data_entry['points'].append(point)
                data_entry['qdr'].append(np.rad2deg(np.arctan2(point[0], point[1])))
                data_entry['dist'].append(np.sqrt(point[0]**2 + point[1]**2))
            geofence_data[i].append(data_entry)

    # Write data to json file
    geofence_data_json = json.dumps(geofence_data, indent=4)
    geofence_json_file = open(loc_output, "w")
    geofence_json_file.write(geofence_data_json)
    geofence_json_file.close()       

def generate_geofences():
    generate_geofence("thesis_tools/data/route.json", "thesis_tools/data/geofence.json")
    generate_geofence("thesis_tools/data/route_wind.json", "thesis_tools/data/geofence_wind.json")