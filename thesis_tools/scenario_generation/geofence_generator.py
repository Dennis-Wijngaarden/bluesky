import numpy as np
import parameters
import json
import scipy as sp 
import random
import os

# Check if there is a aircraft.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/aircraft.json")):
    print("first generate traffic using traffic_generator.py")
    exit()

# Check if there is a wind.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/wind.json")):
    print("first generate wind using wind_generator.py")
    exit()

# Check if there is a route.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/route.json")):
    print("first generate wind using route_generator.py")
    exit()

# initialize randomizeer
random.seed()

# Create geofence dictionary (index corresponds to mission number)
geofence_data = []

# Load aircraft data from aircraft.json
aircraft_json = open("thesis_tools/data/aircraft.json", "r")
aircraft_data = json.load(aircraft_json)
aircraft_json.close()

# Load wind data from wind.json
wind_json = open("thesis_tools/data/wind.json", "r")
wind_data = json.load(wind_json)
wind_json.close()

# load route data from route.json
route_json = open("thesis_tools/data/route.json", "r")
route_data = json.load(route_json)
route_json.close()

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

for i in range(parameters.N_missions):
    # Add list to geofence data to store geofence points for mission
    geofence_data.append([])
    for j in range(parameters.N_vehicles):
        # Create for every vehicle a dictionary to store data
        data_entry = {}
        # add geofence points to data_entry
        data_entry['points'] = []

        # Read route points
        points = route_data[i][j]['points']
        points.insert(0, [0., 0.]) # Add origin to points at index 0

        # Create route segments
        route_segments = []
        for k in range(len(points) - 1):
            route_segment = RouteSegment(points[k], points[k + 1])
            route_segments.append(route_segment)

        # Loop over points and evaluate possible ranges of geofence points
        for k in range(len(points)):
            # First point
            if (k == 0): 
                points_right = False
                points_left = False
                for m in np.arange(2, len(points)):
                    # Check if points at right
                    if (np.dot(points[m], route_segments[0].n_right) >= 0.):
                        points_right = True
                    # Check if points at left
                    if (np.dot(points[m], route_segments[0].n_left) >= 0.):
                        points_left = True
                    
                # Set angle limits
                if (points_right and not points_left):
                    angle_min = route_segments[0].trk + np.pi
                    angle_max = route_segments[0].trk + 1.5 * np.pi
                elif (points_right and points_left):
                    angle_min = route_segments[0].trk + 0.5 * np.pi
                    angle_max = route_segments[0].trk + 1.5 * np.pi
                else: # if only points left
                    angle_min = route_segments[0].trk + 0.5 * np.pi
                    angle_max = route_segments[0].trk + np.pi 

            # Last point
            elif (k == (len(points) - 1 )):
                points_right = False
                points_left = False
                for m in range(len(points) - 2):
                    rel_point = np.array(points[m]) - np.array(points[k])
                    # Check if points at right
                    if (np.dot(rel_point, route_segments[-1].n_right) <= 0.):
                        points_right = True
                    # Check if points at left 
                    if (np.dot(rel_point, route_segments[-1].n_left) <= 0.):
                        points_left= True
                
                # Set angle limits
                if (points_right and not points_left):
                    angle_min = route_segments[-1].trk - 0.5 * np.pi
                    angle_max = route_segments[-1].trk
                elif (points_right and points_left):
                    angle_min = route_segments[-1].trk - 0.5 * np.pi
                    angle_max = route_segments[-1].trk + 0.5 * np.pi
                else: # if only points left
                    angle_min = route_segments[-1].trk 
                    angle_max = route_segments[-1].trk + 0.5 * np.pi 

            # other points
            else:
                # Evaluate from leg
                points_right = False
                points_left = False
                for m in range(len(points)):
                    if ((m == k) or (m == (k - 1))):
                        pass
                    rel_point = np.array(points[m]) - np.array(points[k])
                    # Check if points at right
                    if (np.dot(rel_point, route_segments[k - 1].n_right) <= 0.):
                        points_right = True
                    # Check if points at left
                    if (np.dot(rel_point, route_segments[k - 1].n_left) <= 0.):
                        points_left = True
                
                # Set angle limits
                if (points_right and not points_left):
                    angle_min_from = route_segments[k - 1].trk - 0.5 * np.pi
                    angle_max_from = route_segments[k - 1].trk
                elif (points_right and points_left):
                    angle_min_from = route_segments[k - 1].trk - 0.5 * np.pi
                    angle_max_from = route_segments[k - 1].trk + 0.5 * np.pi
                else: # if only points left
                    angle_min_from = route_segments[k - 1].trk 
                    angle_max_from = route_segments[k - 1].trk + 0.5 * np.pi 
                angles_from = np.array([angle_min_from, angle_max_from])

                # Evaluate to leg
                points_right = False
                points_left = False
                for m in range(len(points)):
                    if ((m == k) or (m == (k + 1))):
                        pass
                    rel_point = np.array(points[m]) - np.array(points[k])
                    # Check if points at right
                    if (np.dot(rel_point, route_segments[k].n_right) >= 0.):
                        points_right = True
                    # Check if points at left
                    if (np.dot(rel_point, route_segments[k].n_left) >= 0.):
                        points_left = True
                
                # Set angle limits
                if (points_right and not points_left):
                    angle_min_to = route_segments[k].trk + np.pi
                    angle_max_to = route_segments[k].trk + 1.5 * np.pi
                elif (points_right and points_left):
                    angle_min_to = route_segments[k].trk + 0.5 * np.pi
                    angle_max_to = route_segments[k].trk + 1.5 * np.pi
                else: # if only points left
                    angle_min_to = route_segments[k].trk + 0.5 * np.pi
                    angle_max_to = route_segments[k].trk + np.pi 
                angles_to = np.array([angle_min_to, angle_max_to])

                angle_limits = intersect_angle_ranges(angles_from, angles_to)
                angle_min = angle_limits[0]
                angle_max = angle_limits[1]
            
            angle = random.uniform(angle_min, angle_max)
            radius = random.uniform(0., parameters.max_gf_dist)
            dx = radius * np.sin(angle) # dx from points[k]
            dy = radius * np.cos(angle) # dy from points[k]
            point = [points[k][0] + dx, points[k][1] + dy]
            data_entry['points'].append(point)
            geofence_data[i].append(data_entry)

# Write data to json file
location = "thesis_tools/data"
geofence_data_json = json.dumps(geofence_data, indent=4)
geofence_json_file = open(location + "/geofence.json", "w")
geofence_json_file.write(geofence_data_json)
geofence_json_file.close()
                