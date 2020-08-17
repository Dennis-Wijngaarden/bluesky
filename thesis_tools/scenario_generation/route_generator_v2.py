import numpy as np
import parameters
import json
import random
import os

# Function declaration
def airspeed_to_groundspeed(airspeed, trk, windspeed, wind_dir):
    wind_ne = np.array([-windspeed * np.cos(wind_dir), -windspeed * np.sin(wind_dir)])

    # Calculate decrab angle (right positivive) and groundspeed
    u_trk = np.array([np.cos(trk), np.sin(trk)]) # (north, east) Unit vector of track
    u_right = np.array([np.cos(trk + 0.5 * np.pi), np.sin(trk + 0.5 * np.pi)]) # (north, east) Unit vector right perpenidcular to track
    wind_right = np.dot(wind_ne, u_right) # Wind component along u_right
    decrab_angle = np.arcsin(-wind_right / airspeed)# to right positive
    wind_trk = np.dot(u_trk, wind_ne) # Wind along track axis (u_trk)
    groundspeed = airspeed * np.cos(decrab_angle) + wind_trk # Speed of vehicle along track
    
    return groundspeed

def generate_route():

    # Check if there is a aircraft.json, otherwise stop
    if (not os.path.isfile("thesis_tools/data/aircraft.json")):
        print("first generate traffic using traffic_generator.py")
        exit()

    # Check if there is a scenario.json, otherwise stop
    if (not os.path.isfile("thesis_tools/data/scenario.json")):
        print("first generate wind using scenario_generator.py")
        exit()

    # Initialize randomizer
    random.seed()

    # Generate route dictionary
    route_data = []

    # Load scenario data from scenario.json
    scenario_json = open("thesis_tools/data/scenario.json", "r")
    scenario_data = json.load(scenario_json)
    scenario_json.close()

    for i in range(parameters.N_missions):
        route_data.append([])

        # normal vector generation of gf0 segment outwards
        gf0_side = scenario_data[i]["gf0_side"]
        if (gf0_side == "Left"):
            course_rad = np.deg2rad(scenario_data[i]['trk0']) - 0.5 * np.pi
        else:
            course_rad = np.deg2rad(scenario_data[i]['trk0']) + 0.5 * np.pi
        n_gf0 = np.array([np.sin(course_rad), np.cos(course_rad)])

        # a and u paramter of r = a + lambda u equation for gf0 segment
        a_gf0 = n_gf0 * parameters.gf_margin_dist
        u_gf0 = np.array([np.sin(np.deg2rad(scenario_data[i]['trk0'])), np.cos(np.deg2rad(scenario_data[i]['trk0']))])

        # Assume 2 vehicles
        for j in range(2):
            data_entry = {}
            data_entry['points'] = []
            data_entry['trk'] = []
            data_entry['qdr'] = []
            data_entry['dist'] = []
            
            route_done = False
            accumulated_leg_time = 0.
            if (j==0):
                trk_rad = np.deg2rad(scenario_data[i]['trk0'])
                start_point = scenario_data[i]['start_pos0_nowind']
                airspeed = scenario_data[i]['spd0']
            else:
                trk_rad = np.deg2rad(scenario_data[i]['trk1'])
                start_point = scenario_data[i]['start_pos1_nowind']
                airspeed = scenario_data[i]['spd1']
            
            k = 0

            while not route_done:
                if k == 0:
                    previous_coordinate = start_point
                else:
                    previous_coordinate = np.array(data_entry['points'][-1])

                valid_leg = False

                while not valid_leg:
                    if (k != 0):
                        trk_rad = random.uniform(0, 2. * np.pi)

                    # Check if trk is converging to predefined GF segment:
                    u_trk = np.array([np.sin(trk_rad), np.cos(trk_rad)])
                    if ((np.dot(n_gf0, u_trk) > 0) and not( k == 0 )):
                        # Caluclate maximum leg time from distance if converging
                        max_dist_gf0 = calc_dist_to_crossing(a_gf0, u_gf0, previous_coordinate, u_trk)
                        max_time_gf0 = max_dist_gf0 / airspeed

                        if (max_time_gf0 > parameters.t_leg_min):
                            leg_time = random.uniform(parameters.t_leg_min, min(max_time_gf0, parameters.t_leg_max))
                            valid_leg = True
                    
                    else:
                        leg_time = random.uniform(parameters.t_leg_min, parameters.t_leg_max)
                        valid_leg = True
                
                # Create waypoint
                leg_distance = leg_time * airspeed
                relative_coordinate = [leg_distance * np.sin(trk_rad), leg_distance * np.cos(trk_rad)] # [east, north]
                absolute_coordinate = [previous_coordinate[0] + relative_coordinate[0], previous_coordinate[1] + relative_coordinate[1]]

                data_entry['trk'].append(trk_rad)
                data_entry['points'].append(absolute_coordinate)

                # Check if accumulated time bigger than mission time
                accumulated_leg_time += leg_time
                if (accumulated_leg_time > (2. * parameters.t_mission)):
                    route_done = True
                k += 1
            
            # Add qdr and dist by looping through points
            for k in range(len(data_entry['points'])):
                data_entry['qdr'].append(np.rad2deg(np.arctan2(data_entry['points'][k][0], data_entry['points'][k][1]))) 
                data_entry['dist'].append(np.sqrt(data_entry['points'][k][0]**2 + data_entry['points'][k][1]**2))

            route_data[i].append(data_entry)

    # Write data to json files
    location = "thesis_tools/data"
    route_data_json = json.dumps(route_data, indent=4)
    route_json_file = open(location + "/route.json", "w")
    route_json_file.write(route_data_json)
    route_json_file.close()

def calc_dist_to_crossing(a_gf, u_gf, a_loc, u_loc):
    # solve using linalg
    a_gf_x = a_gf[0]
    a_gf_y = a_gf[1]

    u_gf_x = u_gf[0]
    u_gf_y = u_gf[1]

    a_loc_x = a_loc[0]
    a_loc_y = a_loc[1]

    u_loc_x = u_loc[0]
    u_loc_y = u_loc[1]

    a = np.array([[u_loc_x, -u_gf_x], [u_loc_y, -u_gf_y]])
    b = np.array([a_gf_x - a_loc_x, a_gf_y - a_loc_y])
    lambdas = np.linalg.solve(a, b)
    return lambdas[0]