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

    # Check if there is a wind.json, otherwise stop
    if (not os.path.isfile("thesis_tools/data/wind.json")):
        print("first generate wind using wind_generator.py")
        exit()

    # Check if there is a scenario.json, otherwise stop
    if (not os.path.isfile("thesis_tools/data/scenario.json")):
        print("first generate wind using scenario_generator.py")
        exit()

    # Initialize randomizer
    random.seed()

    # Generate route dictionaries, alo one with wind included
    route_data_wind_calm = []
    route_data_wind = []

    # Load wind data from wind.json
    wind_json = open("thesis_tools/data/wind.json", "r")
    wind_data = json.load(wind_json)
    wind_json.close()

    # Load scenario data from scenario.json
    scenario_json = open("thesis_tools/data/scenario.json", "r")
    scenario_data = json.load(scenario_json)
    scenario_json.close()

    for i in range(parameters.N_missions):
        route_data_wind_calm.append([])
        route_data_wind.append([])
        
        airspeed = [scenario_data[i]['spd0'], scenario_data[i]['spd1']]

        windspeed = wind_data[i]['speed'] # Wind soeed [m/s]
        wind_dir = np.deg2rad(wind_data[i]['direction']) # Wind direction where the wind is coming from

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

        # assume 2 vehicles:
        for j in range(2):
            data_entry = {} # Constant for windy and wind calm conditions
            data_entry_wind_calm = {} # Specific route parameters for wind calm conditions
            data_entry_wind_calm['points'] = []
            data_entry_wind = {}  # Specific route paramaters for windy conditions
            data_entry_wind['points'] = []
            data_entry['trk'] = []
            if (j==0):
                trk_rad = np.deg2rad(scenario_data[i]['trk0'])
                data_entry_wind_calm['points'].append([0., 0.])
                data_entry_wind['points'].append([0., 0.])
            else:
                trk_rad = np.deg2rad(scenario_data[i]['trk0'] + scenario_data[i]['d_psi'])
                # Determine start point of conflicting vehicle for wind and wind calm conditions
                # without wind (gs = airspeed)
                v_rel_e = scenario_data[i]['spd0'] * np.sin(np.deg2rad(scenario_data[i]['trk0'])) - scenario_data[i]['spd1'] * np.sin(trk_rad)
                v_rel_n = scenario_data[i]['spd0'] * np.cos(np.deg2rad(scenario_data[i]['trk0'])) - scenario_data[i]['spd1'] * np.cos(trk_rad)
                v_rel = np.sqrt(v_rel_n**2 + v_rel_e **2)

                d_rel_cpa = (parameters.t_la + parameters.t_extra) * v_rel + (0. if scenario_data[i]['dist_cpa'] > parameters.R_PZ else np.sqrt(parameters.R_PZ**2 - scenario_data[i]['dist_cpa']**2))
                dist = np.sqrt(d_rel_cpa**2 + scenario_data[i]['dist_cpa']**2)

                rd = d_rel_cpa / dist
                rx = scenario_data[i]['dist_cpa'] / dist
                qdr_rad = np.arctan2(-rx * v_rel_n + rd * v_rel_e, rd * v_rel_n + rx * v_rel_e)
                data_entry_wind_calm['points'].append([np.sin(qdr_rad) * dist, np.cos(qdr_rad) * dist])
                
                # with wind
                v_own = airspeed_to_groundspeed(scenario_data[i]['spd1'], trk_rad, windspeed, wind_dir)
                v_int = airspeed_to_groundspeed(scenario_data[i]['spd0'], np.deg2rad(scenario_data[i]['trk0']), windspeed, wind_dir)
                v_rel_e = v_int * np.sin(np.deg2rad(scenario_data[i]['trk0'])) - v_own * np.sin(trk_rad)
                v_rel_n = v_int * np.cos(np.deg2rad(scenario_data[i]['trk0'])) - v_own * np.cos(trk_rad)
                v_rel = np.sqrt(v_rel_n**2 + v_rel_e **2)

                d_rel_cpa = (parameters.t_la + parameters.t_extra) * v_rel + (0. if scenario_data[i]['dist_cpa'] > parameters.R_PZ else np.sqrt(parameters.R_PZ**2 - scenario_data[i]['dist_cpa']**2))
                dist = np.sqrt(d_rel_cpa**2 + scenario_data[i]['dist_cpa']**2)

                rd = d_rel_cpa / dist
                rx = scenario_data[i]['dist_cpa'] / dist
                qdr_rad = np.arctan2(-rx * v_rel_n + rd * v_rel_e, rd * v_rel_n + rx * v_rel_e)
                data_entry_wind['points'].append([np.sin(qdr_rad) * dist, np.cos(qdr_rad) * dist])

            # Declare other data entries
            data_entry_wind_calm['gs'] = []
            data_entry_wind['gs'] = []
            data_entry_wind_calm['qdr'] = [] # deg
            data_entry_wind['qdr'] = [] # deg
            data_entry_wind_calm['dist'] = [] 
            data_entry_wind['dist'] = []

            route_done = False
            accumulted_leg_time = 0.
            k = 0
            while (not route_done):
                # Create new target
                previous_coordinate_wind_calm = np.array(data_entry_wind_calm['points'][-1])
                previous_coordinate_wind = np.array(data_entry_wind['points'][-1])

                valid_leg = False
                while (not valid_leg):
                    if k != 0:
                        trk_rad = random.uniform(0, 2. * np.pi)

                    gs_wind = airspeed_to_groundspeed(airspeed[j], trk_rad, windspeed, wind_dir)

                    # Check if trk is converging to predefined GF segment:
                    u_trk = np.array([np.sin(trk_rad), np.cos(trk_rad)])
                    if ((np.dot(n_gf0, u_trk) > 0) and not( k == 0 and j == 0)):
                        # Caluclate maximum leg time from distance if converging
                        # for wind calm conditions
                        max_dist_gf0 = calc_dist_to_crossing(a_gf0, u_gf0, previous_coordinate_wind_calm, u_trk)
                        max_time_gf0 = max_dist_gf0 / airspeed[j]

                        # for windy conditions
                        max_dist_gf0_wind = calc_dist_to_crossing(a_gf0, u_gf0, previous_coordinate_wind, u_trk)
                        max_time_gf0_wind = max_dist_gf0_wind / gs_wind

                        # maximum possible leg time such that gf0 segment will not be crossed
                        max_leg_time_gf0 = min(max_time_gf0, max_time_gf0_wind)

                        if (max_leg_time_gf0 > parameters.t_leg_min):
                            if (max_leg_time_gf0 > parameters.t_leg_max):
                                leg_time = random.uniform(parameters.t_leg_min, parameters.t_leg_max)
                            else:
                                leg_time = random.uniform(parameters.t_leg_min, max_leg_time_gf0)
                            valid_leg = True

                    else:
                        leg_time = random.uniform(parameters.t_leg_min, parameters.t_leg_max)
                        valid_leg = True
                        
                # Create waypoint
                if (j == 0):
                    leg_distance = leg_time * max(gs_wind, scenario_data[i]['spd0'])
                else:
                    leg_distance = leg_time * max(gs_wind, scenario_data[i]['spd1'])
                relative_coordinate = [leg_distance * np.sin(trk_rad), leg_distance * np.cos(trk_rad)] # [east, north]
                absolute_coordinate_wind_calm = [previous_coordinate_wind_calm[0] + relative_coordinate[0], previous_coordinate_wind_calm[1] + relative_coordinate[1]]
                absolute_coordinate_wind = [previous_coordinate_wind[0] + relative_coordinate[0], previous_coordinate_wind[1] + relative_coordinate[1]]

                data_entry['trk'].append(trk_rad)
                data_entry_wind_calm['gs'].append(airspeed[j])
                data_entry_wind['gs'].append(gs_wind)
                data_entry_wind_calm['points'].append(absolute_coordinate_wind_calm)
                data_entry_wind['points'].append(absolute_coordinate_wind)

                # Check if accumulated time bigger than mission time
                accumulted_leg_time += leg_time
                if (accumulted_leg_time > parameters.t_mission):
                    route_done = True
                k += 1

            # Add qdr and dist by looping through points
            for k in range(len(data_entry_wind['points'])):
                data_entry_wind_calm['qdr'].append(np.rad2deg(np.arctan2(data_entry_wind_calm['points'][k][0], data_entry_wind_calm['points'][k][1])))
                data_entry_wind['qdr'].append(np.rad2deg(np.arctan2(data_entry_wind['points'][k][0], data_entry_wind['points'][k][1])))
                data_entry_wind_calm['dist'].append(np.sqrt(data_entry_wind_calm['points'][k][0]**2 + data_entry_wind_calm['points'][k][1]**2))
                data_entry_wind['dist'].append(np.sqrt(data_entry_wind['points'][k][0]**2 + data_entry_wind['points'][k][1]**2))
            data_entry_wind_calm.update(data_entry)
            data_entry_wind.update(data_entry)
            route_data_wind_calm[i].append(data_entry_wind_calm)
            route_data_wind[i].append(data_entry_wind)

    # Write data to json files
    location = "thesis_tools/data"
    route_data_json = json.dumps(route_data_wind_calm, indent=4)
    route_json_file = open(location + "/route.json", "w")
    route_json_file.write(route_data_json)
    route_json_file.close()

    location = "thesis_tools/data"
    route_data_wind_json = json.dumps(route_data_wind, indent=4)
    route_wind_json_file = open(location + "/route_wind.json", "w")
    route_wind_json_file.write(route_data_wind_json)
    route_wind_json_file.close()

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