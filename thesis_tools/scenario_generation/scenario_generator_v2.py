import numpy as np
import parameters
import random
import json
import os

import route_generator as rg

def calc_start_position(gs0, gs1, trk0, trk1, d_cpa, t_cpa):
    v0 = gs0 * np.matrix([[np.sin(trk0)], [np.cos(trk0)]])
    v1 = gs1 * np.matrix([[np.sin(trk1)], [np.cos(trk1)]])

    v_rel = v0 - v1
    v_rel_abs = np.sqrt(np.array(v_rel)[0][0]**2 + np.array(v_rel)[1][0]**2)
    
    # create normal vector with randomized direction (Left, Right)
    random.seed()
    if (random.random() < 0.5):
        n_v_rel = np.matrix([[-np.array(v_rel)[1][0]], [np.array(v_rel)[0][0]]]) / v_rel_abs
    else:
        n_v_rel = np.matrix([[np.array(v_rel)[1][0]], [-np.array(v_rel)[0][0]]]) / v_rel_abs

    d_rel = d_cpa * n_v_rel

    # solve lambdas for line equations of first leg
    lambdas = np.linalg.inv(np.matrix([ [np.sin(trk0), -np.sin(trk1)], \
                                        [np.cos(trk0), -np.cos(trk1)]])) * d_rel
    
    cpa_pos0 = np.array(lambdas)[0][0] * np.matrix([[np.sin(trk0)], [np.cos(trk0)]])
    cpa_pos1 = np.array(lambdas)[1][0] * np.matrix([[np.sin(trk1)], [np.cos(trk1)]])

    # Now calculate start positions t_cpa seconds back in time

    start_pos0 = np.array(cpa_pos0) - t_cpa * np.array(v0)
    start_pos1 = np.array(cpa_pos1) - t_cpa * np.array(v1)
    return start_pos0, start_pos1

def get_relative_distance(point0, point1):
    return np.sqrt(point0[0][0]**2 + point1[1][0]**2)

def trk_to_hdg(trk, airspeed, windspeed, wind_dir):
    wind_ne = np.array([-windspeed * np.cos(wind_dir), -windspeed * np.sin(wind_dir)])

    # Calculate decrab angle (right positivive)
    u_right = np.array([np.cos(trk + 0.5 * np.pi), np.sin(trk + 0.5 * np.pi)]) # (north, east) Unit vector right perpenidcular to track
    wind_right = np.dot(wind_ne, u_right) # Wind component along u_right
    decrab_angle = np.arcsin(-wind_right / airspeed)# to right positive

    # so the hdg is the trk + decrab_angle
    hdg = (trk + decrab_angle) % (2. * np.pi)

    return hdg

def generate_scenario():

    # Check if there is a aircraft.json, otherwise stop
    if (not os.path.isfile("thesis_tools/data/aircraft.json")):
        print("first generate traffic using traffic_generator.py")
        exit()

    # Initialize randomizer
    random.seed()

    # create dictionary with scenario data
    scenario_data = []

    # create dictionary for wind data
    wind_data = []

    # Load aircraft data from aircraft.json
    aircraft_json = open("thesis_tools/data/aircraft.json", "r")
    aircraft_data = json.load(aircraft_json)
    aircraft_json.close()

    for i in range(parameters.N_missions):
        valid = False
        wind_entry = {}
        data_entry = {}

        while not valid:
            # Determine v_max_wind
            v_max_wind = None
            for j in range(len(aircraft_data[i])):
                if (j == 0):
                    v_max_wind = aircraft_data[i][j]['v_max']
                else:
                    v_max_wind = min(v_max_wind, aircraft_data[i][j]['v_max'])
            
            wind_speed = random.uniform(0., v_max_wind)
            wind_direction = random.uniform(0., 2.*np.pi)

            # Generate random variables for UAV0 and UAV1
            spd0 = None
            if (aircraft_data[i][0]['type'] == 'RC'):
                spd0 = random.uniform(max(parameters.min_vel_RC, wind_speed), aircraft_data[i][0]['v_max'])
            else:
                spd0 = random.uniform(max(aircraft_data[i][0]['v_min'], wind_speed), aircraft_data[i][0]['v_max'])
            
            spd1 = None
            if (aircraft_data[i][1]['type'] == 'RC'):
                spd1 = random.uniform(max(parameters.min_vel_RC, wind_speed), aircraft_data[i][1]['v_max'])
            else:
                spd1 = random.uniform(max(aircraft_data[i][1]['v_min'], wind_speed), aircraft_data[i][1]['v_max'])

            # Define side of fiest geofence segment w.r.t. UAV0
            gf0_side = random.choice(["Left", "Right"])

            # Define trks of first legs for UAV0 and UAV1
            trk0 = random.uniform(0., 360.)
            d_psi = None
            if (gf0_side == "Left"):
                d_psi = random.uniform(-180., 0.)
            else: # Right
                d_psi = random.uniform(0., 180.)
            trk1 = trk0 + d_psi

            # Define random d_cpa for conflict
            d_cpa = random.uniform(0., parameters.R_PZ)
            t_cpa = 30. 

            # Calculate groundspeed for windy conditions
            gs0 = rg.airspeed_to_groundspeed(spd0, np.deg2rad(trk0), wind_speed, wind_direction)
            gs1 = rg.airspeed_to_groundspeed(spd1, np.deg2rad(trk1), wind_speed, wind_direction)

            # Calculate start positions for wind calm conditions
            start_pos0_nowind, start_pos1_nowind = calc_start_position(spd0, spd1, np.deg2rad(trk0), np.deg2rad(trk1), d_cpa, t_cpa)
            
            # Calculate start positions for windy conditions
            start_pos0_wind, start_pos1_wind = calc_start_position(gs0, gs1, np.deg2rad(trk0), np.deg2rad(trk1), d_cpa, t_cpa)

            ###################
            # Valiidty Checks #
            ###################

            # Check if abs(d_rel) at start of scenario is bigger than R_PZ
            dist_nowind = get_relative_distance(start_pos0_nowind, start_pos1_nowind)
            dist_wind = get_relative_distance(start_pos0_wind, start_pos1_wind)

            if ((dist_nowind > parameters.R_PZ) and (dist_wind > parameters.R_PZ)):
                valid = True

        data_entry['spd0'] = spd0
        data_entry['spd0_wind'] = gs0
        data_entry['turn_rad0'] = abs(spd0**2 / (np.tan(np.deg2rad(parameters.max_bank_angle)) * 9.80665))
        data_entry['trk0'] = trk0
        data_entry['hdg0_wind'] = trk_to_hdg(np.deg2rad(trk0), spd0, wind_speed, wind_direction)
        data_entry['gf0_side'] = gf0_side
        data_entry['spd1'] = spd1
        data_entry['spd1_wind'] = gs1
        data_entry['turn_rad1'] = abs(spd1**2 / (np.tan(np.deg2rad(parameters.max_bank_angle)) * 9.80665))
        data_entry['trk1'] = trk1
        data_entry['hdg1_wind'] = trk_to_hdg(np.deg2rad(trk1), spd1, wind_speed, wind_direction)
        data_entry['dist_cpa'] = d_cpa
        data_entry['start_pos0_nowind'] = [start_pos0_nowind[0][0], start_pos0_nowind[1][0]]
        data_entry['start_pos1_nowind'] = [start_pos1_nowind[0][0], start_pos1_nowind[1][0]]
        data_entry['start_qdr0_nowind'] = np.rad2deg(np.arctan2(start_pos0_nowind[0][0], start_pos0_nowind[1][0]))
        data_entry['start_qdr1_nowind'] = np.rad2deg(np.arctan2(start_pos1_nowind[0][0], start_pos1_nowind[1][0]))
        data_entry['start_dist0_nowind'] = np.sqrt(start_pos0_nowind[0][0]**2 + start_pos0_nowind[1][0]**2)
        data_entry['start_dist1_nowind'] = np.sqrt(start_pos1_nowind[0][0]**2 + start_pos1_nowind[1][0]**2)
        data_entry['start_pos0_wind'] = [start_pos0_wind[0][0], start_pos0_wind[1][0]]
        data_entry['start_pos1_wind'] = [start_pos1_wind[0][0], start_pos1_wind[1][0]]
        data_entry['start_qdr0_wind'] = np.rad2deg(np.arctan2(start_pos0_wind[0][0], start_pos0_wind[1][0]))
        data_entry['start_qdr1_wind'] = np.rad2deg(np.arctan2(start_pos1_wind[0][0], start_pos1_wind[1][0]))
        data_entry['start_dist0_wind'] = np.sqrt(start_pos0_wind[0][0]**2 + start_pos0_wind[1][0]**2)
        data_entry['start_dist1_wind'] = np.sqrt(start_pos1_wind[0][0]**2 + start_pos1_wind[1][0]**2)
        scenario_data.append(data_entry)

        wind_entry['speed'] = wind_speed
        wind_entry['direction'] = np.rad2deg(wind_direction)
        wind_data.append(wind_entry)

    #print(scenario_data)
    # Write data to json files
    location = "thesis_tools/data"
    scenario_data_json = json.dumps(scenario_data, indent=4)
    scenario_json_file = open(location + "/scenario.json", "w")
    scenario_json_file.write(scenario_data_json)
    scenario_json_file.close()

    wind_data_json = json.dumps(wind_data, indent=4)
    wind_json_file = open(location + "/wind.json", "w")
    wind_json_file.write(wind_data_json)
    wind_json_file.close()