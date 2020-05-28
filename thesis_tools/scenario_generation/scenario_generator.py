import numpy as np
import parameters
import random
import json
import os

import route_generator as rg

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

    # Load aircraft data from aircraft.json
    aircraft_json = open("thesis_tools/data/aircraft.json", "r")
    aircraft_data = json.load(aircraft_json)
    aircraft_json.close()

    wind_data = []
    for i in range(parameters.N_missions):
        valid_wind = False
        wind_entry = {}
        data_entry = {}

        # Determine v_max_wind
        v_max_wind = None
        for j in range(len(aircraft_data[i])):
            if (j == 0):
                v_max_wind = aircraft_data[i][j]['v_max']
            else:
                v_max_wind = min(v_max_wind, aircraft_data[i][j]['v_max'])
        
        while not valid_wind:
            # Generate some random windspeed and direction
            wind_speed = random.uniform(0, v_max_wind) # Wind speed [m/s]
            wind_direction = random.uniform(0., 2. * np.pi) # Direction from which the wind is blowing

            # generate some random variables for UAV0
            # It is assumed that 2 aircraft are involved in a mission
            spd0 = None
            if (aircraft_data[i][0]['type'] == 'RC'):
                spd0 = random.uniform(max(parameters.min_vel_RC, wind_speed), aircraft_data[i][0]['v_max'])
            else:
                spd0 = random.uniform(max(aircraft_data[i][0]['v_min'], wind_speed), aircraft_data[i][0]['v_max'])

            # Make sure that rime of leg 0 is between t_leg_min and t_leg max of UAV0
            valid_trk0 = False
            counter0 = 0
            while (not valid_trk0):
                counter0 += 1
                trk0 = random.uniform(0., 360.)
                # Check gs with wind
                gs_wind = rg.airspeed_to_groundspeed(spd0, np.deg2rad(trk0), wind_speed, wind_direction)
                if (gs_wind > spd0):
                    # if max time without wind, does windy condition comply with min time?
                    trk_dist_max = spd0 * parameters.t_leg_max
                    time_min = trk_dist_max / gs_wind
                else:
                    # if max time with wind, does wind calm condition comply with min time?
                    trk_dist_max = gs_wind * parameters.t_leg_max
                    time_min = trk_dist_max / spd0
                if (time_min > parameters.t_leg_min):
                    valid_trk0 = True
                if (counter0 > 100):
                    break
            
            # Generate wind again when trk0 is not valid
            if not valid_trk0:
                continue

            # Define side of fiest geofence segment w.r.t. UAV0
            gf0_side = random.choice(["Left", "Right"])

            # generate some random variables for UAV1
            spd1 = None
            if (aircraft_data[i][1]['type'] == 'RC'):
                spd1 = random.uniform(max(parameters.min_vel_RC, wind_speed), aircraft_data[i][1]['v_max'])
            else:
                spd1 = random.uniform(max(aircraft_data[i][1]['v_min'], wind_speed), aircraft_data[i][1]['v_max'])

            # Make sure that rime of leg 0 is between t_leg_min and t_leg max of UAV1
            valid_trk1 = False
            counter1 = 0
            while (not valid_trk1):
                counter1 += 1
                if (gf0_side == "Left"):
                    d_psi = random.uniform(-180., 0.)
                else: # Right
                    d_psi = random.uniform(0., 180.)
                trk1 = trk0 + d_psi

                # Check gs with wind
                gs_wind = rg.airspeed_to_groundspeed(spd1, np.deg2rad(trk1), wind_speed, wind_direction)
                if (gs_wind > spd1):
                    # if max time without wind, does windy condition comply with min time?
                    trk_dist_max = spd1 * parameters.t_leg_max
                    time_min = trk_dist_max / gs_wind
                else:
                    # if max time with wind, does wind calm condition comply with min time?
                    trk_dist_max = gs_wind * parameters.t_leg_max
                    time_min = trk_dist_max / spd1
                if (time_min > parameters.t_leg_min):
                    valid_trk1 = True
                if (counter1 > 100):
                    break
            
            # Generate wind again when trk1 is not valid
            if not valid_trk1:
                continue
            
            # loop completed -> wind valid
            valid_wind = True
        
        data_entry['spd0'] = spd0
        data_entry['turn_rad0'] = abs(spd0**2 / (np.tan(np.deg2rad(parameters.max_bank_angle)) * 9.80665))
        data_entry['trk0'] = trk0
        data_entry['hdg0_wind'] = np.rad2deg(trk_to_hdg(np.deg2rad(trk0), spd0, wind_speed, wind_direction))
        data_entry['gf0_side'] = gf0_side
        data_entry['spd1'] = spd1
        data_entry['turn_rad1'] = abs(spd1**2 / (np.tan(np.deg2rad(parameters.max_bank_angle)) * 9.80665))
        data_entry['d_psi'] = d_psi
        data_entry['dist_cpa'] = random.uniform(0., parameters.R_PZ)
        scenario_data.append(data_entry)

        wind_entry['speed'] = wind_speed
        wind_entry['direction'] = np.rad2deg(wind_direction)
        wind_data.append(wind_entry)

        

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