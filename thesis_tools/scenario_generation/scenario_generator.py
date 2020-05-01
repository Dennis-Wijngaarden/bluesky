import numpy as np
import parameters
import random
import json
import os

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

    # Check if there is a wind.json, otherwise stop
    if (not os.path.isfile("thesis_tools/data/wind.json")):
        print("first generate wind using wind_generator.py")
        exit()

    # Initialize randomizer
    random.seed()

    # create dictionary with scenario data
    scenario_data = []

    # Load aircraft data from aircraft.json
    aircraft_json = open("thesis_tools/data/aircraft.json", "r")
    aircraft_data = json.load(aircraft_json)
    aircraft_json.close()

    # Load wind data from wind.json
    wind_json = open("thesis_tools/data/wind.json", "r")
    wind_data = json.load(wind_json)
    wind_json.close()

    for i in range(parameters.N_missions):
        data_entry = {}
        # It is assumed that 2 aircraft are involved in a mission
        spd0 = None
        if (aircraft_data[i][0]['type'] == 'RC'):
            spd0 = random.uniform(max(parameters.min_vel_RC, wind_data[i]['speed']), aircraft_data[i][0]['v_max'])
        else:
            spd0 = random.uniform(max(aircraft_data[i][0]['v_min'], wind_data[i]['speed']), aircraft_data[i][0]['v_max'])
        data_entry['spd0'] = spd0
        data_entry['turn_rad0'] = abs(spd0**2 / (np.tan(np.deg2rad(parameters.max_bank_angle)) * 9.80665))

        trk0 = random.uniform(0., 360.)
        data_entry['trk0'] = trk0

        hdg0_wind = np.rad2deg(trk_to_hdg(np.deg2rad(trk0), spd0, wind_data[i]['speed'], np.deg2rad(wind_data[i]['direction'])))
        data_entry['hdg0_wind'] = hdg0_wind

        gf0_side = random.choice(["Left", "Right"])
        data_entry['gf0_side'] = gf0_side

        spd1 = None
        if (aircraft_data[i][1]['type'] == 'RC'):
            spd1 = random.uniform(max(parameters.min_vel_RC, wind_data[i]['speed']), aircraft_data[i][1]['v_max'])
        else:
            spd1 = random.uniform(max(aircraft_data[i][1]['v_min'], wind_data[i]['speed']), aircraft_data[i][1]['v_max'])
        data_entry['spd1'] = spd1
        data_entry['turn_rad1'] = abs(spd1**2 / (np.tan(np.deg2rad(parameters.max_bank_angle)) * 9.80665))

        if (gf0_side == "Left"):
            d_psi = random.uniform(-180., 0.)
        else: # Right
            d_psi = random.uniform(0., 180.)
        data_entry['d_psi'] = d_psi

        dist_cpa = random.uniform(0., parameters.R_PZ)
        data_entry['dist_cpa'] = dist_cpa

        scenario_data.append(data_entry)

    # Write data to json file
    location = "thesis_tools/data"
    scenario_data_json = json.dumps(scenario_data, indent=4)
    scenario_json_file = open(location + "/scenario.json", "w")
    scenario_json_file.write(scenario_data_json)
    scenario_json_file.close()