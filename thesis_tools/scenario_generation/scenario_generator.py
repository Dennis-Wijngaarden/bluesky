import numpy as np
import parameters
import random
import json
import os

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

for i in range(parameters.N_missions):
    data_entry = {}
    # It is assumed that 2 aircraft are involved in a mission
    spd0 = None
    if (aircraft_data[i][0]['type'] == 'RC'):
        spd0 = random.uniform(parameters.min_vel_RC, aircraft_data[i][0]['v_max'])
    else:
        spd0 = random.uniform(aircraft_data[i][0]['v_min'], aircraft_data[i][0]['v_max'])
    data_entry['spd0'] = spd0

    hdg0 = random.uniform(0., 360.)
    data_entry['hdg0'] = hdg0

    spd1 = None
    if (aircraft_data[i][1]['type'] == 'RC'):
        spd1 = random.uniform(parameters.min_vel_RC, aircraft_data[i][1]['v_max'])
    else:
        spd1 = random.uniform(aircraft_data[i][1]['v_min'], aircraft_data[i][1]['v_max'])
    data_entry['spd1'] = spd1

    d_psi = random.uniform(0., 360.)
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