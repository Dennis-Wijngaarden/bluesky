import numpy as np
import parameters
import json
import random
import os

def generate_wind():

    # Check if there is a aircraft.json, otherwise stop
    if (not os.path.isfile("thesis_tools/data/aircraft.json")):
        print("first generate traffic using traffic_generator.py")
        exit()

    # Initialize randomizer
    random.seed()

    # Create wind dictionary (index corresponds to mission number)
    wind_data = []

    # Load aircraft data from aircraft.json
    aircraft_json = open("thesis_tools/data/aircraft.json", "r")
    aircraft_data = json.load(aircraft_json)
    aircraft_json.close()

    for i in range(parameters.N_missions):
        data_entry = {}
        # Determine v_max
        v_max = None
        for j in range(len(aircraft_data[i])):
            if (j == 0):
                v_max = aircraft_data[i][j]['v_max']
            else:
                v_max = min(v_max, aircraft_data[i][j]['v_max'])

        data_entry['speed'] = random.uniform(0, v_max - parameters.min_vel_RC) # Wind speed [m/s] and margin of min_vel_RC
        data_entry['direction'] = random.uniform(0., 360.) # Direction from which the wind is blowing
        wind_data.append(data_entry)

    # Write data to json file

    location = "thesis_tools/data"
    wind_data_json = json.dumps(wind_data, indent=4)
    wind_json_file = open(location + "/wind.json", "w")
    wind_json_file.write(wind_data_json)
    wind_json_file.close()