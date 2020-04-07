import json
import parameters
import random
import os

def generate_traffic():

    # Create empty aircraft dictionary
    aircraft = {}
    # Create empty aircraft data list (each index corresponds to mission index)
    aircraft_data = []

    # Initialize ranomizer
    random.seed()

    # Initialize
    aircraft_name = "__comment"
    aircraft[aircraft_name] = {}
    aircraft[aircraft_name]["name"] = "aircraft full name"
    aircraft[aircraft_name]["n_engines"] = "number of engines"
    aircraft[aircraft_name]["engine_type"] = "number of engines"
    aircraft[aircraft_name]["mtow"] = "max take-of weight (kg)"
    aircraft[aircraft_name]["oew"] = "max take-of weight (kg)"
    aircraft[aircraft_name]["mfc"] = "max fuel capacity (L)"
    aircraft[aircraft_name]["engines"] = "[( engine: power (kW) )]"


    # Loop over number of missions to generate aircraft
    for i in range(parameters.N_missions):
        # Loop over the number of vehicles per mission:
        aircraft_data.append([])
        for j in range(parameters.N_vehicles):
            aircraft_name = "UAV_" + str(i) + "_" + str(j) # UAV_[mission number]_[vehicle number]
            aircraft[aircraft_name] = {}
            aircraft[aircraft_name]["name"] = aircraft_name

            aircraft_type = random.choice(["FW", "RC"]) # Set type for UAV
            aircraft[aircraft_name]["type"] = aircraft_type

            aircraft[aircraft_name]["n_engines"] = 4
            aircraft[aircraft_name]["engine_type"] = "TS"
            aircraft[aircraft_name]["mtow"] = 0.734
            aircraft[aircraft_name]["oew"] = 0.494
            aircraft[aircraft_name]["mfc"] = 0
            aircraft[aircraft_name]["engines"] = [["Motor-E310-1", 0.0669], ["Motor-E310-2", 0.0669], ["Motor-E310-3", 0.0669], ["Motor-E310-4", 0.0669]]
            aircraft[aircraft_name]["envelop"] = {}

            v_min = None
            v_max = None

            if (aircraft_type == "FW"):
                v_max = random.uniform(15., 25.)
                v_min = v_max - random.uniform(5., 10.) 
            else:
                v_min = 0.
                v_max = random.uniform(10., 20.)

            aircraft[aircraft_name]["envelop"]["v_min"] = v_min
            aircraft[aircraft_name]["envelop"]["v_max"] = v_max
            aircraft[aircraft_name]["envelop"]["vs_min"] = -3
            aircraft[aircraft_name]["envelop"]["vs_max"] = 5
            aircraft[aircraft_name]["envelop"]["h_max"] = 5000
            aircraft[aircraft_name]["envelop"]["d_range_max"] = 7

            data_entry = {}
            data_entry['id'] = aircraft_name
            data_entry['type'] = aircraft_type
            data_entry['v_min'] = v_min
            data_entry['v_max'] = v_max
            aircraft_data[i].append(data_entry)

    json_data = json.dumps(aircraft, indent = 4)

    aircraft_file = open("data/performance/OpenAP/rotor/aircraft_experiment.json", "w")   
    aircraft_file.write(json_data)
    aircraft_file.close()

    location = "thesis_tools/data"
    if not os.path.exists(location):
        os.makedirs(location)
    aircraft_data_json = json.dumps(aircraft_data, indent=4)
    aircraft_json_file = open(location + "/aircraft.json", "w")
    aircraft_json_file.write(aircraft_data_json)
    aircraft_json_file.close()