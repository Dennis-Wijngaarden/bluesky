import numpy as np 
import parameters
import json
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

# Check if there is a scenario.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/scenario.json")):
    print("first generate wind using scenario_generator.py")
    exit()

# Initialize randomizer
random.seed()

# Generate route dictionary
route_data = []

# Load aircraft data from aircraft.json
aircraft_json = open("thesis_tools/data/aircraft.json", "r")
aircraft_data = json.load(aircraft_json)
aircraft_json.close()

# Load wind data from wind.json
wind_json = open("thesis_tools/data/wind.json", "r")
wind_data = json.load(wind_json)
wind_json.close()

# Function declaration
def airspeed_to_groundspeed(airspeed, hdg, windspeed, wind_dir):
    wind_ne = np.array([-windspeed * np.cos(wind_dir), -windspeed * np.sin(wind_dir)])

    # Calculate decrab angle (right positivive) and groundspeed
    u_trk = np.array([np.cos(hdg), np.sin(hdg)]) # (north, east) Unit vector of track
    u_right = np.array([np.cos(hdg + 0.5 * np.pi), np.sin(hdg + 0.5 * np.pi)]) # (north, east) Unit vector right perpenidcular to track
    wind_right = np.dot(wind_ne, u_right) # Wind component along u_right
    decrab_angle = - np.sin(wind_right / airspeed)# to right positive
    wind_trk = np.dot(u_trk, wind_ne) # Wind along track axis (u_trk)
    groundspeed = airspeed * np.cos(decrab_angle) + wind_trk # Speed of vehicle along track
    
    return groundspeed

# Load scenario data from scenario.json
scenario_json = open("thesis_tools/data/scenario.json", "r")
scenario_data = json.load(scenario_json)
scenario_json.close()

for i in range(parameters.N_missions):
    route_data.append([])
    # assume 2 vehicles:
    airspeed = [scenario_data[i]['spd0'], scenario_data[i]['spd1']] # Intial airspeeds
    windspeed = wind_data[i]['speed'] # Wind soeed [m/s]
    wind_dir = np.deg2rad(wind_data[i]['direction']) # Wind direction where the wind is coming from
    for j in range(2):
        data_entry = {}
        data_entry['hdg'] = []
        if (j==0):
            data_entry['hdg'].append(np.deg2rad(scenario_data[i]['hdg0']))
        else:
            data_entry['hdg'].append(np.deg2rad(scenario_data[i]['hdg0'] + scenario_data[i]['d_psi']))
        data_entry['points'] = []
        data_entry['gs'] = []
        data_entry['qdr'] = [] # deg
        data_entry['dist'] = []

        route_done = False
        accumulted_leg_time = 0.
        k = 0
        while (not route_done):
            # Create new heading target
            previous_coordinate = (0., 0.)
            if (k != 0):
                data_entry['hdg'].append(random.uniform(0, 2. * np.pi))
                previous_coordinate = data_entry['points'][-1]
                
            # Create waypoint
            data_entry['gs'].append(airspeed_to_groundspeed(airspeed[j], data_entry['hdg'][-1], windspeed, wind_dir))
            leg_time = random.uniform(parameters.t_leg_min, parameters.t_leg_max)
            if (j == 0):
                leg_distance = leg_time * max(data_entry['gs'][-1], scenario_data[i]['spd0'])
            else:
                leg_distance = leg_time * max(data_entry['gs'][-1], scenario_data[i]['spd1'])
            relative_coordinate = (leg_distance * np.cos(data_entry['hdg'][-1]), leg_distance * np.sin(data_entry['hdg'][-1])) # (north, east)
            absolute_coordinate = (previous_coordinate[0] + relative_coordinate[0], previous_coordinate[1] + relative_coordinate[1])
            data_entry['points'].append(absolute_coordinate)
            data_entry['qdr'].append(np.rad2deg(np.arctan2(absolute_coordinate[1], absolute_coordinate[0])))
            data_entry['dist'].append(np.sqrt(absolute_coordinate[0]**2 + absolute_coordinate[1]**2))


            # Check if accumulated time bigger than mission time
            accumulted_leg_time += leg_time
            if (accumulted_leg_time > parameters.t_mission):
                route_done = True
            k += 1

        route_data[i].append(data_entry)

# Write data to json file
location = "thesis_tools/data"
route_data_json = json.dumps(route_data, indent=4)
route_json_file = open(location + "/route.json", "w")
route_json_file.write(route_data_json)
route_json_file.close()