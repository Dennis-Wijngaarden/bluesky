import os
import traffic_generator
import scenario_generator
import wind_generator

traffic_generator.generate_traffic()
scenario_generator.generate_scenario()
wind_generator.generate_wind()
os.system("python thesis_tools/scenario_generation/route_generator.py")
os.system("python thesis_tools/scenario_generation/geofence_generator.py")
